#===============================================================================
# ENTICES: Enceladus Tidally-Induced Circulation and Exchange Simulator
# Authors: Mohit Melwani Daswani & Lucas Fifer
# Last updated: 2026-09-08
#
# Version: 1.0.0
# Release Date: 2025
#
# A computational framework for modeling tidally-driven hydrothermal 
# (hydrotidal?) circulation and water-rock reactions in ocean worlds. Initially 
# developed for Enceladus but generalizable to other planetary bodies.
#
# The model integrates:
# - Two-layer gravity field calculations
# - Tidal stress modeling on orbital timescales
# - Darcy flow in porous media
# - Hydrothermal circulation timescales
# - Water-rock reaction modeling via PHREEQC
#
# Primary references:
# - Fisher et al. (2024) JGR Planets
# - Choblet et al. (2017) Nature Astronomy 
# - Liao et al. (2020) JGR Planets. #update using Bagheri new paper??
# - Zandanel et al. (2021) Icarus
#===============================================================================

#===============================================================================
# Required libraries 
#===============================================================================
library(tidyverse)
library(conflicted)
library(dplyr)
library(patchwork)  # For plot layouts
library(glue)       # For templating
library(R6)

# IOM (organic decomposition) channel definitions live in iom_module.R, kept
# separate so circulation-physics edits here can't touch organics calibration
# and vice versa. Run from the repo root so this relative path resolves, same
# convention as the ENTICES_v3.R / ENTICES_LF_* sourcing in ENTICES_runs.Rmd.
source("iom_module.R")

# Primary rock composition, mol per unit total_rock (i.e. per kg anhydrous rock).
# Single source for all four generators that build primary-mineral inventories;
# they previously each carried their own copy of these numbers.
# Source: RECONSTRUCTION_STATUS.md §1 (final 2026-09-15; Lodders et al. 2025
# Table 4, Scenario R, schreibersite included), wt% x 10 / molar mass.
# Deliberately NOT renormalised: magnesiochromite (0.713 wt%, 0.037079 mol/kg)
# is treated as inert, so it keeps its share of the mass but is not added to
# PHREEQC. (The previous coefficients dropped it and renormalised, inflating
# every other phase by +0.73%.)
# Lawrencite (FeCl2) is intentional (Cl in the rock); its Fe is taken out of
# the metal, so fe = 2.361715 - 0.011250. Total Fe is then 4.8859 mmol/g, which
# matches the reconstruction.
# Schreibersite (Fe3P, 0.047101 mol/kg) is listed here but is wired in
# separately; see the generators.
ROCK_MOL_PER_KG <- c(
  enstatite     = 2.275648,
  forsterite    = 1.606974,
  troilite      = 2.382864,
  albite        = 0.318202,
  diopside      = 0.274529,
  anorthite     = 0.062148,
  kfeldspar     = 0.020515,
  tephroite     = 0.025995,
  fe            = 2.350465,
  ni            = 0.280953,
  lawrencite    = 0.011250,
  schreibersite = 0.047101
)

# Schreibersite (Fe3P) -- kinetic primary, defined once here and spliced into the
# kinetic generators and both restart paths.
# Reaction (RECONSTRUCTION_STATUS.md §4, finalised 2026-09-14):
#     Fe3P + 2.75 O2 + 4 H+ = 3 Fe+2 + HPO4-2 + 1.5 H2O
# Neither database has an Fe3P phase, so the KINETICS entry adds the elements via
# -formula Fe 3 P 1. PHREEQC then speciates them as Fe+2 and HPO4-2, and the
# 11 e- per Fe3P go to COUPLED H2 (5.5 H2 per Fe3P). This is the same mass
# transfer as the reaction above, and it keeps schreibersite's H2 on the same
# books as metal-corrosion H2, as §4 [D] requires (coupled O2/H+, NOT Oxg/Hdg).
# Rate: PLACEHOLDER. Bulk first-order, 1.184 %/day (range 0.59-1.78 %/day) =
# 1.370e-7 s^-1, from §4/§8 item 1b. That value comes from one 1-day H2
# measurement at pH 7, 20 C in Ar-purged DI water. It has no T, pH, W:R or
# surface-area dependence, and there is no SR() term because there is no phase
# to compute one from. Treat any schreibersite-dependent output as provisional.
SCHREIBERSITE_MW <- 198.509
SCHREIBERSITE_RATE_BODY <- paste(
  "\t-start",
  "\t1   REM PLACEHOLDER bulk first-order rate, 1.184 %/day (0.59-1.78), RECONSTRUCTION_STATUS 8(1b)",
  "\t2   REM one 1-day lab measurement, pH 7, 20 C, Ar-purged DI water, no T/pH/SSA dependence",
  "\t10  k = 1.370e-7",
  "\t20  rate = k * M",
  "\t30  IF (M <= 0) THEN rate = 0",
  "\t40  moles = rate * TIME",
  "\t50  SAVE moles",
  "\t-end", sep = "\n")
SCHREIBERSITE_RATE <- paste0("Schreibersite\n", SCHREIBERSITE_RATE_BODY)
schreibersite_kinetics_entry <- function(m0) {
  sprintf("   Schreibersite\n      -formula Fe 3 P 1\n      -m0 %.6e\n      -step_divide 1", m0)
}

#===============================================================================
# Universal Constants and Fluid Properties
#===============================================================================
CONSTANTS <- list(
  # Universal constants
  G = 6.67430e-11,  # Gravitational constant [m³/(kg*s²)]
  
  # Fluid properties (water at 0°C)
  FLUID = list(
    density = 1000,      # [kg/m³]
    viscosity = 1.8e-3,  # [Pa·s]
    heat_capacity = 4186 # [J/(kg·K)]
  )
)
  
#===============================================================================
# PlanetProfile temperature utility
#===============================================================================

# Read a PlanetProfile output text file and return the mean temperature (°C)
# of the porous-rock layer (phase ID 50) from the core top down to
# `layer_thickness` metres depth.  This is used to set SOLUTION temps in
# the PHREEQC template.
#
# Args:
#   profile_file     : path to the PlanetProfile *Profile*.txt output file
#   layer_thickness  : permeable layer thickness [m]  (same as Simulator$layer_thickness)
#
# Returns: mean temperature in °C over the specified layer depth
mean_layer_temp_C <- function(profile_file, layer_thickness) {
  lines      <- readLines(profile_file, warn = FALSE)
  header_end <- grep("P \\(MPa\\)", lines)[1]

  col_names <- c("P_MPa", "T_K", "r_m", "phase_id",
                 "rho", "Cp", "alpha", "g", "phi", "sigma", "k_therm",
                 "VP", "VS", "QS", "KS", "GS", "Ppore",
                 "rhoMatrix", "rhoPore", "MLayer", "VLayer", "Htidal")

  dat <- read.table(
    text      = paste(lines[(header_end + 1):length(lines)], collapse = "\n"),
    header    = FALSE,
    col.names = col_names
  )

  rock <- dat[dat$phase_id == 50, ]
  if (nrow(rock) == 0) stop("No phase-50 rows found in profile file.")

  r_top      <- max(rock$r_m)
  layer_rows <- rock[rock$r_m >= r_top - layer_thickness, ]
  if (nrow(layer_rows) == 0) {
    warning("layer_thickness smaller than one profile step; using single nearest row.")
    layer_rows <- rock[which.max(rock$r_m), ]
  }

  mean(layer_rows$T_K) - 273.15
}

# Same as mean_layer_temp_C but returns mean pressure in MPa over the layer.
# PHREEQC's `pressure` keyword uses atm: divide by 0.101325 to convert.
mean_layer_pressure_MPa <- function(profile_file, layer_thickness) {
  lines      <- readLines(profile_file, warn = FALSE)
  header_end <- grep("P \\(MPa\\)", lines)[1]

  col_names <- c("P_MPa", "T_K", "r_m", "phase_id",
                 "rho", "Cp", "alpha", "g", "phi", "sigma", "k_therm",
                 "VP", "VS", "QS", "KS", "GS", "Ppore",
                 "rhoMatrix", "rhoPore", "MLayer", "VLayer", "Htidal")

  dat <- read.table(
    text      = paste(lines[(header_end + 1):length(lines)], collapse = "\n"),
    header    = FALSE,
    col.names = col_names
  )

  rock <- dat[dat$phase_id == 50, ]
  if (nrow(rock) == 0) stop("No phase-50 rows found in profile file.")

  r_top      <- max(rock$r_m)
  layer_rows <- rock[rock$r_m >= r_top - layer_thickness, ]
  if (nrow(layer_rows) == 0) {
    warning("layer_thickness smaller than one profile step; using single nearest row.")
    layer_rows <- rock[which.max(rock$r_m), ]
  }

  mean(layer_rows$P_MPa)
}

# Returns the pressure at the ocean-rock interface (top of phase 50) in MPa.
# Used to set the ocean SOLUTION pressure in the PHREEQC template.
ocean_floor_pressure_MPa <- function(profile_file) {
  lines      <- readLines(profile_file, warn = FALSE)
  header_end <- grep("P \\(MPa\\)", lines)[1]

  col_names <- c("P_MPa", "T_K", "r_m", "phase_id",
                 "rho", "Cp", "alpha", "g", "phi", "sigma", "k_therm",
                 "VP", "VS", "QS", "KS", "GS", "Ppore",
                 "rhoMatrix", "rhoPore", "MLayer", "VLayer", "Htidal")

  dat <- read.table(
    text      = paste(lines[(header_end + 1):length(lines)], collapse = "\n"),
    header    = FALSE,
    col.names = col_names
  )

  rock <- dat[dat$phase_id == 50, ]
  if (nrow(rock) == 0) stop("No phase-50 rows found in profile file.")

  rock$P_MPa[which.max(rock$r_m)]
}

#===============================================================================
# Configuration Class for Enceladus Parameters
#===============================================================================
EnceladusConfig <- R6::R6Class("EnceladusConfig",
                               public = list(
                                 # Basic dimensional parameters
                                 radius = 252.1e3,        # Surface radius [m]
                                 core_radius = 194.2e3,   # Core radius [m] 
                                 mass = 1.0802e20,        # Total mass [kg] 
                                 
                                 # Material properties
                                 core_density = 2353,     # Core bulk density [kg/m³] #This and above parameters from 4 papers?
                                 core_porosity = 0.32,    # Core porosity [volume fraction] #Default assumed value from Stycinski et al 2023
                                 
                                 # Orbital properties
                                 orbital_period = 118800, # [s]
                                 max_tidal_stress = 1e5,  # [Pa]
                                 
                                 # Optional/derived parameters
                                 hydrosphere_thickness = NA, 
                                 ice_thickness = NA,
                                 ocean_thickness = NA,
                                 
                                 initialize = function(params = list()) {
                                   if(length(params) > 0) {
                                     for(param in names(params)) {
                                       if(param %in% names(self)) {
                                         self[[param]] <- params[[param]]
                                       }
                                     }
                                   }
                                   self$calculate_dimensions()
                                 },
                                 
                                 calculate_dimensions = function() {
                                   if(is.na(self$hydrosphere_thickness)) {
                                     if(is.na(self$ice_thickness) || is.na(self$ocean_thickness)) {
                                       self$hydrosphere_thickness <- self$radius - self$core_radius
                                       message("Using default hydrosphere thickness from radii difference: ", 
                                               sprintf("%.1f km", self$hydrosphere_thickness/1000))
                                     } else {
                                       self$hydrosphere_thickness <- self$ice_thickness + self$ocean_thickness
                                       message("Calculated hydrosphere thickness from ice and ocean components: ",
                                               sprintf("%.1f km", self$hydrosphere_thickness/1000))
                                     }
                                   }
                                   
                                   if(is.na(self$ice_thickness) && is.na(self$ocean_thickness)) {
                                     default_ice_fraction <- 0.36
                                     self$ice_thickness <- self$hydrosphere_thickness * default_ice_fraction
                                     self$ocean_thickness <- self$hydrosphere_thickness * (1 - default_ice_fraction)
                                     message(sprintf("Split hydrosphere using default ice fraction (%.2f)", 
                                                     default_ice_fraction))
                                   } else if(is.na(self$ice_thickness)) {
                                     self$ice_thickness <- self$hydrosphere_thickness - self$ocean_thickness
                                     message("Calculated ice thickness as remainder")
                                   } else if(is.na(self$ocean_thickness)) {
                                     self$ocean_thickness <- self$hydrosphere_thickness - self$ice_thickness
                                     message("Calculated ocean thickness as remainder")
                                   }
                                   
                                   self$validate_dimensions()
                                 },
                                 
                                 validate_dimensions = function() {
                                   if(any(is.na(c(self$hydrosphere_thickness, 
                                                  self$ice_thickness, 
                                                  self$ocean_thickness)))) {
                                     stop("Failed to calculate all dimensions")
                                   }
                                   
                                   abs_tol <- 1  # 1 meter tolerance for floating point comparison
                                   
                                   if(abs(self$hydrosphere_thickness - (self$ice_thickness + self$ocean_thickness)) > abs_tol) {
                                     stop(sprintf("Inconsistent dimensions: hydrosphere (%.1f km) ≠ ice (%.1f km) + ocean (%.1f km)",
                                                  self$hydrosphere_thickness/1000,
                                                  self$ice_thickness/1000,
                                                  self$ocean_thickness/1000))
                                   }
                                   
                                   if(abs(self$radius - (self$core_radius + self$hydrosphere_thickness)) > abs_tol) {
                                     stop(sprintf("Inconsistent dimensions: radius (%.1f km) ≠ core (%.1f km) + hydrosphere (%.1f km)",
                                                  self$radius/1000,
                                                  self$core_radius/1000,
                                                  self$hydrosphere_thickness/1000))
                                   }
                                   
                                   if(any(c(self$ice_thickness, self$ocean_thickness, 
                                            self$hydrosphere_thickness, self$core_radius) <= 0)) {
                                     stop("All dimensions must be positive")
                                   }
                                   
                                   message("\nValidated Enceladus dimensions:")
                                   message(sprintf("Total radius: %.1f km", self$radius/1000))
                                   message(sprintf("Core radius: %.1f km", self$core_radius/1000))
                                   message(sprintf("Hydrosphere: %.1f km", self$hydrosphere_thickness/1000))
                                   message(sprintf("  - Ice shell: %.1f km", self$ice_thickness/1000))
                                   message(sprintf("  - Ocean: %.1f km", self$ocean_thickness/1000))
                                 },
                                 
                                 calculate_gravity_profile = function(step_size = 100) {
                                   depth_m <- seq(0, self$radius, by = step_size)
                                   depth_radius <- self$radius - depth_m
                                   gravity_ms2 <- numeric(length(depth_radius))
                                                                      for(i in seq_along(gravity_ms2)) {
                                     r <- depth_radius[i]
                                     if(r <= self$core_radius) {
                                       gravity_ms2[i] <- (4/3) * pi * CONSTANTS$G * self$core_density * r
                                     } else {
                                       gravity_ms2[i] <- (4/3) * pi * CONSTANTS$G * (
                                         self$core_density * self$core_radius^3 + 
                                           CONSTANTS$FLUID$density * (r^3 - self$core_radius^3) #does it matter that theres no accounting for density difference between ice and water?
                                       ) / r^2
                                     }
                                   }
                                   
                                   return(data.frame(
                                     radius_m = depth_radius,
                                     depth_m = depth_m,
                                     gravity_ms2 = gravity_ms2
                                   ))
                                 },
                                 
                                 calculate_gravity = function(r) { #why is this a different function than calculate_gravity_profile?
                                   if(r <= self$core_radius) {
                                     return((4/3) * pi * CONSTANTS$G * self$core_density * r)
                                   } else {
                                     return((4/3) * pi * CONSTANTS$G * (
                                       self$core_density * self$core_radius^3 + 
                                         CONSTANTS$FLUID$density * (r^3 - self$core_radius^3)
                                     ) / r^2)
                                   }
                                 },
                                 
                                 get_derived_properties = function() {
                                   core_volume <- (4/3) * pi * self$core_radius^3
                                   ocean_volume <- (4/3) * pi * 
                                     ((self$radius - self$ice_thickness)^3 - self$core_radius^3)
                                   
                                   ocean_mass <- ocean_volume * CONSTANTS$FLUID$density
                                   # self$core_density is the Cassini-gravity-derived BULK density (matrix + pore
                                   # fluid combined; e.g. Cadek+16, Hemingway & Mittal 2019, Liao+20/21), related
                                   # to grain density by bulk = (1-porosity)*grain + porosity*fluid_density. So
                                   # solid rock mass = volume*(bulk - porosity*fluid_density), NOT
                                   # volume*bulk*(1-porosity), which double-subtracts the pore-fluid contribution.
                                   core_mass <- core_volume * (self$core_density - self$core_porosity * CONSTANTS$FLUID$density)
                                   
                                   gravity_seafloor <- self$calculate_gravity(self$core_radius)
                                   seafloor_pressure <- CONSTANTS$FLUID$density * gravity_seafloor * self$ocean_thickness #should this not add ice too? or depend on assumptions about isostasy?
                                   
                                   gravity_surface <- self$calculate_gravity(self$radius)
                                   water_rock_ratio <- ocean_mass/core_mass
                                   core_surface_area <- 4 * pi * self$core_radius^2
                                   
                                   message("\nGravity calculations:")
                                   message(sprintf("Surface gravity: %.3f m/s^2", gravity_surface))
                                   message(sprintf("Seafloor gravity: %.3f m/s^2", gravity_seafloor))
                                   
                                   return(list(
                                     core_volume = core_volume,
                                     ocean_volume = ocean_volume,
                                     ocean_mass = ocean_mass,
                                     core_mass = core_mass,
                                     gravity_surface = gravity_surface,
                                     gravity_seafloor = gravity_seafloor,
                                     seafloor_pressure = seafloor_pressure,
                                     water_rock_ratio = water_rock_ratio,
                                     core_surface_area = core_surface_area,
                                     gravity_profile = self$calculate_gravity_profile()
                                   ))
                                 },
                                 
                                 plot_gravity_profile = function() {
                                   profile <- self$calculate_gravity_profile()
                                   
                                   ggplot(profile, aes(x = radius_m/1000, y = gravity_ms2)) +
                                     geom_line() +
                                     geom_vline(xintercept = self$core_radius/1000, 
                                                color = "red", linetype = "dashed") +
                                     labs(x = "Radius (km)",
                                          y = "Gravity (m/s^2)",
                                          title = "Gravity vs. Radius in Enceladus") +
                                     theme_minimal()
                                 }
                               )
)

#===============================================================================
# Utility Functions Module
#===============================================================================
GeometricUtils <- list(
  # Calculate cross-sectional area for flow
  calculate_flow_area = function(layer_thickness, core_radius) {
    return(layer_thickness * core_radius) #confusing - why multiply a depth by a depth?
  },
  
  # Calculate vertical area for volume calculations  
  calculate_vertical_area = function(core_radius) {
    return(pi * core_radius^2) #unsure what vertical area means here
  },
  
  # Calculate depth-dependent porosity
  calculate_mean_porosity = function(surface_porosity, depth_factor = 0.8) { #why depth_factor 0.8?
    base_porosity <- surface_porosity * depth_factor
    return((surface_porosity + base_porosity) / 2)
  }
)

FlowUtils <- list(
  # Calculate Darcy flow
  calculate_darcy_flow = function(k, grad_P, area, porosity,  #this function doesnt get used?
                                  viscosity = CONSTANTS$FLUID$viscosity) {
    q <- -(k/viscosity) * grad_P  # Specific discharge
    Q <- q * area * porosity      # Volumetric flow rate
    return(Q)
  },
  
  # Calculate pressure gradient 
  calculate_pressure_gradient = function(stress, thickness) {
    return(stress/thickness)
  },
  
  # Calculate circulation time
  calculate_circulation_time = function(flow_rate, ocean_mass,
                                        fluid_density = CONSTANTS$FLUID$density) {
    if(flow_rate <= 0 || ocean_mass <= 0) {
      warning("Invalid flow rate or ocean mass")
      return(NA)
    }
    
    t_circ <- ocean_mass / (flow_rate * fluid_density)  # seconds
    t_circ_years <- t_circ / (365.25 * 24 * 3600)      # years
    return(t_circ_years)
  }
)

#===============================================================================
# Core Simulation Module
#===============================================================================
Simulator <- R6::R6Class("HydrothermalSimulator",
                         public = list(
                           # Class fields
                           params = NULL,
                           results = NULL,
                           config = NULL, 
                           
                           # Initialize simulator with parameters and configuration
                           initialize = function(k, layer_thickness, porosity = 0.32, config = NULL) {
                             self$params <- list(
                               k = k,
                               layer_thickness = layer_thickness,
                               porosity = porosity
                             )
                             self$config <- config  # Store configuration
                             if(is.null(self$config)) {
                               self$config <- EnceladusConfig$new()  # Create default if none provided
                             }
                             self$validate_parameters()
                           },
                           
                           # Validate input parameters
                           validate_parameters = function() {
                             if(self$params$k <= 0) stop("Permeability must be positive")
                             if(self$params$layer_thickness <= 0) stop("Layer thickness must be positive")
                             if(self$params$porosity <= 0 || self$params$porosity >= 1) {
                               stop("Porosity must be between 0 and 1")
                             }
                           },
                           
                           # Calculate tidal stress variation
                           calculate_tidal_stress = function(t) {
                             # Debug print
                             print(sprintf("Calculating stress for %d time points", length(t)))
                             print(sprintf("Using orbital period: %f", self$config$orbital_period))
                             print(sprintf("Using max stress: %f", self$config$max_tidal_stress))
                             
                             # Calculate phase in orbital cycle
                             phase <- 2 * pi * t / self$config$orbital_period
                             # Calculate stress - ensure vectorized operation
                             stress <- self$config$max_tidal_stress * sin(phase) #any need to incorporate real/imaginary parts of tidal stress?
                             
                             # Debug print
                             print(sprintf("Calculated %d stress values", length(stress)))
                             
                             return(stress)
                           },
                           
                           
                           # Run full simulation
                           run_simulation = function(timesteps = 1000) {
                             # Debug print
                             print("Starting simulation")
                             print(sprintf("Timesteps: %d", timesteps))
                             
                             # Create time series for one orbital period
                             t <- seq(0, self$config$orbital_period, length.out = timesteps)
                             print(sprintf("Created time vector of length %d", length(t)))
                             
                             # Calculate stress and pressure variations
                             stress <- self$calculate_tidal_stress(t)
                             print(sprintf("Got stress vector of length %d", length(stress)))
                             
                             pressure_gradient <- stress / self$params$layer_thickness
                             print(sprintf("Pressure gradient vector length: %d", length(pressure_gradient)))
                             
                             # Calculate flow parameters
                             flow_area <- 4 * pi * self$config$core_radius^2  # Surface area of core
                             print(sprintf("Flow area: %.2e m²", flow_area))
                             
                             # Create results tibble with explicit lengths
                             results <- tibble(
                               time_s = t,
                               time_hr = t/3600,
                               stress = stress,
                               pressure_gradient = pressure_gradient
                             )
                             print("Created initial results tibble")
                             print(sprintf("Results tibble rows: %d", nrow(results)))
                             
                             # Calculate flow rates
                             flow_velocity <- -(self$params$k/CONSTANTS$FLUID$viscosity) * results$pressure_gradient
                             print(sprintf("Flow velocity vector length: %d", length(flow_velocity)))
                             
                             results <- results %>%
                               mutate(
                                 flow_velocity = flow_velocity,
                                 flow_rate_m3s = abs(flow_velocity) * flow_area * self$params$porosity,
                                 flow_rate_Ls = flow_rate_m3s * 1000
                               )
                             print("Added flow rates")
                             
                             # Calculate fluid volumes
                             dt <- diff(t)[1]  # Time step
                             print(sprintf("Time step: %.2e s", dt))
                             
                             results <- results %>%
                               mutate(
                                 fluid_volume = flow_rate_m3s * dt,
                                 fluid_volume_L = fluid_volume * 1000,
                                 cumulative_fluid_m3 = cumsum(fluid_volume)
                               )
                             print("Added volumes")
                             
                             # Calculate rock reactions
                             vertical_area <- GeometricUtils$calculate_vertical_area(self$config$core_radius)
                             print(sprintf("Vertical area: %.2e m²", vertical_area))
                             
                             permeable_volume <- vertical_area * self$params$layer_thickness * (1 - self$params$porosity)
                             print(sprintf("Permeable volume: %.2e m³", permeable_volume))
                             
                             results <- results %>%
                               mutate(
                                 rock_mass_reacted = fluid_volume * self$config$core_density * 
                                   (fluid_volume / permeable_volume),
                                 cumulative_rock = cumsum(rock_mass_reacted)
                               )
                             print("Added rock reactions")
                             
                             self$results <- results
                             print("Stored results")
                             print(sprintf("Final results tibble rows: %d", nrow(self$results)))
                             
                             return(results)
                           },
                           
                           # Generate summary metrics
                           get_summary_metrics = function() {
                             if (is.null(self$results)) stop("Must run simulation first")
                             
                             circ_time <- FlowUtils$calculate_circulation_time(
                               flow_rate = mean(self$results$flow_rate_m3s),
                               ocean_mass = self$config$get_derived_properties()$ocean_mass
                             )
                             orbital_period_years <- self$config$orbital_period / (365.25 * 24 * 3600)
                             n_orb <- circ_time / orbital_period_years  # orbital cycles per ocean circulation

                             r_core       <- self$config$core_radius
                             layer_thick  <- self$params$layer_thickness
                             layer_vol    <- (4/3) * pi * (r_core^3 - (r_core - layer_thick)^3)
                             porewater_mass_est <- 1000 * self$params$porosity * layer_vol

                             fluid_mass        <- sum(self$results$fluid_volume) * CONSTANTS$FLUID$density
                             eff_fluid_mass    <- min(fluid_mass, porewater_mass_est)
                             corr_circ_time    <- circ_time * (fluid_mass / eff_fluid_mass)
                             corr_n_orb        <- corr_circ_time / orbital_period_years

                             metrics <- list(
                               max_flow_rate_Ls = max(self$results$flow_rate_Ls),
                               mean_flow_rate_Ls = mean(self$results$flow_rate_Ls),
                               total_fluid_volume_m3 = sum(self$results$fluid_volume),
                               exchanging_fluid_mass = fluid_mass,
                               total_rock_reacted = sum(self$results$rock_mass_reacted),
                               circulation_time_years = circ_time,
                               complete_circulation_time = (n_orb * log(n_orb) + 0.577 * n_orb) * orbital_period_years,
                               corrected_circulation_time = corr_circ_time,
                               corrected_complete_circulation_time = (corr_n_orb * log(corr_n_orb) + 0.577 * corr_n_orb) * orbital_period_years
                             )
                             
                             print("Debug metrics:")
                             print(sprintf("Flow rate: %.2e m³/s", mean(self$results$flow_rate_m3s)))
                             print(sprintf("Ocean mass: %.2e kg", self$config$get_derived_properties()$ocean_mass))
                             
                             return(metrics)
                           },
                           
                           # Generate plots
                           generate_plots = function() {
                             if (is.null(self$results)) stop("Must run simulation first")
                             
                             p1 <- ggplot(self$results, aes(x = time_hr, y = stress)) +
                               geom_line() +
                               labs(title = "Tidal Stress", x = "Time (hours)", y = "Stress (Pa)")
                             
                             p2 <- ggplot(self$results, aes(x = time_hr, y = flow_rate_Ls)) +
                               geom_line() +
                               labs(title = "Flow Rate", x = "Time (hours)", y = "Flow Rate (L/s)")
                             
                             p3 <- ggplot(self$results, aes(x = time_hr, y = cumulative_rock)) +
                               geom_line() +
                               labs(title = "Cumulative Rock Reacted", 
                                    x = "Time (hours)", y = "Rock Mass (kg)")
                             
                             return(p1 + p2 + p3)
                           }
                         )
)

#===============================================================================
# PHREEQC Integration Module
#===============================================================================
PhreeqcIntegrator <- R6::R6Class("PhreeqcIntegrator",
                                 public = list(
                                   simulator = NULL,
                                   phreeqc_params = NULL,
                                   organic_wt_percent = 0,
                                   iom_pools = NULL,
                                   profile_file = "EnceladusProfile_Seawater_10.0ppt_Tb272.4578K_PorousRock.txt",
                                   
                                   initialize = function(simulator) {
                                     if (!inherits(simulator, "HydrothermalSimulator")) {
                                       stop("Must provide HydrothermalSimulator instance")
                                     }
                                     self$simulator <- simulator
                                     self$calculate_phreeqc_parameters()
                                   },
                                   
                                   # Add organic_wt_percent parameter to the PhreeqcIntegrator
                                   calculate_phreeqc_parameters = function(organic_wt_percent = 0) {
                                     if (is.null(self$simulator$results)) {
                                       stop("Must run simulation before calculating PHREEQC parameters")
                                     }
                                     
                                     # Get simulator metrics and properties
                                     metrics <- self$simulator$get_summary_metrics()
                                     props <- self$simulator$config$get_derived_properties()
                                     
                                     # In your calculate_phreeqc_parameters method, replace the grouping section with:
                                     
                                     # Calculate layer fraction as volume ratio (layer shell / full core)
                                     # This ensures mineral moles and porewater_mass_norm share the same
                                     # reference mass (core_mass), giving a physically correct W/R ratio.
                                     r_core_gf     <- self$simulator$config$core_radius
                                     layer_thick_gf <- self$simulator$params$layer_thickness
                                     layer_vol      <- (4/3) * pi * (r_core_gf^3 - (r_core_gf - layer_thick_gf)^3)
                                     core_vol       <- (4/3) * pi * r_core_gf^3
                                     layer_fraction <- layer_vol / core_vol

                                     # Use empirical bounds from your analysis (also as volume fractions)
                                     min_circulation_time <- 0.1        # Fastest: k=1e-08, 10m layer
                                     max_circulation_time <- 2.45e9     # Slowest: k=1e-14, full core
                                     min_layer_vol      <- (4/3) * pi * (r_core_gf^3 - (r_core_gf - 10)^3)
                                     min_layer_fraction <- min_layer_vol / core_vol   # 10m layer volume fraction
                                     max_layer_fraction <- 1.0                        # Full core
                                     
                                     # Adaptive grouping targeting 1-1000 year time steps
                                     layer_norm <- pmax(0, pmin(1, (layer_fraction - min_layer_fraction) / (max_layer_fraction - min_layer_fraction)))
                                     time_norm <- pmax(0, pmin(1, (metrics$circulation_time_years - min_circulation_time) / (max_circulation_time - min_circulation_time)))
                                     
                                     # Target time step: 1 year (small/fast) to 1000 years (large/slow)
                                     target_time_step_years <- 1 + (layer_norm * time_norm * 999)
                                     
                                     # Calculate required grouping factor
                                     grouping_factor <- max(1, floor(metrics$circulation_time_years / target_time_step_years))
                                     
                                     # Cap grouping for extreme cases
                                     grouping_factor <- min(grouping_factor, 100000)
                                     
                                     # Debug information
                                     print("Debug adaptive grouping:")
                                     print(sprintf("Layer fraction: %.5f (normalized: %.3f)", layer_fraction, layer_norm))
                                     print(sprintf("Circulation time: %.1f years (normalized: %.6f)", metrics$circulation_time_years, time_norm))
                                     print(sprintf("Target time step: %.1f years", target_time_step_years))
                                     print(sprintf("Final grouping factor: %d", grouping_factor))
                                     print(sprintf("Actual time step will be: %.1f years", metrics$circulation_time_years / grouping_factor))
                                     
                                     # Example outcomes for your test cases:
                                     # k=1e-08, 10m:     layer_norm≈0, time_norm≈0     → 1-year steps
                                     # k=1e-14, full:    layer_norm≈1, time_norm≈1     → 1000-year steps  
                                     # k=1e-12, 600m:    layer_norm≈0.003, time_norm≈mid → intermediate steps
                                     
                                     print("Debug PhreeqcIntegrator calculations:")
                                     print(sprintf("Rock mass reacted: %.2e kg", metrics$total_rock_reacted))
                                     print(sprintf("Core mass: %.2e kg", props$core_mass))
                                     print(sprintf("Ocean mass: %.2e kg", props$ocean_mass))
                                     print(sprintf("Layer fraction: %.3f", layer_fraction))
                                     print(sprintf("Circulation time: %.2f years", metrics$circulation_time_years))
                                     
                                     # Solid rock mass of the permeable layer (same density assumption as core)
                                     layer_mass <- layer_fraction * props$core_mass

                                     # All PHREEQC quantities normalised to layer_mass so that mineral m0 = 2.8833
                                     # (not 2.8833 * layer_fraction) and water masses are ~0.2 kg for every layer size.
                                     # W/R ratio is unchanged: dividing both numerator and denominator by layer_mass cancels.
                                     rock_per_cycle <- metrics$total_rock_reacted / layer_mass

                                     # rock_per_year scaled so total_reaction_time = layer_fraction * circ_time (unchanged)
                                     if(!is.null(metrics$circulation_time_years) && length(metrics$circulation_time_years) > 0) {
                                       rock_per_year <- props$core_mass / (layer_mass * metrics$circulation_time_years)
                                     } else {
                                       rock_per_year <- NA
                                       warning("Circulation time not calculated, cannot determine rock_per_year")
                                     }

                                     water_mass <- props$ocean_mass / layer_mass

                                     # Store parameters
                                     self$organic_wt_percent <- organic_wt_percent
                                     # Porewater mass: physical mass of water filling the permeable layer pore space
                                     r_core <- self$simulator$config$core_radius
                                     layer_thick <- self$simulator$params$layer_thickness
                                     layer_volume <- (4/3) * pi * (r_core^3 - (r_core - layer_thick)^3)  # m^3
                                     porewater_mass <- 1000 * self$simulator$params$porosity * layer_volume  # kg (rho_water * V_pore)

                                     self$phreeqc_params <- list(
                                       grouping_factor = grouping_factor,
                                       total_rock = 1.0,
                                       layer_mass = layer_mass,
                                       rock_per_cycle = rock_per_cycle,
                                       rock_per_year = rock_per_year,
                                       water_mass = water_mass,
                                       porewater_mass = porewater_mass,
                                       ocean_mass = props$ocean_mass,
                                       core_mass = props$core_mass
                                     )

                                     print("Calculated PHREEQC parameters:")
                                     print(sprintf("rock_per_cycle: %.2e", rock_per_cycle))
                                     print(sprintf("rock_per_year: %.2e", rock_per_year))
                                     print(sprintf("water_mass: %.2f", water_mass))
                                     
                                     return(self$phreeqc_params)
                                   },

                                   # Secondary minerals only stable/relevant above a temperature threshold —
                                   # returns the subset of names whose threshold is met at temp_C, for use in
                                   # a given block's EQUILIBRIUM_PHASES list (each block gated by its own
                                   # solution's temperature: porewater_temp_C or ocean_temp_C).
                                   temp_gated_minerals = function(temp_C) {
                                     thresholds_C <- c(Quartz = 150, Antigorite = 200, Foshagite = 200,
                                                        Hatrurite = 200, Epidote = 300, `Epidote-ord` = 300)
                                     names(thresholds_C)[temp_C >= thresholds_C]
                                   },

                                   # CH4 vs Mtg (decoupled methane) mode for a generated .pqi file.
                                   # Default (CH4_redox_override = NULL): decoupled below 150 C, coupled at/above.
                                   # override = "coupled"/"decoupled" forces the choice regardless of temperature.
                                   # The two modes require different databases (CH4 master species present vs.
                                   # commented out in favor of the independent Mtg species) — see CHANGELOG.md.
                                   resolve_ch4_mode = function(porewater_temp_C, CH4_redox_override = NULL) {
                                     if (!is.null(CH4_redox_override) && !CH4_redox_override %in% c("coupled", "decoupled")) {
                                       stop('CH4_redox_override must be NULL, "coupled", or "decoupled"')
                                     }
                                     decoupled <- if (!is.null(CH4_redox_override)) {
                                       CH4_redox_override == "decoupled"
                                     } else {
                                       porewater_temp_C < 150
                                     }
                                     list(
                                       decoupled     = decoupled,
                                       ch4_aq        = if (decoupled) "Mtg"    else "CH4",
                                       ch4_gas       = if (decoupled) "Mtg(g)" else "CH4(g)",
                                       database      = if (decoupled) "coreclath_CH4uncoupled.dat" else "coreclath_CH4coupled.dat",
                                       hydrate_phase = if (decoupled) "Mtg_hydrate" else "CH4_hydrate"
                                     )
                                   },

                                   # Rewrites literal CH4/CH4(g) tokens to Mtg/Mtg(g) in already-rendered
                                   # PHREEQC input text, for decoupled-mode files. Applied once to the final
                                   # assembled text so every block (static text and vector-built lists alike)
                                   # is covered without threading a placeholder through each one individually.
                                   apply_ch4_mode = function(text, ch4_mode) {
                                     if (!ch4_mode$decoupled) return(text)
                                     text <- gsub("CH4\\(g\\)", "Mtg(g)", text)
                                     text <- gsub("(?<![A-Za-z0-9_])CH4(?![A-Za-z0-9_\\(])", "Mtg", text, perl = TRUE)
                                     text
                                   },

                                   # Clathrate stability at the given porewater conditions, for CH4, CO2, and
                                   # H2S hydrates. For each gas, computes the dissociation pressure (the
                                   # pressure above which the clathrate is thermodynamically stable at
                                   # porewater_temp_C) and compares it against the actual porewater pressure.
                                   # Dissociation-pressure formulas are placeholders — TODO fill in.
                                   resolve_clathrate_stability = function(porewater_temp_C, porewater_pressure_atm) {
                                     porewater_temp_K <- porewater_temp_C+273.15
                                     # ---- CH4 clathrate dissociation pressure (atm) ----
                                     dissoc_pressure_CH4_kPa <- 4.654e-10*exp(.1074*porewater_temp_K)
                                     dissoc_pressure_CH4_atm <- dissoc_pressure_CH4_kPa/101.3
                                     
                                     # ---- CO2 clathrate dissociation pressure (atm) ----
                                     temp_threshold_CO2 <- 283
                                     if (porewater_temp_K <= temp_threshold_CO2) {
                                        dissoc_pressure_CO2_kPa <- 1.149e-12*exp(0.1269*porewater_temp_K)
                                     } else {
                                        dissoc_pressure_CO2_kPa <- 8518.6*porewater_temp_K - 2.407e6
                                     }
                                     dissoc_pressure_CO2_atm <- dissoc_pressure_CO2_kPa/101.3

                                     # ---- H2S clathrate dissociation pressure (atm) ----
                                     temp_threshold_H2S <- 301.6
                                     if (porewater_temp_K <= temp_threshold_H2S) {
                                        dissoc_pressure_H2S_kPa <- 2.558e-11*exp(.1062*porewater_temp_K)
                                     } else {
                                        dissoc_pressure_H2S_kPa <- 1.172e4*porewater_temp_K - 3.536e6
                                     }
                                     dissoc_pressure_H2S_atm <- dissoc_pressure_H2S_kPa/101.3

                                     list(
                                       CH4 = list(
                                         dissociation_pressure_atm = dissoc_pressure_CH4_atm,
                                         P_stable = porewater_pressure_atm >= dissoc_pressure_CH4_atm
                                       ),
                                       CO2 = list(
                                         dissociation_pressure_atm = dissoc_pressure_CO2_atm,
                                         P_stable = porewater_pressure_atm >= dissoc_pressure_CO2_atm,
                                         threshold_temperature_C = temp_threshold_CO2 - 273.15,
                                         T_stable = porewater_temp_K <= temp_threshold_CO2
                                       ),
                                       H2S = list(
                                         dissociation_pressure_atm = dissoc_pressure_H2S_atm,
                                         P_stable = porewater_pressure_atm >= dissoc_pressure_H2S_atm,
                                         threshold_temperature_C = temp_threshold_H2S - 273.15,
                                         T_stable = porewater_temp_K <= temp_threshold_H2S
                                       )
                                     )
                                   },

                                   # Builds the PQI header comment block reporting clathrate stability,
                                   # from resolve_clathrate_stability()'s per-species checks. CH4 has no
                                   # T_stable check (no threshold temperature defined for it), so it gets
                                   # one line instead of two; CO2 and H2S each get a pressure line and a
                                   # temperature line. Ends with a summary line of which species are
                                   # considered stable (P_stable, and T_stable where applicable) vs not.
                                   clathrate_stability_note = function(porewater_temp_C, porewater_pressure_atm, location = "Porewater") {
                                     stability <- self$resolve_clathrate_stability(porewater_temp_C, porewater_pressure_atm)

                                     lines <- unlist(lapply(names(stability), function(sp) {
                                       s <- stability[[sp]]

                                       p_op   <- if (isTRUE(s$P_stable)) ">=" else "<"
                                       p_verb <- if (isTRUE(s$P_stable)) "may be stable" else "are not stable"
                                       p_line <- sprintf(
                                         "# %s pressure %.3g atm %s dissociation pressure %.3g atm, so %s clathrates %s.",
                                         location, porewater_pressure_atm, p_op, s$dissociation_pressure_atm, sp, p_verb
                                       )

                                       if (is.null(s$T_stable) || !isTRUE(s$P_stable)) return(p_line)

                                       t_op   <- if (isTRUE(s$T_stable)) "<" else ">"
                                       t_line <- if (isTRUE(s$T_stable)) {
                                         sprintf(
                                           "# %s temperature %.3g C %s threshold temperature %.3g C, so %s clathrates remain stable.",
                                           location, porewater_temp_C, t_op, s$threshold_temperature_C, sp
                                         )
                                       } else {
                                         sprintf(
                                           "# %s temperature %.3g C %s threshold temperature %.3g C, so small temperature fluctuations will make %s clathrates unstable.",
                                           location, porewater_temp_C, t_op, s$threshold_temperature_C, sp
                                         )
                                       }

                                       c(p_line, t_line)
                                     }))

                                     overall_stable <- vapply(names(stability), function(sp) {
                                       s <- stability[[sp]]
                                       isTRUE(s$P_stable) && (is.null(s$T_stable) || isTRUE(s$T_stable))
                                     }, logical(1))

                                     summary_line <- sprintf(
                                       "# %s: clathrates considered stable: %s. Clathrates considered unstable: %s.",
                                       location,
                                       if (any(overall_stable)) paste(names(stability)[overall_stable], collapse = ", ") else "none",
                                       if (any(!overall_stable)) paste(names(stability)[!overall_stable], collapse = ", ") else "none"
                                     )

                                     paste(c(lines, summary_line), collapse = "\n")
                                   },

                                   # PHASES names (from coreclath_CH4uncoupled.dat / the coupled counterpart)
                                   # of the clathrates that are stable — P_stable, and T_stable where
                                   # applicable — at the given conditions. `ch4_mode` (from resolve_ch4_mode())
                                   # supplies the correct CH4- vs Mtg-hydrate phase name to match whichever
                                   # database this file recommends. Meant to be added to an EQUILIBRIUM_PHASES
                                   # block and the matching SELECTED_OUTPUT -equilibrium_phases list.
                                   stable_clathrate_phases = function(porewater_temp_C, porewater_pressure_atm, ch4_mode) {
                                     stability <- self$resolve_clathrate_stability(porewater_temp_C, porewater_pressure_atm)

                                     phase_names <- c(
                                       CH4 = ch4_mode$hydrate_phase,
                                       CO2 = "CO2_hydrate",
                                       H2S = "H2S_hydrate"
                                     )

                                     overall_stable <- vapply(names(stability), function(sp) {
                                       s <- stability[[sp]]
                                       isTRUE(s$P_stable) && (is.null(s$T_stable) || isTRUE(s$T_stable))
                                     }, logical(1))

                                     unname(phase_names[names(stability)][overall_stable])
                                   },

                                   generate_phreeqc_input_tidal_kinetic = function(output_dir = "phreeqc_inputs", mode = "ocean_circulations",
                                                                     n_circulations = NULL, n_orbital_cycles = NULL,
                                                                     grouping_factor = NULL, exchange_ratio = NULL,
                                                                     porewater_temp = "ocean", substeps = 1, suffix = "",
                                                                     CH4_redox_override = NULL) {

                                     # Auto-compute grouping_factor from exchange_ratio if provided.
                                     # exchange_ratio = target (fluid / porewater) fraction per model cycle.
                                     # orbits_per_step = exchange_ratio / (fluid_per_orbit / porewater_mass)
                                     if (!is.null(exchange_ratio)) {
                                       if (exchange_ratio <= 0) stop("exchange_ratio must be > 0")
                                       if (is.null(n_orbital_cycles)) stop("exchange_ratio requires n_orbital_cycles")
                                       .metrics   <- self$simulator$get_summary_metrics()
                                       .f_over_pw <- .metrics$exchanging_fluid_mass / self$phreeqc_params$porewater_mass
                                       grouping_factor <- max(1L, as.integer(round(exchange_ratio / .f_over_pw)))
                                       mode <- "grouped"
                                     }


                                     if (!mode %in% c("ocean_circulations", "ungrouped", "grouped")) {
                                       stop("mode must be 'ocean_circulations', 'ungrouped', or 'grouped'")
                                     }
                                     if (mode == "ocean_circulations") {
                                       if (is.null(n_circulations)) stop("ocean_circulations mode requires n_circulations")
                                       n  <- n_circulations
                                       gf <- self$phreeqc_params$grouping_factor  # auto-computed
                                     } else if (mode == "ungrouped") {
                                       if (is.null(n_orbital_cycles)) stop("ungrouped mode requires n_orbital_cycles")
                                       n  <- n_orbital_cycles
                                       gf <- NULL
                                     } else {  # grouped
                                       if (is.null(n_orbital_cycles)) stop("grouped mode requires n_orbital_cycles")
                                       if (is.null(grouping_factor))  stop("grouped mode requires grouping_factor")
                                       n  <- n_orbital_cycles
                                       gf <- grouping_factor
                                     }
                                   
                                     # Create output directory if needed
                                     dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

                                     # Filename encodes mode:
                                     #   circ = ocean_circulations (auto grouping + n ocean circulations)
                                     #   orb  = ungrouped (1 orbital period per step)
                                     #   grp  = grouped (user grouping factor + n orbital cycles)
                                     temp_tag <- if (is.numeric(porewater_temp)) sprintf("T%g", porewater_temp) else sprintf("T%s", porewater_temp)
                                     org_tag  <- sprintf("org%g", self$organic_wt_percent)
                                     sub_tag  <- if (substeps != 1) sprintf("_sub%d", substeps) else ""
                                     if (mode == "ocean_circulations") {
                                       filename <- file.path(output_dir, sprintf(
                                         "enceladus_k%.2e_d%g_%s_%s_t%.2e_circ%d%s%s.pqi",
                                         self$simulator$params$k,
                                         self$simulator$params$layer_thickness,
                                         temp_tag, org_tag,
                                         self$phreeqc_params$grouping_factor,
                                         n, sub_tag, suffix
                                       ))
                                     } else if (mode == "ungrouped") {
                                       filename <- file.path(output_dir, sprintf(
                                         "enceladus_k%.2e_d%g_%s_%s_orb%g%s%s.pqi",
                                         self$simulator$params$k,
                                         self$simulator$params$layer_thickness,
                                         temp_tag, org_tag,
                                         n, sub_tag, suffix
                                       ))
                                     } else {  # grouped
                                       filename <- file.path(output_dir, sprintf(
                                         "enceladus_k%.2e_d%g_%s_%s_grp%g_orb%g%s%s.pqi",
                                         self$simulator$params$k,
                                         self$simulator$params$layer_thickness,
                                         temp_tag, org_tag,
                                         gf,
                                         n, sub_tag, suffix
                                       ))
                                     }

                                     # Generate input file content
                                     content <- self$create_phreeqc_template_tidal_kinetic(filename, n, mode, grouping_factor = gf, exchange_ratio = exchange_ratio, porewater_temp = porewater_temp, substeps = substeps, CH4_redox_override = CH4_redox_override)

                                     # Write to file
                                     writeLines(content, filename)

                                     return(filename)
                                   },
                                   
                                   create_phreeqc_template_tidal_kinetic = function(filename, n, mode = "ocean_circulations", grouping_factor = NULL, exchange_ratio = NULL, porewater_temp = "ocean", substeps = 1, CH4_redox_override = NULL) {
                                     metrics <- self$simulator$get_summary_metrics()

                                     # Calculate adjusted cycles
                                     orbital_period <- self$simulator$config$orbital_period

                                     # Calculate total cycles based on rock_per_year
                                     if(!is.na(self$phreeqc_params$rock_per_year) && self$phreeqc_params$rock_per_year > 0) {
                                       total_reaction_time <- self$phreeqc_params$total_rock / self$phreeqc_params$rock_per_year
                                       total_cycles <- ceiling((total_reaction_time * 365.25 * 24 * 3600) / orbital_period)
                                     } else {
                                       # Fallback if rock_per_year isn't available
                                       total_cycles <- 1000
                                     }

                                     if (mode == "ocean_circulations") {
                                       # Bundle many orbital cycles using the auto-computed grouping factor.
                                       # n = number of full ocean circulations.
                                       adjusted_cycles   <- ceiling(total_cycles / self$phreeqc_params$grouping_factor)
                                       organic_cycles    <- adjusted_cycles
                                       kinetic_cycles    <- ceiling(adjusted_cycles * n)
                                       time_step_seconds <- (metrics$circulation_time_years / adjusted_cycles) * 365.25 * 24 * 3600
                                     } else if (mode == "grouped") {
                                       # User-specified grouping: n = n_orbital_cycles, each step covers grouping_factor orbital periods.
                                       kinetic_cycles    <- ceiling(n / grouping_factor)
                                       adjusted_cycles   <- kinetic_cycles
                                       organic_cycles    <- kinetic_cycles
                                       time_step_seconds <- grouping_factor * orbital_period
                                     } else {
                                       # Ungrouped: each KINETICS step = exactly 1 orbital period.
                                       # n = number of individual orbital cycles.
                                       adjusted_cycles   <- n
                                       organic_cycles    <- n
                                       kinetic_cycles    <- n
                                       time_step_seconds <- orbital_period
                                     }

                                     # Unified time_step_years used in TITLE comment
                                     time_step_years      <- time_step_seconds / (365.25 * 24 * 3600)
                                     total_sim_years      <- kinetic_cycles * time_step_years
                                     mode_description <- if (mode == "ocean_circulations") {
                                       sprintf("ocean_circulations (auto grouping factor: %d)", self$phreeqc_params$grouping_factor)
                                     } else if (mode == "grouped") {
                                       sprintf("grouped (%d orbital cycles per step, %d total orbital cycles)", grouping_factor, n)
                                     } else {
                                       "ungrouped (1 orbital period per step)"
                                     }

                                     # Fluid masses
                                     porewater_mass      <- self$phreeqc_params$porewater_mass
                                     ocean_mass_actual   <- self$phreeqc_params$ocean_mass
                                     core_mass           <- self$phreeqc_params$core_mass
                                     exchanging_fluid_mass <- metrics$exchanging_fluid_mass

                                     # Mixing fractions using m_total = m_pore + m_ocean as reference.
                                     # Each cycle: m_fluid leaves porewater → enters ocean,
                                     #             then m_fluid leaves ocean → returns to porewater.
                                      m_total <- porewater_mass + ocean_mass_actual
                                      f_pore  <- porewater_mass      / m_total  # porewater fraction of total
                                      f_ocean <- ocean_mass_actual   / m_total  # ocean fraction of total
                                      # (f_pore + f_ocean = 1 by definition)

                                      # --- Compounded exchange_ratio over grouped tidal cycles ---------------------
                                      # exchanging_fluid_mass is integrated over ONE orbital period, but each
                                      # unrolled cycle block advances time_step_seconds = G orbital periods.
                                      # Partial pore-volume replacement compounds over the G cycles:
                                      #     x_step = 1 - (1 - x_tide)^G   (bounded above by 1)
                                      # This makes grouped and ungrouped runs of equal physical duration
                                      # exchange the same total fluid (to first order).
                                      G_cycles <- time_step_seconds / orbital_period
                                      x_tide   <- min(exchanging_fluid_mass / porewater_mass, 1)  # per tidal cycle
                                      x_step   <- 1 - (1 - x_tide)^G_cycles                     # per unrolled step
                                      f_x      <- x_step * f_pore   # exchanged mass as fraction of m_total
                                      if (x_step > 0.3)
                                        warning(sprintf(paste0(
                                          "x_step = %.2f: >30%% of pore volume replaced per grouped step; ",
                                          "grouping approximation is degrading — consider a smaller grouping factor."),
                                          x_step))
                                     # PHREEQC-normalized water masses (divided by core_mass so they match mineral mole scaling)
                                     porewater_mass_norm <- porewater_mass    / self$phreeqc_params$layer_mass
                                     ocean_mass_norm     <- self$phreeqc_params$water_mass  # already = ocean_mass / layer_mass
                                     
                                     # Get base filename without extension for output
                                     base_name <- tools::file_path_sans_ext(basename(filename))

                                     # Temperatures and pressures for SOLUTION blocks
                                     ocean_temp_C <- 274.6 - 273.15  # ~1.45 °C (ocean-rock interface)
                                     if (!is.null(self$profile_file) && file.exists(self$profile_file)) {
                                       pore_temp_C          <- mean_layer_temp_C(
                                         self$profile_file,
                                         self$simulator$params$layer_thickness
                                       )
                                       pore_pressure_atm    <- mean_layer_pressure_MPa(
                                         self$profile_file,
                                         self$simulator$params$layer_thickness
                                       ) / 0.101325
                                       ocean_pressure_atm   <- ocean_floor_pressure_MPa(
                                         self$profile_file
                                       ) / 0.101325
                                     } else {
                                       if (!is.null(self$profile_file))
                                         warning("Profile file not found: '", self$profile_file, "' — falling back to defaults.")
                                       pore_temp_C        <- 25
                                       pore_pressure_atm  <- 70
                                       ocean_pressure_atm <- 70
                                     }

                                     # Resolve porewater_temp to a concrete value
                                     porewater_temp_C <- if (is.numeric(porewater_temp)) {
                                       porewater_temp
                                     } else if (porewater_temp == "profile") {
                                       pore_temp_C
                                     } else if (porewater_temp == "ocean") {
                                       ocean_temp_C
                                     } else {
                                       stop("porewater_temp must be a number (°C), \"profile\", or \"ocean\"")
                                     }

                                     # CH4 (redox-coupled) vs Mtg (decoupled) — see resolve_ch4_mode().
                                     ch4_mode <- self$resolve_ch4_mode(porewater_temp_C, CH4_redox_override)
                                     ch4_override_note <- if (!is.null(CH4_redox_override)) sprintf(", override=%s", CH4_redox_override) else ""
                                     ch4_db_note_line <- sprintf(
                                       "# Recommended database: %s (CH4 %s; porewater T = %.1f C%s)",
                                       ch4_mode$database,
                                       if (ch4_mode$decoupled) "decoupled -> aqueous/gas species are Mtg/Mtg(g)" else "coupled -> aqueous/gas species are CH4/CH4(g)",
                                       porewater_temp_C, ch4_override_note
                                     )

                                     # Clathrate stability — separate checks for the ocean and porewater
                                     # solutions, since this function gives them distinct temperature/pressure.
                                     # See resolve_clathrate_stability()/clathrate_stability_note().
                                     clathrate_note_lines <- paste0(
                                       self$clathrate_stability_note(ocean_temp_C, ocean_pressure_atm, location = "Ocean"),
                                       "\n\n",
                                       self$clathrate_stability_note(porewater_temp_C, pore_pressure_atm, location = "Porewater")
                                     )

                                     # Temperature-gated secondary minerals, appended per block below.
                                     # Clathrate phases stable at each reservoir's own conditions (see
                                     # resolve_clathrate_stability()) are appended the same way.
                                     extra_pw_minerals    <- c(self$temp_gated_minerals(porewater_temp_C),
                                                               self$stable_clathrate_phases(porewater_temp_C, pore_pressure_atm, ch4_mode))
                                     extra_ocean_minerals <- c(self$temp_gated_minerals(ocean_temp_C),
                                                               self$stable_clathrate_phases(ocean_temp_C, ocean_pressure_atm, ch4_mode))
                                     extra_pw_lines    <- if (length(extra_pw_minerals) > 0)
                                       paste0("   ", extra_pw_minerals, " 0 0", collapse = "\n") else ""
                                     extra_ocean_lines <- if (length(extra_ocean_minerals) > 0)
                                       paste0("   ", extra_ocean_minerals, " 0 0", collapse = "\n") else ""
                                     extra_eq_so_line  <- paste(union(extra_pw_minerals, extra_ocean_minerals), collapse = " ")

                                     # Calculate scaling factor for minerals (reduced by organic content)
                                     mineral_scale_factor <- (100 - self$organic_wt_percent) / 100
                                     
                                     # Calculate total mineral amounts (moles) based on permeable layer fraction.
                                     # Coefficients: ROCK_MOL_PER_KG (top of file). 9 kinetic primaries
                                     # (incl. Schreibersite, placeholder rate) + 3 equilibrium primaries
                                     # (Fe, Ni, Lawrencite). Magnetite is no longer a primary mineral here —
                                     # it remains only as an ordinary secondary/precipitate phase.
                                     enstatite_moles  <- ROCK_MOL_PER_KG[["enstatite"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     forsterite_moles <- ROCK_MOL_PER_KG[["forsterite"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     troilite_moles   <- ROCK_MOL_PER_KG[["troilite"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     albite_moles     <- ROCK_MOL_PER_KG[["albite"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     diopside_moles   <- ROCK_MOL_PER_KG[["diopside"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     anorthite_moles  <- ROCK_MOL_PER_KG[["anorthite"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     kfeldspar_moles  <- ROCK_MOL_PER_KG[["kfeldspar"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     tephroite_moles  <- ROCK_MOL_PER_KG[["tephroite"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     fe_moles         <- ROCK_MOL_PER_KG[["fe"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     ni_moles         <- ROCK_MOL_PER_KG[["ni"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     lawrencite_moles <- ROCK_MOL_PER_KG[["lawrencite"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     schreibersite_moles <- ROCK_MOL_PER_KG[["schreibersite"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                         # ---- IOM as kinetic phases (Miller et al. 2025 / Burnham-Sweeney 1989) --
                                      # Formula unit from gram composition per kg IOM (728 C, 46 H, 159 O,
                                      # 30 N, 37 S; sums to 1000 g/mol, so mol formula units = kg IOM).
                                      # NOTE: v3 comments quoted N_3.284 (46 g N) — gram values used here;
                                      # confirm against source composition.
                                      iom_rates_block    <- ""
                                      iom_kinetics_block <- ""
                                      if (self$organic_wt_percent > 0) {
                                        iom_mass_kg_scaled <- (self$organic_wt_percent / 100) * self$phreeqc_params$total_rock
                                        iom_total_moles    <- iom_mass_kg_scaled * 1.0   # 1 mol per kg by construction
                                        # ---- Channel definitions: iom_module.R, not hardcoded here ----
                                        # WIRING CHANGE (2026-09-26): this used to be a hand-copied single-Ea
                                        # table (A=1e13, IOM_CO2/CH4 at 214/228 kJ/mol, IOM_N at 215 kJ/mol) that
                                        # had drifted three patch-generations behind iom_module.R without either
                                        # side noticing -- see FROM_LUCAS_20260924_review.md for the full
                                        # comparison. iom_default_config() is now the single source of truth:
                                        # distributed-Ea sub-pools for CO2/CH4/CHn (Vitrimat 2018, A=2e15),
                                        # corrected IOM_N (242.8 kJ/mol), and the redox_mode argument for Option
                                        # 3 (decoupled Mtg/Sg formulas, tested against real PHREEQC output --
                                        # see PHREEQC_TEST_RESULTS_2026-09-20.md). Full derivation and history
                                        # are in iom_default_config()'s own docstring; do not duplicate it here,
                                        # it will only go stale again.
                                        #
                                        # This is unchanged: -formula is WHAT THE SOLID LOSES, not the product
                                        # molecule (PROGRESS.md Part 1 explains why at length -- Miller's own O
                                        # balance doesn't close without water contributing oxygen to the
                                        # measured CO2). Don't "simplify" formulas to product molecules.
                                        #
                                        # This code loop is generic over whatever rows iom_default_config()
                                        # returns (currently 29: 7+13+7 CO2/CH4/CHn sub-pools + IOM_N + IOM_S),
                                        # not hardcoded to 4 channels -- verified no PUNCH line-number collision
                                        # with the mineral PUNCH statements elsewhere in this template (those
                                        # top out at line 110; IOM_ punch lines start at 210).
                                        if (is.null(self$iom_pools)) {
                                          self$iom_pools <- iom_default_config()
                                        }

                                        iom_rates_block <- paste0(apply(self$iom_pools, 1, function(p) {
                                          sprintf(paste0(
                                            "%s\n -start\n",
                                            " 1  REM Arrhenius 1st-order kerogen degradation (single-Ea pool)\n",
                                            " 10 k = 10^(%s) * exp(-%s / (8.314 * TK))\n",
                                            " 20 rate = k * M\n",
                                            " 30 if (M <= 0) then rate = 0\n",
                                            " 40 moles = rate * TIME\n",
                                            " 50 SAVE moles\n -end\n"),
                                            p[["name"]], p[["logA"]], p[["Ea_J"]])
                                        }), collapse = "")

                                        iom_kinetics_block <- paste0(apply(self$iom_pools, 1, function(p) {
                                          m0 <- as.numeric(p[["m0_per_kg"]]) * iom_total_moles
                                          sprintf("   %s\n      -formula %s\n      -m0 %.6e\n      -step_divide 1\n",
                                                  p[["name"]], trimws(p[["formula"]]), m0)
                                        }), collapse = "")

                                        print(sprintf("IOM as kinetics: %d channels, %.2e kg IOM scaled",
                                                      nrow(self$iom_pools), iom_mass_kg_scaled))
                                        print(data.frame(channel = self$iom_pools$name,
                                                         m0 = self$iom_pools$m0_per_kg * iom_total_moles,
                                                         Ea_kJ = self$iom_pools$Ea_J / 1e3))
                                      }
                                      iom_so_line <- if (self$organic_wt_percent > 0 && !is.null(self$iom_pools))
                                        paste(self$iom_pools$name, collapse = " ") else ""

                                     template <- glue::glue('
TITLE Enceladus hydrothermal simulation
{ch4_db_note_line}
{clathrate_note_lines}

KNOBS
   -logfile true

# k={self$simulator$params$k}, thickness={self$simulator$params$layer_thickness} m
# Organic content: {self$organic_wt_percent} wt.%
# Mode: {mode_description}
# exchange_ratio (fluid/porewater per model cycle): {if (!is.null(exchange_ratio)) sprintf("%.4g", exchange_ratio) else "N/A"}
# Simulation parameters:
# Layer thickness: {self$simulator$params$layer_thickness} m ({sprintf("%.1f", self$simulator$params$layer_thickness/self$simulator$config$core_radius*100)}% of core)
# Layer mass: {sprintf("%.2e", self$phreeqc_params$layer_mass)} kg  (total_rock normalised to 1.0)
# Reaction time: {sprintf("%.2e", total_sim_years)} years
# Total cycles: {adjusted_cycles}
# Time step: {sprintf("%.2e", time_step_years)} years per step
# Mean flow rate: {sprintf("%.2e", metrics$mean_flow_rate_Ls)} L/s
# Circulation time: {sprintf("%.2e", metrics$circulation_time_years)} years
# Complete circulation time: {sprintf("%.2e", metrics$complete_circulation_time)} years
# Corrected circulation time (f_fluid capped at f_pore): {sprintf("%.2e", metrics$corrected_circulation_time)} years
# Corrected complete circulation time: {sprintf("%.2e", metrics$corrected_complete_circulation_time)} years
# Flowing fluid mass (per orbital cycle): {sprintf("%.2e", metrics$exchanging_fluid_mass)} kg
# Porewater mass: {sprintf("%.2e", self$phreeqc_params$porewater_mass)} kg
# Ocean mass: {sprintf("%.2e", self$phreeqc_params$ocean_mass)} kg

SOLUTION 1 Enceladus porewater
    temp      {round(porewater_temp_C, 2)}
    pH        7
    pe        4
    redox     pe
    units     mol/kgw
    density   1
    water     {porewater_mass_norm}
    pressure  {round(ocean_pressure_atm, 1)}#{round(pore_pressure_atm, 1)}

SOLUTION 2 Enceladus ocean
    temp      {round(ocean_temp_C, 2)}
    pH        7
    pe        4
    redox     pe
    units     mol/kgw
    density   1
    water     {ocean_mass_norm}
    pressure  {round(ocean_pressure_atm, 1)}

REACTION_TEMPERATURE 1
    {round(porewater_temp_C, 2)}

REACTION_TEMPERATURE 2
    {round(ocean_temp_C, 2)}

# Define rate equations for each mineral
RATES
Forsterite
	-start
	1   REM Ref PK04
	10  kacid = 10^(-6.85) * exp(-67.2e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^0.47
	20  kneut = 10^(-10.64) * exp(-79.0e3/8.314 * (1/TK-1/298.15))
	21  SSA = 0.1
	22  mw = 140.6715
	40  k = (kacid + kneut) * SSA * mw * M
	50  IF SR("Forsterite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Forsterite"))
	60  moles = rate * TIME
	70  SAVE moles
	-end
	
	Troilite # PK04 has no separate troilite entry; hexagonal pyrrhotite parameters used as proxy (same FeS chemistry)
        -start
        1   REM Acid mechanism only (PK04)
        2   REM Hexagonal pyrrhotite formulation (reaction orders H+ = -0.090, Fe3+ = 0.356) used for troilite
        3   REM Rate depends on H+ and Fe3+ activities
        25  kacid = 10^(-6.79) * exp(-63.0e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^-0.090 * ACT("Fe+3")^0.356
        26  SSA = 5
        27  mw = 87.913
        30  k = kacid * SSA * mw * M
        40  IF SR("Troilite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Troilite"))
        50  moles = rate * TIME
        55  IF moles < 0 THEN moles = 0
        60  SAVE moles
        -end

	Diopside
	-start
	1   REM Ref PK04
	10  kacid = 10^(-6.36) * exp(-96.1e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^0.71
	20  kneut = 10^(-11.11) * exp(-40.6e3/8.314 * (1/TK-1/298.15))
	21  SSA = 0.1
	22  mw = 216.55
	40  k = (kacid + kneut) * SSA * mw * M
	50  IF SR("Diopside") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Diopside"))
	60  moles = rate * TIME
	70  SAVE moles
	-end

	K-Feldspar
	-start
	1   REM Ref PK04
	10  kacid = 10^(-10.06) * exp(-51.7e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^0.5
	20  kneut = 10^(-12.41) * exp(-38.0e3/8.314 * (1/TK-1/298.15))
	30  kbase = 10^(-21.20) * exp(-94.1e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^-0.823
	31  SSA = 5
	32  mw = 278.33
	40  k = (kacid + kneut + kbase) * SSA * mw * M
	50  IF SR("K-Feldspar") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("K-Feldspar"))
	60  moles = rate * TIME
	70  SAVE moles
	-end

	Tephroite # TODO: rate constants not yet supplied (Mn-olivine) — k held at 0 (no dissolution) until provided
	-start
	1   REM TODO awaiting kacid/kneut from user
	10  kacid = 0
	20  kneut = 0
	21  SSA = 0.1
	22  mw = 201.96
	40  k = (kacid + kneut) * SSA * mw * M
	50  IF SR("Tephroite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Tephroite"))
	60  moles = rate * TIME
	70  SAVE moles
	-end

{SCHREIBERSITE_RATE}

	Enstatite
	-start
	1   REM Ref PK04
	10  kacid = 10^(-9.02) * exp(-80.0e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^0.6
	20  kneut = 10^(-12.72) * exp(-80.0e3/8.314 * (1/TK-1/298.15))
	21  SSA = 0.1
	22  mw = 100.3725
	40  k = (kacid + kneut) * SSA * mw * M
	50  IF SR("Enstatite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Enstatite"))
	60  moles = rate * TIME
	70  SAVE moles
	-end

    Anorthite # Ca-endmember of plagioclase
        -start
        1   REM Ref PK04
        10  kacid = 10^(-3.50) * exp(-16.6e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^1.411
        20  kneut = 10^(-9.12) * exp(-17.8e3/8.314 * (1/TK-1/298.15))
        21  SSA = 5
        22  mw = 278.164
        40  k = (kacid + kneut) * SSA * mw * M
        50  IF SR("Anorthite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Anorthite"))
        60  moles = rate * TIME
        70  SAVE moles
        -end
        
    Albite # Na-endmember of plagioclase
        -start
        1   REM 3 mechanisms: acid, neutral, base (PK04)
        2   REM Chemical affinity parameters p and q for albite are 0.760 and 90.0 respectively
        3   REM (Alekseyev et al., 1997), but their use in modeling should be limited to conditions 
        4   REM near the experimental conditions under which they were obtained, 300 Â°C and pH = 9.
        10  kacid = 10^(-10.16) * exp(-65.0e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^0.457
        20  kneut = 10^(-12.56) * exp(-69.8e3/8.314 * (1/TK-1/298.15))
        30  kbase = 10^(-15.60) * exp(-71.0e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^-0.572
        31  SSA = 5
        32  mw = 262.1798
        40  k = (kacid + kneut + kbase) * SSA * mw * M
        50  IF SR("Albite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Albite"))
        60  moles = rate * TIME
        70  SAVE moles
        -end

{iom_rates_block}

# Define kinetic mineral phases with initial amounts
KINETICS 1
   Enstatite
      -m0 {enstatite_moles}
      -step_divide 1
   Forsterite
      -m0 {forsterite_moles}
      -step_divide 1
   Troilite
      -m0 {troilite_moles}
      -step_divide 1
   Albite
      -m0 {albite_moles}
      -step_divide 1
   Diopside
      -m0 {diopside_moles}
      -step_divide 1
   Anorthite
      -m0 {anorthite_moles}
      -step_divide 1
   K-Feldspar
      -m0 {kfeldspar_moles}
      -step_divide 1
   Tephroite
      -m0 {tephroite_moles}
      -step_divide 1
{schreibersite_kinetics_entry(schreibersite_moles)}
   {iom_kinetics_block}
   -steps {time_step_seconds} in {substeps} steps
 
   -cvode true
   -bad_step_max 5000

    
EQUILIBRIUM_PHASES 1
# --- Primary minerals (equilibrium; dissolve_only so they can only deplete) ---
   Fe                   0  {fe_moles}          dissolve_only
   Ni                   0  {ni_moles}          dissolve_only
   Lawrencite           0  {lawrencite_moles}  dissolve_only
#   Alabandite 0 0
#   Alanine 0 0
#   Alum-K 0 0
#   Alunite 0 0
   Analcime 0 0
   Anhydrite 0 0
   Aragonite 0 0
#   Arcanite 0 0
#   Artinite 0 0
   Bassanite 0 0
   Beidellite-Ca 0 0
   Beidellite-Fe 0 0
#   Beidellite-K 0 0
   Beidellite-Mg 0 0
   Beidellite-Na 0 0
   Boehmite 0 0
   Brucite 0 0
   C2H4(g) 0 0
   C2H6(g) 0 0
   C3H8(g) 0 0
   Calcite 0 0
#   Celadonite 0 0
   CH4(g) 0 0
   Chalcedony 0 0
   Chamosite 0 0
   Chrysotile 0 0
   Citric_Acid 0 0
   Clinochlore-14A 0 0
   Clinochlore-7A 0 0
   Clinoptilolite-Ca 0 0
#   Clinoptilolite-K 0 0
   Clinoptilolite-Na 0 0
   Clinozoisite 0 0
   CO(g) 0 0
   CO2(g) 0 0
   Cronstedtite-7A 0 0
   Daphnite-14A 0 0
   Daphnite-7A 0 0
   Dawsonite 0 0
   Diaspore 0 0
   Dolomite 0 0
   Ettringite 0 0
   Fe(OH)2 0 0
   Fe(OH)3 0 0
   Fe2(SO4)3 0 0
   FeSO4 0 0
   Gibbsite 0 0
   Glycine 0 0
   Goethite 0 0
   Greenalite 0 0
   Gypsum 0 0
   Gyrolite 0 0
   H2(g) 0 0
#   H2O(g) 0 0
   H2S(g) 0 0
#   Halite 0 0
#   Hausmannite 0 0
   Hematite 0 0
   Huntite 0 0
   Hydromagnesite 0 0
#   Hydroxyapatite 0 0
#   Ice 0 0
#   Jarosite 0 0
#   KAl(SO4)2 0 0
   Kaolinite 0 0
#   KerogenC128 0 0
#   KerogenC292 0 0
#   KerogenC406 0 0
#   KerogenC415 0 0
#   KerogenC515 0 0
   Lizardite 0 0
   Melanterite 0 0
#   MgOHCl 0 0
   MgSO4 0 0
   Minnesotaite 0 0
   Mirabilite 0 0
#   Mn(OH)2(am) 0 0
#   MnCl2:2H2O 0 0
#   MnCl2:4H2O 0 0
#   MnCl2:H2O 0 0
#   MnSO4 0 0
#   Monohydrocalcite 0 0
   Montmor-Ca 0 0
#   Montmor-K 0 0
   Montmor-Mg 0 0
   Montmor-Na 0 0
   N2(g) 0 0
#   Na2CO3 0 0
#   Na2CO3:7H2O 0 0
#   Nahcolite 0 0
#   Natron 0 0
#   Nesquehonite 0 0
   NH3(g) 0 0
#   NH4Cl 0 0
#   NH4HCO3 0 0
#   Niter 0 0
#   NO(g) 0 0
#   NO2(g) 0 0
   Nontronite-Ca 0 0
#   Nontronite-K 0 0
   Nontronite-Mg 0 0
   Nontronite-Na 0 0
#   O2(g) 0 0
   Portlandite 0 0
#   Pyridine 0 0
   Pyrite 0 0
#   Rhodochrosite 0 0
   Saponite-Fe-Ca 0 0
   Saponite-Fe-Fe 0 0
#   Saponite-Fe-K 0 0
   Saponite-Fe-Mg 0 0
   Saponite-Fe-Na 0 0
   Saponite-Mg-Ca 0 0
   Saponite-Mg-Fe 0 0
#   Saponite-Mg-K 0 0
   Saponite-Mg-Mg 0 0
   Saponite-Mg-Na 0 0
   Sepiolite 0 0
#   Siderite 0 0
   SiO2(am) 0 0
#   Smectite-high-Fe-Mg 0 0
#   Smectite-low-Fe-Mg 0 0
   SO2(g) 0 0
#   Stilbite 0 0
#   Strengite 0 0
#   Sylvite 0 0
   Thenardite 0 0
#   Thermonatrite 0 0
   Tobermorite-11A 0 0
{extra_pw_lines}

# Seafloor equilibrium phases — same candidate list as porewater (EQUILIBRIUM_PHASES 1).
# Applied to the combined ocean+fluid (SOLUTION 4) in simulation 3 each cycle.
EQUILIBRIUM_PHASES 2
#   Alabandite 0 0
#   Alanine 0 0
#   Alum-K 0 0
#   Alunite 0 0
   Analcime 0 0
   Anhydrite 0 0
   Aragonite 0 0
#   Arcanite 0 0
#   Artinite 0 0
   Bassanite 0 0
   Beidellite-Ca 0 0
   Beidellite-Fe 0 0
#   Beidellite-K 0 0
   Beidellite-Mg 0 0
   Beidellite-Na 0 0
   Boehmite 0 0
   Brucite 0 0
   C2H4(g) 0 0
   C2H6(g) 0 0
   C3H8(g) 0 0
   Calcite 0 0
#   Celadonite 0 0
   CH4(g) 0 0
   Chalcedony 0 0
   Chamosite 0 0
   Chrysotile 0 0
   Citric_Acid 0 0
   Clinochlore-14A 0 0
   Clinochlore-7A 0 0
   Clinoptilolite-Ca 0 0
#   Clinoptilolite-K 0 0
   Clinoptilolite-Na 0 0
   Clinozoisite 0 0
   CO(g) 0 0
   CO2(g) 0 0
   Cronstedtite-7A 0 0
   Daphnite-14A 0 0
   Daphnite-7A 0 0
   Dawsonite 0 0
   Diaspore 0 0
   Dolomite 0 0
   Ettringite 0 0
   Fe(OH)2 0 0
   Fe(OH)3 0 0
   Fe2(SO4)3 0 0
   FeSO4 0 0
   Gibbsite 0 0
   Glycine 0 0
   Goethite 0 0
   Greenalite 0 0
   Gypsum 0 0
   Gyrolite 0 0
   H2(g) 0 0
#   H2O(g) 0 0
   H2S(g) 0 0
#   Halite 0 0
#   Hausmannite 0 0
   Hematite 0 0
   Huntite 0 0
   Hydromagnesite 0 0
#   Hydroxyapatite 0 0
#   Ice 0 0
#   Jarosite 0 0
#   KAl(SO4)2 0 0
   Kaolinite 0 0
#   KerogenC128 0 0
#   KerogenC292 0 0
#   KerogenC406 0 0
#   KerogenC415 0 0
#   KerogenC515 0 0
   Lizardite 0 0
   Magnetite 0 0
   Melanterite 0 0
#   MgOHCl 0 0
   MgSO4 0 0
   Minnesotaite 0 0
   Mirabilite 0 0
#   Mn(OH)2(am) 0 0
#   MnCl2:2H2O 0 0
#   MnCl2:4H2O 0 0
#   MnCl2:H2O 0 0
#   MnSO4 0 0
#   Monohydrocalcite 0 0
   Montmor-Ca 0 0
#   Montmor-K 0 0
   Montmor-Mg 0 0
   Montmor-Na 0 0
   N2(g) 0 0
#   Na2CO3 0 0
#   Na2CO3:7H2O 0 0
#   Nahcolite 0 0
#   Natron 0 0
#   Nesquehonite 0 0
   NH3(g) 0 0
#   NH4Cl 0 0
#   NH4HCO3 0 0
#   Niter 0 0
#   NO(g) 0 0
#   NO2(g) 0 0
   Nontronite-Ca 0 0
#   Nontronite-K 0 0
   Nontronite-Mg 0 0
   Nontronite-Na 0 0
#   O2(g) 0 0
   Portlandite 0 0
#   Pyridine 0 0
   Pyrite 0 0
#   Rhodochrosite 0 0
   Saponite-Fe-Ca 0 0
   Saponite-Fe-Fe 0 0
#   Saponite-Fe-K 0 0
   Saponite-Fe-Mg 0 0
   Saponite-Fe-Na 0 0
   Saponite-Mg-Ca 0 0
   Saponite-Mg-Fe 0 0
#   Saponite-Mg-K 0 0
   Saponite-Mg-Mg 0 0
   Saponite-Mg-Na 0 0
   Sepiolite 0 0
#   Siderite 0 0
   SiO2(am) 0 0
#   Smectite-high-Fe-Mg 0 0
#   Smectite-low-Fe-Mg 0 0
   SO2(g) 0 0
#   Stilbite 0 0
#   Strengite 0 0
#   Sylvite 0 0
   Thenardite 0 0
#   Thermonatrite 0 0
   Tobermorite-11A 0 0
{extra_ocean_lines}

SELECTED_OUTPUT
    -file {base_name}.csv
    -reset false
    -reaction true
    -state True
    -solution True
    -step                 true
    -ph                   true
    -pe                   true
    -alkalinity           true
    -ionic_strength       true
    -water                true
    -totals              Al C Ca Fe H #K 
                         Mg #Mn 
                         N Na O P S Si #Cl
    -molalities          OH- H+ C2H4 C2H6 CH4 CO HCO3- CO3-2 CO2 
                         CH3COO- HCOO- HCN Ca+2 Cl- Fe+2 FeOH+ Fe+3 H2 #K+ 
                         Mg+2 #Mn+2 Mn+3 
                         NH4+ NH3 NH4CO3- N2 NO2- 
                         NO3- Na+ H2PO4- HPO4-2 PO4-3 HS- H2S SO3-2 HSO3- 
                         SO2 SO4-2 SiO2 HSiO3-
    -activities          OH- H+ C2H4 C2H6 CH4 CO HCO3- CO3-2 CO2 
                         CH3COO- HCOO- HCN Ca+2 Cl- Fe+2 FeOH+ Fe+3 H2 #K+ 
                         Mg+2 #Mn+2 Mn+3 
                         NH4+ NH3 NH4CO3- N2 NO2- 
                         NO3- Na+ H2PO4- HPO4-2 PO4-3 HS- H2S SO3-2 HSO3- 
                         SO2 SO4-2 SiO2 HSiO3-
    -kinetics            Enstatite Forsterite Troilite Albite Diopside Anorthite K-Feldspar Tephroite Schreibersite {iom_so_line}
    -saturation_indices  Enstatite Forsterite Troilite Albite Diopside Anorthite K-Feldspar Tephroite
    -equilibrium_phases  Fe Ni Lawrencite Analcime Anhydrite Aragonite Bassanite
                         Beidellite-Ca Beidellite-Fe Beidellite-Mg Beidellite-Na
                         Boehmite Brucite
                         C2H4(g) C2H6(g) C3H8(g) Calcite CH4(g) Chalcedony Chamosite Chrysotile
                         Citric_Acid Clinochlore-14A Clinochlore-7A Clinoptilolite-Ca Clinoptilolite-Na
                         Clinozoisite CO(g) CO2(g) Cronstedtite-7A
                         Daphnite-14A Daphnite-7A Dawsonite Diaspore Dolomite
                         Ettringite Fe(OH)2 Fe(OH)3 Fe2(SO4)3 FeSO4
                         Gibbsite Glycine Goethite Greenalite Gypsum Gyrolite
                         H2(g) H2S(g) Hematite Huntite Hydromagnesite Kaolinite
                         Lizardite Magnetite Melanterite MgSO4 Minnesotaite Mirabilite
                         Montmor-Ca Montmor-Mg Montmor-Na
                         N2(g) NH3(g) Nontronite-Ca Nontronite-Mg Nontronite-Na
                         Portlandite Pyrite
                         Saponite-Fe-Ca Saponite-Fe-Fe Saponite-Fe-Mg Saponite-Fe-Na
                         Saponite-Mg-Ca Saponite-Mg-Fe Saponite-Mg-Mg Saponite-Mg-Na
                         Sepiolite SiO2(am) SO2(g) Thenardite Tobermorite-11A
                         {extra_eq_so_line}
   -gases                CO2(g) H2(g) CH4(g) CO(g) H2O(g) H2S(g) N2(g) 
                         NH3(g) NO(g) NO2(g) O2(g) SO2(g)  

    
# USER_PUNCH is defined per-cycle below (with exact cumulative time embedded)
    ')

                                     # Pre-compute USER_PUNCH lines that don't change between cycles.
                                     # Total_Initial_g = full initial primary rock mass: 8 kinetic minerals
                                     # (tracked via KIN() each cycle) plus the 3 equilibrium primaries (Fe,
                                     # Ni, Lawrencite; fixed at their t=0 moles, since -equilibrium_phases
                                     # is deliberately not queried here — this is an initial-mass constant).
                                     total_initial_expr <- sprintf(
                                       "%.6f*100.3725 + %.6f*140.6715 + %.6f*87.913 + %.6f*262.1798 + %.6f*216.55 + %.6f*278.164 + %.6f*278.33 + %.6f*201.96 + %.6f*55.845 + %.6f*58.693 + %.6f*126.75 + %.6f*198.509",
                                       enstatite_moles, forsterite_moles, troilite_moles, albite_moles,
                                       diopside_moles, anorthite_moles, kfeldspar_moles, tephroite_moles,
                                       fe_moles, ni_moles, lawrencite_moles, schreibersite_moles
                                     )

                                     # ---- PATCH D1: punch IOM channel state so the restart can recover it ----
                                     # Without these columns the CSV has no record of remaining IOM, and
                                     # generate_phreeqc_input_tidal_kinetic_restart() cannot carry organics
                                     # over (it rebuilds KINETICS from the CSV). Generated from iom_pools so
                                     # a different channel set stays consistent automatically.
                                     iom_punch_headings <- ""
                                     iom_punch_lines    <- ""
                                     if (self$organic_wt_percent > 0 && !is.null(self$iom_pools)) {
                                       iom_punch_headings <- paste0(" ",
                                         paste(sprintf("%s_mol", self$iom_pools$name), collapse = " "))
                                       iom_punch_lines <- paste0(
                                         vapply(seq_len(nrow(self$iom_pools)), function(j) {
                                           sprintf("   %d PUNCH KIN(\"%s\")\n",
                                                   200 + 10 * j, self$iom_pools$name[j])
                                         }, character(1)), collapse = "")
                                     }

                                     # Generate unrolled cycle blocks.
                                     # Each cycle = 5 simulations (each terminated by END):
                                     #   1. React:    porewater + kinetics + eq_phases → SOLUTION 1 (reacted porewater)
                                     #   2. Extract:  x_step of SOLUTION 1 → SOLUTION 3 (flowing fluid, fractional mass=f_x)
                                     #   3. Combine:  SOLUTION 2 + SOLUTION 3 → SOLUTION 4 (combined ocean, fractional mass=f_ocean+f_fluid)
                                     #   4. Resample: 1-x_step of SOLUTION 1
                                     #              + f_x / (f_ocean + f_x)) of SOLUTION 4 → SOLUTION 1 (fractional mass=f_pore)
                                     #   5. Update ocean: f_ocean/(f_ocean+f_fluid) of SOLUTION 4 → SOLUTION 2 (fractional mass=f_ocean)
                                     cycle_blocks <- paste(sapply(seq_len(kinetic_cycles), function(i) {
                                       cumulative_time <- i * time_step_years
                                       paste0(
                                         sprintf("\n# --- Cycle %d of %d (cumulative time: %.4e years) ---\n", i, kinetic_cycles, cumulative_time),
                                         # Simulation 1: react porewater
                                         "USE solution 1\n",
                                         "USE kinetics 1\n",
                                         "USE equilibrium_phases 1\n",
                                         "USE reaction none\n",
                                         "USE reaction_temperature 1\n",
                                         "SAVE solution 1\n",
                                         "SAVE equilibrium_phases 1\n",
                                         "USER_PUNCH\n",
                                         "   -headings Time_Years Enst_remain_g Forst_remain_g Troil_remain_g Alb_remain_g Diop_remain_g Anorth_remain_g Kfs_remain_g Teph_remain_g Schr_remain_g Total_Initial_g",
                                         iom_punch_headings, "\n",
                                         sprintf("   10 PUNCH %.6f  # Cumulative time (years)\n", cumulative_time),
                                         "   20 PUNCH KIN(\"Enstatite\") * 100.3725\n",
                                         "   30 PUNCH KIN(\"Forsterite\") * 140.6715\n",
                                         "   40 PUNCH KIN(\"Troilite\") * 87.913\n",
                                         "   50 PUNCH KIN(\"Albite\") * 262.1798\n",
                                         "   60 PUNCH KIN(\"Diopside\") * 216.55\n",
                                         "   70 PUNCH KIN(\"Anorthite\") * 278.164\n",
                                         "   80 PUNCH KIN(\"K-Feldspar\") * 278.33\n",
                                         "   85 PUNCH KIN(\"Tephroite\") * 201.96\n",
                                         "   90 PUNCH KIN(\"Schreibersite\") * 198.509\n",
                                         sprintf("   100 total_initial = %s\n", total_initial_expr),
                                         "   110 PUNCH total_initial\n",
                                         iom_punch_lines,
                                         "END\n\n",
                                         # Simulation 2: ocean absorbs flowing fluid directly + seafloor equilibration
                                         # SOLUTION 4: f_fluid/f_pore of porewater + ocean, fractional mass = f_ocean + f_fluid
                                         "USE kinetics none\n",
                                         "USE equilibrium_phases 2\n",
                                         "USE reaction none\n",
                                         "USE reaction_temperature 2\n",
                                         "MIX 4\n",
                                          sprintf("  1  %.6e\n", x_step),
                                          "  2  1.0\n",
                                         "SAVE solution 4\n",
                                         "SAVE equilibrium_phases 2\n",
                                         "END\n\n",
                                         # Simulation 3: resample porewater — replace flowing fluid with combined ocean water
                                         # (f_pore-f_x)/f_pore of reacted porewater + f_x/(f_ocean+f_x) of combined ocean
                                         # fractional mass = f_pore
                                         "USE kinetics none\n",
                                         "USE equilibrium_phases none\n",
                                         "USE reaction none\n",
                                         "USE reaction_temperature 1\n",
                                         "MIX 1\n",
                                          sprintf("  1  %.6e\n", 1 - x_step),
                                          sprintf("  4  %.6e\n", f_x / (f_ocean + f_x)),
                                         "SAVE solution 1\n",
                                         "END\n\n",
                                         # Simulation 4: update ocean from combined ocean, fractional mass = f_ocean
                                         "USE kinetics none\n",
                                         "USE equilibrium_phases none\n",
                                         "USE reaction none\n",
                                         "MIX 2\n",
                                          sprintf("  4  %.6e\n", f_ocean / (f_ocean + f_x)),
                                         "SAVE solution 2\n",
                                         "END\n",
                                         # Early-exit: delay STOP by one full tidal cycle after Fayalite depletion.
                                         # Line 20 is the REACT/MIX discriminator: KIN("Pyrrhotite") > 0 only
                                         # during REACT (USE kinetics 1); it returns 0 in all MIX simulations
                                         # (USE kinetics none), so MIX steps are never interrupted.
                                         # First REACT with Fayalite = 0: PUT(1,1) flag, continue normally.
                                         # Second REACT with Fayalite = 0: EXISTS(1) → STOP.
                                         if (i < kinetic_cycles)
                                           paste0(
                                             "USER_PRINT\n",
                                             "   10 IF KIN(\"Fayalite\") > 0 THEN GOTO 30\n",
                                             "   20 IF KIN(\"Pyrrhotite\") <= 0 THEN GOTO 30\n",
                                             "   25 IF EXISTS(1) THEN STOP\n",
                                             "   27 PUT(1, 1)\n",
                                             "   30 REM\n"
                                           )
                                         else ""
                                       )
                                     }), collapse = "\n")

                                     return(self$apply_ch4_mode(paste0(template, "\n", cycle_blocks), ch4_mode))
                                   },

                                   # ---------------------------------------------------------------
                                   # Equilibrium-porewater variant: primary minerals dissolve via
                                   # EQUILIBRIUM_PHASES instead of kinetics.  Everything else
                                   # (mixing, ocean, seafloor EQ2) is identical.
                                   # ---------------------------------------------------------------
                                   generate_phreeqc_input_tidal_equil = function(output_dir = "phreeqc_inputs",
                                                                           n_orbital_cycles = NULL,
                                                                           mode = NULL,
                                                                           n_circulations = NULL,
                                                                           grouping_factor = NULL,
                                                                           suffix = "") {
                                     if (is.null(n_orbital_cycles))
                                       stop("generate_phreeqc_input_tidal_equil requires n_orbital_cycles (mode/n_circulations/grouping_factor are unused)")
                                     dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
                                     filename <- file.path(output_dir, sprintf(
                                       "enceladus_k%.2e_d%g_orb%g_equil%s.pqi",
                                       self$simulator$params$k,
                                       self$simulator$params$layer_thickness,
                                       n_orbital_cycles, suffix
                                     ))
                                     content <- self$create_phreeqc_template_tidal_equil(filename, n_orbital_cycles)
                                     writeLines(content, filename)
                                     return(filename)
                                   },

                                   create_phreeqc_template_tidal_equil = function(filename, n_orbital_cycles) {
                                     metrics         <- self$simulator$get_summary_metrics()
                                     orbital_period  <- self$simulator$config$orbital_period
                                     time_step_years <- orbital_period / (365.25 * 24 * 3600)
                                     kinetic_cycles  <- n_orbital_cycles
                                     total_sim_years <- kinetic_cycles * time_step_years

                                     # Fluid masses and mixing fractions (identical to kinetic version)
                                     porewater_mass      <- self$phreeqc_params$porewater_mass
                                     ocean_mass_actual   <- self$phreeqc_params$ocean_mass
                                     exchanging_fluid_mass <- metrics$exchanging_fluid_mass
                                     m_total  <- porewater_mass + ocean_mass_actual
                                     f_fluid  <- min(exchanging_fluid_mass / m_total,
                                                     porewater_mass / m_total)
                                     f_pore   <- porewater_mass    / m_total
                                     f_ocean  <- ocean_mass_actual / m_total
                                     porewater_mass_norm <- porewater_mass / self$phreeqc_params$layer_mass
                                     ocean_mass_norm     <- self$phreeqc_params$water_mass

                                     base_name <- tools::file_path_sans_ext(basename(filename))

                                     # Temperatures and pressures
                                     ocean_temp_C <- 274.6 - 273.15
                                     if (!is.null(self$profile_file) && file.exists(self$profile_file)) {
                                       pore_temp_C        <- mean_layer_temp_C(
                                         self$profile_file,
                                         self$simulator$params$layer_thickness
                                       )
                                       pore_pressure_atm  <- mean_layer_pressure_MPa(
                                         self$profile_file,
                                         self$simulator$params$layer_thickness
                                       ) / 0.101325
                                       ocean_pressure_atm <- ocean_floor_pressure_MPa(
                                         self$profile_file
                                       ) / 0.101325
                                     } else {
                                       if (!is.null(self$profile_file))
                                         warning("Profile file not found: '", self$profile_file, "' — falling back to defaults.")
                                       pore_temp_C        <- 25
                                       pore_pressure_atm  <- 70
                                       ocean_pressure_atm <- 70
                                     }

                                     # Mineral moles (ROCK_MOL_PER_KG, shared with the other generators; all 11
                                     # are equilibrium here, no kinetics at all). Magnetite is no longer primary.
                                     # SCHREIBERSITE IS NOT INCLUDED in this equilibrium-only generator. The
                                     # kinetic generators add it via KINETICS -formula, but an equilibrium
                                     # phase needs a log K, and neither database has an Fe3P phase. So this
                                     # generator's rock lacks 0.047 mmol P/g and ~0.26 mmol H2/g of reducing
                                     # capacity relative to the kinetic ones. Options (PI decision): add the
                                     # Fe3P elements as an instantaneous REACTION (full dissolution), or add a
                                     # Schreibersite phase with a sourced log K.
                                     mineral_scale_factor <- (100 - self$organic_wt_percent) / 100
                                     enstatite_moles  <- ROCK_MOL_PER_KG[["enstatite"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     forsterite_moles <- ROCK_MOL_PER_KG[["forsterite"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     troilite_moles   <- ROCK_MOL_PER_KG[["troilite"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     albite_moles     <- ROCK_MOL_PER_KG[["albite"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     diopside_moles   <- ROCK_MOL_PER_KG[["diopside"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     anorthite_moles  <- ROCK_MOL_PER_KG[["anorthite"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     kfeldspar_moles  <- ROCK_MOL_PER_KG[["kfeldspar"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     tephroite_moles  <- ROCK_MOL_PER_KG[["tephroite"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     fe_moles         <- ROCK_MOL_PER_KG[["fe"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     ni_moles         <- ROCK_MOL_PER_KG[["ni"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     lawrencite_moles <- ROCK_MOL_PER_KG[["lawrencite"]] * self$phreeqc_params$total_rock * mineral_scale_factor

                                     total_initial_expr <- sprintf(
                                       "%.6f*100.3725 + %.6f*140.6715 + %.6f*87.913 + %.6f*262.1798 + %.6f*216.55 + %.6f*278.164 + %.6f*278.33 + %.6f*201.96 + %.6f*55.845 + %.6f*58.693 + %.6f*126.75",
                                       enstatite_moles, forsterite_moles, troilite_moles, albite_moles,
                                       diopside_moles, anorthite_moles, kfeldspar_moles, tephroite_moles,
                                       fe_moles, ni_moles, lawrencite_moles
                                     )

                                     template <- glue::glue('
TITLE Enceladus hydrothermal simulation — equilibrium porewater

KNOBS
   -logfile true

# k={self$simulator$params$k}, thickness={self$simulator$params$layer_thickness} m
# Mode: equilibrium (primary minerals in EQUILIBRIUM_PHASES, no kinetics)
# Simulation parameters:
# Layer thickness: {self$simulator$params$layer_thickness} m ({sprintf("%.1f", self$simulator$params$layer_thickness/self$simulator$config$core_radius*100)}% of core)
# Layer mass: {sprintf("%.2e", self$phreeqc_params$layer_mass)} kg  (total_rock normalised to 1.0)
# Reaction time: {sprintf("%.2e", total_sim_years)} years
# Total cycles: {kinetic_cycles}
# Time step: {sprintf("%.2e", time_step_years)} years per step
# Mean flow rate: {sprintf("%.2e", metrics$mean_flow_rate_Ls)} L/s
# Circulation time: {sprintf("%.2e", metrics$circulation_time_years)} years
# Complete circulation time: {sprintf("%.2e", metrics$complete_circulation_time)} years
# Corrected circulation time (f_fluid capped at f_pore): {sprintf("%.2e", metrics$corrected_circulation_time)} years
# Corrected complete circulation time: {sprintf("%.2e", metrics$corrected_complete_circulation_time)} years
# Flowing fluid mass (per orbital cycle): {sprintf("%.2e", metrics$exchanging_fluid_mass)} kg
# Porewater mass: {sprintf("%.2e", self$phreeqc_params$porewater_mass)} kg
# Ocean mass: {sprintf("%.2e", self$phreeqc_params$ocean_mass)} kg

SOLUTION 1 Enceladus porewater
    temp      {round(ocean_temp_C, 2)}#{round(pore_temp_C, 2)}
    pH        7
    pe        4
    redox     pe
    units     mol/kgw
    density   1
    water     {porewater_mass_norm}
    pressure  {round(ocean_pressure_atm, 1)}#{round(pore_pressure_atm, 1)}

SOLUTION 2 Enceladus ocean
    temp      {round(ocean_temp_C, 2)}
    pH        7
    pe        4
    redox     pe
    units     mol/kgw
    density   1
    water     {ocean_mass_norm}
    pressure  {round(ocean_pressure_atm, 1)}

REACTION_TEMPERATURE 1
    {round(ocean_temp_C, 2)}#{round(pore_temp_C, 2)}

REACTION_TEMPERATURE 2
    {round(ocean_temp_C, 2)}

EQUILIBRIUM_PHASES 1
# --- Primary minerals (dissolve to equilibrium each cycle) ---
   Enstatite            0  {enstatite_moles}   dissolve_only
   Forsterite           0  {forsterite_moles}  dissolve_only
   Troilite             0  {troilite_moles}    dissolve_only
   Albite               0  {albite_moles}      dissolve_only
   Diopside             0  {diopside_moles}    dissolve_only
   Anorthite            0  {anorthite_moles}   dissolve_only
   K-Feldspar           0  {kfeldspar_moles}   dissolve_only
   Tephroite            0  {tephroite_moles}   dissolve_only
   Fe                   0  {fe_moles}          dissolve_only
   Ni                   0  {ni_moles}          dissolve_only
   Lawrencite           0  {lawrencite_moles}  dissolve_only
# --- Secondary minerals (precipitate if supersaturated) ---
#   Alabandite 0 0
#   Alanine 0 0
#   Alum-K 0 0
#   Alunite 0 0
   Analcime 0 0
   Anhydrite 0 0
   Aragonite 0 0
#   Arcanite 0 0
#   Artinite 0 0
   Bassanite 0 0
   Beidellite-Ca 0 0
   Beidellite-Fe 0 0
#   Beidellite-K 0 0
   Beidellite-Mg 0 0
   Beidellite-Na 0 0
   Boehmite 0 0
   Brucite 0 0
   C2H4(g) 0 0
   C2H6(g) 0 0
   C3H8(g) 0 0
   Calcite 0 0
#   Celadonite 0 0
   CH4(g) 0 0
   Chalcedony 0 0
   Chamosite 0 0
   Chrysotile 0 0
   Citric_Acid 0 0
   Clinochlore-14A 0 0
   Clinochlore-7A 0 0
   Clinoptilolite-Ca 0 0
#   Clinoptilolite-K 0 0
   Clinoptilolite-Na 0 0
   Clinozoisite 0 0
   CO(g) 0 0
   CO2(g) 0 0
   Cronstedtite-7A 0 0
   Daphnite-14A 0 0
   Daphnite-7A 0 0
   Dawsonite 0 0
   Diaspore 0 0
   Dolomite 0 0
   Ettringite 0 0
   Fe(OH)2 0 0
   Fe(OH)3 0 0
   Fe2(SO4)3 0 0
   FeSO4 0 0
   Gibbsite 0 0
   Glycine 0 0
   Goethite 0 0
   Greenalite 0 0
   Gypsum 0 0
   Gyrolite 0 0
   H2(g) 0 0
#   H2O(g) 0 0
   H2S(g) 0 0
#   Halite 0 0
#   Hausmannite 0 0
   Hematite 0 0
   Huntite 0 0
   Hydromagnesite 0 0
#   Hydroxyapatite 0 0
#   Ice 0 0
#   Jarosite 0 0
#   KAl(SO4)2 0 0
   Kaolinite 0 0
#   KerogenC128 0 0
#   KerogenC292 0 0
#   KerogenC406 0 0
#   KerogenC415 0 0
#   KerogenC515 0 0
   Lizardite 0 0
   Melanterite 0 0
#   MgOHCl 0 0
   MgSO4 0 0
   Minnesotaite 0 0
   Mirabilite 0 0
#   Mn(OH)2(am) 0 0
#   MnCl2:2H2O 0 0
#   MnCl2:4H2O 0 0
#   MnCl2:H2O 0 0
#   MnSO4 0 0
#   Monohydrocalcite 0 0
   Montmor-Ca 0 0
#   Montmor-K 0 0
   Montmor-Mg 0 0
   Montmor-Na 0 0
   N2(g) 0 0
#   Na2CO3 0 0
#   Na2CO3:7H2O 0 0
#   Nahcolite 0 0
#   Natron 0 0
#   Nesquehonite 0 0
   NH3(g) 0 0
#   NH4Cl 0 0
#   NH4HCO3 0 0
#   Niter 0 0
#   NO(g) 0 0
#   NO2(g) 0 0
   Nontronite-Ca 0 0
#   Nontronite-K 0 0
   Nontronite-Mg 0 0
   Nontronite-Na 0 0
#   O2(g) 0 0
   Portlandite 0 0
#   Pyridine 0 0
   Pyrite 0 0
#   Rhodochrosite 0 0
   Saponite-Fe-Ca 0 0
   Saponite-Fe-Fe 0 0
#   Saponite-Fe-K 0 0
   Saponite-Fe-Mg 0 0
   Saponite-Fe-Na 0 0
   Saponite-Mg-Ca 0 0
   Saponite-Mg-Fe 0 0
#   Saponite-Mg-K 0 0
   Saponite-Mg-Mg 0 0
   Saponite-Mg-Na 0 0
   Sepiolite 0 0
#   Siderite 0 0
   SiO2(am) 0 0
#   Smectite-high-Fe-Mg 0 0
#   Smectite-low-Fe-Mg 0 0
   SO2(g) 0 0
#   Stilbite 0 0
#   Strengite 0 0
#   Sylvite 0 0
   Thenardite 0 0
#   Thermonatrite 0 0
   Tobermorite-11A 0 0

EQUILIBRIUM_PHASES 2
#   Alabandite 0 0
#   Alanine 0 0
#   Alum-K 0 0
#   Alunite 0 0
   Analcime 0 0
   Anhydrite 0 0
   Aragonite 0 0
#   Arcanite 0 0
#   Artinite 0 0
   Bassanite 0 0
   Beidellite-Ca 0 0
   Beidellite-Fe 0 0
#   Beidellite-K 0 0
   Beidellite-Mg 0 0
   Beidellite-Na 0 0
   Boehmite 0 0
   Brucite 0 0
   C2H4(g) 0 0
   C2H6(g) 0 0
   C3H8(g) 0 0
   Calcite 0 0
#   Celadonite 0 0
   CH4(g) 0 0
   Chalcedony 0 0
   Chamosite 0 0
   Chrysotile 0 0
   Citric_Acid 0 0
   Clinochlore-14A 0 0
   Clinochlore-7A 0 0
   Clinoptilolite-Ca 0 0
#   Clinoptilolite-K 0 0
   Clinoptilolite-Na 0 0
   Clinozoisite 0 0
   CO(g) 0 0
   CO2(g) 0 0
   Cronstedtite-7A 0 0
   Daphnite-14A 0 0
   Daphnite-7A 0 0
   Dawsonite 0 0
   Diaspore 0 0
   Dolomite 0 0
   Ettringite 0 0
   Fe(OH)2 0 0
   Fe(OH)3 0 0
   Fe2(SO4)3 0 0
   FeSO4 0 0
   Gibbsite 0 0
   Glycine 0 0
   Goethite 0 0
   Greenalite 0 0
   Gypsum 0 0
   Gyrolite 0 0
   H2(g) 0 0
#   H2O(g) 0 0
   H2S(g) 0 0
#   Halite 0 0
#   Hausmannite 0 0
   Hematite 0 0
   Huntite 0 0
   Hydromagnesite 0 0
#   Hydroxyapatite 0 0
#   Ice 0 0
#   Jarosite 0 0
#   KAl(SO4)2 0 0
   Kaolinite 0 0
#   KerogenC128 0 0
#   KerogenC292 0 0
#   KerogenC406 0 0
#   KerogenC415 0 0
#   KerogenC515 0 0
   Lizardite 0 0
   Magnetite 0 0
   Melanterite 0 0
#   MgOHCl 0 0
   MgSO4 0 0
   Minnesotaite 0 0
   Mirabilite 0 0
#   Mn(OH)2(am) 0 0
#   MnCl2:2H2O 0 0
#   MnCl2:4H2O 0 0
#   MnCl2:H2O 0 0
#   MnSO4 0 0
#   Monohydrocalcite 0 0
   Montmor-Ca 0 0
#   Montmor-K 0 0
   Montmor-Mg 0 0
   Montmor-Na 0 0
   N2(g) 0 0
#   Na2CO3 0 0
#   Na2CO3:7H2O 0 0
#   Nahcolite 0 0
#   Natron 0 0
#   Nesquehonite 0 0
   NH3(g) 0 0
#   NH4Cl 0 0
#   NH4HCO3 0 0
#   Niter 0 0
#   NO(g) 0 0
#   NO2(g) 0 0
   Nontronite-Ca 0 0
#   Nontronite-K 0 0
   Nontronite-Mg 0 0
   Nontronite-Na 0 0
#   O2(g) 0 0
   Portlandite 0 0
#   Pyridine 0 0
   Pyrite 0 0
#   Rhodochrosite 0 0
   Saponite-Fe-Ca 0 0
   Saponite-Fe-Fe 0 0
#   Saponite-Fe-K 0 0
   Saponite-Fe-Mg 0 0
   Saponite-Fe-Na 0 0
   Saponite-Mg-Ca 0 0
   Saponite-Mg-Fe 0 0
#   Saponite-Mg-K 0 0
   Saponite-Mg-Mg 0 0
   Saponite-Mg-Na 0 0
   Sepiolite 0 0
#   Siderite 0 0
   SiO2(am) 0 0
#   Smectite-high-Fe-Mg 0 0
#   Smectite-low-Fe-Mg 0 0
   SO2(g) 0 0
#   Stilbite 0 0
#   Strengite 0 0
#   Sylvite 0 0
   Thenardite 0 0
#   Thermonatrite 0 0
   Tobermorite-11A 0 0

SELECTED_OUTPUT
    -file {base_name}.csv
    -reset false
    -reaction true
    -state True
    -solution True
    -step                 true
    -ph                   true
    -pe                   true
    -alkalinity           true
    -ionic_strength       true
    -water                true
    -totals              Al C Ca Fe H Mg N Na O P S Si
    -molalities          OH- H+ C2H4 C2H6 CH4 CO HCO3- CO3-2 CO2
                         CH3COO- HCOO- HCN Ca+2 Cl- Fe+2 FeOH+ Fe+3 H2
                         Mg+2 NH4+ NH3 NH4CO3- N2 NO2-
                         NO3- Na+ H2PO4- HPO4-2 PO4-3 HS- H2S SO3-2 HSO3-
                         SO2 SO4-2 SiO2 HSiO3-
    -activities          OH- H+ C2H4 C2H6 CH4 CO HCO3- CO3-2 CO2
                         CH3COO- HCOO- HCN Ca+2 Cl- Fe+2 FeOH+ Fe+3 H2
                         Mg+2 NH4+ NH3 NH4CO3- N2 NO2-
                         NO3- Na+ H2PO4- HPO4-2 PO4-3 HS- H2S SO3-2 HSO3-
                         SO2 SO4-2 SiO2 HSiO3-
    -equilibrium_phases  Enstatite Forsterite Troilite Albite Diopside Anorthite K-Feldspar Tephroite
                         Fe Ni Lawrencite Analcime Anhydrite Aragonite Bassanite
                         Beidellite-Ca Beidellite-Fe Beidellite-Mg Beidellite-Na
                         Boehmite Brucite
                         C2H4(g) C2H6(g) C3H8(g) Calcite CH4(g) Chalcedony Chamosite Chrysotile
                         Citric_Acid Clinochlore-14A Clinochlore-7A Clinoptilolite-Ca Clinoptilolite-Na
                         Clinozoisite CO(g) CO2(g) Cronstedtite-7A
                         Daphnite-14A Daphnite-7A Dawsonite Diaspore Dolomite
                         Ettringite Fe(OH)2 Fe(OH)3 Fe2(SO4)3 FeSO4
                         Gibbsite Glycine Goethite Greenalite Gypsum Gyrolite
                         H2(g) H2S(g) Hematite Huntite Hydromagnesite Kaolinite
                         Lizardite Magnetite Melanterite MgSO4 Minnesotaite Mirabilite
                         Montmor-Ca Montmor-Mg Montmor-Na
                         N2(g) NH3(g) Nontronite-Ca Nontronite-Mg Nontronite-Na
                         Portlandite Pyrite
                         Saponite-Fe-Ca Saponite-Fe-Fe Saponite-Fe-Mg Saponite-Fe-Na
                         Saponite-Mg-Ca Saponite-Mg-Fe Saponite-Mg-Mg Saponite-Mg-Na
                         Sepiolite SiO2(am) SO2(g) Thenardite Tobermorite-11A
    -saturation_indices  Enstatite Forsterite Troilite Albite Diopside Anorthite K-Feldspar Tephroite Fe Ni Lawrencite
    -gases               CO2(g) H2(g) CH4(g) CO(g) H2O(g) H2S(g) N2(g)
                         NH3(g) NO(g) NO2(g) O2(g) SO2(g)

# USER_PUNCH is defined per-cycle below (with exact cumulative time embedded)
    ')

                                     cycle_blocks <- paste(sapply(seq_len(kinetic_cycles), function(i) {
                                       cumulative_time <- i * time_step_years
                                       paste0(
                                         sprintf("\n# --- Cycle %d of %d (cumulative time: %.4e years) ---\n", i, kinetic_cycles, cumulative_time),
                                         # Simulation 1: equilibrate porewater with EQ phases (primary + secondary)
                                         "USE solution 1\n",
                                         "USE kinetics none\n",
                                         "USE equilibrium_phases 1\n",
                                         "USE reaction none\n",
                                         "USE reaction_temperature 1\n",
                                         "SAVE solution 1\n",
                                         "SAVE equilibrium_phases 1\n",
                                         "USER_PUNCH\n",
                                         "   -headings Time_Years Enst_remain_g Forst_remain_g Troil_remain_g Alb_remain_g Diop_remain_g Anorth_remain_g Kfs_remain_g Teph_remain_g Total_Initial_g\n",
                                         sprintf("   10 PUNCH %.6f\n", cumulative_time),
                                         "   20 PUNCH EQUI(\"Enstatite\") * 100.3725\n",
                                         "   30 PUNCH EQUI(\"Forsterite\") * 140.6715\n",
                                         "   40 PUNCH EQUI(\"Troilite\") * 87.913\n",
                                         "   50 PUNCH EQUI(\"Albite\") * 262.1798\n",
                                         "   60 PUNCH EQUI(\"Diopside\") * 216.55\n",
                                         "   70 PUNCH EQUI(\"Anorthite\") * 278.164\n",
                                         "   80 PUNCH EQUI(\"K-Feldspar\") * 278.33\n",
                                         "   85 PUNCH EQUI(\"Tephroite\") * 201.96\n",
                                         sprintf("   100 total_initial = %s\n", total_initial_expr),
                                         "   110 PUNCH total_initial\n",
                                         "END\n\n",
                                         # Simulation 2: ocean absorbs flowing fluid directly + seafloor equilibration
                                         # SOLUTION 4: f_fluid/f_pore of porewater + ocean, mass = f_ocean + f_fluid
                                         "USE kinetics none\n",
                                         "USE equilibrium_phases 2\n",
                                         "USE reaction none\n",
                                         "USE reaction_temperature 2\n",
                                         "MIX 4\n",
                                         sprintf("  1  %.6e\n", f_fluid / f_pore),
                                         "  2  1.0\n",
                                         "SAVE solution 4\n",
                                         "SAVE equilibrium_phases 2\n",
                                         "END\n\n",
                                         # Simulation 3: resample porewater
                                         "USE kinetics none\n",
                                         "USE equilibrium_phases none\n",
                                         "USE reaction none\n",
                                         "USE reaction_temperature 1\n",
                                         "MIX 1\n",
                                         sprintf("  1  %.6e\n", (f_pore - f_fluid) / f_pore),
                                         sprintf("  4  %.6e\n", f_fluid / (f_ocean + f_fluid)),
                                         "SAVE solution 1\n",
                                         "END\n\n",
                                         # Simulation 4: update ocean
                                         "USE kinetics none\n",
                                         "USE equilibrium_phases none\n",
                                         "USE reaction none\n",
                                         "MIX 2\n",
                                         sprintf("  4  %.6e\n", f_ocean / (f_ocean + f_fluid)),
                                         "SAVE solution 2\n",
                                         "END\n"
                                       )
                                     }), collapse = "\n")

                                     return(paste0(template, "\n", cycle_blocks))
                                   },

                                   # Single-step bulk equilibrium: entire hydrosphere (ocean + porewater)
                                   # equilibrated against the full permeable rock layer in one PHREEQC step.
                                   # No permeability or fluid-flow parameters are used.
                                   generate_phreeqc_input_bulk_equil = function(output_dir = "phreeqc_inputs",
                                                                                        suffix = "") {
                                     dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

                                     filename <- file.path(output_dir, sprintf(
                                       "enceladus_d%g_hydrosphere_equil%s.pqi",
                                       self$simulator$params$layer_thickness, suffix
                                     ))

                                     # Combined hydrosphere mass, normalised to layer mass
                                     porewater_mass_norm <- self$phreeqc_params$porewater_mass / self$phreeqc_params$layer_mass
                                     ocean_mass_norm     <- self$phreeqc_params$water_mass   # already = ocean_mass / layer_mass
                                     combined_mass_norm  <- porewater_mass_norm + ocean_mass_norm

                                     # Temperature and pressure (ocean floor; no pore-gradient averaging needed)
                                     ocean_temp_C <- 274.6 - 273.15
                                     if (!is.null(self$profile_file) && file.exists(self$profile_file)) {
                                       ocean_pressure_atm <- ocean_floor_pressure_MPa(self$profile_file) / 0.101325
                                     } else {
                                       if (!is.null(self$profile_file))
                                         warning("Profile file not found: '", self$profile_file, "' — falling back to defaults.")
                                       ocean_pressure_atm <- 70
                                     }

                                     # Primary mineral moles (ROCK_MOL_PER_KG, shared with the other generators;
                                     # all 11 are equilibrium here since this function has no kinetics at all).
                                     # Magnetite is no longer primary — ordinary secondary phase only.
                                     # SCHREIBERSITE IS NOT INCLUDED here: no Fe3P phase/log K in either
                                     # database, so it can't be an equilibrium primary. See the matching
                                     # note in generate_phreeqc_input_tidal_equil.
                                     mineral_scale_factor <- (100 - self$organic_wt_percent) / 100
                                     enstatite_moles  <- ROCK_MOL_PER_KG[["enstatite"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     forsterite_moles <- ROCK_MOL_PER_KG[["forsterite"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     troilite_moles   <- ROCK_MOL_PER_KG[["troilite"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     albite_moles     <- ROCK_MOL_PER_KG[["albite"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     diopside_moles   <- ROCK_MOL_PER_KG[["diopside"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     anorthite_moles  <- ROCK_MOL_PER_KG[["anorthite"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     kfeldspar_moles  <- ROCK_MOL_PER_KG[["kfeldspar"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     tephroite_moles  <- ROCK_MOL_PER_KG[["tephroite"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     fe_moles         <- ROCK_MOL_PER_KG[["fe"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     ni_moles         <- ROCK_MOL_PER_KG[["ni"]] * self$phreeqc_params$total_rock * mineral_scale_factor
                                     lawrencite_moles <- ROCK_MOL_PER_KG[["lawrencite"]] * self$phreeqc_params$total_rock * mineral_scale_factor

                                     total_initial_expr <- sprintf(
                                       "%.6f*100.3725 + %.6f*140.6715 + %.6f*87.913 + %.6f*262.1798 + %.6f*216.55 + %.6f*278.164 + %.6f*278.33 + %.6f*201.96 + %.6f*55.845 + %.6f*58.693 + %.6f*126.75",
                                       enstatite_moles, forsterite_moles, troilite_moles, albite_moles,
                                       diopside_moles, anorthite_moles, kfeldspar_moles, tephroite_moles,
                                       fe_moles, ni_moles, lawrencite_moles
                                     )

                                     base_name <- tools::file_path_sans_ext(basename(filename))

                                     content <- glue::glue('
TITLE Enceladus bulk hydrosphere equilibrium

KNOBS
   -logfile true

# thickness={self$simulator$params$layer_thickness} m
# Mode: single-step equilibrium — full hydrosphere (ocean + porewater) vs permeable rock layer
# Layer mass: {sprintf("%.2e", self$phreeqc_params$layer_mass)} kg  (total_rock normalised to 1.0)
# Porewater mass: {sprintf("%.2e", self$phreeqc_params$porewater_mass)} kg
# Ocean mass: {sprintf("%.2e", self$phreeqc_params$ocean_mass)} kg
# Combined hydrosphere (normalised): {sprintf("%.6f", combined_mass_norm)}
# Temperature: {round(ocean_temp_C, 2)} C (ocean floor)

SOLUTION 1 Enceladus hydrosphere (ocean + porewater combined)
    temp      {round(ocean_temp_C, 2)}
    pH        7
    pe        4
    redox     pe
    units     mol/kgw
    density   1
    water     {combined_mass_norm}
    pressure  {round(ocean_pressure_atm, 1)}

EQUILIBRIUM_PHASES 1
# --- Primary minerals (dissolve to equilibrium) ---
   Enstatite            0  {enstatite_moles}   dissolve_only
   Forsterite           0  {forsterite_moles}  dissolve_only
   Troilite             0  {troilite_moles}    dissolve_only
   Albite               0  {albite_moles}      dissolve_only
   Diopside             0  {diopside_moles}    dissolve_only
   Anorthite            0  {anorthite_moles}   dissolve_only
   K-Feldspar           0  {kfeldspar_moles}   dissolve_only
   Tephroite            0  {tephroite_moles}   dissolve_only
   Fe                   0  {fe_moles}          dissolve_only
   Ni                   0  {ni_moles}          dissolve_only
   Lawrencite           0  {lawrencite_moles}  dissolve_only
# --- Secondary minerals (precipitate if supersaturated) ---
   Analcime 0 0
   Anhydrite 0 0
   Aragonite 0 0
   Bassanite 0 0
   Beidellite-Ca 0 0
   Beidellite-Fe 0 0
   Beidellite-Mg 0 0
   Beidellite-Na 0 0
   Boehmite 0 0
   Brucite 0 0
   C2H4(g) 0 0
   C2H6(g) 0 0
   C3H8(g) 0 0
   Calcite 0 0
   CH4(g) 0 0
   Chalcedony 0 0
   Chamosite 0 0
   Chrysotile 0 0
   Citric_Acid 0 0
   Clinochlore-14A 0 0
   Clinochlore-7A 0 0
   Clinoptilolite-Ca 0 0
   Clinoptilolite-Na 0 0
   Clinozoisite 0 0
   CO(g) 0 0
   CO2(g) 0 0
   Cronstedtite-7A 0 0
   Daphnite-14A 0 0
   Daphnite-7A 0 0
   Dawsonite 0 0
   Diaspore 0 0
   Dolomite 0 0
   Ettringite 0 0
   Fe(OH)2 0 0
   Fe(OH)3 0 0
   Fe2(SO4)3 0 0
   FeSO4 0 0
   Gibbsite 0 0
   Glycine 0 0
   Goethite 0 0
   Greenalite 0 0
   Gypsum 0 0
   Gyrolite 0 0
   H2(g) 0 0
   H2S(g) 0 0
   Hematite 0 0
   Huntite 0 0
   Hydromagnesite 0 0
   Kaolinite 0 0
   Lizardite 0 0
   Magnetite 0 0
   Melanterite 0 0
   MgSO4 0 0
   Minnesotaite 0 0
   Mirabilite 0 0
   Montmor-Ca 0 0
   Montmor-Mg 0 0
   Montmor-Na 0 0
   N2(g) 0 0
   NH3(g) 0 0
   Nontronite-Ca 0 0
   Nontronite-Mg 0 0
   Nontronite-Na 0 0
   Portlandite 0 0
   Pyrite 0 0
   Saponite-Fe-Ca 0 0
   Saponite-Fe-Fe 0 0
   Saponite-Fe-Mg 0 0
   Saponite-Fe-Na 0 0
   Saponite-Mg-Ca 0 0
   Saponite-Mg-Fe 0 0
   Saponite-Mg-Mg 0 0
   Saponite-Mg-Na 0 0
   Sepiolite 0 0
   SiO2(am) 0 0
   SO2(g) 0 0
   Thenardite 0 0
   Tobermorite-11A 0 0

SELECTED_OUTPUT
    -file {base_name}.csv
    -reset false
    -state true
    -solution true
    -step true
    -ph true
    -pe true
    -alkalinity true
    -ionic_strength true
    -water true
    -totals              Al C Ca Fe H Mg N Na O P S Si
    -molalities          OH- H+ C2H4 C2H6 CH4 CO HCO3- CO3-2 CO2
                         CH3COO- HCOO- HCN Ca+2 Cl- Fe+2 FeOH+ Fe+3 H2
                         Mg+2 NH4+ NH3 NH4CO3- N2 NO2-
                         NO3- Na+ H2PO4- HPO4-2 PO4-3 HS- H2S SO3-2 HSO3-
                         SO2 SO4-2 SiO2 HSiO3-
    -activities          OH- H+ C2H4 C2H6 CH4 CO HCO3- CO3-2 CO2
                         CH3COO- HCOO- HCN Ca+2 Cl- Fe+2 FeOH+ Fe+3 H2
                         Mg+2 NH4+ NH3 NH4CO3- N2 NO2-
                         NO3- Na+ H2PO4- HPO4-2 PO4-3 HS- H2S SO3-2 HSO3-
                         SO2 SO4-2 SiO2 HSiO3-
    -equilibrium_phases  Enstatite Forsterite Troilite Albite Diopside Anorthite K-Feldspar Tephroite
                         Fe Ni Lawrencite
                         Analcime Anhydrite Aragonite Bassanite Beidellite-Ca Beidellite-Fe Beidellite-Mg Beidellite-Na
                         Boehmite Brucite C2H4(g) C2H6(g) C3H8(g) Calcite CH4(g) Chalcedony Chamosite Chrysotile
                         Citric_Acid Clinochlore-14A Clinochlore-7A Clinoptilolite-Ca Clinoptilolite-Na
                         Clinozoisite CO(g) CO2(g) Cronstedtite-7A
                         Daphnite-14A Daphnite-7A Dawsonite Diaspore Dolomite
                         Ettringite Fe(OH)2 Fe(OH)3 Fe2(SO4)3 FeSO4 Gibbsite Glycine Goethite Greenalite Gypsum Gyrolite
                         H2(g) H2S(g) Hematite Huntite Hydromagnesite Kaolinite Lizardite Magnetite
                         Melanterite MgSO4 Minnesotaite Mirabilite
                         Montmor-Ca Montmor-Mg Montmor-Na N2(g) NH3(g) Nontronite-Ca Nontronite-Mg Nontronite-Na
                         Portlandite Pyrite
                         Saponite-Fe-Ca Saponite-Fe-Fe Saponite-Fe-Mg Saponite-Fe-Na
                         Saponite-Mg-Ca Saponite-Mg-Fe Saponite-Mg-Mg Saponite-Mg-Na
                         Sepiolite SiO2(am) SO2(g) Thenardite Tobermorite-11A
    -saturation_indices  Enstatite Forsterite Troilite Albite Diopside Anorthite K-Feldspar Tephroite Fe Ni Lawrencite
    -gases               CO2(g) H2(g) CH4(g) CO(g) H2O(g) H2S(g) N2(g)
                         NH3(g) NO(g) NO2(g) O2(g) SO2(g)

USER_PUNCH
   -headings Time_Years Enst_remain_g Forst_remain_g Troil_remain_g Alb_remain_g Diop_remain_g Anorth_remain_g Kfs_remain_g Teph_remain_g Total_Initial_g
   10 PUNCH 0.0
   20 PUNCH EQUI("Enstatite") * 100.3725
   30 PUNCH EQUI("Forsterite") * 140.6715
   40 PUNCH EQUI("Troilite") * 87.913
   50 PUNCH EQUI("Albite") * 262.1798
   60 PUNCH EQUI("Diopside") * 216.55
   70 PUNCH EQUI("Anorthite") * 278.164
   80 PUNCH EQUI("K-Feldspar") * 278.33
   85 PUNCH EQUI("Tephroite") * 201.96
   100 total_initial = {total_initial_expr}
   110 PUNCH total_initial

END
')

                                     writeLines(content, filename)
                                     return(filename)
                                   },

                                   generate_phreeqc_input_bulk_kinetic = function(
                                     output_dir          = "phreeqc_inputs",
                                     n_orbital_cycles    = NULL,
                                     mode                = c("porewater", "hydrosphere"),
                                     water_rock_ratio    = NULL,
                                     porewater_temp      = "ocean",
                                     suffix              = "",
                                     CH4_redox_override  = NULL
                                   ) {
                                     mode <- match.arg(mode)
                                     if (mode == "porewater" && (is.null(water_rock_ratio) || water_rock_ratio <= 0))
                                       stop("water_rock_ratio must be a positive number when mode = 'porewater'")
                                     if (is.null(n_orbital_cycles) || n_orbital_cycles <= 0)
                                       stop("n_orbital_cycles must be a positive integer > 0")

                                     dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

                                     temp_tag <- if (is.numeric(porewater_temp)) sprintf("T%g", porewater_temp) else sprintf("T%s", porewater_temp)
                                     org_tag  <- sprintf("org%g", self$organic_wt_percent)

                                     filename <- if (mode == "porewater") {
                                       file.path(output_dir, sprintf(
                                         "enceladus_wr%.4g_%s_%s_orb%g_porewater_kinetic%s.pqi",
                                         water_rock_ratio, temp_tag, org_tag, n_orbital_cycles, suffix
                                       ))
                                     } else {
                                       file.path(output_dir, sprintf(
                                         "enceladus_d%g_%s_%s_orb%g_hydrosphere_kinetic%s.pqi",
                                         self$simulator$params$layer_thickness, temp_tag, org_tag, n_orbital_cycles, suffix
                                       ))
                                     }

                                     # Temperature and pressure (ocean floor)
                                     ocean_temp_C <- 274.6 - 273.15
                                     if (!is.null(self$profile_file) && file.exists(self$profile_file)) {
                                       pore_temp_C        <- mean_layer_temp_C(self$profile_file, self$simulator$params$layer_thickness)
                                       pore_pressure_atm  <- mean_layer_pressure_MPa(self$profile_file, self$simulator$params$layer_thickness) / 0.101325
                                       ocean_pressure_atm <- ocean_floor_pressure_MPa(self$profile_file) / 0.101325
                                     } else {
                                       if (!is.null(self$profile_file))
                                         warning("Profile file not found: '", self$profile_file, "' — falling back to defaults.")
                                       pore_temp_C        <- 25
                                       pore_pressure_atm  <- 70
                                       ocean_pressure_atm <- 70
                                     }

                                     # Resolve porewater_temp to a concrete value — this is the reaction
                                     # temperature for the single combined solution used here.
                                     porewater_temp_C <- if (is.numeric(porewater_temp)) porewater_temp else
                                       if (porewater_temp == "profile") pore_temp_C else ocean_temp_C

                                     # CH4 (redox-coupled) vs Mtg (decoupled) — see resolve_ch4_mode().
                                     ch4_mode <- self$resolve_ch4_mode(porewater_temp_C, CH4_redox_override)
                                     ch4_override_note <- if (!is.null(CH4_redox_override)) sprintf(", override=%s", CH4_redox_override) else ""
                                     ch4_db_note_line <- sprintf(
                                       "# Recommended database: %s (CH4 %s; porewater T = %.1f C%s)",
                                       ch4_mode$database,
                                       if (ch4_mode$decoupled) "decoupled -> aqueous/gas species are Mtg/Mtg(g)" else "coupled -> aqueous/gas species are CH4/CH4(g)",
                                       porewater_temp_C, ch4_override_note
                                     )

                                     # Clathrate stability at porewater conditions — see resolve_clathrate_stability().
                                     # Single combined solution here uses ocean_pressure_atm as its reaction
                                     # pressure (pore_pressure_atm is computed above but not otherwise used
                                     # in this function), so the stability check matches that.
                                     clathrate_note_lines <- self$clathrate_stability_note(porewater_temp_C, ocean_pressure_atm)

                                     # Temperature-gated secondary minerals (single solution, so a single
                                     # temperature — porewater_temp_C — gates the whole block). Stable
                                     # clathrate phases (checked against ocean_pressure_atm, matching this
                                     # function's own reaction pressure) are appended the same way.
                                     extra_ocean_minerals <- c(self$temp_gated_minerals(porewater_temp_C),
                                                               self$stable_clathrate_phases(porewater_temp_C, ocean_pressure_atm, ch4_mode))
                                     extra_ocean_lines <- if (length(extra_ocean_minerals) > 0)
                                       paste0("   ", extra_ocean_minerals, " 0 0", collapse = "\n") else ""
                                     extra_eq_so_line <- paste(extra_ocean_minerals, collapse = " ")

                                     # ---- mode-specific mass computations ----
                                     if (mode == "hydrosphere") {
                                       porewater_mass_norm <- self$phreeqc_params$porewater_mass / self$phreeqc_params$layer_mass
                                       ocean_mass_norm     <- self$phreeqc_params$water_mass
                                       combined_mass_norm  <- porewater_mass_norm + ocean_mass_norm
                                     }
                                     .rock_kg    <- if (mode == "hydrosphere") self$phreeqc_params$total_rock else 1.0
                                     .water_mass <- if (mode == "porewater") water_rock_ratio else combined_mass_norm

                                     # Primary mineral moles. Coefficients (mol per unit .rock_kg): ROCK_MOL_PER_KG
                                     # (top of file). 9 kinetic primaries (incl. Schreibersite, placeholder
                                     # rate) + 3 equilibrium primaries (Fe, Ni, Lawrencite). Magnetite is no
                                     # longer primary — ordinary secondary phase only.
                                     mineral_scale_factor <- (100 - self$organic_wt_percent) / 100
                                     enstatite_moles  <- ROCK_MOL_PER_KG[["enstatite"]] * .rock_kg * mineral_scale_factor
                                     forsterite_moles <- ROCK_MOL_PER_KG[["forsterite"]] * .rock_kg * mineral_scale_factor
                                     troilite_moles   <- ROCK_MOL_PER_KG[["troilite"]] * .rock_kg * mineral_scale_factor
                                     albite_moles     <- ROCK_MOL_PER_KG[["albite"]] * .rock_kg * mineral_scale_factor
                                     diopside_moles   <- ROCK_MOL_PER_KG[["diopside"]] * .rock_kg * mineral_scale_factor
                                     anorthite_moles  <- ROCK_MOL_PER_KG[["anorthite"]] * .rock_kg * mineral_scale_factor
                                     kfeldspar_moles  <- ROCK_MOL_PER_KG[["kfeldspar"]] * .rock_kg * mineral_scale_factor
                                     tephroite_moles  <- ROCK_MOL_PER_KG[["tephroite"]] * .rock_kg * mineral_scale_factor
                                     fe_moles         <- ROCK_MOL_PER_KG[["fe"]] * .rock_kg * mineral_scale_factor
                                     ni_moles         <- ROCK_MOL_PER_KG[["ni"]] * .rock_kg * mineral_scale_factor
                                     lawrencite_moles <- ROCK_MOL_PER_KG[["lawrencite"]] * .rock_kg * mineral_scale_factor
                                     schreibersite_moles <- ROCK_MOL_PER_KG[["schreibersite"]] * .rock_kg * mineral_scale_factor

                                     orbital_period     <- self$simulator$config$orbital_period
                                     total_time_seconds <- n_orbital_cycles * orbital_period
                                     total_time_years   <- total_time_seconds / (365.25 * 24 * 3600)

                                     # Total_Initial_g = full initial primary rock mass: 8 kinetic minerals
                                     # (tracked via KIN() each cycle) plus the 3 equilibrium primaries.
                                     total_initial_expr <- sprintf(
                                       "%.6f*100.3725 + %.6f*140.6715 + %.6f*87.913 + %.6f*262.1798 + %.6f*216.55 + %.6f*278.164 + %.6f*278.33 + %.6f*201.96 + %.6f*55.845 + %.6f*58.693 + %.6f*126.75 + %.6f*198.509",
                                       enstatite_moles, forsterite_moles, troilite_moles, albite_moles,
                                       diopside_moles, anorthite_moles, kfeldspar_moles, tephroite_moles,
                                       fe_moles, ni_moles, lawrencite_moles, schreibersite_moles
                                     )

                                     # ---- mode-specific labels and header ----
                                     .pqi_title  <- if (mode == "porewater") "Enceladus porewater kinetic" else
                                                      "Enceladus bulk hydrosphere kinetic"
                                     .soln_label <- if (mode == "porewater") "Enceladus porewater" else
                                                      "Enceladus hydrosphere (ocean + porewater combined)"
                                     .header_lines <- if (mode == "porewater") {
                                       paste0(
                                         "# Mode: kinetic — porewater only vs permeable rock layer (no ocean)\n",
                                         "# water_rock_ratio: ", water_rock_ratio, "\n",
                                         "# n_orbital_cycles: ", n_orbital_cycles, "\n",
                                         "# Total simulation time: ", sprintf("%.4e", total_time_years),
                                         " years  (", sprintf("%.6e", total_time_seconds), " s)\n",
                                         "# Rock normalised to 1.0 kg; water = ", water_rock_ratio, " kg (porewater only)\n",
                                         "# Temperature: ", round(porewater_temp_C, 2), " C"
                                       )
                                     } else {
                                       paste0(
                                         "# thickness=", self$simulator$params$layer_thickness, " m\n",
                                         "# Mode: kinetic — full hydrosphere (ocean + porewater) vs permeable rock layer\n",
                                         "# n_orbital_cycles: ", n_orbital_cycles, "\n",
                                         "# Total simulation time: ", sprintf("%.4e", total_time_years),
                                         " years  (", sprintf("%.6e", total_time_seconds), " s)\n",
                                         "# Layer mass: ", sprintf("%.2e", self$phreeqc_params$layer_mass),
                                         " kg  (total_rock normalised to 1.0)\n",
                                         "# Porewater mass: ", sprintf("%.2e", self$phreeqc_params$porewater_mass), " kg\n",
                                         "# Ocean mass: ", sprintf("%.2e", self$phreeqc_params$ocean_mass), " kg\n",
                                         "# Combined hydrosphere (normalised): ",
                                         sprintf("%.6f", combined_mass_norm), "\n",
                                         "# Temperature: ", round(porewater_temp_C, 2), " C"
                                       )
                                     }

                                     base_name <- tools::file_path_sans_ext(basename(filename))

                                     content <- glue::glue('
TITLE {.pqi_title}
{ch4_db_note_line}
{clathrate_note_lines}

KNOBS
   -logfile true

{.header_lines}

SOLUTION 1 {.soln_label}
    temp      {round(porewater_temp_C, 2)}
    pH        7
    pe        4
    redox     pe
    units     mol/kgw
    density   1
    water     {.water_mass}
    pressure  {round(ocean_pressure_atm, 1)}

INCREMENTAL_REACTIONS true

RATES
Forsterite
\t-start
\t1   REM Ref PK04
\t10  kacid = 10^(-6.85) * exp(-67.2e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^0.47
\t20  kneut = 10^(-10.64) * exp(-79.0e3/8.314 * (1/TK-1/298.15))
\t21  SSA = 0.1
\t22  mw = 140.6715
\t40  k = (kacid + kneut) * SSA * mw * M
\t50  IF SR("Forsterite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Forsterite"))
\t60  moles = rate * TIME
\t70  SAVE moles
\t-end

Troilite # PK04 has no separate troilite entry; hexagonal pyrrhotite parameters used as proxy (same FeS chemistry)
\t-start
\t1   REM Acid mechanism only (PK04, hexagonal pyrrhotite formulation)
\t25  kacid = 10^(-6.79) * exp(-63.0e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^-0.090 * ACT("Fe+3")^0.356
\t26  SSA = 5
\t27  mw = 87.913
\t30  k = kacid * SSA * mw * M
\t40  IF SR("Troilite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Troilite"))
\t50  moles = rate * TIME
\t55  IF moles < 0 THEN moles = 0
\t60  SAVE moles
\t-end

Diopside
\t-start
\t1   REM Ref PK04
\t10  kacid = 10^(-6.36) * exp(-96.1e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^0.71
\t20  kneut = 10^(-11.11) * exp(-40.6e3/8.314 * (1/TK-1/298.15))
\t21  SSA = 0.1
\t22  mw = 216.55
\t40  k = (kacid + kneut) * SSA * mw * M
\t50  IF SR("Diopside") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Diopside"))
\t60  moles = rate * TIME
\t70  SAVE moles
\t-end

K-Feldspar
\t-start
\t1   REM Ref PK04
\t10  kacid = 10^(-10.06) * exp(-51.7e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^0.5
\t20  kneut = 10^(-12.41) * exp(-38.0e3/8.314 * (1/TK-1/298.15))
\t30  kbase = 10^(-21.20) * exp(-94.1e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^-0.823
\t31  SSA = 5
\t32  mw = 278.33
\t40  k = (kacid + kneut + kbase) * SSA * mw * M
\t50  IF SR("K-Feldspar") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("K-Feldspar"))
\t60  moles = rate * TIME
\t70  SAVE moles
\t-end

Tephroite # TODO: rate constants not yet supplied (Mn-olivine) — k held at 0 (no dissolution) until provided
\t-start
\t1   REM TODO awaiting kacid/kneut from user
\t10  kacid = 0
\t20  kneut = 0
\t21  SSA = 0.1
\t22  mw = 201.96
\t40  k = (kacid + kneut) * SSA * mw * M
\t50  IF SR("Tephroite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Tephroite"))
\t60  moles = rate * TIME
\t70  SAVE moles
\t-end

{SCHREIBERSITE_RATE}

Enstatite
\t-start
\t1   REM Ref PK04
\t10  kacid = 10^(-9.02) * exp(-80.0e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^0.6
\t20  kneut = 10^(-12.72) * exp(-80.0e3/8.314 * (1/TK-1/298.15))
\t21  SSA = 0.1
\t22  mw = 100.3725
\t40  k = (kacid + kneut) * SSA * mw * M
\t50  IF SR("Enstatite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Enstatite"))
\t60  moles = rate * TIME
\t70  SAVE moles
\t-end

Anorthite
\t-start
\t1   REM Ref PK04
\t10  kacid = 10^(-3.50) * exp(-16.6e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^1.411
\t20  kneut = 10^(-9.12) * exp(-17.8e3/8.314 * (1/TK-1/298.15))
\t21  SSA = 5
\t22  mw = 278.164
\t40  k = (kacid + kneut) * SSA * mw * M
\t50  IF SR("Anorthite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Anorthite"))
\t60  moles = rate * TIME
\t70  SAVE moles
\t-end

Albite
\t-start
\t1   REM Ref PK04
\t10  kacid = 10^(-10.16) * exp(-65.0e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^0.457
\t20  kneut = 10^(-12.56) * exp(-69.8e3/8.314 * (1/TK-1/298.15))
\t30  kbase = 10^(-15.60) * exp(-71.0e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^-0.572
\t31  SSA = 5
\t32  mw = 262.1798
\t40  k = (kacid + kneut + kbase) * SSA * mw * M
\t50  IF SR("Albite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Albite"))
\t60  moles = rate * TIME
\t70  SAVE moles
\t-end

KINETICS 1
   Enstatite
      -m0 {enstatite_moles}
      -step_divide 1
   Forsterite
      -m0 {forsterite_moles}
      -step_divide 1
   Troilite
      -m0 {troilite_moles}
      -step_divide 1
   Albite
      -m0 {albite_moles}
      -step_divide 1
   Diopside
      -m0 {diopside_moles}
      -step_divide 1
   Anorthite
      -m0 {anorthite_moles}
      -step_divide 1
   K-Feldspar
      -m0 {kfeldspar_moles}
      -step_divide 1
   Tephroite
      -m0 {tephroite_moles}
      -step_divide 1
{schreibersite_kinetics_entry(schreibersite_moles)}
   -steps {sprintf("%.6e", total_time_seconds)} in {n_orbital_cycles} steps
   -cvode true
   -bad_step_max 5000

EQUILIBRIUM_PHASES 1
# --- Primary minerals (equilibrium; dissolve_only so they can only deplete) ---
   Fe                   0  {fe_moles}          dissolve_only
   Ni                   0  {ni_moles}          dissolve_only
   Lawrencite           0  {lawrencite_moles}  dissolve_only
# --- Secondary minerals (precipitate if supersaturated) ---
   Analcime 0 0
   Anhydrite 0 0
   Aragonite 0 0
   Bassanite 0 0
   Beidellite-Ca 0 0
   Beidellite-Fe 0 0
   Beidellite-Mg 0 0
   Beidellite-Na 0 0
   Boehmite 0 0
   Brucite 0 0
   C2H4(g) 0 0
   C2H6(g) 0 0
   C3H8(g) 0 0
   Calcite 0 0
   CH4(g) 0 0
   Chalcedony 0 0
   Chamosite 0 0
   Chrysotile 0 0
   Citric_Acid 0 0
   Clinochlore-14A 0 0
   Clinochlore-7A 0 0
   Clinoptilolite-Ca 0 0
   Clinoptilolite-Na 0 0
   Clinozoisite 0 0
   CO(g) 0 0
   CO2(g) 0 0
   Cronstedtite-7A 0 0
   Daphnite-14A 0 0
   Daphnite-7A 0 0
   Dawsonite 0 0
   Diaspore 0 0
   Dolomite 0 0
   Ettringite 0 0
   Fe(OH)2 0 0
   Fe(OH)3 0 0
   Fe2(SO4)3 0 0
   FeSO4 0 0
   Gibbsite 0 0
   Glycine 0 0
   Goethite 0 0
   Greenalite 0 0
   Gypsum 0 0
   Gyrolite 0 0
   H2(g) 0 0
   H2S(g) 0 0
   Hematite 0 0
   Huntite 0 0
   Hydromagnesite 0 0
   Kaolinite 0 0
   Lizardite 0 0
   Magnetite 0 0
   Melanterite 0 0
   MgSO4 0 0
   Minnesotaite 0 0
   Mirabilite 0 0
   Montmor-Ca 0 0
   Montmor-Mg 0 0
   Montmor-Na 0 0
   N2(g) 0 0
   NH3(g) 0 0
   Nontronite-Ca 0 0
   Nontronite-Mg 0 0
   Nontronite-Na 0 0
   Portlandite 0 0
   Pyrite 0 0
   Saponite-Fe-Ca 0 0
   Saponite-Fe-Fe 0 0
   Saponite-Fe-Mg 0 0
   Saponite-Fe-Na 0 0
   Saponite-Mg-Ca 0 0
   Saponite-Mg-Fe 0 0
   Saponite-Mg-Mg 0 0
   Saponite-Mg-Na 0 0
   Sepiolite 0 0
   SiO2(am) 0 0
   SO2(g) 0 0
   Thenardite 0 0
   Tobermorite-11A 0 0
{extra_ocean_lines}

SELECTED_OUTPUT
    -file {base_name}.csv
    -reset false
    -state true
    -solution true
    -step true
    -ph true
    -pe true
    -alkalinity true
    -ionic_strength true
    -water true
    -totals              Al C Ca Fe H Mg N Na O P S Si
    -molalities          OH- H+ C2H4 C2H6 CH4 CO HCO3- CO3-2 CO2
                         CH3COO- HCOO- HCN Ca+2 Cl- Fe+2 FeOH+ Fe+3 H2
                         Mg+2 NH4+ NH3 NH4CO3- N2 NO2-
                         NO3- Na+ H2PO4- HPO4-2 PO4-3 HS- H2S SO3-2 HSO3-
                         SO2 SO4-2 SiO2 HSiO3-
    -activities          OH- H+ C2H4 C2H6 CH4 CO HCO3- CO3-2 CO2
                         CH3COO- HCOO- HCN Ca+2 Cl- Fe+2 FeOH+ Fe+3 H2
                         Mg+2 NH4+ NH3 NH4CO3- N2 NO2-
                         NO3- Na+ H2PO4- HPO4-2 PO4-3 HS- H2S SO3-2 HSO3-
                         SO2 SO4-2 SiO2 HSiO3-
    -kinetics            Enstatite Forsterite Troilite Albite Diopside Anorthite K-Feldspar Tephroite Schreibersite
    -saturation_indices  Enstatite Forsterite Troilite Albite Diopside Anorthite K-Feldspar Tephroite
    -equilibrium_phases  Fe Ni Lawrencite Analcime Anhydrite Aragonite Bassanite
                         Beidellite-Ca Beidellite-Fe Beidellite-Mg Beidellite-Na
                         Boehmite Brucite
                         C2H4(g) C2H6(g) C3H8(g) Calcite CH4(g) Chalcedony Chamosite Chrysotile
                         Citric_Acid Clinochlore-14A Clinochlore-7A Clinoptilolite-Ca Clinoptilolite-Na
                         Clinozoisite CO(g) CO2(g) Cronstedtite-7A
                         Daphnite-14A Daphnite-7A Dawsonite Diaspore Dolomite
                         Ettringite Fe(OH)2 Fe(OH)3 Fe2(SO4)3 FeSO4
                         Gibbsite Glycine Goethite Greenalite Gypsum Gyrolite
                         H2(g) H2S(g) Hematite Huntite Hydromagnesite Kaolinite
                         Lizardite Magnetite Melanterite MgSO4 Minnesotaite Mirabilite
                         Montmor-Ca Montmor-Mg Montmor-Na
                         N2(g) NH3(g) Nontronite-Ca Nontronite-Mg Nontronite-Na
                         Portlandite Pyrite
                         Saponite-Fe-Ca Saponite-Fe-Fe Saponite-Fe-Mg Saponite-Fe-Na
                         Saponite-Mg-Ca Saponite-Mg-Fe Saponite-Mg-Mg Saponite-Mg-Na
                         Sepiolite SiO2(am) SO2(g) Thenardite Tobermorite-11A
                         {extra_eq_so_line}
    -gases               CO2(g) H2(g) CH4(g) CO(g) H2O(g) H2S(g) N2(g)
                         NH3(g) NO(g) NO2(g) O2(g) SO2(g)

USER_PUNCH
   -headings Time_Years Enst_remain_g Forst_remain_g Troil_remain_g Alb_remain_g Diop_remain_g Anorth_remain_g Kfs_remain_g Teph_remain_g Schr_remain_g Total_Initial_g
   10 PUNCH TOTAL_TIME / (365.25 * 24 * 3600)
   20 PUNCH KIN("Enstatite") * 100.3725
   30 PUNCH KIN("Forsterite") * 140.6715
   40 PUNCH KIN("Troilite") * 87.913
   50 PUNCH KIN("Albite") * 262.1798
   60 PUNCH KIN("Diopside") * 216.55
   70 PUNCH KIN("Anorthite") * 278.164
   80 PUNCH KIN("K-Feldspar") * 278.33
   85 PUNCH KIN("Tephroite") * 201.96
   90 PUNCH KIN("Schreibersite") * 198.509
   100 total_initial = {total_initial_expr}
   110 PUNCH total_initial

END
')

                                     content <- self$apply_ch4_mode(content, ch4_mode)
                                     writeLines(content, filename)
                                     return(filename)
                                   },

                                   generate_phreeqc_input_bulk_kinetic_restart = function(
                                     output_dir          = "phreeqc_inputs",
                                     prev_csv            = NULL,
                                     n_orbital_cycles    = NULL,
                                     orbits_per_step     = 1,
                                     equil_minerals      = character(0),
                                     mode                = c("porewater", "hydrosphere"),
                                     porewater_temp      = "ocean",
                                     suffix              = "",
                                     CH4_redox_override  = NULL
                                   ) {
                                     mode <- match.arg(mode)
                                     if (is.null(prev_csv) || !file.exists(prev_csv))
                                       stop("prev_csv must be a valid path to an existing CSV file")
                                     if (is.null(n_orbital_cycles) || n_orbital_cycles <= 0)
                                       stop("n_orbital_cycles must be a positive integer > 0")
                                     if (orbits_per_step < 1 || n_orbital_cycles %% orbits_per_step != 0)
                                       stop("orbits_per_step must be >= 1 and must evenly divide n_orbital_cycles")

                                     dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

                                     # Read last valid row from previous CSV
                                     prev_data <- read.delim(prev_csv, sep = "\t", header = TRUE, strip.white = TRUE)
                                     prev_data <- prev_data[prev_data$step != -99, ]
                                     if (nrow(prev_data) == 0)
                                       stop("No valid (non-sentinel) rows found in: ", prev_csv)
                                     last_row <- tail(prev_data, 1)

                                     # Solution state from last row
                                     ph_val   <- as.numeric(last_row[["pH"]])
                                     pe_val   <- as.numeric(last_row[["pe"]])
                                     water_kg <- as.numeric(last_row[["mass_H2O"]])

                                     temp_tag <- if (is.numeric(porewater_temp)) sprintf("T%g", porewater_temp) else sprintf("T%s", porewater_temp)
                                     org_tag  <- sprintf("org%g", self$organic_wt_percent)

                                     filename <- if (mode == "porewater") {
                                       file.path(output_dir, sprintf(
                                         "enceladus_wr%.4g_%s_%s_orb%g_step%g_porewater_kinetic_restart%s.pqi",
                                         water_kg, temp_tag, org_tag, n_orbital_cycles, orbits_per_step, suffix
                                       ))
                                     } else {
                                       file.path(output_dir, sprintf(
                                         "enceladus_d%g_%s_%s_orb%g_step%g_hydrosphere_kinetic_restart%s.pqi",
                                         self$simulator$params$layer_thickness, temp_tag, org_tag, n_orbital_cycles,
                                         orbits_per_step, suffix
                                       ))
                                     }

                                     elements  <- c("Al", "C", "Ca", "Fe", "Mg", "N", "Na", "P", "S", "Si")
                                     elem_vals <- sapply(elements, function(el) {
                                       val <- suppressWarnings(as.numeric(last_row[[el]]))
                                       if (!is.na(val) && val > 0) val else NA_real_
                                     })
                                     elem_lines <- paste(
                                       Filter(Negate(is.null), lapply(seq_along(elements), function(i) {
                                         if (is.na(elem_vals[i])) return(NULL)
                                         sprintf("    %-4s  %.6e", elements[i], elem_vals[i])
                                       })),
                                       collapse = "\n"
                                     )

                                     # Temperature and pressure (ocean floor)
                                     ocean_temp_C <- 274.6 - 273.15
                                     if (!is.null(self$profile_file) && file.exists(self$profile_file)) {
                                       pore_temp_C        <- mean_layer_temp_C(self$profile_file, self$simulator$params$layer_thickness)
                                       pore_pressure_atm  <- mean_layer_pressure_MPa(self$profile_file, self$simulator$params$layer_thickness) / 0.101325
                                       ocean_pressure_atm <- ocean_floor_pressure_MPa(self$profile_file) / 0.101325
                                     } else {
                                       if (!is.null(self$profile_file))
                                         warning("Profile file not found: '", self$profile_file, "' - using default 70 atm.")
                                       pore_temp_C        <- 25
                                       pore_pressure_atm  <- 70
                                       ocean_pressure_atm <- 70
                                     }

                                     # Resolve porewater_temp to a concrete value — this is the reaction
                                     # temperature for the single combined solution used here.
                                     porewater_temp_C <- if (is.numeric(porewater_temp)) porewater_temp else
                                       if (porewater_temp == "profile") pore_temp_C else ocean_temp_C

                                     # CH4 (redox-coupled) vs Mtg (decoupled) — see resolve_ch4_mode().
                                     ch4_mode <- self$resolve_ch4_mode(porewater_temp_C, CH4_redox_override)
                                     ch4_override_note <- if (!is.null(CH4_redox_override)) sprintf(", override=%s", CH4_redox_override) else ""
                                     ch4_db_note_line <- sprintf(
                                       "# Recommended database: %s (CH4 %s; porewater T = %.1f C%s)",
                                       ch4_mode$database,
                                       if (ch4_mode$decoupled) "decoupled -> aqueous/gas species are Mtg/Mtg(g)" else "coupled -> aqueous/gas species are CH4/CH4(g)",
                                       porewater_temp_C, ch4_override_note
                                     )

                                     # Clathrate stability at porewater conditions — see resolve_clathrate_stability().
                                     # Single combined solution here uses ocean_pressure_atm as its reaction
                                     # pressure (pore_pressure_atm is computed above but not otherwise used
                                     # in this function), so the stability check matches that.
                                     clathrate_note_lines <- self$clathrate_stability_note(porewater_temp_C, ocean_pressure_atm)

                                     # Primary mineral classification. Fe/Ni/Lawrencite are handled
                                     # separately below (always_equil_primary) — always equilibrium,
                                     # recovered via get_equil_moles() (their native PHREEQC
                                     # -equilibrium_phases output column), never kinetic-switchable.
                                     all_primary      <- c("Enstatite", "Forsterite", "Troilite", "Albite",
                                                           "Diopside", "Anorthite", "K-Feldspar", "Tephroite", "Schreibersite")
                                     if ("Schreibersite" %in% equil_minerals)
                                       stop("Schreibersite cannot be an equilibrium primary: neither database has an Fe3P phase (no log K). It is kinetic only.")
                                     equil_primary    <- intersect(equil_minerals, all_primary)
                                     kinetic_minerals <- setdiff(all_primary, equil_primary)
                                     always_equil_primary <- c("Fe", "Ni", "Lawrencite")

                                     molar_masses <- c(Enstatite = 100.3725, Forsterite = 140.6715, Troilite = 87.913,
                                                       Albite = 262.1798, Diopside = 216.55, Anorthite = 278.164,
                                                       `K-Feldspar` = 278.33, Tephroite = 201.96,
                                                      Schreibersite = SCHREIBERSITE_MW)

                                     # USER_PUNCH column names for each primary mineral (grams)
                                     remain_g_cols <- c(
                                       Enstatite  = "Enst_remain_g",  Forsterite   = "Forst_remain_g",
                                       Troilite   = "Troil_remain_g", Albite       = "Alb_remain_g",
                                       Diopside   = "Diop_remain_g",  Anorthite    = "Anorth_remain_g",
                                       `K-Feldspar` = "Kfs_remain_g", Tephroite    = "Teph_remain_g",
                                       Schreibersite = "Schr_remain_g"
                                     )

                                     get_k_moles <- function(mineral) {
                                       col <- paste0("k_", mineral)
                                       val <- if (col %in% names(last_row)) suppressWarnings(as.numeric(last_row[[col]])) else NA_real_
                                       if (!is.na(val) && val >= 0) return(val)
                                       # Fallback: mineral was equilibrium in prev run — read *_remain_g column
                                       g_col <- remain_g_cols[mineral]
                                       if (!is.na(g_col) && g_col %in% names(last_row)) {
                                         g_val <- suppressWarnings(as.numeric(last_row[[g_col]]))
                                         if (!is.na(g_val) && g_val >= 0) return(g_val / molar_masses[mineral])
                                       }
                                       0
                                     }

                                     # ---- PATCH D2: recover IOM channels from prev_csv for restart ----
                                     iom_rates_restart          <- ""
                                     iom_kinetics_restart       <- ""
                                     iom_kin_so                 <- ""
                                     iom_punch_headings_restart <- ""
                                     iom_punch_lines_restart    <- ""

                                     if (self$organic_wt_percent > 0) {
                                       # See the main-path comment (iom_default_config()) for why this is no
                                       # longer a hardcoded table -- WIRING CHANGE 2026-09-26.
                                       if (is.null(self$iom_pools)) {
                                         self$iom_pools <- iom_default_config()
                                       }

                                       # Guard: IOM_*_mol columns must exist — written by PATCH D1.
                                       # Silently restarting with zero organics is the failure mode most
                                       # likely to survive undetected into a published run.
                                       iom_mol_cols <- paste0(self$iom_pools$name, "_mol")
                                       missing_iom  <- iom_mol_cols[!iom_mol_cols %in% names(last_row)]
                                       if (length(missing_iom) > 0)
                                         stop(paste0(
                                           "Restart error: organic_wt_percent > 0 but IOM columns missing from prev_csv:\n  ",
                                           paste(missing_iom, collapse = ", "), "\n",
                                           "The CSV must come from a run with PATCH D1 applied ",
                                           "(IOM channels punched to output)."
                                         ))

                                       get_iom_moles <- function(channel) {
                                         val <- suppressWarnings(as.numeric(last_row[[paste0(channel, "_mol")]]))
                                         if (!is.na(val) && val >= 0) val else 0
                                       }

                                       iom_rates_restart <- paste0(apply(self$iom_pools, 1, function(p) {
                                         sprintf(paste0(
                                           "%s\n -start\n",
                                           " 1  REM Arrhenius 1st-order kerogen degradation (single-Ea pool)\n",
                                           " 10 k = 10^(%s) * exp(-%s / (8.314 * TK))\n",
                                           " 20 rate = k * M\n",
                                           " 30 if (M <= 0) then rate = 0\n",
                                           " 40 moles = rate * TIME\n",
                                           " 50 SAVE moles\n -end\n"),
                                           p[["name"]], p[["logA"]], p[["Ea_J"]])
                                       }), collapse = "")

                                       iom_kinetics_restart <- paste0(apply(self$iom_pools, 1, function(p) {
                                         m0 <- get_iom_moles(p[["name"]])
                                         sprintf("   %s\n      -formula %s\n      -m0 %.6e\n      -step_divide 1",
                                                 p[["name"]], trimws(p[["formula"]]), m0)
                                       }), collapse = "\n")

                                       iom_kin_so <- paste(self$iom_pools$name, collapse = " ")

                                       iom_punch_headings_restart <- paste0(" ",
                                         paste(sprintf("%s_mol", self$iom_pools$name), collapse = " "))
                                       iom_punch_lines_restart <- paste0(
                                         vapply(seq_len(nrow(self$iom_pools)), function(j) {
                                           sprintf("   %d PUNCH KIN(\"%s\")\n",
                                                   200 + 10 * j, self$iom_pools$name[j])
                                         }, character(1)), collapse = "")
                                     }

                                     orbital_period     <- self$simulator$config$orbital_period
                                     total_time_seconds <- n_orbital_cycles * orbital_period
                                     total_time_years   <- total_time_seconds / (365.25 * 24 * 3600)
                                     n_steps            <- n_orbital_cycles / orbits_per_step

                                     secondary_minerals <- c(
                                       "Analcime", "Anhydrite", "Aragonite", "Bassanite",
                                       "Beidellite-Ca", "Beidellite-Fe", "Beidellite-Mg", "Beidellite-Na",
                                       "Boehmite", "Brucite",
                                       "C2H4(g)", "C2H6(g)", "C3H8(g)", "Calcite", "CH4(g)",
                                       "Chalcedony", "Chamosite", "Chrysotile", "Citric_Acid",
                                       "Clinochlore-14A", "Clinochlore-7A", "Clinoptilolite-Ca", "Clinoptilolite-Na",
                                       "Clinozoisite", "CO(g)", "CO2(g)", "Cronstedtite-7A",
                                       "Daphnite-14A", "Daphnite-7A", "Dawsonite", "Diaspore", "Dolomite",
                                       "Ettringite", "Fe(OH)2", "Fe(OH)3", "Fe2(SO4)3", "FeSO4",
                                       "Gibbsite", "Glycine", "Goethite", "Greenalite", "Gypsum", "Gyrolite",
                                       "H2(g)", "H2S(g)", "Hematite", "Huntite", "Hydromagnesite",
                                       "Kaolinite", "Lizardite", "Melanterite", "MgSO4", "Minnesotaite", "Mirabilite",
                                       "Montmor-Ca", "Montmor-Mg", "Montmor-Na",
                                       "N2(g)", "NH3(g)", "Nontronite-Ca", "Nontronite-Mg", "Nontronite-Na",
                                       "Portlandite", "Pyrite",
                                       "Saponite-Fe-Ca", "Saponite-Fe-Fe", "Saponite-Fe-Mg", "Saponite-Fe-Na",
                                       "Saponite-Mg-Ca", "Saponite-Mg-Fe", "Saponite-Mg-Mg", "Saponite-Mg-Na",
                                       "Sepiolite", "SiO2(am)", "SO2(g)", "Thenardite", "Tobermorite-11A"
                                     )
                                     # Single solution (porewater_temp_C) gates the whole assemblage here.
                                     # Clathrate phases stable against ocean_pressure_atm (this function's own
                                     # reaction pressure) are appended the same way.
                                     secondary_minerals <- c(secondary_minerals, self$temp_gated_minerals(porewater_temp_C),
                                                             self$stable_clathrate_phases(porewater_temp_C, ocean_pressure_atm, ch4_mode))

                                     get_equil_moles <- function(mineral) {
                                       col <- make.names(mineral)
                                       # "Fe" collides with the -totals Fe (dissolved iron) column — PHREEQC
                                       # writes both under the literal header "Fe", and read.delim's default
                                       # check.names=TRUE renames the second occurrence to "Fe.1". Since
                                       # -equilibrium_phases always comes after -totals in these templates,
                                       # take the LAST matching column so the mineral (not the solution
                                       # total) is what's recovered.
                                       matches <- names(last_row)[names(last_row) == col |
                                                                   grepl(paste0("^", col, "\\.[0-9]+$"), names(last_row))]
                                       val <- if (length(matches) > 0) suppressWarnings(as.numeric(last_row[[tail(matches, 1)]])) else NA_real_
                                       if (is.na(val) || val < 0) 0 else val
                                     }

                                     # ---- RATES block ----
                                     .rate_laws <- list(
                                       Schreibersite = SCHREIBERSITE_RATE_BODY,
                                       Forsterite = paste(
                                         '\t-start',
                                         '\t1   REM Ref PK04',
                                         '\t10  kacid = 10^(-6.85) * exp(-67.2e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^0.47',
                                         '\t20  kneut = 10^(-10.64) * exp(-79.0e3/8.314 * (1/TK-1/298.15))',
                                         '\t21  SSA = 0.1',
                                         '\t22  mw = 140.6715',
                                         '\t40  k = (kacid + kneut) * SSA * mw * M',
                                         '\t50  IF SR("Forsterite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Forsterite"))',
                                         '\t60  moles = rate * TIME',
                                         '\t70  SAVE moles',
                                         '\t-end', sep = "\n"),
                                       Troilite = paste(
                                         '\t-start',
                                         '\t1   REM PK04 has no separate troilite entry; hexagonal pyrrhotite parameters used as proxy',
                                         '\t25  kacid = 10^(-6.79) * exp(-63.0e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^-0.090 * ACT("Fe+3")^0.356',
                                         '\t26  SSA = 5',
                                         '\t27  mw = 87.913',
                                         '\t30  k = kacid * SSA * mw * M',
                                         '\t40  IF SR("Troilite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Troilite"))',
                                         '\t50  moles = rate * TIME',
                                         '\t55  IF moles < 0 THEN moles = 0',
                                         '\t60  SAVE moles',
                                         '\t-end', sep = "\n"),
                                       Diopside = paste(
                                         '\t-start',
                                         '\t1   REM Ref PK04',
                                         '\t10  kacid = 10^(-6.36) * exp(-96.1e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^0.71',
                                         '\t20  kneut = 10^(-11.11) * exp(-40.6e3/8.314 * (1/TK-1/298.15))',
                                         '\t21  SSA = 0.1',
                                         '\t22  mw = 216.55',
                                         '\t40  k = (kacid + kneut) * SSA * mw * M',
                                         '\t50  IF SR("Diopside") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Diopside"))',
                                         '\t60  moles = rate * TIME',
                                         '\t70  SAVE moles',
                                         '\t-end', sep = "\n"),
                                       `K-Feldspar` = paste(
                                         '\t-start',
                                         '\t1   REM Ref PK04',
                                         '\t10  kacid = 10^(-10.06) * exp(-51.7e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^0.5',
                                         '\t20  kneut = 10^(-12.41) * exp(-38.0e3/8.314 * (1/TK-1/298.15))',
                                         '\t30  kbase = 10^(-21.20) * exp(-94.1e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^-0.823',
                                         '\t31  SSA = 5',
                                         '\t32  mw = 278.33',
                                         '\t40  k = (kacid + kneut + kbase) * SSA * mw * M',
                                         '\t50  IF SR("K-Feldspar") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("K-Feldspar"))',
                                         '\t60  moles = rate * TIME',
                                         '\t70  SAVE moles',
                                         '\t-end', sep = "\n"),
                                       Tephroite = paste(
                                         '\t-start',
                                         '\t1   REM TODO awaiting kacid/kneut from user (Mn-olivine, no rate constants supplied yet)',
                                         '\t10  kacid = 0',
                                         '\t20  kneut = 0',
                                         '\t21  SSA = 0.1',
                                         '\t22  mw = 201.96',
                                         '\t40  k = (kacid + kneut) * SSA * mw * M',
                                         '\t50  IF SR("Tephroite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Tephroite"))',
                                         '\t60  moles = rate * TIME',
                                         '\t70  SAVE moles',
                                         '\t-end', sep = "\n"),
                                       Enstatite = paste(
                                         '\t-start',
                                         '\t1   REM Ref PK04',
                                         '\t10  kacid = 10^(-9.02) * exp(-80.0e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^0.6',
                                         '\t20  kneut = 10^(-12.72) * exp(-80.0e3/8.314 * (1/TK-1/298.15))',
                                         '\t21  SSA = 0.1',
                                         '\t22  mw = 100.3725',
                                         '\t40  k = (kacid + kneut) * SSA * mw * M',
                                         '\t50  IF SR("Enstatite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Enstatite"))',
                                         '\t60  moles = rate * TIME',
                                         '\t70  SAVE moles',
                                         '\t-end', sep = "\n"),
                                       Anorthite = paste(
                                         '\t-start',
                                         '\t1   REM Ref PK04',
                                         '\t10  kacid = 10^(-3.50) * exp(-16.6e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^1.411',
                                         '\t20  kneut = 10^(-9.12) * exp(-17.8e3/8.314 * (1/TK-1/298.15))',
                                         '\t21  SSA = 5',
                                         '\t22  mw = 278.164',
                                         '\t40  k = (kacid + kneut) * SSA * mw * M',
                                         '\t50  IF SR("Anorthite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Anorthite"))',
                                         '\t60  moles = rate * TIME',
                                         '\t70  SAVE moles',
                                         '\t-end', sep = "\n"),
                                       Albite = paste(
                                         '\t-start',
                                         '\t1   REM Ref PK04',
                                         '\t10  kacid = 10^(-10.16) * exp(-65.0e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^0.457',
                                         '\t20  kneut = 10^(-12.56) * exp(-69.8e3/8.314 * (1/TK-1/298.15))',
                                         '\t30  kbase = 10^(-15.60) * exp(-71.0e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^-0.572',
                                         '\t31  SSA = 5',
                                         '\t32  mw = 262.1798',
                                         '\t40  k = (kacid + kneut + kbase) * SSA * mw * M',
                                         '\t50  IF SR("Albite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Albite"))',
                                         '\t60  moles = rate * TIME',
                                         '\t70  SAVE moles',
                                         '\t-end', sep = "\n")
                                     )

                                     .primary_rate_defs <- if (length(kinetic_minerals) > 0)
                                       paste(sapply(kinetic_minerals, function(m)
                                         paste0(m, "\n", .rate_laws[[m]])), collapse = "\n\n")
                                     else ""
                                     .all_rate_defs <- paste(
                                       Filter(nzchar, c(.primary_rate_defs, iom_rates_restart)),
                                       collapse = "\n\n")
                                     rates_block <- if (nzchar(.all_rate_defs))
                                       paste0("RATES\n", .all_rate_defs) else ""

                                     # ---- KINETICS block ----
                                     .primary_kin_entries <- if (length(kinetic_minerals) > 0)
                                       paste(sapply(kinetic_minerals, function(m)
                                         if (m == "Schreibersite") schreibersite_kinetics_entry(get_k_moles(m)) else sprintf("   %s\n      -m0 %.6e\n      -step_divide 1", m, get_k_moles(m))),
                                         collapse = "\n")
                                     else ""
                                     .all_kin_entries <- paste(
                                       Filter(nzchar, c(.primary_kin_entries, iom_kinetics_restart)),
                                       collapse = "\n")
                                     kinetics_block <- if (nzchar(.all_kin_entries))
                                       paste0("KINETICS 1\n", .all_kin_entries, "\n",
                                              sprintf("   -steps %.6e in %d steps\n", total_time_seconds, n_steps),
                                              "   -cvode true\n",
                                              "   -bad_step_max 5000")
                                     else ""

                                     # ---- EQUILIBRIUM_PHASES block ----
                                     # Fe/Ni/Lawrencite are always equilibrium (never in all_primary, so
                                     # never kinetic-switchable) — recovered via get_equil_moles() since
                                     # they're ordinary EQUILIBRIUM_PHASES entries whose remaining moles
                                     # PHREEQC writes natively (same mechanism as secondary_minerals).
                                     always_equil_prim_lines <- paste(
                                       sapply(always_equil_primary, function(m)
                                         sprintf("   %-20s 0 %.6e  dissolve_only", m, get_equil_moles(m))),
                                       collapse = "\n"
                                     )
                                     equil_prim_lines <- paste(
                                       c("# --- Primary minerals (equilibrium) ---",
                                         always_equil_prim_lines,
                                         if (length(equil_primary) > 0)
                                           sapply(equil_primary, function(m)
                                             sprintf("   %-20s 0 %.6e  dissolve_only", m, get_k_moles(m)))
                                         else character(0)),
                                       collapse = "\n")

                                     secondary_eq_lines <- paste(
                                       c("# --- Secondary minerals ---",
                                         sapply(secondary_minerals, function(m)
                                           sprintf("   %-20s 0 %.6e", m, get_equil_moles(m)))),
                                       collapse = "\n"
                                     )

                                     equil_phases_block <- paste(
                                       Filter(nzchar, c(equil_prim_lines, secondary_eq_lines)),
                                       collapse = "\n"
                                     )

                                     # ---- SELECTED_OUTPUT dynamic lines ----
                                     .all_kin_so <- c(kinetic_minerals,
                                                       if (nzchar(iom_kin_so)) strsplit(iom_kin_so, " ")[[1]]
                                                       else character(0))
                                     kin_so_line <- if (length(.all_kin_so) > 0)
                                       paste("    -kinetics           ", paste(.all_kin_so, collapse = " "))
                                     else ""

                                     .all_eq_so  <- c(always_equil_primary, equil_primary, secondary_minerals)
                                     .eq_chunks  <- split(.all_eq_so, ceiling(seq_along(.all_eq_so) / 7))
                                     equil_so_line <- paste(
                                       c(paste0("    -equilibrium_phases  ", paste(.eq_chunks[[1]], collapse = " ")),
                                         if (length(.eq_chunks) > 1)
                                           sapply(.eq_chunks[-1], function(ch)
                                             paste0("                         ", paste(ch, collapse = " ")))),
                                       collapse = "\n"
                                     )

                                     # ---- USER_PUNCH ----
                                     total_initial_expr <- paste(
                                       sapply(all_primary, function(m)
                                         sprintf("%.6f*%.4f", get_k_moles(m), molar_masses[[m]])),
                                       collapse = " + "
                                     )

                                     punch_lines <- paste(
                                       sapply(seq_along(all_primary), function(i) {
                                         m  <- all_primary[i]
                                         mm <- molar_masses[[m]]
                                         fn <- if (m %in% equil_primary) sprintf('EQUI("%s")', m) else sprintf('KIN("%s")', m)
                                         sprintf("   %d PUNCH %s * %.4f", 10 + i * 10, fn, mm)
                                       }),
                                       collapse = "\n"
                                     )

                                     # ---- mode-specific labels ----
                                     .pqi_title    <- if (mode == "porewater") "Enceladus porewater kinetic (restart)" else
                                                        "Enceladus bulk hydrosphere kinetic (restart)"
                                     .soln_label   <- if (mode == "porewater") "Enceladus porewater (restarted)" else
                                                        "Enceladus hydrosphere (restarted)"
                                     .thick_line   <- if (mode == "hydrosphere")
                                       sprintf("# thickness=%g m", self$simulator$params$layer_thickness) else ""
                                     .mode_comment <- if (mode == "porewater")
                                       "porewater kinetic restart from previous CSV state" else
                                       "kinetic restart from previous CSV state"

                                     base_name <- tools::file_path_sans_ext(basename(filename))

                                     content <- glue::glue('
TITLE {.pqi_title}
{ch4_db_note_line}
{clathrate_note_lines}

KNOBS
   -logfile true

{.thick_line}
# Mode: {.mode_comment}
# prev_csv: {prev_csv}
# equil_minerals: {paste(equil_primary, collapse = ", ")}
# n_orbital_cycles: {n_orbital_cycles}
# orbits_per_step: {orbits_per_step}  (step duration = {sprintf("%.4e", orbits_per_step * orbital_period)} s)
# n_steps: {n_steps}
# Total simulation time: {sprintf("%.4e", total_time_years)} years  ({sprintf("%.6e", total_time_seconds)} s)

SOLUTION 1 {.soln_label}
    temp      {round(porewater_temp_C, 2)}
    pH        {sprintf("%.4f", ph_val)}
    pe        {sprintf("%.4f", pe_val)}
    redox     pe
    units     mol/kgw
    density   1
    water     {sprintf("%.6f", water_kg)}
    pressure  {round(ocean_pressure_atm, 1)}
{elem_lines}

INCREMENTAL_REACTIONS true

{rates_block}

{kinetics_block}

EQUILIBRIUM_PHASES 1
{equil_phases_block}

SELECTED_OUTPUT
    -file {base_name}.csv
    -reset false
    -state true
    -solution true
    -step true
    -ph true
    -pe true
    -alkalinity true
    -ionic_strength true
    -water true
    -totals              Al C Ca Fe H Mg N Na O P S Si
    -molalities          OH- H+ C2H4 C2H6 CH4 CO HCO3- CO3-2 CO2
                         CH3COO- HCOO- HCN Ca+2 Cl- Fe+2 FeOH+ Fe+3 H2
                         Mg+2 NH4+ NH3 NH4CO3- N2 NO2-
                         NO3- Na+ H2PO4- HPO4-2 PO4-3 HS- H2S SO3-2 HSO3-
                         SO2 SO4-2 SiO2 HSiO3-
    -activities          OH- H+ C2H4 C2H6 CH4 CO HCO3- CO3-2 CO2
                         CH3COO- HCOO- HCN Ca+2 Cl- Fe+2 FeOH+ Fe+3 H2
                         Mg+2 NH4+ NH3 NH4CO3- N2 NO2-
                         NO3- Na+ H2PO4- HPO4-2 PO4-3 HS- H2S SO3-2 HSO3-
                         SO2 SO4-2 SiO2 HSiO3-
{kin_so_line}
{equil_so_line}
    -saturation_indices  {paste(setdiff(all_primary, "Schreibersite"), collapse = " ")}
    -gases               CO2(g) H2(g) CH4(g) CO(g) H2O(g) H2S(g) N2(g)
                         NH3(g) NO(g) NO2(g) O2(g) SO2(g)

USER_PUNCH
   -headings Time_Years Enst_remain_g Forst_remain_g Troil_remain_g Alb_remain_g Diop_remain_g Anorth_remain_g Kfs_remain_g Teph_remain_g Schr_remain_g Total_Initial_g{iom_punch_headings_restart}
   10 PUNCH TOTAL_TIME / (365.25 * 24 * 3600)
{punch_lines}
   {10 + (length(all_primary) + 1) * 10} total_initial = {total_initial_expr}
   {10 + (length(all_primary) + 2) * 10} PUNCH total_initial
{iom_punch_lines_restart}
END
')

                                     content <- self$apply_ch4_mode(content, ch4_mode)
                                     writeLines(content, filename)
                                     return(filename)
                                   },

                                   # ---------------------------------------------------------------
                                   # Tidal kinetic restart: continues a tidal-exchange kinetic run
                                   # from the end of a previous run.  Reads solution 1 (porewater),
                                   # solution 2 (ocean), kinetics, and both EP slots from the CSV.
                                   # Primary minerals listed in equil_minerals are moved from KINETICS
                                   # to EQUILIBRIUM_PHASES 1 (same as hydrosphere restart).
                                   # ---------------------------------------------------------------
                                   generate_phreeqc_input_tidal_kinetic_restart = function(
                                     output_dir          = "phreeqc_inputs",
                                     prev_csv            = NULL,
                                     n_orbital_cycles    = NULL,
                                     substeps            = 1,
                                     equil_minerals      = character(0),
                                     mode                = "ungrouped",
                                     grouping_factor     = NULL,
                                     exchange_ratio      = NULL,
                                     porewater_temp      = "ocean",
                                     depletion_mineral   = NULL,
                                     suffix              = "",
                                     CH4_redox_override  = NULL
                                   ) {
                                     if (is.null(prev_csv) || !file.exists(prev_csv))
                                       stop("prev_csv must be a valid path to an existing CSV file")
                                     if (is.null(n_orbital_cycles) || n_orbital_cycles <= 0)
                                       stop("n_orbital_cycles must be a positive integer")

                                     dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

                                     sub_tag  <- if (substeps != 1) sprintf("_sub%d", substeps) else ""
                                     temp_tag <- if (is.numeric(porewater_temp)) sprintf("T%g", porewater_temp) else sprintf("T%s", porewater_temp)
                                     org_tag  <- sprintf("org%g", self$organic_wt_percent)
                                     filename <- file.path(output_dir, if (mode == "grouped" && !is.null(grouping_factor))
                                       sprintf("enceladus_k%.2e_d%g_%s_%s_grp%g_orb%g%s_tidal_kinetic_restart%s.pqi",
                                               self$simulator$params$k,
                                               self$simulator$params$layer_thickness,
                                               temp_tag, org_tag,
                                               grouping_factor, n_orbital_cycles, sub_tag, suffix)
                                     else
                                       sprintf("enceladus_k%.2e_d%g_%s_%s_orb%g%s_tidal_kinetic_restart%s.pqi",
                                               self$simulator$params$k,
                                               self$simulator$params$layer_thickness,
                                               temp_tag, org_tag,
                                               n_orbital_cycles, sub_tag, suffix)
                                     )

                                     # ---- Read CSV and locate last complete tidal cycle ----
                                     prev_data <- read.delim(prev_csv, sep = "\t", header = TRUE,
                                                             strip.white = TRUE, check.names = FALSE)
                                     prev_data <- prev_data[trimws(prev_data[["step"]]) != "-99", ]
                                     if (nrow(prev_data) == 0)
                                       stop("No valid (non-sentinel) rows found in: ", prev_csv)

                                     # Anchor on last soln=2 row (end of a complete MIX2 step)
                                     soln_col  <- trimws(prev_data[["soln"]])
                                     mix2_rows <- which(soln_col == "2")
                                     if (length(mix2_rows) == 0)
                                       stop("No soln=2 rows found — cannot identify a complete tidal cycle in: ", prev_csv)

                                     idx_mix2  <- tail(mix2_rows, 1)
                                     idx_mix1  <- idx_mix2 - 1
                                     idx_mix4  <- idx_mix2 - 2
                                     idx_react <- idx_mix2 - 3

                                     if (idx_react < 1)
                                       stop("Not enough rows before last soln=2 to reconstruct full cycle state")

                                     row_mix2  <- prev_data[idx_mix2, ]   # ocean chemistry  (MIX2)
                                     row_mix1  <- prev_data[idx_mix1, ]   # porewater chem   (MIX1)
                                     row_mix4  <- prev_data[idx_mix4, ]   # EP2 state        (MIX4)
                                     row_react <- prev_data[idx_react, ]  # kinetics + EP1   (REACT)

                                     # Sanity-check soln fields
                                     expected <- c(idx_mix1 = "1", idx_mix4 = "4", idx_react = "1")
                                     for (nm in names(expected)) {
                                       obs <- trimws(prev_data[get(nm), "soln"])
                                       if (obs != expected[[nm]])
                                         warning(sprintf("Expected soln=%s at %s (row %d), got %s — CSV structure may differ",
                                                         expected[[nm]], nm, get(nm), obs))
                                     }

                                     # Cumulative time offset from previous run
                                     restart_time_years <- suppressWarnings(as.numeric(row_mix2[["Time_Years"]]))
                                     if (is.na(restart_time_years))
                                       restart_time_years <- suppressWarnings(as.numeric(row_react[["Time_Years"]]))
                                     if (is.na(restart_time_years)) restart_time_years <- 0

                                     # ---- Orbital / cycle parameters (same as non-restart template) ----
                                     metrics        <- self$simulator$get_summary_metrics()
                                     orbital_period <- self$simulator$config$orbital_period

                                     gf <- if (!is.null(grouping_factor)) grouping_factor else 1
                                     if (mode == "ocean_circulations") {
                                       adjusted_cycles   <- ceiling(
                                         (ceiling((self$phreeqc_params$total_rock / self$phreeqc_params$rock_per_year *
                                                     365.25 * 24 * 3600) / orbital_period)) /
                                           self$phreeqc_params$grouping_factor)
                                       kinetic_cycles    <- ceiling(adjusted_cycles * n_orbital_cycles)
                                       time_step_seconds <- (metrics$circulation_time_years / adjusted_cycles) *
                                         365.25 * 24 * 3600
                                     } else if (mode == "grouped") {
                                       kinetic_cycles    <- ceiling(n_orbital_cycles / gf)
                                       time_step_seconds <- gf * orbital_period
                                     } else {
                                       kinetic_cycles    <- n_orbital_cycles
                                       time_step_seconds <- orbital_period
                                     }
                                     time_step_years <- time_step_seconds / (365.25 * 24 * 3600)

                                     # Mixing fractions
                                     porewater_mass  <- self$phreeqc_params$porewater_mass
                                     ocean_mass      <- self$phreeqc_params$ocean_mass
                                     m_total         <- porewater_mass + ocean_mass
                                     f_pore          <- porewater_mass / m_total
                                     f_ocean         <- ocean_mass     / m_total
                                     G_cycles        <- time_step_seconds / orbital_period
                                     x_tide          <- min(metrics$exchanging_fluid_mass / porewater_mass, 1)
                                     x_step          <- 1 - (1 - x_tide)^G_cycles
                                     f_x             <- x_step * f_pore

                                     porewater_mass_norm <- porewater_mass / self$phreeqc_params$layer_mass
                                     ocean_mass_norm     <- self$phreeqc_params$water_mass

                                     # ---- Temperature / pressure ----
                                     ocean_temp_C <- 274.6 - 273.15
                                     if (!is.null(self$profile_file) && file.exists(self$profile_file)) {
                                       pore_temp_C        <- mean_layer_temp_C(self$profile_file, self$simulator$params$layer_thickness)
                                       pore_pressure_atm  <- mean_layer_pressure_MPa(self$profile_file, self$simulator$params$layer_thickness) / 0.101325
                                       ocean_pressure_atm <- ocean_floor_pressure_MPa(self$profile_file) / 0.101325
                                     } else {
                                       pore_temp_C        <- 25
                                       pore_pressure_atm  <- 70
                                       ocean_pressure_atm <- 70
                                     }
                                     porewater_temp_C <- if (is.numeric(porewater_temp)) porewater_temp else
                                       if (porewater_temp == "profile") pore_temp_C else ocean_temp_C

                                     # CH4 (redox-coupled) vs Mtg (decoupled) — see resolve_ch4_mode().
                                     ch4_mode <- self$resolve_ch4_mode(porewater_temp_C, CH4_redox_override)
                                     ch4_override_note <- if (!is.null(CH4_redox_override)) sprintf(", override=%s", CH4_redox_override) else ""
                                     ch4_db_note_line <- sprintf(
                                       "# Recommended database: %s (CH4 %s; porewater T = %.1f C%s)",
                                       ch4_mode$database,
                                       if (ch4_mode$decoupled) "decoupled -> aqueous/gas species are Mtg/Mtg(g)" else "coupled -> aqueous/gas species are CH4/CH4(g)",
                                       porewater_temp_C, ch4_override_note
                                     )

                                     # Clathrate stability — separate checks for the ocean and porewater
                                     # solutions, since this function gives them distinct temperature/pressure.
                                     # See resolve_clathrate_stability()/clathrate_stability_note().
                                     clathrate_note_lines <- paste0(
                                       self$clathrate_stability_note(ocean_temp_C, ocean_pressure_atm, location = "Ocean"),
                                       "\n\n",
                                       self$clathrate_stability_note(porewater_temp_C, pore_pressure_atm, location = "Porewater")
                                     )

                                     # ---- Primary mineral classification ----
                                     # Fe/Ni/Lawrencite are handled separately (always_equil_primary, below) —
                                     # always equilibrium, never kinetic-switchable, porewater (EP1) only.
                                     all_primary      <- c("Enstatite", "Forsterite", "Troilite", "Albite",
                                                           "Diopside", "Anorthite", "K-Feldspar", "Tephroite", "Schreibersite")
                                     if ("Schreibersite" %in% equil_minerals)
                                       stop("Schreibersite cannot be an equilibrium primary: neither database has an Fe3P phase (no log K). It is kinetic only.")
                                     equil_primary    <- intersect(equil_minerals, all_primary)
                                     kinetic_minerals <- setdiff(all_primary, equil_primary)
                                     always_equil_primary <- c("Fe", "Ni", "Lawrencite")

                                     molar_masses <- c(Enstatite = 100.3725, Forsterite = 140.6715, Troilite = 87.913,
                                                       Albite = 262.1798, Diopside = 216.55, Anorthite = 278.164,
                                                       `K-Feldspar` = 278.33, Tephroite = 201.96,
                                                      Schreibersite = SCHREIBERSITE_MW)

                                     # Helper: read kinetic moles from REACT row
                                     get_k_moles <- function(mineral) {
                                       col <- paste0("k_", mineral)
                                       val <- if (col %in% names(row_react)) suppressWarnings(as.numeric(row_react[[col]])) else NA_real_
                                       if (!is.na(val) && val >= 0) return(val)
                                       # Fallback: use _remain_g column
                                       g_cols <- c(Enstatite="Enst_remain_g", Forsterite="Forst_remain_g",
                                                   Troilite="Troil_remain_g", Albite="Alb_remain_g",
                                                   Diopside="Diop_remain_g", Anorthite="Anorth_remain_g",
                                                   `K-Feldspar`="Kfs_remain_g", Tephroite="Teph_remain_g", Schreibersite="Schr_remain_g")
                                       g_col <- g_cols[mineral]
                                       if (!is.na(g_col) && g_col %in% names(row_react)) {
                                         g_val <- suppressWarnings(as.numeric(row_react[[g_col]]))
                                         if (!is.na(g_val) && g_val >= 0) return(g_val / molar_masses[mineral])
                                       }
                                       0
                                     }

                                     # Helper: read equilibrium phase moles from a given row. "Fe" collides
                                     # with the -totals Fe (dissolved iron) column — PHREEQC writes both
                                     # under the literal header "Fe", so with check.names = FALSE both stay
                                     # named "Fe" (no automatic .1 suffix); take the LAST match, since
                                     # -equilibrium_phases always comes after -totals in these templates.
                                     get_ep_moles <- function(mineral, row) {
                                       col <- make.names(mineral)
                                       matches <- which(make.names(names(row)) == col)
                                       val <- if (length(matches) > 0)
                                         suppressWarnings(as.numeric(row[[tail(matches, 1)]])) else NA_real_
                                       if (is.na(val) || val < 0) 0 else val
                                     }

                                     # Helper: read element totals from solution row
                                     get_chem <- function(row, elements) {
                                       sapply(elements, function(el) {
                                         val <- suppressWarnings(as.numeric(row[[el]]))
                                         if (!is.na(val) && val > 0) val else NA_real_
                                       })
                                     }

                                     # ---- Solution chemistry ----
                                     elements <- c("Al", "C", "Ca", "Fe", "Mg", "N", "Na", "P", "S", "Si")

                                     make_soln_lines <- function(row) {
                                       vals <- get_chem(row, elements)
                                       paste(Filter(Negate(is.null), lapply(seq_along(elements), function(i) {
                                         if (is.na(vals[i])) return(NULL)
                                         sprintf("    %-4s  %.6e", elements[i], vals[i])
                                       })), collapse = "\n")
                                     }

                                     soln1_elem_lines <- make_soln_lines(row_mix1)
                                     soln2_elem_lines <- make_soln_lines(row_mix2)
                                     soln1_pH   <- as.numeric(row_mix1[["pH"]])
                                     soln1_pe   <- as.numeric(row_mix1[["pe"]])
                                     soln1_h2o  <- as.numeric(row_mix1[["mass_H2O"]])
                                     soln2_pH   <- as.numeric(row_mix2[["pH"]])
                                     soln2_pe   <- as.numeric(row_mix2[["pe"]])
                                     soln2_h2o  <- as.numeric(row_mix2[["mass_H2O"]])

                                     # ---- Secondary minerals list ----
                                     secondary_minerals <- c(
                                       "Analcime", "Anhydrite", "Aragonite", "Bassanite",
                                       "Beidellite-Ca", "Beidellite-Fe", "Beidellite-Mg", "Beidellite-Na",
                                       "Boehmite", "Brucite",
                                       "C2H4(g)", "C2H6(g)", "C3H8(g)", "Calcite", "CH4(g)",
                                       "Chalcedony", "Chamosite", "Chrysotile", "Citric_Acid",
                                       "Clinochlore-14A", "Clinochlore-7A", "Clinoptilolite-Ca", "Clinoptilolite-Na",
                                       "Clinozoisite", "CO(g)", "CO2(g)", "Cronstedtite-7A",
                                       "Daphnite-14A", "Daphnite-7A", "Dawsonite", "Diaspore", "Dolomite",
                                       "Ettringite", "Fe(OH)2", "Fe(OH)3", "Fe2(SO4)3", "FeSO4",
                                       "Gibbsite", "Glycine", "Goethite", "Greenalite", "Gypsum", "Gyrolite",
                                       "H2(g)", "H2S(g)", "Hematite", "Huntite", "Hydromagnesite",
                                       "Kaolinite", "Lizardite", "Melanterite", "MgSO4", "Minnesotaite", "Mirabilite",
                                       "Montmor-Ca", "Montmor-Mg", "Montmor-Na",
                                       "N2(g)", "NH3(g)", "Nontronite-Ca", "Nontronite-Mg", "Nontronite-Na",
                                       "Portlandite", "Pyrite",
                                       "Saponite-Fe-Ca", "Saponite-Fe-Fe", "Saponite-Fe-Mg", "Saponite-Fe-Na",
                                       "Saponite-Mg-Ca", "Saponite-Mg-Fe", "Saponite-Mg-Mg", "Saponite-Mg-Na",
                                       "Sepiolite", "SiO2(am)", "SO2(g)", "Thenardite", "Tobermorite-11A"
                                     )
                                     # Two solutions here, gated by their own temperatures: porewater_temp_C
                                     # for EQUILIBRIUM_PHASES 1, ocean_temp_C for EQUILIBRIUM_PHASES 2.
                                     # Clathrate phases stable at each reservoir's own conditions are
                                     # appended the same way.
                                     extra_pw_minerals    <- c(self$temp_gated_minerals(porewater_temp_C),
                                                               self$stable_clathrate_phases(porewater_temp_C, pore_pressure_atm, ch4_mode))
                                     extra_ocean_minerals <- c(self$temp_gated_minerals(ocean_temp_C),
                                                               self$stable_clathrate_phases(ocean_temp_C, ocean_pressure_atm, ch4_mode))

                                     # ---- RATES block ----
                                     .rate_laws <- list(
                                       Schreibersite = SCHREIBERSITE_RATE_BODY,
                                       Forsterite = paste(
                                         '\t-start',
                                         '\t1   REM Ref PK04',
                                         '\t10  kacid = 10^(-6.85) * exp(-67.2e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^0.47',
                                         '\t20  kneut = 10^(-10.64) * exp(-79.0e3/8.314 * (1/TK-1/298.15))',
                                         '\t21  SSA = 0.1',
                                         '\t22  mw = 140.6715',
                                         '\t40  k = (kacid + kneut) * SSA * mw * M',
                                         '\t50  IF SR("Forsterite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Forsterite"))',
                                         '\t60  moles = rate * TIME',
                                         '\t70  SAVE moles',
                                         '\t-end', sep = "\n"),
                                       Troilite = paste(
                                         '\t-start',
                                         '\t1   REM PK04 has no separate troilite entry; hexagonal pyrrhotite parameters used as proxy',
                                         '\t25  kacid = 10^(-6.79) * exp(-63.0e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^-0.090 * ACT("Fe+3")^0.356',
                                         '\t26  SSA = 5',
                                         '\t27  mw = 87.913',
                                         '\t30  k = kacid * SSA * mw * M',
                                         '\t40  IF SR("Troilite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Troilite"))',
                                         '\t50  moles = rate * TIME',
                                         '\t55  IF moles < 0 THEN moles = 0',
                                         '\t60  SAVE moles',
                                         '\t-end', sep = "\n"),
                                       Diopside = paste(
                                         '\t-start',
                                         '\t1   REM Ref PK04',
                                         '\t10  kacid = 10^(-6.36) * exp(-96.1e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^0.71',
                                         '\t20  kneut = 10^(-11.11) * exp(-40.6e3/8.314 * (1/TK-1/298.15))',
                                         '\t21  SSA = 0.1',
                                         '\t22  mw = 216.55',
                                         '\t40  k = (kacid + kneut) * SSA * mw * M',
                                         '\t50  IF SR("Diopside") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Diopside"))',
                                         '\t60  moles = rate * TIME',
                                         '\t70  SAVE moles',
                                         '\t-end', sep = "\n"),
                                       `K-Feldspar` = paste(
                                         '\t-start',
                                         '\t1   REM Ref PK04',
                                         '\t10  kacid = 10^(-10.06) * exp(-51.7e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^0.5',
                                         '\t20  kneut = 10^(-12.41) * exp(-38.0e3/8.314 * (1/TK-1/298.15))',
                                         '\t30  kbase = 10^(-21.20) * exp(-94.1e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^-0.823',
                                         '\t31  SSA = 5',
                                         '\t32  mw = 278.33',
                                         '\t40  k = (kacid + kneut + kbase) * SSA * mw * M',
                                         '\t50  IF SR("K-Feldspar") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("K-Feldspar"))',
                                         '\t60  moles = rate * TIME',
                                         '\t70  SAVE moles',
                                         '\t-end', sep = "\n"),
                                       Tephroite = paste(
                                         '\t-start',
                                         '\t1   REM TODO awaiting kacid/kneut from user (Mn-olivine, no rate constants supplied yet)',
                                         '\t10  kacid = 0',
                                         '\t20  kneut = 0',
                                         '\t21  SSA = 0.1',
                                         '\t22  mw = 201.96',
                                         '\t40  k = (kacid + kneut) * SSA * mw * M',
                                         '\t50  IF SR("Tephroite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Tephroite"))',
                                         '\t60  moles = rate * TIME',
                                         '\t70  SAVE moles',
                                         '\t-end', sep = "\n"),
                                       Enstatite = paste(
                                         '\t-start',
                                         '\t1   REM Ref PK04',
                                         '\t10  kacid = 10^(-9.02) * exp(-80.0e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^0.6',
                                         '\t20  kneut = 10^(-12.72) * exp(-80.0e3/8.314 * (1/TK-1/298.15))',
                                         '\t21  SSA = 0.1',
                                         '\t22  mw = 100.3725',
                                         '\t40  k = (kacid + kneut) * SSA * mw * M',
                                         '\t50  IF SR("Enstatite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Enstatite"))',
                                         '\t60  moles = rate * TIME',
                                         '\t70  SAVE moles',
                                         '\t-end', sep = "\n"),
                                       Anorthite = paste(
                                         '\t-start',
                                         '\t1   REM Ref PK04',
                                         '\t10  kacid = 10^(-3.50) * exp(-16.6e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^1.411',
                                         '\t20  kneut = 10^(-9.12) * exp(-17.8e3/8.314 * (1/TK-1/298.15))',
                                         '\t21  SSA = 5',
                                         '\t22  mw = 278.164',
                                         '\t40  k = (kacid + kneut) * SSA * mw * M',
                                         '\t50  IF SR("Anorthite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Anorthite"))',
                                         '\t60  moles = rate * TIME',
                                         '\t70  SAVE moles',
                                         '\t-end', sep = "\n"),
                                       Albite = paste(
                                         '\t-start',
                                         '\t1   REM Ref PK04',
                                         '\t10  kacid = 10^(-10.16) * exp(-65.0e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^0.457',
                                         '\t20  kneut = 10^(-12.56) * exp(-69.8e3/8.314 * (1/TK-1/298.15))',
                                         '\t30  kbase = 10^(-15.60) * exp(-71.0e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^-0.572',
                                         '\t31  SSA = 5',
                                         '\t32  mw = 262.1798',
                                         '\t40  k = (kacid + kneut + kbase) * SSA * mw * M',
                                         '\t50  IF SR("Albite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Albite"))',
                                         '\t60  moles = rate * TIME',
                                         '\t70  SAVE moles',
                                         '\t-end', sep = "\n")
                                     )

                                     # ---- IOM kinetic phases (PATCH C/D2): Miller-calibrated channels,
                                     #      remaining moles recovered from IOM_*_mol columns (PATCH D1). ----
                                     iom_rates_block          <- ""
                                     iom_kinetics_block       <- ""
                                     iom_so_line              <- ""
                                     iom_punch_headings_tidal <- ""
                                     iom_punch_lines_tidal    <- ""

                                     if (self$organic_wt_percent > 0) {
                                       # See the main-path comment (iom_default_config()) for why this is no
                                       # longer a hardcoded table -- WIRING CHANGE 2026-09-26.
                                       if (is.null(self$iom_pools)) {
                                         self$iom_pools <- iom_default_config()
                                       }

                                       # Guard: IOM_*_mol columns must be present — written by PATCH D1.
                                       iom_mol_cols <- paste0(self$iom_pools$name, "_mol")
                                       missing_iom  <- iom_mol_cols[!iom_mol_cols %in% names(row_react)]
                                       if (length(missing_iom) > 0)
                                         stop(paste0(
                                           "Restart error: organic_wt_percent > 0 but IOM columns missing from prev_csv:\n  ",
                                           paste(missing_iom, collapse = ", "), "\n",
                                           "The CSV must come from a run with PATCH D1 applied ",
                                           "(IOM channels punched to output)."
                                         ))

                                       get_iom_moles <- function(channel) {
                                         val <- suppressWarnings(as.numeric(row_react[[paste0(channel, "_mol")]]))
                                         if (!is.na(val) && val >= 0) val else 0
                                       }

                                       iom_rates_block <- paste0(apply(self$iom_pools, 1, function(p) {
                                         sprintf(paste0(
                                           "%s\n -start\n",
                                           " 1  REM Arrhenius 1st-order kerogen degradation (single-Ea pool)\n",
                                           " 10 k = 10^(%s) * exp(-%s / (8.314 * TK))\n",
                                           " 20 rate = k * M\n",
                                           " 30 if (M <= 0) then rate = 0\n",
                                           " 40 moles = rate * TIME\n",
                                           " 50 SAVE moles\n -end\n"),
                                           p[["name"]], p[["logA"]], p[["Ea_J"]])
                                       }), collapse = "")

                                       iom_kinetics_block <- paste0(apply(self$iom_pools, 1, function(p) {
                                         m_remaining <- get_iom_moles(p[["name"]])
                                         sprintf("   %s\n      -formula %s\n      -m0 %.6e\n      -step_divide 1\n",
                                                 p[["name"]], trimws(p[["formula"]]), m_remaining)
                                       }), collapse = "")

                                       iom_so_line <- paste(self$iom_pools$name, collapse = " ")

                                       iom_punch_headings_tidal <- paste0(" ",
                                         paste(sprintf("%s_mol", self$iom_pools$name), collapse = " "))
                                       iom_punch_lines_tidal <- paste0(
                                         vapply(seq_len(nrow(self$iom_pools)), function(j) {
                                           sprintf("   %d PUNCH KIN(\"%s\")\n",
                                                   200 + 10 * j, self$iom_pools$name[j])
                                         }, character(1)), collapse = "")

                                       message(sprintf("IOM restart: %d channels recovered from prev_csv (PATCH D1 columns)",
                                                       nrow(self$iom_pools)))
                                     }

                                     rates_block <- if (length(kinetic_minerals) > 0)
                                       paste0("RATES\n",
                                              paste(sapply(kinetic_minerals, function(m)
                                                paste0(m, "\n", .rate_laws[[m]])),
                                                collapse = "\n\n"))
                                     else ""
                                     if (nzchar(iom_rates_block)) {
                                       rates_block <- if (nzchar(rates_block))
                                         paste0(rates_block, "\n\n", iom_rates_block)
                                       else
                                         paste0("RATES\n", iom_rates_block)
                                     }

                                     # ---- KINETICS 1 block ----
                                     kinetics_block <- if (length(kinetic_minerals) > 0 || nzchar(iom_kinetics_block)) {
                                       min_entries <- if (length(kinetic_minerals) > 0)
                                         paste(sapply(kinetic_minerals, function(m)
                                           if (m == "Schreibersite") schreibersite_kinetics_entry(get_k_moles(m)) else sprintf("   %s\n      -m0 %.6e\n      -step_divide 1", m, get_k_moles(m))),
                                           collapse = "\n")
                                       else ""
                                       all_entries <- paste(Filter(nzchar, c(min_entries, iom_kinetics_block)),
                                                            collapse = "\n")
                                       paste0("KINETICS 1\n", all_entries, "\n",
                                              sprintf("   -steps %.6e in %d steps\n", time_step_seconds, substeps),
                                              "   -cvode true\n",
                                              "   -bad_step_max 5000")
                                     } else ""

                                     # ---- EQUILIBRIUM_PHASES 1 (porewater: EP1 state + equil primaries) ----
                                     # Fe/Ni/Lawrencite always appear here (porewater), recovered via
                                     # get_ep_moles() from their own native EQUILIBRIUM_PHASES output column.
                                     ep1_prim_lines <- paste(
                                       c("# --- Primary minerals (equilibrium) ---",
                                         sapply(always_equil_primary, function(m)
                                           sprintf("   %-20s 0 %.6e  dissolve_only", m, get_ep_moles(m, row_react))),
                                         if (length(equil_primary) > 0)
                                           sapply(equil_primary, function(m)
                                             sprintf("   %-20s 0 %.6e  dissolve_only", m, get_k_moles(m)))
                                         else character(0)),
                                       collapse = "\n")

                                     ep1_sec_lines <- paste(
                                       c("# --- Secondary minerals (porewater, from previous run) ---",
                                         sapply(c(secondary_minerals, extra_pw_minerals), function(m)
                                           sprintf("   %-20s 0 %.6e", m, get_ep_moles(m, row_react)))),
                                       collapse = "\n"
                                     )

                                     ep1_block <- paste0("EQUILIBRIUM_PHASES 1\n",
                                       paste(Filter(nzchar, c(ep1_prim_lines, ep1_sec_lines)), collapse = "\n"))

                                     # ---- EQUILIBRIUM_PHASES 2 (ocean: EP2 state, equil primaries + Magnetite + secondaries) ----
                                     ep2_prim_lines <- if (length(equil_primary) > 0)
                                       paste(c("# --- Primary minerals (equilibrium, ocean) ---",
                                               sapply(equil_primary, function(m)
                                                 sprintf("   %-20s 0 %.6e  dissolve_only", m, get_ep_moles(m, row_mix4)))),
                                             collapse = "\n")
                                     else ""
                                     ep2_minerals <- c(append(secondary_minerals, "Magnetite",
                                                              after = which(secondary_minerals == "Lizardite")),
                                                        extra_ocean_minerals)
                                     ep2_sec_lines <- paste(
                                       c("# --- Secondary minerals (ocean, from previous run) ---",
                                         sapply(ep2_minerals, function(m)
                                           sprintf("   %-20s 0 %.6e", m, get_ep_moles(m, row_mix4)))),
                                       collapse = "\n"
                                     )
                                     ep2_block <- paste0("EQUILIBRIUM_PHASES 2\n",
                                       paste(Filter(nzchar, c(ep2_prim_lines, ep2_sec_lines)), collapse = "\n"))

                                     # ---- SELECTED_OUTPUT dynamic lines ----
                                     .kin_names  <- c(kinetic_minerals,
                                                       if (nzchar(iom_so_line)) self$iom_pools$name else character(0))
                                     kin_so_line <- if (length(.kin_names) > 0)
                                       paste("    -kinetics           ", paste(.kin_names, collapse = " "))
                                     else ""

                                     .all_eq_so  <- c(always_equil_primary, equil_primary,
                                                       append(secondary_minerals, "Magnetite",
                                                              after = which(secondary_minerals == "Lizardite")),
                                                       union(extra_pw_minerals, extra_ocean_minerals))
                                     .eq_chunks  <- split(.all_eq_so, ceiling(seq_along(.all_eq_so) / 7))
                                     equil_so_line <- paste(
                                       c(paste0("    -equilibrium_phases  ", paste(.eq_chunks[[1]], collapse = " ")),
                                         if (length(.eq_chunks) > 1)
                                           sapply(.eq_chunks[-1], function(ch)
                                             paste0("                         ", paste(ch, collapse = " ")))),
                                       collapse = "\n"
                                     )

                                     # ---- USER_PRINT early-exit: adapted to kinetic_minerals ----
                                     # Depletion mineral is user-specified; proxy is the last kinetic mineral
                                     # (longest-lived, used to confirm we are in a REACT simulation).
                                     .iom_names <- if (self$organic_wt_percent > 0 && !is.null(self$iom_pools))
                                       self$iom_pools$name else character(0)
                                     .all_kinetic <- c(kinetic_minerals, .iom_names)
                                     # Provisional default order (smallest initial abundance first, so it's
                                     # the most likely to deplete first; largest last, as the most robust
                                     # REACT/MIX proxy) — unverified against actual depletion behavior for
                                     # this composition; pass depletion_mineral explicitly to override.
                                     depletion_order <- c("K-Feldspar", "Tephroite", "Anorthite", "Diopside",
                                                          "Albite", "Forsterite", "Enstatite", "Troilite",
                                                          "IOM_CO2", "IOM_CH4", "IOM_N", "IOM_S")
                                     if (!is.null(depletion_mineral)) {
                                       if (!depletion_mineral %in% .all_kinetic)
                                         stop(sprintf("depletion_mineral '%s' is not a kinetic phase in this restart (check equil_minerals and organic_wt_percent)", depletion_mineral))
                                     } else {
                                       depletion_mineral <- head(intersect(depletion_order, .all_kinetic), 1)
                                     }
                                     proxy_mineral <- tail(intersect(depletion_order, .all_kinetic), 1)

                                     # ---- USER_PUNCH total_initial expression ----
                                     total_initial_expr <- paste(
                                       sapply(all_primary, function(m)
                                         sprintf("%.6f*%.4f", get_k_moles(m), molar_masses[[m]])),
                                       collapse = " + "
                                     )

                                     # Punch function: KIN() for kinetic minerals, EQUI() for equil_primary
                                     make_punch_line <- function(line_num, mineral) {
                                       mm <- molar_masses[[mineral]]
                                       fn <- if (mineral %in% equil_primary)
                                         sprintf('EQUI("%s")', mineral) else sprintf('KIN("%s")', mineral)
                                       sprintf("   %d PUNCH %s * %.4f", line_num, fn, mm)
                                     }

                                     base_name <- tools::file_path_sans_ext(basename(filename))

                                     # ---- Template header (static part) ----
                                     template <- glue::glue('
TITLE Enceladus tidal kinetic (restart)
{ch4_db_note_line}
{clathrate_note_lines}

KNOBS
   -logfile true

# k={self$simulator$params$k}, thickness={self$simulator$params$layer_thickness} m
# Mode: tidal kinetic restart from previous CSV state
# prev_csv: {prev_csv}
# restart_time_offset: {sprintf("%.6e", restart_time_years)} years
# equil_minerals: {paste(equil_primary, collapse = ", ")}
# n_orbital_cycles: {n_orbital_cycles}
# substeps: {substeps}
# mode: {mode}
# time_step: {sprintf("%.4e", time_step_years)} years

SOLUTION 1 Enceladus porewater (restarted)
    temp      {round(porewater_temp_C, 2)}
    pH        {sprintf("%.4f", soln1_pH)}
    pe        {sprintf("%.4f", soln1_pe)}
    redox     pe
    units     mol/kgw
    density   1
    water     {sprintf("%.6f", soln1_h2o)}
    pressure  {round(ocean_pressure_atm, 1)}
{soln1_elem_lines}

SOLUTION 2 Enceladus ocean (restarted)
    temp      {round(ocean_temp_C, 2)}
    pH        {sprintf("%.4f", soln2_pH)}
    pe        {sprintf("%.4f", soln2_pe)}
    redox     pe
    units     mol/kgw
    density   1
    water     {sprintf("%.6f", soln2_h2o)}
    pressure  {round(ocean_pressure_atm, 1)}
{soln2_elem_lines}

REACTION_TEMPERATURE 1
    {round(porewater_temp_C, 2)}

REACTION_TEMPERATURE 2
    {round(ocean_temp_C, 2)}

{rates_block}

{kinetics_block}

{ep1_block}

{ep2_block}

SELECTED_OUTPUT
    -file {base_name}.csv
    -reset false
    -reaction true
    -state true
    -solution true
    -step true
    -ph true
    -pe true
    -alkalinity true
    -ionic_strength true
    -water true
    -totals              Al C Ca Fe H Mg N Na O P S Si
    -molalities          OH- H+ C2H4 C2H6 CH4 CO HCO3- CO3-2 CO2
                         CH3COO- HCOO- HCN Ca+2 Cl- Fe+2 FeOH+ Fe+3 H2
                         Mg+2 NH4+ NH3 NH4CO3- N2 NO2-
                         NO3- Na+ H2PO4- HPO4-2 PO4-3 HS- H2S SO3-2 HSO3-
                         SO2 SO4-2 SiO2 HSiO3-
    -activities          OH- H+ C2H4 C2H6 CH4 CO HCO3- CO3-2 CO2
                         CH3COO- HCOO- HCN Ca+2 Cl- Fe+2 FeOH+ Fe+3 H2
                         Mg+2 NH4+ NH3 NH4CO3- N2 NO2-
                         NO3- Na+ H2PO4- HPO4-2 PO4-3 HS- H2S SO3-2 HSO3-
                         SO2 SO4-2 SiO2 HSiO3-
{kin_so_line}
{equil_so_line}
    -saturation_indices  {paste(setdiff(all_primary, "Schreibersite"), collapse = " ")}
    -gases               CO2(g) H2(g) CH4(g) CO(g) H2O(g) H2S(g) N2(g)
                         NH3(g) NO(g) NO2(g) O2(g) SO2(g)

# USER_PUNCH is defined per-cycle below (with exact cumulative time embedded)
    ')

                                     # ---- Cycle loop ----
                                     cycle_blocks <- paste(sapply(seq_len(kinetic_cycles), function(i) {
                                       cumulative_time <- restart_time_years + i * time_step_years
                                       paste0(
                                         sprintf("\n# --- Cycle %d of %d (cumulative time: %.4e years) ---\n",
                                                 i, kinetic_cycles, cumulative_time),
                                         "USE solution 1\n",
                                         "USE kinetics 1\n",
                                         "USE equilibrium_phases 1\n",
                                         "USE reaction none\n",
                                         "USE reaction_temperature 1\n",
                                         "SAVE solution 1\n",
                                         "SAVE equilibrium_phases 1\n",
                                         "USER_PUNCH\n",
                                         "   -headings Time_Years Enst_remain_g Forst_remain_g Troil_remain_g Alb_remain_g Diop_remain_g Anorth_remain_g Kfs_remain_g Teph_remain_g Schr_remain_g Total_Initial_g",
                                         iom_punch_headings_tidal, "\n",
                                         sprintf("   10 PUNCH %.6f  # Cumulative time (years)\n", cumulative_time),
                                         paste(sapply(seq_along(all_primary), function(j) {
                                           make_punch_line(10 + j * 10, all_primary[j])
                                         }), collapse = "\n"), "\n",
                                         sprintf("   %d total_initial = %s\n", 10 + (length(all_primary) + 1) * 10, total_initial_expr),
                                         sprintf("   %d PUNCH total_initial\n", 10 + (length(all_primary) + 2) * 10),
                                         iom_punch_lines_tidal,
                                         "END\n\n",
                                         "USE kinetics none\n",
                                         "USE equilibrium_phases 2\n",
                                         "USE reaction none\n",
                                         "USE reaction_temperature 2\n",
                                         "MIX 4\n",
                                         sprintf("  1  %.6e\n", x_step),
                                         "  2  1.0\n",
                                         "SAVE solution 4\n",
                                         "SAVE equilibrium_phases 2\n",
                                         "END\n\n",
                                         "USE kinetics none\n",
                                         "USE equilibrium_phases none\n",
                                         "USE reaction none\n",
                                         "USE reaction_temperature 1\n",
                                         "MIX 1\n",
                                         sprintf("  1  %.6e\n", 1 - x_step),
                                         sprintf("  4  %.6e\n", f_x / (f_ocean + f_x)),
                                         "SAVE solution 1\n",
                                         "END\n\n",
                                         "USE kinetics none\n",
                                         "USE equilibrium_phases none\n",
                                         "USE reaction none\n",
                                         "MIX 2\n",
                                         sprintf("  4  %.6e\n", f_ocean / (f_ocean + f_x)),
                                         "SAVE solution 2\n",
                                         "END\n",
                                         if (i < kinetic_cycles &&
                                             length(depletion_mineral) == 1 && length(proxy_mineral) == 1 &&
                                             depletion_mineral != proxy_mineral)
                                           paste0(
                                             "USER_PRINT\n",
                                             sprintf('   10 IF KIN("%s") > 1e-9 THEN GOTO 30\n', depletion_mineral),
                                             sprintf('   20 IF KIN("%s") > 0 THEN STOP\n', proxy_mineral),
                                             "   30 REM\n"
                                           )
                                         else ""
                                       )
                                     }), collapse = "\n")

                                     content <- self$apply_ch4_mode(paste0(template, "\n", cycle_blocks), ch4_mode)
                                     writeLines(content, filename)
                                     return(filename)
                                   }
                                 )
)

