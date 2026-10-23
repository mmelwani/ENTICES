#===============================================================================
# ENTICES: Enceladus Tidally-Induced Circulation and Exchange Simulator
# Author: Mohit Melwani Daswani
# Last updated: 2025-08-25
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
# - Liao et al. (2020) JGR Planets
# - Zandanel et al. (2021) Icarus
#===============================================================================

#===============================================================================
# Required libraries 
#===============================================================================
library(tidyverse)
library(patchwork)  # For plot layouts
library(glue)       # For templating
library(R6)

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
# Configuration Class for Enceladus Parameters
#===============================================================================
EnceladusConfig <- R6::R6Class("EnceladusConfig",
                               public = list(
                                 # Basic dimensional parameters
                                 radius = 252.1e3,        # Surface radius [m]
                                 core_radius = 194.2e3,   # Core radius [m]
                                 mass = 1.0802e20,        # Total mass [kg]
                                 
                                 # Material properties
                                 core_density = 2353,     # Core bulk density [kg/m³]
                                 core_porosity = 0.32,    # Core porosity [volume fraction]
                                 
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
                                           CONSTANTS$FLUID$density * (r^3 - self$core_radius^3)
                                       ) / r^2
                                     }
                                   }
                                   
                                   return(data.frame(
                                     radius_m = depth_radius,
                                     depth_m = depth_m,
                                     gravity_ms2 = gravity_ms2
                                   ))
                                 },
                                 
                                 calculate_gravity = function(r) {
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
                                   core_mass <- core_volume * self$core_density * (1 - self$core_porosity)
                                   
                                   gravity_seafloor <- self$calculate_gravity(self$core_radius)
                                   seafloor_pressure <- CONSTANTS$FLUID$density * gravity_seafloor * self$ocean_thickness
                                   
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
    return(layer_thickness * core_radius)
  },
  
  # Calculate vertical area for volume calculations  
  calculate_vertical_area = function(core_radius) {
    return(pi * core_radius^2)
  },
  
  # Calculate depth-dependent porosity
  calculate_mean_porosity = function(surface_porosity, depth_factor = 0.8) {
    base_porosity <- surface_porosity * depth_factor
    return((surface_porosity + base_porosity) / 2)
  }
)

FlowUtils <- list(
  # Calculate Darcy flow
  calculate_darcy_flow = function(k, grad_P, area, porosity, 
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
                             stress <- self$config$max_tidal_stress * sin(phase)
                             
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
                             
                             metrics <- list(
                               max_flow_rate_Ls = max(self$results$flow_rate_Ls),
                               mean_flow_rate_Ls = mean(self$results$flow_rate_Ls),
                               total_fluid_volume_m3 = sum(self$results$fluid_volume),
                               total_rock_reacted = sum(self$results$rock_mass_reacted),
                               circulation_time_years = FlowUtils$calculate_circulation_time(
                                 flow_rate = mean(self$results$flow_rate_m3s),
                                 ocean_mass = self$config$get_derived_properties()$ocean_mass
                               )
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
                                     
                                     # Calculate layer fraction
                                     layer_fraction <- self$simulator$params$layer_thickness / self$simulator$config$core_radius
                                     
                                     # Use empirical bounds from your analysis
                                     min_circulation_time <- 0.1        # Fastest: k=1e-08, 10m layer
                                     max_circulation_time <- 2.45e9     # Slowest: k=1e-14, full core  
                                     min_layer_fraction <- 10 / self$simulator$config$core_radius     # 10m layer
                                     max_layer_fraction <- 1.0          # Full core
                                     
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
                                     
                                     # Calculate parameters
                                     rock_per_cycle <- metrics$total_rock_reacted / props$core_mass
                                     
                                     # Only calculate rock_per_year if circulation time exists
                                     if(!is.null(metrics$circulation_time_years) && length(metrics$circulation_time_years) > 0) {
                                       rock_per_year <- 1.0 / metrics$circulation_time_years
                                     } else {
                                       rock_per_year <- NA
                                       warning("Circulation time not calculated, cannot determine rock_per_year")
                                     }
                                     
                                     water_mass <- props$ocean_mass / props$core_mass
                                     
                                     # Store parameters
                                     # Store organic percentage
                                     self$organic_wt_percent <- organic_wt_percent
                                     self$phreeqc_params <- list(
                                       grouping_factor = grouping_factor,
                                       total_rock = layer_fraction,
                                       rock_per_cycle = rock_per_cycle,
                                       rock_per_year = rock_per_year,
                                       water_mass = water_mass
                                     )
                                     
                                     print("Calculated PHREEQC parameters:")
                                     print(sprintf("rock_per_cycle: %.2e", rock_per_cycle))
                                     print(sprintf("rock_per_year: %.2e", rock_per_year))
                                     print(sprintf("water_mass: %.2f", water_mass))
                                     
                                     return(self$phreeqc_params)
                                   },
                                   
                                   generate_phreeqc_input = function(output_dir = "phreeqc_inputs", n_circulations = 1) {
                                     # Create output directory if needed
                                     dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
                                     
                                     # Create filename
                                     filename <- file.path(output_dir, sprintf(
                                       "enceladus_k%.2e_d%g_t%.2e_circ%d.pqi",
                                       self$simulator$params$k,
                                       self$simulator$params$layer_thickness,
                                       self$phreeqc_params$grouping_factor,
                                       n_circulations
                                     ))
                                     
                                     
                                     # Generate input file content
                                     content <- self$create_phreeqc_template(filename, n_circulations)
                                     
                                     # Write to file
                                     writeLines(content, filename)
                                     
                                     return(filename)
                                   },
                                   
                                   create_phreeqc_template = function(filename, n_circulations = 1) {
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
                                     
                                     # Adjust for grouping
                                     adjusted_cycles <- ceiling(total_cycles / self$phreeqc_params$grouping_factor)
                                     
                                     # Extend total simulation time
                                     extended_reaction_time <- total_reaction_time * n_circulations
                                     extended_cycles <- ceiling((extended_reaction_time * 365.25 * 24 * 3600) / orbital_period)
                                     
                                     # But keep organic addition to original timeframe
                                     # Original cycles for organic addition (1 circulation)
                                     organic_cycles <- adjusted_cycles
                                     
                                     # Extended cycles for kinetic simulation
                                     kinetic_cycles <- ceiling(adjusted_cycles * n_circulations)
                                     
                                     # Calculate time step for kinetics (years)
                                     time_step_years <- metrics$circulation_time_years / adjusted_cycles
                                     
                                     # Get base filename without extension for output
                                     base_name <- tools::file_path_sans_ext(basename(filename))
                                     
                                     # Calculate scaling factor for minerals (reduced by organic content)
                                     mineral_scale_factor <- (100 - self$organic_wt_percent) / 100
                                     
                                     # Calculate total mineral amounts (moles) based on permeable layer fraction
                                     forsterite_moles <- 2.8833 * self$phreeqc_params$total_rock * mineral_scale_factor
                                     fayalite_moles <- 0.2212 * self$phreeqc_params$total_rock * mineral_scale_factor
                                     enstatite_moles <- 2.9845 * self$phreeqc_params$total_rock * mineral_scale_factor
                                     ferrosilite_moles <- 0.6569 * self$phreeqc_params$total_rock * mineral_scale_factor
                                     pyrrhotite_moles <- 0.7981 * self$phreeqc_params$total_rock * mineral_scale_factor
                                     anorthite_moles <- 0.2323 * self$phreeqc_params$total_rock * mineral_scale_factor
                                     albite_moles <- 0.0616 * self$phreeqc_params$total_rock * mineral_scale_factor
                                     magnetite_moles <- 0.0527 * self$phreeqc_params$total_rock * mineral_scale_factor
                                     
                                     iom_reaction_block <- ""
                                     if (self$organic_wt_percent > 0) {
                                       # IOM composition per kg: 728g C, 46g H, 159g O, 30g N, 37g S
                                       iom_molar_composition_per_kg <- list(
                                         C = 728 / 12.011,      # 60.63 mol C/kg IOM
                                         H2 = (46 / 1.0079)/2,  # 22.82 mol H2/kg IOM  
                                         O = 159 / 15.994,      # 9.94 mol [(aro)-O-(aro)]/kg IOM
                                         N2 = (30 / 14.007)/2,  # 1.07 mol N2/kg IOM
                                         S = 37 / 32.066        # 1.15 mol [(6)(CB)(CB)S]/kg IOM
                                       )
                                       
                                       # Calculate total IOM mass in the permeable layer (PHREEQC scaled)
                                       iom_mass_kg_scaled <- (self$organic_wt_percent / 100) * self$phreeqc_params$total_rock
                                       
                                       # Calculate moles of each element from IOM (using scaled mass)
                                       iom_C_moles <- iom_molar_composition_per_kg$C * iom_mass_kg_scaled
                                       iom_H2_moles <- iom_molar_composition_per_kg$H2 * iom_mass_kg_scaled
                                       iom_O_moles <- iom_molar_composition_per_kg$O * iom_mass_kg_scaled
                                       iom_N2_moles <- iom_molar_composition_per_kg$N2 * iom_mass_kg_scaled
                                       iom_S_moles <- iom_molar_composition_per_kg$S * iom_mass_kg_scaled
                                       
                                       # Create the reaction block
                                       total_organic_moles <- max(iom_C_moles, iom_H2_moles, iom_O_moles, iom_N2_moles, iom_S_moles)
                                       
                                       iom_reaction_block <- paste0("REACTION 1\n",
                                                                    "   C                 ", sprintf("%.6f", iom_C_moles), "\n",
                                                                    "   H2                ", sprintf("%.6f", iom_H2_moles), "\n", 
                                                                    "   [(aro)-O-(aro)]   ", sprintf("%.6f", iom_O_moles), "\n",
                                                                    "   N2                ", sprintf("%.6f", iom_N2_moles), "\n",
                                                                    "   [(6)(CB)(CB)S]    ", sprintf("%.6f", iom_S_moles), "\n",
                                                                    # Single increment over ALL kinetic steps
                                                                    "   ", sprintf("%.6f", total_organic_moles), " moles in ", kinetic_cycles, " steps\n")
                                       
                                       # Debug output for organic mass calculation
                                       print(sprintf("Debug organic calculation:"))
                                       print(sprintf("Organic mass (PHREEQC scaled): %.2e kg", iom_mass_kg_scaled))
                                     }
                                     
                                     template <- glue::glue('
TITLE Enceladus hydrothermal simulation 
# k={self$simulator$params$k}, thickness={self$simulator$params$layer_thickness} m
# Organic content: {self$organic_wt_percent} wt.%
# Simulation parameters:
# Layer thickness: {self$simulator$params$layer_thickness} m ({sprintf("%.1f", self$simulator$params$layer_thickness/self$simulator$config$core_radius*100)}% of core)
# Available rock mass: {self$phreeqc_params$total_rock} kg
# Reaction time: {metrics$circulation_time_years} years
# Total cycles: {adjusted_cycles} 
# Time step: {sprintf("%.2e", time_step_years)} years per step
# Mean flow rate: {sprintf("%.2e", metrics$mean_flow_rate_Ls)} L/s
# Circulation time: {sprintf("%.2e", metrics$circulation_time_years)} years

SOLUTION 1 Enceladus ocean
    temp      25
    pH        7
    pe        4
    redox     pe
    units     mol/kgw
    Ca        1e-6
    Mg        1e-6  
    Fe        1e-6
    Na        1e-6
    Si        1e-6
    Al        1e-6
    S         1e-6
    N         1e-6
    density   1
    water     {self$phreeqc_params$water_mass}
    pressure  70
    
# Define rate equations for each mineral
RATES
Forsterite
	-start
	1   REM Ref PK04
	10  kacid = 10^(-6.85) * exp(-67.2e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^0.47
	20  kneut = 10^(-10.64) * exp(-79.0e3/8.314 * (1/TK-1/298.15))
	40  k = kacid + kneut
	50  IF SR("Forsterite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Forsterite"))
	60  moles = rate * TIME
	70  SAVE moles
	-end
	
	Fayalite
	-start
	1   REM Ref PK04
	10  kacid = 10^(-4.80) * exp(-94.4e3/8.314 * (1/TK-1/298.15))
	20  kneut = 10^(-12.80) * exp(-94.4e3/8.314 * (1/TK-1/298.15))
	40  k = kacid + kneut
	50  IF SR("Fayalite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Fayalite"))
	60  moles = rate * TIME
	70  SAVE moles
	-end
	
	Enstatite
	-start
	1   REM Ref PK04
	10  kacid = 10^(-9.02) * exp(-80.0e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^0.6
	20  kneut = 10^(-12.72) * exp(-80.0e3/8.314 * (1/TK-1/298.15))
	40  k = kacid + kneut
	50  IF SR("Enstatite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Enstatite"))
	60  moles = rate * TIME
	70  SAVE moles
	-end
	
    Ferrosilite
        -start
        10  kacid = 10^(-8.30) * exp(-47.2e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^0.650
        20  kneut = 10^(-11.70) * exp(-66.1e3/8.314 * (1/TK-1/298.15))
        30  k = kacid + kneut
        40  IF SR("Ferrosilite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Ferrosilite"))
        50  moles = rate * TIME
        60  SAVE moles
        -end

Pyrrhotite # In lieu of Troilite
        -start
        1   REM Acid mechanism only (PK04)
        2   REM Two formulations available: monoclinic and hexagonal pyrrhotite
        3   REM Using hexagonal pyrrhotite parameters (least common form, but troilite is hexagonal)
        4   REM Rate depends on H+, Fe3+, and O2 activities
        5   REM Reaction orders: H+ = -0.597, Fe3+ = 0.355, O2 = not specified in neutral mechanism
        10  REM Monoclinic pyrrhotite: log k = -8.04, E = 50.8 kJ/mol
        15  REM kacid = 10^(-8.04) * exp(-50.8e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^-0.597 * ACT("Fe+3")^0.355
        20  REM Hexagonal pyrrhotite formulation reaction orders H+ = -0.090, Fe3+ = 0.356
        25  kacid = 10^(-6.79) * exp(-63.0e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^-0.090 * ACT("Fe+3")^0.356
        30  k = kacid
        40  IF SR("Pyrrhotite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Pyrrhotite"))
        50  moles = rate * TIME
        60  SAVE moles
        -end

    Anorthite # Ca-endmember of plagioclase
        -start
        1   REM Ref PK04
        10  kacid = 10^(-3.50) * exp(-16.6e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^1.411
        20  kneut = 10^(-9.12) * exp(-17.8e3/8.314 * (1/TK-1/298.15))
        40  k = kacid + kneut
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
        40  k = kacid + kneut + kbase
        50  IF SR("Albite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Albite"))
        60  moles = rate * TIME
        70  SAVE moles
        -end
        
Magnetite
	-start
	1   REM Ref PK04
	10   if (M < 0) then goto 200
	110  kacid = 10^(-8.59) * exp(-18.6e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^0.279
	120  kneut = 10^(-10.78) * exp(-18.6e3/8.314 * (1/TK-1/298.15))
	140  k = kacid + kneut
	150  rate = k * (1 - SR("Magnetite"))
	160  moles = rate * TIME
	200  SAVE moles
	-end
	
    
# Define kinetic mineral phases with initial amounts
KINETICS 1
   Forsterite
      -m0 {forsterite_moles}
      -step_divide 1
   Fayalite
      -m0 {fayalite_moles}
      -step_divide 1
   Enstatite
      -m0 {enstatite_moles}
      -step_divide 1
   Ferrosilite
      -m0 {ferrosilite_moles}
      -step_divide 1
   Pyrrhotite
      -m0 {pyrrhotite_moles}
      -step_divide 1
   Anorthite
      -m0 {anorthite_moles}
      -step_divide 1
   Albite
      -m0 {albite_moles}
      -step_divide 1
   Magnetite
      -m0 {magnetite_moles}
      -step_divide 1
   -step {time_step_years * 365.25 * 24 * 3600}  # Convert years to seconds
   -steps {kinetic_cycles}  # Extended simulation with ocean overturns > 1, when investigating full equilibration timescales
   -cvode true
   -bad_step_max 5000

# Add IOM as elemental composition via REACTION block (if present)
{iom_reaction_block}
    
EQUILIBRIUM_PHASES 1
#   Alabandite 0 0
   Alanine 0 0
#   Alum-K 0 0
#   Alunite 0 0
   Analcime 0 0
   Anhydrite 0 0
   Aragonite 0 0
#   Arcanite 0 0
   Artinite 0 0
   Bassanite 0 0
   Beidellite-Ca 0 0
   Beidellite-Fe 0 0
#   Beidellite-K 0 0
   Beidellite-Mg 0 0
   Beidellite-Na 0 0
   Boehmite 0 0
   Brucite 0 0
   Calcite 0 0
#   Celadonite 0 0
   Citric_Acid 0 0
   Clinoptilolite-Ca 0 0
#   Clinoptilolite-K 0 0
   Clinoptilolite-Na 0 0
   CO(g) 0 0
   CO2(g) 0 0
   Cronstedtite-7A 0 0
   Dawsonite 0 0
   Ettringite 0 0
   Fe(OH)2 0 0
   Fe(OH)3 0 0
   FeSO4 0 0
   Gibbsite 0 0
   Glycine 0 0
   Goethite 0 0
   Greenalite 0 0
   Gypsum 0 0
   H2(g) 0 0
   H2O(g) 0 0
   H2S(g) 0 0
#   Halite 0 0
#   Hausmannite 0 0
   Hematite 0 0
   Huntite 0 0
#   Hydroxyapatite 0 0
   Ice 0 0
#   Jarosite 0 0
#   KAl(SO4)2 0 0
   Kaolinite 0 0
   KerogenC128 0 0
   KerogenC292 0 0
   KerogenC406 0 0
   KerogenC415 0 0
   KerogenC515 0 0
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
   Monohydrocalcite 0 0
   Montmor-Ca 0 0
#   Montmor-K 0 0
   Montmor-Mg 0 0
   Montmor-Na 0 0
   N2(g) 0 0
   Na2CO3 0 0
   Na2CO3:7H2O 0 0
   Nahcolite 0 0
   Natron 0 0
   Nesquehonite 0 0
   NH3(g) 0 0
#   NH4Cl 0 0
   NH4HCO3 0 0
#   Niter 0 0
#   NO(g) 0 0
#   NO2(g) 0 0
   Nontronite-Ca 0 0
#   Nontronite-K 0 0
   Nontronite-Mg 0 0
   Nontronite-Na 0 0
   O2(g) 0 0
   Portlandite 0 0
   Pyridine 0 0
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
   Siderite 0 0
   SiO2(am) 0 0
#   Smectite-high-Fe-Mg 0 0
#   Smectite-low-Fe-Mg 0 0
   SO2(g) 0 0
#   Stilbite 0 0
#   Strengite 0 0
#   Sylvite 0 0
   Thenardite 0 0
   Thermonatrite 0 0
   Tobermorite-11A 0 0
    
SELECTED_OUTPUT 
    -file {base_name}.csv
    -reset false
    -reaction true
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
    -kinetics            Forsterite Fayalite Enstatite Ferrosilite Pyrrhotite Anorthite Albite Magnetite
    -equilibrium_phases  #Alabandite 
                         Alanine 
                         #Alum-K Alunite 
                         Analcime Anhydrite Aragonite #Arcanite 
                         Artinite 
                         Bassanite Beidellite-Ca Beidellite-Fe 
                         #Beidellite-K 
                         Beidellite-Mg Beidellite-Na Boehmite 
                         Brucite Calcite #Celadonite 
                         Citric_Acid 
                         Clinoptilolite-Ca #Clinoptilolite-K 
                         Clinoptilolite-Na CO(g) 
                         CO2(g) 
                         Cronstedtite-7A Dawsonite 
                         Ettringite 
                         Fe(OH)2 Fe(OH)3 FeSO4 Gibbsite 
                         Glycine 
                         Goethite Greenalite Gypsum 
                         H2(g) H2O(g) H2S(g) #Halite 
                         #Hausmannite 
                         Hematite Huntite #Hydroxyapatite 
                         Ice #Jarosite KAl(SO4)2 
                         Kaolinite 
                         KerogenC128 KerogenC292 KerogenC406 KerogenC415 
                         KerogenC515 
                         Melanterite #MgOHCl 
                         MgSO4 
                         Minnesotaite Mirabilite #Mn(OH)2(am) MnCl2:2H2O 
                         #MnCl2:4H2O MnCl2:H2O MnSO4 
                         Monohydrocalcite 
                         Montmor-Ca #Montmor-K 
                         Montmor-Mg Montmor-Na 
                         N2(g) Na2CO3 Na2CO3:7H2O Nahcolite 
                         Natron Nesquehonite NH3(g) #NH4Cl 
                         NH4HCO3 Niter NO(g) NO2(g) 
                         Nontronite-Ca #Nontronite-K 
                         Nontronite-Mg Nontronite-Na 
                         O2(g) Portlandite Pyridine 
                         Pyrite 
                         #Rhodochrosite 
                         Saponite-Fe-Ca Saponite-Fe-Fe #Saponite-Fe-K 
                         Saponite-Fe-Mg Saponite-Fe-Na Saponite-Mg-Ca Saponite-Mg-Fe 
                         #Saponite-Mg-K 
                         Saponite-Mg-Mg Saponite-Mg-Na Sepiolite 
                         Siderite SiO2(am) #Smectite-high-Fe-Mg Smectite-low-Fe-Mg 
                         SO2(g) #Stilbite 
                         #Strengite #Sylvite 
                         Thenardite Thermonatrite 
                         Tobermorite-11A
   -gases                CO2(g) H2(g) CH4(g) CO(g) H2O(g) H2S(g) N2(g) 
                         NH3(g) NO(g) NO2(g) O2(g) SO2(g)  

    
USER_PUNCH
   -headings Time_Years Forst_remain_g Fayal_remain_g Enst_remain_g Ferros_remain_g Pyrr_remain_g Anorth_remain_g Alb_remain_g Mag_remain_g Total_Initial_g
   10 PUNCH STEP_NO * {1/self$phreeqc_params$rock_per_year/adjusted_cycles}  # Time in years
   20 PUNCH KIN("Forsterite") * 140.6715    # Remaining Forsterite mass (g)
   30 PUNCH KIN("Fayalite") * 203.7555      # Remaining Fayalite mass (g)
   40 PUNCH KIN("Enstatite") * 100.3725     # Remaining Enstatite mass (g)
   50 PUNCH KIN("Ferrosilite") * 131.9145   # Remaining Ferrosilite mass (g)
   60 PUNCH KIN("Pyrrhotite") * 87.913      # Remaining Pyrrhotite mass (g)
   70 PUNCH KIN("Anorthite") * 278.164      # Remaining Anorthite mass (g)
   80 PUNCH KIN("Albite") * 262.1798        # Remaining Albite mass (g)
   90 PUNCH KIN("Magnetite") * 231.517      # Remaining Magnetite mass (g)
   100 total_initial = {forsterite_moles}*140.6715 + {fayalite_moles}*203.7555 + {enstatite_moles}*100.3725 + {ferrosilite_moles}*131.9145 + {pyrrhotite_moles}*87.913 + {anorthite_moles}*278.164 + {albite_moles}*262.1798 + {magnetite_moles}*231.517
   110 PUNCH total_initial  # Total initial primary mineral mass (g)


END
    ')
                                     
                                     return(template)
                                   }
                                 )
)

