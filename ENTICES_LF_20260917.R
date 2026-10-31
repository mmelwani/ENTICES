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
# convention as the ENTICES_LF_* sourcing in ENTICES_runs.Rmd.
source("iom_module.R")
# k_primary_mineral() -- reused by generate_phreeqc_input_bulk_kinetic's
# rec_eqb_minerals half-life check, so the auto-equilibrium decision uses the
# exact same rate laws as the RATES blocks below (verified consistent, see
# CHANGELOG 2026-10 Tephroite/Albite/K-Feldspar half-life entries).
source("mineral_dissolution_constants.r")
# phreeqc_tools.R -- PHREEQC_EXPECTED_VERSION (stamped into every .pqi header by
# stamp_pqi_header()) and run_phreeqc(), which warns if the executable differs.
source("phreeqc_tools.R")
# organic_contribution.R -- pair_organic_run(): a run and its organic_wt_percent = 0 blank, and the
# organics' contribution to ocean DIC, CH4, NH4+ and H2 by difference.
source("organic_contribution.R")
# porosity_evolution.R -- run_tidal_evolution(): tidal_kinetic in chunks, with the porosity (and k, alpha,
# skin depth, exchange) updated between chunks from the solid volume (PI 2026-10-29).
source("porosity_evolution.R")

# Primary rock composition, mol per unit total_rock (i.e. per kg anhydrous rock).
# Single source for all four generators that build primary-mineral inventories.
# Source: reconstruction/rock_mol_per_kg.csv (written by
# reconstruction/Primordial_CI_chondrite_mineralogy.Rmd: Scenario R, Lodders et
# al. 2025 Table 4, schreibersite included), rounded to 6 decimals. Regenerated
# 2026-10-23 from the CSV; the previous values came from 3-decimal wt% and were
# off by up to 3e-4 relative. tests/rock_table_check.R asserts the match.
# Deliberately NOT renormalised: magnesiochromite (0.037102 mol/kg) is treated
# as inert, so it keeps its share of the mass but is not added to PHREEQC.
# Lawrencite (FeCl2) is intentional (Cl in the rock); it is not part of the
# reconstruction, and its Fe is taken out of the metal: fe = 2.361630 - 0.011250.
# Total Fe (fe + lawrencite + troilite + 3 x schreibersite) is 4.8859 mmol/g.
# Schreibersite (Fe3P) is listed here but is wired in separately; see the
# generators.
ROCK_MOL_PER_KG <- c(
  enstatite     = 2.275669,
  forsterite    = 1.606976,
  troilite      = 2.383003,
  albite        = 0.318204,
  diopside      = 0.274515,
  anorthite     = 0.062135,
  kfeldspar     = 0.020521,
  tephroite     = 0.025987,
  fe            = 2.350380,
  ni            = 0.280938,
  lawrencite    = 0.011250,
  schreibersite = 0.047093
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
  sprintf("   Schreibersite\n      -formula Fe 3 P 1\n      -m0 %.15e\n      -step_divide 1", m0)
}

# USER_PUNCH columns every generator ends with. Fe metal is punched as Fe_metal / d_Fe_metal rather
# than listed in -equilibrium_phases, where its column was headed "Fe", the same as the -totals Fe
# (dissolved iron) column. The phase is still called Fe in the databases and in EQUILIBRIUM_PHASES
# blocks. Rho_soln is the solution density PHREEQC reports (RHO, g/cm3), the fluid term of
# bulk_density() (phreeqc_tools.R). H_total_mol / O_total_mol are the solution's total H and O in
# moles, water included (TOT("H") x kg water); the -totals H and O columns are the H+ and H2O master
# species' totals, which leave the water out.
SHARED_PUNCH_HEADINGS <- " Fe_metal d_Fe_metal Rho_soln H_total_mol O_total_mol"
SHARED_PUNCH_LINES <- '   900 PUNCH EQUI("Fe")\n   910 PUNCH EQUI_DELTA("Fe")\n   920 PUNCH RHO\n   925 PUNCH TOT("H") * TOT("water"), TOT("O") * TOT("water")\n'
# Optional trace carbon in the initial SOLUTIONs of the kinetic generators (trace_C_molal, mol/kgw;
# NULL = pure water). A numerical seed, not a composition: the science case is pure water, so that
# all ocean carbon can be attributed to the organics. Returns the SOLUTION line, with its leading
# newline, to append to the block's last line ("" for NULL).
# KNOBS around every ocean step (MIX 4) of the tidal generators. At 25 C the ocean step fails to converge
# for scattered porewater fractions ("Numerical method failed on all combinations": 3e-8, 2e-6 and 5e-6
# fail, neighbouring fractions run; 2026-10-29). -diagonal_scale true converges all fractions tested
# (1e-8 to 1e-3) and leaves converging runs unchanged (a 160 C orbit agrees to 1e-8 relative). It is
# switched off after the ocean step: set globally it shifted the porewater step's m_H2S by 0.9%.
# Kozeny-Carman permeability (m2) of a granular medium of grain diameter d (m) and porosity phi (item 3, PI
# 2026-10-29): k = d^2 phi^3 / (180 (1 - phi)^2). The 180 is Carman's (1937) value for packed spheres (36 x
# his constant 5, with the specific surface 6/d); Carman, P.C. (1937), Fluid flow through granular beds,
# Trans. Inst. Chem. Eng. 15, 150-166. A lower bound where fractures carry the flow.
kozeny_carman_k <- function(d, phi) d^2 * phi^3 / (180 * (1 - phi)^2)
# The grain diameter that gives permeability k at porosity phi (the inverse).
kozeny_carman_d <- function(k, phi) sqrt(k * 180 * (1 - phi)^2 / phi^3)
# Geometric specific surface area (m2/g) of spheres of diameter d (m) and density rho (g/cm3): 6 / (rho d).
geometric_ssa_m2_g <- function(d, rho) 6 / (rho * 1e6 * d)

OCEAN_STEP_KNOBS_ON  <- "KNOBS\n    -diagonal_scale true\n"
OCEAN_STEP_KNOBS_OFF <- "KNOBS\n    -diagonal_scale false\n"

trace_C_solution_line <- function(trace_C_molal) {
  if (is.null(trace_C_molal)) return("")
  if (!is.numeric(trace_C_molal) || length(trace_C_molal) != 1 || !isTRUE(trace_C_molal > 0))
    stop("trace_C_molal must be NULL or one positive number (mol/kgw)")
  sprintf("\n    C         %.6e", trace_C_molal)
}
trace_C_header_line <- function(trace_C_molal) {
  sprintf("# trace_C_molal: %s", if (is.null(trace_C_molal)) "NULL (no initial C)" else
    sprintf("%.6e mol/kgw C added to every initial SOLUTION", trace_C_molal))
}
# KNOBS -logfile (the generators' logfile argument). PHREEQC's phreeqc.log records every solver
# iteration: 258-346 MB per 12-orbit production run at 150-160 C, and 11.9 GB from one stalled run.
# The templates write "-logfile false"; logfile = TRUE switches it on (for convergence debugging).
knobs_logfile <- function(content, logfile = FALSE) {
  if (!is.logical(logfile) || length(logfile) != 1 || is.na(logfile)) stop("logfile must be TRUE or FALSE")
  if (!logfile) return(content)
  gsub("\n   -logfile false\n", "\n   -logfile true\n", content)
}

# A restart's record of its parent's initial fluid: the trace_C_molal and initial_ocean header lines
# of the parent .pqi (prev_csv with .csv -> .pqi), re-emitted as "# parent <key>: <value>". A parent
# that is itself a restart carries "# parent <key>:", so the original value passes down a chain.
parent_fluid_header_lines <- function(prev_csv) {
  keys <- c("trace_C_molal", "initial_ocean")
  pqi <- sub("\\.csv$", ".pqi", prev_csv)
  h <- if (file.exists(pqi)) readLines(pqi, n = 200, warn = FALSE) else character(0)
  vals <- vapply(keys, function(k) {
    l <- grep(sprintf("^# (parent )?%s: ", k), h, value = TRUE)
    if (!file.exists(pqi)) sprintf("unknown (parent input not found: %s)", pqi)
    else if (!length(l)) "not stamped in the parent header"
    else sub(sprintf("^# (parent )?%s: ", k), "", l[1])
  }, character(1))
  paste(sprintf("# parent %s: %s", keys, vals), collapse = "\n")
}

# Every element a database defines: its element-level SOLUTION_MASTER_SPECIES, pseudo-elements
# (Mtg, Ntg, Hdg, Oxg, Sg) included, E and Alkalinity not. The three coreclath databases define the
# same 36. Every generator punches all of them (totals_block()), so a restart can carry every element
# whatever the run contained; before 2026-10-28 the totals were a fixed list without K, Mn, Cl, Ni
# and (in the equilibrium generators) Mtg, so restarts dropped those elements.
database_elements <- function(database) {
  g <- names(element_gfw(database))
  g[grepl("^[A-Z][a-z]*$", g) & !g %in% c("E", "Alkalinity")]
}
# The SELECTED_OUTPUT -totals option listing database_elements(database).
totals_block <- function(database) {
  el <- database_elements(database)
  rows <- vapply(split(el, ceiling(seq_along(el) / 16)), paste, "", collapse = " ")
  paste0("    -totals              ", paste(rows, collapse = "\n                         "))
}
# Elements that ENTICES inputs can contain (rock, IOM, initial fluids): a restart warns if its parent
# CSV has no column for one of them (a parent generated before 2026-10-28), since it cannot carry it.
ENTICES_ELEMENTS <- c("Al", "C", "Ca", "Cl", "Fe", "K", "Mg", "Mn", "N", "Na", "Ni", "P", "S", "Si", "Mtg", "Ntg")
# The elements a restart sets in its SOLUTIONs: every database element but H and O (which come from
# the water mass, pH and pe) that the parent CSV has a column for.
restart_elements <- function(database, csv_names, prev_csv) {
  el <- setdiff(database_elements(database), c("H", "O"))
  miss <- setdiff(intersect(ENTICES_ELEMENTS, el), csv_names)
  if (length(miss))
    warning("restart: the parent CSV has no total for ", paste(miss, collapse = ", "), " (", prev_csv,
            "; generated before 2026-10-28): these elements are not carried into the restart")
  intersect(el, csv_names)
}

# A parent run's carbon mode, from the database its .pqi header names ("# Database:", else
# "# Recommended database:"). NA when there is no parent .pqi or it names none of these.
CARBON_MODE_DATABASES <- c(coupled = "coreclath_CH4coupled.dat", ch4_uncoupled = "coreclath_CH4uncoupled.dat",
                           c4_only = "coreclath_C4only.dat")
parent_carbon_mode <- function(prev_csv) {
  pqi <- sub("\\.csv$", ".pqi", prev_csv)
  if (!file.exists(pqi)) return(NA_character_)
  h <- readLines(pqi, n = 200, warn = FALSE)
  for (l in c(grep("^# Database: ", h, value = TRUE), grep("^# Recommended database: ", h, value = TRUE))) {
    hit <- names(CARBON_MODE_DATABASES)[vapply(CARBON_MODE_DATABASES, grepl, logical(1), x = l, fixed = TRUE)]
    if (length(hit)) return(hit[1])
  }
  NA_character_
}
# A restart's carbon mode against its parent's. A mismatch is an error unless the caller chose the
# mode (carbon_mode or CH4_redox_override given: ch4_mode$explicit), and then a warning. A parent
# whose mode can't be read is a warning. Returns the "# parent carbon_mode:" header line.
check_restart_carbon_mode <- function(prev_csv, ch4_mode) {
  parent <- parent_carbon_mode(prev_csv)
  if (is.na(parent)) {
    warning("restart: the parent's carbon mode is unknown (no parent .pqi naming a known database: ",
            sub("\\.csv$", ".pqi", prev_csv), "); not checked")
    return("# parent carbon_mode: unknown (not checked)")
  }
  if (parent == ch4_mode$carbon_mode) return(sprintf("# parent carbon_mode: %s", parent))
  msg <- sprintf('restart: the parent ran carbon_mode "%s", this restart would run "%s"', parent, ch4_mode$carbon_mode)
  if (!ch4_mode$explicit)
    stop(msg, sprintf('. Give carbon_mode = "%s" to continue the parent, or carbon_mode = "%s" to switch.',
                      parent, ch4_mode$carbon_mode))
  warning(msg, "; the parent's solutions are re-speciated in the restart's mode")
  sprintf("# parent carbon_mode: %s (this restart switches to %s)", parent, ch4_mode$carbon_mode)
}

# Initial ocean composition for the tidal generators (applied to the ocean SOLUTION 2 only).
# "pure_water" (default) or a named list of element totals in mol/kgw, e.g. list(C = ..., N = ...),
# for a cometary-ice scenario; names are PHREEQC elements or valence states such as "N(-3)".
# charge_balance: NULL (default, no charge balance), "pH", or one of the initial_ocean names, whose
# total PHREEQC then adjusts to balance SOLUTION 2's charge.
# Returns the SOLUTION lines, each with its leading newline ("" for "pure_water").
initial_ocean_solution_lines <- function(initial_ocean, trace_C_molal = NULL, charge_balance = NULL) {
  if (inherits(initial_ocean, "cometary_ice")) {           # SOLUTION 2 is pure water; the ices enter by REACTION
    if (!is.null(charge_balance)) stop("charge_balance is not used with cometary_ice_ocean()")
    return("")
  }
  check_ocean_charge_balance(initial_ocean, charge_balance)
  if (identical(initial_ocean, "pure_water")) return("")
  nm <- names(initial_ocean)
  if (!is.list(initial_ocean) || length(initial_ocean) == 0 || is.null(nm) || any(!nzchar(nm)) || anyDuplicated(nm))
    stop('initial_ocean must be "pure_water" or a named list of element totals (mol/kgw)')
  bad <- nm[!grepl("^[A-Z][a-z]?(\\([+-]?[0-9]+\\))?$", nm)]
  if (length(bad)) stop("initial_ocean: not an element name: ", paste(bad, collapse = ", "))
  ok <- vapply(initial_ocean, function(v) is.numeric(v) && length(v) == 1 && isTRUE(v > 0), logical(1))
  if (!all(ok)) stop("initial_ocean: each total must be one positive number (mol/kgw): ", paste(nm[!ok], collapse = ", "))
  if (!is.null(trace_C_molal) && any(grepl("^C(\\(|$)", nm)))
    stop("trace_C_molal and initial_ocean both set the ocean's C; give only one")
  paste0(sprintf("\n    %-9s %.6e%s", nm, unlist(initial_ocean),
                 ifelse(nm %in% charge_balance, " charge", "")), collapse = "")
}
check_ocean_charge_balance <- function(initial_ocean, charge_balance) {
  if (is.null(charge_balance) || identical(charge_balance, "pH")) return(invisible(TRUE))
  ok <- is.character(charge_balance) && length(charge_balance) == 1 && is.list(initial_ocean) &&
    charge_balance %in% names(initial_ocean)
  if (!ok) stop('charge_balance must be NULL, "pH", or the name of one of the initial_ocean totals')
  invisible(TRUE)
}
# " charge" for SOLUTION 2's pH line when charge_balance = "pH", else "".
ocean_pH_charge <- function(charge_balance) if (identical(charge_balance, "pH")) " charge" else ""
initial_ocean_header_line <- function(initial_ocean, charge_balance = NULL) {
  cb <- if (is.null(charge_balance)) "" else sprintf("; charge balance on %s", charge_balance)
  if (identical(initial_ocean, "pure_water")) return(paste0("# initial_ocean: pure_water", cb))
  if (inherits(initial_ocean, "cometary_ice"))
    return(sprintf("# initial_ocean: cometary_ice_ocean(v = %g, include_NH3 = %s, include_O2 = %s): %s mol/kgw, reacted into the pure-water ocean at ocean T and P (Miller et al. 2026 Table 5; Rubin et al. 2019)",
                   initial_ocean$v, initial_ocean$include_NH3, initial_ocean$include_O2,
                   paste(sprintf("%s=%.4e", names(initial_ocean$mol_per_kgw), initial_ocean$mol_per_kgw), collapse = ", ")))
  sprintf("# initial_ocean: %s mol/kgw in the ocean SOLUTION 2 only%s",
          paste(sprintf("%s=%.6e", names(initial_ocean), unlist(initial_ocean)), collapse = ", "), cb)
}

# Header note written into every tidal_equil / bulk_equil input.
EQUIL_INORGANIC_NOTE <- paste(
  "# Inorganic limiting case: no IOM and no Schreibersite. organic_wt_percent only scales the",
  "# mineral mass by (100 - organic_wt_percent)/100; no organic carbon, N or S is added.", sep = "\n")

# Index of the output column holding an equilibrium phase's moles. An index, not a name: with
# check.names = FALSE a pre-Fe_metal CSV has two columns named "Fe", and [["Fe"]] takes the first.
# Those older CSVs carry Fe metal in the second "Fe" column, used here so restarts from them still work.
ep_column <- function(phase, col_names) {
  if (phase == "Fe") {
    if ("Fe_metal" %in% col_names) return(match("Fe_metal", col_names))
    hits <- which(col_names == "Fe" | grepl("^Fe\\.[0-9]+$", col_names))
    return(if (length(hits) > 1) tail(hits, 1) else NA_integer_)
  }
  # A phase named like an element (Ni) shares its name with that element's -totals column, which comes
  # first; read.delim's default check.names renames the second "Ni.1". Either way, take the last.
  mn <- make.names(col_names); mp <- make.names(phase)
  hits <- which(mn == mp | startsWith(mn, paste0(mp, ".")) & grepl("^\\.[0-9]+$", substring(mn, nchar(mp) + 1)))
  if (length(hits) > 0) tail(hits, 1) else NA_integer_
}

#===============================================================================
# Gas phases (2026-10-28)
#===============================================================================
# Gases are a fixed-pressure GAS_PHASE in each reservoir, not EQUILIBRIUM_PHASES. An "X(g) 0 0"
# equilibrium phase is an infinite gas reservoir at 1 atm: it capped dissolved gases at their 1-atm
# solubility and kept whatever exsolved, whatever the pressure. A fixed-pressure gas phase (Parkhurst
# & Appelo 2013, GAS_PHASE) forms a bubble only when the partial pressures sum to more than the
# reservoir's pressure, and holds a finite amount. Every gas is Peng-Robinson (Peng & Robinson 1976;
# critical constants from scripts/gas_eos.py), which PHREEQC requires of all components of a phase.
#
# GAS_PHASE 1: porewater, at pore_pressure_atm; GAS_PHASE 2: ocean, at ocean_pressure_atm (tidal
# generators; never emptied: the PI's default, plume loss being unknown). Both start empty. The bulk
# generators have one, GAS_PHASE 1, at their solution's pressure. pore_gas_fate (tidal; default "to_ocean",
# PI 2026-10-29):
#   "retain"   -- the pore bubble stays with the porewater and is carried to the next orbit (MIX
#                 moves solution only). When the bubble is larger than the pore space (Vgas_Vpore > 1)
#                 this is an upper bound on retention, not a physical state.
#   "to_ocean" -- after each porewater step the bubble's moles go to the ocean step (see
#                 gas_transfer_block()) and the porewater starts the next orbit with an empty phase.

# The components: the gases the EQUILIBRIUM_PHASES lists carried before 2026-10-28, with methane as
# Mtg(g) when it is decoupled, and only those the carbon mode's database defines (c4_only has no CO(g)
# or C2-C3 gases). extra: more components (Ntg(g) for an initial ocean holding Ntg).
GAS_COMPONENTS <- c("C2H4(g)", "C2H6(g)", "C3H8(g)", "CH4(g)", "CO(g)", "CO2(g)", "H2(g)", "H2S(g)", "N2(g)",
                    "NH3(g)", "SO2(g)")
gas_components <- function(ch4_mode, extra = character(0)) {
  g <- c(GAS_COMPONENTS, extra)
  if (ch4_mode$decoupled) g[g == "CH4(g)"] <- "Mtg(g)"
  unique(intersect(g, database_names(ch4_mode$database)$phases))
}
# An empty fixed-pressure GAS_PHASE block.
gas_phase_block <- function(n, comps, pressure_atm, temp_C, label) {
  paste(c(sprintf("GAS_PHASE %d %s", n, label),
          "    -fixed_pressure",
          sprintf("    -pressure    %.1f", pressure_atm),
          "    -volume      1e-6      # negligible: every component starts at partial pressure 0, so the phase is empty",
          sprintf("    -temperature %s", format(round(temp_C, 2))),
          sprintf("    %-10s 0", comps)), collapse = "\n")
}
# GAS_PHASE_MODIFY setting a restart's gas phase n to the moles in a parent row (its g_ columns).
gas_restore_block <- function(n, comps, row, scale = 1) {
  # by make.names(), so a row read with or without check.names matches ("g_CO2(g)" / "g_CO2.g.")
  v <- vapply(comps, function(c) { k <- match(make.names(paste0("g_", c)), make.names(names(row)))
    if (!is.na(k)) suppressWarnings(as.numeric(row[[k]])) else NA_real_ }, 0)
  v <- v[!is.na(v) & v > 0] * scale
  if (!length(v)) return(sprintf("# gas phase %d: empty in the parent", n))
  paste(c(sprintf("GAS_PHASE_MODIFY %d   # the parent's gas phase %d", n, n),
          unlist(lapply(names(v), function(c) c(sprintf("    -component %s", c), sprintf("        -moles %.15e", v[[c]]))))),
        collapse = "\n")
}
# SELECTED_OUTPUT -gases: the component moles (g_ columns), with the phase's pressure, total moles
# and volume.
gas_so_line <- function(comps) sprintf("    -gases               %s", paste(comps, collapse = " "))
# pore_gas_fate = "to_ocean". PHREEQC has no keyword that adds one gas phase's moles to another, so
# the transfer goes through the ocean step's solution: the porewater step's USER_PUNCH stores each
# component's moles with PUT (gas_punch(put = TRUE)), and KINETICS 3, used in the ocean step, adds
# them to the mixed ocean solution at a constant rate over 1 s (a constant rate integrates exactly),
# after which they partition between the ocean solution and GAS_PHASE 2 at the ocean's pressure.
# GAS_PHASE 1 is then redefined empty at the start of the next porewater step. Mass is conserved to
# the integrator's precision (tested: < 1e-11 relative); there is never more than one gas reservoir
# per reacting solution. Each component enters by its own formula (CO2, H2, Mtg, ...).
gas_transfer_block <- function(comps) {
  nm <- sprintf("PoreGasIn_%d", seq_along(comps))
  paste(c("RATES",
          unlist(lapply(seq_along(comps), function(j) c(nm[j], " -start",
                        sprintf(" 10 SAVE GET(%d) * TIME   # moles of %s in the pore bubble", 9100 + j, comps[j]), " -end"))),
          "",
          "KINETICS 3   # used only in the ocean step, when pore_gas_fate = to_ocean",
          unlist(lapply(seq_along(comps), function(j) c(sprintf("   %s", nm[j]), sprintf("      -formula %s 1", sub("\\(g\\)$", "", comps[j])),
                                                         "      -m0 1e6"))),
          "   -steps 1"), collapse = "\n")
}
# USER_PUNCH columns for the gas phase of the simulation in use: its volume (L), pressure (atm; 0 when
# there is no bubble) and volume over the pore volume (vpore_L, from the same normalisation as the
# porewater mass, at 1 kg/L). On ocean rows Vgas_Vpore uses the same pore volume. put = TRUE adds the
# PUT lines that store the pore bubble for the to_ocean transfer.
gas_punch <- function(comps, vpore_L, put = FALSE) {
  sum_lines <- sprintf("   %d ng = ng + GAS(\"%s\")\n", 931 + seq_along(comps) - 1, comps)
  put_lines <- if (put) sprintf("   %d PUT(GAS(\"%s\"), %d)\n", 970 + seq_along(comps) - 1, comps, 9100 + seq_along(comps)) else character(0)
  list(headings = " Vgas_L P_gas_atm Vgas_Vpore",
       lines = paste0(c("   930 ng = 0\n", sum_lines, "   950 vg = 0\n", "   951 IF ng > 0 THEN vg = GAS_VM * ng\n",
                        sprintf("   960 PUNCH vg, GAS_P, vg / %.15g\n", vpore_L), put_lines), collapse = ""))
}

# Cometary-ice initial ocean (an option; the default ocean is pure water). Miller et al. (2026),
# Icarus 449, 116961, Table 5, "molar abundance used (volatile depleted)": comet 67P ice from ROSINA
# (Rubin et al. 2019, MNRAS 489, 594-607, the "67P" column below, mol per 100 mol H2O), with every
# ice more volatile than NH3 at 273 K (Fray & Schmitt 2009) depleted by v_NH3 = 0.02. v_NH3 comes
# from Titan's atmospheric N budget (Miller et al. 2026, eqs. 21-24): it is not an Enceladus
# constraint. NH3 (0.67 per 100 H2O, Rubin et al. 2019 Table 4) is not in Miller's Table 5, whose
# "Total N" (2 x N2) excludes it, while their text says v_NH3 applies to NH3 too: include_NH3
# (default FALSE, Table 5 as printed) is pending the PI. include_O2 exists because O2 decides the
# redox outcome. v = 1 is undepleted 67P ice: an unrealistic upper bound (2.6 molal CO2 + 1.7 molal
# CO, far above CO2-clathrate dissociation at 1.45 C).
COMET_67P_PER_100_H2O <- c(CO2 = 4.7, CO = 3.1, O2 = 3.1, H2S = 1.1, CH4 = 0.34, C2H6 = 0.29, N2 = 8.9e-2, NH3 = 0.67)
cometary_ice_ocean <- function(v = 0.02, include_NH3 = FALSE, include_O2 = TRUE) {
  if (!is.numeric(v) || length(v) != 1 || !isTRUE(v > 0 && v <= 1)) stop("v must be one number in (0, 1]")
  sp <- setdiff(names(COMET_67P_PER_100_H2O), c(if (!include_NH3) "NH3", if (!include_O2) "O2"))
  # mol per 100 mol H2O x v -> mol/kgw: / 100 / 0.0180153 kg per mol H2O
  structure(list(mol_per_kgw = COMET_67P_PER_100_H2O[sp] * v / 100 / 0.0180153, v = v,
                 include_NH3 = include_NH3, include_O2 = include_O2), class = "cometary_ice")
}
# REACTION lines for the ices, by carbon mode: CO2, CO, O2, H2S and NH3 are coupled (equilibrium
# decides CO's fate; under c4_only C can only be +4, so CO + H2O -> CO2 + H2);
# CH4 is Mtg and C2H6 is "Mtg 2 H -2" (the electron-count rule, as IOM_CHn) when methane is
# decoupled, CH4 and C2H6 when coupled; N2 is Ntg (kinetically inert at ocean temperature) in every
# mode. Amounts are mol/kgw x the ocean's water mass (kg). "Mtg 2 H -2" reduces ethane's C from -3
# to -4, drawing 2 e- per C2H6 from the coupled system (0.0064 eq at v = 0.02): it oxidises 0.0008
# molal more H2S to sulfate when O2 is present (tests/cometary_ocean_check.R).
cometary_reaction_lines <- function(ice, ch4_mode, water_kg) {
  m <- ice$mol_per_kgw * water_kg; g <- function(n) if (n %in% names(m)) m[[n]] else 0
  dec <- ch4_mode$decoupled
  # CO enters as CO2 minus O: REACTION matches phase names first, case-insensitively, and "CO" is
  # the cobalt phase Co (tested 2026-10-28: it added 0.034 molal Co).
  r <- c(CO2 = g("CO2") + g("CO"), O = -g("CO"), O2 = g("O2"), H2S = g("H2S"), Ntg = g("N2"), NH3 = g("NH3"),
         Mtg = if (dec) g("CH4") + 2 * g("C2H6") else 0, H = if (dec) -2 * g("C2H6") else 0,
         CH4 = if (dec) 0 else g("CH4"), C2H6 = if (dec) 0 else g("C2H6"))
  r <- r[r != 0]
  c(sprintf("    %-5s %.15e", names(r), r), "    1 mol")
}
# The simulation that reacts the ices into the pure-water ocean SOLUTION 2 at ocean temperature and
# pressure, with the ocean's gas phase (which it saves, so a bubble would be kept), and writes the
# resulting ocean to <base>_initial_ocean.csv (SELECTED_OUTPUT 2, switched off afterwards).
cometary_ocean_simulation <- function(ice, ch4_mode, water_kg, ocean_gas_block, base_name, database) {
  comps <- gas_components(ch4_mode, "Ntg(g)")
  paste(c("END",
          sprintf("# Cometary-ice ocean: cometary_ice_ocean(v = %g, include_NH3 = %s, include_O2 = %s), in mol/kgw: %s.",
                  ice$v, ice$include_NH3, ice$include_O2,
                  paste(sprintf("%s %.4e", names(ice$mol_per_kgw), ice$mol_per_kgw), collapse = ", ")),
          "# Entered as a REACTION into pure water at ocean T and P, not as SOLUTION totals: the ice mixture",
          "# is out of redox equilibrium (O2 with CO and H2S), and equilibrium decides the outcome. That is",
          "# an assumption at 1.45 C, where the CO shift and O2 oxidation are kinetically slow.",
          ocean_gas_block,
          "USE solution 2",
          "REACTION 9 Cometary ices", cometary_reaction_lines(ice, ch4_mode, water_kg),
          "USE reaction_temperature 2", "USE reaction_pressure 2",
          "SAVE solution 2", "SAVE gas_phase 2",
          "SELECTED_OUTPUT 2",
          sprintf("    -file %s_initial_ocean.csv", base_name),
          "    -high_precision true", "    -state true", "    -solution true", "    -reaction true", "    -ph true", "    -pe true", "    -alkalinity true", "    -water true",
          totals_block(database),
          "    -molalities CO2 HCO3- CO3-2 Mtg H2 O2 HS- H2S SO4-2 HSO4- Ntg NH3 NH4+ HCOO- CO CH3COO- C2H6",
          gas_so_line(comps),
          "USER_PUNCH 2", "    -headings H_total_mol O_total_mol Vgas_L", "    10 PUNCH TOT(\"H\") * TOT(\"water\"), TOT(\"O\") * TOT(\"water\")",
          "    20 ng = 0", sprintf("    %d ng = ng + GAS(\"%s\")", 20 + seq_along(comps), comps), "    90 PUNCH GAS_VM * ng",
          "END",
          "SELECTED_OUTPUT 2", "    -active false"), collapse = "\n")
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

# Porewater and ocean conditions from the PlanetProfile profile, used by every generator. There is
# no fallback: without a profile there is no defensible pressure. (Before 2026-10-27 the generators
# fell back to 70 atm, and five of the six wrote the ocean-floor pressure, 66.7 atm for the
# Enceladus profile, into the porewater SOLUTION.)
#   pore_temp_C        -- mean temperature of the permeable layer (porewater_temp = "profile")
#   pore_pressure_atm  -- mean pressure of the permeable layer: the porewater SOLUTION, and the
#                         reaction pressure of every porewater step
#   ocean_pressure_atm -- pressure at the ocean floor (top of phase 50): the ocean SOLUTION
MPA_PER_ATM <- 0.101325
profile_conditions <- function(profile_file, layer_thickness) {
  if (is.null(profile_file) || !nzchar(profile_file) || !file.exists(profile_file))
    stop(sprintf(paste0("profile_file not found: '%s'. Set PhreeqcIntegrator$profile_file to a PlanetProfile ",
                        "*Profile*.txt (the Enceladus profile is in profiles/). Pressures come only from ",
                        "the profile; there is no default."),
                 if (is.null(profile_file)) "NULL" else profile_file))
  list(pore_temp_C        = mean_layer_temp_C(profile_file, layer_thickness),
       pore_pressure_atm  = mean_layer_pressure_MPa(profile_file, layer_thickness) / MPA_PER_ATM,
       ocean_pressure_atm = ocean_floor_pressure_MPa(profile_file) / MPA_PER_ATM)
}

# The two pressure lines stamp_pqi_header() writes, and parent_pressures_atm() reads back.
pressure_header_lines <- function(pore_pressure_atm, ocean_pressure_atm) {
  c(sprintf("# Porewater pressure: %.1f atm (profile mean over the reacting box)", pore_pressure_atm),
    sprintf("# Ocean pressure: %.1f atm (profile, ocean floor)", ocean_pressure_atm))
}
# A parent run's pressures (atm), from its .pqi header. NA where there is no parent .pqi or no
# stamp (inputs generated before 2026-10-27).
parent_pressures_atm <- function(prev_csv) {
  pqi <- sub("\\.csv$", ".pqi", prev_csv)
  h <- if (file.exists(pqi)) readLines(pqi, n = 200, warn = FALSE) else character(0)
  get <- function(rx) {
    l <- grep(rx, h, value = TRUE)
    if (length(l)) as.numeric(sub(".*pressure: ([0-9.eE+-]+) atm.*", "\\1", l[1])) else NA_real_
  }
  c(pore = get("^# (parent )?Porewater pressure: "), ocean = get("^# (parent )?Ocean pressure: "))
}
# A restart's pressures against its parent's, compared at the 0.1 atm written into the inputs. Like
# check_restart_carbon_mode(): a mismatch is an error, or a warning with allow_pressure_change =
# TRUE; a parent whose pressures can't be read is a warning. Returns the "# parent pressures:" line.
check_restart_pressures <- function(prev_csv, pore_pressure_atm, ocean_pressure_atm, allow_pressure_change = FALSE) {
  parent <- parent_pressures_atm(prev_csv)
  if (anyNA(parent)) {
    warning("restart: the parent's pressures are not stamped in ", sub("\\.csv$", ".pqi", prev_csv),
            " (inputs generated before 2026-10-27, or no parent .pqi); not checked")
    return("# parent pressures: unknown (not checked)")
  }
  line <- sprintf("# parent pressures: porewater %.1f atm, ocean %.1f atm", parent[["pore"]], parent[["ocean"]])
  same <- sprintf("%.1f", parent) == sprintf("%.1f", c(pore_pressure_atm, ocean_pressure_atm))
  if (all(same)) return(line)
  msg <- sprintf("restart: the parent ran at porewater %.1f / ocean %.1f atm, this restart would run at %.1f / %.1f atm (from profile_file)",
                 parent[["pore"]], parent[["ocean"]], pore_pressure_atm, ocean_pressure_atm)
  if (!allow_pressure_change)
    stop(msg, ". Set profile_file to the parent's profile, or give allow_pressure_change = TRUE to switch.")
  warning(msg, "; continuing (allow_pressure_change = TRUE)")
  paste0(line, sprintf(" (this restart switches to %.1f / %.1f atm)", pore_pressure_atm, ocean_pressure_atm))
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
                                 # radius, mass, core_radius and core_density are the PlanetProfile profile's, rounded (profiles/: R_m 252100,
                                 # M_kg 1.08022e20, RsilMean_m 194243.7, rhoSil_kgm3 2353.19; Styczinski et al. 2023), a structure fitted
                                 # to C/MR2 = 0.335 (Iess et al. 2014). Neither Iess et al. 2014 nor Bagheri et al. 2025 gives these
                                 # values; both bracket them (Bagheri Table 5: core 2160-2400 kg/m3, 183-213 km). README section 7.
                                 core_density = 2353,     # Core bulk density [kg/m³]
                                 core_porosity = 0.32,    # Core porosity [volume fraction] #Default assumed value from Stycinski et al 2023
                                 
                                 # Orbital properties
                                 orbital_period = 118800, # [s]
                                 # Tidal pore-pressure amplitude [Pa]. NULL (default): derived from the strain, dP = alpha eps_v M
                                 # (Simulator$tidal_pressure(); PI 2026-10-30). A number overrides it (the old default was an
                                 # UNSOURCED 1e5 Pa); README section 6.
                                 max_tidal_stress = NULL,
                                 # Ocean tidal loading at the seafloor (PI 2026-10-31): the whole-body surface Love number h2 of
                                 # Bagheri et al. (2025) Table 3, "dissipative core" (|h2| 0.0459 +- 0.0294, lag 2.95 +- 2.69 deg);
                                 # the ocean surface beneath the shell is taken to move with the surface. ocean_load = FALSE drops it.
                                 # Scale on the core's tidal strain eps0 (PI 2026-10-31; N6 sweep). The isolated core's h2 behind eps0
                                 # (0.29 at the placeholders) is 6x Bagheri et al. (2025)'s whole-body h2 (0.046): the bracket is
                                 # strain_scale = 0.046/0.29 = 0.16 to 1 (README section 6). Default 1.
                                 strain_scale = 1,
                                 surface_h2 = 0.0459,
                                 surface_h2_lag_deg = 2.95,
                                 ocean_load = TRUE,
                                 eccentricity = 0.0047,   # Enceladus's orbital eccentricity (the value with which Liao et al. 2020's Fig. 2a is reproduced)

                                 # Poroelastic exchange (README section 6, Exchange). Core rigidity and viscosity: Bagheri et al.
                                 # (2025) Table 5, "dissipation in the core" (mu 4-70 GPa, eta 4.4e11-2.3e12 Pa s); PLACEHOLDERS at
                                 # the geometric means until the PI picks values. (The Biot coefficient is the Simulator's,
                                 # from porosity: Simulator$biot_alpha().)
                                 core_rigidity = 17e9,    # [Pa]
                                 core_viscosity = 1e12,   # [Pa s]
                                 # Bulk moduli of the pore water and the rock grains, for the poroelastic skin depth:
                                 # Liao et al. (2020) Table A1.
                                 water_bulk_modulus = 2.2e9,  # [Pa] K_f
                                 rock_bulk_modulus = 10e9,    # [Pa] K_s
                                 
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
  # Polar tidal strain amplitude of the core, Liao et al. (2020) section 2.2 and Fig. 2 caption:
  #   eps0 = 9/(4 pi) e w^2/(rho G) (5/3) |k2|,  k2 = (3/2) / (1 + 19 mu* / (2 rho g R)),
  #   mu* = mu i w tau / (1 + i w tau), tau = eta / mu (Maxwell), g = (4/3) pi G rho R.
  # Their section 2.2 text and Appendix A write k2 with an extra factor 3/5; the caption's form is the one
  # that reproduces their Fig. 2a contours (checked 2026-10-29, README section 6), so it is used here.
  tidal_strain_amplitude = function(mu, eta, rho, R, e, period) {
    w <- 2 * pi / period; tau <- eta / mu
    mus <- mu * (1i * w * tau) / (1 + 1i * w * tau)
    g <- 4 / 3 * pi * CONSTANTS$G * rho * R
    9 / (4 * pi) * e * w^2 / (rho * CONSTANTS$G) * 5 / 3 * Mod(1.5 / (1 + 19 * mus / (2 * rho * g * R)))
  },

  # Volume mean of the positive part of a degree-2 dilatation field (r/R)^2 P2(cos theta), amplitude 1 at
  # the pole surface, over the shell R - L to R: the angular mean of max(P2, 0), 0.19245 (P2 > 0 for
  # |cos theta| > 1/sqrt(3)), times the radial mean of (r/R)^2, (3/5)(1 - x^5)/(1 - x^3), x = (R - L)/R.
  # The negative part is equal (a degree-2 field has zero mean): this is the fraction of eps0 * V_layer
  # compressed (expelling pore fluid) while an equal volume expands (drawing it in), each half-orbit.
  degree2_dilatation_factor = function(R, L) {
    x <- max(0, (R - L) / R)
    (1 / sqrt(3) - 1 / (3 * sqrt(3))) / 2 * 3 / 5 * (1 - x^5) / (1 - x^3)
  },

  # q is the Darcy flux (volume per unit TOTAL area per unit time), so the volumetric flow through
  # area A is q * A; porosity enters only in the pore velocity q / porosity. (Not called.)
  calculate_darcy_flow = function(k, grad_P, area,
                                  viscosity = CONSTANTS$FLUID$viscosity) {
    q <- -(k/viscosity) * grad_P  # Specific discharge
    Q <- q * area                 # Volumetric flow rate
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
                           exchange = NULL,
                           config = NULL,
                           
                           # Initialize simulator with parameters and configuration
                           # reactive_depth: "skin" (default) makes the reacting box the top L_eff = min(layer_thickness,
                           # delta) of the layer, delta the poroelastic skin depth (skin_depth()); "layer" uses the whole
                           # layer (the behaviour before 2026-10-29). See reactive_thickness().
                           # phi0: initial porosity (default 0.32, the PlanetProfile profile's phiRockMax: ONE fit, Styczinski
                           # et al. 2023; N6 sweep 0.1-0.4). `porosity` is the older name for it and, if given, wins.
                           # phi_c: critical porosity of the critical-porosity model (Nur et al. 1998), default 0.4 (PI 10-30:
                           # Dvorkin & Nur 2000's table, 36-40% for sorted grain packs and sandstones; still a placeholder);
                           # the Biot coefficient is then min(1, phi / phi_c) (biot_alpha()); phi >= phi_c gives alpha = 1, a
                           # suspension rather than a frame.
                           # alpha: a number overrides that; it must be >= the porosity.
                           # k: the permeability (m2), or "kozeny_carman" (item 3, PI 2026-10-29): k = d^2 phi^3 / (180 (1 - phi)^2)
                           # from grain_diameter d (m) and the porosity (kozeny_carman_k()). ssa: "legacy" (default: the
                           # per-mineral surface areas in the RATES blocks, 0.1 and 5 m2/g, UNSOURCED) or "geometric":
                           # 6 / (rho_i d_i) per mineral (rho_i from data/mineral_molar_volumes.csv), times `roughness`
                           # (default 1; BET over geometric area); d_i from grain_diameter_by_mineral (a named vector, m) or
                           # grain_diameter (applied by PhreeqcIntegrator$apply_ssa() to every generated input).
                           initialize = function(k, layer_thickness, porosity = NULL, config = NULL,
                                                 reactive_depth = c("skin", "layer"), phi0 = 0.32, phi_c = 0.4, alpha = NULL,
                                                 grain_diameter = NULL, ssa = c("legacy", "geometric"), roughness = 1,
                                                 grain_diameter_by_mineral = NULL) {
                             phi_ <- if (is.null(porosity)) phi0 else porosity
                             perm <- if (identical(k, "kozeny_carman")) "kozeny_carman" else "given"
                             if (perm == "kozeny_carman") {
                               if (!(is.numeric(grain_diameter) && length(grain_diameter) == 1 && grain_diameter > 0))
                                 stop("k = \"kozeny_carman\" needs grain_diameter (m)")
                               k <- kozeny_carman_k(grain_diameter, phi_)
                             }
                             ssa <- match.arg(ssa)
                             if (ssa == "geometric" && is.null(grain_diameter) && is.null(grain_diameter_by_mineral))
                               stop("ssa = \"geometric\" needs grain_diameter or grain_diameter_by_mineral")
                             self$params <- list(
                               k = k,
                               permeability = perm,
                               grain_diameter = grain_diameter,
                               ssa = ssa,
                               roughness = roughness,
                               grain_diameter_by_mineral = grain_diameter_by_mineral,
                               layer_thickness = layer_thickness,
                               porosity = phi_,
                               reactive_depth = match.arg(reactive_depth),
                               phi_c = phi_c,
                               alpha = alpha
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
                             if (!(is.numeric(self$params$phi_c) && self$params$phi_c > 0 && self$params$phi_c <= 1))
                               stop("phi_c must be in (0, 1]")
                             if (!is.null(self$params$alpha) && (!is.numeric(self$params$alpha) || self$params$alpha > 1 ||
                                                                 self$params$alpha < self$params$porosity))
                               stop(sprintf("alpha = %s: the Biot coefficient must be >= the porosity (%g) and <= 1",
                                            format(self$params$alpha), self$params$porosity))
                           },
                           
                           # Calculate tidal stress variation
                           calculate_tidal_stress = function(t) {
                             # Debug print
                             print(sprintf("Calculating stress for %d time points", length(t)))
                             print(sprintf("Using orbital period: %f", self$config$orbital_period))
                             print(sprintf("Using tidal pressure amplitude: %f Pa", self$tidal_pressure()))
                             
                             # Calculate phase in orbital cycle
                             phase <- 2 * pi * t / self$config$orbital_period
                             # Calculate stress - ensure vectorized operation
                             stress <- self$tidal_pressure() * sin(phase)
                             
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
                             
                             # The tidal pressure perturbation diffuses only a skin depth delta into the rock, so the
                             # gradient acts over L_eff = min(layer_thickness, delta) (skin_depth(); README section 6).
                             pressure_gradient <- stress / self$skin_depth()$L_eff
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
                                 # flow_velocity is the Darcy flux q (per unit total area) and flow_area the
                                 # core's whole surface, so the flow is q * A. Until 2026-10-29 it was also
                                 # multiplied by porosity, which made it 1/porosity (3.1x) too small.
                                 flow_rate_m3s = abs(flow_velocity) * flow_area,
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
                             self$exchange <- self$calculate_exchange()
                             print("Stored results")
                             print(sprintf("Final results tibble rows: %d", nrow(self$results)))
                             
                             return(results)
                           },
                           
                           # Tidal strain of the reacting box. eps0: Liao et al. (2020)'s polar strain amplitude. eps_v: the
                           # volume mean of the compressed part of a degree-2 dilatation field with pole amplitude eps0 over the
                           # box, eps0 * degree2_dilatation_factor(R, L): (mean |P2| over the sphere = 0.385) / 2 = 0.192 times
                           # the radial mean of (r/R)^2 over the box: 0.192 for a thin skin, 0.115 for the full core (reactive
                           # depth "layer", 194.2 km). PLACEHOLDER pending the PI (README section 6, The strain).
                           tidal_strain = function() {
                             cfg <- self$config; R <- cfg$core_radius; L <- self$reactive_thickness()
                             eps0 <- cfg$strain_scale * FlowUtils$tidal_strain_amplitude(cfg$core_rigidity, cfg$core_viscosity, cfg$core_density,
                                                                                    R, cfg$eccentricity, cfg$orbital_period)
                             f <- FlowUtils$degree2_dilatation_factor(R, L)
                             list(eps0 = eps0, factor = f, eps_v = eps0 * f)
                           },

                           # Tidal pore-pressure amplitude (Pa) driving the Darcy flow: the undrained pressure of the strain,
                           # dP = alpha eps_v M with the Biot modulus M = 1/S (Liao et al. 2020 Table A1: (K_u - K_m)/alpha^2;
                           # their Eq. A5 with zero fluid-content change). config$max_tidal_stress, if a number, overrides it.
                           tidal_pressure = function() {
                             if (!is.null(self$config$max_tidal_stress)) return(self$config$max_tidal_stress)
                             self$tidal_pressure_strain() + self$tidal_load_pressure()$mean
                           },
                           tidal_pressure_strain = function() self$biot_alpha() * self$tidal_strain()$eps_v / self$skin_depth()$S,

                           # Ocean tidal loading at the seafloor (PI 2026-10-31): the change of the ocean column's weight,
                           # rho_w g_c (u_ocean - u_core), with u_r = h2 U2 / g (Bagheri et al. 2025 Eq. 2) and the eccentricity
                           # tide's amplitude at the pole, U2 = (3/2) r^2 omega^2 e (their Eq. 1 with P20 = 1, P22 = 0; omega the
                           # mean motion, e the eccentricity). u_ocean: the surface, |h2| and lag from config (Bagheri Table 3,
                           # dissipative core), at R = config$radius; u_core: the core's own h2 = (5/3) k2* (Liao et al. 2020's
                           # uniform Maxwell body, the h2 behind eps0) at the core radius; the difference is taken with both
                           # phases. The pole amplitude is averaged over the box like the strain (degree-2 factor), and that
                           # mean is ADDED to the strain's dP (amplitudes added: an upper bound, the relative phase of the
                           # loading and strain pressures in the matrix not being known). Self-gravity of the tide is left out.
                           tidal_load_pressure = function() {
                             cfg <- self$config
                             if (!isTRUE(cfg$ocean_load)) return(list(pole = 0, mean = 0, u_ocean = 0, u_core = 0, h2_core = 0))
                             w <- 2 * pi / cfg$orbital_period; e <- cfg$eccentricity; Rc <- cfg$core_radius; RE <- cfg$radius
                             U2 <- function(r) 1.5 * r^2 * w^2 * e
                             gE <- cfg$calculate_gravity(RE); gc <- cfg$calculate_gravity(Rc)
                             tau <- cfg$core_viscosity / cfg$core_rigidity
                             mus <- cfg$core_rigidity * (1i * w * tau) / (1 + 1i * w * tau)
                             h2c <- 5 / 3 * 1.5 / (1 + 19 * mus / (2 * cfg$core_density * (4 / 3 * pi * CONSTANTS$G * cfg$core_density * Rc) * Rc))
                             h2s <- cfg$surface_h2 * exp(-1i * cfg$surface_h2_lag_deg * pi / 180)
                             us <- h2s * U2(RE) / gE; uc <- h2c * U2(Rc) / gc
                             pole <- CONSTANTS$FLUID$density * gc * Mod(us - uc)
                             list(pole = pole, mean = pole * self$tidal_strain()$factor, u_ocean = Mod(us), u_core = Mod(uc),
                                  h2_core = Mod(h2c), h2_core_lag_deg = -Arg(h2c) * 180 / pi)
                           },

                           # Exchange per half-orbit, the smallest of three volumes (PI 2026-10-29; README section 6):
                           #   poro:        alpha * eps_v * V_layer, the fluid the tidal volumetric strain expels from the
                           #                matrix (drained limit). eps_v = eps0 * degree2_dilatation_factor(): the
                           #                compressed part of a degree-2 dilatation field whose pole amplitude is eps0,
                           #                Liao et al. (2020)'s polar strain (an assumption, README).
                           #   darcy:       what drains through the rock in half an orbit, the |q| A integral of
                           #                run_simulation() over the orbit, halved.
                           #   pore_volume: phi * V_layer.
                           # x_half = dV / V_pore per half-orbit; per orbit the two halves compound, x_orbit = 1 - (1 - x_half)^2,
                           # and exchanging_fluid_mass = x_orbit * the porewater mass, so x <= 1 by construction.
                           calculate_exchange = function() {
                             cfg <- self$config; R <- cfg$core_radius; L <- self$reactive_thickness()
                             V_layer <- (4/3) * pi * (R^3 - max(0, R - L)^3)   # the reacting box (the skin, by default)
                             V_pore  <- self$params$porosity * V_layer
                             st      <- self$tidal_strain(); eps0 <- st$eps0; eps_v <- st$eps_v
                             dP      <- self$tidal_pressure()
                             ld      <- self$tidal_load_pressure()
                             dP_load <- if (is.null(cfg$max_tidal_stress)) ld$mean else 0
                             S_      <- self$skin_depth()$S
                             # poro: the strain's fluid volume plus the volume the ocean loading stores in the box, S dP_load L A
                             vols    <- c(poro = self$biot_alpha() * eps_v * V_layer + S_ * dP_load * L * 4 * pi * R^2,
                                          darcy = sum(self$results$fluid_volume) / 2,
                                          pore_volume = V_pore)
                             dV      <- min(vols)
                             x_half  <- dV / V_pore
                             x_orbit <- 1 - (1 - x_half)^2
                             sd      <- self$skin_depth()
                             list(eps0 = eps0, eps_v = eps_v, strain_factor = st$factor, dP = dP,
                                  dP_strain = if (is.null(cfg$max_tidal_stress)) dP - dP_load else NA_real_, dP_load = dP_load,
                                  dP_load_pole = ld$pole, u_ocean = ld$u_ocean, u_core = ld$u_core, h2_core = ld$h2_core,
                                  dP_source = if (is.null(cfg$max_tidal_stress)) "strain" else "override",
                                  V_layer = V_layer, V_pore = V_pore,
                                  delta = sd$delta, L_eff = sd$L_eff, storage = sd$S,
                                  # cross-check: the volume the pressure change stores in the skin, S * dP * delta * A
                                  dV_storage = sd$S * dP * sd$L_eff * 4 * pi * R^2,
                                  dV_poro = vols[["poro"]], dV_darcy = vols[["darcy"]], dV_pore_volume = V_pore,
                                  dV = dV, binding = names(vols)[which.min(vols)],
                                  x_half = x_half, x_orbit = x_orbit,
                                  mass_per_orbit = x_orbit * V_pore * CONSTANTS$FLUID$density)
                           },

                           # Biot coefficient. Critical-porosity model (Nur et al. 1998): the frame's bulk modulus falls
                           # linearly to zero at the critical porosity, K_frame = K_grain (1 - phi/phi_c), so
                           # alpha = 1 - K_frame/K_grain = min(1, phi/phi_c); alpha >= phi whenever phi_c <= 1. A number given
                           # as `alpha` overrides it.
                           biot_alpha = function() {
                             if (!is.null(self$params$alpha)) self$params$alpha else min(1, self$params$porosity / self$params$phi_c)
                           },

                           # Thickness (m) of the box that reacts and exchanges: the top min(layer_thickness, delta) under
                           # reactive_depth = "skin" (only the skin's porewater exchanges in a tide), else the whole layer.
                           # It sets the porewater temperature and pressure (profile means over it), the rock mass and pore
                           # volume (PhreeqcIntegrator$calculate_phreeqc_parameters()) and the exchange (calculate_exchange()).
                           reactive_thickness = function() {
                             if (identical(self$params$reactive_depth, "layer")) self$params$layer_thickness else self$skin_depth()$L_eff
                           },

                           # Poroelastic skin depth, Liao et al. (2020) Table A1 and Fig. B4: delta = sqrt(c / omega), with
                           #   c = (k / eta_f) (K_m + 4/3 mu)(K_u - K_m) / (alpha^2 (K_u + 4/3 mu)),
                           #   K_m = (1 - alpha) K_s, K_u = K_m + alpha^2 K_s K_f / (phi K_s + (alpha - phi) K_f).
                           # (K_u - K_m) / alpha^2 is the Biot modulus M = 1/S, S = phi/K_f + (alpha - phi)/K_s the storage, so
                           # c = (k / eta_f) M (K_m + 4/3 mu)/(K_u + 4/3 mu), the diffusivity k/(eta_f S) times a factor <= 1.
                           # mu is the matrix rigidity (core_rigidity), as in their table. L_eff = min(layer_thickness, delta).
                           skin_depth = function() {
                             cfg <- self$config; a <- self$biot_alpha(); phi <- self$params$porosity
                             Ks <- cfg$rock_bulk_modulus; Kf <- cfg$water_bulk_modulus; mu <- cfg$core_rigidity
                             S  <- phi / Kf + (a - phi) / Ks
                             Km <- (1 - a) * Ks; Ku <- Km + a^2 / S
                             c  <- self$params$k / CONSTANTS$FLUID$viscosity / S * (Km + 4/3 * mu) / (Ku + 4/3 * mu)
                             delta <- sqrt(c / (2 * pi / cfg$orbital_period))
                             list(S = S, c = c, delta = delta, L_eff = min(self$params$layer_thickness, delta))
                           },

                           # One header line on the exchange, for the generated inputs.
                           exchange_header_line = function() {
                             x <- self$exchange
                             sprintf(paste0("# Reacting box: top %.4g m of the %.4g m layer (reactive_depth = \"%s\"; skin depth %.4g m)\n",
                                            "# Exchange per half-orbit (min of three, README section 6): %.3e of the pore volume, limited by %s ",
                                            "(poro %.3e, darcy %.3e, pore_volume %.3e m3; skin depth %.3g m; eps0 %.3e, eps_v/eps0 %.4f, ",
                                            "dP %.4g Pa (%s; ocean load %.3g Pa), mu %.3g GPa, eta %.3g Pa s, alpha %.3f, phi %.3g)"),
                                     self$reactive_thickness(), self$params$layer_thickness, self$params$reactive_depth, x$delta,
                                     x$x_half, x$binding, x$dV_poro, x$dV_darcy, x$dV_pore_volume, x$delta, x$eps0, x$strain_factor,
                                     x$dP, x$dP_source, x$dP_load,
                                     self$config$core_rigidity / 1e9, self$config$core_viscosity, self$biot_alpha(), self$params$porosity)
                           },

                           # Generate summary metrics
                           get_summary_metrics = function() {
                             if (is.null(self$results)) stop("Must run simulation first")
                             
                             if (is.null(self$exchange)) self$exchange <- self$calculate_exchange()
                             # The exchanged mass per orbit is the min-of-three exchange (calculate_exchange()), which
                             # is <= the pore volume by construction; the circulation time is the ocean's mass over it.
                             fluid_mass <- self$exchange$mass_per_orbit
                             circ_time  <- FlowUtils$calculate_circulation_time(
                               flow_rate = fluid_mass / CONSTANTS$FLUID$density / self$config$orbital_period,
                               ocean_mass = self$config$get_derived_properties()$ocean_mass
                             )
                             orbital_period_years <- self$config$orbital_period / (365.25 * 24 * 3600)
                             n_orb <- circ_time / orbital_period_years  # orbital cycles per ocean circulation
                             # (The "corrected" times once capped a Darcy volume at the pore volume; the cap is now
                             # one of the three limits, so they equal the plain ones. Kept for the headers.)
                             corr_circ_time <- circ_time
                             corr_n_orb     <- n_orb

                             metrics <- list(
                               max_flow_rate_Ls = max(self$results$flow_rate_Ls),
                               mean_flow_rate_Ls = mean(self$results$flow_rate_Ls),
                               total_fluid_volume_m3 = sum(self$results$fluid_volume),
                               exchanging_fluid_mass = fluid_mass,
                               darcy_fluid_mass_per_orbit = sum(self$results$fluid_volume) * CONSTANTS$FLUID$density,
                               exchange_fraction_per_half_orbit = self$exchange$x_half,
                               exchange_fraction_per_orbit = self$exchange$x_orbit,
                               exchange_limit = self$exchange$binding,
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
                                   # organic wt% the mineral mass is scaled for; NULL = organic_wt_percent
                                   mineral_basis_wt_percent = NULL,
                                   iom_pools = NULL,
                                   profile_file = "profiles/EnceladusProfile_Seawater_10.0ppt_Tb272.4578K_PorousRock.txt",
                                   
                                   # Factor applied to every primary mineral's moles (the rock minus its organics).
                                   mineral_scale_factor = function() {
                                     basis <- if (is.null(self$mineral_basis_wt_percent)) self$organic_wt_percent else self$mineral_basis_wt_percent
                                     (100 - basis) / 100
                                   },
                                   # "" by default; a header line when mineral_basis_wt_percent is set.
                                   mineral_basis_header_line = function() {
                                     if (is.null(self$mineral_basis_wt_percent)) return("")
                                     sprintf("\n# Mineral basis: minerals scaled by (100 - %g)/100 = %.6g, as for a %g wt%% organic run (organic_wt_percent = %g)",
                                             self$mineral_basis_wt_percent, self$mineral_scale_factor(),
                                             self$mineral_basis_wt_percent, self$organic_wt_percent)
                                   },

                                   initialize = function(simulator) {
                                     if (!inherits(simulator, "HydrothermalSimulator")) {
                                       stop("Must provide HydrothermalSimulator instance")
                                     }
                                     self$simulator <- simulator
                                     self$calculate_phreeqc_parameters()
                                   },
                                   
                                   # Add organic_wt_percent parameter to the PhreeqcIntegrator
                                   # mineral_basis_wt_percent: scale the minerals by (100 - this)/100 instead of
                                   # (100 - organic_wt_percent)/100. pair_organic_run() sets it to the organic
                                   # run's value for an organics-absent blank with the same mineral inventory.
                                   calculate_phreeqc_parameters = function(organic_wt_percent = 0, mineral_basis_wt_percent = NULL) {
                                     if (!is.null(mineral_basis_wt_percent) &&
                                         !(is.numeric(mineral_basis_wt_percent) && length(mineral_basis_wt_percent) == 1 &&
                                           isTRUE(mineral_basis_wt_percent >= 0 && mineral_basis_wt_percent < 100)))
                                       stop("mineral_basis_wt_percent must be NULL or one number in [0, 100)")
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
                                     layer_thick_gf <- self$simulator$reactive_thickness()   # the reacting box (skin by default)
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
                                     # Its solid mass at the core's bulk density and the Simulator's porosity (phi0):
                                     # V (rho_core - phi rho_water). Equals layer_fraction * props$core_mass when phi0 is the
                                     # config's core_porosity (0.32).
                                     layer_mass <- layer_vol * (self$simulator$config$core_density -
                                                                self$simulator$params$porosity * CONSTANTS$FLUID$density)

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

                                     # The ocean's mass for the PHREEQC basis: the profile's (Mocean_kg, PI 2026-10-31), so that the
                                     # ocean per kg of the box's rock is the geometric ratio; the config's geometric ocean (at 1000
                                     # kg/m3) is 0.9% lighter and is kept for the Simulator's circulation time.
                                     ocean_mass_phreeqc <- tryCatch(profile_ocean_mass_kg(self$profile_file), error = function(e) NA_real_)
                                     if (!is.finite(ocean_mass_phreeqc)) ocean_mass_phreeqc <- props$ocean_mass
                                     water_mass <- ocean_mass_phreeqc / layer_mass

                                     # Store parameters
                                     self$organic_wt_percent <- organic_wt_percent
                                     self$mineral_basis_wt_percent <- mineral_basis_wt_percent
                                     # Porewater mass: physical mass of water filling the permeable layer pore space
                                     r_core <- self$simulator$config$core_radius
                                     layer_thick <- self$simulator$reactive_thickness()
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
                                       ocean_mass = ocean_mass_phreeqc,
                                       core_mass = props$core_mass
                                     )

                                     print("Calculated PHREEQC parameters:")
                                     print(sprintf("rock_per_cycle: %.2e", rock_per_cycle))
                                     print(sprintf("rock_per_year: %.2e", rock_per_year))
                                     print(sprintf("water_mass: %.2f", water_mass))
                                     
                                     return(self$phreeqc_params)
                                   },

                                   # Surface areas from the grain size (item 3): with the Simulator's ssa = "geometric", each RATES
                                   # entry's "SSA = <value>" line gets 6 / (rho_i d_i) x roughness (m2/g); "legacy" leaves the input
                                   # unchanged. rho_i = molar mass / Vm (data/mineral_molar_volumes.csv); a mineral without a
                                   # density there is an error.
                                   apply_ssa = function(content) {
                                     pr <- self$simulator$params
                                     if (!identical(pr$ssa, "geometric")) return(content)
                                     v <- read_mineral_volumes()
                                     l <- strsplit(content, "\n", fixed = TRUE)[[1]]
                                     in_rates <- FALSE; cur <- NA_character_
                                     for (i in seq_along(l)) {
                                       if (grepl("^RATES\\b", l[i])) { in_rates <- TRUE; cur <- NA_character_; next }
                                       if (in_rates && grepl("^(END|KINETICS|EQUILIBRIUM_PHASES|SOLUTION|SELECTED_OUTPUT|USER_PUNCH|USER_PRINT|USE|MIX|SAVE|KNOBS|GAS_PHASE|REACTION|PHASES|TITLE|INCREMENTAL)\\b", l[i]))
                                         in_rates <- FALSE
                                       if (!in_rates) next
                                       if (grepl("^[A-Za-z]", l[i])) { cur <- sub("[ \t#].*$", "", l[i]); next }
                                       if (!is.na(cur) && grepl("SSA = [0-9.eE+-]+", l[i])) {
                                         j <- match(cur, v$phase)
                                         if (is.na(j)) stop("apply_ssa: no density for ", cur)
                                         rho <- v$molar_mass_g_mol[j] / v$Vm_cm3_mol[j]
                                         d <- if (!is.null(pr$grain_diameter_by_mineral) && cur %in% names(pr$grain_diameter_by_mineral))
                                           pr$grain_diameter_by_mineral[[cur]] else pr$grain_diameter
                                         s <- pr$roughness * geometric_ssa_m2_g(d, rho)
                                         l[i] <- sub("SSA = [0-9.eE+-]+.*$", sprintf("SSA = %.6g  # geometric 6/(rho d) x roughness: rho %.4g g/cm3, d %.4g m, roughness %g (item 3)",
                                                                                   s, rho, d, pr$roughness), l[i])
                                       }
                                     }
                                     paste(l, collapse = "\n")
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

                                   # Insert provenance lines after the TITLE line of a rendered .pqi:
                                   # the PHREEQC version these inputs are written for, and the database
                                   # they need. run_phreeqc() (phreeqc_tools.R) reads the version line back
                                   # and warns if the executable reports a different one.
                                   # With pressures = list(pore =, ocean =) (atm), also the two pressure lines
                                   # that a restart checks against (check_restart_pressures()).
                                   stamp_pqi_header = function(content, database, pressures = NULL) {
                                     lines <- strsplit(paste(content, collapse = "\n"), "\n", fixed = TRUE)[[1]]
                                     dropped <- attr(content, "dropped")
                                     stamp <- c(sprintf("# Target PHREEQC version: %s", PHREEQC_EXPECTED_VERSION),
                                                sprintf("# Database: %s", database),
                                                if (!is.null(pressures)) pressure_header_lines(pressures$pore, pressures$ocean),
                                                if (length(dropped)) sprintf("# Not defined in %s, so removed from this input: %s",
                                                                             basename(database), paste(dropped, collapse = " ")))
                                     at <- if (length(lines) && startsWith(lines[1], "TITLE")) 1 else 0
                                     paste(append(lines, stamp, after = at), collapse = "\n")
                                   },

                                   # Carbon mode for a generated .pqi file, and the database it needs:
                                   #   "coupled"       -- coreclath_CH4coupled.dat; all carbon redox-coupled (CH4/CH4(g)).
                                   #   "ch4_uncoupled" -- coreclath_CH4uncoupled.dat; methane decoupled as Mtg/Mtg(g),
                                   #                      other reduced carbon (C2H6, acetate, formate, ...) coupled.
                                   #   "c4_only"       -- coreclath_C4only.dat; methane is Mtg, and coupled carbon is
                                   #                      C(+4) only (the database comments out the rest).
                                   #   "temperature"   -- the rule that was the default before 2026-10-26:
                                   #                      "ch4_uncoupled" below ch4_cutoff_C (150 C), "coupled" at/above.
                                   # Default (carbon_mode = NULL): "c4_only" at every temperature. Coupled runs
                                   # let equilibrium reduce carbonate with H2, so they are an upper bound on
                                   # abiotic methanation, not the expected case. CH4_redox_override =
                                   # "coupled"/"decoupled" is the older spelling of carbon_mode =
                                   # "coupled"/"ch4_uncoupled"; give one or the other. ch4_cutoff_C is used only
                                   # by "temperature". `decoupled` (methane is Mtg) is TRUE for both
                                   # "ch4_uncoupled" and "c4_only"; `explicit` is TRUE when the caller chose the
                                   # mode (carbon_mode or CH4_redox_override given).
                                   resolve_ch4_mode = function(porewater_temp_C, CH4_redox_override = NULL, carbon_mode = NULL,
                                                               ch4_cutoff_C = 150) {
                                     if (!is.null(CH4_redox_override) && !CH4_redox_override %in% c("coupled", "decoupled")) {
                                       stop('CH4_redox_override must be NULL, "coupled", or "decoupled"')
                                     }
                                     if (!is.null(carbon_mode) && !carbon_mode %in% c("coupled", "ch4_uncoupled", "c4_only", "temperature")) {
                                       stop('carbon_mode must be NULL, "coupled", "ch4_uncoupled", "c4_only", or "temperature"')
                                     }
                                     if (!is.null(carbon_mode) && !is.null(CH4_redox_override)) {
                                       stop("give carbon_mode or CH4_redox_override, not both")
                                     }
                                     explicit <- !is.null(carbon_mode) || !is.null(CH4_redox_override)
                                     # The header note: "temperature" writes what the old default wrote, so its
                                     # inputs match those generated before the default changed.
                                     note <- if (!is.null(CH4_redox_override)) sprintf(", override=%s", CH4_redox_override)
                                             else if (identical(carbon_mode, "temperature")) {
                                               if (ch4_cutoff_C != 150) sprintf(", cutoff=%g C", ch4_cutoff_C) else ""
                                             }
                                             else sprintf(", carbon_mode=%s", if (is.null(carbon_mode)) "c4_only" else carbon_mode)
                                     carbon_mode <- if (!is.null(CH4_redox_override)) {
                                       if (CH4_redox_override == "decoupled") "ch4_uncoupled" else "coupled"
                                     } else if (is.null(carbon_mode)) "c4_only"
                                     else if (carbon_mode == "temperature") {
                                       if (porewater_temp_C < ch4_cutoff_C) "ch4_uncoupled" else "coupled"
                                     } else carbon_mode
                                     decoupled <- carbon_mode != "coupled"
                                     list(
                                       carbon_mode   = carbon_mode,
                                       explicit      = explicit,
                                       decoupled     = decoupled,
                                       ch4_aq        = if (decoupled) "Mtg"    else "CH4",
                                       ch4_gas       = if (decoupled) "Mtg(g)" else "CH4(g)",
                                       database      = switch(carbon_mode,
                                                              coupled       = "databases/coreclath_CH4coupled.dat",
                                                              ch4_uncoupled = "databases/coreclath_CH4uncoupled.dat",
                                                              c4_only       = "databases/coreclath_C4only.dat"),
                                       hydrate_phase = if (decoupled) "Mtg_hydrate" else "CH4_hydrate",
                                       description   = switch(carbon_mode,
                                                              coupled       = "coupled -> aqueous/gas species are CH4/CH4(g)",
                                                              ch4_uncoupled = "decoupled -> aqueous/gas species are Mtg/Mtg(g)",
                                                              c4_only       = "decoupled -> aqueous/gas species are Mtg/Mtg(g); coupled C is C(+4) only"),
                                       note          = note
                                     )
                                   },

                                   # IOM channel config matching the carbon mode's database. With the uncoupled
                                   # database (no coupled CH4), IOM_CH4 releases its carbon as Mtg
                                   # (redox_mode C = "decoupled": "Mtg 1 H 0.435"); with C4only, IOM_CHn too
                                   # (C = "c4_only", see iom_default_config()); otherwise everything is
                                   # coupled. N and S stay coupled. Rebuilt on every call so a config from
                                   # one temperature is never reused at another; a config the caller set
                                   # in self$iom_pools (without the "auto_ch4_mode" attribute) is kept.
                                   resolve_iom_pools = function(ch4_mode) {
                                     if (!is.null(self$iom_pools) && !isTRUE(attr(self$iom_pools, "auto_ch4_mode")))
                                       return(invisible(self$iom_pools))
                                     cfg <- iom_default_config(redox_mode = c(C = switch(ch4_mode$carbon_mode, coupled = "coupled",
                                                                                         ch4_uncoupled = "decoupled", c4_only = "c4_only"),
                                                                              N = "coupled", S = "coupled"))
                                     attr(cfg, "auto_ch4_mode") <- TRUE
                                     self$iom_pools <- cfg
                                     invisible(cfg)
                                   },

                                   # Rewrites literal CH4/CH4(g) tokens to Mtg/Mtg(g) in already-rendered
                                   # PHREEQC input text, for decoupled-mode files. Applied once to the final
                                   # assembled text so every block (static text and vector-built lists alike)
                                   # is covered without threading a placeholder through each one individually.
                                   # With carbon_mode "c4_only" it then removes every phase and species the
                                   # C4only database does not define but its source, coreclath_CH4uncoupled.dat,
                                   # does (the reduced-carbon entries the templates list), via
                                   # restrict_to_database() (phreeqc_tools.R). stamp_pqi_header() records what
                                   # was removed.
                                   apply_ch4_mode = function(text, ch4_mode) {
                                     if (!ch4_mode$decoupled) return(text)
                                     text <- gsub("CH4\\(g\\)", "Mtg(g)", text)
                                     text <- gsub("(?<![A-Za-z0-9_])CH4(?![A-Za-z0-9_\\(])", "Mtg", text, perl = TRUE)
                                     if (identical(ch4_mode$carbon_mode, "c4_only")) {
                                       text <- restrict_to_database(text, ch4_mode$database, base = "databases/coreclath_CH4uncoupled.dat")
                                       if (length(attr(text, "dropped_nonzero")))
                                         warning("carbon_mode c4_only: removed phases that had a nonzero amount: ",
                                                 paste(attr(text, "dropped_nonzero"), collapse = ", "))
                                     }
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

                                   # Writes the main ENTICES input: porewater reacting with the kinetic rock (silicates, troilite, schreibersite,
                                   # IOM channels) and equilibrium phases, exchanging fluid with the ocean every step (REACT, MIX 4 ocean +
                                   # seafloor phases, MIX 1, MIX 2). Porewater and ocean each have a fixed-pressure gas phase; pressures come
                                   # from profile_file. Inputs only; run them with run_phreeqc().
                                   #
                                   # @param output_dir Directory for the .pqi (and, when run there, its .csv/.out). Default "phreeqc_inputs".
                                   # @param mode "ocean_circulations" (default: auto grouping, n_circulations), "ungrouped" (1 orbit per step, n_orbital_cycles) or "grouped" (grouping_factor orbits per step; exchange compounds as 1 - (1 - x_tide)^G).
                                   # @param n_circulations Number of ocean circulations (mode "ocean_circulations").
                                   # @param n_orbital_cycles Number of orbits (1.375 d each).
                                   # @param grouping_factor Orbits per step (mode "grouped").
                                   # @param exchange_ratio Target exchanged fluid / porewater per step; sets grouping_factor and mode "grouped". Default NULL.
                                   # @param porewater_temp Porewater reaction temperature: a number (deg C), "profile" (the profile's layer mean) or "ocean" (1.45 C, the profile's ocean floor). Default "ocean". Production: 150 and 160 C.
                                   # @param substeps KINETICS substeps per step. Default 1.
                                   # @param suffix Appended to the file name. Default "".
                                   # @param CH4_redox_override Older spelling of carbon_mode: "coupled" or "decoupled" (= "ch4_uncoupled"). Default NULL; give this or carbon_mode, not both.
                                   # @param carbon_mode "c4_only" (NULL, the default: coupled carbon is C(+4) only, decided by the PI 2026-10-27), "coupled" (an upper bound on abiotic methanation), "ch4_uncoupled" (diagnostic), or "temperature" (the pre-10-26 rule). Picks the database.
                                   # @param ch4_cutoff_C Cutoff (deg C) for carbon_mode = "temperature" only. Default 150.
                                   # @param trace_C_molal Numerical C seed in every initial SOLUTION, mol/kgw. Default NULL; production 1e-14 (the smallest tested to converge at 25 C).
                                   # @param initial_ocean "pure_water" (default, for attribution), a named list of element totals (mol/kgw), or cometary_ice_ocean() (Miller et al. 2026 Table 5).
                                   # @param charge_balance NULL (default), "pH", or an initial_ocean element to charge-balance the ocean SOLUTION 2 on. Not used with cometary_ice_ocean().
                                   # @param pore_gas_fate "to_ocean" (default, PI 2026-10-29: the pore bubble moves to the ocean step each orbit; without plume loss an upper bound on ocean H2) or "retain" (it stays with the porewater; an upper bound on retention, ~2.2 pore volumes of gas). Carbon attribution is the same in both.
                                   # @return The path of the .pqi written.
                                   # @examples
                                   #   sim <- Simulator$new(k = 1e-7, layer_thickness = 194200, config = EnceladusConfig$new()); sim$run_simulation()
                                   #   ph <- PhreeqcIntegrator$new(sim); ph$calculate_phreeqc_parameters(organic_wt_percent = 15)
                                   #   f <- ph$generate_phreeqc_input_tidal_kinetic(mode = "ungrouped", n_orbital_cycles = 2, porewater_temp = 25, trace_C_molal = 1e-14)
                                   generate_phreeqc_input_tidal_kinetic = function(output_dir = "phreeqc_inputs", mode = "ocean_circulations",
                                                                     n_circulations = NULL, n_orbital_cycles = NULL,
                                                                     grouping_factor = NULL, exchange_ratio = NULL,
                                                                     porewater_temp = "ocean", substeps = 1, suffix = "",
                                                                     CH4_redox_override = NULL, carbon_mode = NULL, ch4_cutoff_C = 150, trace_C_molal = NULL,
                                                                     initial_ocean = "pure_water", charge_balance = NULL,
                                                                     pore_gas_fate = c("to_ocean", "retain"), logfile = FALSE) {
                                     pore_gas_fate <- match.arg(pore_gas_fate)

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
                                     content <- self$create_phreeqc_template_tidal_kinetic(filename, n, mode, grouping_factor = gf, exchange_ratio = exchange_ratio, porewater_temp = porewater_temp, substeps = substeps, CH4_redox_override = CH4_redox_override, carbon_mode = carbon_mode, ch4_cutoff_C = ch4_cutoff_C,
                                                                                trace_C_molal = trace_C_molal, initial_ocean = initial_ocean, charge_balance = charge_balance,
                                                                                pore_gas_fate = pore_gas_fate)

                                     # Write to file
                                     writeLines(knobs_logfile(self$apply_ssa(content), logfile), filename)

                                     return(filename)
                                   },
                                   
                                   create_phreeqc_template_tidal_kinetic = function(filename, n, mode = "ocean_circulations", grouping_factor = NULL, exchange_ratio = NULL, porewater_temp = "ocean", substeps = 1, CH4_redox_override = NULL, carbon_mode = NULL, ch4_cutoff_C = 150, trace_C_molal = NULL, initial_ocean = "pure_water", charge_balance = NULL, pore_gas_fate = "to_ocean") {
                                     trace_C_line <- trace_C_solution_line(trace_C_molal)
                                     ocean_lines  <- paste0(trace_C_line, initial_ocean_solution_lines(initial_ocean, trace_C_molal, charge_balance))
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

                                     # Temperatures and pressures for SOLUTION blocks. The ocean temperature is fixed.
                                     ocean_temp_C <- 274.6 - 273.15  # the profile's ocean floor, top of the rock layer (274.566 K), rounded
                                     # Porewater and ocean conditions from the profile; stops if profile_file is missing.
                                     pc <- profile_conditions(self$profile_file, self$simulator$reactive_thickness())
                                     pore_temp_C        <- pc$pore_temp_C
                                     pore_pressure_atm  <- pc$pore_pressure_atm
                                     ocean_pressure_atm <- pc$ocean_pressure_atm

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
                                     ch4_mode <- self$resolve_ch4_mode(porewater_temp_C, CH4_redox_override, carbon_mode, ch4_cutoff_C)
                                     ch4_override_note <- ch4_mode$note
                                     # Gas phases (see gas_components()): the porewater one, the ocean one, the
                                     # punch columns, and the empty porewater phase to_ocean starts each orbit with.
                                     cometary <- inherits(initial_ocean, "cometary_ice")
                                     # Tested 2026-10-28: a cometary ocean with a carbon-free porewater (no seed, no IOM)
                                     # left the first porewater step running for > 20 min without converging; with
                                     # trace_C_molal = 1e-14 it completes.
                                     if (cometary && is.null(trace_C_molal) && self$organic_wt_percent == 0)
                                       warning("cometary_ice_ocean() with a carbon-free porewater: give trace_C_molal = 1e-14 (without it the first porewater step did not converge in tests)")
                                     gas_comps <- gas_components(ch4_mode, if (cometary) "Ntg(g)" else character(0))
                                     gp <- gas_punch(gas_comps, porewater_mass_norm, put = pore_gas_fate == "to_ocean")
                                     pore_gas_block <- gas_phase_block(1, gas_comps, pore_pressure_atm, porewater_temp_C,
                                                                       "Porewater gas (fixed pressure, starts empty)")
                                     ocean_gas_block <- gas_phase_block(2, gas_comps, ocean_pressure_atm, ocean_temp_C,
                                                                        "Ocean gas (fixed pressure, starts empty, never removed)")
                                     # A cometary-ice ocean: the ices react into the pure-water ocean in a simulation of their
                                     # own, before any reactant is defined (cometary_ocean_simulation()), which also defines
                                     # and saves GAS_PHASE 2; otherwise GAS_PHASE 2 is defined empty with the others.
                                     cometary_block <- if (cometary) paste0(cometary_ocean_simulation(initial_ocean, ch4_mode, ocean_mass_norm,
                                                                            ocean_gas_block, base_name, ch4_mode$database), "\n") else ""
                                     ocean_main_block <- if (cometary) "# GAS_PHASE 2 (ocean): defined, and saved with the cometary ices, above" else ocean_gas_block
                                     ch4_db_note_line <- sprintf(
                                       "# Recommended database: %s (CH4 %s; porewater T = %.1f C%s)",
                                       ch4_mode$database,
                                       ch4_mode$description,
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
                                     mineral_scale_factor <- self$mineral_scale_factor()
                                     
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
                                        self$resolve_iom_pools(ch4_mode)

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
{trace_C_header_line(trace_C_molal)}{self$mineral_basis_header_line()}
{initial_ocean_header_line(initial_ocean, charge_balance)}
{ch4_db_note_line}
# pore_gas_fate: {pore_gas_fate}
# Gas phases (fixed pressure, Peng-Robinson): {paste(gas_comps, collapse = " ")}
{clathrate_note_lines}

KNOBS
   -logfile false

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
# Exchanged fluid mass (per orbital cycle): {sprintf("%.2e", metrics$exchanging_fluid_mass)} kg
{self$simulator$exchange_header_line()}
# Porewater mass: {sprintf("%.2e", self$phreeqc_params$porewater_mass)} kg
# Ocean mass: {sprintf("%.2e", self$phreeqc_params$ocean_mass)} kg

# Initial fluids are defined at 25 C, where pe 4 adds no oxidant or reductant (O(0) and H(0)
# negligible), and are brought to the run temperature by REACTION_TEMPERATURE. Defined at the
# run temperature, pe 4 adds O2: 4.9 molal O(0) at 250 C, and does not converge at 300 C.
SOLUTION 1 Enceladus porewater
    temp      25
    pH        7
    pe        4
    redox     pe
    units     mol/kgw
    density   1
    water     {porewater_mass_norm}
    pressure  {round(pore_pressure_atm, 1)}{trace_C_line}

SOLUTION 2 Enceladus ocean
    temp      25
    pH        7{ocean_pH_charge(charge_balance)}
    pe        4
    redox     pe
    units     mol/kgw
    density   1
    water     {ocean_mass_norm}
    pressure  {round(ocean_pressure_atm, 1)}{ocean_lines}

REACTION_TEMPERATURE 1
    {round(porewater_temp_C, 2)}

REACTION_TEMPERATURE 2
    {round(ocean_temp_C, 2)}

# Reaction pressures, used by every step like the reaction temperatures. MIX averages the
# solution pressures by mass, so without these the porewater would drift towards the ocean
# pressure after each exchange.
REACTION_PRESSURE 1
    {round(pore_pressure_atm, 1)}

REACTION_PRESSURE 2
    {round(ocean_pressure_atm, 1)}

{cometary_block}{pore_gas_block}

{ocean_main_block}
{if (pore_gas_fate == "to_ocean") paste0("\n", gas_transfer_block(gas_comps)) else ""}

# Define rate equations for each mineral
RATES
Forsterite
	-start
	1   REM Ref PK04
	10  kacid = 10^(-6.85) * exp(-67.2e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^0.47
	20  kneut = 10^(-10.64) * exp(-79.0e3/8.314 * (1/TK-1/298.15))
	21  SSA = 0.1  # UNSOURCED surface area (m2/g): README section 6, decision D7
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
  26  SSA = 5  # UNSOURCED surface area (m2/g): README section 6, decision D7
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
  21  SSA = 0.1  # UNSOURCED surface area (m2/g): README section 6, decision D7
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
  31  SSA = 5  # UNSOURCED surface area (m2/g): README section 6, decision D7
  32  mw = 278.33
  40  k = (kacid + kneut + kbase) * SSA * mw * M
  50  IF SR("K-Feldspar") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("K-Feldspar"))
  60  moles = rate * TIME
  70  SAVE moles
  -end

Tephroite
  -start
  1   REM Pincus et al. (2026) ACS Earth Space Chem 10, 1174-1184, Table S1,
  2   REM citing Casey et al. (1993) GCA 57, 785-793. log kH+ = -5.4 mol
  3   REM m-2 s-1 (measured, 25C, pH-2 reference), n = 0.47. Cross-checked
  4   REM 2026-09-15 by direct regression of Casey\'s own raw rate-vs-pH data
  5   REM (25C, synthetic tephroite): independently derived value agrees to
  6   REM 0.14 log units. See RECONSTRUCTION_STATUS.md S4/S6 for full derivation.
  7   REM Pincus\'s tabulated -5.4 is the RATE EVALUATED AT pH 2 (their own SI
  8   REM convention), not the pre-exponential this PHASES-block form needs.
  9   REM Pre-exponential = log_kH+ + n*2 (same relationship that reproduces
  10  REM Core11\'s own Fayalite entry: -4.80 + 1.0*2 = -6.8, Pincus\'s fayalite
  11  REM value -- checked here, not assumed). For tephroite: -5.4 + 0.47*2 =
  12  REM -4.46. Hand-checked: 10^(-4.46)*ACT(H+)^0.47 at pH2 = 10^(-5.4), exact.
  13  REM
  14  REM Ea: Casey et al. (1993) Fig. 4 gives Eexp = 54, 26, 22 and 4 kJ/mol at
  15  REM pH 2.5, 3.5, 4.2 and 5.15, and < 4 at pH 6.2 (a slope not distinct from 0).
  16  REM Taking 0 at pH 6.2, a least-squares fit Ea = a ln(pH) + b to these 5 points
  17  REM gives a = -60.11, b = 106.17 (not the appendix rate table: see PROGRESS).
  18  REM The fit is not forced through (6.2, 0): it reaches 0 at pH 5.85 and is
  19  REM -3.5 at pH 6.2. Ea >= 0 is imposed, so Ea = 0 for pH > 5.85.
  22  REM
  23  REM ANOXIC CONDITIONS ONLY. Casey\'s experiments were run rigorously
  24  REM O2-free (N2-degassed, Fe-filing/hot-Cu scrubbed, <6e-6 mol/L O2);
  25  REM Mn2+ -> Mn3+/4+ oxidation "severely retards" dissolution per the
  26  REM source. Do not reuse this rate law under any oxidizing scenario
  27  REM without revisiting.
  28  REM
  29  REM No neutral/base branch: source data at pH>7 are internally
  30  REM inconsistent between methods by up to 1.8x (STATUS.md S4/S6) and
  31  REM not fittable. Neutral rate at pH > 7 assumed to be same as pH=7. 
  32  REM This is generally consistent with other minerals that lack base
  33  REM dissolution mechanisms. 

  34  REM pH >= 7: acid branch evaluated at a(H+) = 1e-7 (flat at the pH-7 rate). No neutral term is summed.
  35  aH = ACT("H+")
  36  IF aH < 1e-7 THEN aH = 1e-7
  37  Ea = -60.11*LOG(-LOG10(aH)) + 106.17
  38  IF Ea < 0 THEN Ea = 0
  40  kacid = 10^(-4.46) * exp(-Ea*1e3/8.314 * (1/TK-1/298.15)) * aH^0.47
  46  SSA = 0.1  # UNSOURCED surface area (m2/g): README section 6, decision D7
  47  mw = 201.96
  50  k = kacid * SSA * mw * M
  60  IF SR("Tephroite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Tephroite"))
  70  moles = rate * TIME
  80  SAVE moles
  -end

{SCHREIBERSITE_RATE}

Enstatite
  -start
  1   REM Ref PK04
  10  kacid = 10^(-9.02) * exp(-80.0e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^0.6
  20  kneut = 10^(-12.72) * exp(-80.0e3/8.314 * (1/TK-1/298.15))
  21  SSA = 0.1  # UNSOURCED surface area (m2/g): README section 6, decision D7
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
  21  SSA = 5  # UNSOURCED surface area (m2/g): README section 6, decision D7
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
  31  SSA = 5  # UNSOURCED surface area (m2/g): README section 6, decision D7
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
# --- Primary minerals (equilibrium; dissolve_only so they can only deplete). Fe and Ni corrode completely
#     in the first step: UNSOURCED instantaneous corrosion, README section 6, decision D6 ---
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
   Calcite 0 0
#   Celadonite 0 0
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
#   H2O(g) 0 0
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
   Mesolite 0 0
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
   Natrolite 0 0
#   Na2CO3 0 0
#   Na2CO3:7H2O 0 0
#   Nahcolite 0 0
#   Natron 0 0
#   Nesquehonite 0 0
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
   Calcite 0 0
#   Celadonite 0 0
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
#   H2O(g) 0 0
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
   Mesolite 0 0
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
   Natrolite 0 0
#   Na2CO3 0 0
#   Na2CO3:7H2O 0 0
#   Nahcolite 0 0
#   Natron 0 0
#   Nesquehonite 0 0
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
    -high_precision true
    -reaction true
    -state True
    -solution True
    -step                 true
    -ph                   true
    -pe                   true
    -alkalinity           true
    -ionic_strength       true
    -water                true
{totals_block(ch4_mode$database)}
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
    -equilibrium_phases  Ni Lawrencite Analcime Anhydrite Aragonite Bassanite
                         Beidellite-Ca Beidellite-Fe Beidellite-Mg Beidellite-Na
                         Boehmite Brucite
                         Calcite Chalcedony Chamosite Chrysotile
                         Citric_Acid Clinochlore-14A Clinochlore-7A Clinoptilolite-Ca Clinoptilolite-Na
                         Clinozoisite Cronstedtite-7A
                         Daphnite-14A Daphnite-7A Dawsonite Diaspore Dolomite
                         Ettringite Fe(OH)2 Fe(OH)3 Fe2(SO4)3 FeSO4
                         Gibbsite Glycine Goethite Greenalite Gypsum Gyrolite
                         Hematite Huntite Hydromagnesite Kaolinite
                         Lizardite Magnetite Melanterite Mesolite MgSO4 Minnesotaite Mirabilite
                         Montmor-Ca Montmor-Mg Montmor-Na
                         Natrolite Nontronite-Ca Nontronite-Mg Nontronite-Na
                         Portlandite Pyrite
                         Saponite-Fe-Ca Saponite-Fe-Fe Saponite-Fe-Mg Saponite-Fe-Na
                         Saponite-Mg-Ca Saponite-Mg-Fe Saponite-Mg-Mg Saponite-Mg-Na
                         Sepiolite SiO2(am) Thenardite Tobermorite-11A
                         {extra_eq_so_line}
{gas_so_line(gas_comps)}

# -high_precision (above) also sets the mole-balance convergence criterion to 1e-12; with a gas phase
# that made the first porewater step (the metal pulse) run for more than 15 minutes. This KNOBS, read
# after it, puts back the default 1e-8 (the last setting read wins, PHREEQC read.cpp), and keeps the
# 12-digit output.
KNOBS
    -convergence_tolerance 1e-8

    
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
                                         # to_ocean: the porewater gas phase is redefined empty at the start of each orbit
                                         if (pore_gas_fate == "to_ocean" && i > 1) paste0(pore_gas_block, "\n") else "",
                                         "USE solution 1\n",
                                         "USE kinetics 1\n",
                                         "USE equilibrium_phases 1\n",
                                         "USE gas_phase 1\n",
                                         "USE reaction none\n",
                                         "USE reaction_temperature 1\n",
                                         "USE reaction_pressure 1\n",
                                         "SAVE solution 1\n",
                                         "SAVE equilibrium_phases 1\n",
                                         if (pore_gas_fate == "retain") "SAVE gas_phase 1\n" else "",
                                         "USER_PUNCH\n",
                                         "   -headings Time_Years Enst_remain_g Forst_remain_g Troil_remain_g Alb_remain_g Diop_remain_g Anorth_remain_g Kfs_remain_g Teph_remain_g Schr_remain_g Total_Initial_g",
                                         iom_punch_headings, SHARED_PUNCH_HEADINGS, gp$headings, "\n",
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
                                         SHARED_PUNCH_LINES, gp$lines,
                                         "END\n\n",
                                         # Simulation 2: ocean absorbs flowing fluid directly + seafloor equilibration
                                         # SOLUTION 4: f_fluid/f_pore of porewater + ocean, fractional mass = f_ocean + f_fluid
                                         OCEAN_STEP_KNOBS_ON,
                                         if (pore_gas_fate == "to_ocean") "USE kinetics 3\n" else "USE kinetics none\n",  # 3: the pore bubble
                                         "USE equilibrium_phases 2\n",
                                         "USE gas_phase 2\n",
                                         "USE reaction none\n",
                                         "USE reaction_temperature 2\n",
                                         "USE reaction_pressure 2\n",
                                         "MIX 4\n",
                                          sprintf("  1  %.15e\n", x_step),
                                          "  2  1.0\n",
                                         "SAVE solution 4\n",
                                         "SAVE equilibrium_phases 2\n",
                                         "SAVE gas_phase 2\n",
                                         "END\n\n",
                                         OCEAN_STEP_KNOBS_OFF,
                                         # Simulation 3: resample porewater — replace flowing fluid with combined ocean water
                                         # (f_pore-f_x)/f_pore of reacted porewater + f_x/(f_ocean+f_x) of combined ocean
                                         # fractional mass = f_pore
                                         "USE kinetics none\n",
                                         "USE equilibrium_phases none\n",
                                         "USE reaction none\n",
                                         "USE reaction_temperature 1\n",
                                         "USE reaction_pressure 1\n",
                                         "MIX 1\n",
                                          sprintf("  1  %.15e\n", 1 - x_step),
                                          sprintf("  4  %.15e\n", f_x / (f_ocean + f_x)),
                                         "SAVE solution 1\n",
                                         "END\n\n",
                                         # Simulation 4: update ocean from combined ocean, fractional mass = f_ocean
                                         "USE kinetics none\n",
                                         "USE equilibrium_phases none\n",
                                         "USE reaction none\n",
                                         "MIX 2\n",
                                          sprintf("  4  %.15e\n", f_ocean / (f_ocean + f_x)),
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

                                     return(self$stamp_pqi_header(self$apply_ch4_mode(paste0(template, "\n", cycle_blocks), ch4_mode),
                                                                  ch4_mode$database, list(pore = pore_pressure_atm, ocean = ocean_pressure_atm)))
                                   },

                                   # ---------------------------------------------------------------
                                   # Equilibrium-porewater variant: primary minerals dissolve via
                                   # EQUILIBRIUM_PHASES instead of kinetics.  Everything else
                                   # (mixing, ocean, seafloor EQ2) is identical.
                                   # An inorganic limiting case: no IOM and no Schreibersite (see
                                   # EQUIL_INORGANIC_NOTE, written into the .pqi header).
                                   # ---------------------------------------------------------------
                                   # As tidal_kinetic, with every primary mineral as a dissolve_only equilibrium phase: an inorganic limiting
                                   # case (no IOM, no schreibersite).
                                   #
                                   # @param output_dir Directory for the .pqi (and, when run there, its .csv/.out). Default "phreeqc_inputs".
                                   # @param n_orbital_cycles Number of orbits (1.375 d each).
                                   # @param porewater_temp Porewater reaction temperature: a number (deg C), "profile" (the profile's layer mean) or "ocean" (1.45 C, the profile's ocean floor). Default "ocean". Production: 150 and 160 C.
                                   # @param suffix Appended to the file name. Default "".
                                   # @param CH4_redox_override Older spelling of carbon_mode: "coupled" or "decoupled" (= "ch4_uncoupled"). Default NULL; give this or carbon_mode, not both.
                                   # @param carbon_mode "c4_only" (NULL, the default: coupled carbon is C(+4) only, decided by the PI 2026-10-27), "coupled" (an upper bound on abiotic methanation), "ch4_uncoupled" (diagnostic), or "temperature" (the pre-10-26 rule). Picks the database.
                                   # @param ch4_cutoff_C Cutoff (deg C) for carbon_mode = "temperature" only. Default 150.
                                   # @param initial_ocean "pure_water" (default, for attribution), a named list of element totals (mol/kgw), or cometary_ice_ocean() (Miller et al. 2026 Table 5).
                                   # @param charge_balance NULL (default), "pH", or an initial_ocean element to charge-balance the ocean SOLUTION 2 on. Not used with cometary_ice_ocean().
                                   # @param pore_gas_fate "to_ocean" (default, PI 2026-10-29: the pore bubble moves to the ocean step each orbit; without plume loss an upper bound on ocean H2) or "retain" (it stays with the porewater; an upper bound on retention, ~2.2 pore volumes of gas). Carbon attribution is the same in both.
                                   # @param trace_C_molal Numerical C seed in both initial SOLUTIONs, mol/kgw. Default NULL; give 1e-14 with a cometary ocean.
                                   # @return The path of the .pqi written.
                                   # @examples
                                   #   sim <- Simulator$new(k = 1e-7, layer_thickness = 194200, config = EnceladusConfig$new()); sim$run_simulation()
                                   #   ph <- PhreeqcIntegrator$new(sim); ph$calculate_phreeqc_parameters(organic_wt_percent = 15)
                                   #   f <- ph$generate_phreeqc_input_tidal_equil(n_orbital_cycles = 2, porewater_temp = 100)
                                   generate_phreeqc_input_tidal_equil = function(output_dir = "phreeqc_inputs",
                                                                           n_orbital_cycles = NULL,
                                                                           mode = NULL,
                                                                           n_circulations = NULL,
                                                                           grouping_factor = NULL,
                                                                           suffix = "",
                                                                           porewater_temp = "ocean",
                                                                           CH4_redox_override = NULL,
                                                                           carbon_mode = NULL, ch4_cutoff_C = 150,
                                                                           initial_ocean = "pure_water", charge_balance = NULL,
                                                                           pore_gas_fate = c("to_ocean", "retain"), trace_C_molal = NULL,
                                                                           logfile = FALSE) {
                                     pore_gas_fate <- match.arg(pore_gas_fate)
                                     if (is.null(n_orbital_cycles))
                                       stop("generate_phreeqc_input_tidal_equil requires n_orbital_cycles (mode/n_circulations/grouping_factor are unused)")
                                     dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
                                     temp_tag <- if (is.numeric(porewater_temp)) sprintf("T%g", porewater_temp) else sprintf("T%s", porewater_temp)
                                     filename <- file.path(output_dir, sprintf(
                                       "enceladus_k%.2e_d%g_%s_orb%g_equil%s.pqi",
                                       self$simulator$params$k,
                                       self$simulator$params$layer_thickness,
                                       temp_tag, n_orbital_cycles, suffix
                                     ))
                                     content <- self$create_phreeqc_template_tidal_equil(filename, n_orbital_cycles, porewater_temp = porewater_temp,
                                                                                         CH4_redox_override = CH4_redox_override, carbon_mode = carbon_mode, ch4_cutoff_C = ch4_cutoff_C,
                                                                                         initial_ocean = initial_ocean, charge_balance = charge_balance,
                                                                                         pore_gas_fate = pore_gas_fate, trace_C_molal = trace_C_molal)
                                     writeLines(knobs_logfile(self$apply_ssa(content), logfile), filename)
                                     return(filename)
                                   },

                                   create_phreeqc_template_tidal_equil = function(filename, n_orbital_cycles, porewater_temp = "ocean",
                                                                                  CH4_redox_override = NULL,
                                                                                  carbon_mode = NULL, ch4_cutoff_C = 150,
                                                                                  initial_ocean = "pure_water", charge_balance = NULL,
                                                                                  pore_gas_fate = "to_ocean", trace_C_molal = NULL) {
                                     trace_C_line <- trace_C_solution_line(trace_C_molal)
                                     ocean_lines  <- paste0(trace_C_line, initial_ocean_solution_lines(initial_ocean, trace_C_molal, charge_balance))
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
                                     ocean_temp_C <- 274.6 - 273.15  # the profile's ocean floor, top of the rock layer (274.566 K), rounded
                                     # Porewater and ocean conditions from the profile; stops if profile_file is missing.
                                     pc <- profile_conditions(self$profile_file, self$simulator$reactive_thickness())
                                     pore_temp_C        <- pc$pore_temp_C
                                     pore_pressure_atm  <- pc$pore_pressure_atm
                                     ocean_pressure_atm <- pc$ocean_pressure_atm

                                     # Resolve porewater_temp as in create_phreeqc_template_tidal_kinetic()
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
                                     ch4_mode <- self$resolve_ch4_mode(porewater_temp_C, CH4_redox_override, carbon_mode, ch4_cutoff_C)
                                     ch4_override_note <- ch4_mode$note
                                     # Gas phases (see gas_components()): the porewater one, the ocean one, the
                                     # punch columns, and the empty porewater phase to_ocean starts each orbit with.
                                     cometary <- inherits(initial_ocean, "cometary_ice")
                                     # A cometary ocean with a carbon-free porewater failed ("Cl has not converged",
                                     # 2026-10-28, every temperature); with trace_C_molal = 1e-14 it runs.
                                     if (cometary && is.null(trace_C_molal))
                                       warning("cometary_ice_ocean() in tidal_equil: give trace_C_molal = 1e-14 (without it the run failed to converge in tests)")
                                     gas_comps <- gas_components(ch4_mode, if (cometary) "Ntg(g)" else character(0))
                                     gp <- gas_punch(gas_comps, porewater_mass_norm, put = pore_gas_fate == "to_ocean")
                                     pore_gas_block <- gas_phase_block(1, gas_comps, pore_pressure_atm, porewater_temp_C,
                                                                       "Porewater gas (fixed pressure, starts empty)")
                                     ocean_gas_block <- gas_phase_block(2, gas_comps, ocean_pressure_atm, ocean_temp_C,
                                                                        "Ocean gas (fixed pressure, starts empty, never removed)")
                                     # A cometary-ice ocean: the ices react into the pure-water ocean in a simulation of their
                                     # own, before any reactant is defined (cometary_ocean_simulation()), which also defines
                                     # and saves GAS_PHASE 2; otherwise GAS_PHASE 2 is defined empty with the others.
                                     cometary_block <- if (cometary) paste0(cometary_ocean_simulation(initial_ocean, ch4_mode, ocean_mass_norm,
                                                                            ocean_gas_block, base_name, ch4_mode$database), "\n") else ""
                                     ocean_main_block <- if (cometary) "# GAS_PHASE 2 (ocean): defined, and saved with the cometary ices, above" else ocean_gas_block
                                     ch4_db_note_line <- sprintf(
                                       "# Recommended database: %s (CH4 %s; porewater T = %.1f C%s)",
                                       ch4_mode$database,
                                       ch4_mode$description,
                                       porewater_temp_C, ch4_override_note
                                     )

                                     # Temperature-gated secondary minerals per block, as in tidal_kinetic. Each
                                     # block also gets the clathrate phases stable at its reservoir's own
                                     # conditions, gated as in tidal_kinetic (stable_clathrate_phases()).
                                     extra_pw_minerals    <- c(self$temp_gated_minerals(porewater_temp_C),
                                                               self$stable_clathrate_phases(porewater_temp_C, pore_pressure_atm, ch4_mode))
                                     extra_ocean_minerals <- c(self$temp_gated_minerals(ocean_temp_C),
                                                               self$stable_clathrate_phases(ocean_temp_C, ocean_pressure_atm, ch4_mode))
                                     clathrate_note_lines <- paste0(
                                       self$clathrate_stability_note(ocean_temp_C, ocean_pressure_atm, location = "Ocean"),
                                       "\n\n",
                                       self$clathrate_stability_note(porewater_temp_C, pore_pressure_atm, location = "Porewater"))
                                     extra_pw_lines    <- if (length(extra_pw_minerals) > 0)
                                       paste0("   ", extra_pw_minerals, " 0 0", collapse = "\n") else ""
                                     extra_ocean_lines <- if (length(extra_ocean_minerals) > 0)
                                       paste0("   ", extra_ocean_minerals, " 0 0", collapse = "\n") else ""
                                     extra_eq_so_line  <- paste(union(extra_pw_minerals, extra_ocean_minerals), collapse = " ")

                                     # Mineral moles (ROCK_MOL_PER_KG, shared with the other generators; all 11
                                     # are equilibrium here, no kinetics at all). Magnetite is no longer primary.
                                     # SCHREIBERSITE IS NOT INCLUDED in this equilibrium-only generator. The
                                     # kinetic generators add it via KINETICS -formula, but an equilibrium
                                     # phase needs a log K, and neither database has an Fe3P phase. So this
                                     # generator's rock lacks 0.047 mmol P/g and ~0.26 mmol H2/g of reducing
                                     # capacity relative to the kinetic ones. Options (PI decision): add the
                                     # Fe3P elements as an instantaneous REACTION (full dissolution), or add a
                                     # Schreibersite phase with a sourced log K.
                                     mineral_scale_factor <- self$mineral_scale_factor()
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
{ch4_db_note_line}
# pore_gas_fate: {pore_gas_fate}
# Gas phases (fixed pressure, Peng-Robinson): {paste(gas_comps, collapse = " ")}
{clathrate_note_lines}
{EQUIL_INORGANIC_NOTE}{self$mineral_basis_header_line()}
{initial_ocean_header_line(initial_ocean, charge_balance)}
{trace_C_header_line(trace_C_molal)}

KNOBS
   -logfile false

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
# Exchanged fluid mass (per orbital cycle): {sprintf("%.2e", metrics$exchanging_fluid_mass)} kg
{self$simulator$exchange_header_line()}
# Porewater mass: {sprintf("%.2e", self$phreeqc_params$porewater_mass)} kg
# Ocean mass: {sprintf("%.2e", self$phreeqc_params$ocean_mass)} kg

# Initial fluids are defined at 25 C, where pe 4 adds no oxidant or reductant (O(0) and H(0)
# negligible), and are brought to the run temperature by REACTION_TEMPERATURE. Defined at the
# run temperature, pe 4 adds O2: 4.9 molal O(0) at 250 C, and does not converge at 300 C.
SOLUTION 1 Enceladus porewater
    temp      25
    pH        7
    pe        4
    redox     pe
    units     mol/kgw
    density   1
    water     {porewater_mass_norm}
    pressure  {round(pore_pressure_atm, 1)}{trace_C_line}

SOLUTION 2 Enceladus ocean
    temp      25
    pH        7{ocean_pH_charge(charge_balance)}
    pe        4
    redox     pe
    units     mol/kgw
    density   1
    water     {ocean_mass_norm}
    pressure  {round(ocean_pressure_atm, 1)}{ocean_lines}

REACTION_TEMPERATURE 1
    {round(porewater_temp_C, 2)}

REACTION_TEMPERATURE 2
    {round(ocean_temp_C, 2)}

# Reaction pressures, used by every step like the reaction temperatures. MIX averages the
# solution pressures by mass, so without these the porewater would drift towards the ocean
# pressure after each exchange.
REACTION_PRESSURE 1
    {round(pore_pressure_atm, 1)}

REACTION_PRESSURE 2
    {round(ocean_pressure_atm, 1)}

{cometary_block}{pore_gas_block}

{ocean_main_block}
{if (pore_gas_fate == "to_ocean") paste0("\n", gas_transfer_block(gas_comps)) else ""}

EQUILIBRIUM_PHASES 1
# --- Primary minerals (dissolve to equilibrium each cycle). Fe and Ni corrode completely in the first step: UNSOURCED instantaneous corrosion, README section 6, decision D6 ---
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
   Calcite 0 0
#   Celadonite 0 0
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
#   H2O(g) 0 0
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
   Mesolite 0 0
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
   Natrolite 0 0
#   Na2CO3 0 0
#   Na2CO3:7H2O 0 0
#   Nahcolite 0 0
#   Natron 0 0
#   Nesquehonite 0 0
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
#   Stilbite 0 0
#   Strengite 0 0
#   Sylvite 0 0
   Thenardite 0 0
#   Thermonatrite 0 0
   Tobermorite-11A 0 0
{extra_pw_lines}

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
   Calcite 0 0
#   Celadonite 0 0
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
#   H2O(g) 0 0
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
   Mesolite 0 0
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
   Natrolite 0 0
#   Na2CO3 0 0
#   Na2CO3:7H2O 0 0
#   Nahcolite 0 0
#   Natron 0 0
#   Nesquehonite 0 0
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
    -high_precision true
    -reaction true
    -state True
    -solution True
    -step                 true
    -ph                   true
    -pe                   true
    -alkalinity           true
    -ionic_strength       true
    -water                true
{totals_block(ch4_mode$database)}
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
                         Ni Lawrencite Analcime Anhydrite Aragonite Bassanite
                         Beidellite-Ca Beidellite-Fe Beidellite-Mg Beidellite-Na
                         Boehmite Brucite
                         Calcite Chalcedony Chamosite Chrysotile
                         Citric_Acid Clinochlore-14A Clinochlore-7A Clinoptilolite-Ca Clinoptilolite-Na
                         Clinozoisite Cronstedtite-7A
                         Daphnite-14A Daphnite-7A Dawsonite Diaspore Dolomite
                         Ettringite Fe(OH)2 Fe(OH)3 Fe2(SO4)3 FeSO4
                         Gibbsite Glycine Goethite Greenalite Gypsum Gyrolite
                         Hematite Huntite Hydromagnesite Kaolinite
                         Lizardite Magnetite Melanterite Mesolite MgSO4 Minnesotaite Mirabilite
                         Montmor-Ca Montmor-Mg Montmor-Na
                         Natrolite Nontronite-Ca Nontronite-Mg Nontronite-Na
                         Portlandite Pyrite
                         Saponite-Fe-Ca Saponite-Fe-Fe Saponite-Fe-Mg Saponite-Fe-Na
                         Saponite-Mg-Ca Saponite-Mg-Fe Saponite-Mg-Mg Saponite-Mg-Na
                         Sepiolite SiO2(am) Thenardite Tobermorite-11A
                         {extra_eq_so_line}
    -saturation_indices  Enstatite Forsterite Troilite Albite Diopside Anorthite K-Feldspar Tephroite Fe Ni Lawrencite
{gas_so_line(gas_comps)}

# -high_precision (above) also sets the mole-balance convergence criterion to 1e-12; with a gas phase
# that made the first porewater step (the metal pulse) run for more than 15 minutes. This KNOBS, read
# after it, puts back the default 1e-8 (the last setting read wins, PHREEQC read.cpp), and keeps the
# 12-digit output.
KNOBS
    -convergence_tolerance 1e-8

# USER_PUNCH is defined per-cycle below (with exact cumulative time embedded)
    ')

                                     cycle_blocks <- paste(sapply(seq_len(kinetic_cycles), function(i) {
                                       cumulative_time <- i * time_step_years
                                       paste0(
                                         sprintf("\n# --- Cycle %d of %d (cumulative time: %.4e years) ---\n", i, kinetic_cycles, cumulative_time),
                                         # Simulation 1: equilibrate porewater with EQ phases (primary + secondary)
                                         # to_ocean: the porewater gas phase is redefined empty at the start of each orbit
                                         if (pore_gas_fate == "to_ocean" && i > 1) paste0(pore_gas_block, "\n") else "",
                                         "USE solution 1\n",
                                         "USE kinetics none\n",
                                         "USE equilibrium_phases 1\n",
                                         "USE gas_phase 1\n",
                                         "USE reaction none\n",
                                         "USE reaction_temperature 1\n",
                                         "USE reaction_pressure 1\n",
                                         "SAVE solution 1\n",
                                         "SAVE equilibrium_phases 1\n",
                                         if (pore_gas_fate == "retain") "SAVE gas_phase 1\n" else "",
                                         "USER_PUNCH\n",
                                         "   -headings Time_Years Enst_remain_g Forst_remain_g Troil_remain_g Alb_remain_g Diop_remain_g Anorth_remain_g Kfs_remain_g Teph_remain_g Total_Initial_g", SHARED_PUNCH_HEADINGS, gp$headings, "\n",
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
                                         SHARED_PUNCH_LINES, gp$lines,
                                         "END\n\n",
                                         # Simulation 2: ocean absorbs flowing fluid directly + seafloor equilibration
                                         # SOLUTION 4: f_fluid/f_pore of porewater + ocean, mass = f_ocean + f_fluid
                                         OCEAN_STEP_KNOBS_ON,
                                         if (pore_gas_fate == "to_ocean") "USE kinetics 3\n" else "USE kinetics none\n",  # 3: the pore bubble
                                         "USE equilibrium_phases 2\n",
                                         "USE gas_phase 2\n",
                                         "USE reaction none\n",
                                         "USE reaction_temperature 2\n",
                                         "USE reaction_pressure 2\n",
                                         "MIX 4\n",
                                         sprintf("  1  %.15e\n", f_fluid / f_pore),
                                         "  2  1.0\n",
                                         "SAVE solution 4\n",
                                         "SAVE equilibrium_phases 2\n",
                                         "SAVE gas_phase 2\n",
                                         "END\n\n",
                                         OCEAN_STEP_KNOBS_OFF,
                                         # Simulation 3: resample porewater
                                         "USE kinetics none\n",
                                         "USE equilibrium_phases none\n",
                                         "USE reaction none\n",
                                         "USE reaction_temperature 1\n",
                                         "USE reaction_pressure 1\n",
                                         "MIX 1\n",
                                         sprintf("  1  %.15e\n", (f_pore - f_fluid) / f_pore),
                                         sprintf("  4  %.15e\n", f_fluid / (f_ocean + f_fluid)),
                                         "SAVE solution 1\n",
                                         "END\n\n",
                                         # Simulation 4: update ocean
                                         "USE kinetics none\n",
                                         "USE equilibrium_phases none\n",
                                         "USE reaction none\n",
                                         "MIX 2\n",
                                         sprintf("  4  %.15e\n", f_ocean / (f_ocean + f_fluid)),
                                         "SAVE solution 2\n",
                                         "END\n"
                                       )
                                     }), collapse = "\n")

                                     return(self$stamp_pqi_header(self$apply_ch4_mode(paste0(template, "\n", cycle_blocks), ch4_mode),
                                                                  ch4_mode$database, list(pore = pore_pressure_atm, ocean = ocean_pressure_atm)))
                                   },

                                   # Single-step bulk equilibrium: entire hydrosphere (ocean + porewater)
                                   # equilibrated against the full permeable rock layer in one PHREEQC step.
                                   # No permeability or fluid-flow parameters are used.
                                   # An inorganic limiting case: no IOM and no Schreibersite (see
                                   # EQUIL_INORGANIC_NOTE, written into the .pqi header).
                                   # Single-step equilibrium of one solution with the whole rock (inorganic limiting case: no IOM, no
                                   # schreibersite), with one fixed-pressure gas phase.
                                   #
                                   # @param output_dir Directory for the .pqi (and, when run there, its .csv/.out). Default "phreeqc_inputs".
                                   # @param mode "porewater" (porewater only against 1 kg rock, at the porewater pressure) or "hydrosphere" (ocean + porewater combined, at the ocean-floor pressure).
                                   # @param water_rock_ratio kg water per kg rock (mode "porewater"; required there).
                                   # @param water_temp As porewater_temp, for the single solution.
                                   # @param suffix Appended to the file name. Default "".
                                   # @param CH4_redox_override Older spelling of carbon_mode: "coupled" or "decoupled" (= "ch4_uncoupled"). Default NULL; give this or carbon_mode, not both.
                                   # @param carbon_mode "c4_only" (NULL, the default: coupled carbon is C(+4) only, decided by the PI 2026-10-27), "coupled" (an upper bound on abiotic methanation), "ch4_uncoupled" (diagnostic), or "temperature" (the pre-10-26 rule). Picks the database.
                                   # @param ch4_cutoff_C Cutoff (deg C) for carbon_mode = "temperature" only. Default 150.
                                   # @return The path of the .pqi written.
                                   # @examples
                                   #   sim <- Simulator$new(k = 1e-7, layer_thickness = 194200, config = EnceladusConfig$new()); sim$run_simulation()
                                   #   ph <- PhreeqcIntegrator$new(sim); ph$calculate_phreeqc_parameters(organic_wt_percent = 15)
                                   #   f <- ph$generate_phreeqc_input_bulk_equil(mode = "porewater", water_rock_ratio = 1, water_temp = 100)
                                   generate_phreeqc_input_bulk_equil = function(output_dir = "phreeqc_inputs",
                                                                                 mode             = c("hydrosphere", "porewater"),
                                                                                 water_rock_ratio = NULL,
                                                                                 water_temp       = "ocean",
                                                                                 suffix = "",
                                                                                 CH4_redox_override = NULL, carbon_mode = NULL, ch4_cutoff_C = 150,
                                                                                 logfile = FALSE) {
                                     mode <- match.arg(mode)
                                     if (mode == "porewater" && (is.null(water_rock_ratio) || water_rock_ratio <= 0))
                                       stop("water_rock_ratio must be a positive number when mode = 'porewater'")

                                     dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

                                     temp_tag <- if (is.numeric(water_temp)) sprintf("T%g", water_temp) else sprintf("T%s", water_temp)

                                     filename <- if (mode == "porewater") {
                                       file.path(output_dir, sprintf(
                                         "enceladus_wr%.4g_%s_porewater_equil%s.pqi",
                                         water_rock_ratio, temp_tag, suffix
                                       ))
                                     } else {
                                       file.path(output_dir, sprintf(
                                         "enceladus_d%g_%s_hydrosphere_equil%s.pqi",
                                         self$simulator$params$layer_thickness, temp_tag, suffix
                                       ))
                                     }

                                     # Combined hydrosphere mass, normalised to layer mass (hydrosphere mode only)
                                     if (mode == "hydrosphere") {
                                       porewater_mass_norm <- self$phreeqc_params$porewater_mass / self$phreeqc_params$layer_mass
                                       ocean_mass_norm     <- self$phreeqc_params$water_mass   # already = ocean_mass / layer_mass
                                       combined_mass_norm  <- porewater_mass_norm + ocean_mass_norm
                                     }
                                     .rock_kg    <- if (mode == "hydrosphere") self$phreeqc_params$total_rock else 1.0
                                     .water_mass <- if (mode == "porewater") water_rock_ratio else combined_mass_norm

                                     # Temperature and pressure (ocean floor, or profile-averaged over the
                                     # layer if water_temp = "profile")
                                     ocean_temp_C <- 274.6 - 273.15  # the profile's ocean floor, top of the rock layer (274.566 K), rounded
                                     # Porewater and ocean conditions from the profile; stops if profile_file is missing.
                                     pc <- profile_conditions(self$profile_file, self$simulator$reactive_thickness())
                                     pore_temp_C        <- pc$pore_temp_C
                                     pore_pressure_atm  <- pc$pore_pressure_atm
                                     ocean_pressure_atm <- pc$ocean_pressure_atm

                                     # Resolve water_temp to a concrete value — this is the reaction
                                     # temperature for the single combined solution used here.
                                     water_temp_C <- if (is.numeric(water_temp)) {
                                       water_temp
                                     } else if (water_temp == "profile") {
                                       pore_temp_C
                                     } else if (water_temp == "ocean") {
                                       ocean_temp_C
                                     } else {
                                       stop('water_temp must be a number (°C), "profile", or "ocean"')
                                     }

                                     # CH4 (redox-coupled) vs Mtg (decoupled) and the database — see
                                     # resolve_ch4_mode(); water_temp_C plays the role of porewater T.
                                     ch4_mode <- self$resolve_ch4_mode(water_temp_C, CH4_redox_override, carbon_mode, ch4_cutoff_C)
                                     ch4_override_note <- ch4_mode$note
                                     # The gas phase (see gas_components()).
                                     gas_comps <- gas_components(ch4_mode)
                                     ch4_db_note_line <- sprintf(
                                       "# Recommended database: %s (CH4 %s; water T = %.1f C%s)",
                                       ch4_mode$database,
                                       ch4_mode$description,
                                       water_temp_C, ch4_override_note
                                     )

                                     # Temperature-gated secondary minerals, as in bulk_kinetic (single
                                     # solution, so water_temp_C gates the whole block), plus the clathrate
                                     # phases stable at the solution's temperature and pressure, gated as in
                                     # bulk_kinetic.
                                     # Reaction pressure: porewater mode reacts porewater only, so it uses the
                                     # layer-mean porewater pressure; hydrosphere mode keeps the ocean-floor value.
                                     .solution_pressure_atm <- if (mode == "porewater") pore_pressure_atm else ocean_pressure_atm
                                     extra_gated_minerals <- c(self$temp_gated_minerals(water_temp_C),
                                                               self$stable_clathrate_phases(water_temp_C, .solution_pressure_atm,
                                                                                            ch4_mode))
                                     clathrate_note_lines <- self$clathrate_stability_note(water_temp_C, .solution_pressure_atm, location = "Solution")
                                     extra_gated_lines <- if (length(extra_gated_minerals) > 0)
                                       paste0("   ", extra_gated_minerals, " 0 0", collapse = "\n") else ""
                                     extra_gated_so_line <- paste(extra_gated_minerals, collapse = " ")

                                     # Primary mineral moles (ROCK_MOL_PER_KG, shared with the other generators;
                                     # all 11 are equilibrium here since this function has no kinetics at all).
                                     # Magnetite is no longer primary — ordinary secondary phase only.
                                     # SCHREIBERSITE IS NOT INCLUDED here: no Fe3P phase/log K in either
                                     # database, so it can't be an equilibrium primary. See the matching
                                     # note in generate_phreeqc_input_tidal_equil.
                                     mineral_scale_factor <- self$mineral_scale_factor()
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

                                     total_initial_expr <- sprintf(
                                       "%.6f*100.3725 + %.6f*140.6715 + %.6f*87.913 + %.6f*262.1798 + %.6f*216.55 + %.6f*278.164 + %.6f*278.33 + %.6f*201.96 + %.6f*55.845 + %.6f*58.693 + %.6f*126.75",
                                       enstatite_moles, forsterite_moles, troilite_moles, albite_moles,
                                       diopside_moles, anorthite_moles, kfeldspar_moles, tephroite_moles,
                                       fe_moles, ni_moles, lawrencite_moles
                                     )

                                     # ---- mode-specific labels and header ----
                                     .pqi_title  <- if (mode == "porewater") "Enceladus porewater equilibrium" else
                                                      "Enceladus bulk hydrosphere equilibrium"
                                     .soln_label <- if (mode == "porewater") "Enceladus porewater" else
                                                      "Enceladus hydrosphere (ocean + porewater combined)"
                                     .header_lines <- if (mode == "porewater") {
                                       paste0(
                                         "# Mode: single-step equilibrium — porewater only vs permeable rock layer (no ocean)\n",
                                         "# water_rock_ratio: ", water_rock_ratio, "\n",
                                         "# Rock normalised to 1.0 kg; water = ", water_rock_ratio, " kg (porewater only)\n",
                                         "# Temperature: ", round(water_temp_C, 2), " C"
                                       )
                                     } else {
                                       paste0(
                                         "# thickness=", self$simulator$params$layer_thickness, " m\n",
                                         "# Mode: single-step equilibrium — full hydrosphere (ocean + porewater) vs permeable rock layer\n",
                                         "# Layer mass: ", sprintf("%.2e", self$phreeqc_params$layer_mass),
                                         " kg  (total_rock normalised to 1.0)\n",
                                         "# Porewater mass: ", sprintf("%.2e", self$phreeqc_params$porewater_mass), " kg\n",
                                         "# Ocean mass: ", sprintf("%.2e", self$phreeqc_params$ocean_mass), " kg\n",
                                         "# Combined hydrosphere (normalised): ",
                                         sprintf("%.6f", combined_mass_norm), "\n",
                                         "# Temperature: ", round(water_temp_C, 2), " C"
                                       )
                                     }

                                     base_name <- tools::file_path_sans_ext(basename(filename))

                                     content <- glue::glue('
TITLE {.pqi_title}
{ch4_db_note_line}
# Gas phase (fixed pressure, Peng-Robinson): {paste(gas_comps, collapse = " ")}
{clathrate_note_lines}
{EQUIL_INORGANIC_NOTE}{self$mineral_basis_header_line()}

KNOBS
   -logfile false

{.header_lines}

# Initial fluids are defined at 25 C, where pe 4 adds no oxidant or reductant (O(0) and H(0)
# negligible), and are brought to the run temperature by REACTION_TEMPERATURE. Defined at the
# run temperature, pe 4 adds O2: 4.9 molal O(0) at 250 C, and does not converge at 300 C.
SOLUTION 1 {.soln_label}
    temp      25
    pH        7
    pe        4
    redox     pe
    units     mol/kgw
    density   1
    water     {.water_mass}
    pressure  {round(.solution_pressure_atm, 1)}

REACTION_TEMPERATURE 1
    {round(water_temp_C, 2)}

{gas_phase_block(1, gas_comps, .solution_pressure_atm, water_temp_C, "Gas (fixed pressure, starts empty)")}

EQUILIBRIUM_PHASES 1
# --- Primary minerals (dissolve to equilibrium). Fe and Ni corrode completely in the first step: UNSOURCED instantaneous corrosion, README section 6, decision D6 ---
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
   Calcite 0 0
   Chalcedony 0 0
   Chamosite 0 0
   Chrysotile 0 0
   Citric_Acid 0 0
   Clinochlore-14A 0 0
   Clinochlore-7A 0 0
   Clinoptilolite-Ca 0 0
   Clinoptilolite-Na 0 0
   Clinozoisite 0 0
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
   Hematite 0 0
   Huntite 0 0
   Hydromagnesite 0 0
   Kaolinite 0 0
   Lizardite 0 0
   Magnetite 0 0
   Melanterite 0 0
   Mesolite 0 0
   MgSO4 0 0
   Minnesotaite 0 0
   Mirabilite 0 0
   Montmor-Ca 0 0
   Montmor-Mg 0 0
   Montmor-Na 0 0
   Natrolite 0 0
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
   Thenardite 0 0
   Tobermorite-11A 0 0
{extra_gated_lines}

SELECTED_OUTPUT
    -file {base_name}.csv
    -reset false
    -high_precision true
    -state true
    -solution true
    -step true
    -ph true
    -pe true
    -alkalinity true
    -ionic_strength true
    -water true
{totals_block(ch4_mode$database)}
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
                         Ni Lawrencite
                         Analcime Anhydrite Aragonite Bassanite Beidellite-Ca Beidellite-Fe Beidellite-Mg Beidellite-Na
                         Boehmite Brucite Calcite Chalcedony Chamosite Chrysotile
                         Citric_Acid Clinochlore-14A Clinochlore-7A Clinoptilolite-Ca Clinoptilolite-Na
                         Clinozoisite Cronstedtite-7A
                         Daphnite-14A Daphnite-7A Dawsonite Diaspore Dolomite
                         Ettringite Fe(OH)2 Fe(OH)3 Fe2(SO4)3 FeSO4 Gibbsite Glycine Goethite Greenalite Gypsum Gyrolite
                         Hematite Huntite Hydromagnesite Kaolinite Lizardite Magnetite
                         Melanterite Mesolite MgSO4 Minnesotaite Mirabilite
                         Montmor-Ca Montmor-Mg Montmor-Na Natrolite Nontronite-Ca Nontronite-Mg Nontronite-Na
                         Portlandite Pyrite
                         Saponite-Fe-Ca Saponite-Fe-Fe Saponite-Fe-Mg Saponite-Fe-Na
                         Saponite-Mg-Ca Saponite-Mg-Fe Saponite-Mg-Mg Saponite-Mg-Na
                         Sepiolite SiO2(am) Thenardite Tobermorite-11A
                         {extra_gated_so_line}
    -saturation_indices  Enstatite Forsterite Troilite Albite Diopside Anorthite K-Feldspar Tephroite Fe Ni Lawrencite
{gas_so_line(gas_comps)}

# -high_precision (above) also sets the mole-balance convergence criterion to 1e-12; with a gas phase
# that made the first porewater step (the metal pulse) run for more than 15 minutes. This KNOBS, read
# after it, puts back the default 1e-8 (the last setting read wins, PHREEQC read.cpp), and keeps the
# 12-digit output.
KNOBS
    -convergence_tolerance 1e-8

USER_PUNCH
   -headings Time_Years Enst_remain_g Forst_remain_g Troil_remain_g Alb_remain_g Diop_remain_g Anorth_remain_g Kfs_remain_g Teph_remain_g Total_Initial_g{SHARED_PUNCH_HEADINGS}{gas_punch(gas_comps, .water_mass)$headings}
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
{SHARED_PUNCH_LINES}{gas_punch(gas_comps, .water_mass)$lines}
END
')

                                     content <- self$stamp_pqi_header(self$apply_ch4_mode(content, ch4_mode), ch4_mode$database,
                                                                      list(pore = pore_pressure_atm, ocean = ocean_pressure_atm))
                                     writeLines(knobs_logfile(self$apply_ssa(content), logfile), filename)
                                     return(filename)
                                   },

                                   # One solution reacting kinetically with the rock (as tidal_kinetic, no ocean) for n_orbital_cycles orbits,
                                   # with one fixed-pressure gas phase.
                                   #
                                   # @param output_dir Directory for the .pqi (and, when run there, its .csv/.out). Default "phreeqc_inputs".
                                   # @param n_orbital_cycles Number of orbits (1.375 d each).
                                   # @param mode "porewater" (porewater only against 1 kg rock, at the porewater pressure) or "hydrosphere" (ocean + porewater combined, at the ocean-floor pressure).
                                   # @param water_rock_ratio kg water per kg rock (mode "porewater"; required there).
                                   # @param porewater_temp Porewater reaction temperature: a number (deg C), "profile" (the profile's layer mean) or "ocean" (1.45 C, the profile's ocean floor). Default "ocean". Production: 150 and 160 C.
                                   # @param suffix Appended to the file name. Default "".
                                   # @param CH4_redox_override Older spelling of carbon_mode: "coupled" or "decoupled" (= "ch4_uncoupled"). Default NULL; give this or carbon_mode, not both.
                                   # @param carbon_mode "c4_only" (NULL, the default: coupled carbon is C(+4) only, decided by the PI 2026-10-27), "coupled" (an upper bound on abiotic methanation), "ch4_uncoupled" (diagnostic), or "temperature" (the pre-10-26 rule). Picks the database.
                                   # @param ch4_cutoff_C Cutoff (deg C) for carbon_mode = "temperature" only. Default 150.
                                   # @param rec_eqb_minerals Promote primaries whose half-life (x rec_eqb_margin) is under one orbit to equilibrium. Default FALSE (a modelling choice; not needed for convergence).
                                   # @param rec_eqb_margin Safety margin for rec_eqb_minerals. Default 2.
                                   # @param trace_C_molal Numerical C seed in every initial SOLUTION, mol/kgw. Default NULL; production 1e-14 (the smallest tested to converge at 25 C).
                                   # @return The path of the .pqi written.
                                   # @examples
                                   #   sim <- Simulator$new(k = 1e-7, layer_thickness = 194200, config = EnceladusConfig$new()); sim$run_simulation()
                                   #   ph <- PhreeqcIntegrator$new(sim); ph$calculate_phreeqc_parameters(organic_wt_percent = 15)
                                   #   f <- ph$generate_phreeqc_input_bulk_kinetic(n_orbital_cycles = 2, mode = "porewater", water_rock_ratio = 1, porewater_temp = 100)
                                   generate_phreeqc_input_bulk_kinetic = function(
                                     output_dir          = "phreeqc_inputs",
                                     n_orbital_cycles    = NULL,
                                     mode                = c("porewater", "hydrosphere"),
                                     water_rock_ratio    = NULL,
                                     porewater_temp      = "ocean",
                                     suffix              = "",
                                     CH4_redox_override  = NULL,
                                     carbon_mode         = NULL,
                                     ch4_cutoff_C        = 150,
                                     rec_eqb_minerals    = FALSE,
                                     rec_eqb_margin      = 2,
                                     trace_C_molal       = NULL,
                                     logfile             = FALSE
                                   ) {
                                     mode <- match.arg(mode)
                                     trace_C_line <- trace_C_solution_line(trace_C_molal)
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

                                     # Temperature and pressure
                                     ocean_temp_C <- 274.6 - 273.15  # the profile's ocean floor, top of the rock layer (274.566 K), rounded
                                     # Porewater and ocean conditions from the profile; stops if profile_file is missing.
                                     pc <- profile_conditions(self$profile_file, self$simulator$reactive_thickness())
                                     pore_temp_C        <- pc$pore_temp_C
                                     pore_pressure_atm  <- pc$pore_pressure_atm
                                     ocean_pressure_atm <- pc$ocean_pressure_atm

                                     # Resolve porewater_temp to a concrete value — this is the reaction
                                     # temperature for the single combined solution used here.
                                     porewater_temp_C <- if (is.numeric(porewater_temp)) porewater_temp else
                                       if (porewater_temp == "profile") pore_temp_C else ocean_temp_C

                                     # CH4 (redox-coupled) vs Mtg (decoupled) — see resolve_ch4_mode().
                                     ch4_mode <- self$resolve_ch4_mode(porewater_temp_C, CH4_redox_override, carbon_mode, ch4_cutoff_C)
                                     ch4_override_note <- ch4_mode$note
                                     # The gas phase (see gas_components()).
                                     gas_comps <- gas_components(ch4_mode)
                                     ch4_db_note_line <- sprintf(
                                       "# Recommended database: %s (CH4 %s; porewater T = %.1f C%s)",
                                       ch4_mode$database,
                                       ch4_mode$description,
                                       porewater_temp_C, ch4_override_note
                                     )

                                     # Reaction pressure of the single solution, as in bulk_equil: porewater mode
                                     # reacts porewater only, so it uses the layer-mean porewater pressure;
                                     # hydrosphere mode (ocean + porewater combined) keeps the ocean-floor value.
                                     .solution_pressure_atm <- if (mode == "porewater") pore_pressure_atm else ocean_pressure_atm
                                     # Clathrate stability at that temperature and pressure — see resolve_clathrate_stability().
                                     clathrate_note_lines <- self$clathrate_stability_note(porewater_temp_C, .solution_pressure_atm)

                                     # Temperature-gated secondary minerals (single solution, so a single
                                     # temperature — porewater_temp_C — gates the whole block). Stable
                                     # clathrate phases (checked against .solution_pressure_atm) are
                                     # appended the same way.
                                     extra_ocean_minerals <- c(self$temp_gated_minerals(porewater_temp_C),
                                                               self$stable_clathrate_phases(porewater_temp_C, .solution_pressure_atm, ch4_mode))
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
                                     mineral_scale_factor <- self$mineral_scale_factor()
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

                                     # ---- rec_eqb_minerals: auto-promote fast-dissolving primaries to equilibrium ----
                                     # Off by default (rec_eqb_minerals = FALSE). When on, a primary whose
                                     # kinetic half-life, times a safety margin (rec_eqb_margin), is shorter
                                     # than one orbital period is moved to dissolve_only equilibrium.
                                     # Reference pH is chosen to be conservative (it under-promotes) for each
                                     # rate law: pH 11 for acid-branch minerals, whose rate is slowest at high
                                     # pH, and pH 9 for the base-branch minerals (Albite, K-Feldspar), whose
                                     # rate is slowest at the low end of the alkaline range seen in runs
                                     # (pH ~9-11.5). At pH 11 the base branch is ~14x (Albite) and ~44x
                                     # (K-Feldspar) faster than at pH 9, so a single pH-11 reference
                                     # over-promoted them. rec_eqb_margin = 2 means promote only if double the
                                     # half-life still fits inside one orbit. Promotion is a modelling choice,
                                     # not a convergence fix: the 10-06 hang was the unguarded tephroite rate
                                     # (CHANGELOG 2026-10-19). Uses k_primary_mineral()
                                     # (mineral_dissolution_constants.r), the same rate laws as the RATES
                                     # blocks below. Schreibersite is excluded: it cannot be an equilibrium
                                     # primary (no Fe3P phase/log K in the database).
                                     all_primary <- c("Enstatite", "Forsterite", "Troilite", "Albite",
                                                      "Diopside", "Anorthite", "K-Feldspar", "Tephroite")
                                     orbital_period <- self$simulator$config$orbital_period
                                     if (rec_eqb_minerals) {
                                       base_branch_minerals <- c("Albite", "K-Feldspar")
                                       ref_pH <- ifelse(all_primary %in% base_branch_minerals, 9, 11)
                                       names(ref_pH) <- all_primary
                                       half_life_s <- sapply(all_primary, function(m) {
                                         k <- k_primary_mineral(m, H = 10^(-ref_pH[[m]]), TK = porewater_temp_C + 273.15)
                                         if (k <= 0) Inf else log(2) / k
                                       })
                                       equil_primary <- names(half_life_s)[rec_eqb_margin * half_life_s < orbital_period]
                                       print(sprintf("rec_eqb_minerals: half-life @ pH 9 (base-branch) / 11 (acid-branch), T=%.1fC, x%g margin vs 1 orbit (%.0f s)",
                                                     porewater_temp_C, rec_eqb_margin, orbital_period))
                                       print(data.frame(mineral = all_primary, ref_pH = unname(ref_pH), half_life_yr = half_life_s / (365.25*24*3600),
                                                        promoted_to_equilibrium = all_primary %in% equil_primary))
                                     } else {
                                       equil_primary <- character(0)
                                     }
                                     kinetic_minerals <- setdiff(all_primary, equil_primary)
                                     rec_eqb_note <- if (!rec_eqb_minerals) {
                                       "# rec_eqb_minerals=FALSE: all primaries kinetic regardless of half-life"
                                     } else if (length(equil_primary) > 0) {
                                       sprintf("# rec_eqb_minerals=TRUE (margin=%gx): promoted to equilibrium (half-life at pH 9 base-branch / 11 acid-branch, T=%.1fC, x%g < 1 orbital period [%.0f s]): %s",
                                               rec_eqb_margin, porewater_temp_C, rec_eqb_margin, orbital_period, paste(equil_primary, collapse = ", "))
                                     } else {
                                       sprintf("# rec_eqb_minerals=TRUE (margin=%gx): no primaries promoted (all half-lives at pH 9 base-branch / 11 acid-branch, T=%.1fC, x%g exceed 1 orbital period [%.0f s])",
                                               rec_eqb_margin, porewater_temp_C, rec_eqb_margin, orbital_period)
                                     }

                                     molar_masses <- c(Enstatite = 100.3725, Forsterite = 140.6715, Troilite = 87.913,
                                                       Albite = 262.1798, Diopside = 216.55, Anorthite = 278.164,
                                                       `K-Feldspar` = 278.33, Tephroite = 201.96)
                                     primary_moles <- c(Enstatite = enstatite_moles, Forsterite = forsterite_moles,
                                                        Troilite = troilite_moles, Albite = albite_moles,
                                                        Diopside = diopside_moles, Anorthite = anorthite_moles,
                                                        `K-Feldspar` = kfeldspar_moles, Tephroite = tephroite_moles)

                                     # ---- RATES block (only for kinetic_minerals; Schreibersite/IOM separate) ----
                                     .rate_laws <- list(
                                       Forsterite = paste(
                                         '\t-start',
                                         '\t1   REM Ref PK04',
                                         '\t10  kacid = 10^(-6.85) * exp(-67.2e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^0.47',
                                         '\t20  kneut = 10^(-10.64) * exp(-79.0e3/8.314 * (1/TK-1/298.15))',
                                         '\t21  SSA = 0.1  # UNSOURCED surface area (m2/g): README section 6, decision D7',
                                         '\t22  mw = 140.6715',
                                         '\t40  k = (kacid + kneut) * SSA * mw * M',
                                         '\t50  IF SR("Forsterite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Forsterite"))',
                                         '\t60  moles = rate * TIME',
                                         '\t70  SAVE moles',
                                         '\t-end', sep = "\n"),
                                       Troilite = paste(
                                         '\t-start',
                                         '\t1   REM Acid mechanism only (PK04, hexagonal pyrrhotite formulation)',
                                         '\t25  kacid = 10^(-6.79) * exp(-63.0e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^-0.090 * ACT("Fe+3")^0.356',
                                         '\t26  SSA = 5  # UNSOURCED surface area (m2/g): README section 6, decision D7',
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
                                         '\t21  SSA = 0.1  # UNSOURCED surface area (m2/g): README section 6, decision D7',
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
                                         '\t31  SSA = 5  # UNSOURCED surface area (m2/g): README section 6, decision D7',
                                         '\t32  mw = 278.33',
                                         '\t40  k = (kacid + kneut + kbase) * SSA * mw * M',
                                         '\t50  IF SR("K-Feldspar") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("K-Feldspar"))',
                                         '\t60  moles = rate * TIME',
                                         '\t70  SAVE moles',
                                         '\t-end', sep = "\n"),
                                       Tephroite = paste(
                                         '\t-start',
                                         '\t1   REM Pincus et al. (2026) ACS Earth Space Chem 10, 1174-1184, Table S1,',
                                         '\t2   REM citing Casey et al. (1993) GCA 57, 785-793. See generate_phreeqc_input_tidal_kinetic',
                                         '\t3   REM for full derivation notes -- identical rate law, ported here verbatim.',
                                         '\t34  REM pH >= 7: acid branch evaluated at a(H+) = 1e-7 (flat at the pH-7 rate). No neutral term is summed.',
                                         '\t35  aH = ACT("H+")',
                                         '\t36  IF aH < 1e-7 THEN aH = 1e-7',
                                         '\t37  Ea = -60.11*LOG(-LOG10(aH)) + 106.17',
                                         '\t38  IF Ea < 0 THEN Ea = 0',
                                         '\t40  kacid = 10^(-4.46) * exp(-Ea*1e3/8.314 * (1/TK-1/298.15)) * aH^0.47',
                                         '\t46  SSA = 0.1  # UNSOURCED surface area (m2/g): README section 6, decision D7',
                                         '\t47  mw = 201.96',
                                         '\t50  k = kacid * SSA * mw * M',
                                         '\t60  IF SR("Tephroite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Tephroite"))',
                                         '\t70  moles = rate * TIME',
                                         '\t80  SAVE moles',
                                         '\t-end', sep = "\n"),
                                       Enstatite = paste(
                                         '\t-start',
                                         '\t1   REM Ref PK04',
                                         '\t10  kacid = 10^(-9.02) * exp(-80.0e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^0.6',
                                         '\t20  kneut = 10^(-12.72) * exp(-80.0e3/8.314 * (1/TK-1/298.15))',
                                         '\t21  SSA = 0.1  # UNSOURCED surface area (m2/g): README section 6, decision D7',
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
                                         '\t21  SSA = 5  # UNSOURCED surface area (m2/g): README section 6, decision D7',
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
                                         '\t31  SSA = 5  # UNSOURCED surface area (m2/g): README section 6, decision D7',
                                         '\t32  mw = 262.1798',
                                         '\t40  k = (kacid + kneut + kbase) * SSA * mw * M',
                                         '\t50  IF SR("Albite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Albite"))',
                                         '\t60  moles = rate * TIME',
                                         '\t70  SAVE moles',
                                         '\t-end', sep = "\n")
                                     )

                                     .primary_rate_defs <- if (length(kinetic_minerals) > 0)
                                       paste(sapply(kinetic_minerals, function(m) paste0(m, "\n", .rate_laws[[m]])), collapse = "\n\n")
                                     else ""

                                     # ---- KINETICS entries (only for kinetic_minerals) ----
                                     .primary_kin_entries <- if (length(kinetic_minerals) > 0)
                                       paste(sapply(kinetic_minerals, function(m)
                                         sprintf("   %s\n      -m0 %.6e\n      -step_divide 1", m, primary_moles[[m]])),
                                         collapse = "\n")
                                     else ""

                                     # ---- EQUILIBRIUM_PHASES primary lines: Fe/Ni/Lawrencite (always) + equil_primary (promoted) ----
                                     equil_prim_lines <- paste(
                                       c("# --- Primary minerals (equilibrium; dissolve_only so they can only deplete). Fe and Ni corrode completely
#     in the first step: UNSOURCED instantaneous corrosion, README section 6, decision D6 ---",
                                         sprintf("   Fe                   0  %.6e          dissolve_only", fe_moles),
                                         sprintf("   Ni                   0  %.6e          dissolve_only", ni_moles),
                                         sprintf("   Lawrencite           0  %.6e          dissolve_only", lawrencite_moles),
                                         if (length(equil_primary) > 0)
                                           sapply(equil_primary, function(m)
                                             sprintf("   %-20s 0  %.6e  dissolve_only", m, primary_moles[[m]]))
                                         else character(0)),
                                       collapse = "\n"
                                     )

                                     # ---- SELECTED_OUTPUT dynamic lines ----
                                     kin_so_line <- paste(c(kinetic_minerals, "Schreibersite"), collapse = " ")
                                     equil_so_line <- paste(c("Ni", "Lawrencite", equil_primary), collapse = " ")

                                     # ---- USER_PUNCH: EQUI() for promoted minerals, KIN() otherwise ----
                                     punch_lines <- paste(
                                       sapply(seq_along(all_primary), function(i) {
                                         m  <- all_primary[i]
                                         fn <- if (m %in% equil_primary) sprintf('EQUI("%s")', m) else sprintf('KIN("%s")', m)
                                         sprintf("   %d PUNCH %s * %.4f", 10 + i * 10, fn, molar_masses[[m]])
                                       }),
                                       collapse = "\n"
                                     )

                                     # ---- IOM as kinetic phases (see generate_phreeqc_input_tidal_kinetic
                                     # for full derivation notes; identical wiring, generic over iom_pools) --
                                     iom_rates_block    <- ""
                                     iom_kinetics_block <- ""
                                     if (self$organic_wt_percent > 0) {
                                       iom_mass_kg_scaled <- (self$organic_wt_percent / 100) * .rock_kg
                                       iom_total_moles    <- iom_mass_kg_scaled * 1.0   # 1 mol per kg by construction
                                       self$resolve_iom_pools(ch4_mode)

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

                                     # ---- PATCH D1: punch IOM channel state so the restart can recover it ----
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
{trace_C_header_line(trace_C_molal)}{self$mineral_basis_header_line()}
{ch4_db_note_line}
# Gas phase (fixed pressure, Peng-Robinson): {paste(gas_comps, collapse = " ")}
{clathrate_note_lines}

KNOBS
   -logfile false

{.header_lines}
{rec_eqb_note}

# Initial fluids are defined at 25 C, where pe 4 adds no oxidant or reductant (O(0) and H(0)
# negligible), and are brought to the run temperature by REACTION_TEMPERATURE. Defined at the
# run temperature, pe 4 adds O2: 4.9 molal O(0) at 250 C, and does not converge at 300 C.
SOLUTION 1 {.soln_label}
    temp      25
    pH        7
    pe        4
    redox     pe
    units     mol/kgw
    density   1
    water     {.water_mass}
    pressure  {round(.solution_pressure_atm, 1)}{trace_C_line}

REACTION_TEMPERATURE 1
    {round(porewater_temp_C, 2)}

{gas_phase_block(1, gas_comps, .solution_pressure_atm, porewater_temp_C, "Gas (fixed pressure, starts empty)")}

INCREMENTAL_REACTIONS true

RATES
{.primary_rate_defs}

{SCHREIBERSITE_RATE}

{iom_rates_block}

KINETICS 1
{.primary_kin_entries}
{schreibersite_kinetics_entry(schreibersite_moles)}
   {iom_kinetics_block}
   -steps {sprintf("%.6e", total_time_seconds)} in {n_orbital_cycles} steps
   -cvode true
   -bad_step_max 5000

EQUILIBRIUM_PHASES 1
{equil_prim_lines}
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
   Calcite 0 0
   Chalcedony 0 0
   Chamosite 0 0
   Chrysotile 0 0
   Citric_Acid 0 0
   Clinochlore-14A 0 0
   Clinochlore-7A 0 0
   Clinoptilolite-Ca 0 0
   Clinoptilolite-Na 0 0
   Clinozoisite 0 0
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
   Hematite 0 0
   Huntite 0 0
   Hydromagnesite 0 0
   Kaolinite 0 0
   Lizardite 0 0
   Magnetite 0 0
   Melanterite 0 0
   Mesolite 0 0
   MgSO4 0 0
   Minnesotaite 0 0
   Mirabilite 0 0
   Montmor-Ca 0 0
   Montmor-Mg 0 0
   Montmor-Na 0 0
   Natrolite 0 0
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
   Thenardite 0 0
   Tobermorite-11A 0 0
{extra_ocean_lines}

SELECTED_OUTPUT
    -file {base_name}.csv
    -reset false
    -high_precision true
    -state true
    -solution true
    -step true
    -ph true
    -pe true
    -alkalinity true
    -ionic_strength true
    -water true
{totals_block(ch4_mode$database)}
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
    -kinetics            {kin_so_line} {iom_so_line}
    -saturation_indices  Enstatite Forsterite Troilite Albite Diopside Anorthite K-Feldspar Tephroite
    -equilibrium_phases  {equil_so_line} Analcime Anhydrite Aragonite Bassanite
                         Beidellite-Ca Beidellite-Fe Beidellite-Mg Beidellite-Na
                         Boehmite Brucite
                         Calcite Chalcedony Chamosite Chrysotile
                         Citric_Acid Clinochlore-14A Clinochlore-7A Clinoptilolite-Ca Clinoptilolite-Na
                         Clinozoisite Cronstedtite-7A
                         Daphnite-14A Daphnite-7A Dawsonite Diaspore Dolomite
                         Ettringite Fe(OH)2 Fe(OH)3 Fe2(SO4)3 FeSO4
                         Gibbsite Glycine Goethite Greenalite Gypsum Gyrolite
                         Hematite Huntite Hydromagnesite Kaolinite
                         Lizardite Magnetite Melanterite Mesolite MgSO4 Minnesotaite Mirabilite
                         Montmor-Ca Montmor-Mg Montmor-Na
                         Natrolite Nontronite-Ca Nontronite-Mg Nontronite-Na
                         Portlandite Pyrite
                         Saponite-Fe-Ca Saponite-Fe-Fe Saponite-Fe-Mg Saponite-Fe-Na
                         Saponite-Mg-Ca Saponite-Mg-Fe Saponite-Mg-Mg Saponite-Mg-Na
                         Sepiolite SiO2(am) Thenardite Tobermorite-11A
                         {extra_eq_so_line}
{gas_so_line(gas_comps)}

# -high_precision (above) also sets the mole-balance convergence criterion to 1e-12; with a gas phase
# that made the first porewater step (the metal pulse) run for more than 15 minutes. This KNOBS, read
# after it, puts back the default 1e-8 (the last setting read wins, PHREEQC read.cpp), and keeps the
# 12-digit output.
KNOBS
    -convergence_tolerance 1e-8

USER_PUNCH
   -headings Time_Years Enst_remain_g Forst_remain_g Troil_remain_g Alb_remain_g Diop_remain_g Anorth_remain_g Kfs_remain_g Teph_remain_g Schr_remain_g Total_Initial_g{iom_punch_headings}{SHARED_PUNCH_HEADINGS}{gas_punch(gas_comps, .water_mass)$headings}
   10 PUNCH TOTAL_TIME / (365.25 * 24 * 3600)
{punch_lines}
   95 PUNCH KIN("Schreibersite") * 198.509
   100 total_initial = {total_initial_expr}
   110 PUNCH total_initial
{iom_punch_lines}
{SHARED_PUNCH_LINES}{gas_punch(gas_comps, .water_mass)$lines}
END
')

                                     content <- self$stamp_pqi_header(self$apply_ch4_mode(content, ch4_mode), ch4_mode$database,
                                                                      list(pore = pore_pressure_atm, ocean = ocean_pressure_atm))
                                     writeLines(knobs_logfile(self$apply_ssa(content), logfile), filename)
                                     return(filename)
                                   },

                                   # Continues a bulk_kinetic run from the last row of its CSV: every element the parent punched (15 digits),
                                   # its equilibrium phases, kinetic reactants, IOM channels and gas phase.
                                   #
                                   # @param output_dir Directory for the .pqi (and, when run there, its .csv/.out). Default "phreeqc_inputs".
                                   # @param prev_csv The parent run's CSV; its .pqi must be beside it (for the carbon-mode, pressure and fluid checks).
                                   # @param n_orbital_cycles Number of orbits (1.375 d each).
                                   # @param orbits_per_step Orbits per KINETICS step; must divide n_orbital_cycles. Default 1.
                                   # @param equil_minerals Primaries to move from KINETICS to EQUILIBRIUM_PHASES (dissolve_only). Default none; Schreibersite cannot be one.
                                   # @param mode "porewater" (porewater only against 1 kg rock, at the porewater pressure) or "hydrosphere" (ocean + porewater combined, at the ocean-floor pressure).
                                   # @param porewater_temp Porewater reaction temperature: a number (deg C), "profile" (the profile's layer mean) or "ocean" (1.45 C, the profile's ocean floor). Default "ocean". Production: 150 and 160 C.
                                   # @param suffix Appended to the file name. Default "".
                                   # @param CH4_redox_override Older spelling of carbon_mode: "coupled" or "decoupled" (= "ch4_uncoupled"). Default NULL; give this or carbon_mode, not both.
                                   # @param carbon_mode "c4_only" (NULL, the default: coupled carbon is C(+4) only, decided by the PI 2026-10-27), "coupled" (an upper bound on abiotic methanation), "ch4_uncoupled" (diagnostic), or "temperature" (the pre-10-26 rule). Picks the database.
                                   # @param ch4_cutoff_C Cutoff (deg C) for carbon_mode = "temperature" only. Default 150.
                                   # @param allow_pressure_change TRUE lets a restart run at pressures other than its parent's (with a warning). Default FALSE (an error).
                                   # @return The path of the .pqi written.
                                   # @examples
                                   #   f <- ph$generate_phreeqc_input_bulk_kinetic_restart(prev_csv = "phreeqc_inputs/<parent>.csv", n_orbital_cycles = 2, porewater_temp = 100)
                                   generate_phreeqc_input_bulk_kinetic_restart = function(
                                     output_dir          = "phreeqc_inputs",
                                     prev_csv            = NULL,
                                     n_orbital_cycles    = NULL,
                                     orbits_per_step     = 1,
                                     equil_minerals      = character(0),
                                     mode                = c("porewater", "hydrosphere"),
                                     porewater_temp      = "ocean",
                                     suffix              = "",
                                     CH4_redox_override  = NULL,
                                     carbon_mode         = NULL,
                                     ch4_cutoff_C        = 150,
                                     allow_pressure_change = FALSE,
                                     logfile             = FALSE
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

                                     # Every element the parent punched (restart_elements(); the three
                                     # databases define the same elements, so the coupled one is used here,
                                     # before this restart's carbon mode is resolved).
                                     elements  <- restart_elements(file.path("databases", CARBON_MODE_DATABASES[["coupled"]]),
                                                                   names(last_row), prev_csv)
                                     elem_vals <- sapply(elements, function(el) {
                                       val <- suppressWarnings(as.numeric(last_row[[el]]))
                                       if (!is.na(val) && val > 0) val else NA_real_
                                     })
                                     elem_lines <- paste(
                                       Filter(Negate(is.null), lapply(seq_along(elements), function(i) {
                                         if (is.na(elem_vals[i])) return(NULL)
                                         sprintf("    %-4s  %.15e", elements[i], elem_vals[i])
                                       })),
                                       collapse = "\n"
                                     )

                                     # Temperature and pressure
                                     ocean_temp_C <- 274.6 - 273.15  # the profile's ocean floor, top of the rock layer (274.566 K), rounded
                                     # Porewater and ocean conditions from the profile; stops if profile_file is missing.
                                     pc <- profile_conditions(self$profile_file, self$simulator$reactive_thickness())
                                     pore_temp_C        <- pc$pore_temp_C
                                     pore_pressure_atm  <- pc$pore_pressure_atm
                                     ocean_pressure_atm <- pc$ocean_pressure_atm

                                     # Resolve porewater_temp to a concrete value — this is the reaction
                                     # temperature for the single combined solution used here.
                                     porewater_temp_C <- if (is.numeric(porewater_temp)) porewater_temp else
                                       if (porewater_temp == "profile") pore_temp_C else ocean_temp_C

                                     # CH4 (redox-coupled) vs Mtg (decoupled) — see resolve_ch4_mode().
                                     ch4_mode <- self$resolve_ch4_mode(porewater_temp_C, CH4_redox_override, carbon_mode, ch4_cutoff_C)
                                     parent_mode_line <- check_restart_carbon_mode(prev_csv, ch4_mode)
                                     parent_pressure_line <- check_restart_pressures(prev_csv, pore_pressure_atm, ocean_pressure_atm,
                                                                                     allow_pressure_change)
                                     ch4_override_note <- ch4_mode$note
                                     # The gas phase (see gas_components()).
                                     gas_comps <- gas_components(ch4_mode)
                                     ch4_db_note_line <- sprintf(
                                       "# Recommended database: %s (CH4 %s; porewater T = %.1f C%s)",
                                       ch4_mode$database,
                                       ch4_mode$description,
                                       porewater_temp_C, ch4_override_note
                                     )

                                     # Reaction pressure of the single solution, as in bulk_kinetic: the layer-mean
                                     # porewater pressure in porewater mode, the ocean-floor value in hydrosphere mode.
                                     .solution_pressure_atm <- if (mode == "porewater") pore_pressure_atm else ocean_pressure_atm
                                     # Clathrate stability at that temperature and pressure — see resolve_clathrate_stability().
                                     clathrate_note_lines <- self$clathrate_stability_note(porewater_temp_C, .solution_pressure_atm)

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

                                     # Floor for near-zero moles: a kinetic mineral that's dissolved down to
                                     # floating-point dust (e.g. 1e-21 mol) is chemically depleted, but writing
                                     # that literal value into a dissolve_only EQUILIBRIUM_PHASES entry still
                                     # makes PHREEQC try to satisfy its saturation condition with that
                                     # infinitesimal mass -- which can fail to converge if the solution is far
                                     # from saturation for that phase (seen under a low water_rock_ratio restart,
                                     # SI(Albite) = -4.95 but m0 = 1.06e-21: PHREEQC couldn't reconcile "dissolve
                                     # toward SR=1" with "only this much mass exists"; clamping to exactly 0
                                     # removes the phase from the system entirely and converged cleanly). Any
                                     # real mineral mass in this system is many orders of magnitude above this
                                     # threshold, so it only ever catches floating-point dust, not genuine
                                     # small-but-real remaining amounts.
                                     ZERO_MOLES_FLOOR <- 1e-15
                                     get_k_moles <- function(mineral) {
                                       col <- paste0("k_", mineral)
                                       val <- if (col %in% names(last_row)) suppressWarnings(as.numeric(last_row[[col]])) else NA_real_
                                       if (!is.na(val) && val >= 0) return(if (val < ZERO_MOLES_FLOOR) 0 else val)
                                       # Fallback: mineral was equilibrium in prev run — read *_remain_g column
                                       g_col <- remain_g_cols[mineral]
                                       if (!is.na(g_col) && g_col %in% names(last_row)) {
                                         g_val <- suppressWarnings(as.numeric(last_row[[g_col]]))
                                         if (!is.na(g_val) && g_val >= 0) {
                                           moles <- g_val / molar_masses[mineral]
                                           return(if (moles < ZERO_MOLES_FLOOR) 0 else moles)
                                         }
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
                                       self$resolve_iom_pools(ch4_mode)

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
                                         sprintf("   %s\n      -formula %s\n      -m0 %.15e\n      -step_divide 1",
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
                                       "Calcite", 
                                       "Chalcedony", "Chamosite", "Chrysotile", "Citric_Acid",
                                       "Clinochlore-14A", "Clinochlore-7A", "Clinoptilolite-Ca", "Clinoptilolite-Na",
                                       "Clinozoisite", "Cronstedtite-7A",
                                       "Daphnite-14A", "Daphnite-7A", "Dawsonite", "Diaspore", "Dolomite",
                                       "Ettringite", "Fe(OH)2", "Fe(OH)3", "Fe2(SO4)3", "FeSO4",
                                       "Gibbsite", "Glycine", "Goethite", "Greenalite", "Gypsum", "Gyrolite",
                                       "Hematite", "Huntite", "Hydromagnesite",
                                       "Kaolinite", "Lizardite", "Magnetite", "Melanterite", "Mesolite", "MgSO4", "Minnesotaite", "Mirabilite",
                                       "Montmor-Ca", "Montmor-Mg", "Montmor-Na",
                                       "Natrolite", "Nontronite-Ca", "Nontronite-Mg", "Nontronite-Na",
                                       "Portlandite", "Pyrite",
                                       "Saponite-Fe-Ca", "Saponite-Fe-Fe", "Saponite-Fe-Mg", "Saponite-Fe-Na",
                                       "Saponite-Mg-Ca", "Saponite-Mg-Fe", "Saponite-Mg-Mg", "Saponite-Mg-Na",
                                       "Sepiolite", "SiO2(am)", "Thenardite", "Tobermorite-11A"
                                     )
                                     # Single solution (porewater_temp_C) gates the whole assemblage here.
                                     # Clathrate phases stable against .solution_pressure_atm (this function's
                                     # own reaction pressure) are appended the same way.
                                     secondary_minerals <- c(secondary_minerals, self$temp_gated_minerals(porewater_temp_C),
                                                             self$stable_clathrate_phases(porewater_temp_C, .solution_pressure_atm, ch4_mode))

                                     get_equil_moles <- function(mineral) {
                                       col <- ep_column(mineral, names(last_row))
                                       val <- if (!is.na(col)) suppressWarnings(as.numeric(last_row[[col]])) else NA_real_
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
                                         '\t21  SSA = 0.1  # UNSOURCED surface area (m2/g): README section 6, decision D7',
                                         '\t22  mw = 140.6715',
                                         '\t40  k = (kacid + kneut) * SSA * mw * M',
                                         '\t50  IF SR("Forsterite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Forsterite"))',
                                         '\t60  moles = rate * TIME',
                                         '\t70  SAVE moles',
                                         '\t-end', sep = "\n"),
                                       Troilite = paste(
                                         '\t-start',
                                         '\t1   REM PK04 has no separate troilite entry, hexagonal pyrrhotite parameters used as proxy',
                                         '\t25  kacid = 10^(-6.79) * exp(-63.0e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^-0.090 * ACT("Fe+3")^0.356',
                                         '\t26  SSA = 5  # UNSOURCED surface area (m2/g): README section 6, decision D7',
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
                                         '\t21  SSA = 0.1  # UNSOURCED surface area (m2/g): README section 6, decision D7',
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
                                         '\t31  SSA = 5  # UNSOURCED surface area (m2/g): README section 6, decision D7',
                                         '\t32  mw = 278.33',
                                         '\t40  k = (kacid + kneut + kbase) * SSA * mw * M',
                                         '\t50  IF SR("K-Feldspar") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("K-Feldspar"))',
                                         '\t60  moles = rate * TIME',
                                         '\t70  SAVE moles',
                                         '\t-end', sep = "\n"),
                                       Tephroite = paste(
                                         '\t-start',
                                         '\t1   REM Pincus et al. (2026) ACS Earth Space Chem 10, 1174-1184, Table S1,',
                                         '\t2   REM citing Casey et al. (1993) GCA 57, 785-793. See generate_phreeqc_input_tidal_kinetic',
                                         '\t3   REM for full derivation notes -- identical rate law, ported here verbatim.',
                                         '\t34  REM pH >= 7: acid branch evaluated at a(H+) = 1e-7 (flat at the pH-7 rate). No neutral term is summed.',
                                         '\t35  aH = ACT("H+")',
                                         '\t36  IF aH < 1e-7 THEN aH = 1e-7',
                                         '\t37  Ea = -60.11*LOG(-LOG10(aH)) + 106.17',
                                         '\t38  IF Ea < 0 THEN Ea = 0',
                                         '\t40  kacid = 10^(-4.46) * exp(-Ea*1e3/8.314 * (1/TK-1/298.15)) * aH^0.47',
                                         '\t46  SSA = 0.1  # UNSOURCED surface area (m2/g): README section 6, decision D7',
                                         '\t47  mw = 201.96',
                                         '\t50  k = kacid * SSA * mw * M',
                                         '\t60  IF SR("Tephroite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Tephroite"))',
                                         '\t70  moles = rate * TIME',
                                         '\t80  SAVE moles',
                                         '\t-end', sep = "\n"),
                                       Enstatite = paste(
                                         '\t-start',
                                         '\t1   REM Ref PK04',
                                         '\t10  kacid = 10^(-9.02) * exp(-80.0e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^0.6',
                                         '\t20  kneut = 10^(-12.72) * exp(-80.0e3/8.314 * (1/TK-1/298.15))',
                                         '\t21  SSA = 0.1  # UNSOURCED surface area (m2/g): README section 6, decision D7',
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
                                         '\t21  SSA = 5  # UNSOURCED surface area (m2/g): README section 6, decision D7',
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
                                         '\t31  SSA = 5  # UNSOURCED surface area (m2/g): README section 6, decision D7',
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
                                         if (m == "Schreibersite") schreibersite_kinetics_entry(get_k_moles(m)) else sprintf("   %s\n      -m0 %.15e\n      -step_divide 1", m, get_k_moles(m))),
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
                                         sprintf("   %-20s 0 %.15e  dissolve_only", m, get_equil_moles(m))),
                                       collapse = "\n"
                                     )
                                     equil_prim_lines <- paste(
                                       c("# --- Primary minerals (equilibrium). Fe and Ni corrode completely in the first step: UNSOURCED instantaneous corrosion, README section 6, decision D6 ---",
                                         always_equil_prim_lines,
                                         if (length(equil_primary) > 0)
                                           sapply(equil_primary, function(m)
                                             sprintf("   %-20s 0 %.15e  dissolve_only", m, get_k_moles(m)))
                                         else character(0)),
                                       collapse = "\n")

                                     secondary_eq_lines <- paste(
                                       c("# --- Secondary minerals ---",
                                         sapply(secondary_minerals, function(m)
                                           sprintf("   %-20s 0 %.15e", m, get_equil_moles(m)))),
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

                                     .all_eq_so  <- c(setdiff(always_equil_primary, "Fe"), equil_primary, secondary_minerals)
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
# Gas phase (fixed pressure, Peng-Robinson): {paste(gas_comps, collapse = " ")}
{clathrate_note_lines}

KNOBS
   -logfile false

{.thick_line}
# Mode: {.mode_comment}
# prev_csv: {prev_csv}
{parent_fluid_header_lines(prev_csv)}
{parent_mode_line}
{parent_pressure_line}
# equil_minerals: {paste(equil_primary, collapse = ", ")}
# n_orbital_cycles: {n_orbital_cycles}
# orbits_per_step: {orbits_per_step}  (step duration = {sprintf("%.4e", orbits_per_step * orbital_period)} s)
# n_steps: {n_steps}
# Total simulation time: {sprintf("%.4e", total_time_years)} years  ({sprintf("%.6e", total_time_seconds)} s)

SOLUTION 1 {.soln_label}
    temp      {round(porewater_temp_C, 2)}
    pH        {sprintf("%.15g", ph_val)}
    pe        {sprintf("%.15g", pe_val)}
    redox     pe
    units     mol/kgw
    density   1
    water     {sprintf("%.15g", water_kg)}
    pressure  {round(.solution_pressure_atm, 1)}
{elem_lines}

{gas_phase_block(1, gas_comps, .solution_pressure_atm, porewater_temp_C, "Gas (fixed pressure)")}
{gas_restore_block(1, gas_comps, last_row)}

INCREMENTAL_REACTIONS true

{rates_block}

{kinetics_block}

EQUILIBRIUM_PHASES 1
{equil_phases_block}

SELECTED_OUTPUT
    -file {base_name}.csv
    -reset false
    -high_precision true
    -state true
    -solution true
    -step true
    -ph true
    -pe true
    -alkalinity true
    -ionic_strength true
    -water true
{totals_block(ch4_mode$database)}
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
{gas_so_line(gas_comps)}

# -high_precision (above) also sets the mole-balance convergence criterion to 1e-12; with a gas phase
# that made the first porewater step (the metal pulse) run for more than 15 minutes. This KNOBS, read
# after it, puts back the default 1e-8 (the last setting read wins, PHREEQC read.cpp), and keeps the
# 12-digit output.
KNOBS
    -convergence_tolerance 1e-8

USER_PUNCH
   -headings Time_Years Enst_remain_g Forst_remain_g Troil_remain_g Alb_remain_g Diop_remain_g Anorth_remain_g Kfs_remain_g Teph_remain_g Schr_remain_g Total_Initial_g{iom_punch_headings_restart}{SHARED_PUNCH_HEADINGS}{gas_punch(gas_comps, water_kg)$headings}
   10 PUNCH TOTAL_TIME / (365.25 * 24 * 3600)
{punch_lines}
   {10 + (length(all_primary) + 1) * 10} total_initial = {total_initial_expr}
   {10 + (length(all_primary) + 2) * 10} PUNCH total_initial
{iom_punch_lines_restart}
{SHARED_PUNCH_LINES}{gas_punch(gas_comps, water_kg)$lines}
END
')

                                     content <- self$stamp_pqi_header(self$apply_ch4_mode(content, ch4_mode), ch4_mode$database,
                                                                      list(pore = pore_pressure_atm, ocean = ocean_pressure_atm))
                                     writeLines(knobs_logfile(self$apply_ssa(content), logfile), filename)
                                     return(filename)
                                   },

                                   # ---------------------------------------------------------------
                                   # Tidal kinetic restart: continues a tidal-exchange kinetic run
                                   # from the end of a previous run.  Reads solution 1 (porewater),
                                   # solution 2 (ocean), kinetics, and both EP slots from the CSV.
                                   # Primary minerals listed in equil_minerals are moved from KINETICS
                                   # to EQUILIBRIUM_PHASES 1 (same as hydrosphere restart).
                                   # ---------------------------------------------------------------
                                   # Continues a tidal_kinetic run from its last complete orbit (REACT, MIX 4, MIX 1, MIX 2 rows): both
                                   # solutions (every element, 15 digits), both equilibrium-phase sets, kinetic reactants, IOM channels and both
                                   # gas phases. Checks the parent's carbon mode and pressures.
                                   #
                                   # @param output_dir Directory for the .pqi (and, when run there, its .csv/.out). Default "phreeqc_inputs".
                                   # @param prev_csv The parent run's CSV; its .pqi must be beside it (for the carbon-mode, pressure and fluid checks).
                                   # @param n_orbital_cycles Number of orbits (1.375 d each).
                                   # @param substeps KINETICS substeps per step. Default 1.
                                   # @param equil_minerals Primaries to move from KINETICS to EQUILIBRIUM_PHASES (dissolve_only). Default none; Schreibersite cannot be one.
                                   # @param mode "ocean_circulations" (default: auto grouping, n_circulations), "ungrouped" (1 orbit per step, n_orbital_cycles) or "grouped" (grouping_factor orbits per step; exchange compounds as 1 - (1 - x_tide)^G).
                                   # @param grouping_factor Orbits per step (mode "grouped").
                                   # @param exchange_ratio Target exchanged fluid / porewater per step; sets grouping_factor and mode "grouped". Default NULL.
                                   # @param porewater_temp Porewater reaction temperature: a number (deg C), "profile" (the profile's layer mean) or "ocean" (1.45 C, the profile's ocean floor). Default "ocean". Production: 150 and 160 C.
                                   # @param depletion_mineral Kinetic phase whose depletion stops the run early (USER_PRINT). Default: the fastest in depletion_order.
                                   # @param refill_kg Water (kg, on the run's normalised basis) moved from the ocean into the porewater
                                   #   before the first orbit, by MIX (negative: from the porewater to the ocean), so that the porewater fills
                                   #   the pore volume again (run_tidal_evolution(); PI 2026-10-30). Default 0. The MIX rows are not punched.
                                   # @param ocean_scale Factor on the parent's ocean (its water, the ocean's equilibrium phases and gas
                                   #   phase 2) on the run's normalised basis: when the reacting box changes between chunks (the skin depth
                                   #   moves with porosity), the real ocean is unchanged but the ocean per kg of the box's rock is not; the
                                   #   driver passes old box rock mass / new box rock mass (PI 2026-10-31). Default 1.
                                   # @param stop_on_depletion TRUE (default): stop the run (USER_PRINT STOP, which PHREEQC reports as 2 errors, status 3) at the
                                   #   end of an orbit once depletion_mineral is exhausted. FALSE runs every orbit (run_tidal_evolution() uses it).
                                   # @param suffix Appended to the file name. Default "".
                                   # @param CH4_redox_override Older spelling of carbon_mode: "coupled" or "decoupled" (= "ch4_uncoupled"). Default NULL; give this or carbon_mode, not both.
                                   # @param carbon_mode "c4_only" (NULL, the default: coupled carbon is C(+4) only, decided by the PI 2026-10-27), "coupled" (an upper bound on abiotic methanation), "ch4_uncoupled" (diagnostic), or "temperature" (the pre-10-26 rule). Picks the database.
                                   # @param ch4_cutoff_C Cutoff (deg C) for carbon_mode = "temperature" only. Default 150.
                                   # @param allow_pressure_change TRUE lets a restart run at pressures other than its parent's (with a warning). Default FALSE (an error).
                                   # @param pore_gas_fate "to_ocean" (default, PI 2026-10-29: the pore bubble moves to the ocean step each orbit; without plume loss an upper bound on ocean H2) or "retain" (it stays with the porewater; an upper bound on retention, ~2.2 pore volumes of gas). Carbon attribution is the same in both.
                                   # @return The path of the .pqi written.
                                   # @examples
                                   #   f <- ph$generate_phreeqc_input_tidal_kinetic_restart(prev_csv = "phreeqc_inputs/<parent>.csv", n_orbital_cycles = 2,
                                   #                                                         porewater_temp = 160, depletion_mineral = "Troilite")
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
                                     CH4_redox_override  = NULL,
                                     carbon_mode         = NULL,
                                     ch4_cutoff_C        = 150,
                                     allow_pressure_change = FALSE,
                                     pore_gas_fate       = c("to_ocean", "retain"),
                                     logfile             = FALSE,
                                     stop_on_depletion   = TRUE,
                                     refill_kg           = 0,
                                     ocean_scale         = 1
                                   ) {
                                     pore_gas_fate <- match.arg(pore_gas_fate)
                                     if (is.null(prev_csv) || !file.exists(prev_csv))
                                       stop("prev_csv must be a valid path to an existing CSV file")
                                     if (is.null(n_orbital_cycles) || n_orbital_cycles <= 0)
                                       stop("n_orbital_cycles must be a positive integer")
                                     # The parent's pore_gas_fate, from its header; a different one here is a warning.
                                     .ph <- sub("\\.csv$", ".pqi", prev_csv)
                                     .pf <- if (file.exists(.ph)) grep("^# (parent )?pore_gas_fate: ", readLines(.ph, n = 200, warn = FALSE), value = TRUE) else character(0)
                                     .pf <- if (length(.pf)) sub("^# (parent )?pore_gas_fate: ", "", .pf[1]) else "unknown (generated before 2026-10-28: gases were equilibrium phases)"
                                     if (.pf %in% c("retain", "to_ocean") && .pf != pore_gas_fate)
                                       warning(sprintf("restart: the parent ran pore_gas_fate = \"%s\", this restart runs \"%s\"", .pf, pore_gas_fate))
                                     parent_fate_line <- sprintf("# parent pore_gas_fate: %s", .pf)

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
                                     ocean_temp_C <- 274.6 - 273.15  # the profile's ocean floor, top of the rock layer (274.566 K), rounded
                                     # Porewater and ocean conditions from the profile; stops if profile_file is missing.
                                     pc <- profile_conditions(self$profile_file, self$simulator$reactive_thickness())
                                     pore_temp_C        <- pc$pore_temp_C
                                     pore_pressure_atm  <- pc$pore_pressure_atm
                                     ocean_pressure_atm <- pc$ocean_pressure_atm
                                     porewater_temp_C <- if (is.numeric(porewater_temp)) porewater_temp else
                                       if (porewater_temp == "profile") pore_temp_C else ocean_temp_C

                                     # CH4 (redox-coupled) vs Mtg (decoupled) — see resolve_ch4_mode().
                                     ch4_mode <- self$resolve_ch4_mode(porewater_temp_C, CH4_redox_override, carbon_mode, ch4_cutoff_C)
                                     parent_mode_line <- check_restart_carbon_mode(prev_csv, ch4_mode)
                                     parent_pressure_line <- check_restart_pressures(prev_csv, pore_pressure_atm, ocean_pressure_atm,
                                                                                     allow_pressure_change)
                                     ch4_override_note <- ch4_mode$note
                                     # Gas phases (see gas_components()): the porewater one, the ocean one, the
                                     # punch columns, and the empty porewater phase to_ocean starts each orbit with.
                                     # Ntg(g) when the parent carried it (a cometary-ice ocean)
                                     gas_comps <- gas_components(ch4_mode, if ("g_Ntg(g)" %in% names(row_mix4)) "Ntg(g)" else character(0))
                                     gp <- gas_punch(gas_comps, porewater_mass_norm, put = pore_gas_fate == "to_ocean")
                                     pore_gas_block <- gas_phase_block(1, gas_comps, pore_pressure_atm, porewater_temp_C,
                                                                       "Porewater gas (fixed pressure, starts empty)")
                                     ocean_gas_block <- gas_phase_block(2, gas_comps, ocean_pressure_atm, ocean_temp_C,
                                                                        "Ocean gas (fixed pressure, starts empty, never removed)")
                                     ch4_db_note_line <- sprintf(
                                       "# Recommended database: %s (CH4 %s; porewater T = %.1f C%s)",
                                       ch4_mode$database,
                                       ch4_mode$description,
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
                                     # ZERO_MOLES_FLOOR: precautionary, mirroring the fix in
                                     # bulk_kinetic_restart -- a dissolve_only EQUILIBRIUM_PHASES entry at a
                                     # literal floating-point-dust amount (e.g. 1e-21 mol) can make PHREEQC try
                                     # to satisfy that phase's saturation condition with essentially no mass
                                     # to work with, which failed to converge under a low water_rock_ratio
                                     # restart (bulk_kinetic_restart, SI far from 0 but m0 ~1e-21). Not yet
                                     # observed here, but the same dissolve_only pattern applies, so clamping
                                     # defensively; any real remaining mass is far above this threshold.
                                     ZERO_MOLES_FLOOR <- 1e-15
                                     get_k_moles <- function(mineral) {
                                       col <- paste0("k_", mineral)
                                       val <- if (col %in% names(row_react)) suppressWarnings(as.numeric(row_react[[col]])) else NA_real_
                                       if (!is.na(val) && val >= 0) return(if (val < ZERO_MOLES_FLOOR) 0 else val)
                                       # Fallback: use _remain_g column
                                       g_cols <- c(Enstatite="Enst_remain_g", Forsterite="Forst_remain_g",
                                                   Troilite="Troil_remain_g", Albite="Alb_remain_g",
                                                   Diopside="Diop_remain_g", Anorthite="Anorth_remain_g",
                                                   `K-Feldspar`="Kfs_remain_g", Tephroite="Teph_remain_g", Schreibersite="Schr_remain_g")
                                       g_col <- g_cols[mineral]
                                       if (!is.na(g_col) && g_col %in% names(row_react)) {
                                         g_val <- suppressWarnings(as.numeric(row_react[[g_col]]))
                                         if (!is.na(g_val) && g_val >= 0) {
                                           moles <- g_val / molar_masses[mineral]
                                           return(if (moles < ZERO_MOLES_FLOOR) 0 else moles)
                                         }
                                       }
                                       0
                                     }

                                     # Helper: read equilibrium phase moles from a given row (column via ep_column()).
                                     get_ep_moles <- function(mineral, row) {
                                       col <- ep_column(mineral, names(row))
                                       val <- if (!is.na(col)) suppressWarnings(as.numeric(row[[col]])) else NA_real_
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
                                     # Every element the parent punched (restart_elements()).
                                     elements <- restart_elements(ch4_mode$database, names(row_mix1), prev_csv)

                                     make_soln_lines <- function(row) {
                                       vals <- get_chem(row, elements)
                                       paste(Filter(Negate(is.null), lapply(seq_along(elements), function(i) {
                                         if (is.na(vals[i])) return(NULL)
                                         sprintf("    %-4s  %.15e", elements[i], vals[i])
                                       })), collapse = "\n")
                                     }

                                     soln1_elem_lines <- make_soln_lines(row_mix1)
                                     soln2_elem_lines <- make_soln_lines(row_mix2)
                                     soln1_pH   <- as.numeric(row_mix1[["pH"]])
                                     soln1_pe   <- as.numeric(row_mix1[["pe"]])
                                     soln1_h2o  <- as.numeric(row_mix1[["mass_H2O"]])
                                     soln2_pH   <- as.numeric(row_mix2[["pH"]])
                                     soln2_pe   <- as.numeric(row_mix2[["pe"]])
                                     if (!is.numeric(ocean_scale) || length(ocean_scale) != 1 || !(ocean_scale > 0)) stop("ocean_scale must be one positive number")
                                     soln2_h2o  <- as.numeric(row_mix2[["mass_H2O"]]) * ocean_scale
                                     # Refill (PI 2026-10-30): move refill_kg of water ocean -> porewater (negative: porewater ->
                                     # ocean) by MIX before the first orbit. MIX scales whole solutions, so the moved water carries
                                     # its solutes and every element is conserved between the two. Its simulations run before the
                                     # SELECTED_OUTPUT is defined, so they add no rows (and this restart has no i_soln rows).
                                     if (!is.numeric(refill_kg) || length(refill_kg) != 1 || !is.finite(refill_kg))
                                       stop("refill_kg must be one finite number (kg)")
                                     refill_block <- if (refill_kg > 0) {
                                       if (refill_kg >= soln2_h2o) stop(sprintf("refill_kg %.4g >= the ocean's water %.4g kg", refill_kg, soln2_h2o))
                                       paste0("END\n# Refill: ", sprintf("%.6e", refill_kg), " kg of ocean water into the porewater\n",
                                              "MIX 1\n  1  1.0\n", sprintf("  2  %.15e\n", refill_kg / soln2_h2o), "SAVE solution 1\nEND\n",
                                              "MIX 2\n", sprintf("  2  %.15e\n", 1 - refill_kg / soln2_h2o), "SAVE solution 2\nEND\n")
                                     } else if (refill_kg < 0) {
                                       if (-refill_kg >= soln1_h2o) stop(sprintf("-refill_kg %.4g >= the porewater's water %.4g kg", -refill_kg, soln1_h2o))
                                       paste0("END\n# Refill: ", sprintf("%.6e", -refill_kg), " kg of porewater into the ocean\n",
                                              "MIX 2\n  2  1.0\n", sprintf("  1  %.15e\n", -refill_kg / soln1_h2o), "SAVE solution 2\nEND\n",
                                              "MIX 1\n", sprintf("  1  %.15e\n", 1 + refill_kg / soln1_h2o), "SAVE solution 1\nEND\n")
                                     } else ""
                                     refill_header_line <- sprintf("# refill_kg: %.6e (ocean -> porewater before the first orbit; README section 6); ocean_scale: %.15g", refill_kg, ocean_scale)

                                     # ---- Secondary minerals list ----
                                     secondary_minerals <- c(
                                       "Analcime", "Anhydrite", "Aragonite", "Bassanite",
                                       "Beidellite-Ca", "Beidellite-Fe", "Beidellite-Mg", "Beidellite-Na",
                                       "Boehmite", "Brucite",
                                       "Calcite", 
                                       "Chalcedony", "Chamosite", "Chrysotile", "Citric_Acid",
                                       "Clinochlore-14A", "Clinochlore-7A", "Clinoptilolite-Ca", "Clinoptilolite-Na",
                                       "Clinozoisite", "Cronstedtite-7A",
                                       "Daphnite-14A", "Daphnite-7A", "Dawsonite", "Diaspore", "Dolomite",
                                       "Ettringite", "Fe(OH)2", "Fe(OH)3", "Fe2(SO4)3", "FeSO4",
                                       "Gibbsite", "Glycine", "Goethite", "Greenalite", "Gypsum", "Gyrolite",
                                       "Hematite", "Huntite", "Hydromagnesite",
                                       "Kaolinite", "Lizardite", "Magnetite", "Melanterite", "Mesolite", "MgSO4", "Minnesotaite", "Mirabilite",
                                       "Montmor-Ca", "Montmor-Mg", "Montmor-Na",
                                       "Natrolite", "Nontronite-Ca", "Nontronite-Mg", "Nontronite-Na",
                                       "Portlandite", "Pyrite",
                                       "Saponite-Fe-Ca", "Saponite-Fe-Fe", "Saponite-Fe-Mg", "Saponite-Fe-Na",
                                       "Saponite-Mg-Ca", "Saponite-Mg-Fe", "Saponite-Mg-Mg", "Saponite-Mg-Na",
                                       "Sepiolite", "SiO2(am)", "Thenardite", "Tobermorite-11A"
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
                                         '\t21  SSA = 0.1  # UNSOURCED surface area (m2/g): README section 6, decision D7',
                                         '\t22  mw = 140.6715',
                                         '\t40  k = (kacid + kneut) * SSA * mw * M',
                                         '\t50  IF SR("Forsterite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Forsterite"))',
                                         '\t60  moles = rate * TIME',
                                         '\t70  SAVE moles',
                                         '\t-end', sep = "\n"),
                                       Troilite = paste(
                                         '\t-start',
                                         '\t1   REM PK04 has no separate troilite entry, hexagonal pyrrhotite parameters used as proxy',
                                         '\t25  kacid = 10^(-6.79) * exp(-63.0e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^-0.090 * ACT("Fe+3")^0.356',
                                         '\t26  SSA = 5  # UNSOURCED surface area (m2/g): README section 6, decision D7',
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
                                         '\t21  SSA = 0.1  # UNSOURCED surface area (m2/g): README section 6, decision D7',
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
                                         '\t31  SSA = 5  # UNSOURCED surface area (m2/g): README section 6, decision D7',
                                         '\t32  mw = 278.33',
                                         '\t40  k = (kacid + kneut + kbase) * SSA * mw * M',
                                         '\t50  IF SR("K-Feldspar") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("K-Feldspar"))',
                                         '\t60  moles = rate * TIME',
                                         '\t70  SAVE moles',
                                         '\t-end', sep = "\n"),
                                       Tephroite = paste(
                                         '\t-start',
                                         '\t1   REM Pincus et al. (2026) ACS Earth Space Chem 10, 1174-1184, Table S1,',
                                         '\t2   REM citing Casey et al. (1993) GCA 57, 785-793. See generate_phreeqc_input_tidal_kinetic',
                                         '\t3   REM for full derivation notes -- identical rate law, ported here verbatim.',
                                         '\t34  REM pH >= 7: acid branch evaluated at a(H+) = 1e-7 (flat at the pH-7 rate). No neutral term is summed.',
                                         '\t35  aH = ACT("H+")',
                                         '\t36  IF aH < 1e-7 THEN aH = 1e-7',
                                         '\t37  Ea = -60.11*LOG(-LOG10(aH)) + 106.17',
                                         '\t38  IF Ea < 0 THEN Ea = 0',
                                         '\t40  kacid = 10^(-4.46) * exp(-Ea*1e3/8.314 * (1/TK-1/298.15)) * aH^0.47',
                                         '\t46  SSA = 0.1  # UNSOURCED surface area (m2/g): README section 6, decision D7',
                                         '\t47  mw = 201.96',
                                         '\t50  k = kacid * SSA * mw * M',
                                         '\t60  IF SR("Tephroite") > 1 THEN rate = 0 ELSE rate = k * (1 - SR("Tephroite"))',
                                         '\t70  moles = rate * TIME',
                                         '\t80  SAVE moles',
                                         '\t-end', sep = "\n"),
                                       Enstatite = paste(
                                         '\t-start',
                                         '\t1   REM Ref PK04',
                                         '\t10  kacid = 10^(-9.02) * exp(-80.0e3/8.314 * (1/TK-1/298.15)) * ACT("H+")^0.6',
                                         '\t20  kneut = 10^(-12.72) * exp(-80.0e3/8.314 * (1/TK-1/298.15))',
                                         '\t21  SSA = 0.1  # UNSOURCED surface area (m2/g): README section 6, decision D7',
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
                                         '\t21  SSA = 5  # UNSOURCED surface area (m2/g): README section 6, decision D7',
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
                                         '\t31  SSA = 5  # UNSOURCED surface area (m2/g): README section 6, decision D7',
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
                                       self$resolve_iom_pools(ch4_mode)

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
                                         sprintf("   %s\n      -formula %s\n      -m0 %.15e\n      -step_divide 1\n",
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
                                           if (m == "Schreibersite") schreibersite_kinetics_entry(get_k_moles(m)) else sprintf("   %s\n      -m0 %.15e\n      -step_divide 1", m, get_k_moles(m))),
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
                                       c("# --- Primary minerals (equilibrium). Fe and Ni corrode completely in the first step: UNSOURCED instantaneous corrosion, README section 6, decision D6 ---",
                                         sapply(always_equil_primary, function(m)
                                           sprintf("   %-20s 0 %.15e  dissolve_only", m, get_ep_moles(m, row_react))),
                                         if (length(equil_primary) > 0)
                                           sapply(equil_primary, function(m)
                                             sprintf("   %-20s 0 %.15e  dissolve_only", m, get_k_moles(m)))
                                         else character(0)),
                                       collapse = "\n")

                                     ep1_sec_lines <- paste(
                                       c("# --- Secondary minerals (porewater, from previous run) ---",
                                         sapply(c(secondary_minerals, extra_pw_minerals), function(m)
                                           sprintf("   %-20s 0 %.15e", m, get_ep_moles(m, row_react)))),
                                       collapse = "\n"
                                     )

                                     ep1_block <- paste0("EQUILIBRIUM_PHASES 1\n",
                                       paste(Filter(nzchar, c(ep1_prim_lines, ep1_sec_lines)), collapse = "\n"))

                                     # ---- EQUILIBRIUM_PHASES 2 (ocean: EP2 state, equil primaries + secondaries) ----
                                     ep2_prim_lines <- if (length(equil_primary) > 0)
                                       paste(c("# --- Primary minerals (equilibrium, ocean) ---",
                                               sapply(equil_primary, function(m)
                                                 sprintf("   %-20s 0 %.15e  dissolve_only", m, ocean_scale * get_ep_moles(m, row_mix4)))),
                                             collapse = "\n")
                                     else ""
                                     ep2_minerals <- c(secondary_minerals, extra_ocean_minerals)
                                     ep2_sec_lines <- paste(
                                       c("# --- Secondary minerals (ocean, from previous run) ---",
                                         sapply(ep2_minerals, function(m)
                                           sprintf("   %-20s 0 %.15e", m, ocean_scale * get_ep_moles(m, row_mix4)))),
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

                                     # Magnetite is already in secondary_minerals (2026-10-07); appending it
                                     # again here listed it twice in -equilibrium_phases.
                                     .all_eq_so  <- c(setdiff(always_equil_primary, "Fe"), equil_primary, secondary_minerals,
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
                                     # Default order: primary minerals fastest-to-slowest by pseudo-first-
                                     # order rate constant (kacid+kneut[+kbase])*SSA*mw at 150 C, pH 7,
                                     # SR=0 (see CHANGELOG 2026-10-01 rate ranking), with Troilite forced
                                     # last regardless of rank — its rate depends on ACT(Fe3+), which swings
                                     # it from fastest to slowest primary across a plausible activity range,
                                     # so it's deliberately not trusted as an early depletion/proxy pick.
                                     # Tephroite deliberately OMITTED for now: its rate law is still a
                                     # kacid=kneut=0 placeholder (never depletes), pending a real rate —
                                     # insert it once that's supplied. Omitting it from this list doesn't
                                     # remove it as a kinetic phase, it just can't be auto-selected as
                                     # depletion_mineral/proxy_mineral until a real rate makes that meaningful;
                                     # explicit depletion_mineral="Tephroite" still works regardless.
                                     # IOM channels follow, fastest-to-slowest by k at 150 C (same ranking).
                                     # Pass depletion_mineral explicitly to override any of this.
                                     depletion_order <- c("Albite", "Diopside", "Anorthite", "Forsterite",
                                                          "Schreibersite", "K-Feldspar", "Enstatite", "Troilite",
                                                          "IOM_CO2_44", "IOM_CO2_46", "IOM_CO2_48", "IOM_CO2_50",
                                                          "IOM_CHn_50", "IOM_CO2_52", "IOM_CHn_52", "IOM_CH4_52",
                                                          "IOM_CO2_54", "IOM_CHn_54", "IOM_CH4_54", "IOM_N", "IOM_S",
                                                          "IOM_CO2_56", "IOM_CHn_56", "IOM_CH4_56",
                                                          "IOM_CHn_58", "IOM_CH4_58", "IOM_CHn_60", "IOM_CH4_60",
                                                          "IOM_CHn_62", "IOM_CH4_62", "IOM_CH4_64", "IOM_CH4_66",
                                                          "IOM_CH4_68", "IOM_CH4_70", "IOM_CH4_72", "IOM_CH4_74",
                                                          "IOM_CH4_76")
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
# pore_gas_fate: {pore_gas_fate}
# Gas phases (fixed pressure, Peng-Robinson): {paste(gas_comps, collapse = " ")}
{parent_fate_line}
{clathrate_note_lines}

KNOBS
   -logfile false

# k={self$simulator$params$k}, thickness={self$simulator$params$layer_thickness} m
# Mode: tidal kinetic restart from previous CSV state
{refill_header_line}
{self$simulator$exchange_header_line()}
# prev_csv: {prev_csv}
{parent_fluid_header_lines(prev_csv)}
{parent_mode_line}
{parent_pressure_line}
# restart_time_offset: {sprintf("%.6e", restart_time_years)} years
# equil_minerals: {paste(equil_primary, collapse = ", ")}
# n_orbital_cycles: {n_orbital_cycles}
# substeps: {substeps}
# mode: {mode}
# time_step: {sprintf("%.4e", time_step_years)} years

SOLUTION 1 Enceladus porewater (restarted)
    temp      {round(porewater_temp_C, 2)}
    pH        {sprintf("%.15g", soln1_pH)}
    pe        {sprintf("%.15g", soln1_pe)}
    redox     pe
    units     mol/kgw
    density   1
    water     {sprintf("%.15g", soln1_h2o)}
    pressure  {round(pore_pressure_atm, 1)}
{soln1_elem_lines}

SOLUTION 2 Enceladus ocean (restarted)
    temp      {round(ocean_temp_C, 2)}
    pH        {sprintf("%.15g", soln2_pH)}
    pe        {sprintf("%.15g", soln2_pe)}
    redox     pe
    units     mol/kgw
    density   1
    water     {sprintf("%.15g", soln2_h2o)}
    pressure  {round(ocean_pressure_atm, 1)}
{soln2_elem_lines}
{refill_block}

REACTION_TEMPERATURE 1
    {round(porewater_temp_C, 2)}

REACTION_TEMPERATURE 2
    {round(ocean_temp_C, 2)}

# Reaction pressures, used by every step like the reaction temperatures. MIX averages the
# solution pressures by mass, so without these the porewater would drift towards the ocean
# pressure after each exchange.
REACTION_PRESSURE 1
    {round(pore_pressure_atm, 1)}

REACTION_PRESSURE 2
    {round(ocean_pressure_atm, 1)}

{pore_gas_block}
{if (pore_gas_fate == "retain") gas_restore_block(1, gas_comps, row_react) else "# pore_gas_fate = to_ocean: the porewater gas phase starts empty"}

{ocean_gas_block}
{gas_restore_block(2, gas_comps, row_mix4, ocean_scale)}
{if (pore_gas_fate == "to_ocean") paste0("\n", gas_transfer_block(gas_comps)) else ""}

{rates_block}

{kinetics_block}

{ep1_block}

{ep2_block}

SELECTED_OUTPUT
    -file {base_name}.csv
    -reset false
    -high_precision true
    -reaction true
    -state true
    -solution true
    -step true
    -ph true
    -pe true
    -alkalinity true
    -ionic_strength true
    -water true
{totals_block(ch4_mode$database)}
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
{gas_so_line(gas_comps)}

# -high_precision (above) also sets the mole-balance convergence criterion to 1e-12; with a gas phase
# that made the first porewater step (the metal pulse) run for more than 15 minutes. This KNOBS, read
# after it, puts back the default 1e-8 (the last setting read wins, PHREEQC read.cpp), and keeps the
# 12-digit output.
KNOBS
    -convergence_tolerance 1e-8

# USER_PUNCH is defined per-cycle below (with exact cumulative time embedded)
    ')

                                     # ---- Cycle loop ----
                                     cycle_blocks <- paste(sapply(seq_len(kinetic_cycles), function(i) {
                                       cumulative_time <- restart_time_years + i * time_step_years
                                       paste0(
                                         sprintf("\n# --- Cycle %d of %d (cumulative time: %.4e years) ---\n",
                                                 i, kinetic_cycles, cumulative_time),
                                         # to_ocean: the porewater gas phase is redefined empty at the start of each orbit
                                         if (pore_gas_fate == "to_ocean" && i > 1) paste0(pore_gas_block, "\n") else "",
                                         "USE solution 1\n",
                                         "USE kinetics 1\n",
                                         "USE equilibrium_phases 1\n",
                                         "USE gas_phase 1\n",
                                         "USE reaction none\n",
                                         "USE reaction_temperature 1\n",
                                         "USE reaction_pressure 1\n",
                                         "SAVE solution 1\n",
                                         "SAVE equilibrium_phases 1\n",
                                         if (pore_gas_fate == "retain") "SAVE gas_phase 1\n" else "",
                                         "USER_PUNCH\n",
                                         "   -headings Time_Years Enst_remain_g Forst_remain_g Troil_remain_g Alb_remain_g Diop_remain_g Anorth_remain_g Kfs_remain_g Teph_remain_g Schr_remain_g Total_Initial_g",
                                         iom_punch_headings_tidal, SHARED_PUNCH_HEADINGS, gp$headings, "\n",
                                         sprintf("   10 PUNCH %.6f  # Cumulative time (years)\n", cumulative_time),
                                         paste(sapply(seq_along(all_primary), function(j) {
                                           make_punch_line(10 + j * 10, all_primary[j])
                                         }), collapse = "\n"), "\n",
                                         sprintf("   %d total_initial = %s\n", 10 + (length(all_primary) + 1) * 10, total_initial_expr),
                                         sprintf("   %d PUNCH total_initial\n", 10 + (length(all_primary) + 2) * 10),
                                         iom_punch_lines_tidal,
                                         SHARED_PUNCH_LINES, gp$lines,
                                         "END\n\n",
                                         OCEAN_STEP_KNOBS_ON,
                                         if (pore_gas_fate == "to_ocean") "USE kinetics 3\n" else "USE kinetics none\n",  # 3: the pore bubble
                                         "USE equilibrium_phases 2\n",
                                         "USE gas_phase 2\n",
                                         "USE reaction none\n",
                                         "USE reaction_temperature 2\n",
                                         "USE reaction_pressure 2\n",
                                         "MIX 4\n",
                                         sprintf("  1  %.15e\n", x_step),
                                         "  2  1.0\n",
                                         "SAVE solution 4\n",
                                         "SAVE equilibrium_phases 2\n",
                                         "SAVE gas_phase 2\n",
                                         "END\n\n",
                                         OCEAN_STEP_KNOBS_OFF,
                                         "USE kinetics none\n",
                                         "USE equilibrium_phases none\n",
                                         "USE reaction none\n",
                                         "USE reaction_temperature 1\n",
                                         "USE reaction_pressure 1\n",
                                         "MIX 1\n",
                                         sprintf("  1  %.15e\n", 1 - x_step),
                                         sprintf("  4  %.15e\n", f_x / (f_ocean + f_x)),
                                         "SAVE solution 1\n",
                                         "END\n\n",
                                         "USE kinetics none\n",
                                         "USE equilibrium_phases none\n",
                                         "USE reaction none\n",
                                         "MIX 2\n",
                                         sprintf("  4  %.15e\n", f_ocean / (f_ocean + f_x)),
                                         "SAVE solution 2\n",
                                         "END\n",
                                         if (stop_on_depletion && i < kinetic_cycles &&
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

                                     content <- self$stamp_pqi_header(self$apply_ch4_mode(paste0(template, "\n", cycle_blocks), ch4_mode),
                                                                      ch4_mode$database, list(pore = pore_pressure_atm, ocean = ocean_pressure_atm))
                                     writeLines(knobs_logfile(self$apply_ssa(content), logfile), filename)
                                     return(filename)
                                   }
                                 )
)

