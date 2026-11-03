# The reference case (PI 2026-11-01/02). Every other case inherits from this file and overrides a few fields.
# Run: Rscript run_case.R cases/reference.R  ->  runs/reference/ (PRODUCTION_GUIDE.md).
# Each field: units; allowed range; the section of docs/MODEL.md (or README.md) that justifies it. Do not add fields: an unknown field is
# an error (production_kit.R, validate_case()).
case <- list(
  name        = "reference",   # the run directory runs/<name>; letters, digits, _ . + -
  description = "Reference: Kozeny-Carman k (d = 1.59 cm), phi0 0.32, eta 1e12 Pa s dissipative core, 15 wt% IOM, pure-water ocean",
  label       = "",            # free text carried to RUNINFO.txt and summary.csv (e.g. "not at steady state by design")

  # --- geometry and temperature ---
  layer_thickness_m = 194200,    # m; (0, 194200]; the permeable layer, 194.2 km = the whole core. docs/MODEL.md (The reacting box)
  reactive_depth    = "skin",    # "skin" | "layer"; the box that reacts: the top min(layer, skin depth). docs/MODEL.md (The reacting box is the skin)
  porewater_temp    = "profile", # "profile" (the skin's own profile mean) | "ocean" | a number 0-300 degC (a bound, e.g. 160). docs/MODEL.md

  # --- rock and permeability ---
  phi0     = 0.32,               # -; (0, phi_c); initial porosity (PlanetProfile phiRockMax). docs/MODEL.md (Porosity and the Biot coefficient)
  phi_c    = 0.4,                # -; (phi0, 1]; critical porosity: alpha = phi/phi_c, K_d = K_grain (1 - phi/phi_c). docs/MODEL.md
  phi_min  = 0.01,               # -; (0, phi0); the porosity at which the rock counts as clogged. docs/MODEL.md (Porosity evolution)
  grain_diameter_m = 0.0159,     # m; 1e-7 to 0.1; Kozeny-Carman k = d^2 phi^3 / (180 (1 - phi)^2), following phi(t). docs/MODEL.md (Permeability and surface area from grain size)
  k_m2     = NULL,               # m2; 1e-20 to 1e-6; a fixed k instead of grain_diameter_m (give one of the two)
  k_floor_m2 = NULL,             # m2; 1e-20 to 1e-8, or NULL (none: the rock clogs and the run ends); fracture-maintained k (Kelemen & Hirth 2012). docs/MODEL.md (Clogging and k_floor)

  # --- tides ---
  core_viscosity_Pa_s = 1e12,    # Pa s; within the scenario's table: dissipative_core 4.4e11-2.3e12, dissipative_shell 3e13-1e22. docs/MODEL.md (The strain)
  tidal_scenario = "dissipative_core",  # "dissipative_core" | "dissipative_shell" (Bagheri et al. 2025). docs/MODEL.md (The strain)

  # --- chemistry ---
  iom_wt_percent = 15,           # wt%; 0-30; IOM in the rock (organic run; the blank has the same minerals, no IOM). docs/MODEL.md (IOM composition)
  carbon_mode    = "c4_only",    # "c4_only" | "coupled" (upper bound on methanation) | "ch4_uncoupled" | "temperature". docs/MODEL.md (Carbon mode)
  initial_ocean  = "pure_water", # "pure_water" | "cometary" (cometary_ice_ocean()). docs/MODEL.md (Initial ocean; The cometary-ice ocean)
  cometary_f     = 1,            # -; (0, 1]; fraction of the cometary volatiles retained (initial_ocean = "cometary" only). docs/MODEL.md
  pore_gas_fate  = "to_ocean",   # "to_ocean" | "retain" | "retain_to_pore_volume"; the pore bubble's fate. docs/MODEL.md (Gas)
  trace_C_molal  = 1e-14,        # mol/kgw; 0-1e-6; numerical carbon seed. docs/MODEL.md (Trace C seed)
  phosphate_phases = "apatite",  # "apatite" (Hydroxyapatite + precursors) | "precursor" (apatite kinetically inhibited); no "none". docs/MODEL.md (Phosphate phases)

  # --- time and end conditions (the first reached ends the run) ---
  t_end_yr = 1e8,                # yr; (0, 4.5e9]; the backstop (PI 2026-11-02). README 3 (Production runs)
  ss_tol   = 0.01,               # -; (0, 0.5], or NULL (no steady-state stop); ocean C(+4), pH, Na within ss_tol over one circulation time. README 3

  # --- accuracy (set by validation, PI 2026-11-02: the largest values within 1% of much finer steps; do not
  #     loosen: docs/MODEL.md, Step limits) ---
  phi_tol    = 0.003,            # -; largest relative porosity change per chunk (0.01 moved ocean C(+4) by 2%)
  x_step_max = 0.03,             # -; largest fraction of the pore water exchanged per PHREEQC step (0.1 moved Na by 2%)
  steps_per_chunk = 4,           # -; PHREEQC steps per chunk

  generator = "tidal_kinetic",   # the only generator for long runs
  blank     = TRUE               # TRUE: also run the same-rock blank (no IOM), for organic minus blank
)
