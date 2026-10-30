# Porosity that evolves with the chemistry (exchange-physics E, PI 2026-10-29). Sourced by
# ENTICES_LF_20260917.R after phreeqc_tools.R.
#
# Rigid frame, no compaction: phi(t) = 1 - V_solid(t) / V_total, with V_total = V_solid(0) / (1 - phi0) fixed
# and V_solid from bulk_density() (minerals, the remaining IOM at the Delarue density, clathrates). The run
# is split into chunks of tidal_kinetic (first) and tidal_kinetic_restart (the rest); between chunks, on the
# R side, phi sets a new Simulator: k (if k_of_phi is given, e.g. Kozeny-Carman), alpha (phi/phi_c), the
# storage, the skin depth, the reacting box, the strain-limited volume and the exchange fraction.
#
# Chunk length: a chunk whose porosity changes by more than phi_tol (relative) is rerun with half the orbits
# (down to 1 orbit, which is accepted and flagged); the next chunk's length is scaled from the last one's
# rate, n * phi_tol / |dphi/phi|, between 1 and max_chunk orbits.
#
# Clogging: if V_solid reaches V_total (phi <= phi_min), phi is capped at phi_min and `clogged` is set.
#
# Water accounting: PHREEQC's porewater mass (it loses water to hydration and gains it by exchange) is
# compared with the pore water phi(t) implies, phi * V_total * 1 g/cm3. PHREEQC's mass is authoritative for
# the chemistry (it is what conserves H and O); phi from the solids is authoritative for the transport
# parameters. Both are reported (water_phreeqc_kg, water_from_phi_kg) on the run's normalised basis.
#
# porosity_evolution = FALSE (default until tested; then the PI decides) runs the same chunks but keeps
# phi0 throughout, reporting the phi the solids imply.

# One Simulator for porosity phi and permeability k, with everything else from the template simulator.
evolution_simulator <- function(template, phi, k) {
  s <- Simulator$new(k = k, layer_thickness = template$params$layer_thickness, config = template$config,
                     reactive_depth = template$params$reactive_depth, phi0 = phi,
                     phi_c = template$params$phi_c, alpha = template$params$alpha)
  invisible(capture.output(s$run_simulation()))
  s
}

# template: a Simulator (k, layer, phi0, phi_c, config, reactive_depth). organic_wt_percent: as for
# calculate_phreeqc_parameters(). n_orbits: total orbits. Further arguments go to both generators
# (porewater_temp, pore_gas_fate, carbon_mode, ...); trace_C_molal only to the first.
# Returns list(chunks = one row per accepted chunk, files = the chunks' .pqi).
run_tidal_evolution <- function(template, organic_wt_percent, n_orbits, output_dir, porosity_evolution = FALSE,
                                phi_tol = 0.01, phi_min = 0.01, k_of_phi = NULL, first_chunk = 1, max_chunk = 64,
                                trace_C_molal = 1e-14, ...) {
  quiet <- function(e) { invisible(capture.output(r <- suppressMessages(e))); r }
  phi0 <- template$params$porosity; phi <- phi0
  k <- if (is.null(k_of_phi)) template$params$k else k_of_phi(phi)
  done <- 0; n <- first_chunk; prev_csv <- NULL; origin <- NULL; V_total <- NA; clogged <- FALSE
  rows <- list(); files <- character(0); chunk <- 0
  while (done < n_orbits) {
    n <- min(n, n_orbits - done)
    sim <- evolution_simulator(template, phi, k)
    ph <- quiet(PhreeqcIntegrator$new(sim)); quiet(ph$calculate_phreeqc_parameters(organic_wt_percent = organic_wt_percent))
    dir <- file.path(output_dir, sprintf("chunk%03d_orb%g-%g", chunk + 1, done + 1, done + n))
    unlink(dir, recursive = TRUE)
    f <- suppressWarnings(quiet(if (is.null(prev_csv))
      ph$generate_phreeqc_input_tidal_kinetic(output_dir = dir, mode = "ungrouped", n_orbital_cycles = n, trace_C_molal = trace_C_molal, ...)
    else
      ph$generate_phreeqc_input_tidal_kinetic_restart(output_dir = dir, prev_csv = prev_csv, n_orbital_cycles = n,
                                                      allow_pressure_change = TRUE, stop_on_depletion = FALSE, ...)))
    db <- sub("^# Database: ", "", grep("^# Database: ", readLines(f, n = 30), value = TRUE)[1])
    r <- run_phreeqc(f, db, output = sub("\\.pqi$", ".out", f), workdir = dirname(f), warn_version = FALSE)
    if (r$status != 0 || r$errors > 0) stop(sprintf("chunk %d (%s): PHREEQC status %s, %d errors", chunk + 1, basename(f), r$status, r$errors))
    csv <- sub("\\.pqi$", ".csv", f)
    bd <- suppressWarnings(bulk_density(csv, porosity = phi, origin_pqi = if (is.null(origin)) f else origin))$summary
    if (is.null(origin)) V_total <- bd$solid_volume_cm3[1] / (1 - phi0)
    V_end <- bd$solid_volume_cm3[nrow(bd)]
    phi_raw <- 1 - V_end / V_total
    rel <- abs(phi_raw - phi) / phi
    if (porosity_evolution && rel > phi_tol && n > 1) { n <- max(1, floor(n / 2)); next }   # too big a step: redo
    chunk <- chunk + 1; files <- c(files, f)
    if (is.null(origin)) origin <- f
    x <- read.delim(csv, check.names = FALSE, strip.white = TRUE); names(x) <- trimws(names(x))
    rr <- which(trimws(x$state) == "react" & trimws(x$soln) == "1")
    w_ph <- as.numeric(x$mass_H2O[rr[length(rr)]])          # the porewater after the last exchange (MIX 1)
    xch <- sim$exchange
    rows[[chunk]] <- data.frame(chunk = chunk, orbit_from = done + 1, orbit_to = done + n, n_orbits = n,
      phi_start = phi, phi_end_solids = phi_raw, dphi_rel = rel, over_tol = rel > phi_tol, k = k, alpha = sim$biot_alpha(),
      delta_m = xch$delta, L_eff_m = sim$reactive_thickness(), x_half = xch$x_half, binding = xch$binding,
      circulation_yr = quiet(sim$get_summary_metrics())$circulation_time_years, clogged = clogged,
      V_solid_cm3 = V_end, V_total_cm3 = V_total, water_phreeqc_kg = w_ph, water_from_phi_kg = max(phi_raw, 0) * V_total / 1000,
      seconds = r$seconds)
    done <- done + n; prev_csv <- csv
    if (porosity_evolution) {
      clogged <- phi_raw <= phi_min
      phi <- max(phi_min, phi_raw)
      if (!is.null(k_of_phi)) k <- k_of_phi(phi)
    }
    n <- if (rel > 0) max(1, min(max_chunk, floor(n * phi_tol / rel))) else min(max_chunk, 2 * n)
  }
  list(chunks = do.call(rbind, rows), files = files)
}
