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
# Water accounting and refill (PI 2026-10-30): hydration consumes porewater. With refill = TRUE (the default
# when porosity evolves), between chunks the porewater is brought back to the pore volume phi(t) implies,
# target = m_pw(0) * phi(t) / phi0 (m_pw(0): the first chunk's initial porewater, so the target is on the
# run's own normalised basis and the refill is 0 at t = 0), by moving ocean water into it
# (tidal_kinetic_restart(refill_kg =), a MIX that conserves every element; negative if the pores shrink
# below the water they hold). Reported per chunk: refill_in_kg (moved before the chunk), refill_kg_per_orbit
# (that over the previous chunk's orbits, i.e. the hydration-driven flux), exchange_kg_per_orbit (the tidal
# exchange, x_orbit * porewater), and the porewater and ocean water at the chunk's end. PHREEQC's masses are
# authoritative for the chemistry; phi from the solids for transport (water_from_phi_kg is phi V_total for
# reference).
#
# The ocean per kg of the box's rock (PI 2026-10-31): the box's rock mass changes when the skin depth moves,
# while the real ocean does not, so each restart scales the parent's ocean (water, ocean minerals, ocean gas)
# by old box rock mass / new box rock mass (tidal_kinetic_restart(ocean_scale =)); the ocean per kg of rock
# then stays the profile's ocean / the box's rock (PhreeqcIntegrator uses the profile's Mocean_kg).
#
# Drift (PI 2026-10-31): each chunk's whole-system element totals (porewater, ocean, rock, IOM, pore and ocean
# gas) are compared with its start (t = 0 for the first chunk; the parent's end with the ocean scaled, for a
# restart; the refill only moves water between the two); the differences accumulate over chunks, and
# drift_max is the largest element's cumulative difference relative to its current whole-system total
# (drift_element names it). A warning is given above 1e-4. The production tolerance is unchanged.
#
# porosity_evolution: TRUE for production runs (default, PI 10-30, once the refill passed its checks);
# FALSE keeps phi0 (and no refill) and reports the phi the solids imply; the quick self-tests use both.

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
run_tidal_evolution <- function(template, organic_wt_percent, n_orbits, output_dir, porosity_evolution = TRUE,
                                phi_tol = 0.01, phi_min = 0.01, k_of_phi = NULL, first_chunk = 1, max_chunk = 64,
                                trace_C_molal = 1e-14, refill = porosity_evolution, ...) {
  quiet <- function(e) { invisible(capture.output(r <- suppressMessages(e))); r }
  phi0 <- template$params$porosity; phi <- phi0
  k <- if (is.null(k_of_phi)) template$params$k else k_of_phi(phi)
  done <- 0; n <- first_chunk; prev_csv <- NULL; origin <- NULL; V_total <- NA; clogged <- FALSE
  m_pw0 <- NA; refill_next <- 0; n_prev <- NA; lm_prev <- NA; pq0 <- NULL; parts_prev <- NULL; drift <- NULL
  rows <- list(); files <- character(0); chunk <- 0
  while (done < n_orbits) {
    n <- min(n, n_orbits - done)
    sim <- evolution_simulator(template, phi, k)
    ph <- quiet(PhreeqcIntegrator$new(sim)); quiet(ph$calculate_phreeqc_parameters(organic_wt_percent = organic_wt_percent))
    ocean_scale <- if (is.na(lm_prev)) 1 else lm_prev / ph$phreeqc_params$layer_mass
    dir <- file.path(output_dir, sprintf("chunk%03d_orb%g-%g", chunk + 1, done + 1, done + n))
    unlink(dir, recursive = TRUE)
    f <- suppressWarnings(quiet(if (is.null(prev_csv))
      ph$generate_phreeqc_input_tidal_kinetic(output_dir = dir, mode = "ungrouped", n_orbital_cycles = n, trace_C_molal = trace_C_molal, ...)
    else
      ph$generate_phreeqc_input_tidal_kinetic_restart(output_dir = dir, prev_csv = prev_csv, n_orbital_cycles = n,
                                                      allow_pressure_change = TRUE, stop_on_depletion = FALSE,
                                                      refill_kg = refill_next, ocean_scale = ocean_scale, ...)))
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
    ro <- which(trimws(x$state) == "react" & trimws(x$soln) == "2")
    w_oc <- as.numeric(x$mass_H2O[ro[length(ro)]])          # the ocean (MIX 2)
    if (is.na(m_pw0)) m_pw0 <- as.numeric(x$mass_H2O[which(trimws(x$state) == "i_soln" & trimws(x$soln) == "1")[1]])
    # drift: this chunk's end against its start
    if (is.null(pq0)) pq0 <- pqi_initial_solids(f)
    fate <- sub("^# pore_gas_fate: ", "", grep("^# pore_gas_fate: ", readLines(f, n = 80), value = TRUE)[1])
    parts <- whole_system_parts(x, pq0, fate)
    start <- if (is.null(parts_prev)) whole_system_t0(x, pq0) else parts_prev$pore + ocean_scale * parts_prev$ocean
    end <- parts$pore + parts$ocean
    drift <- (if (is.null(drift)) 0 else drift) + (end - start)
    rel_d <- ifelse(abs(end) > 0, abs(drift) / abs(end), 0)
    if (max(rel_d) > 1e-4) warning(sprintf("run_tidal_evolution: chunk %d: cumulative drift of %s is %.2e relative", chunk, names(rel_d)[which.max(rel_d)], max(rel_d)))
    parts_prev <- parts
    xch <- sim$exchange
    rows[[chunk]] <- data.frame(chunk = chunk, orbit_from = done + 1, orbit_to = done + n, n_orbits = n,
      phi_start = phi, phi_end_solids = phi_raw, dphi_rel = rel, over_tol = rel > phi_tol, k = k, alpha = sim$biot_alpha(),
      delta_m = xch$delta, L_eff_m = sim$reactive_thickness(), x_half = xch$x_half, binding = xch$binding,
      circulation_yr = quiet(sim$get_summary_metrics())$circulation_time_years, clogged = clogged,
      V_solid_cm3 = V_end, V_total_cm3 = V_total, water_phreeqc_kg = w_ph, water_ocean_kg = w_oc,
      water_from_phi_kg = max(phi_raw, 0) * V_total / 1000, water_target_kg = m_pw0 * max(phi_min, phi_raw) / phi0,
      refill_in_kg = refill_next, refill_kg_per_orbit = if (is.na(n_prev)) 0 else refill_next / n_prev,
      exchange_kg_per_orbit = xch$x_orbit * ph$phreeqc_params$porewater_mass / ph$phreeqc_params$layer_mass,
      layer_mass_kg = ph$phreeqc_params$layer_mass, ocean_scale = ocean_scale,
      ocean_per_rock = w_oc, ocean_per_rock_geometric = ph$phreeqc_params$ocean_mass / ph$phreeqc_params$layer_mass,
      drift_max = max(rel_d), drift_element = names(rel_d)[which.max(rel_d)],
      seconds = r$seconds)
    done <- done + n; prev_csv <- csv
    if (porosity_evolution) {
      clogged <- phi_raw <= phi_min
      phi <- max(phi_min, phi_raw)
      if (!is.null(k_of_phi)) k <- k_of_phi(phi)
    }
    refill_next <- if (refill) m_pw0 * phi / phi0 - w_ph else 0
    n_prev <- n; lm_prev <- ph$phreeqc_params$layer_mass
    n <- if (rel > 0) max(1, min(max_chunk, floor(n * phi_tol / rel))) else min(max_chunk, 2 * n)
  }
  list(chunks = do.call(rbind, rows), files = files)
}

# The ocean mass in a PlanetProfile profile's header (Mocean_kg).
profile_ocean_mass_kg <- function(profile_file) {
  h <- readLines(profile_file, n = 81, warn = FALSE)
  as.numeric(sub(".*=\\s*", "", grep("^\\s*Mocean_kg\\s*=", h, value = TRUE)[1]))
}

# The water the refill drew from the ocean (PI 2026-10-31), per chunk and cumulatively: on the run's normalised
# basis (kg per kg of the reacting box's rock), scaled to the whole core by the box's rock mass (the skin is a
# shell over the whole core, so layer_mass_kg is already whole-core; rebuilt from phi_start and k for runs made
# before that column), and as a fraction of the profile's ocean mass. The ocean's concentration factor from this
# water loss is M_ocean / (M_ocean - cumulative); the model's own ocean, water_ocean_kg x layer_mass_kg (real kg,
# at the chunk's end), gives M_ocean / that alongside (the two agree once restarts scale the ocean, 10-31).
refill_summary <- function(chunks, template = NULL, organic_wt_percent = 15,
                           profile_file = "profiles/EnceladusProfile_Seawater_10.0ppt_Tb272.4578K_PorousRock.txt") {
  if (!"layer_mass_kg" %in% names(chunks)) {
    if (is.null(template)) stop("refill_summary: give the template Simulator to rebuild layer_mass_kg")
    chunks$layer_mass_kg <- vapply(seq_len(nrow(chunks)), function(i) {
      s <- evolution_simulator(template, chunks$phi_start[i], chunks$k[i])
      invisible(capture.output(ph <- PhreeqcIntegrator$new(s), p <- ph$calculate_phreeqc_parameters(organic_wt_percent = organic_wt_percent)))
      p$layer_mass }, 0)
  }
  M <- profile_ocean_mass_kg(profile_file)
  real <- chunks$refill_in_kg * chunks$layer_mass_kg
  data.frame(chunk = chunks$chunk, orbit_from = chunks$orbit_from, refill_norm_kg = chunks$refill_in_kg,
             refill_core_kg = real, cumulative_core_kg = cumsum(real), cumulative_fraction_of_ocean = cumsum(real) / M,
             ocean_concentration_factor = M / (M - cumsum(real)),
             model_ocean_kg = chunks$water_ocean_kg * chunks$layer_mass_kg,
             model_ocean_concentration_factor = M / (chunks$water_ocean_kg * chunks$layer_mass_kg),
             layer_mass_kg = chunks$layer_mass_kg)
}

# Whole-system element totals of a tidal run's last orbit (rows REACT, MIX 4, MIX 1, MIX 2), split into the
# porewater side (porewater after MIX 1; rock, IOM and, under "retain", pore gas of the REACT row) and the ocean
# side (ocean after MIX 2; ocean minerals and gas of the MIX 4 row; under "to_ocean" the pore bubble is there).
whole_system_parts <- function(x, pq, fate) {
  x <- x[trimws(x$state) == "react", , drop = FALSE]
  el <- intersect(pq$totals, names(x)); v <- utils::read.csv(MINERAL_VOLUMES_CSV, stringsAsFactors = FALSE)
  sol <- function(i) { s <- stats::setNames(vapply(el, function(e) as.numeric(x[i, match(e, names(x))]) * as.numeric(x$mass_H2O[i]), 0), el)
    s[c("H", "O")] <- c(as.numeric(x$H_total_mol[i]), as.numeric(x$O_total_mol[i])); s }
  tot <- function(i) system_totals(x, i, el, pq$iom, v)[el]
  gas <- function(i) { g <- stats::setNames(rep(0, length(el)), el)
    for (k in grep("^g_", names(x), value = TRUE)) { fe <- formula_elements(sub("\\(g\\)$", "", sub("^g_", "", k)))
      for (e in intersect(names(fe), el)) g[e] <- g[e] + as.numeric(x[i, k]) * fe[[e]] }
    g }
  i <- nrow(x) - 3
  list(pore = sol(i + 2) + tot(i) - sol(i) - (if (identical(fate, "to_ocean")) gas(i) else 0),
       ocean = sol(i + 3) + tot(i + 1) - sol(i + 1))
}
# The t = 0 whole-system totals of a first chunk: its two initial solutions, the .pqi's minerals and IOM.
whole_system_t0 <- function(x, pq) {
  ini <- x[trimws(x$state) == "i_soln", , drop = FALSE]
  el <- intersect(pq$totals, names(x)); v <- utils::read.csv(MINERAL_VOLUMES_CSV, stringsAsFactors = FALSE)
  s <- stats::setNames(rep(0, length(el)), el)
  for (i in 1:2) { s <- s + stats::setNames(vapply(el, function(e) as.numeric(ini[i, match(e, names(ini))]) * as.numeric(ini$mass_H2O[i]), 0), el)
    s[c("H", "O")] <- s[c("H", "O")] + c(as.numeric(ini$H_total_mol[i]), as.numeric(ini$O_total_mol[i])) }
  for (p in names(pq$minerals)) { fe <- formula_elements(v$formula[match(p, v$phase)])
    for (e in intersect(names(fe), el)) s[e] <- s[e] + pq$minerals[[p]] * fe[[e]] }
  for (j in seq_len(nrow(pq$iom))) { fe <- formula_elements(pq$iom$formula[j])
    for (e in intersect(names(fe), el)) s[e] <- s[e] + pq$iom$m0[j] * fe[[e]] }
  s
}
