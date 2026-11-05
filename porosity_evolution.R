# Porosity that evolves with the chemistry (exchange-physics E, PI 2026-10-29). Sourced by
# ENTICES.R after phreeqc_tools.R.
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
# The core's drained frame moduli follow phi(t) (PI 2026-11-01): each chunk's Simulator has K_d = K_grain (1 -
# phi/phi_c), which sets the layered strain (K_d_Pa and mu_d_Pa columns). If phi reaches phi_c (alpha = 1, no
# frame) the run is flagged and stopped: a warning, and the chunks so far are returned with
# stopped = "phi_c".
#
# Long runs: grouped chunks with a geometric time step (production kit, item 1, PI 2026-11-01). With
# grouped = TRUE (the default when t_end_years is given) a chunk is steps_per_chunk PHREEQC steps of G orbits
# each (tidal_kinetic's grouped mode: the exchange compounds as x_step = 1 - (1 - x_orbit)^G), so the output
# is log-spaced in time. G starts at 1 and grows by g_growth per accepted chunk, limited by accuracy, not
# convenience:
#   - x_step <= x_step_max (0.1): G <= ln(1 - x_step_max) / ln(1 - x_orbit), from the chunk's own exchange;
#     this bounds the operator-splitting error of reacting the porewater for G orbits between exchanges;
#   - the porosity change per chunk within phi_tol: G <= phi_tol / (|dphi/phi| per orbit x steps_per_chunk),
#     from the last chunk's rate; a chunk that still exceeds phi_tol is redone at G/2;
#   - g_max (default Inf). There is no fixed cap: the 1e5-orbit cap of calculate_phreeqc_parameters()'s
#     automatic grouping (mode "ocean_circulations") is not used here, because the two criteria above already
#     bound the step and would be hidden by it (at the reference, x_step <= 0.1 allows ~1e6 orbits).
# Near the end the step shrinks so the run lands on the end time. End conditions, the first reached:
# t_end_years (default 4.5e9 in run_case.R); a steady state, ss_tol (relative change of the ocean's C(+4)
# (carbon_by_oxidation_state()), pH and Na over one ocean circulation time all below ss_tol); clogging
# (phi <= phi_min); phi reaching phi_c; max_chunks. `stopped` names it. Without t_end_years and grouped,
# n_orbits ends the run and chunks are ungrouped, as before.
#
# porosity_evolution: TRUE for production runs (default, PI 10-30, once the refill passed its checks);
# FALSE keeps phi0 (and no refill) and reports the phi the solids imply; the quick self-tests use both.

# One Simulator for porosity phi and permeability k, with everything else from the template simulator.
evolution_simulator <- function(template, phi, k, iom_vol_frac = template$params$iom_vol_frac) {
  s <- Simulator$new(k = k, layer_thickness = template$params$layer_thickness, config = template$config,
                     reactive_depth = template$params$reactive_depth, phi0 = phi,
                     phi_c = template$params$phi_c, alpha = template$params$alpha,
                     grain_diameter = template$params$grain_diameter, ssa = template$params$ssa,
                     roughness = template$params$roughness, grain_diameter_by_mineral = template$params$grain_diameter_by_mineral,
                     iom_wt_percent = template$params$iom_wt_percent, iom_vol_frac = iom_vol_frac)
  invisible(capture.output(s$run_simulation()))
  s
}

# template: a Simulator (k, layer, phi0, phi_c, config, reactive_depth). organic_wt_percent: as for
# calculate_phreeqc_parameters() (mineral_basis_wt_percent too: the blank's). n_orbits: total orbits. Further
# arguments go to both generators (porewater_temp, pore_gas_fate, carbon_mode, ...); trace_C_molal and first_args
# (e.g. list(initial_ocean = cometary_ice_ocean())) only to the first.
# Returns list(chunks = one row per accepted chunk, files = the chunks' .pqi, stopped, clog_time_yr).
#
# k_floor (PI 2026-11-02; NULL = none, the rigid-frame default): fracture-maintained permeability from
# reaction-driven cracking (Kelemen & Hirth 2012, EPSL 345-348, 81-89: the solid volume increase of hydration
# and carbonation can fracture the rock and keep or raise its permeability). With k_floor, k = max(k(phi),
# k_floor), and where the solids would clog the rock (phi <= phi_min) the run goes on with phi held at phi_min
# and k = k_floor (the further solid volume taken up by the cracking, not by the pores) instead of stopping.
# clog_time_yr: the first chunk end at which phi <= phi_min (NA if never), with or without k_floor.
#
# Resumable (production kit, item 2): after every accepted chunk the loop's state goes to state_file
# (<output_dir>/state.rds: time, phi, k, G, the chunk index, the drift, and what the next chunk needs) and the
# chunk table to <output_dir>/chunks.csv. resume = TRUE (default) with a state file present continues after
# its last chunk: the next chunk's directory is rebuilt from scratch, so a chunk killed half-way is redone and
# none is duplicated. A finished run returns at once, unless it ended at t_end or max_chunks and these were
# raised. The settings that shape the trajectory must match the state's (an error otherwise).
# log_file: one line per chunk (and per redo) appended there.
#
# follow (the same-rock blank, item 3): the output_dir of a leader run (the organic case). Each chunk takes the
# leader's phi, k, IOM volume fraction, G and orbits for the same chunk, waiting for the leader to finish it, so
# the two runs share transport, exchange and output times and differ only in their chemistry; it ends where the
# leader ended. No step control, steady-state or clogging test of its own.
run_tidal_evolution <- function(template, organic_wt_percent, n_orbits = NULL, output_dir, porosity_evolution = TRUE,
                                phi_tol = 0.01, phi_min = 0.01, k_of_phi = NULL, first_chunk = 1, max_chunk = 64,
                                trace_C_molal = 1e-14, refill = porosity_evolution,
                                t_end_years = NULL, grouped = !is.null(t_end_years), steps_per_chunk = 4, g_growth = 2,
                                x_step_max = 0.1, g_max = Inf, ss_tol = NULL, max_chunks = Inf,
                                k_floor = NULL, mineral_basis_wt_percent = NULL, resume = TRUE,
                                state_file = file.path(output_dir, "state.rds"), log_file = NULL, follow = NULL,
                                follow_poll_s = 5, first_args = list(), ...) {
  quiet <- function(e) { invisible(capture.output(r <- suppressMessages(e))); r }
  logf <- function(fmt, ...) if (!is.null(log_file)) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S "), sprintf(fmt, ...), "\n", sep = "", file = log_file, append = TRUE)
  phi0 <- template$params$porosity; phi <- phi0
  yr_per_orbit <- template$config$orbital_period / (365.25 * 24 * 3600)
  if (is.null(n_orbits) && is.null(t_end_years)) stop("run_tidal_evolution: give n_orbits or t_end_years")
  n_tot <- if (!is.null(t_end_years)) ceiling(t_end_years / yr_per_orbit) else n_orbits
  S <- steps_per_chunk; G_try <- 1
  # a Kozeny-Carman template follows phi(t) through k unless k_of_phi is given (item 3)
  if (is.null(k_of_phi) && identical(template$params$permeability, "kozeny_carman"))
    k_of_phi <- function(phi) kozeny_carman_k(template$params$grain_diameter, phi)
  k_at <- function(phi) { k <- if (is.null(k_of_phi)) template$params$k else k_of_phi(phi); if (is.null(k_floor)) k else max(k, k_floor) }
  k <- k_at(phi)
  done <- 0; n <- first_chunk; prev_csv <- NULL; origin <- NULL; V_total <- NA; clogged <- FALSE
  m_pw0 <- NA; w_ph_prev <- NA; n_prev <- NA; lm_prev <- NA; pq0 <- NULL; parts_prev <- NULL; drift <- NULL
  rows <- list(); files <- character(0); chunk <- 0; stopped <- NA_character_; clog_time <- NA_real_
  # the IOM's share of the solid volume for the frame moduli: t = 0 from the template, then from each chunk's
  # remaining IOM at the Delarue density (bulk_density()); the minerals' volume is the reconstruction's
  iom_vf <- template$params$iom_vol_frac

  # ---- checkpoint and resume ----
  env <- environment()
  st_vars <- c("phi", "k", "done", "n", "prev_csv", "origin", "V_total", "clogged", "m_pw0", "w_ph_prev", "n_prev",
               "lm_prev", "pq0", "parts_prev", "drift", "rows", "files", "chunk", "stopped", "clog_time", "iom_vf", "G_try")
  settings <- list(organic_wt_percent = organic_wt_percent, mineral_basis_wt_percent = mineral_basis_wt_percent,
                   porosity_evolution = porosity_evolution, phi_tol = phi_tol, phi_min = phi_min, grouped = grouped,
                   steps_per_chunk = S, g_growth = g_growth, x_step_max = x_step_max, g_max = g_max, k_floor = k_floor,
                   trace_C_molal = trace_C_molal, refill = refill, follow = !is.null(follow), generator_args = list(...),
                   first_args = first_args)
  rel_path <- function(p) if (is.null(p)) NULL else file.path(basename(dirname(p)), basename(p))
  abs_path <- function(p) if (is.null(p)) NULL else file.path(output_dir, p)
  save_state <- function() {
    s <- mget(st_vars, envir = env)
    s$prev_csv <- rel_path(prev_csv); s$origin <- rel_path(origin); s$files <- vapply(files, rel_path, "", USE.NAMES = FALSE)
    s$settings <- settings; s$time_yr <- done * yr_per_orbit; s$G <- if (length(rows)) rows[[length(rows)]]$G else NA
    s$saved <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    tmp <- paste0(state_file, ".tmp"); saveRDS(s, tmp)
    if (!file.rename(tmp, state_file)) { file.copy(tmp, state_file, overwrite = TRUE); unlink(tmp) }
    utils::write.csv(do.call(rbind, rows), file.path(output_dir, "chunks.csv"), row.names = FALSE)
  }
  result <- function() list(chunks = do.call(rbind, rows), files = files, stopped = stopped, clog_time_yr = clog_time)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  if (resume && file.exists(state_file)) {
    s <- readRDS(state_file)
    same <- function(a, b) isTRUE(all.equal(a, b))
    bad <- names(settings)[!vapply(names(settings), function(nm) same(settings[[nm]], s$settings[[nm]]), TRUE)]
    if (length(bad)) stop(sprintf("run_tidal_evolution: %s holds a run made with different settings (%s); give a new output_dir, or delete the old run to start again",
                                  state_file, paste(bad, collapse = ", ")))
    list2env(s[st_vars], envir = env)
    prev_csv <- abs_path(s$prev_csv); origin <- abs_path(s$origin); files <- vapply(s$files, abs_path, "", USE.NAMES = FALSE)
    if (stopped %in% c("t_end", "max_chunks")) stopped <- NA_character_   # the end was raised: go on
    logf("resume after chunk %d (t = %.6g yr, phi %.5g, k %.4g)%s", chunk, done * yr_per_orbit, phi, k,
         if (is.na(stopped)) "" else paste0("; already ended: ", stopped))
    if (!is.na(stopped)) return(result())
  } else logf("start: %s", output_dir)
  pv <- phreeqc_version_check()   # before chunk 1: a version mismatch warns now, not after hours
  logf("PHREEQC %s%s", pv$version, if (pv$ok) "" else sprintf(" (WARNING: differs from %s)", PHREEQC_EXPECTED_VERSION))
  leader_row <- function(i) {
    lf <- file.path(follow, "state.rds"); dead <- liveness_watch(follow, file.path(follow, "lock"))
    repeat {
      L <- if (file.exists(lf)) tryCatch(readRDS(lf), error = function(e) NULL) else NULL
      if (!is.null(L) && length(L$rows) >= i) return(L$rows[[i]])
      # A leader that ended at an earlier t_end is being extended (t_end raised): wait for its next rows
      if (!is.null(L) && !is.na(L$stopped) && !(L$stopped == "t_end" && L$done < n_tot)) return(structure(list(), stopped = L$stopped))
      if (file.exists(file.path(follow, "FAILED"))) stop("run_tidal_evolution: the leader run failed (", follow, ")")
      if (dead()) stop("run_tidal_evolution: the leader run (", follow, ") is not running and has not finished; rerun the case to resume both")
      Sys.sleep(follow_poll_s)
    }
  }

  while (done < n_tot && chunk < max_chunks) {
    if (!is.null(follow)) {
      lr <- leader_row(chunk + 1)
      if (!length(lr)) { stopped <- attr(lr, "stopped"); break }
      phi <- lr$phi_start; k <- lr$k; iom_vf <- lr$iom_vol_frac; G <- lr$G; n <- lr$n_orbits
    }
    if (phi >= template$params$phi_c) {
      stopped <- "phi_c"
      warning(sprintf("run_tidal_evolution: after chunk %d (orbit %g) phi = %.4g has reached phi_c = %.4g: no frame (alpha = 1, K_d = 0); stopped",
                      chunk, done, phi, template$params$phi_c))
      break
    }
    sim <- evolution_simulator(template, phi, k, iom_vf)
    if (is.null(follow)) {
      if (grouped) {
        x_orb <- sim$exchange$x_orbit
        G_x <- if (x_orb > 0) max(1, floor(log(1 - x_step_max) / log(1 - x_orb))) else Inf
        rem <- n_tot - done
        G <- max(1, min(G_try, G_x, g_max, floor(rem / S)))
        n <- if (rem < S) rem else S * G
      } else {
        G <- 1; n <- min(n, n_tot - done)
      }
    }
    ph <- quiet(PhreeqcIntegrator$new(sim))
    quiet(ph$calculate_phreeqc_parameters(organic_wt_percent = organic_wt_percent, mineral_basis_wt_percent = mineral_basis_wt_percent))
    ocean_scale <- if (is.na(lm_prev)) 1 else lm_prev / ph$phreeqc_params$layer_mass
    # refill (PI 10-30): bring the porewater back to the pore volume phi implies, from the ocean
    refill_next <- if (refill && !is.na(w_ph_prev)) m_pw0 * phi / phi0 - w_ph_prev else 0
    # short chunk directories: the generators' file names are long, and Windows refuses paths over 260
    # characters (PHREEQC then cannot open its DUMP file); the orbits are in the chunk table and the .pqi header
    dir <- file.path(output_dir, if (grouped || !is.null(follow)) sprintf("c%04d", chunk + 1) else sprintf("chunk%03d_orb%g-%g", chunk + 1, done + 1, done + n))
    if (.Platform$OS.type == "windows" && nchar(normalizePath(dir, mustWork = FALSE)) > 120)
      stop(sprintf("run_tidal_evolution: the output path %s is %d characters long; with the generated file names (~130 characters) it would pass Windows' 260-character limit. Use a shorter output directory (e.g. clone the repository nearer the drive root).",
                   normalizePath(dir, mustWork = FALSE), nchar(normalizePath(dir, mustWork = FALSE))))
    gm <- if (G > 1) list(mode = "grouped", grouping_factor = G) else list(mode = "ungrouped")
    unlink(dir, recursive = TRUE)
    f <- suppressWarnings(quiet(if (is.null(prev_csv))
      do.call(ph$generate_phreeqc_input_tidal_kinetic, c(list(output_dir = dir, n_orbital_cycles = n, trace_C_molal = trace_C_molal), gm, first_args, list(...)))
    else
      do.call(ph$generate_phreeqc_input_tidal_kinetic_restart, c(list(output_dir = dir, prev_csv = prev_csv, n_orbital_cycles = n,
                                                      allow_pressure_change = TRUE, stop_on_depletion = FALSE,
                                                      refill_kg = refill_next, ocean_scale = ocean_scale), gm, list(...)))))
    db <- sub("^# Database: ", "", grep("^# Database: ", readLines(f, n = 30), value = TRUE)[1])
    r <- run_phreeqc(f, db, output = sub("\\.pqi$", ".out", f), workdir = dirname(f), warn_version = FALSE)
    if (r$status != 0 || r$errors > 0) stop(sprintf("chunk %d (%s): PHREEQC status %s, %d errors", chunk + 1, basename(f), r$status, r$errors))
    csv <- sub("\\.pqi$", ".csv", f)
    bd <- suppressWarnings(bulk_density(csv, porosity = phi, origin_pqi = if (is.null(origin)) f else origin))$summary
    if (is.null(origin)) V_total <- bd$solid_volume_cm3[1] / (1 - phi0)
    V_end <- bd$solid_volume_cm3[nrow(bd)]
    phi_raw <- 1 - V_end / V_total
    # the porosity change the step control sees: with k_floor, phi is held at phi_min once the rock clogs
    rel <- abs((if (is.null(k_floor)) phi_raw else max(phi_raw, phi_min)) - phi) / phi
    if (is.null(follow) && porosity_evolution && rel > phi_tol && grouped && G > 1) {   # too big a step: redo
      logf("chunk %d: redo at G %g (porosity change %.3g > phi_tol %g at G %g)", chunk + 1, max(1, floor(G / 2)), rel, phi_tol, G)
      G_try <- max(1, floor(G / 2)); next
    }
    if (is.null(follow) && porosity_evolution && rel > phi_tol && !grouped && n > 1) { n <- max(1, floor(n / 2)); next }
    chunk <- chunk + 1; files <- c(files, f)
    if (is.null(origin)) origin <- f
    x <- read.delim(csv, check.names = FALSE, strip.white = TRUE); names(x) <- trimws(names(x))
    rr <- which(trimws(x$state) == "react" & trimws(x$soln) == "1")
    w_ph <- as.numeric(x$mass_H2O[rr[length(rr)]])          # the porewater after the last exchange (MIX 1)
    ro <- which(trimws(x$state) == "react" & trimws(x$soln) == "2")
    w_oc <- as.numeric(x$mass_H2O[ro[length(ro)]])          # the ocean (MIX 2)
    # the ocean at the chunk's end, for the steady-state test: pH, Na and C(+4) by actual oxidation state
    io <- ro[length(ro)]
    cb <- carbon_by_oxidation_state(csv)$by_state
    oc_C4 <- sum(cb$C_molal[cb$block == io & cb$ox == 4 & !cb$decoupled])
    if (is.na(m_pw0)) m_pw0 <- as.numeric(x$mass_H2O[which(trimws(x$state) == "i_soln" & trimws(x$soln) == "1")[1]])
    if (is.na(m_pw0)) {   # no t = 0 rows (a cometary initial ocean): the porewater SOLUTION's water in the .pqi
      p <- readLines(f); s1 <- grep("^SOLUTION 1\\b", p)[1]
      m_pw0 <- as.numeric(sub("^\\s*water\\s+", "", grep("^\\s*water\\s", p[s1:(s1 + 15)], value = TRUE)[1]))
    }
    # drift: this chunk's end against its start
    if (is.null(pq0)) pq0 <- pqi_initial_solids(f)
    fate <- sub("^# pore_gas_fate: ", "", grep("^# pore_gas_fate: ", readLines(f, n = 80), value = TRUE)[1])
    parts <- whole_system_parts(x, pq0, fate)
    start <- if (is.null(parts_prev)) whole_system_t0(x, pq0) else parts_prev$pore + ocean_scale * parts_prev$ocean
    end <- parts$pore + parts$ocean
    # a cometary initial ocean is made in a separate simulation, so the run's CSV has no t = 0 solutions: the
    # drift is then counted from the end of the first chunk
    if (anyNA(start)) { logf("chunk %d: no t = 0 solutions in the CSV (cometary initial ocean): drift counted from this chunk's end", chunk); start <- end }
    drift <- (if (is.null(drift)) 0 else drift) + (end - start)
    rel_d <- ifelse(abs(end) > 0, abs(drift) / abs(end), 0)
    if (max(rel_d) > 1e-4) warning(sprintf("run_tidal_evolution: chunk %d: cumulative drift of %s is %.2e relative", chunk, names(rel_d)[which.max(rel_d)], max(rel_d)))
    parts_prev <- parts
    xch <- sim$exchange
    pc <- profile_conditions(ph$profile_file, sim$reactive_thickness())   # the skin's profile mean T and P
    rows[[chunk]] <- data.frame(chunk = chunk, orbit_from = done + 1, orbit_to = done + n, n_orbits = n,
      phi_start = phi, phi_end_solids = phi_raw, dphi_rel = rel, over_tol = rel > phi_tol, k = k, alpha = sim$biot_alpha(),
      delta_m = xch$delta, L_eff_m = sim$reactive_thickness(), x_half = xch$x_half, binding = xch$binding,
      circulation_yr = quiet(sim$get_summary_metrics())$circulation_time_years, clogged = clogged,
      V_solid_cm3 = V_end, V_total_cm3 = V_total, water_phreeqc_kg = w_ph, water_ocean_kg = w_oc,
      water_from_phi_kg = max(phi_raw, 0) * V_total / 1000, water_target_kg = m_pw0 * max(phi_min, phi_raw) / phi0,
      refill_in_kg = refill_next, refill_kg_per_orbit = if (is.na(n_prev)) 0 else refill_next / n_prev,
      exchange_kg_per_orbit = xch$x_orbit * ph$phreeqc_params$porewater_mass / ph$phreeqc_params$layer_mass,
      layer_mass_kg = ph$phreeqc_params$layer_mass, ocean_scale = ocean_scale,
      K_d_Pa = sim$frame_moduli()$K_d, mu_d_Pa = sim$frame_moduli()$mu_d, eps_v = xch$eps_v,
      iom_vol_frac = sim$iom_volume_fraction(),
      ocean_per_rock = w_oc, ocean_per_rock_geometric = ph$phreeqc_params$ocean_mass / ph$phreeqc_params$layer_mass,
      drift_max = max(rel_d), drift_element = names(rel_d)[which.max(rel_d)],
      G = G, steps = if (G > 1) ceiling(n / G) else n, x_step = 1 - (1 - xch$x_orbit)^G,
      time_end_yr = (done + n) * yr_per_orbit, ocean_pH = as.numeric(x$pH[io]), ocean_Na = as.numeric(x$Na[io]),
      ocean_C4 = oc_C4, ocean_H2 = as.numeric(x$m_H2[io]), steady = FALSE,
      k_floor_on = !is.null(k_floor) && clogged,
      rho_solid = bd$rho_solid[nrow(bd)], rho_bulk = bd$rho_bulk[nrow(bd)], rho_bulk_in_bagheri = bd$core_density_ok[nrow(bd)],
      skin_T_C = pc$pore_temp_C, skin_P_atm = pc$pore_pressure_atm, phreeqc_version = r$version, seconds = r$seconds)
    done <- done + n; prev_csv <- csv; w_ph_prev <- w_ph
    # steady state: the ocean's C(+4), pH and Na each within ss_tol (relative) of their values one circulation
    # time earlier (the last chunk end at or before t - tau)
    if (!is.null(ss_tol) && is.null(follow)) {
      tnow <- rows[[chunk]]$time_end_yr; tau <- rows[[chunk]]$circulation_yr
      past <- Filter(function(r) r$time_end_yr <= tnow - tau, rows)
      if (length(past)) {
        p <- past[[length(past)]]; q <- rows[[chunk]]
        rc <- c(abs(q$ocean_C4 - p$ocean_C4) / abs(q$ocean_C4), abs(q$ocean_pH - p$ocean_pH) / abs(q$ocean_pH),
                abs(q$ocean_Na - p$ocean_Na) / abs(q$ocean_Na))
        if (all(is.finite(rc)) && all(rc < ss_tol)) { rows[[chunk]]$steady <- TRUE; stopped <- "steady state" }
      }
    }
    if (porosity_evolution && is.null(follow)) {
      clogged <- phi_raw <= phi_min
      if (clogged && is.na(clog_time)) clog_time <- rows[[chunk]]$time_end_yr
      phi <- max(phi_min, phi_raw)
      k <- if (clogged && !is.null(k_floor)) k_floor else k_at(phi)
    }
    nb <- nrow(bd); iom_vf <- iom_volume_fraction(template$params$iom_wt_percent, iom_g = bd$iom_g[nb], iom_density = bd$iom_density[nb])
    n_prev <- n; lm_prev <- ph$phreeqc_params$layer_mass
    if (is.null(follow)) {
      if (grouped) {
        # (with porosity_evolution = FALSE phi is frozen, so only x_step limits G)
        G_phi <- if (porosity_evolution && rel > 0) floor(phi_tol / (rel / n * S)) else Inf
        G_try <- max(1, min(G * g_growth, G_phi))
      } else n <- if (rel > 0) max(1, min(max_chunk, floor(n * phi_tol / rel))) else min(max_chunk, 2 * n)
      if (grouped && clogged && is.null(k_floor)) stopped <- "clogging"
    }
    if (is.na(stopped) && !(done < n_tot && chunk < max_chunks))
      stopped <- if (done >= n_tot) "t_end" else "max_chunks"
    logf("chunk %d: orbits %.0f-%.0f, G %g, t = %.6g yr, phi %.5g -> %.5g, k %.4g, x_step %.3g, drift %.2g (%s), %.1f s%s",
         chunk, done - n + 1, done, G, done * yr_per_orbit, rows[[chunk]]$phi_start, phi_raw, rows[[chunk]]$k,
         rows[[chunk]]$x_step, max(rel_d), names(rel_d)[which.max(rel_d)], r$seconds,
         paste0(if (clogged) "; clogged" else "", if (rows[[chunk]]$k_floor_on) " (k at k_floor)" else "",
                if (!is.na(stopped)) paste0("; ended: ", stopped) else ""))
    save_state()
    if (!is.na(stopped)) break
  }
  if (is.na(stopped)) stopped <- if (done >= n_tot) "t_end" else if (chunk >= max_chunks) "max_chunks" else NA_character_
  if (length(rows) && !identical(readRDS(state_file)$stopped, stopped)) save_state()
  logf("end: %s after %d chunks, t = %.6g yr", stopped, chunk, done * yr_per_orbit)
  result()
}

# Whether process `pid` is running: TRUE, FALSE (the process table answered and it is not there), or NA (no
# answer: tasklist failed or printed nothing, as it can on a loaded machine; tried 3 times).
pid_state <- function(pid) {
  if (length(pid) != 1 || is.na(pid)) return(FALSE)
  if (.Platform$OS.type != "windows") return(isTRUE(tools::pskill(pid, 0L)))
  for (i in 1:3) {
    o <- suppressWarnings(system2("tasklist", c("/FI", shQuote(sprintf("PID eq %d", pid)), "/NH", "/FO", "CSV"), stdout = TRUE, stderr = FALSE))
    if (any(grepl(sprintf('"%d"', pid), o, fixed = TRUE))) return(TRUE)
    if (is.null(attr(o, "status")) && any(nzchar(trimws(o)))) return(FALSE)   # answered (the "no tasks" line): not running
    Sys.sleep(1)
  }
  NA
}
# The run locks of run_case.R: an unanswered check counts as running, so a busy machine never frees a live lock.
pid_alive <- function(pid) !identical(pid_state(pid), FALSE)

# Seconds since anything in run directory `dir` last changed: its log, console log, state and chunk table, and the
# files of its newest chunk directory (PHREEQC's .out grows during a chunk). Inf if there is nothing.
run_activity_age <- function(dir) {
  f <- file.path(dir, c("run.log", "console.log", "state.rds", "chunks.csv", "lock", "run_case.console.log"))
  cd <- sort(list.files(dir, "^c[0-9]{4}$", full.names = TRUE), decreasing = TRUE)[1]
  if (!is.na(cd)) f <- c(f, list.files(cd, full.names = TRUE))
  m <- file.mtime(f[file.exists(f)])
  if (!length(m)) Inf else as.numeric(difftime(Sys.time(), max(m), units = "secs"))
}

# A watchdog for a run that should be alive (PI 2026-11-02: a busy shared machine must never stop a healthy run).
# Each call is one check; it returns TRUE only when the run has looked stopped on n_checks consecutive checks
# spanning at least span_s seconds. A check looks stopped only if the file activity under `dirs` is older than
# stale_s (the primary signal) and the process table says the lock's process is gone (tasklist answered "no";
# no answer counts as alive). A dead run is reported ~stale_s after its last write.
liveness_watch <- function(dirs, lock, stale_s = 600, n_checks = 6, span_s = 60) {
  n <- 0; first <- NULL
  function() {
    age <- min(vapply(dirs, run_activity_age, 0))
    pid <- if (file.exists(lock)) suppressWarnings(as.integer(readLines(lock, n = 1, warn = FALSE))[1]) else NA
    if (age <= stale_s || !identical(pid_state(pid), FALSE)) { n <<- 0; first <<- NULL; return(FALSE) }
    n <<- n + 1; if (is.null(first)) first <<- Sys.time()
    n >= n_checks && as.numeric(difftime(Sys.time(), first, units = "secs")) >= span_s
  }
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
