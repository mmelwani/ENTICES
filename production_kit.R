# production_kit.R -- the production kit (PI 2026-11-01/02): case files, one-command runs of a case and its
# same-rock blank, resumable unattended runs, and the per-case summary. Used by run_case.R, run_sweep.R and
# make_figures.R; sourced after ENTICES_LF_20260917.R. PRODUCTION_GUIDE.md is the user's guide; README section 3
# (Production runs) the reference.
#
# A case is an R file in cases/ that assigns a list `case` (cases/reference.R documents every field). A case
# may set `inherits = "<name>"` (a file cases/<name>.R; for cases/sweep/*.R, "reference" is cases/reference.R)
# and override fields: the parent's list is loaded first and modifyList()ed with the child's.

# ---- case files ----

KIT_CASE_DIR <- "cases"

# The parent of an inheriting case: cases/<name>.R, or a sibling file <dir>/<name>.R.
case_parent_path <- function(name, child_path) {
  cand <- c(file.path(dirname(child_path), paste0(name, ".R")), file.path(KIT_CASE_DIR, paste0(name, ".R")))
  hit <- cand[file.exists(cand)]
  if (!length(hit)) stop(sprintf("case %s: inherits = \"%s\", but no %s exists", child_path, name, paste(unique(cand), collapse = " or ")))
  hit[1]
}

# Load a case file (with inheritance) into a list; `name` defaults to the file name.
load_case <- function(path, depth = 0) {
  if (!file.exists(path)) stop("no case file ", path)
  if (depth > 10) stop("case inheritance deeper than 10 levels (a loop?) at ", path)
  e <- new.env(parent = globalenv())
  sys.source(path, envir = e)
  if (!exists("case", envir = e, inherits = FALSE) || !is.list(e$case)) stop(sprintf("case file %s must assign a list called `case`", path))
  child <- e$case
  out <- if (!is.null(child$inherits)) {
    parent <- load_case(case_parent_path(child$inherits, path), depth + 1)
    parent$inherits <- NULL; parent$name <- NULL; parent$label <- ""; parent$description <- NULL
    utils::modifyList(parent, child, keep.null = TRUE)
  } else child
  if (depth == 0) {
    if (is.null(out$name)) out$name <- tools::file_path_sans_ext(basename(path))
    out$case_file <- path
  }
  out
}

# Every field a case may set, with the reference value (cases/reference.R repeats these with comments).
KIT_CASE_DEFAULTS <- list(
  name = NULL, description = "", label = "", inherits = NULL,
  layer_thickness_m = 194200, reactive_depth = "skin", porewater_temp = "profile",
  phi0 = 0.32, phi_c = 0.4, phi_min = 0.01, grain_diameter_m = 0.0159, k_m2 = NULL, k_floor_m2 = NULL,
  core_viscosity_Pa_s = 1e12, tidal_scenario = "dissipative_core",
  iom_wt_percent = 15, carbon_mode = "c4_only", initial_ocean = "pure_water", cometary_f = 1,
  pore_gas_fate = "to_ocean", trace_C_molal = 1e-14, phosphate_phases = "apatite",
  t_end_yr = 1e8, ss_tol = 0.01,
  phi_tol = 0.003, x_step_max = 0.03, steps_per_chunk = 4,
  generator = "tidal_kinetic", blank = TRUE
)

# Check a case and fill in defaults. All problems are collected and reported together, in plain words.
validate_case <- function(case) {
  unknown <- setdiff(names(case), c(names(KIT_CASE_DEFAULTS), "case_file"))
  p <- character(0)
  if (length(unknown)) p <- c(p, sprintf("unknown field(s): %s (see cases/reference.R for the list)", paste(unknown, collapse = ", ")))
  for (nm in names(KIT_CASE_DEFAULTS)) if (!nm %in% names(case)) case[nm] <- list(KIT_CASE_DEFAULTS[[nm]])
  for (nm in c("ss_tol", "k_floor_m2", "k_m2", "grain_diameter_m", "trace_C_molal"))   # NA means "none", as NULL
    if (length(case[[nm]]) == 1 && is.na(case[[nm]])) case[nm] <- list(NULL)
  num <- function(f, lo, hi, null_ok = FALSE, lo_open = FALSE) {
    v <- case[[f]]
    if (is.null(v) || (length(v) == 1 && is.na(v))) { if (!null_ok) p <<- c(p, sprintf("%s is missing", f)); return(invisible()) }
    if (!(is.numeric(v) && length(v) == 1 && is.finite(v)) || v > hi || (if (lo_open) v <= lo else v < lo))
      p <<- c(p, sprintf("%s = %s is outside the allowed range %s%g, %g]", f, format(v), if (lo_open) "(" else "[", lo, hi))
  }
  one_of <- function(f, allowed) if (!(is.character(case[[f]]) && length(case[[f]]) == 1 && case[[f]] %in% allowed))
    p <<- c(p, sprintf("%s = %s: must be one of %s", f, format(case[[f]]), paste(sprintf("\"%s\"", allowed), collapse = ", ")))
  if (!(is.character(case$name) && length(case$name) == 1 && grepl("^[A-Za-z0-9_.+-]+$", case$name)))
    p <- c(p, "name must be letters, digits and _ . + - only (it becomes the directory runs/<name>)")
  num("layer_thickness_m", 0, 194200, lo_open = TRUE)
  one_of("reactive_depth", c("skin", "layer"))
  if (is.numeric(case$porewater_temp)) num("porewater_temp", 0, 300) else one_of("porewater_temp", c("profile", "ocean"))
  num("phi0", 0, 1, lo_open = TRUE); num("phi_c", 0, 1, lo_open = TRUE); num("phi_min", 0, 0.2, lo_open = TRUE)
  if (is.numeric(case$phi0) && is.numeric(case$phi_c) && case$phi0 >= case$phi_c)
    p <- c(p, sprintf("phi0 = %g must be below phi_c = %g (at phi_c the rock has no frame)", case$phi0, case$phi_c))
  if (is.numeric(case$phi0) && is.numeric(case$phi_min) && case$phi_min >= case$phi0)
    p <- c(p, sprintf("phi_min = %g must be below phi0 = %g", case$phi_min, case$phi0))
  if (is.null(case$grain_diameter_m) == is.null(case$k_m2))
    p <- c(p, "give exactly one of grain_diameter_m (Kozeny-Carman k, following the porosity) and k_m2 (a fixed k)")
  num("grain_diameter_m", 1e-7, 0.1, null_ok = TRUE); num("k_m2", 1e-20, 1e-6, null_ok = TRUE)
  num("k_floor_m2", 1e-20, 1e-8, null_ok = TRUE)
  one_of("tidal_scenario", c("dissipative_core", "dissipative_shell"))
  if (is.numeric(case$core_viscosity_Pa_s) && is.character(case$tidal_scenario)) {
    ok <- tryCatch({ FlowUtils$layered_strain(case$tidal_scenario, case$core_viscosity_Pa_s); TRUE }, error = function(e) conditionMessage(e))
    if (!isTRUE(ok)) p <- c(p, sprintf("core_viscosity_Pa_s: %s", ok))
  } else p <- c(p, "core_viscosity_Pa_s must be a number")
  num("iom_wt_percent", 0, 30)
  if (isTRUE(case$blank) && is.numeric(case$iom_wt_percent) && case$iom_wt_percent == 0)
    p <- c(p, "iom_wt_percent = 0 with blank = TRUE: there are no organics to compare with a blank")
  one_of("carbon_mode", c("c4_only", "coupled", "ch4_uncoupled", "temperature"))
  one_of("initial_ocean", c("pure_water", "cometary"))
  if (identical(case$initial_ocean, "cometary")) num("cometary_f", 0, 1, lo_open = TRUE)
  one_of("pore_gas_fate", c("to_ocean", "retain", "retain_to_pore_volume"))
  one_of("phosphate_phases", c("apatite", "precursor"))
  num("trace_C_molal", 0, 1e-6, null_ok = TRUE)
  num("t_end_yr", 0, 4.5e9, lo_open = TRUE)
  num("ss_tol", 0, 0.5, null_ok = TRUE, lo_open = TRUE)
  num("phi_tol", 0, 0.2, lo_open = TRUE); num("x_step_max", 0, 1, lo_open = TRUE)
  num("steps_per_chunk", 1, 64)
  if (!identical(case$generator, "tidal_kinetic")) {
    known <- identical(case$generator, "tidal_equil") && identical(case$carbon_mode, "c4_only") &&
      identical(case$pore_gas_fate, "to_ocean") && identical(case$initial_ocean, "cometary")
    p <- c(p, if (known) "generator = \"tidal_equil\" with c4_only, to_ocean and a cometary ocean is the known PHREEQC 3.8.6 crash (tests/known_failures/README.md); long runs use tidal_kinetic"
              else sprintf("generator = %s: long runs use \"tidal_kinetic\" only", format(case$generator)))
  }
  if (!is.logical(case$blank) || length(case$blank) != 1) p <- c(p, "blank must be TRUE or FALSE")
  if (length(p)) stop(sprintf("case %s is not valid:\n  - %s", if (is.null(case$case_file)) case$name else case$case_file,
                              paste(p, collapse = "\n  - ")), call. = FALSE)
  case
}

# The template Simulator of a case.
case_template <- function(case) {
  cfg <- EnceladusConfig$new()
  cfg$core_viscosity <- case$core_viscosity_Pa_s
  cfg$tidal_scenario <- case$tidal_scenario
  Simulator$new(k = if (is.null(case$k_m2)) "kozeny_carman" else case$k_m2, grain_diameter = case$grain_diameter_m,
                layer_thickness = case$layer_thickness_m, config = cfg, reactive_depth = case$reactive_depth,
                phi0 = case$phi0, phi_c = case$phi_c, iom_wt_percent = case$iom_wt_percent)
}

# The run_tidal_evolution() arguments of one variant ("organic" or "blank") of a case.
case_run_args <- function(case, variant, run_dir, max_chunks = Inf) {
  out <- file.path(run_dir, variant)
  first <- if (identical(case$initial_ocean, "cometary")) list(initial_ocean = cometary_ice_ocean(f = case$cometary_f)) else list()
  list(template = case_template(case),
       organic_wt_percent = if (variant == "organic") case$iom_wt_percent else 0,
       mineral_basis_wt_percent = if (variant == "blank") case$iom_wt_percent else NULL,
       t_end_years = case$t_end_yr, output_dir = out, grouped = TRUE,
       phi_tol = case$phi_tol, phi_min = case$phi_min, x_step_max = case$x_step_max, steps_per_chunk = case$steps_per_chunk,
       ss_tol = case$ss_tol, k_floor = case$k_floor_m2, trace_C_molal = case$trace_C_molal, max_chunks = max_chunks,
       log_file = file.path(out, "run.log"), follow = if (variant == "blank") file.path(run_dir, "organic") else NULL,
       first_args = first, porewater_temp = case$porewater_temp, carbon_mode = case$carbon_mode, pore_gas_fate = case$pore_gas_fate,
       phosphate_phases = case$phosphate_phases)
}

kit_log <- function(file, fmt, ...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S "), sprintf(fmt, ...), "\n", sep = "", file = file, append = TRUE)

# The PHREEQC executable and the version(s) its runs reported (run_phreeqc(); the chunk table's
# phreeqc_version column). PHREEQC has no non-interactive --version (it waits for input), so before a run only
# the path is known.
phreeqc_version_string <- function(chunks = NULL) {
  exe <- tryCatch(phreeqc_exe(), error = function(e) "PHREEQC not found")
  v <- if (!is.null(chunks) && "phreeqc_version" %in% names(chunks)) unique(stats::na.omit(as.character(chunks$phreeqc_version))) else character(0)
  paste0(exe, if (length(v)) paste0(" (version ", paste(v, collapse = ", "), " as run)") else "")
}

# Run one variant in this process: a lock with this PID (another live process holding it is an error), the log,
# warnings into the log, FAILED with the message on an error. Resumes from the variant's state file.
run_variant <- function(case, variant, run_dir, max_chunks = Inf) {
  out <- file.path(run_dir, variant); dir.create(out, recursive = TRUE, showWarnings = FALSE)
  lock <- file.path(out, "lock"); log <- file.path(out, "run.log")
  if (file.exists(lock)) {
    pid <- suppressWarnings(as.integer(readLines(lock, n = 1, warn = FALSE)))
    if (!identical(pid, Sys.getpid()) && pid_alive(pid)) stop(sprintf("%s is already being run by process %d", out, pid), call. = FALSE)
  }
  writeLines(as.character(Sys.getpid()), lock); unlink(file.path(out, "FAILED"))
  on.exit(unlink(lock), add = TRUE)
  kit_log(log, "%s / %s: %s; PID %d; R %s; %s", case$name, variant, code_version(), Sys.getpid(), getRversion(), phreeqc_version_string())
  t0 <- Sys.time()
  res <- withCallingHandlers(
    tryCatch(do.call(run_tidal_evolution, case_run_args(case, variant, run_dir, max_chunks)),
             error = function(e) {
               kit_log(log, "ERROR: %s", conditionMessage(e))
               writeLines(conditionMessage(e), file.path(out, "FAILED"))
               stop(e)
             }),
    warning = function(w) { kit_log(log, "WARNING: %s", conditionMessage(w)); invokeRestart("muffleWarning") })
  kit_log(log, "%s / %s: finished (%s), %.0f s in this session", case$name, variant, res$stopped, as.numeric(difftime(Sys.time(), t0, units = "secs")))
  invisible(res)
}

# The status of a variant directory: "done" (its state has ended), "failed", "running" (its lock's process is
# alive), or "stopped" (none of these: never started, or killed).
variant_status <- function(out) {
  if (file.exists(file.path(out, "FAILED"))) return("failed")
  st <- file.path(out, "state.rds")
  if (file.exists(st)) { s <- tryCatch(readRDS(st), error = function(e) NULL); if (!is.null(s) && !is.na(s$stopped)) return("done") }
  lk <- file.path(out, "lock")
  if (file.exists(lk) && pid_alive(suppressWarnings(as.integer(readLines(lk, n = 1, warn = FALSE))))) return("running")
  "stopped"
}

# The fields that must match for runs/<name> to be resumed (end conditions and labels may change).
case_identity <- function(case) case[setdiff(names(case), c("t_end_yr", "ss_tol", "label", "description", "case_file"))]

# Run a case: validate, then the organic run and its blank (in two processes when parallel), then summary.csv,
# summary_units.csv and RUNINFO.txt. Resumes whatever runs/<name> already holds. runs_dir: where run dirs go.
run_case <- function(path, parallel = NA, max_chunks = Inf, runs_dir = "runs", name = NULL, rscript = file.path(R.home("bin"), "Rscript")) {
  case <- load_case(path); if (!is.null(name)) case$name <- name
  case <- validate_case(case)
  run_dir <- file.path(runs_dir, case$name); dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  if (.Platform$OS.type == "windows" && nchar(normalizePath(run_dir, mustWork = FALSE)) > 100)
    stop(sprintf("the run directory %s is %d characters long; Windows limits paths to 260 characters and the generated file names need ~150. Clone the repository nearer the drive root (e.g. C:/ENTICES).",
                 normalizePath(run_dir, mustWork = FALSE), nchar(normalizePath(run_dir, mustWork = FALSE))), call. = FALSE)
  saved <- file.path(run_dir, "case.rds")
  if (file.exists(saved)) {
    old <- case_identity(readRDS(saved)); new <- case_identity(case)
    changed <- union(setdiff(names(new), names(old)), Filter(function(nm) !isTRUE(all.equal(old[[nm]], new[[nm]])), intersect(names(new), names(old))))
    if (length(changed)) stop(sprintf("%s holds a run of a different case (changed: %s). Give the case a new name, or move the old run away.",
                                      run_dir, paste(changed, collapse = ", ")), call. = FALSE)
  }
  saveRDS(case, saved)
  file.copy(case$case_file, file.path(run_dir, "case_file.R"), overwrite = TRUE)
  lock <- file.path(run_dir, "lock")
  if (file.exists(lock)) {
    pid <- suppressWarnings(as.integer(readLines(lock, n = 1, warn = FALSE)))
    if (pid_alive(pid) && pid != Sys.getpid()) stop(sprintf("%s is already being run by process %d", run_dir, pid), call. = FALSE)
  }
  writeLines(as.character(Sys.getpid()), lock); on.exit(unlink(lock), add = TRUE)
  variants <- c("organic", if (isTRUE(case$blank)) "blank")
  if (is.na(parallel)) parallel <- length(variants) == 2 && isTRUE(parallel::detectCores() >= 2)
  t0 <- Sys.time()
  if (parallel && length(variants) == 2) {
    for (v in variants) {
      if (variant_status(file.path(run_dir, v)) %in% c("done", "running")) next
      dir.create(file.path(run_dir, v), showWarnings = FALSE)
      unlink(file.path(run_dir, v, "FAILED"))
      system2(rscript, c("run_case.R", "--variant", v, "--runs-dir", shQuote(runs_dir), "--name", case$name,
                         if (is.finite(max_chunks)) c("--max-chunks", max_chunks), shQuote(case$case_file)),
              wait = FALSE, stdout = file.path(run_dir, v, "console.log"), stderr = file.path(run_dir, v, "console.log"))
    }
    started <- Sys.time()
    repeat {
      st <- vapply(variants, function(v) variant_status(file.path(run_dir, v)), "")
      if (any(st == "failed")) {
        v <- variants[st == "failed"][1]
        stop(sprintf("the %s run failed: %s (log: %s)", v, paste(readLines(file.path(run_dir, v, "FAILED"), warn = FALSE), collapse = " "),
                     file.path(run_dir, v, "run.log")), call. = FALSE)
      }
      if (all(st == "done")) break
      if (any(st == "stopped") && difftime(Sys.time(), started, units = "secs") > 300)
        stop(sprintf("the %s run is not running and has not finished; see %s", variants[st == "stopped"][1],
                     file.path(run_dir, variants[st == "stopped"][1], "console.log")), call. = FALSE)
      Sys.sleep(5)
    }
  } else for (v in variants) if (variant_status(file.path(run_dir, v)) != "done") run_variant(case, v, run_dir, max_chunks)
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  s <- build_summary(run_dir, case)
  utils::write.csv(s, file.path(run_dir, "summary.csv"), row.names = FALSE)
  utils::write.csv(summary_units(names(s)), file.path(run_dir, "summary_units.csv"), row.names = FALSE)
  write_runinfo(run_dir, case, wall)
  invisible(s)
}

# ---- summary.csv ----

KIT_OCEAN_SPECIES <- list(pH = "pH", pe = "pe", DIC_C4 = NA, Mtg = "Mtg", H2 = "m_H2", NH4_NH3 = c("m_NH4+", "m_NH3"),
                          SO4 = "m_SO4-2", S_minus2 = c("m_HS-", "m_H2S"), Na = "Na", Cl = "Cl", K = "K", Ca = "Ca",
                          Mg = "Mg", Si = "Si", P = "P", HPO4 = "m_HPO4-2", PO4 = "m_PO4-3", ionic_strength = "mu")
KIT_PHOSPHATE_PHASES <- c("Hydroxyapatite", "Whitlockite", "Brushite", "OCP", "Bobbierite", "Newberyite", "Vivianite", "Struvite")
KIT_GASES <- c("H2(g)", "CO2(g)", "Mtg(g)", "CH4(g)", "H2S(g)", "N2(g)", "NH3(g)", "SO2(g)", "Ntg(g)")

# One row per chunk end of one variant: the ocean, the gases, the H2 budget (organic only when `h2` is TRUE).
variant_rows <- function(out, chunks, case, h2 = TRUE) {
  st <- readRDS(file.path(out, "state.rds")); files <- file.path(out, st$files)
  cfg <- case_template(case)$config
  fate <- case$pore_gas_fate; gen_cum <- 0
  do.call(rbind, lapply(seq_along(files), function(i) {
    csv <- sub("\\.pqi$", ".csv", files[i])
    x <- utils::read.delim(csv, check.names = FALSE, strip.white = TRUE); names(x) <- trimws(names(x))
    r <- x[trimws(x$state) == "react", , drop = FALSE]; j <- nrow(r) - 3   # REACT, MIX 4, MIX 1, MIX 2 of the last step
    oc <- r[j + 3, , drop = FALSE]
    val <- function(cols) { cols <- intersect(cols, names(oc)); if (!length(cols)) NA_real_ else sum(vapply(cols, function(cc) as.numeric(oc[[cc]][1]), 0)) }
    cb <- carbon_by_oxidation_state(csv)$by_state
    io <- which(trimws(x$state) == "react")[j + 3]
    row <- data.frame(chunk = i)
    for (s in names(KIT_OCEAN_SPECIES)) row[[paste0("ocean_", s)]] <- if (s == "DIC_C4")
      sum(cb$C_molal[cb$block == io & cb$ox == 4 & !cb$decoupled]) else val(KIT_OCEAN_SPECIES[[s]])
    # the porewater's phosphate (after the last exchange, MIX 1) and the phosphate phases held (REACT: pores; MIX 4: ocean)
    pw <- r[j + 2, , drop = FALSE]
    pcol <- c(P = "P", HPO4 = "m_HPO4-2", PO4 = "m_PO4-3")
    for (s in names(pcol)) row[[paste0("pore_", s)]] <- if (pcol[[s]] %in% names(pw)) as.numeric(pw[[pcol[[s]]]][1]) else NA_real_
    ep <- function(k, p) { w <- which(names(r) == p); if (length(w)) as.numeric(r[k, max(w)]) else 0 }
    for (p in KIT_PHOSPHATE_PHASES) { row[[paste0("pore_", p, "_mol")]] <- ep(j, p); row[[paste0("ocean_", p, "_mol")]] <- ep(j + 1, p) }
    gas <- function(k, g) { cc <- paste0("g_", g); if (cc %in% names(r)) as.numeric(r[[cc]][k]) else 0 }
    pg <- if (identical(fate, "retain_to_pore_volume")) {
      y <- utils::read.delim(sub("\\.csv$", "_poregas.csv", csv), check.names = FALSE, strip.white = TRUE); names(y) <- trimws(names(y)); y } else NULL
    for (g in KIT_GASES) {
      row[[paste0("ocean_gas_", sub("\\(g\\)", "", g), "_mol")]] <- gas(j + 1, g)
      row[[paste0("pore_gas_", sub("\\(g\\)", "", g), "_mol")]] <- switch(fate, retain = gas(j, g),
        retain_to_pore_volume = { cc <- paste0("g_", g); if (cc %in% names(pg)) as.numeric(pg[[cc]][nrow(pg)]) else 0 }, 0)
    }
    if (h2) {
      b <- utils::tail(h2_budget(csv), 1)
      gen_cum <<- gen_cum + b$H2_generated
      row$h2_generated_chunk <- b$H2_generated; row$h2_generated_cum <- gen_cum
      row$h2_present <- b$H2_present; row$h2_ocean_dissolved <- b$ocean_dissolved; row$h2_ocean_gas <- b$ocean_gas
      row$h2_porewater_dissolved <- b$porewater_dissolved; row$h2_pore_gas <- b$pore_gas
      row$h2_other_redox_cum <- gen_cum - b$H2_present
      # whole-core values: x the rock of the whole core (scale_to_core(), with this chunk's porosity and skin)
      sim <- list(config = cfg, params = list(porosity = chunks$phi_start[i], layer_thickness = case$layer_thickness_m),
                  reactive_thickness = function() chunks$L_eff_m[i])
      sc <- scale_to_core(data.frame(orbit = i, Time_Years = NA, h2_generated_cum = gen_cum, h2_present = b$H2_present,
                                     h2_ocean_dissolved = b$ocean_dissolved, h2_ocean_gas = b$ocean_gas,
                                     h2_porewater_dissolved = b$porewater_dissolved, h2_pore_gas = b$pore_gas), sim)
      for (cc in grep("_core_mol$", names(sc), value = TRUE)) row[[cc]] <- sc[[cc]]
      row$rock_box_kg <- attr(sc, "rock_kg")[["box"]]; row$rock_core_kg <- attr(sc, "rock_kg")[["core"]]
    }
    row
  }))
}

# summary.csv of a run directory: one row per output time (the end of each chunk, log-spaced in time).
build_summary <- function(run_dir, case) {
  oc <- utils::read.csv(file.path(run_dir, "organic", "chunks.csv"), stringsAsFactors = FALSE)
  so <- readRDS(file.path(run_dir, "organic", "state.rds"))
  o <- variant_rows(file.path(run_dir, "organic"), oc, case, h2 = TRUE)
  out <- data.frame(case = case$name, time_yr = oc$time_end_yr, orbit = oc$orbit_to, chunk = oc$chunk, G = oc$G)
  sp <- paste0("ocean_", names(KIT_OCEAN_SPECIES))
  out <- cbind(out, o[sp])
  if (isTRUE(case$blank) && file.exists(file.path(run_dir, "blank", "state.rds"))) {
    bc <- utils::read.csv(file.path(run_dir, "blank", "chunks.csv"), stringsAsFactors = FALSE)
    b <- variant_rows(file.path(run_dir, "blank"), bc, case, h2 = FALSE)
    n <- min(nrow(o), nrow(b))
    if (n < nrow(o)) warning(sprintf("build_summary: the blank has %d of the organic run's %d chunks; the rest have no blank", nrow(b), nrow(o)))
    for (s in sp) {
      bl <- c(b[[s]][seq_len(n)], rep(NA_real_, nrow(o) - n))
      out[[paste0("blank_", s)]] <- bl
      out[[paste0("org_minus_blank_", s)]] <- out[[s]] - bl
    }
  }
  gascols <- grep("^(ocean|pore)_gas_.*_mol$", names(o), value = TRUE)
  gascols <- gascols[vapply(gascols, function(cc) any(o[[cc]] != 0, na.rm = TRUE), TRUE)]   # gases that appear
  phcols <- paste0(rep(c("pore_", "ocean_"), each = length(KIT_PHOSPHATE_PHASES)), KIT_PHOSPHATE_PHASES, "_mol")
  phcols <- phcols[vapply(phcols, function(cc) any(o[[cc]] != 0, na.rm = TRUE), TRUE)]          # phases that formed
  out <- cbind(out, o[c("pore_P", "pore_HPO4", "pore_PO4")], o[phcols], o[gascols], o[grep("^h2_|^rock_", names(o), value = TRUE)])
  out$phi <- pmax(oc$phi_end_solids, case$phi_min); out$phi_solids <- oc$phi_end_solids; out$phi_start <- oc$phi_start
  out$k_m2 <- oc$k; out$delta_m <- oc$delta_m; out$L_eff_m <- oc$L_eff_m
  out$skin_T_C <- oc$skin_T_C; out$skin_P_atm <- oc$skin_P_atm
  out$porewater_T_C <- if (is.numeric(case$porewater_temp)) case$porewater_temp else if (case$porewater_temp == "profile") oc$skin_T_C else 274.6 - 273.15
  out$alpha <- oc$alpha; out$K_d_Pa <- oc$K_d_Pa; out$circulation_yr <- oc$circulation_yr
  out$x_orbit <- 1 - (1 - oc$x_half)^2; out$x_step <- oc$x_step
  out$exchange_kg_per_orbit <- oc$exchange_kg_per_orbit; out$refill_kg_per_orbit <- oc$refill_kg_per_orbit
  out$rho_solid <- oc$rho_solid; out$rho_bulk <- oc$rho_bulk; out$rho_bulk_in_bagheri <- oc$rho_bulk_in_bagheri
  out$clogged <- oc$phi_end_solids <= case$phi_min; out$k_floor_on <- oc$k_floor_on
  out$drift_max <- oc$drift_max; out$drift_element <- oc$drift_element; out$steady <- oc$steady
  out$end_reason <- so$stopped; out$clog_time_yr <- so$clog_time
  out$label <- case$label
  out
}

# summary_units.csv: unit and meaning of every summary column.
summary_units <- function(cols) {
  u <- list(
    case = c("-", "case name"), time_yr = c("yr", "model time at the end of the chunk"), orbit = c("-", "orbits since t = 0"),
    chunk = c("-", "chunk index"), G = c("orbits", "orbits per PHREEQC step in the chunk"),
    ocean_pH = c("-", "ocean pH (end of chunk)"), ocean_pe = c("-", "ocean pe"),
    ocean_DIC_C4 = c("mol/kgw", "ocean dissolved inorganic carbon: C(+4) by actual oxidation state (carbon_by_oxidation_state())"),
    ocean_Mtg = c("mol/kgw", "ocean decoupled methane (Mtg total)"), ocean_H2 = c("mol/kgw", "ocean dissolved H2"),
    ocean_NH4_NH3 = c("mol/kgw", "ocean NH4+ + NH3"), ocean_SO4 = c("mol/kgw", "ocean SO4-2 (free species)"),
    ocean_S_minus2 = c("mol/kgw", "ocean HS- + H2S"), ocean_Na = c("mol/kgw", "ocean total Na"), ocean_Cl = c("mol/kgw", "ocean total Cl"),
    ocean_K = c("mol/kgw", "ocean total K"), ocean_Ca = c("mol/kgw", "ocean total Ca"), ocean_Mg = c("mol/kgw", "ocean total Mg"),
    ocean_Si = c("mol/kgw", "ocean total Si"), ocean_P = c("mol/kgw", "ocean total P"), ocean_ionic_strength = c("mol/kgw", "ocean ionic strength"),
    ocean_HPO4 = c("mol/kgw", "ocean HPO4-2 (free species)"), ocean_PO4 = c("mol/kgw", "ocean PO4-3 (free species)"),
    pore_P = c("mol/kgw", "porewater total P (after the last exchange)"), pore_HPO4 = c("mol/kgw", "porewater HPO4-2 (free species)"),
    pore_PO4 = c("mol/kgw", "porewater PO4-3 (free species)"),
    h2_generated_chunk = c("mol/kg box rock", "H2 generated in the chunk by electron balance (h2_budget())"),
    h2_generated_cum = c("mol/kg box rock", "H2 generated since t = 0 (sum of the chunks; each on its own box's rock)"),
    h2_present = c("mol/kg box rock", "H2 present: ocean dissolved + ocean gas + porewater + pore gas"),
    h2_ocean_dissolved = c("mol/kg box rock", "H2 dissolved in the ocean"), h2_ocean_gas = c("mol/kg box rock", "H2 in the ocean gas phase"),
    h2_porewater_dissolved = c("mol/kg box rock", "H2 dissolved in the porewater"), h2_pore_gas = c("mol/kg box rock", "H2 in the pore gas"),
    h2_other_redox_cum = c("mol/kg box rock", "H2 generated minus present: taken up by the other redox couples"),
    rock_box_kg = c("kg", "rock in the reacting box (the skin)"), rock_core_kg = c("kg", "rock in the whole core"),
    phi = c("-", "porosity at the chunk's end (held at phi_min once clogged)"), phi_solids = c("-", "porosity the solids imply (may go below phi_min)"),
    phi_start = c("-", "porosity during the chunk"), k_m2 = c("m2", "permeability during the chunk"),
    delta_m = c("m", "poroelastic skin depth"), L_eff_m = c("m", "thickness of the reacting box, min(layer, skin depth)"),
    skin_T_C = c("degC", "profile mean temperature over the reacting box"), skin_P_atm = c("atm", "profile mean pressure over the reacting box"),
    porewater_T_C = c("degC", "porewater reaction temperature used"), alpha = c("-", "Biot coefficient"),
    K_d_Pa = c("Pa", "drained frame bulk modulus"), circulation_yr = c("yr", "ocean circulation time (ocean mass / exchange per orbit)"),
    x_orbit = c("-", "fraction of the pore water exchanged per orbit"), x_step = c("-", "fraction exchanged per PHREEQC step (G orbits)"),
    exchange_kg_per_orbit = c("kg/kg box rock", "tidal exchange per orbit"), refill_kg_per_orbit = c("kg/kg box rock", "ocean water drawn into the pores per orbit by hydration"),
    rho_solid = c("g/cm3", "solid density"), rho_bulk = c("g/cm3", "bulk density"),
    rho_bulk_in_bagheri = c("-", "rho_bulk within Bagheri et al. (2025)'s core density 2.16-2.40 g/cm3"),
    clogged = c("-", "phi_solids <= phi_min at the chunk's end"), k_floor_on = c("-", "k held at k_floor during the chunk (clogged, with k_floor)"),
    drift_max = c("-", "largest cumulative element drift, relative"), drift_element = c("-", "the element with drift_max"),
    steady = c("-", "steady state reached at this chunk (ocean C(+4), pH, Na within ss_tol over one circulation time)"),
    end_reason = c("-", "what ended the run: steady state, t_end, clogging, phi_c (or max_chunks for a test run)"),
    clog_time_yr = c("yr", "first chunk end with phi_solids <= phi_min (NA: never)"), label = c("-", "the case's label"))
  d <- do.call(rbind, lapply(cols, function(cc) {
    base <- sub("^(blank_|org_minus_blank_)", "", cc)
    if (startsWith(cc, "ocean_gas_") || startsWith(cc, "pore_gas_"))
      return(data.frame(column = cc, unit = "mol per kg box rock", description = sprintf("%s gas phase: moles of %s", sub("_gas_.*", "", cc), sub(".*_gas_(.*)_mol$", "\\1", cc))))
    if (grepl("^(pore|ocean)_.+_mol$", cc) && sub("^(pore|ocean)_(.+)_mol$", "\2", cc) %in% KIT_PHOSPHATE_PHASES)
      return(data.frame(column = cc, unit = "mol per kg box rock", description = sprintf("%s held in the %s (phosphate phase)", sub("^(pore|ocean)_(.+)_mol$", "\2", cc),
                                                                                        if (startsWith(cc, "pore_")) "pores (REACT step)" else "ocean (ocean step)")))
    if (grepl("_core_mol$", cc))
      return(data.frame(column = cc, unit = "mol", description = paste(u[[sub("_core_mol$", "", cc)]][2], "- whole core (scale_to_core())")))
    v <- u[[base]]
    if (is.null(v)) v <- c("", "")
    pre <- if (startsWith(cc, "blank_")) "blank: " else if (startsWith(cc, "org_minus_blank_")) "organic minus blank: " else ""
    data.frame(column = cc, unit = if (base %in% c("ocean_pH", "ocean_pe") && pre != "") "-" else v[1], description = paste0(pre, v[2]))
  }))
  d
}

# RUNINFO.txt: what was run, with what, where, and how it ended.
write_runinfo <- function(run_dir, case, wall_s) {
  st <- function(v) { f <- file.path(run_dir, v, "state.rds"); if (file.exists(f)) readRDS(f) else NULL }
  o <- st("organic"); b <- st("blank")
  si <- Sys.info()
  ch <- utils::read.csv(file.path(run_dir, "organic", "chunks.csv"), stringsAsFactors = FALSE)
  lines <- c(
    sprintf("case: %s (%s)", case$name, case$case_file),
    if (nzchar(case$label)) sprintf("label: %s", case$label),
    if (nzchar(case$description)) sprintf("description: %s", case$description),
    sprintf("code: %s", code_version()),
    sprintf("PHREEQC: %s", phreeqc_version_string(ch)),
    sprintf("database: %s", sub("^# Database: ", "", grep("^# Database: ", readLines(file.path(run_dir, "organic", o$files[1]), n = 30), value = TRUE)[1])),
    sprintf("machine: %s (%s %s), %s, R %s", si[["nodename"]], si[["sysname"]], si[["release"]], si[["user"]], getRversion()),
    sprintf("written: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("wall time of this invocation: %.0f s; PHREEQC time over all chunks (organic): %.0f s", wall_s, sum(ch$seconds)),
    sprintf("ENDED BY: %s at t = %.6g yr (%d chunks, %d PHREEQC steps)", o$stopped, o$time_yr, o$chunk, sum(ch$steps)),
    sprintf("clog time: %s", if (is.na(o$clog_time)) "never clogged" else sprintf("%.6g yr", o$clog_time)),
    if (!is.null(b)) sprintf("blank: ended by %s at t = %.6g yr (%d chunks)", b$stopped, b$time_yr, b$chunk),
    "", "case (all fields, after inheritance and defaults):",
    vapply(names(case), function(nm) sprintf("  %s = %s", nm, paste(format(case[[nm]]), collapse = " ")), ""))
  writeLines(lines, file.path(run_dir, "RUNINFO.txt"))
}
