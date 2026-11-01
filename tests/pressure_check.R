# Porewater and ocean pressures from the PlanetProfile profile (Stage 1, 2026-10-27).
#
#  1. profile_conditions() on profiles/EnceladusProfile_...: porewater 261.2 atm (mean over the
#     194.2 km layer), ocean 66.7 atm (ocean floor).
#  2. Every generator writes the porewater pressure into the porewater SOLUTION (the ocean-floor
#     pressure into a bulk "hydrosphere" solution and the ocean SOLUTION), sets REACTION_PRESSURE
#     1/2 in the tidal generators, and stamps both pressures in the header.
#  3. A missing profile is an error in every generator (no fallback pressure).
#  4. Restarts compare their pressures with the parent's header: a match passes, a mismatch is an
#     error (a warning with allow_pressure_change = TRUE), an unstamped parent is a warning.
#  5. KNOBS -logfile is false unless logfile = TRUE, in every generator (2026-10-29).
#
# Generates inputs only; nothing is run in PHREEQC.
# Run from the repo root:  Rscript tests/pressure_check.R

suppressPackageStartupMessages({ library(tidyverse); library(glue); library(R6) })
invisible(capture.output(source("ENTICES_LF_20260917.R")))

failures <- 0
check <- function(ok, what) {
  cat(if (isTRUE(ok)) "[PASS]" else "[FAIL]", what, "\n")
  if (!isTRUE(ok)) failures <<- failures + 1
}
quiet <- function(expr) { invisible(capture.output(r <- suppressMessages(expr))); r }

# reactive_depth = "layer": these checks are of the profile mean over the whole 194.2 km layer (section 6
# checks the default skin box).
sim <- Simulator$new(k = 1e-07, layer_thickness = 194200, config = EnceladusConfig$new(), reactive_depth = "layer")
quiet(sim$run_simulation())
ph <- quiet(PhreeqcIntegrator$new(sim))
invisible(quiet(ph$calculate_phreeqc_parameters(organic_wt_percent = 0)))
profile <- "profiles/EnceladusProfile_Seawater_10.0ppt_Tb272.4578K_PorousRock.txt"
check(identical(ph$profile_file, profile) && file.exists(profile), "default profile_file is the copy in profiles/")

# 1. Profile pressures.
pc <- profile_conditions(profile, 194200)
check(sprintf("%.1f", pc$pore_pressure_atm) == "261.2" && sprintf("%.1f", pc$ocean_pressure_atm) == "66.7",
      sprintf("profile: porewater %.2f atm, ocean %.2f atm", pc$pore_pressure_atm, pc$ocean_pressure_atm))

# 2. What each generator writes.
od <- tempfile("pressure_check_")
pressure_lines <- function(f) {
  l <- readLines(f)
  sol <- function(n) {
    i <- grep(sprintf("^SOLUTION %d ", n), l)
    if (!length(i)) return(NA_character_)
    p <- grep("^\\s+pressure\\s", l[i:min(length(l), i + 12)], value = TRUE)[1]
    sub("^\\s+pressure\\s+(\\S+).*", "\\1", p)
  }
  rp <- function(n) { i <- grep(sprintf("^REACTION_PRESSURE %d$", n), l); if (length(i)) trimws(l[i + 1]) else NA_character_ }
  list(s1 = sol(1), s2 = sol(2), rp1 = rp(1), rp2 = rp(2),
       stamp = all(c("# Porewater pressure: 261.2 atm (profile mean over the reacting box)",
                     "# Ocean pressure: 66.7 atm (profile, ocean floor)") %in% l[1:10]),
       use_rp = sum(l == "USE reaction_pressure 1"), use_rp2 = sum(l == "USE reaction_pressure 2"))
}
gens <- list(
  tidal_kinetic   = list(fn = "generate_phreeqc_input_tidal_kinetic", args = list(mode = "ungrouped", n_orbital_cycles = 2, porewater_temp = 100),
                         s1 = "261.2", s2 = "66.7", tidal = TRUE),
  tidal_equil     = list(fn = "generate_phreeqc_input_tidal_equil", args = list(n_orbital_cycles = 2, porewater_temp = 100),
                         s1 = "261.2", s2 = "66.7", tidal = TRUE),
  bulk_kinetic_pw = list(fn = "generate_phreeqc_input_bulk_kinetic", args = list(n_orbital_cycles = 2, mode = "porewater", water_rock_ratio = 1, porewater_temp = 100),
                         s1 = "261.2", s2 = NA_character_, tidal = FALSE),
  bulk_kinetic_hy = list(fn = "generate_phreeqc_input_bulk_kinetic", args = list(n_orbital_cycles = 2, mode = "hydrosphere", porewater_temp = 100),
                         s1 = "66.7", s2 = NA_character_, tidal = FALSE),
  bulk_equil_pw   = list(fn = "generate_phreeqc_input_bulk_equil", args = list(mode = "porewater", water_rock_ratio = 1, water_temp = 100),
                         s1 = "261.2", s2 = NA_character_, tidal = FALSE),
  bulk_equil_hy   = list(fn = "generate_phreeqc_input_bulk_equil", args = list(mode = "hydrosphere", water_temp = 100),
                         s1 = "66.7", s2 = NA_character_, tidal = FALSE))
files <- list()
for (nm in names(gens)) {
  g <- gens[[nm]]
  f <- quiet(do.call(ph[[g$fn]], c(list(output_dir = file.path(od, nm)), g$args)))
  files[[nm]] <- f
  p <- pressure_lines(f)
  ok <- identical(p$s1, g$s1) && identical(p$s2, g$s2) && p$stamp &&
    if (g$tidal) identical(p$rp1, "261.2") && identical(p$rp2, "66.7") && p$use_rp == 4 && p$use_rp2 == 2
    else is.na(p$rp1) && is.na(p$rp2)
  check(ok, sprintf("%-15s SOLUTION 1 %s, SOLUTION 2 %s, REACTION_PRESSURE %s/%s, header stamped: %s",
                    nm, p$s1, p$s2, p$rp1, p$rp2, p$stamp))
}

# 3. No profile, no input.
ph$profile_file <- file.path(od, "no_such_profile.txt")
for (nm in names(gens)) {
  g <- gens[[nm]]
  e <- tryCatch({ quiet(do.call(ph[[g$fn]], c(list(output_dir = file.path(od, "missing")), g$args))); "no error" },
                error = function(e) conditionMessage(e))
  check(grepl("profile_file not found", e, fixed = TRUE), sprintf("%-15s stops without a profile", nm))
}
ph$profile_file <- profile

# 4. Restarts against a stub parent (a .pqi header and the CSV rows a restart reads).
stub_parent <- function(dir, header, tidal) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  csv <- file.path(dir, "parent.csv")
  rows <- if (tidal) data.frame(state = "react", soln = c(1, 4, 1, 2), step = 1, pH = 9, pe = -8, mass_H2O = c(0.157, 0.497, 0.157, 0.340), Time_Years = 0.004)
          else data.frame(state = "react", soln = 1, step = 1, pH = 9, pe = -8, mass_H2O = 1, Time_Years = 0.004)
  for (el in ENTICES_ELEMENTS) rows[[el]] <- 1e-6
  write.table(rows, csv, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(c("TITLE stub parent", "# Target PHREEQC version: 3.8.6", "# Database: databases/coreclath_C4only.dat", header), sub("\\.csv$", ".pqi", csv))
  csv
}
restart <- function(tidal, csv, ...) {
  fn <- if (tidal) "generate_phreeqc_input_tidal_kinetic_restart" else "generate_phreeqc_input_bulk_kinetic_restart"
  args <- c(list(output_dir = file.path(dirname(csv), "out"), prev_csv = csv, n_orbital_cycles = 1, porewater_temp = 100), list(...))
  tryCatch(withCallingHandlers({ f <- quiet(do.call(ph[[fn]], args)); list(file = f, warn = NULL, err = NULL) },
                               warning = function(w) { if (grepl("pressure", conditionMessage(w))) assign("w_msg", conditionMessage(w), envir = globalenv()); invokeRestart("muffleWarning") }),
           error = function(e) list(file = NULL, err = conditionMessage(e)))
}
for (tidal in c(TRUE, FALSE)) {
  kind <- if (tidal) "tidal_kinetic_restart" else "bulk_kinetic_restart"
  good <- pressure_header_lines(261.2, 66.7)
  r <- restart(tidal, stub_parent(file.path(od, kind, "match"), good, tidal))
  p <- if (!is.null(r$file)) pressure_lines(r$file) else list()
  check(!is.null(r$file) && "# parent pressures: porewater 261.2 atm, ocean 66.7 atm" %in% readLines(r$file) &&
          identical(p$s1, "261.2") && isTRUE(p$stamp) &&
          if (tidal) identical(p$s2, "66.7") && identical(p$rp1, "261.2") && identical(p$rp2, "66.7") else is.na(p$s2),
        sprintf("%s: parent at the same pressures passes; SOLUTION 1 %s, SOLUTION 2 %s, REACTION_PRESSURE %s/%s",
                kind, p$s1, p$s2, p$rp1, p$rp2))
  bad <- stub_parent(file.path(od, kind, "mismatch"), pressure_header_lines(66.7, 66.7), tidal)
  r <- restart(tidal, bad)
  check(is.null(r$file) && grepl("parent ran at porewater 66.7", r$err, fixed = TRUE), sprintf("%s: a pressure mismatch is an error", kind))
  w_msg <- NULL
  r <- restart(tidal, bad, allow_pressure_change = TRUE)
  check(!is.null(r$file) && grepl("allow_pressure_change = TRUE", w_msg, fixed = TRUE) &&
          any(grepl("this restart switches to 261.2 / 66.7 atm", readLines(r$file), fixed = TRUE)),
        sprintf("%s: allow_pressure_change = TRUE turns the mismatch into a warning", kind))
  w_msg <- NULL
  r <- restart(tidal, stub_parent(file.path(od, kind, "unstamped"), character(0), tidal))
  check(!is.null(r$file) && grepl("not stamped", w_msg, fixed = TRUE) &&
          "# parent pressures: unknown (not checked)" %in% readLines(r$file),
        sprintf("%s: an unstamped parent is a warning", kind))
}

# 5. KNOBS -logfile: false by default in every generator (restarts included), true with logfile = TRUE.
logfile_line <- function(f) grep("^   -logfile ", readLines(f), value = TRUE)
for (nm in names(gens)) {
  g <- gens[[nm]]
  on <- quiet(do.call(ph[[g$fn]], c(list(output_dir = file.path(od, "logfile", nm), logfile = TRUE), g$args)))
  check(identical(logfile_line(files[[nm]]), "   -logfile false") && identical(logfile_line(on), "   -logfile true"),
        sprintf("%-15s -logfile false by default, true with logfile = TRUE", nm))
}
for (tidal in c(TRUE, FALSE)) {
  kind <- if (tidal) "tidal_kinetic_restart" else "bulk_kinetic_restart"
  csv <- stub_parent(file.path(od, "logfile", kind), pressure_header_lines(261.2, 66.7), tidal)
  off <- restart(tidal, csv)$file
  on <- restart(tidal, csv, logfile = TRUE, suffix = "_log")$file
  check(identical(logfile_line(off), "   -logfile false") && identical(logfile_line(on), "   -logfile true"),
        sprintf("%-15s -logfile false by default, true with logfile = TRUE", kind))
}

# 6. The default skin box (reactive_depth = "skin"): at k = 1e-7 the box is the top skin depth delta (73.8 km
#    with alpha = phi/phi_c), and the porewater pressure is the profile mean over it; the
#    ocean pressure is unchanged.
sk <- Simulator$new(k = 1e-07, layer_thickness = 194200, config = EnceladusConfig$new()); quiet(sk$run_simulation())
phs <- quiet(PhreeqcIntegrator$new(sk)); invisible(quiet(phs$calculate_phreeqc_parameters(organic_wt_percent = 15)))
fs <- quiet(phs$generate_phreeqc_input_tidal_kinetic(output_dir = file.path(od, "skin"), mode = "ungrouped", n_orbital_cycles = 1, porewater_temp = 100))
ps <- pressure_lines(fs)
p_skin <- sprintf("%.1f", quiet(profile_conditions(profile, sk$skin_depth()$delta))$pore_pressure_atm)
check(sk$reactive_thickness() == sk$skin_depth()$delta && sk$reactive_thickness() < 194200 && identical(ps$s1, p_skin) &&
        identical(ps$s2, "66.7") && any(grepl(sprintf("# Reacting box: top %.4g m of the 1.942e+05 m layer", sk$reactive_thickness()), readLines(fs), fixed = TRUE)),
      sprintf("skin box: top %.4g m (= delta); SOLUTION 1 %s atm (profile mean over it: %s), SOLUTION 2 %s atm", sk$reactive_thickness(), ps$s1, p_skin, ps$s2))

unlink(od, recursive = TRUE)
cat(sprintf("\nSUMMARY: %d failure(s).\n", failures))
if (failures > 0) quit(status = 1)
