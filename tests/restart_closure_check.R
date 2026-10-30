# Restarts close: a run of 2 orbits and a 1-orbit run restarted for 1 more orbit agree on every
# element's system total to <= 1e-9 relative.
#
# tidal_kinetic + tidal_kinetic_restart, and bulk_kinetic + bulk_kinetic_restart (porewater mode), at
# 160 C (TEMP_C to change it) with 15 wt% IOM. The rows compared are the second orbit's four steps
# (tidal: REACT, MIX 4, MIX 1, MIX 2) or its last row (bulk). A row's system total of an element is
# system_totals() (phreeqc_tools.R): solution and its water, equilibrium phases, kinetic minerals and
# IOM channels. An element passes if it agrees to 1e-9 relative, or to within the solver's own
# mole-balance criterion: 1e-12 x the element's total in that orbit's REACT row (KNOBS
# -convergence_tolerance, which -high_precision sets to 1e-12; Parkhurst & Appelo 2013, p. 112).
# That covers trace elements pinned by a mineral (dissolved Ni at ~1e-16 molal against 0.24 mol of
# Ni metal), whose dissolved amount the solver fixes only to that tolerance.
#
# Before 2026-10-28 restarts carried only Al C Ca Fe Mg N Na P S Si (and Mtg) at 4-5 significant
# figures: K, Mn, Cl and Ni were lost, and every other element differed at ~1e-5.
#
# Runs PHREEQC (a few minutes). Run from the repo root:
#   Rscript tests/restart_closure_check.R
# TEMP_C=<C> sets the temperature; KEEP=<dir> keeps the runs in <dir>.

TEMP_C <- as.numeric(Sys.getenv("TEMP_C", "160"))
suppressPackageStartupMessages({ library(tidyverse); library(glue); library(R6) })
invisible(capture.output(source("ENTICES_LF_20260917.R")))

failures <- 0
check <- function(ok, what) {
  cat(if (isTRUE(ok)) "[PASS]" else "[FAIL]", what, "\n")
  if (!isTRUE(ok)) failures <<- failures + 1
}
quiet <- function(expr) { invisible(capture.output(r <- suppressWarnings(suppressMessages(expr)))); r }

sim <- Simulator$new(k = 1e-07, layer_thickness = 194200, config = EnceladusConfig$new())
quiet(sim$run_simulation())
ph <- quiet(PhreeqcIntegrator$new(sim))
invisible(quiet(ph$calculate_phreeqc_parameters(organic_wt_percent = 15)))
od <- normalizePath(if (nzchar(Sys.getenv("KEEP"))) Sys.getenv("KEEP") else tempfile("restart_closure_"), mustWork = FALSE)
dir.create(od, showWarnings = FALSE, recursive = TRUE)

run <- function(pqi) {
  db <- sub("^# Database: ", "", grep("^# Database: ", readLines(pqi, n = 20), value = TRUE)[1])
  r <- run_phreeqc(pqi, db, output = sub("\\.pqi$", ".out", pqi), workdir = dirname(pqi), warn_version = FALSE)
  if (r$errors > 0) stop(basename(pqi), ": ", r$errors, " PHREEQC errors")
  csv <- sub("\\.pqi$", ".csv", pqi)
  x <- utils::read.delim(csv, check.names = FALSE, strip.white = TRUE); names(x) <- trimws(names(x))
  list(x = x[trimws(x$state) == "react", ], csv = csv)
}
totals <- function(run, i, pq) system_totals(run$x, i, intersect(pq$totals, names(run$x)), pq$iom)
dissolved <- function(run, i, pq) {
  el <- setdiff(intersect(pq$totals, names(run$x)), c("H", "O"))
  stats::setNames(vapply(el, function(e) as.numeric(run$x[i, match(e, names(run$x))]) * as.numeric(run$x$mass_H2O[i]), 0), el)
}
# Worst element: relative difference, and whether it is within 1e-9 relative or 1e-12 x scale.
judge <- function(a, b, scale) {
  e <- union(names(a), names(b)); a <- a[e]; b <- b[e]; a[is.na(a)] <- 0; b[is.na(b)] <- 0; s <- scale[e]; s[is.na(s)] <- 0
  d <- abs(a - b); rel <- ifelse(pmax(abs(a), abs(b)) > 0, d / pmax(abs(a), abs(b)), 0)
  ok <- rel <= 1e-9 | d <= 1e-12 * s
  list(ok = all(ok), worst = e[which.max(rel)], rel = max(rel), fails = e[!ok])
}
compare <- function(A, ia, R, ir, pq, scale, label) {
  j <- judge(totals(A, ia, pq), totals(R, ir, pq), scale); d <- judge(dissolved(A, ia, pq), dissolved(R, ir, pq), scale)
  check(j$ok, sprintf("%s: element totals max rel. diff. %.1e (%s)%s; dissolved alone %.1e (%s)", label, j$rel, j$worst,
                      if (length(j$fails)) paste0(", beyond tolerance: ", paste(j$fails, collapse = " ")) else ", all within tolerance",
                      d$rel, d$worst))
}

# Tidal.
tk <- function(n, dir) quiet(ph$generate_phreeqc_input_tidal_kinetic(output_dir = file.path(od, dir), mode = "ungrouped",
                                                                     n_orbital_cycles = n, porewater_temp = TEMP_C, trace_C_molal = 1e-14))
A <- run(tk(2, "tidal_A")); P <- run(tk(1, "tidal_P"))
R <- run(quiet(ph$generate_phreeqc_input_tidal_kinetic_restart(output_dir = file.path(od, "tidal_R"), prev_csv = P$csv,
                                                              n_orbital_cycles = 1, porewater_temp = TEMP_C, depletion_mineral = "Troilite")))
pq <- pqi_initial_solids(sub("\\.csv$", ".pqi", A$csv))
steps <- c("REACT", "MIX 4", "MIX 1", "MIX 2")
scale <- totals(A, 5, pq)                          # orbit 2 REACT row: rock + porewater
for (j in 1:4) compare(A, 4 + j, R, j, pq, scale, sprintf("tidal orbit 2, %-6s", steps[j]))

# Refill (exchange-physics H, 2026-10-30): the same restart with refill_kg = +0.02 kg (ocean -> porewater) and
# -0.01 kg (porewater -> ocean). The whole system at the end of the orbit (porewater after MIX 1 + ocean after
# MIX 2 + the rock and pore gas of the REACT row + the ocean minerals and gas of MIX 4, less the pore gas moved
# to the ocean) is the same as without the refill, and the water has moved.
sol <- function(run, i) { s <- dissolved(run, i, pq); c(s, H = as.numeric(run$x$H_total_mol[i]), O = as.numeric(run$x$O_total_mol[i])) }
gas <- function(run, i) { el <- intersect(pq$totals, names(run$x)); g <- stats::setNames(rep(0, length(el)), el)
  for (k in grep("^g_", names(run$x), value = TRUE)) { f <- formula_elements(sub("\\(g\\)$", "", sub("^g_", "", k)))
    for (e in names(f)) g[e] <- g[e] + as.numeric(run$x[i, k]) * f[[e]] }
  g }
whole <- function(run) { a <- sol(run, 3); b <- sol(run, 4); e <- names(totals(run, 1, pq))
  s <- totals(run, 1, pq) - sol(run, 1)[e] + totals(run, 2, pq) - sol(run, 2)[e] - gas(run, 1)[e]
  s + a[e] + b[e] }
for (rk in c(0.02, -0.01)) {
  R1 <- run(quiet(ph$generate_phreeqc_input_tidal_kinetic_restart(output_dir = file.path(od, sprintf("tidal_refill_%g", rk)), prev_csv = P$csv,
                                                                 n_orbital_cycles = 1, porewater_temp = TEMP_C, depletion_mineral = "Troilite",
                                                                 refill_kg = rk)))
  j <- judge(whole(R), whole(R1), scale)
  dpw <- as.numeric(R1$x$mass_H2O[3]) - as.numeric(R$x$mass_H2O[3]); doc <- as.numeric(R1$x$mass_H2O[4]) - as.numeric(R$x$mass_H2O[4])
  check(j$ok && nrow(R1$x) == nrow(R$x) && abs(dpw / rk - 1) < 0.05 && abs(doc / -rk - 1) < 0.05,
        sprintf("refill %+g kg: whole-system totals max rel. diff. %.1e (%s)%s; porewater %+.4g kg, ocean %+.4g kg after the orbit; %d rows",
                rk, j$rel, j$worst, if (length(j$fails)) paste0(", beyond tolerance: ", paste(j$fails, collapse = " ")) else "",
                dpw, doc, nrow(R1$x)))
}

# Bulk (porewater mode).
bk <- function(n, dir) quiet(ph$generate_phreeqc_input_bulk_kinetic(output_dir = file.path(od, dir), n_orbital_cycles = n, mode = "porewater",
                                                                    water_rock_ratio = 1, porewater_temp = TEMP_C, trace_C_molal = 1e-14))
A <- run(bk(2, "bulk_A")); P <- run(bk(1, "bulk_P"))
R <- run(quiet(ph$generate_phreeqc_input_bulk_kinetic_restart(output_dir = file.path(od, "bulk_R"), prev_csv = P$csv,
                                                             n_orbital_cycles = 1, mode = "porewater", porewater_temp = TEMP_C)))
pq <- pqi_initial_solids(sub("\\.csv$", ".pqi", A$csv))
compare(A, nrow(A$x), R, nrow(R$x), pq, totals(A, nrow(A$x), pq), "bulk orbit 2, last row ")

if (nzchar(Sys.getenv("KEEP"))) cat("kept:", od, "\n") else unlink(od, recursive = TRUE)
cat(sprintf("\nSUMMARY: %d failure(s).\n", failures))
if (failures > 0) quit(status = 1)
