# Gas phases conserve mass, and the two pore_gas_fate options differ only through the bubble.
#
#  1. Mass balance. tidal_kinetic, 160 C, 15 wt% IOM, 2 orbits, for pore_gas_fate "retain" and
#     "to_ocean": every element's whole-system total at the end of each orbit equals its t = 0 total
#     to <= 1e-9 relative, or within the run's mole-balance criterion, 1e-8 x the element's whole-system
#     total (PI 2026-10-31). The end-of-orbit system is the porewater (MIX 1 row) and ocean (MIX 2 row)
#     solutions, plus the solids, IOM and pore gas of the REACT row and the seafloor solids and ocean
#     gas of the MIX 4 row; under to_ocean the pore bubble is counted once, in the ocean. t = 0 is the
#     two initial solutions plus the .pqi's kinetic m0 and equilibrium-phase amounts.
#  2. No-bubble identity. The same at 25 C with the Fe and Ni metal set to 0 (no H2 pulse), where no
#     bubble forms: "retain" and "to_ocean" give the same output. Every element's system total in
#     every row agrees to 1e-9 relative or within the run's mole-balance criterion (1e-8 x the
#     element's whole-system total), and every numeric column to 1e-6 relative or 1e-12 absolute
#     (molal; PI 2026-10-29: with the smaller exchange of exchange-physics the ocean's H2 is 2.6e-8
#     molal and differs by 4e-14, 1.6e-6 relative). The ocean step of
#     "to_ocean" integrates the (zero) transfer over 1 s, so PHREEQC solves the same equilibrium from
#     a different start, and species differ at its 1e-8 convergence criterion (measured 2026-10-28:
#     5e-8 relative in dissolved H2, only in the ocean step).
#
# Runs PHREEQC (about 15 minutes). Run from the repo root:  Rscript tests/gas_mass_balance_check.R
# KEEP=<dir> keeps the runs.

suppressPackageStartupMessages({ library(tidyverse); library(glue); library(R6) })
invisible(capture.output(source("ENTICES_LF_20260917.R")))

failures <- 0
check <- function(ok, what) {
  cat(if (isTRUE(ok)) "[PASS]" else "[FAIL]", what, "\n")
  if (!isTRUE(ok)) failures <<- failures + 1
}
quiet <- function(expr) { invisible(capture.output(r <- suppressWarnings(suppressMessages(expr)))); r }
od <- normalizePath(if (nzchar(Sys.getenv("KEEP"))) Sys.getenv("KEEP") else tempfile("gas_mass_balance_"), mustWork = FALSE)
dir.create(od, showWarnings = FALSE, recursive = TRUE)

sim <- Simulator$new(k = 1e-07, layer_thickness = 194200, config = EnceladusConfig$new())
quiet(sim$run_simulation())
ph <- quiet(PhreeqcIntegrator$new(sim))
invisible(quiet(ph$calculate_phreeqc_parameters(organic_wt_percent = 15)))

run <- function(pqi) {
  db <- sub("^# Database: ", "", grep("^# Database: ", readLines(pqi, n = 20), value = TRUE)[1])
  r <- run_phreeqc(pqi, db, output = sub("\\.pqi$", ".out", pqi), workdir = dirname(pqi), warn_version = FALSE)
  if (r$errors > 0) stop(basename(pqi), ": ", r$errors, " PHREEQC errors")
  csv <- sub("\\.pqi$", ".csv", pqi)
  x <- utils::read.delim(csv, check.names = FALSE, strip.white = TRUE); names(x) <- trimws(names(x)); x$state <- trimws(x$state)
  list(x = x, pq = pqi_initial_solids(pqi), seconds = r$seconds)
}
gen <- function(dir, T, fate) quiet(ph$generate_phreeqc_input_tidal_kinetic(output_dir = file.path(od, dir), mode = "ungrouped",
                                    n_orbital_cycles = 2, porewater_temp = T, pore_gas_fate = fate, trace_C_molal = 1e-14))

# Whole-system element totals: t = 0, and the end of each orbit.
orbit_totals <- function(R, fate) {
  x <- R$x; pq <- R$pq; el <- intersect(pq$totals, names(x)); v <- utils::read.csv(MINERAL_VOLUMES_CSV, stringsAsFactors = FALSE)
  ini <- x[x$state == "i_soln", ]; r <- x[x$state == "react", ]
  tot <- function(d, i) system_totals(d, i, el, pq$iom, v)
  sol <- function(d, i) {                                      # the solution alone: dissolved + water
    s <- stats::setNames(vapply(el, function(e) as.numeric(d[i, match(e, names(d))]) * as.numeric(d$mass_H2O[i]), 0), el)
    s[c("H", "O")] <- c(as.numeric(d$H_total_mol[i]), as.numeric(d$O_total_mol[i])); s }
  gas <- function(d, i) { g <- stats::setNames(rep(0, length(el)), el)
    for (k in grep("^g_", names(d), value = TRUE)) { f <- formula_elements(sub("\\(g\\)$", "", sub("^g_", "", k)))
      for (e in names(f)) g[e] <- g[e] + as.numeric(d[i, k]) * f[[e]] }
    g }
  t0 <- sol(ini, 1) + sol(ini, 2)
  for (p in names(pq$minerals)) { f <- formula_elements(v$formula[match(p, v$phase)])
    for (e in names(f)) t0[e] <- t0[e] + pq$minerals[[p]] * f[[e]] }
  for (j in seq_len(nrow(pq$iom))) { f <- formula_elements(pq$iom$formula[j]); for (e in names(f)) t0[e] <- t0[e] + pq$iom$m0[j] * f[[e]] }
  out <- list(t0 = t0)
  for (o in seq_len(nrow(r) %/% 4)) {
    i <- 4 * (o - 1)
    s <- sol(r, i + 3) + sol(r, i + 4) + (tot(r, i + 1) - sol(r, i + 1)) + (tot(r, i + 2) - sol(r, i + 2))
    if (fate == "to_ocean") s <- s - gas(r, i + 1)              # moved to the ocean, already in MIX 4's gas and solution
    out[[sprintf("orbit %d", o)]] <- s
  }
  do.call(rbind, out)
}

# 1. Mass balance.
runs <- list()
for (fate in c("retain", "to_ocean")) {
  R <- runs[[fate]] <- run(gen(paste0("mb_", fate), 160, fate))
  m <- orbit_totals(R, fate)
  for (o in 2:nrow(m)) {
    rel <- abs(m[o, ] - m[1, ]) / pmax(abs(m[1, ]), 1e-300); rel[m[1, ] == 0 & m[o, ] == 0] <- 0
    key <- c("H", "O", "C", "Mtg", "S", "N", "Fe")
    # within 1e-9 relative, or within the run's mole-balance criterion, 1e-8 x the element's whole-system total
    # (PI 2026-10-31: with phi_c = 0.4 the first REACT step converges Al only to 1.0e-9 mol of 0.39, inside it)
    ok_el <- rel <= 1e-9 | abs(m[o, ] - m[1, ]) <= 1e-8 * abs(m[1, ])
    check(all(ok_el), sprintf("160 C %-8s %s: worst element %s %.1e; H %.1e O %.1e C %.1e Mtg %.1e S %.1e N %.1e Fe %.1e",
                                    fate, rownames(m)[o], names(rel)[which.max(rel)], max(rel), rel["H"], rel["O"], rel["C"],
                                    rel["Mtg"], rel["S"], rel["N"], rel["Fe"]))
  }
  rr <- R$x[R$x$state == "react", ]
  cat(sprintf("   %s: %.0f s; orbit 1 pore bubble %.4g L (Vgas/Vpore %.3g), ocean gas %.4g L\n", fate, R$seconds,
              rr$Vgas_L[1], rr$Vgas_Vpore[1], rr$Vgas_L[2]))
}

# 2. No-bubble identity: 25 C, Fe and Ni metal removed.
nometal <- function(pqi) {
  l <- readLines(pqi)
  l <- sub("^(   (Fe|Ni) +0 +)[0-9.eE+-]+( +dissolve_only)", "\\10\\3", l)
  writeLines(l, pqi); pqi
}
A <- run(nometal(gen("nb_retain", 25, "retain"))); B <- run(nometal(gen("nb_to_ocean", 25, "to_ocean")))
num <- names(A$x)[vapply(A$x, is.numeric, logical(1))]
# a value passes if within 1e-6 relative or 1e-12 absolute (molal); d is each column's largest relative
# difference among the values that are not within 1e-12 absolute
d <- vapply(num, function(c) { a <- A$x[[c]]; b <- B$x[[c]]; k <- abs(a - b) > 1e-12
  if (!any(k)) 0 else max(abs(a - b)[k] / pmax(abs(a), abs(b))[k]) }, 0)
d_all <- vapply(num, function(c) { a <- A$x[[c]]; b <- B$x[[c]]; k <- abs(a) > 1e-20 | abs(b) > 1e-20
  if (!any(k)) 0 else max(abs(a - b)[k] / pmax(abs(a), abs(b))[k]) }, 0)
el <- intersect(A$pq$totals, names(A$x)); ra <- which(A$x$state == "react")
# an element passes if within 1e-9 relative, or within 1e-8 (the run's mole-balance convergence
# criterion, KNOBS -convergence_tolerance) x its whole-system total in the first REACT row
scale <- system_totals(A$x, ra[1], el, A$pq$iom)
worst <- c(rel = 0); bad <- character(0)
for (i in ra) { a <- system_totals(A$x, i, el, A$pq$iom); b <- system_totals(B$x, i, el, B$pq$iom)
  dd <- abs(a - b); rel <- ifelse(pmax(abs(a), abs(b)) > 0, dd / pmax(abs(a), abs(b)), 0)
  bad <- union(bad, names(rel)[!(rel <= 1e-9 | dd <= 1e-8 * scale[names(rel)])])
  if (max(rel) > worst[["rel"]]) worst <- c(rel = max(rel), row = i, el = names(rel)[which.max(rel)]) }
te <- length(bad) == 0
vg <- max(c(A$x$Vgas_L, B$x$Vgas_L))
check(vg == 0 && te && max(d) <= 1e-6,
      sprintf("25 C, no metal: largest gas volume %g L; retain vs to_ocean: element totals within tolerance%s (largest relative %s, %s, row %s); %d columns within 1e-6 relative or 1e-12 absolute (largest relative beyond 1e-12 absolute %.1e, %s; largest relative overall %.1e, %s)",
              vg, if (te) "" else paste0(" EXCEPT ", paste(bad, collapse = " ")), signif(as.numeric(worst[["rel"]]), 2), worst["el"], worst["row"],
              length(d), max(d), names(d)[which.max(d)], max(d_all), names(d_all)[which.max(d_all)]))

if (nzchar(Sys.getenv("KEEP"))) cat("kept:", od, "\n") else unlink(od, recursive = TRUE)
cat(sprintf("\nSUMMARY: %d failure(s).\n", failures))
if (failures > 0) quit(status = 1)
