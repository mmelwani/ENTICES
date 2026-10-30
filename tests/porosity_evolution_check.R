# run_tidal_evolution() (porosity_evolution.R; exchange-physics E, 2026-10-29): the bookkeeping, on a
# short 25 C run (fast; the chemistry is the production case's at 25 C).
#
#  1. porosity_evolution = TRUE: each chunk starts at the porosity the previous one ended with (capped at
#     phi_min); k = k_of_phi(phi) (Kozeny-Carman here); alpha = phi/phi_c; V_total is fixed (rigid frame);
#     a chunk over phi_tol is only accepted at 1 orbit; the orbits add up.
#  2. porosity_evolution = FALSE: phi stays phi0 and k stays k0 in every chunk.
#  3. phi = 1 - V_solid / V_total, recomputed here from the chunk's bulk_density() (with origin_pqi).
#
# Runs PHREEQC (a few minutes). Run from the repo root:  Rscript tests/porosity_evolution_check.R

suppressPackageStartupMessages({ library(tidyverse); library(glue); library(R6) })
invisible(capture.output(source("ENTICES_LF_20260917.R")))

failures <- 0
check <- function(ok, what) {
  cat(if (isTRUE(ok)) "[PASS]" else "[FAIL]", what, "\n")
  if (!isTRUE(ok)) failures <<- failures + 1
}
od <- tempfile("porosity_evolution_")
kc <- function(phi) (0.01594)^2 * phi^3 / (180 * (1 - phi)^2)        # Kozeny-Carman, k(0.32) = 1e-7
tmpl <- Simulator$new(k = 1e-7, layer_thickness = 194200, config = EnceladusConfig$new(), phi0 = 0.32)

on <- run_tidal_evolution(tmpl, organic_wt_percent = 15, n_orbits = 4, output_dir = file.path(od, "on"),
                          porosity_evolution = TRUE, k_of_phi = kc, porewater_temp = 25)$chunks
n <- nrow(on)
chain <- n == 1 || all(abs(on$phi_start[-1] - pmax(0.01, on$phi_end_solids[-n])) < 1e-12)
check(on$phi_start[1] == 0.32 && chain && all(abs(on$k / kc(on$phi_start) - 1) < 1e-12) &&
        all(abs(on$alpha - pmin(1, on$phi_start / tmpl$params$phi_c)) < 1e-12) && length(unique(on$V_total_cm3)) == 1 &&
        all(!on$over_tol | on$n_orbits == 1) && sum(on$n_orbits) == 4 && max(on$orbit_to) == 4,
      sprintf("evolution on: %d chunks (%s orbits); phi %s; k %s",
              n, paste(on$n_orbits, collapse = "+"), paste(sprintf("%.4f", c(on$phi_start, on$phi_end_solids[n])), collapse = " -> "),
              paste(sprintf("%.3g", on$k), collapse = ", ")))

off <- run_tidal_evolution(tmpl, organic_wt_percent = 15, n_orbits = 2, output_dir = file.path(od, "off"),
                           porosity_evolution = FALSE, k_of_phi = kc, porewater_temp = 25)$chunks
check(all(off$phi_start == 0.32) && all(off$k == kc(0.32)) && sum(off$n_orbits) == 2,
      sprintf("evolution off: phi %s, k %s in %d chunks (solids imply phi %.4f at the end)",
              paste(unique(off$phi_start), collapse = ","), paste(unique(sprintf("%.3g", off$k)), collapse = ","), nrow(off),
              off$phi_end_solids[nrow(off)]))

# 3. phi from the last chunk's solids, recomputed.
f <- list.files(file.path(od, "on"), "\\.pqi$", recursive = TRUE, full.names = TRUE)
first <- f[grepl("chunk001", f)]; last <- f[grepl(sprintf("chunk%03d", n), f)]
bd <- suppressWarnings(bulk_density(sub("\\.pqi$", ".csv", last), porosity = on$phi_start[n], origin_pqi = first))$summary
phi_re <- 1 - bd$solid_volume_cm3[nrow(bd)] / on$V_total_cm3[n]
check(abs(phi_re - on$phi_end_solids[n]) < 1e-12, sprintf("phi from the solids: %.6f (recomputed %.6f)", on$phi_end_solids[n], phi_re))

unlink(od, recursive = TRUE)
cat(sprintf("\nSUMMARY: %d failure(s).\n", failures))
if (failures > 0) quit(status = 1)
