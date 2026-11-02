# Grouped time steps against one orbit per step (production kit, item 1, PI 2026-11-01): the reference case
# (cases/reference.R's physics: Kozeny-Carman k with d = 1.59 cm, phi0 0.32, the skin at its profile
# temperature, 15 wt% IOM, to_ocean, trace C 1e-14) run for N_ORBITS orbits by run_tidal_evolution() twice,
# ungrouped (chunks of up to 64 single-orbit steps) and grouped (4 steps of G orbits per chunk, G doubling,
# x_step <= 0.1, porosity change within phi_tol). The ocean's C(+4), pH, H2 and Na at the end agree within 1%.
#
# Long: ~1-2 h at N_ORBITS = 1000 (the ungrouped run is most of it). Run from the repo root:
#   Rscript tests/grouped_validation_check.R          (N_ORBITS=<n> to change the length; KEEP=<dir> keeps the runs)
N_ORBITS <- as.numeric(Sys.getenv("N_ORBITS", "1000"))
suppressPackageStartupMessages({ library(tidyverse); library(glue); library(R6) })
invisible(capture.output(source("ENTICES_LF_20260917.R")))

failures <- 0
check <- function(ok, what) {
  cat(if (isTRUE(ok)) "[PASS]" else "[FAIL]", what, "\n")
  if (!isTRUE(ok)) failures <<- failures + 1
}
od <- normalizePath(if (nzchar(Sys.getenv("KEEP"))) Sys.getenv("KEEP") else tempfile("grouped_validation_"), mustWork = FALSE)
tmpl <- Simulator$new(k = "kozeny_carman", grain_diameter = 0.0159, layer_thickness = 194200, config = EnceladusConfig$new(), phi0 = 0.32)
run <- function(grouped, dir) {
  t0 <- Sys.time()
  r <- suppressWarnings(run_tidal_evolution(tmpl, organic_wt_percent = 15, n_orbits = N_ORBITS, output_dir = file.path(od, dir),
                                            grouped = grouped, porewater_temp = "profile", trace_C_molal = 1e-14))
  r$wall <- as.numeric(difftime(Sys.time(), t0, units = "secs")); r
}
U <- run(FALSE, "ungrouped"); G <- run(TRUE, "grouped")
u <- U$chunks[nrow(U$chunks), ]; g <- G$chunks[nrow(G$chunks), ]
rel <- c(C4 = abs(g$ocean_C4 / u$ocean_C4 - 1), pH = abs(g$ocean_pH / u$ocean_pH - 1), H2 = abs(g$ocean_H2 / u$ocean_H2 - 1), Na = abs(g$ocean_Na / u$ocean_Na - 1))
cat(sprintf("ungrouped: %d chunks, %d steps, %.0f s; grouped: %d chunks, %d steps, %.0f s, G up to %g (x_step %.3g)\n",
            nrow(U$chunks), sum(U$chunks$steps), U$wall, nrow(G$chunks), sum(G$chunks$steps), G$wall, max(G$chunks$G), max(G$chunks$x_step)))
check(u$orbit_to == N_ORBITS && g$orbit_to == N_ORBITS && all(rel < 0.01),
      sprintf("orbit %g, ocean grouped / ungrouped: C(+4) %.4g / %.4g, pH %.4f / %.4f, H2 %.4g / %.4g, Na %.4g / %.4g (worst %.2g%%, %s)",
              N_ORBITS, g$ocean_C4, u$ocean_C4, g$ocean_pH, u$ocean_pH, g$ocean_H2, u$ocean_H2, g$ocean_Na, u$ocean_Na,
              100 * max(rel), names(rel)[which.max(rel)]))
if (nzchar(Sys.getenv("KEEP"))) { saveRDS(list(U = U, G = G), file.path(od, "validation.rds")); cat("kept:", od, "\n") } else unlink(od, recursive = TRUE)
cat(sprintf("\nSUMMARY: %d failure(s).\n", failures))
if (failures > 0) quit(status = 1)
