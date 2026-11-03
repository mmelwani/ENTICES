# Step-limit validation (production kit, PI 2026-11-02; docs/MODEL.md, Step limits). From the repo root:
#   Rscript scripts/step_limit_validation.R <repo> <out> <case: ref|k11> <test: x|phi> <value> <t_end_yr>
# x test: porosity held fixed, x_step_max = value; phi test: porosity evolving, phi_tol = value, x_step_max 0.1.
# ref: Kozeny-Carman d = 1.59 cm; k11: d = 1.59e-4 m (k ~ 1e-11). Prints the ocean C(+4), pH, H2, Na at the end.
a <- commandArgs(TRUE); setwd(a[1]); out <- a[2]; cs <- a[3]; test <- a[4]; v <- as.numeric(a[5]); te <- as.numeric(a[6])
suppressPackageStartupMessages({ library(tidyverse); library(glue); library(R6) })
invisible(capture.output(source("ENTICES.R")))
d <- if (cs == "ref") 0.0159 else 1.59e-4
tmpl <- Simulator$new(k = "kozeny_carman", grain_diameter = d, layer_thickness = 194200, config = EnceladusConfig$new(), phi0 = 0.32)
t0 <- Sys.time()
r <- suppressWarnings(if (test == "x")
  run_tidal_evolution(tmpl, organic_wt_percent = 15, t_end_years = te, output_dir = out, porosity_evolution = FALSE,
                      x_step_max = v, porewater_temp = "profile", trace_C_molal = 1e-14)
else
  run_tidal_evolution(tmpl, organic_wt_percent = 15, t_end_years = te, output_dir = out, phi_tol = v, x_step_max = 0.1,
                      porewater_temp = "profile", trace_C_molal = 1e-14))
r$wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
saveRDS(r, file.path(out, "val.rds"))
ch <- r$chunks; n <- nrow(ch)
cat(sprintf("%s %s %g: %s, %d chunks, %d steps, %.0f s, t %.4g yr, G max %g, x_step max %.3g, phi %.4g; ocean C4 %.6g pH %.6f H2 %.6g Na %.6g\n",
            cs, test, v, r$stopped, n, sum(ch$steps), r$wall, ch$time_end_yr[n], max(ch$G), max(ch$x_step), ch$phi_end_solids[n],
            ch$ocean_C4[n], ch$ocean_pH[n], ch$ocean_H2[n], ch$ocean_Na[n]))
