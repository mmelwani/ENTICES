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

# Refill (H): before each chunk after the first, the porewater is brought to m_pw(0) phi / phi0 by moving ocean
# water in; so each chunk's refill equals the previous chunk's target less its porewater, and after the refill
# the porewater (before the next chunk's orbits) is on target: checked by the next chunk's first REACT row
# water being near the target (hydration in that REACT step moves it a little).
check(on$refill_in_kg[1] == 0 && n > 1 && all(abs(on$refill_in_kg[-1] - (on$water_target_kg[-n] - on$water_phreeqc_kg[-n])) < 1e-12) &&
        all(on$refill_in_kg[-1] > 0),
      sprintf("refill: %s kg before chunks 2..%d (exchange per orbit %s kg)", paste(sprintf("%.3g", on$refill_in_kg[-1]), collapse = ", "), n,
              paste(sprintf("%.2g", on$exchange_kg_per_orbit), collapse = ", ")))

# Drift and the ocean per kg of rock (PI 2026-10-31): the cumulative whole-system drift stays within 1e-6
# relative; each restart scales the ocean by old box rock / new box rock (ocean_scale), and the ocean per kg
# of rock at a chunk's start is the profile's ocean / the box's rock (before that chunk's refill, below).
os_ok <- n == 1 || all(abs(on$ocean_scale[-1] - on$layer_mass_kg[-n] / on$layer_mass_kg[-1]) < 1e-12)
check(max(on$drift_max) < 1e-6 && on$ocean_scale[1] == 1 && os_ok,
      sprintf("drift: max %.1e (%s); ocean_scale %s", max(on$drift_max), on$drift_element[which.max(on$drift_max)],
              paste(sprintf("%.4f", on$ocean_scale), collapse = ", ")))

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

# 4. The frame moduli follow the chunk (PI 2026-11-01): K_d = K_grain(IOM fraction) (1 - phi/phi_c) with the IOM's
#    share of the solid volume at t = 0 from 15 wt%, then from each chunk's remaining IOM (bulk_density()); at
#    phi >= phi_c the run is flagged and stopped.
kd_ok <- all(vapply(seq_len(n), function(i) abs(on$K_d_Pa[i] / (grain_elastic_moduli(on$iom_vol_frac[i])$K * (1 - on$phi_start[i] / 0.4)) - 1) < 1e-12, TRUE))
bd_prev <- if (n > 1) suppressWarnings(bulk_density(sub("\\.pqi$", ".csv", f[grepl(sprintf("chunk%03d", n - 1), f)]), porosity = on$phi_start[n - 1], origin_pqi = first))$summary else NULL
vf_ok <- abs(on$iom_vol_frac[1] / iom_volume_fraction(15) - 1) < 1e-12 &&
  (n == 1 || abs(on$iom_vol_frac[n] / iom_volume_fraction(15, iom_g = bd_prev$iom_g[nrow(bd_prev)], iom_density = bd_prev$iom_density[nrow(bd_prev)]) - 1) < 1e-12)
st <- suppressWarnings(run_tidal_evolution(Simulator$new(k = 1e-7, layer_thickness = 194200, config = EnceladusConfig$new(), phi0 = 0.4),
                                           organic_wt_percent = 15, n_orbits = 1, output_dir = file.path(od, "pc"), k_of_phi = kc, porewater_temp = 25))
check(kd_ok && vf_ok && identical(st$stopped, "phi_c") && is.null(st$chunks),
      sprintf("frame moduli per chunk: IOM fraction %s, K_d %s GPa; phi0 = phi_c: stopped (%s)",
              paste(sprintf("%.4f", on$iom_vol_frac), collapse = ", "), paste(sprintf("%.4g", on$K_d_Pa / 1e9), collapse = ", "), st$stopped))

unlink(od, recursive = TRUE)
cat(sprintf("\nSUMMARY: %d failure(s).\n", failures))
if (failures > 0) quit(status = 1)
