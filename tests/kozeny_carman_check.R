# Permeability and surface area from grain size (item 3, PI 2026-10-29).
#
#  1. kozeny_carman_k() reproduces the PI's table (m2, within rounding):
#        phi   d=1um   d=30um  d=1mm   d=2cm
#        0.2   7e-17   6e-14   7e-11   3e-8
#        0.3   3e-16   3e-13   3e-10   1e-7
#  2. Simulator$new(k = "kozeny_carman", grain_diameter = d) sets k from d and the porosity; a number for k
#     still works (old runs); "kozeny_carman" without grain_diameter is an error.
#  3. ssa = "geometric": every mineral's SSA line in a generated tidal_kinetic input is 6/(rho d) x roughness
#     (rho = M/Vm from data/mineral_molar_volumes.csv), per-mineral diameters override; "legacy" leaves the
#     0.1 / 5 m2/g lines.
#
# Generates inputs only. Run from the repo root:  Rscript tests/kozeny_carman_check.R

suppressPackageStartupMessages({ library(tidyverse); library(glue); library(R6) })
invisible(capture.output(source("ENTICES_LF_20260917.R")))
failures <- 0
check <- function(ok, what) {
  cat(if (isTRUE(ok)) "[PASS]" else "[FAIL]", what, "\n")
  if (!isTRUE(ok)) failures <<- failures + 1
}
quiet <- function(e) { invisible(capture.output(r <- suppressWarnings(suppressMessages(e)))); r }

# 1. The table.
d <- c(1e-6, 30e-6, 1e-3, 2e-2)
tab <- rbind(`0.2` = c(7e-17, 6e-14, 7e-11, 3e-8), `0.3` = c(3e-16, 3e-13, 3e-10, 1e-7))
got <- rbind(kozeny_carman_k(d, 0.2), kozeny_carman_k(d, 0.3))
check(all(signif(got, 1) == tab), sprintf("Kozeny-Carman table: phi 0.2: %s; phi 0.3: %s",
                                           paste(sprintf("%.2g", got[1, ]), collapse = " "), paste(sprintf("%.2g", got[2, ]), collapse = " ")))
check(abs(kozeny_carman_k(kozeny_carman_d(1e-13, 0.32), 0.32) / 1e-13 - 1) < 1e-12,
      sprintf("inverse: phi 0.32 needs d = %.3g-%.3g m for k 1e-14-1e-13 (Choblet), %.3g m for 1e-7",
              kozeny_carman_d(1e-14, 0.32), kozeny_carman_d(1e-13, 0.32), kozeny_carman_d(1e-7, 0.32)))

# 2. The Simulator.
s <- Simulator$new(k = "kozeny_carman", layer_thickness = 194200, grain_diameter = 1e-5)
e <- tryCatch({ Simulator$new(k = "kozeny_carman", layer_thickness = 194200); "no error" }, error = function(e) conditionMessage(e))
check(abs(s$params$k / kozeny_carman_k(1e-5, 0.32) - 1) < 1e-12 && s$params$permeability == "kozeny_carman" &&
        Simulator$new(k = 1e-7, layer_thickness = 194200)$params$k == 1e-7 && grepl("needs grain_diameter", e),
      sprintf("k = \"kozeny_carman\", d 10 um, phi 0.32: k %.3g m2; numeric k kept; no grain_diameter errors", s$params$k))

# 3. Geometric SSA in a generated input.
od <- tempfile("kc_")
gen <- function(sim, dir) quiet({ sim$run_simulation(); ph <- PhreeqcIntegrator$new(sim); ph$calculate_phreeqc_parameters(organic_wt_percent = 15)
  ph$generate_phreeqc_input_tidal_kinetic(output_dir = file.path(od, dir), mode = "ungrouped", n_orbital_cycles = 1, porewater_temp = 25, trace_C_molal = 1e-14) })
ssa_lines <- function(f) { l <- readLines(f); i <- grep("SSA = [0-9.eE+-]+", l)
  nm <- vapply(i, function(j) { k <- max(grep("^[A-Za-z]", l[seq_len(j)])); sub("[ \t#].*$", "", l[k]) }, "")
  stats::setNames(as.numeric(sub(".*SSA = ([0-9.eE+-]+).*", "\\1", l[i])), nm) }
v <- read_mineral_volumes(); rho <- function(m) v$molar_mass_g_mol[v$phase == m] / v$Vm_cm3_mol[v$phase == m]
leg <- ssa_lines(gen(Simulator$new(k = 1e-7, layer_thickness = 194200), "legacy"))
geo <- ssa_lines(gen(Simulator$new(k = "kozeny_carman", layer_thickness = 194200, grain_diameter = 1e-4, ssa = "geometric", roughness = 3,
                                   grain_diameter_by_mineral = c(Albite = 1e-5)), "geometric"))
want <- vapply(names(geo), function(m) 3 * geometric_ssa_m2_g(if (m == "Albite") 1e-5 else 1e-4, rho(m)), 0)
check(all(leg %in% c(0.1, 5)) && length(geo) == length(leg) && all(abs(geo / want - 1) < 1e-5),
      sprintf("SSA: legacy %s; geometric (d 0.1 mm, albite 10 um, roughness 3) %s m2/g",
              paste(sprintf("%s %g", names(leg), leg), collapse = ", "), paste(sprintf("%s %.3g", names(geo), geo), collapse = ", ")))
unlink(od, recursive = TRUE)
cat(sprintf("\nSUMMARY: %d failure(s).\n", failures))
if (failures > 0) quit(status = 1)
