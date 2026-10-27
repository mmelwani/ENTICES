# The initial porewater: pure water, no oxidant or reductant, converging at every run temperature.
#
# Generators define the initial SOLUTIONs at 25 C with pe 4 and bring them to the run temperature by
# REACTION_TEMPERATURE (2026-10-27). Defined at the run temperature instead, pe 4 added 4.9 molal
# O(0) at 250 C and did not converge at 300 C, at any pressure.
#
# For tidal_kinetic, bulk_kinetic and bulk_equil at 25, 160 and 300 C, this takes the porewater
# SOLUTION 1 block and its REACTION_TEMPERATURE / REACTION_PRESSURE from the generated input, runs
# that heating step alone in PHREEQC, and checks: 0 errors, the run temperature and the profile
# porewater pressure, and O(0) and H(0) below 1e-8 molal. (Measured: O(0) 2e-14 molal at 160 C and
# 3e-9 at 300 C, i.e. ~5e-10 mol in the 0.157 kg porewater, against ~2.6 mol H2 from the metal.)
#
# Run from the repo root:  Rscript tests/initial_fluid_check.R

suppressPackageStartupMessages({ library(tidyverse); library(glue); library(R6) })
invisible(capture.output(source("ENTICES_LF_20260917.R")))

failures <- 0
check <- function(ok, what) {
  cat(if (isTRUE(ok)) "[PASS]" else "[FAIL]", what, "\n")
  if (!isTRUE(ok)) failures <<- failures + 1
}
quiet <- function(expr) { invisible(capture.output(r <- suppressMessages(expr))); r }

sim <- Simulator$new(k = 1e-07, layer_thickness = 194200, config = EnceladusConfig$new())
quiet(sim$run_simulation())
ph <- quiet(PhreeqcIntegrator$new(sim))
invisible(quiet(ph$calculate_phreeqc_parameters(organic_wt_percent = 0)))
pore_p <- sprintf("%.1f", profile_conditions(ph$profile_file, 194200)$pore_pressure_atm)

# The keyword block starting at a line matching `rx` (to the next blank line).
block <- function(l, rx) {
  i <- grep(rx, l)[1]
  if (is.na(i)) return(character(0))
  j <- i; while (j < length(l) && nzchar(trimws(l[j + 1]))) j <- j + 1
  l[i:j]
}
od <- tempfile("initial_fluid_")
gens <- list(
  tidal_kinetic = function(T) ph$generate_phreeqc_input_tidal_kinetic(output_dir = od, mode = "ungrouped", n_orbital_cycles = 1, porewater_temp = T),
  bulk_kinetic  = function(T) ph$generate_phreeqc_input_bulk_kinetic(output_dir = od, n_orbital_cycles = 1, mode = "porewater", water_rock_ratio = 1, porewater_temp = T),
  bulk_equil    = function(T) ph$generate_phreeqc_input_bulk_equil(output_dir = od, mode = "porewater", water_rock_ratio = 1, water_temp = T))
for (nm in names(gens)) for (T in c(25, 160, 300)) {
  l <- readLines(quiet(gens[[nm]](T)))
  db <- sub("^# Database: ", "", grep("^# Database: ", l, value = TRUE)[1])
  sol <- block(l, "^SOLUTION 1 ")
  rt  <- block(l, "^REACTION_TEMPERATURE 1$")
  rp  <- block(l, "^REACTION_PRESSURE 1$")
  inp <- tempfile(fileext = ".pqi"); out <- tempfile(fileext = ".out")
  writeLines(c(sol, "END", "USE solution 1", rt, rp, "END"), inp)
  r <- run_phreeqc(inp, db, output = out, warn_version = FALSE)
  o <- r$output_lines
  last <- function(rx) { v <- grep(rx, o, value = TRUE); if (length(v)) tail(v, 1) else NA_character_ }
  num <- function(x) suppressWarnings(as.numeric(sub("^\\S+\\s+(\\S+).*", "\\1", x)))
  tc <- as.numeric(sub(".*=\\s*", "", last("Temperature \\(.+C\\)\\s+=")))
  pa <- sub(".*=\\s*", "", last("Pressure \\(atm\\)\\s+="))
  o0 <- num(last("^O\\(0\\)\\s")); h0 <- num(last("^H\\(0\\)\\s"))
  o0 <- if (is.na(o0)) 0 else o0; h0 <- if (is.na(h0)) 0 else h0
  check(r$errors == 0 && identical(sol[2], "    temp      25") && isTRUE(abs(tc - T) < 0.01) &&
          isTRUE(abs(as.numeric(pa) - as.numeric(pore_p)) < 0.05) && o0 < 1e-8 && h0 < 1e-8,
        sprintf("%-13s %3d C: %d errors, at %.2f C and %s atm, O(0) %.2e, H(0) %.2e molal", nm, T, r$errors, tc, trimws(pa), o0, h0))
}

unlink(od, recursive = TRUE)
cat(sprintf("\nSUMMARY: %d failure(s).\n", failures))
if (failures > 0) quit(status = 1)
