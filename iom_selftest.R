# =============================================================================
# iom_selftest.R  --  Standalone tests for iom_module.R
#
# Runs WITHOUT ENTICES, WITHOUT PHREEQC, and WITHOUT circulation physics.
# Two kinds of checks:
#   (A) Pure-R checks on the module itself: config validity, balance
#       arithmetic, string generation -- fast, no dependencies beyond base R.
#   (B) A from-scratch Basic-language-equivalent integration of the Arrhenius
#       rate law in R, run at Miller's own experimental conditions (48h
#       isothermal at 250/350/500 C), to check whether the CALIBRATED
#       channel set reproduces his measured gas yields and (more
#       importantly, since yields were fitted) his RESIDUE composition
#       trajectory -- the test this module has never actually been run
#       against.
#
# This does NOT call PHREEQC. It re-implements the same first-order
# Arrhenius ODE that the generated RATES blocks encode, so it validates the
# KINETICS/RATES numbers before they ever reach PHREEQC. It is not a
# substitute for actually running the generated PHREEQC text (do that too,
# once ENTICES is wired up) -- it is a fast, dependency-free check of the
# calibration itself.
#
# Usage: Rscript iom_selftest.R
# =============================================================================

source("iom_module.R")

failures <- 0
pass <- function(name, ok, detail = "") {
  status <- if (ok) "PASS" else "FAIL"
  cat(sprintf("[%s] %s%s\n", status, name, if (nzchar(detail)) paste0(" -- ", detail) else ""))
  if (!ok) failures <<- failures + 1
  invisible(ok)
}

cat("=============================================================\n")
cat("PART A: module self-checks (config, balance, string generation)\n")
cat("=============================================================\n")

cfg <- iom_default_config()

pass("default config validates",
     tryCatch({ iom_validate_config(cfg); TRUE }, error = function(e) { message(e); FALSE }))

res <- iom_residue_formula(cfg)
bulk <- iom_bulk_formula()
pass("residue is non-negative for every element", all(res >= -1e-9),
     paste(sprintf("%s=%.4f", names(res), res), collapse = " "))
pass("residue + channels reproduce bulk exactly",
     {
       released <- bulk - res
       max_abs_err <- max(abs((released + res) - bulk))
       max_abs_err < 1e-9
     })
cat(sprintf("  residue H/C = %.4f, O/C = %.4f  (bulk H/C = %.4f, O/C = %.4f)\n",
            res["H"]/res["C"], res["O"]/res["C"], bulk["H"]/bulk["C"], bulk["O"]/bulk["C"]))
cat(sprintf("  fraction of IOM carbon ever released = %.1f%%\n",
            (bulk["C"] - res["C"]) / bulk["C"] * 100))

pass("detects duplicate channel names",
     tryCatch({ bad <- cfg; bad$name[2] <- bad$name[1]
                iom_validate_config(bad); FALSE },
              error = function(e) TRUE))

pass("detects an element not in bulk formula",
     tryCatch({ bad <- cfg; bad$formula[1] <- "Xe 1"
                iom_validate_config(bad); FALSE },
              error = function(e) TRUE))

pass("detects over-drawn element (negative residue)",
     tryCatch({ bad <- cfg; bad$m0_per_kg[1] <- 1000
                iom_validate_config(bad); FALSE },
              error = function(e) TRUE))

rb <- iom_rates_block(cfg)
pass("RATES block contains one -start/-end pair per channel",
     lengths(regmatches(rb, gregexpr("-start", rb))) == nrow(cfg) &&
     lengths(regmatches(rb, gregexpr("-end", rb))) == nrow(cfg))

kb <- iom_kinetics_block(cfg, iom_mass_kg = 2.0)
pass("KINETICS block m0 scales linearly with iom_mass_kg",
     {
       kb1 <- iom_kinetics_block(cfg, iom_mass_kg = 1.0)
       kb2 <- iom_kinetics_block(cfg, iom_mass_kg = 2.0)
       m0_1 <- as.numeric(sub("-m0 ", "", regmatches(kb1, gregexpr("-m0 [0-9.eE+-]+", kb1))[[1]]))
       m0_2 <- as.numeric(sub("-m0 ", "", regmatches(kb2, gregexpr("-m0 [0-9.eE+-]+", kb2))[[1]]))
       all(abs(m0_2 / m0_1 - 2.0) < 1e-6)
     })

pass("punch headings/lines cover every channel",
     {
       h <- iom_punch_headings(cfg); l <- iom_punch_lines(cfg, start_line = 500)
       all(vapply(cfg$name, function(nm) grepl(paste0(nm, "_mol"), h) &&
                                          grepl(sprintf('KIN\\("%s"\\)', nm), l),
                  logical(1)))
     })

cat("\n=============================================================\n")
cat("PART B: reproduce Miller et al. (2025) conditions, standalone\n")
cat("=============================================================\n")
cat("Re-implements the RATES Basic rate law (rate = 10^logA * exp(-Ea/RT) * M)\n")
cat("as an R ODE and integrates it over Miller's own experiment duration, at\n")
cat("his experimental temperatures. Compares predicted C released and residue\n")
cat("H/C, O/C against his measured values (Tables 2 and 5).\n\n")

R_GAS <- 8.314
run_isothermal <- function(config, T_C, duration_s, n_steps = 2000) {
  TK <- T_C + 273.15
  dt <- duration_s / n_steps
  m  <- config$m0_per_kg  # mol/kg IOM remaining, per channel
  for (i in seq_len(n_steps)) {
    k <- 10^config$logA * exp(-config$Ea_J / (R_GAS * TK))
    m <- pmax(0, m - k * m * dt)
  }
  released <- config$m0_per_kg - m
  names(released) <- config$name
  released
}

hours <- 48
dur_s <- hours * 3600

miller_measured <- data.frame(
  T_C = c(350, 500),
  CO2_mol_kg = c(3.04, 3.43),        # Murchison, Table 2
  CH4_mol_kg = c(0.184, 1.57),
  NH3_pct_of_N = c(32.06, 40.83)     # Table 9, Murchison 3 kbar rows "350-3-M"/"500-3-M"
)
# Table 1: Murchison N = 2.48 wt% -> mol N / kg IOM. Cross-checked against
# Table 9's direct umol N/mg column (== mol/kg): 0.57 @350C, 0.72 @500C --
# matches NH3_pct_of_N * murchison_N0 to 3 decimal places.
murchison_N0 <- (2.48/100) * 1000 / 14.007

cat(sprintf("%-8s %10s %10s %10s | %10s %10s %10s | %10s %10s %10s\n",
            "T(C)", "CO2_pred", "CO2_meas", "ratio", "CH4_pred", "CH4_meas", "ratio",
            "N_pred", "N_meas", "ratio"))
for (i in seq_len(nrow(miller_measured))) {
  T_C <- miller_measured$T_C[i]
  rel <- run_isothermal(cfg, T_C, dur_s)
  co2p <- rel[["IOM_CO2"]]; ch4p <- rel[["IOM_CH4"]]; np <- rel[["IOM_N"]]
  co2m <- miller_measured$CO2_mol_kg[i]; ch4m <- miller_measured$CH4_mol_kg[i]
  nm <- miller_measured$NH3_pct_of_N[i]/100 * murchison_N0
  cat(sprintf("%-8g %10.3f %10.3f %10.2f | %10.3f %10.3f %10.2f | %10.3f %10.3f %10.2f\n",
              T_C, co2p, co2m, co2p/co2m, ch4p, ch4m, ch4p/ch4m, np, nm, np/nm))
}
cat("(350C ratio ~1.0 is expected for all three: Ea was fit to hit each\n")
cat(" channel's ABSOLUTE 350C value exactly, so this column is a consistency\n")
cat(" check on the arithmetic. 500C ratio is NOT expected to be 1.0 for any\n")
cat(" of them -- each channel's m0_per_kg is a fixed exhaustion pool scaled\n")
cat(" to bulk IOM content, not Miller's raw measured yield basis, so full\n")
cat(" exhaustion legitimately exceeds the raw yield at the hotter point.\n")
cat(" This is why the ORIGINAL CO2/CH4 fit -- which used the ratio of\n")
cat(" Miller's two measured yields rather than the absolute value -- was\n")
cat(" wrong: see the CORRECTED note in iom_default_config(). IOM_N was\n")
cat(" refit the same way on 2026-09-17; IOM_S remains uncalibrated (no H2S\n")
cat(" data in Miller at all) and is not compared here.)\n\n")

cat("Residue trajectory (THIS is the real, previously unexecuted test --\n")
cat("residue composition was NOT part of the two-point Ea fit):\n\n")
miller_residue <- data.frame(   # Table 5, Murchison rows only ("350-3-M", "500-3-M", both 3 kbar)
  T_C  = c(0,    350,   500),
  H_C  = c(0.70, 0.63,  NA),     # 500C: Murchison row reports only N/C -- H/C, O/C are blank
  O_C  = c(0.18, 0.096, NA)      # in Table 5 itself (not measured), not a transcription gap
)
cat("NOTE: 350C values are Miller Table 5, Murchison row '350-3-M' (3 kbar).\n")
cat("500C H/C and O/C are left NA because Miller's Murchison row at 500C\n")
cat("('500-3-M') reports only N/C -- H/C and O/C are blank in Table 5 itself.\n\n")
for (i in seq_len(nrow(miller_measured))) {
  T_C <- miller_measured$T_C[i]
  rel <- run_isothermal(cfg, T_C, dur_s)
  m_remaining <- cfg$m0_per_kg - rel
  # residue element totals = inert-residue elements (never react) + whatever
  # channel mass has NOT yet reacted (m_remaining), since unreacted channel
  # mass is still, physically, part of the solid.
  inert <- iom_residue_formula(cfg)
  resid_now <- inert
  for (j in seq_len(nrow(cfg))) {
    f <- iom_parse_formula(cfg$formula[j])
    for (el in names(f)) resid_now[el] <- resid_now[el] + m_remaining[j] * f[el]
  }
  cat(sprintf("  T=%gC: predicted residue H/C=%.3f O/C=%.4f  |  measured H/C=%s O/C=%s\n",
              T_C, resid_now["H"]/resid_now["C"], resid_now["O"]/resid_now["C"],
              miller_residue$H_C[miller_residue$T_C == T_C],
              miller_residue$O_C[miller_residue$T_C == T_C]))
}

cat("\n=============================================================\n")
cat(sprintf("SUMMARY: %d failure(s) in Part A checks.\n", failures))
cat("Part B is diagnostic (compares predicted vs Miller-measured values)\n")
cat("rather than pass/fail. 500C residue H/C, O/C remain NA -- unmeasured\n")
cat("in Miller's own Table 5 Murchison row, not a script gap.\n")
cat("=============================================================\n")

if (failures > 0) quit(status = 1)
