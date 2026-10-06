# =============================================================================
# iom_selftest.R  --  Standalone tests for iom_module.R
#
# Runs WITHOUT ENTICES, WITHOUT PHREEQC, and WITHOUT circulation physics.
# Three kinds of checks:
#   (A) Pure-R checks on the module itself: config validity, balance
#       arithmetic, string generation -- fast, no dependencies beyond base R.
#   (B) A from-scratch Basic-language-equivalent integration of the Arrhenius
#       rate law in R, run at Miller's own experimental conditions (48h
#       isothermal at 250/350/500 C), to check whether the CALIBRATED
#       channel set reproduces his measured gas yields and (more
#       importantly, since yields were fitted) his RESIDUE composition
#       trajectory -- the test this module has never actually been run
#       against.
#   (C) Structural checks on redox_mode (Option 3, coupled/decoupled redox):
#       confirms the mechanics of switching formulas work and bad input is
#       rejected. Deliberately does NOT and CANNOT check the thing
#       PROGRESS.md Part 3 actually cares about (whether decoupled mode
#       conserves elements in real PHREEQC output) -- that is blocked on
#       PHREEQC access, not something this script fakes. See Part C's own
#       header and iom_default_config()'s docstring.
#
# This does NOT call PHREEQC. It re-implements the same first-order
# Arrhenius ODE that the generated RATES blocks encode, so it validates the
# KINETICS/RATES numbers before they ever reach PHREEQC. It is not a
# substitute for actually running the generated PHREEQC text (do that too,
# once ENTICES is wired up, and once PHREEQC access is available in this
# environment) -- it is a fast, dependency-free check of the calibration
# itself.
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

# CO2/CH4/CHn are now several sub-pool rows each (IOM_CO2_44, IOM_CO2_46, ...)
# sharing one -formula; sum by family to get the quantity Miller actually
# measured. Works for both multi-pool families and the still-single-pool
# IOM_N/IOM_S (falls back to exact-name match).
family_sum <- function(rel, prefix) {
  sum(rel[names(rel) == prefix | startsWith(names(rel), paste0(prefix, "_"))])
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

# BASIS NOTE (found 2026-09-18 while re-deriving IOM_N's Ea, see
# iom_default_config()'s docstring for the full story): cfg's m0_per_kg is on
# the ENTICES-normalized bulk basis (bulk() = mol per 1000 g of C+H+O+N+S
# only), NOT the same basis as Miller's raw per-kg-of-real-Murchison-IOM
# measurements -- real Murchison IOM is only ~82% organic elements by mass
# (Table 1's wt% sum to 82.1%; the rest is ash/S/P), so ENTICES's kg is a
# renormalization, not the same kg. Comparing an ENTICES-basis prediction
# directly to Miller's raw value is dimensionally wrong; divide by the
# element-specific scale factor first. This bug was previously latent in
# EVERY single-Ea fit (CO2/CH4/N) -- it happened not to matter for 350C
# (the sole free Ea parameter absorbed it) but produced the ~20-22% "500C
# overshoot" that earlier sessions attributed (incorrectly) to "m0 being a
# fixed exhaustion pool." The distributed CO2/CH4/CHn fit avoids the bug by
# construction (shift fit from the ratio of two Murchison-basis numbers,
# basis-independent; m0 scaled up only afterward) -- so only N's comparison
# below still needs the explicit correction.
c_scale <- iom_bulk_formula()[["C"]] / 49.871   # Murchison C, Table 1
n_scale <- iom_bulk_formula()[["N"]] / murchison_N0
cat(sprintf("Basis scale factors (ENTICES-normalized / real-Murchison): C=%.4f  N=%.4f\n\n",
            c_scale, n_scale))

cat(sprintf("%-8s %10s %10s %10s | %10s %10s %10s | %10s %10s %10s | %10s\n",
            "T(C)", "CO2_pred", "CO2_meas", "ratio", "CH4_pred", "CH4_meas", "ratio",
            "N_pred", "N_meas", "ratio", "CHn_pred"))
released_by_T <- list()
for (i in seq_len(nrow(miller_measured))) {
  T_C <- miller_measured$T_C[i]
  rel <- run_isothermal(cfg, T_C, dur_s)
  released_by_T[[as.character(T_C)]] <- rel
  co2p <- family_sum(rel, "IOM_CO2") / c_scale; ch4p <- family_sum(rel, "IOM_CH4") / c_scale
  np <- family_sum(rel, "IOM_N") / n_scale; chnp <- family_sum(rel, "IOM_CHn") / c_scale
  co2m <- miller_measured$CO2_mol_kg[i]; ch4m <- miller_measured$CH4_mol_kg[i]
  nm <- miller_measured$NH3_pct_of_N[i]/100 * murchison_N0
  cat(sprintf("%-8g %10.3f %10.3f %10.2f | %10.3f %10.3f %10.2f | %10.3f %10.3f %10.2f | %10.3f\n",
              T_C, co2p, co2m, co2p/co2m, ch4p, ch4m, ch4p/ch4m, np, nm, np/nm, chnp))
}
cat("(All *_pred values above are ENTICES-basis model output divided by the\n")
cat(" relevant basis scale factor, i.e. converted BACK to Miller's raw\n")
cat(" Murchison basis for a fair comparison -- see the basis note above.\n")
cat(" CO2/CH4 are distributed-Ea sub-pools fit to BOTH 350C and 500C as\n")
cat(" exact targets, so both ratio columns should read ~1.0. IOM_N is still\n")
cat(" single-Ea (no distribution data exists for it), fit only at 350C, so\n")
cat(" its 500C ratio is a genuine PREDICTION. CHn has no Miller target at\n")
cat(" all (GC doesn't report C2+) -- printed for transparency only.)\n\n")

cat("Asserted benchmark checks (tolerances explained inline):\n")
# CO2/CH4 350C AND 500C are now both exact fit targets (2 unknowns -- shift,
# m0 -- solved against these 2 equations), so both should match tightly once
# converted back to Miller's raw basis. 1% tolerance covers the ~0.3-0.5%
# numerical error of the 2000-step explicit-Euler ODE relative to the
# closed-form fit.
rel350 <- released_by_T[["350"]]; rel500 <- released_by_T[["500"]]
m350 <- miller_measured[miller_measured$T_C == 350, ]
m500 <- miller_measured[miller_measured$T_C == 500, ]
nm350 <- m350$NH3_pct_of_N / 100 * murchison_N0
nm500 <- m500$NH3_pct_of_N / 100 * murchison_N0
for (ch in c("CO2", "CH4")) {
  fam <- paste0("IOM_", ch)
  pass(sprintf("350C %s matches Miller's absolute measured yield (fit target, 1%% tol)", ch),
       abs((family_sum(rel350, fam)/c_scale) / m350[[sprintf("%s_mol_kg", ch)]] - 1) < 0.01)
  pass(sprintf("500C %s matches Miller's absolute measured yield (fit target, 1%% tol)", ch),
       abs((family_sum(rel500, fam)/c_scale) / m500[[sprintf("%s_mol_kg", ch)]] - 1) < 0.01)
}
pass("350C N matches Miller's absolute measured yield (fit target, 1% tol)",
     abs((family_sum(rel350, "IOM_N")/n_scale) / nm350 - 1) < 0.01)
# 500C N is a genuine prediction (single-Ea, no distribution) -- expected to
# sit close to full exhaustion of its own m0 pool (Murchison-basis-equivalent
# ~1.0 vs Miller, since m0 was itself calibrated so 500C approximates full
# exhaustion). Loose band, informational rather than a strict target.
n500_ratio <- (family_sum(rel500, "IOM_N")/n_scale) / nm500
pass("500C N ratio is a sane prediction (informational, 0.95-1.10 loose band)",
     n500_ratio > 0.95 && n500_ratio < 1.10, sprintf("ratio=%.3f", n500_ratio))
cat("\n")

# ---------------------------------------------------------------------------
# Independent shape validation (NOT fit to any of our data): the UNSHIFTED
# Vitrimat-2018 CO2 distribution should reproduce the RATIO between Miller's
# HC113 250C and 500C measured yields. This is the strongest evidence the
# module has for using this distribution's shape at all -- see
# PATCH_F_vitrimat2018.md sec 3. Deliberately reimplemented here from the
# raw bins/weights, independent of cfg (cfg's CO2 rows are shifted+fit to
# Murchison, not HC113, and are the wrong thing to test this against).
cat("Independent validation: unshifted Vitrimat-2018 CO2 shape vs HC113\n")
cat("(not fit to any of our data -- a genuine out-of-sample check):\n")
co2_bins <- c(44, 46, 48, 50, 52, 54, 56)
co2_w    <- c(10, 15, 15, 15, 15, 15, 15) / 100
frac_released_shape <- function(bins, w, shift_kcal, T_C, t, A) {
  TK <- T_C + 273.15
  Ea_J <- (bins + shift_kcal) * 4184
  k <- A * exp(-Ea_J / (R_GAS * TK))
  sum(w * (1 - exp(-k * t)))
}
f250 <- frac_released_shape(co2_bins, co2_w, 0, 250, dur_s, 2e15)
f500 <- frac_released_shape(co2_bins, co2_w, 0, 500, dur_s, 2e15)
pred_ratio <- f250 / f500
hc113_ratio <- 1.158 / 2.55
cat(sprintf("  predicted 250/500 ratio = %.4f, HC113 measured = %.4f, agreement = %.3f\n",
            pred_ratio, hc113_ratio, pred_ratio / hc113_ratio))
pass("Unshifted Vitrimat-2018 CO2 shape reproduces HC113 250/500 ratio (2% tol)",
     abs(pred_ratio / hc113_ratio - 1) < 0.02)
cat("\n")

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
  h_c_meas <- miller_residue$H_C[miller_residue$T_C == T_C]
  o_c_meas <- miller_residue$O_C[miller_residue$T_C == T_C]
  h_c_pred <- resid_now["H"]/resid_now["C"]
  o_c_pred <- resid_now["O"]/resid_now["C"]
  cat(sprintf("  T=%gC: predicted residue H/C=%.3f O/C=%.4f  |  measured H/C=%s O/C=%s\n",
              T_C, h_c_pred, o_c_pred, h_c_meas, o_c_meas))
  # Loose (25%) sanity bound only -- NOT a tight calibration target. The
  # residue was never part of the two-point Ea fit (that's the whole point
  # of testing it), and the model is known to be missing a low-Ea decay
  # tail (see handover notes), which biases residue H/C, O/C high. This
  # check exists to catch a GROSS error (e.g. a sign flip or unit mixup in
  # iom_residue_formula/formula parsing), not to enforce a match that the
  # model isn't expected to hit yet.
  if (!is.na(h_c_meas)) {
    pass(sprintf("%gC residue H/C within loose 25%% sanity bound", T_C),
         abs(h_c_pred / h_c_meas - 1) < 0.25,
         sprintf("pred=%.3f meas=%.3f", h_c_pred, h_c_meas))
  }
  if (!is.na(o_c_meas)) {
    pass(sprintf("%gC residue O/C within loose 25%% sanity bound", T_C),
         abs(o_c_pred / o_c_meas - 1) < 0.25,
         sprintf("pred=%.4f meas=%.4f", o_c_pred, o_c_meas))
  }
}

cat("\n=============================================================\n")
cat("PART C: redox_mode (Option 3) structural checks\n")
cat("=============================================================\n")
cat("Structure checks on the redox_mode switch. The FORMULAS themselves were\n")
cat("tested in PHREEQC on 2026-09-20 -- see PHREEQC_TEST_RESULTS_2026-09-20.md\n")
cat("for the measured numbers (Mtg 1 H 0.435 reproduces the coupled aqueous\n")
cat("state exactly; Sg 1 H 1 matches the coupled H release, while the\n")
cat("previously-proposed Sg 1 H -1 was wrong by 2 mol H per mol S). What is\n")
cat("checked HERE is only that the switch wires those strings up correctly.\n\n")

pass("default redox_mode (all coupled) is unchanged from calling with no args",
     identical(iom_default_config(), iom_default_config(c(C = "coupled", N = "coupled", S = "coupled"))))

cfg_dS <- iom_default_config(c(C = "coupled", N = "coupled", S = "decoupled"))
pass("decoupled S uses the PHREEQC-verified 'Sg 1 H 1' formula",
     cfg_dS$formula[cfg_dS$name == "IOM_S"] == "Sg 1 H 1")
pass("decoupled S leaves every other channel's formula untouched",
     all(cfg_dS$formula[cfg_dS$name != "IOM_S"] == cfg$formula[cfg$name != "IOM_S"]))

cfg_dC <- iom_default_config(c(C = "decoupled", N = "coupled", S = "coupled"))
pass("decoupled C switches IOM_CH4 sub-pools to 'Mtg 1 H 0.435'",
     all(cfg_dC$formula[startsWith(cfg_dC$name, "IOM_CH4")] == "Mtg 1 H 0.435"))
# IOM_CO2/IOM_CHn deliberately keep their coupled formulas: their decoupling
# is the caller entering H2 as Hdg in the solution, which no -formula can
# express. Verified in PHREEQC (carbon reduction suppressed ~26,000x).
pass("decoupled C leaves IOM_CO2 and IOM_CHn formulas unchanged (Hdg is solution-level)",
     all(cfg_dC$formula[startsWith(cfg_dC$name, "IOM_CO2")] ==
         cfg$formula[startsWith(cfg$name, "IOM_CO2")]) &&
     all(cfg_dC$formula[startsWith(cfg_dC$name, "IOM_CHn")] ==
         cfg$formula[startsWith(cfg$name, "IOM_CHn")]))
pass("decoupled C attaches the Hdg reminder attribute",
     !is.na(attr(cfg_dC, "carbon_decoupling_note")) &&
     grepl("Hdg", attr(cfg_dC, "carbon_decoupling_note")))

# Decoupling N is a decision made AGAINST, not an unimplemented gap.
pass("requesting decoupled N is refused, with the reason",
     tryCatch({ iom_default_config(c(C = "coupled", N = "decoupled", S = "coupled")); FALSE },
              error = function(e) grepl("Ntg|deliberately", conditionMessage(e))))

pass("invalid redox_mode value is rejected",
     tryCatch({ iom_default_config(c(C = "coupled", N = "sort-of", S = "coupled")); FALSE },
              error = function(e) TRUE))

# NOT a bug: validate_config()'s balance check expects real element symbols,
# so it rejects "Mtg"/"Sg". Left as an honest failure -- see the docstring.
pass("decoupled formulas are correctly REJECTED by validate_config (known limitation)",
     tryCatch({ iom_validate_config(cfg_dS); FALSE },
              error = function(e) grepl("Sg", conditionMessage(e))))
cat("\n")

cat("\n=============================================================\n")
cat(sprintf("SUMMARY: %d failure(s) across all asserted checks.\n", failures))
cat("Gas-yield ratios, the independent HC113 shape check, and residue\n")
cat("sanity bounds above are all asserted (pass/fail), not just printed.\n")
cat("500C residue H/C, O/C remain unmeasured in Miller's own Table 5\n")
cat("Murchison row -- not checked. IOM_CHn and IOM_S have no Miller target\n")
cat("at all and are reported for transparency only, never asserted.\n")
cat("=============================================================\n")

if (failures > 0) quit(status = 1)
