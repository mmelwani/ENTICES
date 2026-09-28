setwd("C:/Users/melwa/GitHub/ENTICES")
source("iom_module.R")
cfg <- iom_default_config()
R_GAS <- 8.314
dur_s <- 48*3600

run_isothermal <- function(config, T_C, duration_s) {
  TK <- T_C + 273.15
  k <- 10^config$logA * exp(-config$Ea_J/(R_GAS*TK))
  rel <- config$m0_per_kg * (1 - exp(-k*duration_s))   # closed form, exact
  names(rel) <- config$name
  rel
}

bulk <- iom_bulk_formula()
inert <- iom_residue_formula(cfg)

residue_at <- function(T_C) {
  rel <- run_isothermal(cfg, T_C, dur_s)
  m_rem <- cfg$m0_per_kg - rel
  r <- inert
  for (j in seq_len(nrow(cfg))) {
    f <- iom_parse_formula(cfg$formula[j])
    for (el in names(f)) r[el] <- r[el] + m_rem[j]*f[el]
  }
  r
}

cat("=== STARTING COMPOSITION: model vs Murchison ===\n")
bulk_HC <- bulk[["H"]]/bulk[["C"]]; bulk_OC <- bulk[["O"]]/bulk[["C"]]
cat(sprintf("ENTICES bulk formula : H/C = %.4f   O/C = %.4f\n", bulk_HC, bulk_OC))
cat("Miller Table 1 Murch. : H/C = 0.70     O/C = 0.18\n")
cat(sprintf("  -> model starts %.1f%% too H-rich, %.1f%% too O-poor\n",
            (bulk_HC/0.70-1)*100, (bulk_OC/0.18-1)*100))

cat("\n=== 350 C RESIDUE: absolute vs dimensionless R-ratio ===\n")
r350 <- residue_at(350)
p_HC <- r350[["H"]]/r350[["C"]]; p_OC <- r350[["O"]]/r350[["C"]]
m_HC <- 0.63; m_OC <- 0.096            # Table 5, Murchison 350-3-M
m_RHC <- 0.90; m_ROC <- 0.53           # Table 5, same row, RX/C columns
cat(sprintf("ABSOLUTE  H/C: pred %.4f vs meas %.3f  -> %+.1f%%\n", p_HC, m_HC, (p_HC/m_HC-1)*100))
cat(sprintf("ABSOLUTE  O/C: pred %.4f vs meas %.3f  -> %+.1f%%\n", p_OC, m_OC, (p_OC/m_OC-1)*100))
pred_RHC <- p_HC/bulk_HC; pred_ROC <- p_OC/bulk_OC
cat(sprintf("R-RATIO  RH/C: pred %.4f vs meas %.2f (+/-0.12) -> %+.1f%%\n",
            pred_RHC, m_RHC, (pred_RHC/m_RHC-1)*100))
cat(sprintf("R-RATIO  RO/C: pred %.4f vs meas %.2f (+/-0.06) -> %+.1f%%\n",
            pred_ROC, m_ROC, (pred_ROC/m_ROC-1)*100))

cat("\n=== IS THE ABSOLUTE OVERSHOOT JUST THE STARTING-COMPOSITION OFFSET? ===\n")
cat(sprintf("ratio of predicted/measured FINAL   H/C = %.4f\n", p_HC/m_HC))
cat(sprintf("ratio of model/measured  STARTING   H/C = %.4f\n", bulk_HC/0.70))
cat(sprintf("  difference between the two        = %+.2f%%\n",
            ((p_HC/m_HC)/(bulk_HC/0.70)-1)*100))
