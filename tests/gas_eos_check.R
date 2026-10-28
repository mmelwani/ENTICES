# Peng-Robinson gases (scripts/gas_eos.py), and how much PR changes H2 at porewater conditions.
#
#  1. Every gas the generators put in a GAS_PHASE has -T_c, -P_c and -Omega in all three databases
#     (PHREEQC uses Peng-Robinson only if every component of a gas phase has all three).
#  2. H2 at 160 C, 261.2 atm, into 0.157 kg of pure water through an empty fixed-pressure gas phase:
#     Peng-Robinson (the databases) against ideal (a copy with the restored constants commented out).
#     Prints the saturation molality of dissolved H2 and the bubble's molar volume for each, and
#     checks that both run, that PR's fugacity coefficient is > 1 at these conditions, and that the
#     two differ by less than 15% (measured 2026-10-28: +8.3% dissolved, +8.7% molar volume).
#
# Run from the repo root:  Rscript tests/gas_eos_check.R

suppressPackageStartupMessages({ library(tidyverse); library(glue); library(R6) })
invisible(capture.output(source("ENTICES_LF_20260917.R")))

failures <- 0
check <- function(ok, what) {
  cat(if (isTRUE(ok)) "[PASS]" else "[FAIL]", what, "\n")
  if (!isTRUE(ok)) failures <<- failures + 1
}

# 1. Critical constants.
gas_consts <- function(db) {
  l <- sub("#.*", "", readLines(db, warn = FALSE, encoding = "UTF-8"))
  out <- list(); cur <- NULL; sec <- ""
  for (s in l) {
    if (trimws(s) %in% c("PHASES", "RATES", "SOLUTION_SPECIES", "END")) { sec <- trimws(s); cur <- NULL; next }
    if (sec != "PHASES") next
    if (grepl("^[A-Za-z]", s) && !grepl("=", s)) { cur <- sub("\\s.*", "", s); out[[cur]] <- character(0); next }
    if (!is.null(cur)) out[[cur]] <- c(out[[cur]], regmatches(s, gregexpr("-(T_c|P_c|Omega)", s))[[1]])
  }
  out
}
for (db in list.files("databases", "\\.dat$", full.names = TRUE)) {
  mode <- names(CARBON_MODE_DATABASES)[match(basename(db), CARBON_MODE_DATABASES)]
  comps <- c(gas_components(list(decoupled = mode != "coupled", database = db)), "Hdg(g)", "Ntg(g)", "Oxg(g)", "H2Sg(g)")
  k <- gas_consts(db)
  ok <- vapply(intersect(comps, names(k)), function(g) all(c("-T_c", "-P_c", "-Omega") %in% k[[g]]), logical(1))
  check(all(ok), sprintf("%s: Peng-Robinson constants for %s", basename(db), paste(names(ok), collapse = " ")))
}

# 2. H2: Peng-Robinson against ideal.
od <- tempfile("gas_eos_"); dir.create(od)
pr <- "databases/coreclath_C4only.dat"
ideal <- file.path(od, "ideal.dat")
writeLines(sub("^(.*)  # PR constants restored 2026-10-28 \\(scripts/gas_eos.py\\)$", "#\\1", readLines(pr, warn = FALSE)), ideal)
inp <- file.path(od, "h2.pqi")
writeLines(c("SOLUTION 1", " temp 160", " pressure 261.2", " water 0.157",
             "GAS_PHASE 1", " -fixed_pressure", " -pressure 261.2", " -volume 1e-6", " -temperature 160", " H2(g) 0",
             "REACTION 1", " H2 1", " 2 mol",
             "SELECTED_OUTPUT", " -file h2.sel", " -high_precision true", " -molalities H2", " -gases H2(g)",
             "USER_PUNCH", " -headings phi_H2 Vm_gas", " 10 PUNCH PR_PHI(\"H2(g)\"), GAS_VM", "END"), inp)
res <- lapply(c(pr = pr, ideal = ideal), function(db) {
  r <- run_phreeqc(inp, db, workdir = od, warn_version = FALSE)
  x <- utils::read.delim(file.path(od, "h2.sel"), check.names = FALSE, strip.white = TRUE); names(x) <- trimws(names(x))
  c(errors = r$errors, m_H2 = x$m_H2[nrow(x)], Vm = x$Vm_gas[nrow(x)], phi = x$phi_H2[nrow(x)], V = x$volume[nrow(x)])
})
dm <- res$pr[["m_H2"]] / res$ideal[["m_H2"]] - 1; dv <- res$pr[["Vm"]] / res$ideal[["Vm"]] - 1
check(res$pr[["errors"]] == 0 && res$ideal[["errors"]] == 0 && res$pr[["phi"]] > 1 && abs(dm) < 0.15 && abs(dv) < 0.15,
      sprintf("H2, 160 C, 261.2 atm: dissolved at saturation PR %.4f vs ideal %.4f molal (%+.1f%%); bubble molar volume PR %.4f vs ideal %.4f L/mol (%+.1f%%); PR phi %.3f",
              res$pr[["m_H2"]], res$ideal[["m_H2"]], 100 * dm, res$pr[["Vm"]], res$ideal[["Vm"]], 100 * dv, res$pr[["phi"]]))
unlink(od, recursive = TRUE)
cat(sprintf("\nSUMMARY: %d failure(s).\n", failures))
if (failures > 0) quit(status = 1)
