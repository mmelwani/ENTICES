# Asserts that ENTICES's hard-coded ROCK_MOL_PER_KG matches the reconstruction's
# machine-readable output, reconstruction/rock_mol_per_kg.csv.
#
# Run from the repo root:  Rscript tests/rock_table_check.R
#
# Tolerance: a phase passes if it is within 1e-4 relative OR 1e-6 mol/kg absolute.
# ENTICES-side modelling steps applied before comparing:
#   - Fe metal: ENTICES takes lawrencite's Fe out of the metal, so the comparison
#     is ROCK_MOL_PER_KG$fe + ROCK_MOL_PER_KG$lawrencite vs the CSV's Iron.
#   - Lawrencite is ENTICES's own addition (not in the reconstruction): not checked.
#   - Inert phases (role == "inert", i.e. magnesiochromite) must NOT appear in
#     ROCK_MOL_PER_KG; they stay in the mass basis only.

suppressPackageStartupMessages({ library(tidyverse); library(glue); library(R6) })
source("ENTICES_LF_20260917.R")

REL_TOL <- 1e-4
ABS_TOL <- 1e-6

csv <- read.csv("reconstruction/rock_mol_per_kg.csv", comment.char = "#",
                stringsAsFactors = FALSE)

# reconstruction phase name -> ROCK_MOL_PER_KG key
key_map <- c(Forsterite = "forsterite", Enstatite = "enstatite", Diopside = "diopside",
             Tephroite = "tephroite", Anorthite = "anorthite", Albite = "albite",
             `K-Feldspar` = "kfeldspar", Troilite = "troilite", Schreibersite = "schreibersite",
             Iron = "fe", Nickel = "ni")

failures <- 0
check <- function(ok, msg) {
  cat(if (ok) "[PASS]" else "[FAIL]", msg, "\n")
  if (!ok) failures <<- failures + 1
}

for (i in seq_len(nrow(csv))) {
  ph <- csv$phase[i]
  if (csv$role[i] == "inert") {
    check(!(tolower(ph) %in% names(ROCK_MOL_PER_KG)),
          sprintf("%-16s inert: absent from ROCK_MOL_PER_KG", ph))
    next
  }
  key <- key_map[[ph]]
  if (is.null(key) || !(key %in% names(ROCK_MOL_PER_KG))) {
    check(FALSE, sprintf("%-16s no ROCK_MOL_PER_KG entry", ph)); next
  }
  entices <- ROCK_MOL_PER_KG[[key]]
  if (key == "fe") entices <- entices + ROCK_MOL_PER_KG[["lawrencite"]]
  recon <- csv$mol_per_kg[i]
  abs_d <- abs(entices - recon); rel_d <- abs_d / recon
  check(rel_d <= REL_TOL || abs_d <= ABS_TOL,
        sprintf("%-16s ENTICES %.6f  reconstruction %.6f  rel %.1e  abs %.1e",
                ph, entices, recon, rel_d, abs_d))
}

extra <- setdiff(names(ROCK_MOL_PER_KG), c(unname(key_map), "lawrencite"))
check(length(extra) == 0,
      sprintf("no ROCK_MOL_PER_KG entries without a reconstruction row%s",
              if (length(extra)) paste0(" (extra: ", paste(extra, collapse = ", "), ")") else ""))

cat(sprintf("\nSUMMARY: %d failure(s).\n", failures))
if (failures > 0) quit(status = 1)
