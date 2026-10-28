# bulk_density() / solid_density() (phreeqc_tools.R) and data/mineral_molar_volumes.csv.
#
#  1. Every phase the generators list (production = "yes") has a molar volume.
#  2. Mass and volume close on a hand-built case: 1 mol forsterite + 2 mol quartz + 10 g IOM.
#  3. t = 0 against the reconstruction: the Scenario R rock (reconstruction/rock_mol_per_kg.csv)
#     through solid_density(), against the density the reconstruction's own phase library implies
#     (reconstruction/Primordial_CI_chondrite_mineralogy.Rmd, section 4, phase_lib Vm). Neither the
#     Rmd nor RECONSTRUCTION_STATUS.md states a single rock density, so this is the comparison.
#     The two differ only through Schreibersite's Vm (this table: 27.56 cm3/mol from mindat unit-cell
#     data; phase_lib: M/7.44 g/cm3 from the Handbook of Mineralogy) and the molar masses (the
#     database's element weights vs the reconstruction's); the test requires agreement within 0.1%.
#
# Run from the repo root:  Rscript tests/bulk_density_check.R

source("phreeqc_tools.R")

failures <- 0
check <- function(ok, what) {
  cat(if (isTRUE(ok)) "[PASS]" else "[FAIL]", what, "\n")
  if (!isTRUE(ok)) failures <<- failures + 1
}
vt <- utils::read.csv(MINERAL_VOLUMES_CSV, stringsAsFactors = FALSE)
v <- read_mineral_volumes()

# 1. Production phases.
prod <- vt$phase[vt$production == "yes"]
check(all(!is.na(vt$Vm_cm3_mol[vt$production == "yes"])),
      sprintf("all %d production phases have a molar volume", length(prod)))

# 2. Hand-built case.
s <- solid_density(c(Forsterite = 1, Quartz = 2), iom_g = 10, iom_density = 1.4, volumes = v)
fo <- v[v$phase == "Forsterite", ]; q <- v[v$phase == "Quartz", ]
m <- fo$molar_mass_g_mol + 2 * q$molar_mass_g_mol + 10
vol <- fo$Vm_cm3_mol + 2 * q$Vm_cm3_mol + 10 / 1.4
check(abs(s$mass_g - m) < 1e-9 * m && abs(s$volume_cm3 - vol) < 1e-9 * vol && abs(s$rho_solid - m / vol) < 1e-12 &&
        abs(sum(s$breakdown$mass_frac) - 1) < 1e-12 && abs(sum(s$breakdown$vol_frac) - 1) < 1e-12,
      sprintf("forsterite + 2 quartz + 10 g IOM: %.4f g / %.4f cm3 = %.5f g/cm3; fractions sum to 1", s$mass_g, s$volume_cm3, s$rho_solid))
e <- tryCatch({ solid_density(c(Unobtainium = 1), volumes = v); "no error" }, error = function(e) conditionMessage(e))
check(grepl("no molar volume", e), "a phase with no molar volume is an error, not dropped")

# 3. t = 0 against the reconstruction's phase library.
rock <- utils::read.csv("reconstruction/rock_mol_per_kg.csv", comment.char = "#", stringsAsFactors = FALSE)
rmd <- paste(readLines("reconstruction/Primordial_CI_chondrite_mineralogy.Rmd"), collapse = "\n")
lib <- sub(".*phase_lib <- data.frame\\((.*?)\\n  Vm_src.*", "\\1", rmd)
ph  <- eval(parse(text = sub(".*phase = (c\\(.*?\\)),\\s*\\n\\s*formula.*", "\\1", lib)))
vm_txt <- sub(".*Vm = (c\\(.*?formula_mass\\(\"Fe3P\"\\) / 7.44\\)).*", "\\1", lib)
vm  <- eval(parse(text = sub('formula_mass("Fe3P")', as.character(rock$mw_g_mol[rock$phase == "Schreibersite"]), vm_txt, fixed = TRUE)))
lib_vm <- stats::setNames(vm, ph)
name <- c(Iron = "Fe", Nickel = "Ni")
mol <- stats::setNames(rock$mol_per_kg, ifelse(rock$phase %in% names(name), name[rock$phase], rock$phase))
ours <- solid_density(mol, volumes = v)$rho_solid
theirs <- sum(rock$mol_per_kg * rock$mw_g_mol) / sum(rock$mol_per_kg * lib_vm[rock$phase])
check(abs(ours / theirs - 1) < 1e-3,
      sprintf("Scenario R rock: %.4f g/cm3 (this table) vs %.4f g/cm3 (reconstruction phase_lib), %+.3f%%",
              ours, theirs, 100 * (ours / theirs - 1)))

cat(sprintf("\nSUMMARY: %d failure(s).\n", failures))
if (failures > 0) quit(status = 1)
