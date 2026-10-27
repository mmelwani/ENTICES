# databases/coreclath_C4only.dat (carbon_mode "c4_only") and restrict_to_database().
#
#  1. The header documents exactly the entries the database comments out.
#  2. Every carbon species it still defines is C(+4) or decoupled Mtg.
#  3. A reducing carbonate solution (pe -12, 25 and 100 C) speciates as C(+4)
#     only, and its carbon all stays in solution.
#  4. restrict_to_database() removes the reduced-carbon names from a generated
#     input (tests/N1_leak_T100.pqi) and from a small input that then runs.
#
# Run from the repo root:  Rscript tests/c4only_db_check.R

source("phreeqc_tools.R")

failures <- 0
check <- function(ok, what) {
  cat(if (isTRUE(ok)) "[PASS]" else "[FAIL]", what, "\n")
  if (!isTRUE(ok)) failures <<- failures + 1
}
db  <- "databases/coreclath_C4only.dat"
src <- "databases/coreclath_CH4uncoupled.dat"

# 1. Header vs. what is actually gone.
u <- database_names(src); c4 <- database_names(db)
gone <- c(setdiff(u$phases, c4$phases), setdiff(u$species, c4$species))
hdr <- grep("^#   (master|species|phase) ", readLines(db, n = 200), value = TRUE)
hdr_names <- sub("^#   \\S+\\s+(\\S+).*", "\\1", hdr)
hdr_kind  <- sub("^#   (\\S+).*", "\\1", hdr)
check(setequal(hdr_names[hdr_kind != "master"], gone),
      sprintf("header lists the %d species/phases removed (%d species, %d phases)", length(gone),
              length(setdiff(u$species, c4$species)), length(setdiff(u$phases, c4$phases))))
check(setequal(hdr_names[hdr_kind == "master"], c("C(-3)", "C(-2)", "C(+2)")), "master species C(-3), C(-2), C(+2) removed")
check(length(setdiff(c4$phases, u$phases)) + length(setdiff(c4$species, u$species)) == 0, "nothing added")
check(all(c("Mtg", "HCO3-", "CO3-2", "CO2") %in% c4$species) && all(c("Mtg(g)", "Mtg_hydrate", "CO2(g)", "CO2_hydrate", "Calcite") %in% c4$phases),
      "Mtg, Mtg(g), Mtg_hydrate and the carbonate system are kept")

# 2. Remaining carbon species: classified, and C(+4) unless Mtg.
c_sp <- c4$species[grepl("C(?![a-z])", c4$species, perl = TRUE) | c4$species == "Mtg"]
m <- lapply(c_sp, carbon_moiety)
check(!any(vapply(m, is.null, logical(1))), sprintf("all %d remaining carbon species classified", length(c_sp)))
ox <- vapply(m, function(x) if (is.null(x)) NA_real_ else x$ox, 0)
check(all(ox[c_sp != "Mtg"] == 4), sprintf("all non-Mtg carbon species are C(+4): %s", paste(c_sp[c_sp != "Mtg"], collapse = " ")))

# 3. Reducing carbonate solution.
for (T in c(25, 100)) {
  inp <- tempfile(fileext = ".pqi"); out <- tempfile(fileext = ".out")
  writeLines(c("SOLUTION 1", sprintf(" temp %d", T), " pH 10", " pe -12", " units mol/kgw",
               " C 1e-2", " Na 1e-2 charge", " N 1e-3", " S 1e-3", " Mtg 1e-4", "END"), inp)
  r <- run_phreeqc(inp, db, output = out, warn_version = FALSE)
  st <- carbon_by_oxidation_state(out)$by_state
  cpl <- st[!st$decoupled, ]
  check(r$errors == 0 && all(cpl$ox == 4) && abs(sum(cpl$C_molal) - 1e-2) <= 5e-4 * 1e-2 &&
          abs(sum(st$C_molal[st$decoupled]) - 1e-4) <= 5e-4 * 1e-4,
        sprintf("%3d C, pe -12: coupled C %.5e molal, all C(+4); Mtg %.3e", T, sum(cpl$C_molal), sum(st$C_molal[st$decoupled])))
}

# 4. restrict_to_database() on generated and synthetic inputs.
n1 <- readLines("tests/N1_leak_T100.pqi")
r <- restrict_to_database(n1, db, base = src)
expect <- c("C2H4(g)", "C2H6(g)", "C3H8(g)", "Citric_Acid", "CO(g)", "Glycine", "C2H4", "C2H6", "CO", "CH3COO-", "HCOO-", "HCN")
check(setequal(attr(r, "dropped"), expect), sprintf("N1 input: removed %s", paste(attr(r, "dropped"), collapse = " ")))
check(all(attr(r, "dropped") %in% gone), "N1 input: everything removed is a removed database entry")
check(length(attr(r, "dropped_nonzero")) == 0, "N1 input: no removed EP entry had a nonzero amount")
rl <- strsplit(r, "\n", fixed = TRUE)[[1]]
check(sum(grepl("#", rl, fixed = TRUE)) == sum(grepl("#", n1, fixed = TRUE)), "N1 input: every comment kept")
check(identical(restrict_to_database(n1, src)[1], paste(n1, collapse = "\n")), "N1 input: unchanged against its own database")
check(setequal(attr(restrict_to_database(n1, db), "dropped"), expect), "N1 input: same removals without base")

syn <- c("SOLUTION 1", " temp 25", " pH 9", " pe -8", " C 1e-3", " Na 1e-3 charge",
         "EQUILIBRIUM_PHASES 1", "   Calcite 0 0", "   C2H6(g) 0 0.5", "   -force_equality true", "   Glycine 0 0", "#   Alanine 0 0",
         "GAS_PHASE 1", "   -fixed_volume", "   CO2(g) 0", "   CO(g) 0", "   Mtg(g) 0",
         "SELECTED_OUTPUT 1", "   -file sel.csv", "   -molalities HCO3- C2H6 CO3-2 #CH3COO-", "               HCOO- CH3OH",
         "   -saturation_indices Calcite C2H6(g)", "   -gases CO2(g) CO(g) Mtg(g)", "END")
rs <- restrict_to_database(syn, db, base = src)
rsl <- strsplit(rs, "\n", fixed = TRUE)[[1]]
check(identical(rsl[grep("-molalities", rsl)], "   -molalities HCO3- CO3-2 #CH3COO-") && !any(grepl("^\\s+HCOO- CH3OH", rsl)),
      "SELECTED_OUTPUT list edited, comment kept, emptied continuation line removed")
check(!any(grepl("C2H6\\(g\\)|force_equality|Glycine", rsl)) && any(grepl("#   Alanine", rsl)),
      "EP entry removed with its option line; commented EP line kept")
check(identical(attr(rs, "dropped_nonzero"), "C2H6(g) 0.5"), "EP entry with a nonzero amount reported")
typo <- c("SELECTED_OUTPUT 1", "   -saturation_indices Calcite C2H6(g) Calcte")
check(identical(attr(restrict_to_database(typo, db, base = src), "dropped"), "C2H6(g)"), "with base, a misspelt name is left for PHREEQC to reject")
check(setequal(attr(restrict_to_database(typo, db), "dropped"), c("C2H6(g)", "Calcte")), "without base, every name the database lacks is removed")
inp <- tempfile(fileext = ".pqi"); writeLines(rs, inp)
check(run_phreeqc(inp, db, warn_version = FALSE)$errors == 0, "restricted synthetic input runs with 0 errors")
inp0 <- tempfile(fileext = ".pqi"); writeLines(syn, inp0)
check(run_phreeqc(inp0, db, warn_version = FALSE)$errors > 0, "unrestricted synthetic input fails with the C4only database")

# 5. scripts/make_c4only.py regenerates the committed database exactly (CR characters
#    ignored: git may check the file out with either line ending).
rd_lf <- function(f) gsub("\r", "", rawToChar(readBin(f, "raw", file.size(f))), fixed = TRUE)
regen <- NA
for (py in unique(Sys.which(c("python3", "python", "py")))) {
  if (!nzchar(py)) next
  tmp <- tempfile(fileext = ".dat")
  suppressWarnings(system2(py, c("scripts/make_c4only.py", src, tmp), stdout = FALSE, stderr = FALSE))
  if (file.exists(tmp)) { regen <- tmp; break }
}
if (is.na(regen)) {
  cat("[SKIP] no working python on the PATH: regeneration of coreclath_C4only.dat not checked\n")
} else {
  check(identical(rd_lf(regen), rd_lf(db)), "scripts/make_c4only.py regenerates databases/coreclath_C4only.dat exactly")
}

cat(sprintf("\nSUMMARY: %d failure(s).\n", failures))
if (failures > 0) quit(status = 1)
