# Every database in databases/ must load in PHREEQC with 0 errors.
# (2e13e34 shipped a coreclath_CH4uncoupled.dat that failed to load; this
# would have caught it.)
#
# Run from the repo root:  Rscript tests/db_load_check.R

source("phreeqc_tools.R")

input <- tempfile(fileext = ".pqi")
writeLines(c("SOLUTION 1", " temp 25", "END"), input)

failures <- 0
dbs <- list.files("databases", pattern = "[.]dat$", full.names = TRUE)
if (length(dbs) == 0) { cat("[FAIL] no .dat files found in databases/\n"); failures <- 1 }

for (db in dbs) {
  r <- run_phreeqc(input, db, warn_version = FALSE)
  ok <- r$status == 0 && r$errors == 0 && length(r$output_lines) > 0
  cat(if (ok) "[PASS]" else "[FAIL]",
      sprintf("%-32s errors %d  warnings %d  (PHREEQC %s)", basename(db), r$errors, r$warnings, r$version), "\n")
  if (!ok) {
    failures <- failures + 1
    cat(paste("   ", head(grep("^ERROR", r$output_lines, value = TRUE), 5)), sep = "\n")
  }
  if (!identical(r$version, PHREEQC_EXPECTED_VERSION))
    cat(sprintf("   note: tested with PHREEQC %s, expected %s\n", r$version, PHREEQC_EXPECTED_VERSION))
}

cat(sprintf("\nSUMMARY: %d failure(s).\n", failures))
if (failures > 0) quit(status = 1)
