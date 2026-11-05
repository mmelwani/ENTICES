# Run the whole test suite (production kit, item 6). From the repo root:
#
#   Rscript tests/run_all.R            (everything, 20 tests: a few hours)
#   Rscript tests/run_all.R --quick    (skip the long ones: a few minutes)
#
# Each test runs in its own R process; its output goes to tests/output/run_all/<test>.log. At the end a table
# lists PASS or FAIL and the time per test, and the script exits with status 1 if any failed. Passing looks like
# "ALL n TESTS PASSED". If one fails, send its log (PRODUCTION_GUIDE.md, If something fails).
args <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args
# name, minutes (rough; measured on a shared machine 2026-11-02), long?
TESTS <- list(
  list("db_load_check.R", 1, FALSE), list("rock_table_check.R", 0.2, FALSE), list("iom_selftest.R", 0.5, FALSE),
  list("carbon_rebin_check.R", 0.2, FALSE), list("c4only_db_check.R", 0.5, FALSE), list("gas_eos_check.R", 0.5, FALSE),
  list("bulk_density_check.R", 0.5, FALSE), list("kozeny_carman_check.R", 1, FALSE), list("pressure_check.R", 1, FALSE),
  list("initial_fluid_check.R", 2, FALSE), list("cometary_ocean_check.R", 2, FALSE), list("exchange_check.R", 3, FALSE),
  list("porosity_evolution_check.R", 3, FALSE), list("liveness_check.R", 1, FALSE),
  list("restart_closure_check.R", 20, TRUE), list("gas_mass_balance_check.R", 20, TRUE),
  list("pipeline_check.R", 15, TRUE), list("resume_check.R", 20, TRUE), list("grouped_validation_check.R", 90, TRUE),
  list("paths_check.R", 30, TRUE))
if (!file.exists("ENTICES.R")) stop("run from the repository root: Rscript tests/run_all.R", call. = FALSE)
run <- Filter(function(t) !(quick && t[[3]]), TESTS)
logdir <- file.path("tests", "output", "run_all"); dir.create(logdir, recursive = TRUE, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), "Rscript")
cat(sprintf("Running %d tests%s (logs in %s)\n", length(run), if (quick) " (--quick: long tests skipped)" else "", logdir))
res <- do.call(rbind, lapply(run, function(t) {
  log <- file.path(logdir, sub("\\.R$", ".log", t[[1]]))
  cat(sprintf("  %-30s (~%g min) ... ", t[[1]], t[[2]])); flush.console()
  t0 <- Sys.time()
  st <- suppressWarnings(system2(rscript, file.path("tests", t[[1]]), stdout = log, stderr = log))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  out <- readLines(log, warn = FALSE)
  pass <- identical(as.integer(st), 0L) && any(grepl("SUMMARY: 0 failure", out))
  cat(if (pass) "PASS" else "FAIL", sprintf("(%.0f s)\n", secs))
  data.frame(test = t[[1]], result = if (pass) "PASS" else "FAIL", seconds = round(secs),
             summary = { s <- grep("SUMMARY", out, value = TRUE); if (length(s)) trimws(s[length(s)]) else sprintf("no SUMMARY line (exit status %s)", st) })
}))
cat("\n"); options(width = 250); print(res, row.names = FALSE)
nf <- sum(res$result == "FAIL")
cat(if (nf == 0) sprintf("\nALL %d TESTS PASSED\n", nrow(res)) else sprintf("\n%d OF %d TESTS FAILED: see the logs in %s\n", nf, nrow(res), logdir))
if (nf > 0) quit(status = 1)
