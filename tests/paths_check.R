# Paths handed to PHREEQC (the macOS out-of-memory failure, 2026-11-05): run_phreeqc() must make every
# path absolute before it changes to the working directory, because normalizePath(mustWork = FALSE)
# leaves a path that doesn't exist yet relative on macOS and Linux. PHREEQC then could not open its
# output and re-prompted for a file name in a loop. Checks: a relative output with a different workdir;
# a missing output directory stops at once with a clear message; the version check; and run_case.R with
# a relative runs directory and case path. ~4 min. From the repo root:  Rscript tests/paths_check.R
suppressPackageStartupMessages({ library(tidyverse); library(glue); library(R6) })
invisible(capture.output(source("ENTICES.R")))
failures <- 0
check <- function(ok, what) { cat(if (isTRUE(ok)) "[PASS]" else "[FAIL]", what, "\n"); if (!isTRUE(ok)) failures <<- failures + 1 }
db <- normalizePath(file.path("databases", "coreclath_C4only.dat"))
root <- getwd()

# 1. relative input and output, run from a different working directory
base <- tempfile("paths_"); dir.create(file.path(base, "rel"), recursive = TRUE)
old <- setwd(base)
writeLines(c("SOLUTION 1", "SELECTED_OUTPUT", " -file b.csv", "END"), file.path("rel", "b.pqi"))
r <- run_phreeqc(file.path("rel", "b.pqi"), db, output = file.path("rel", "b.out"), workdir = "rel", warn_version = FALSE)
check(r$status == 0 && r$errors == 0 && file.exists(file.path("rel", "b.out")) && length(r$output_lines) > 10 &&
        file.exists(file.path("rel", "b.csv")),
      "relative input, output and workdir: PHREEQC writes the output where asked, and SELECTED_OUTPUT beside it")
r <- run_phreeqc(file.path("rel", "b.pqi"), db, output = file.path("rel", "c.out"), workdir = tempfile("wd_"), warn_version = FALSE)
check(r$status == 0 && file.exists(file.path("rel", "c.out")), "relative output with a different (temporary) workdir")

# 2. a missing output directory stops before PHREEQC starts
t0 <- Sys.time()
e <- tryCatch(run_phreeqc(file.path("rel", "b.pqi"), db, output = file.path("nodir", "b.out"), warn_version = FALSE),
              error = function(e) conditionMessage(e))
check(is.character(e) && grepl("output directory .* does not exist", e) && difftime(Sys.time(), t0, units = "secs") < 5,
      sprintf("missing output directory: stops at once (%s)", if (is.character(e)) e else "no error"))
setwd(old)

# 3. the version check run before chunk 1
v <- phreeqc_version_check()
check(isTRUE(v$ok) && identical(v$version, PHREEQC_EXPECTED_VERSION), sprintf("phreeqc_version_check(): PHREEQC %s", v$version))

# 4. run_case.R with a relative runs directory and a relative case path
runs <- file.path("runs", "_paths_test"); unlink(runs, recursive = TRUE)
st <- system2(file.path(R.home("bin"), "Rscript"), c("run_case.R", "--runs-dir", runs, file.path("tests", "cases", "short_reference.R")),
              stdout = TRUE, stderr = TRUE)
s <- file.path(runs, "short_reference", "summary.csv")
check(is.null(attr(st, "status")) && file.exists(s) && nrow(read.csv(s)) > 0,
      "run_case.R with a relative runs directory and case path finishes and writes summary.csv")
lg <- file.path(runs, "short_reference", "organic", "run.log")
check(file.exists(lg) && any(grepl(paste0("PHREEQC ", PHREEQC_EXPECTED_VERSION), readLines(lg))), "the run log records the PHREEQC version checked before chunk 1")
if (!nzchar(Sys.getenv("KEEP"))) unlink(runs, recursive = TRUE)
unlink(base, recursive = TRUE)
cat(sprintf("\nSUMMARY: %d failure(s).\n", failures))
if (failures > 0) quit(status = 1)
