# Run every case in cases/sweep/ (production kit; PRODUCTION_GUIDE.md), resuming any that are unfinished.
# From the repo root:
#
#   Rscript run_sweep.R                 (2 cases at a time; each case runs its organic and blank in parallel)
#   Rscript run_sweep.R --parallel 4    (4 cases at a time: 8 PHREEQC processes)
#
# Options: --parallel N (cases at a time, default 2); --dir DIR (default cases/sweep); --runs-dir DIR (default
# runs). Each case's own log is runs/<name>/organic/run.log (and blank/); this script writes runs/sweep.log.
args <- commandArgs(trailingOnly = TRUE)
opt <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i)) default else args[i + 1] }
n_par <- as.integer(opt("--parallel", 2)); case_dir <- opt("--dir", "cases/sweep"); runs_dir <- opt("--runs-dir", "runs")
suppressPackageStartupMessages({ library(tidyverse); library(glue); library(R6) })
invisible(capture.output(source("ENTICES_LF_20260917.R")))
source("production_kit.R")
files <- sort(list.files(case_dir, "\\.R$", full.names = TRUE))
if (!length(files)) stop("no case files in ", case_dir, call. = FALSE)
for (f in files) validate_case(load_case(f))                   # every case valid before anything runs
dir.create(runs_dir, showWarnings = FALSE)
log <- file.path(runs_dir, "sweep.log")
kit_log(log, "sweep of %d cases in %s, %d at a time; %s", length(files), case_dir, n_par, code_version())
rscript <- file.path(R.home("bin"), "Rscript")
done_case <- function(f) { nm <- validate_case(load_case(f))$name; file.exists(file.path(runs_dir, nm, "summary.csv")) &&
  all(vapply(c("organic", "blank"), function(v) !dir.exists(file.path(runs_dir, nm, v)) || variant_status(file.path(runs_dir, nm, v)) == "done", TRUE)) }
todo <- files[!vapply(files, done_case, TRUE)]
kit_log(log, "%d finished already, %d to run", length(files) - length(todo), length(todo))
running <- list()
while (length(todo) || length(running)) {
  while (length(todo) && length(running) < n_par) {
    f <- todo[1]; todo <- todo[-1]; nm <- validate_case(load_case(f))$name
    dir.create(file.path(runs_dir, nm), showWarnings = FALSE, recursive = TRUE)
    out <- file.path(runs_dir, nm, "run_case.console.log")
    system2(rscript, c("run_case.R", "--runs-dir", shQuote(runs_dir), shQuote(f)), wait = FALSE, stdout = out, stderr = out)
    running[[f]] <- list(name = nm, started = Sys.time())
    kit_log(log, "started %s (%s)", nm, f)
    Sys.sleep(10)
  }
  Sys.sleep(30)
  for (f in names(running)) {
    nm <- running[[f]]$name; lk <- file.path(runs_dir, nm, "lock")
    alive <- file.exists(lk) && pid_alive(suppressWarnings(as.integer(readLines(lk, n = 1, warn = FALSE))))
    if (!alive && difftime(Sys.time(), running[[f]]$started, units = "secs") > 120) {
      ok <- done_case(f)
      kit_log(log, "%s %s%s", if (ok) "finished" else "STOPPED WITHOUT FINISHING:", nm,
              if (ok) "" else sprintf(" (see %s)", file.path(runs_dir, nm, "run_case.console.log")))
      running[[f]] <- NULL
    }
  }
}
kit_log(log, "sweep done")
cat(readLines(log), sep = "\n")
