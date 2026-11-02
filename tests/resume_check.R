# Resuming an interrupted run (production kit, item 2): run_case.R on a short case is killed mid-chunk (the
# whole process tree: the case, its organic and blank processes and PHREEQC), then rerun. The rerun must resume
# from the last complete chunk (the log says so), duplicate no chunk, and end with the same summary.csv as an
# uninterrupted run of the same case (1e-4 relative or 1e-6 absolute). ~10 min on a free machine. From the repo
# root:  Rscript tests/resume_check.R      (KEEP=1 keeps runs/_resume_test)
suppressPackageStartupMessages({ library(tidyverse); library(glue); library(R6) })
invisible(capture.output(source("ENTICES_LF_20260917.R")))
source("production_kit.R")
failures <- 0
check <- function(ok, what) { cat(if (isTRUE(ok)) "[PASS]" else "[FAIL]", what, "\n"); if (!isTRUE(ok)) failures <<- failures + 1 }
runs <- file.path("runs", "_resume_test"); unlink(runs, recursive = TRUE)
case <- "tests/cases/short_reference.R"; rscript <- file.path(R.home("bin"), "Rscript")
kill_tree <- function(pid) if (.Platform$OS.type == "windows") system2("taskkill", c("/F", "/T", "/PID", pid), stdout = FALSE, stderr = FALSE) else
  system(sprintf("pkill -KILL -P %d; kill -KILL %d", pid, pid))
# 1. uninterrupted
st <- system2(rscript, c("run_case.R", "--runs-dir", runs, "--name", "whole", case), stdout = TRUE, stderr = TRUE)
check(is.null(attr(st, "status")) && file.exists(file.path(runs, "whole", "summary.csv")), "the uninterrupted run finished")
# 2. interrupted: start in the background, kill once chunk 3 of the organic run is under way
d <- file.path(runs, "killed")
system2(rscript, c("run_case.R", "--runs-dir", runs, "--name", "killed", case), wait = FALSE, stdout = FALSE, stderr = FALSE)
t0 <- Sys.time(); killed_at <- NA
repeat {
  s <- if (file.exists(file.path(d, "organic", "state.rds"))) tryCatch(readRDS(file.path(d, "organic", "state.rds")), error = function(e) NULL) else NULL
  if (!is.null(s) && s$chunk >= 2 && dir.exists(file.path(d, "organic", sprintf("c%04d", s$chunk + 1)))) {
    Sys.sleep(2)
    pid <- as.integer(readLines(file.path(d, "lock"), n = 1)); kill_tree(pid)
    for (v in c("organic", "blank")) { lk <- file.path(d, v, "lock"); if (file.exists(lk)) kill_tree(as.integer(readLines(lk, n = 1))) }
    killed_at <- readRDS(file.path(d, "organic", "state.rds"))$chunk; break
  }
  if (difftime(Sys.time(), t0, units = "mins") > 30) break
  Sys.sleep(1)
}
Sys.sleep(3)
check(!is.na(killed_at) && variant_status(file.path(d, "organic")) == "stopped" && !file.exists(file.path(d, "summary.csv")),
      sprintf("the run was killed mid-chunk (organic state at chunk %s, chunk %s under way)", killed_at, killed_at + 1))
# 3. rerun the same command: resumes
st <- system2(rscript, c("run_case.R", "--runs-dir", runs, "--name", "killed", case), stdout = TRUE, stderr = TRUE)
check(is.null(attr(st, "status")) && file.exists(file.path(d, "summary.csv")), "the rerun finished")
lg <- readLines(file.path(d, "organic", "run.log"))
check(any(grepl(sprintf("resume after chunk %d", killed_at), lg)), sprintf("the organic log shows the resume after chunk %s", killed_at))
oc <- read.csv(file.path(d, "organic", "chunks.csv")); wc <- read.csv(file.path(runs, "whole", "organic", "chunks.csv"))
done_lines <- grep("^.{20}chunk [0-9]+: orbits", lg, value = TRUE)
nums <- as.integer(sub("^.{20}chunk ([0-9]+):.*", "\\1", done_lines))
check(identical(oc$chunk, seq_len(nrow(oc))) && nrow(oc) == nrow(wc) && !anyDuplicated(nums),
      sprintf("no chunk duplicated or lost (%d chunks, as uninterrupted; each logged once)", nrow(oc)))
a <- read.csv(file.path(runs, "whole", "summary.csv")); b <- read.csv(file.path(d, "summary.csv"))
num <- names(a)[vapply(a, is.numeric, TRUE)]; num <- setdiff(num, "seconds")
bad <- num[!vapply(num, function(cc) all(abs(a[[cc]] - b[[cc]]) <= pmax(1e-6, 1e-4 * abs(a[[cc]])) | (is.na(a[[cc]]) & is.na(b[[cc]]))), TRUE)]
check(nrow(a) == nrow(b) && !length(bad), sprintf("summary.csv equals the uninterrupted run's (%d columns compared)%s", length(num),
                                                  if (length(bad)) paste(": differ in", paste(bad, collapse = ", ")) else ""))
if (!nzchar(Sys.getenv("KEEP"))) unlink(runs, recursive = TRUE)
cat(sprintf("\nSUMMARY: %d failure(s).\n", failures))
if (failures > 0) quit(status = 1)
