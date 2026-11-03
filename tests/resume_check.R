# Resuming an interrupted run (production kit, item 2): run_case.R on a short case is killed mid-chunk (the
# whole process tree: the case, its organic and blank processes and PHREEQC), then rerun. The rerun must resume
# from the last complete chunk (the log says so), duplicate no chunk, and end with the same summary.csv as an
# uninterrupted run of the same case (1e-4 relative or 1e-6 absolute). Then the staged runs of PRODUCTION_GUIDE.md:
# raising t_end extends a finished run (organic and blank in lockstep), lowering it keeps the run whole, and a run
# that ended at a steady state stays done. ~20 min on a free machine. From the repo
# root:  Rscript tests/resume_check.R      (KEEP=1 keeps runs/_resume_test)
suppressPackageStartupMessages({ library(tidyverse); library(glue); library(R6) })
invisible(capture.output(source("ENTICES.R")))
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
# 4. extend the finished uninterrupted run: the same case with t_end raised from 0.1 to 0.2 yr, rerun under the same
#    name (the staged runs of PRODUCTION_GUIDE.md). It must resume after the last chunk, keep the earlier chunks
#    unchanged, and end at the new t_end, for the organic run and its blank, with the same chunk boundaries.
ext <- file.path(runs, "short_reference.R")
writeLines('case <- list(inherits = "reference", name = "short_reference", t_end_yr = 0.2)', ext)
w <- file.path(runs, "whole"); before <- lapply(c(organic = "organic", blank = "blank"), function(v) read.csv(file.path(w, v, "chunks.csv")))
st <- system2(rscript, c("run_case.R", "--runs-dir", runs, "--name", "whole", ext), stdout = TRUE, stderr = TRUE)
check(is.null(attr(st, "status")) && file.exists(file.path(w, "summary.csv")), "raising t_end on the finished run: the rerun finished")
for (v in c("organic", "blank")) {
  n0 <- nrow(before[[v]]); after <- read.csv(file.path(w, v, "chunks.csv")); lg <- readLines(file.path(w, v, "run.log"))
  same <- isTRUE(all.equal(after[seq_len(n0), setdiff(names(after), "seconds")], before[[v]][, setdiff(names(after), "seconds")]))
  check(any(grepl(sprintf("resume after chunk %d", n0), lg)) && nrow(after) > n0 && same && identical(after$chunk, seq_len(nrow(after))) &&
        abs(tail(after$time_end_yr, 1) - 0.2) < 0.01,
        sprintf("%s: extended, not restarted (resumed after chunk %d, %d -> %d chunks, earlier chunks unchanged, t = %.4g -> %.4g yr)",
                v, n0, n0, nrow(after), tail(before[[v]]$time_end_yr, 1), tail(after$time_end_yr, 1)))
}
s <- read.csv(file.path(w, "summary.csv"))
check(all(s$end_reason == "t_end") && abs(max(s$time_yr) - 0.2) < 0.01, sprintf("the extended summary.csv reaches %.4g yr, ended by t_end", max(s$time_yr)))
oc <- read.csv(file.path(w, "organic", "chunks.csv")); bc <- read.csv(file.path(w, "blank", "chunks.csv"))
check(identical(oc$orbit_to, bc$orbit_to) && isTRUE(all.equal(oc$time_end_yr, bc$time_end_yr)) && isTRUE(all.equal(oc$phi_start, bc$phi_start)),
      sprintf("the blank extended in lockstep with the organic run (%d chunks, the same chunk boundaries and phi)", nrow(oc)))
# 5. lowering t_end below what the run has reached (0.2 -> 0.05 yr): skipped with a note, neither truncated nor rerun
low <- file.path(runs, "low", "short_reference.R"); dir.create(dirname(low), showWarnings = FALSE)
writeLines('case <- list(inherits = "reference", name = "short_reference", t_end_yr = 0.05)', low)
st <- system2(rscript, c("run_case.R", "--runs-dir", runs, "--name", "whole", low), stdout = TRUE, stderr = TRUE)
oc2 <- read.csv(file.path(w, "organic", "chunks.csv")); s2 <- read.csv(file.path(w, "summary.csv"))
check(is.null(attr(st, "status")) && sum(grepl("already at t = .* past t_end_yr = 0.05: kept as it is", st)) == 2 &&
      isTRUE(all.equal(oc2[, names(oc2) != "seconds"], oc[, names(oc) != "seconds"])) && abs(max(s2$time_yr) - max(oc$time_end_yr)) < 1e-9,
      "lowering t_end: both variants skipped as already past it; the chunks and summary.csv are kept whole")
# 6. a run that ended at a steady state stays done when t_end is raised (the state of the resumed run, relabelled)
k <- file.path(runs, "killed")
for (v in c("organic", "blank")) { f <- file.path(k, v, "state.rds"); x <- readRDS(f); x$stopped <- "steady state"; saveRDS(x, f) }
n_k <- nrow(read.csv(file.path(k, "organic", "chunks.csv")))
st <- system2(rscript, c("run_case.R", "--runs-dir", runs, "--name", "killed", ext), stdout = TRUE, stderr = TRUE)
check(is.null(attr(st, "status")) && sum(grepl("ended by steady state .*raising t_end_yr won't reopen it", st)) == 2 &&
      nrow(read.csv(file.path(k, "organic", "chunks.csv"))) == n_k,
      sprintf("a run ended at a steady state stays done when t_end is raised (skipped with a note; %d chunks)", n_k))
if (!nzchar(Sys.getenv("KEEP"))) unlink(runs, recursive = TRUE)
cat(sprintf("\nSUMMARY: %d failure(s).\n", failures))
if (failures > 0) quit(status = 1)
