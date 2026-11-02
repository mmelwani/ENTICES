# The run watchdog and the staged-run status (PI 2026-11-02), without running PHREEQC: liveness_watch() gives up on a
# run only when its files are stale AND the process table says its process is gone, on consecutive checks over a
# minimum span; an unanswered process-table check counts as alive. variant_status() reopens only runs that ended at
# an earlier t_end. ~1 min. From the repo root:  Rscript tests/liveness_check.R
suppressPackageStartupMessages({ library(tidyverse); library(glue); library(R6) })
invisible(capture.output(source("ENTICES_LF_20260917.R")))
source("production_kit.R")
failures <- 0
check <- function(ok, what) { cat(if (isTRUE(ok)) "[PASS]" else "[FAIL]", what, "\n"); if (!isTRUE(ok)) failures <<- failures + 1 }
d <- tempfile("live_"); dir.create(file.path(d, "c0003"), recursive = TRUE)
files <- file.path(d, c("run.log", "state.rds", "c0003/run.out")); for (f in files) writeLines("x", f)
lock <- file.path(d, "lock")
dead_pid <- as.integer(system2(file.path(R.home("bin"), "Rscript"), c("-e", shQuote("cat(Sys.getpid())")), stdout = TRUE))
age_files <- function(s) Sys.setFileTime(c(files, lock), Sys.time() - s)

check(isTRUE(pid_state(Sys.getpid())) && identical(pid_state(dead_pid), FALSE) && pid_alive(Sys.getpid()) && !pid_alive(dead_pid),
      sprintf("pid_state: this process TRUE, a finished one (%d) FALSE", dead_pid))
writeLines(as.character(dead_pid), lock)
w <- liveness_watch(d, lock, stale_s = 600, n_checks = 3, span_s = 2)
r <- vapply(1:5, function(i) { Sys.sleep(0.5); w() }, TRUE)
check(!any(r), "fresh files, process gone: never given up on (the files are the primary signal)")
age_files(3600); w <- liveness_watch(d, lock, stale_s = 600, n_checks = 3, span_s = 2)
r <- vapply(1:6, function(i) { x <- w(); Sys.sleep(1); x }, TRUE)
check(!any(r[1:2]) && r[length(r)], sprintf("stale files, process gone: given up on only after 3 checks over 2 s (%s)", paste(r, collapse = " ")))
w <- liveness_watch(d, lock, stale_s = 600, n_checks = 3, span_s = 2)
r1 <- c(w(), w()); writeLines("y", files[3]); r2 <- w(); age_files(3600); Sys.sleep(1)
r3 <- vapply(1:2, function(i) { Sys.sleep(1); w() }, TRUE)
check(!any(c(r1, r2, r3)), "a write in the chunk directory resets the count")
writeLines(as.character(Sys.getpid()), lock); age_files(3600); w <- liveness_watch(d, lock, stale_s = 600, n_checks = 3, span_s = 2)
r <- vapply(1:5, function(i) { Sys.sleep(0.5); w() }, TRUE)
check(!any(r), "stale files, process alive (a long PHREEQC step): never given up on")
writeLines(as.character(dead_pid), lock); age_files(3600)
real <- pid_state; pid_state <- function(pid) NA                  # tasklist gives no answer (a loaded machine)
w <- liveness_watch(d, lock, stale_s = 600, n_checks = 3, span_s = 2)
r <- vapply(1:5, function(i) { Sys.sleep(0.5); w() }, TRUE); pid_state <- real
check(!any(r), "stale files, process table unanswered: counted as alive, never given up on")

# variant_status() and the notes, on a state that ended at t_end = 0.1 yr (time_yr 0.1016)
v <- file.path(d, "organic"); dir.create(v)
saveRDS(list(stopped = "t_end", time_yr = 0.1016), file.path(v, "state.rds"))
check(variant_status(v, 0.2) == "stopped" && variant_status(v, 0.1) == "done" && variant_status(v) == "done" && variant_status(v, 0.05) == "done",
      "ended at t_end: reopens when t_end is raised (0.2), done at 0.1 and below")
check(grepl("already at t = 0.1016 yr, past t_end_yr = 0.05: kept as it is", variant_done_note(v, 0.05)), paste("lowered:", variant_done_note(v, 0.05)))
for (why in c("steady state", "clogging", "phi_c")) {
  saveRDS(list(stopped = why, time_yr = 0.1016), file.path(v, "state.rds"))
  check(variant_status(v, 1e8) == "done" && grepl("won't reopen it", variant_done_note(v, 1e8)), sprintf("ended by %s: stays done at t_end 1e8 (%s)", why, variant_done_note(v, 1e8)))
}
unlink(d, recursive = TRUE)
cat(sprintf("\nSUMMARY: %d failure(s).\n", failures))
if (failures > 0) quit(status = 1)
