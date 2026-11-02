# The production pipeline end to end (production kit, item 4): run_case.R on two short cases (a few chunks each,
# organic and blank in parallel), then every figure. Checks the run directories, that summary.csv has the
# listed columns and summary_units.csv covers them, that the blank shares the organic run's times, and that
# every figure renders (PDF and PNG). ~10 min on a free machine. From the repo root:
#   Rscript tests/pipeline_check.R      (KEEP=1 keeps runs/_pipeline_test)
suppressPackageStartupMessages({ library(tidyverse); library(glue); library(R6) })
invisible(capture.output(source("ENTICES_LF_20260917.R")))
source("production_kit.R"); source("make_figures.R")
failures <- 0
check <- function(ok, what) { cat(if (isTRUE(ok)) "[PASS]" else "[FAIL]", what, "\n"); if (!isTRUE(ok)) failures <<- failures + 1 }
runs <- file.path("runs", "_pipeline_test"); unlink(runs, recursive = TRUE)
cases <- c("tests/cases/short_reference.R", "tests/cases/short_k_1e-11.R")
rscript <- file.path(R.home("bin"), "Rscript")
for (f in cases) {
  st <- system2(rscript, c("run_case.R", "--runs-dir", runs, f), stdout = TRUE, stderr = TRUE)
  nm <- validate_case(load_case(f))$name; d <- file.path(runs, nm)
  check(is.null(attr(st, "status")) && all(file.exists(file.path(d, c("summary.csv", "summary_units.csv", "RUNINFO.txt", "organic/run.log", "blank/run.log", "organic/state.rds", "blank/state.rds")))),
        sprintf("%s: run_case.R finished and wrote summary.csv, summary_units.csv, RUNINFO.txt, logs and state files", nm))
  if (!file.exists(file.path(d, "summary.csv"))) { cat(tail(st, 20), sep = "\n"); next }
  s <- read.csv(file.path(d, "summary.csv")); u <- read.csv(file.path(d, "summary_units.csv"))
  need <- c("time_yr", "orbit", paste0("ocean_", names(KIT_OCEAN_SPECIES)), paste0("org_minus_blank_ocean_", names(KIT_OCEAN_SPECIES)),
            "ocean_gas_H2_mol", "h2_generated_cum", "h2_present", "h2_generated_cum_core_mol", "phi", "k_m2", "delta_m", "skin_T_C", "skin_P_atm",
            "alpha", "K_d_Pa", "exchange_kg_per_orbit", "refill_kg_per_orbit", "rho_solid", "rho_bulk", "rho_bulk_in_bagheri", "clogged",
            "drift_max", "steady", "end_reason", "clog_time_yr")
  check(all(need %in% names(s)), sprintf("%s: summary.csv has the listed columns (%d rows, %d columns)%s", nm, nrow(s), ncol(s),
                                         if (all(need %in% names(s))) "" else paste(": missing", paste(setdiff(need, names(s)), collapse = ", "))))
  check(setequal(u$column, names(s)) && all(nzchar(u$description)), sprintf("%s: summary_units.csv describes every column", nm))
  oc <- read.csv(file.path(d, "organic", "chunks.csv")); bc <- read.csv(file.path(d, "blank", "chunks.csv"))
  check(nrow(oc) >= 3 && identical(oc$orbit_to, bc$orbit_to) && isTRUE(all.equal(oc$phi_start, bc$phi_start)) && isTRUE(all.equal(oc$k, bc$k)),
        sprintf("%s: the blank followed the organic run (%d chunks, same orbits, phi and k)", nm, nrow(oc)))
  ri <- readLines(file.path(d, "RUNINFO.txt"))
  check(any(grepl("^ENDED BY: t_end", ri)) && any(grepl("^code: ENTICES commit", ri)) && all(s$end_reason == "t_end"),
        sprintf("%s: RUNINFO and summary.csv state what ended the run (t_end) and the commit", nm))
  pq <- readLines(file.path(d, "organic", readRDS(file.path(d, "organic", "state.rds"))$files[1]), n = 12)
  check(any(grepl("^# Code: ENTICES commit [0-9a-f]{40}", pq)), sprintf("%s: the .pqi header carries the commit", nm))
}
S <- read_summaries(runs)
S$case[S$case == "short_reference"] <- "reference"   # the tornado figure needs a reference
out <- tempfile("figures_"); f <- tryCatch(make_all_figures(S, out = out), error = function(e) { cat(conditionMessage(e), "\n"); character(0) })
check(length(f) == 12 && all(file.exists(f)) && all(file.size(f) > 2000), sprintf("all six figures render as PDF and PNG (%d files)", length(f)))
ov <- data.frame(panel = c("pH", "DIC_C4"), ymin = c(8, NA), ymax = c(9, NA), x = c(NA, 0.05), y = c(NA, 1e-14), label = c("band", "point"))
g <- tryCatch(save_figure(fig_ocean(S, species = c("pH", "DIC_C4"), overlay = ov), "overlay", out), error = function(e) { cat(conditionMessage(e), "\n"); character(0) })
check(length(g) == 2 && all(file.exists(g)), "an overlay of bands and points renders")
if (!nzchar(Sys.getenv("KEEP"))) unlink(runs, recursive = TRUE)
unlink(out, recursive = TRUE)
cat(sprintf("\nSUMMARY: %d failure(s).\n", failures))
if (failures > 0) quit(status = 1)
