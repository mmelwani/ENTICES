# Run one case and its same-rock blank (production kit; PRODUCTION_GUIDE.md). From the repo root:
#
#   Rscript run_case.R cases/reference.R
#
# Writes runs/<name>/organic/, runs/<name>/blank/, runs/<name>/summary.csv, summary_units.csv and RUNINFO.txt.
# Rerunning the same command resumes an interrupted run from its last complete chunk.
#
# Options: --sequential (run the organic case, then the blank, in this process); --max-chunks N (stop after N
# chunks: for tests); --runs-dir DIR (default runs); --name NAME (override the case's name). Used internally:
# --variant organic|blank (run one variant in this process).
args <- commandArgs(trailingOnly = TRUE)
opt <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i)) default else args[i + 1] }
flags <- c("--max-chunks", "--runs-dir", "--name", "--variant")
pos <- args[!(args %in% c(flags, "--sequential")) & !(seq_along(args) - 1) %in% match(flags, args)]
if (length(pos) != 1) stop("usage: Rscript run_case.R [--sequential] [--max-chunks N] [--runs-dir DIR] [--name NAME] cases/<case>.R", call. = FALSE)
suppressPackageStartupMessages({ library(tidyverse); library(glue); library(R6) })
invisible(capture.output(source("ENTICES.R")))
source("production_kit.R")
max_chunks <- as.numeric(opt("--max-chunks", Inf)); runs_dir <- opt("--runs-dir", "runs")
variant <- opt("--variant")
if (!is.null(variant)) {
  case <- load_case(pos); if (!is.null(opt("--name"))) case$name <- opt("--name")
  case <- validate_case(case)
  run_variant(case, variant, file.path(runs_dir, case$name), max_chunks)
} else {
  t0 <- Sys.time()
  s <- run_case(pos, parallel = if ("--sequential" %in% args) FALSE else NA, max_chunks = max_chunks, runs_dir = runs_dir, name = opt("--name"))
  cat(sprintf("%s: %s at t = %.4g yr after %d chunks (%.0f s). Outputs in %s\n", s$case[1], s$end_reason[1], max(s$time_yr), nrow(s),
              as.numeric(difftime(Sys.time(), t0, units = "secs")), file.path(runs_dir, s$case[1])))
}
