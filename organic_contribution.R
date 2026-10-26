# organic_contribution.R -- the organics' contribution to ocean chemistry, by difference.
#
# pair_organic_run() generates an input from a PhreeqcIntegrator at its own organic_wt_percent
# and an otherwise identical one at organic_wt_percent = 0, runs both with run_phreeqc(), and
# returns organic_contribution() of the two CSVs.
#
# The blank has the same mineral inventory as the organic run, with the organics simply absent: it
# is built with mineral_basis_wt_percent = the run's organic_wt_percent, so its minerals are scaled
# by the same (100 - organic_wt_percent)/100 rather than rescaled up to a full rock.
#
# "Ocean" is the last SOLUTION 2 row for the tidal generators (the updated ocean of the final
# cycle) and the last SOLUTION 1 row for the bulk generators, which have a single solution.
# Quantities are molalities (mol/kgw) summed over the listed SELECTED_OUTPUT species. DIC counts
# only the C(4) species written by -molalities, so carbonate complexes of Ca, Mg, Fe, Na are not
# included; CH4 is m_CH4, or m_Mtg in decoupled runs.

ORGANIC_CONTRIBUTION_SPECIES <- list(
  DIC   = c("m_CO2", "m_HCO3-", "m_CO3-2", "m_NH4CO3-"),
  CH4   = c("m_CH4", "m_Mtg"),
  `NH4+` = "m_NH4+",
  H2    = "m_H2",
  C_total = "C"
)

# The ocean row of one ENTICES SELECTED_OUTPUT CSV, as a named numeric vector of the columns used.
ocean_row <- function(csv) {
  d <- suppressWarnings(read.delim(csv, sep = "\t", check.names = FALSE, strip.white = TRUE))
  names(d) <- trimws(names(d))
  ocean_soln <- if (any(d$soln == 2)) 2 else 1
  rows <- d[d$soln == ocean_soln & trimws(d$state) != "i_soln", , drop = FALSE]
  if (!nrow(rows)) stop("no reacted ocean row in ", csv)
  cols <- intersect(unlist(ORGANIC_CONTRIBUTION_SPECIES), names(d))
  vapply(cols, function(cl) as.numeric(rows[[cl]][nrow(rows)]), numeric(1))
}

# Organics' contribution by difference: run minus blank, per quantity, in mol/kgw.
organic_contribution <- function(run_csv, blank_csv) {
  run <- ocean_row(run_csv); blank <- ocean_row(blank_csv)
  total <- function(v, sp) { sp <- intersect(sp, names(v)); if (length(sp)) sum(v[sp]) else NA_real_ }
  out <- data.frame(
    quantity      = names(ORGANIC_CONTRIBUTION_SPECIES),
    with_organics = vapply(ORGANIC_CONTRIBUTION_SPECIES, function(sp) total(run, sp), numeric(1)),
    blank         = vapply(ORGANIC_CONTRIBUTION_SPECIES, function(sp) total(blank, sp), numeric(1)),
    row.names = NULL)
  out$organic_contribution <- out$with_organics - out$blank
  out
}

# Generates and runs the pair, writing each into its own subdirectory of output_dir ("run" and
# "org0"). `generator` is a generator name without the generate_phreeqc_input_ prefix, e.g.
# "tidal_kinetic"; `...` are its arguments other than output_dir. Run from the repo root so the
# stamped database path resolves.
pair_organic_run <- function(ph, generator, ..., output_dir) {
  if (!isTRUE(ph$organic_wt_percent > 0))
    stop("ph has organic_wt_percent = 0: there is nothing to pair it with")
  blank <- ph$clone(deep = TRUE)
  basis <- if (is.null(ph$mineral_basis_wt_percent)) ph$organic_wt_percent else ph$mineral_basis_wt_percent
  invisible(capture.output(blank$calculate_phreeqc_parameters(organic_wt_percent = 0, mineral_basis_wt_percent = basis)))
  run_one <- function(obj, sub) {
    f <- obj[[paste0("generate_phreeqc_input_", generator)]](output_dir = file.path(output_dir, sub), ...)
    db <- trimws(sub("^# Database:", "", grep("^# Database:", readLines(f, n = 20), value = TRUE)[1]))
    r <- run_phreeqc(f, db, output = sub("\\.pqi$", ".out", f), workdir = dirname(f))
    if (r$status != 0 || r$errors > 0)
      stop(sprintf("%s: PHREEQC status %s, %d errors", basename(f), r$status, r$errors))
    list(input = f, csv = file.path(dirname(f), paste0(tools::file_path_sans_ext(basename(f)), ".csv")),
         seconds = r$seconds)
  }
  run <- run_one(ph, "run"); org0 <- run_one(blank, "org0")
  list(contribution = organic_contribution(run$csv, org0$csv), run = run, blank = org0,
       organic_wt_percent = ph$organic_wt_percent, mineral_basis_wt_percent = basis)
}
