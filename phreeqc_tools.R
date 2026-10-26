# phreeqc_tools.R -- locate and run the PHREEQC batch executable from R.
# Also carbon_by_oxidation_state(): carbon in PHREEQC output re-binned by the
# actual oxidation state of C in each species (see the section at the end).
#
# ENTICES targets PHREEQC 3.8.6 (PHREEQC_EXPECTED_VERSION). Results from other
# versions can differ in minor phases and trace species (see
# docs/reviews/UPDATES_CLEANUP_AND_RUNS_for_review_2026-10-22.md, section 3.3),
# so run_phreeqc() warns when the executable reports a different version.
#
# Executable lookup order: the PHREEQC_EXE environment variable, the default
# Windows install path, then `phreeqc` on the PATH.

PHREEQC_EXPECTED_VERSION <- "3.8.6"

phreeqc_exe <- function() {
  env <- Sys.getenv("PHREEQC_EXE")
  if (nzchar(env)) return(env)
  win <- "C:/Program Files/USGS/phreeqc-3.8.6-17100-x64/bin/Release/phreeqc.exe"
  if (file.exists(win)) return(win)
  on_path <- Sys.which("phreeqc")
  if (nzchar(on_path)) return(unname(on_path))
  stop("PHREEQC executable not found: set the PHREEQC_EXE environment variable")
}

# Runs PHREEQC on `input` with `database`, inside a temporary working
# directory so phreeqc.log and SELECTED_OUTPUT files don't land in the repo.
# Returns list(status, version, output_lines, errors, warnings, seconds).
run_phreeqc <- function(input, database, output = tempfile(fileext = ".out"),
                        workdir = tempfile("phreeqc_run_"), warn_version = TRUE) {
  dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
  input <- normalizePath(input, mustWork = TRUE)
  database <- normalizePath(database, mustWork = TRUE)
  output <- normalizePath(output, mustWork = FALSE)
  old <- setwd(workdir); on.exit(setwd(old), add = TRUE)
  t0 <- Sys.time()
  console <- suppressWarnings(system2(phreeqc_exe(), c(shQuote(input), shQuote(output), shQuote(database)),
                                      stdout = TRUE, stderr = TRUE))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  status <- attr(console, "status"); if (is.null(status)) status <- 0L
  ver <- regmatches(console, regexpr("PHREEQC-[0-9.]+[0-9]", console))
  ver <- if (length(ver)) sub("PHREEQC-", "", ver[1], fixed = TRUE) else NA_character_
  # The version an input was generated for: the '# Target PHREEQC version:' line
  # that stamp_pqi_header() writes, else PHREEQC_EXPECTED_VERSION.
  head_lines <- readLines(input, n = 20, warn = FALSE)
  stamped <- grep("^# Target PHREEQC version: ", head_lines, value = TRUE)
  target <- if (length(stamped)) trimws(sub("# Target PHREEQC version: ", "", stamped[1], fixed = TRUE))
            else PHREEQC_EXPECTED_VERSION
  if (warn_version && !identical(ver, target))
    warning(sprintf("PHREEQC version %s differs from %s, the version %s was generated for",
                    ver, target, basename(input)))
  out <- if (file.exists(output)) readLines(output, warn = FALSE) else character(0)
  list(status = status, version = ver, output_lines = out,
       errors = sum(startsWith(out, "ERROR")), warnings = sum(startsWith(out, "WARNING")),
       seconds = secs)
}

# ---- Carbon by actual oxidation state -------------------------------------
#
# PHREEQC groups each species under the master species its reaction is written
# from, and in the coreclath databases formate, acetate, methanol, propane and
# the other organics are written from HCO3- + O2, so PHREEQC's "C(4)" total
# includes them whatever their redox state. carbon_by_oxidation_state() re-bins
# carbon by the oxidation state of C in each species (mean over its C atoms),
# using CARBON_MOIETIES. Decoupled methane (Mtg) is counted as C(-4), flagged
# decoupled = TRUE.

# Carbon moiety: regex on the species name, C atoms, mean C oxidation state
# (H +1, O -2, N -3, S -2). First match wins; a "(X)n" ligand multiplies.
CARBON_MOIETIES <- data.frame(
  pattern = c("^Mtg$", "C4H8N2O3", "C3H7NO2", "CH3COCOO", "C3H5O2|C3H6O2", "C3H8",
              "C2H6NO2|C2H5NO2|C2H4NO2", "C2H5NO$", "CH3CH2OH", "CH3COO", "C2H6", "C2H4",
              "CH3NH2|CH3NH3", "CH3OH", "CH2O$", "HCOO", "SCN", "CN", "CO\\(NH2\\)2", "^CO$",
              "CO3|HCO3|^CO2$"),
  moiety  = c("Mtg", "glycylglycine", "alanine", "pyruvate", "propanoate", "propane",
              "glycine", "acetamide", "ethanol", "acetate", "ethane", "ethene",
              "methylamine", "methanol", "formaldehyde", "formate", "thiocyanate", "cyanide", "urea", "CO",
              "carbonate"),
  nC      = c(1, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1),
  ox      = c(-4, 1, 0, 2/3, -2/3, -8/3, 1, 0, -2, 0, -3, -2, -2, -2, 0, 2, 4, 2, 4, 2, 4),
  stringsAsFactors = FALSE)

carbon_moiety <- function(species) {
  hit <- which(vapply(CARBON_MOIETIES$pattern, function(p) grepl(p, species), logical(1)))[1]
  if (is.na(hit)) return(NULL)
  m <- CARBON_MOIETIES[hit, ]
  k <- regmatches(species, regexec("\\((CH3COO|HCOO|CN|CO3|HCO3)\\)([0-9]+)", species))[[1]]
  if (length(k)) m$nC <- m$nC * as.numeric(k[3])
  m
}

# file: a PHREEQC output file (every "Distribution of species" block, all carbon
# species) or a SELECTED_OUTPUT file (every row, but only the m_ species it
# punches, plus the Mtg total if punched). Returns list(by_state, by_species),
# keyed by block (output: the solution or mix it used; selected output: row,
# state, soln). C_molal = nC * molality.
carbon_by_oxidation_state <- function(file) {
  lines <- readLines(file, warn = FALSE)
  rows <- if (any(grepl("Distribution of species", lines, fixed = TRUE))) .carbon_from_output(lines)
          else .carbon_from_selected_output(file)
  unk <- unique(rows$species[is.na(rows$ox)])
  if (length(unk)) warning("carbon_by_oxidation_state: unclassified carbon species: ", paste(unk, collapse = ", "))
  rows$C_molal <- rows$nC * rows$molality
  by_state <- aggregate(C_molal ~ block + context + ox + decoupled, data = rows[!is.na(rows$ox), ], FUN = sum)
  by_state <- by_state[order(by_state$block, by_state$ox), ]
  rownames(by_state) <- NULL
  list(by_state = by_state, by_species = rows)
}

.carbon_rows <- function(block, context, species, molality) {
  ms <- lapply(species, carbon_moiety)
  data.frame(block = rep(block, length(species)), context = rep(context, length(species)),
             species = species, molality = molality,
             moiety = vapply(ms, function(m) if (is.null(m)) NA_character_ else m$moiety, ""),
             nC = vapply(ms, function(m) if (is.null(m)) NA_real_ else m$nC, 0),
             ox = vapply(ms, function(m) if (is.null(m)) NA_real_ else m$ox, 0),
             decoupled = species == "Mtg", stringsAsFactors = FALSE)
}

.carbon_from_output <- function(lines) {
  out <- list(); ctx <- "?"; i <- 1; b <- 0
  while (i <= length(lines)) {
    l <- lines[i]
    if (grepl("^Initial solution [0-9]+\\.", l)) ctx <- sub("^Initial solution ([0-9]+)\\..*", "initial solution \\1", l)
    if (grepl("^Using (solution|mix) [0-9]+\\.", l)) ctx <- sub("^Using (solution|mix) ([0-9]+)\\..*", "\\1 \\2", l)
    if (grepl("Distribution of species", l, fixed = TRUE)) {
      b <- b + 1; grp <- NA_character_; sp <- character(0); mol <- numeric(0); i <- i + 1
      while (i <= length(lines) && !grepl("Saturation indices", lines[i], fixed = TRUE)) {
        t <- lines[i]
        if (grepl("^\\S", t)) grp <- if (grepl("^(C\\([^)]*\\)|Mtg)\\s+[-0-9.e+]+\\s*$", t)) sub("\\s.*", "", t) else NA_character_
        else if (!is.na(grp) && grepl("^\\s{3}\\S+\\s+[-0-9.e+]+\\s", t)) {
          f <- strsplit(trimws(t), "\\s+")[[1]]; sp <- c(sp, f[1]); mol <- c(mol, as.numeric(f[2]))
        }
        i <- i + 1
      }
      out[[b]] <- .carbon_rows(b, ctx, sp, mol)
    }
    i <- i + 1
  }
  do.call(rbind, out)
}

.carbon_from_selected_output <- function(file) {
  x <- read.delim(file, check.names = FALSE, strip.white = TRUE); names(x) <- trimws(names(x))
  sp <- sub("^m_", "", grep("^m_", names(x), value = TRUE))
  sp <- sp[vapply(sp, function(s) !is.null(carbon_moiety(s)), logical(1))]
  use_mtg_total <- "Mtg" %in% names(x) && !"Mtg" %in% sp       # -totals Mtg, when m_Mtg isn't punched
  do.call(rbind, lapply(seq_len(nrow(x)), function(r) {
    ctx <- sprintf("row %d: %s soln %s", r, trimws(x$state[r]), x$soln[r])
    s <- c(sp, if (use_mtg_total) "Mtg"); m <- c(vapply(sp, function(z) as.numeric(x[[paste0("m_", z)]][r]), 0),
                                                 if (use_mtg_total) as.numeric(x$Mtg[r]))
    .carbon_rows(r, ctx, s, m)
  }))
}
