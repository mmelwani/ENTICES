# phreeqc_tools.R -- locate and run the PHREEQC batch executable from R.
# Also carbon_by_oxidation_state(): carbon in PHREEQC output re-binned by the
# actual oxidation state of C in each species (see the section at the end).
# And database_names() / restrict_to_database(): the PHASES and species a
# database defines, and an input text cut down to them (see the database section).
# And bulk_density(): rock and bulk density through a run, from the mineral and
# IOM moles in its selected output (see the rock density section), and
# system_totals(): each element summed over solution, solids, IOM and gas (at the end).
#
# ENTICES targets PHREEQC 3.8.6 (PHREEQC_EXPECTED_VERSION). Results from other
# versions can differ in minor phases and trace species (see
# docs/reviews/UPDATES_CLEANUP_AND_RUNS_for_review_2026-10-22.md, section 3.3),
# so run_phreeqc() warns when the executable reports a different version.
#
# Executable lookup order: the PHREEQC_EXE environment variable, the default
# Windows install path, then `phreeqc` on the PATH.

PHREEQC_EXPECTED_VERSION <- "3.8.6"

# The code's version for the .pqi headers and RUNINFO (production kit, item 5): "ENTICES commit <sha>", with
# " (dirty: uncommitted changes)" if tracked files differ from it, or "ENTICES commit unknown (not a git
# checkout)". Run from the repo root. Cached for the session.
code_version <- function() {
  if (!is.null(CODE_VERSION$v)) return(CODE_VERSION$v)
  git <- function(...) tryCatch(suppressWarnings(system2("git", c("--no-optional-locks", ...), stdout = TRUE, stderr = FALSE)),
                                error = function(e) character(0))
  sha <- git("rev-parse", "HEAD")
  v <- if (length(sha) == 1 && grepl("^[0-9a-f]{40}$", sha))
    paste0("ENTICES commit ", sha, if (length(git("status", "--porcelain", "--untracked-files=no"))) " (dirty: uncommitted changes)" else "")
  else "ENTICES commit unknown (not a git checkout)"
  CODE_VERSION$v <- v
  v
}
CODE_VERSION <- new.env()

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
  # Every path PHREEQC gets must be absolute before setwd(workdir). normalizePath(mustWork = FALSE) leaves
  # a path that doesn't exist yet unchanged on macOS and Linux (Windows makes it absolute), so a relative
  # output path pointed to the wrong place after setwd() and PHREEQC sat at its file-name prompt. Build
  # the output from its directory, which must exist.
  if (!dir.exists(dirname(output)))
    stop(sprintf("run_phreeqc: the output directory %s does not exist", dirname(output)), call. = FALSE)
  output <- file.path(normalizePath(dirname(output), mustWork = TRUE), basename(output))
  # PHREEQC prompts for a file it cannot open, and at end of input it re-prompts in a tight loop (measured
  # on Windows: 50,661 prompts in 20 s, all captured below), which is how a run can grow until the OS kills
  # it. So make the prompt impossible: the input and database exist (checked above), and the output is
  # created here, so PHREEQC can open all three.
  if (!file.create(output, showWarnings = FALSE))
    stop(sprintf("run_phreeqc: cannot create the output file %s", output), call. = FALSE)
  workdir <- normalizePath(workdir, mustWork = TRUE)
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

# The executable's version, from a one-line run, checked against PHREEQC_EXPECTED_VERSION. The kit
# calls it before the first chunk, so a mismatch shows before hours of running rather than in the
# results. Returns list(version, ok); warns when they differ.
phreeqc_version_check <- function(database = file.path("databases", "coreclath_C4only.dat")) {
  d <- tempfile("phreeqc_version_"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  inp <- file.path(d, "version.pqi"); writeLines(c("SOLUTION 1", "END"), inp)
  r <- run_phreeqc(inp, database, output = file.path(d, "version.out"), workdir = d, warn_version = FALSE)
  ok <- identical(r$version, PHREEQC_EXPECTED_VERSION)
  if (!ok) warning(sprintf("PHREEQC version %s differs from %s, the version ENTICES inputs are generated for: minor phases and trace species may differ",
                           r$version, PHREEQC_EXPECTED_VERSION), call. = FALSE)
  list(version = r$version, ok = ok)
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
  pattern = c("^Mtg$", "C6H8O7", "C4H8N2O3", "C3H7NO2", "CH3COCOO", "C3H5O2|C3H6O2", "C3H8",
              "C2H6NO2|C2H5NO2|C2H4NO2", "C2H5NO$", "CH3CH2OH", "CH3COO", "C2H6", "C2H4",
              "CH3NH2|CH3NH3", "CH3OH", "CH2O$", "HCOO", "SCN", "CN", "CO\\(NH2\\)2", "^CO$",
              "CO3|HCO3|^CO2$"),
  moiety  = c("Mtg", "citric acid", "glycylglycine", "alanine", "pyruvate", "propanoate", "propane",
              "glycine", "acetamide", "ethanol", "acetate", "ethane", "ethene",
              "methylamine", "methanol", "formaldehyde", "formate", "thiocyanate", "cyanide", "urea", "CO",
              "carbonate"),
  nC      = c(1, 6, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1),
  ox      = c(-4, 1, 1, 0, 2/3, -2/3, -8/3, 1, 0, -2, 0, -3, -2, -2, -2, 0, 2, 4, 2, 4, 2, 4),
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

# ---- Names a database defines ----------------------------------------------
#
# coreclath_C4only.dat comments out reduced-carbon species and phases that the
# generators' templates still list. restrict_to_database() removes those names
# from a rendered input, so it runs with the database it recommends.

PHREEQC_DB_KEYWORDS <- c("SOLUTION_MASTER_SPECIES", "SOLUTION_SPECIES", "PHASES", "RATES",
                         "EXCHANGE_MASTER_SPECIES", "EXCHANGE_SPECIES", "SURFACE_MASTER_SPECIES",
                         "SURFACE_SPECIES", "LLNL_AQUEOUS_MODEL_PARAMETERS", "NAMED_EXPRESSIONS",
                         "CALCULATE_VALUES", "END")

# list(phases, species): PHASES names and SOLUTION_SPECIES product species that
# `database` defines (commented-out entries excluded).
database_names <- function(database) {
  l <- sub("#.*", "", readLines(database, warn = FALSE, encoding = "UTF-8"))
  first <- sub("\\s.*", "", l)
  sec <- rep(NA_character_, length(l)); cur <- NA_character_
  for (i in seq_along(l)) {
    if (first[i] %in% PHREEQC_DB_KEYWORDS && !nzchar(trimws(sub("^\\S+", "", l[i])))) cur <- first[i]
    sec[i] <- cur
  }
  col0 <- grepl("^\\S", l) & !first %in% PHREEQC_DB_KEYWORDS
  sp <- l[col0 & sec %in% "SOLUTION_SPECIES" & grepl("=", l, fixed = TRUE)]
  list(phases  = unique(first[col0 & sec %in% "PHASES" & !grepl("=", l, fixed = TRUE)]),
       species = unique(sub("\\s.*", "", trimws(sub("^[^=]*=", "", sp)))))
}

# Input keywords whose blocks or list options name database phases or species.
.PQI_KEYWORD_RX <- "^(SOLUTION|SOLUTION_SPREAD|MIX|REACTION|REACTION_TEMPERATURE|REACTION_PRESSURE|KINETICS|RATES|EQUILIBRIUM_PHASES|GAS_PHASE|SELECTED_OUTPUT|USER_PUNCH|USER_PRINT|INCREMENTAL_REACTIONS|SAVE|USE|END|TITLE|PRINT|KNOBS|EXCHANGE|SURFACE|SOLID_SOLUTIONS|DELETE|RUN_CELLS|COPY|DATABASE)(\\s|$)"
.SO_PHASE_OPTS   <- c("-equilibrium_phases", "-saturation_indices", "-si", "-gases")
.SO_SPECIES_OPTS <- c("-molalities", "-m", "-activities", "-a")

# Removes from PHREEQC input `text` every name `database` does not define (with
# `base`, only names `base` defines and `database` does not, so anything else
# undefined, a typo say, still fails in PHREEQC):
# EQUILIBRIUM_PHASES and GAS_PHASE entries (an EP entry with its option lines),
# and names in SELECTED_OUTPUT -equilibrium_phases / -saturation_indices /
# -gases (phases) and -molalities / -activities (species) lists. Comments are
# kept. Returns the text with attr "dropped" (the names removed) and
# "dropped_nonzero" (EP/GAS_PHASE entries removed with a nonzero amount).
restrict_to_database <- function(text, database, base = NULL) {
  def <- database_names(database)
  if (!is.null(base)) {
    b <- database_names(base)
    gone <- list(phases = setdiff(b$phases, def$phases), species = setdiff(b$species, def$species))
  }
  undefined <- function(x, kind) if (is.null(base)) !x %in% def[[kind]] else x %in% gone[[kind]]
  lines <- strsplit(paste(text, collapse = "\n"), "\n", fixed = TRUE)[[1]]
  keep <- rep(TRUE, length(lines)); dropped <- character(0); nonzero <- character(0)
  block <- ""; so_list <- NULL; in_dropped_ep <- FALSE
  for (i in seq_along(lines)) {
    l <- lines[i]
    if (grepl(.PQI_KEYWORD_RX, l)) { block <- sub("\\s.*", "", l); so_list <- NULL; in_dropped_ep <- FALSE; next }
    code <- sub("#.*", "", l); cmt <- if (grepl("#", l, fixed = TRUE)) sub("^[^#]*", "", l) else ""
    toks <- strsplit(trimws(code), "\\s+")[[1]]; toks <- toks[nzchar(toks)]
    if (!length(toks)) next
    if (block %in% c("EQUILIBRIUM_PHASES", "GAS_PHASE")) {
      if (startsWith(toks[1], "-")) { if (block == "EQUILIBRIUM_PHASES" && in_dropped_ep) keep[i] <- FALSE; next }
      in_dropped_ep <- FALSE
      if (undefined(toks[1], "phases")) {
        keep[i] <- FALSE; dropped <- c(dropped, toks[1]); in_dropped_ep <- TRUE
        amt <- suppressWarnings(as.numeric(toks[if (block == "GAS_PHASE") 2 else 3]))
        if (!is.na(amt) && amt != 0) nonzero <- c(nonzero, sprintf("%s %g", toks[1], amt))
      }
      next
    }
    if (block != "SELECTED_OUTPUT") next
    lead <- regmatches(code, regexpr("^\\s*", code))
    if (startsWith(toks[1], "-")) {
      so_list <- if (toks[1] %in% .SO_PHASE_OPTS) "phases" else if (toks[1] %in% .SO_SPECIES_OPTS) "species" else NULL
      if (is.null(so_list)) next
      opt <- toks[1]; names_ <- toks[-1]
    } else {
      if (is.null(so_list)) next
      opt <- NULL; names_ <- toks
    }
    out_ <- names_[undefined(names_, so_list) & !tolower(names_) %in% c("true", "false")]
    if (!length(out_)) next
    dropped <- c(dropped, out_)
    kept <- setdiff(names_, out_)
    if (is.null(opt) && !length(kept) && !nzchar(cmt)) { keep[i] <- FALSE; next }
    lines[i] <- paste0(lead, paste(c(opt, kept), collapse = " "), if (nzchar(cmt)) paste0(" ", cmt) else "")
  }
  out <- paste(lines[keep], collapse = "\n")
  attr(out, "dropped") <- unique(dropped)
  attr(out, "dropped_nonzero") <- nonzero
  out
}

# ---- Rock density ------------------------------------------------------------
#
# PHREEQC uses a mineral's -Vm only to correct log K for pressure; it reports no solid volumes or
# densities (RHO is the solution's). bulk_density() computes them in R from the mineral and IOM
# moles in a run's SELECTED_OUTPUT, with the molar volumes in data/mineral_molar_volumes.csv
# (written by scripts/export_mineral_volumes.py; 25 C, 1 bar values, no thermal or compressibility
# correction: leaving it out costs ~1% for silicates at 160 C and 26 MPa, less than the porosity
# and IOM-density uncertainty).
#
#   rho_solid = sum(n_i M_i) / sum(n_i V_i) over every mineral with nonzero moles (equilibrium
#               phases, kinetic primaries, schreibersite, Fe/Ni metal, clathrates) plus the
#               remaining IOM, at the density Delarue et al. (2026) Eq. 10 gives for its H/C
#   rho_bulk  = (1 - porosity) rho_solid + porosity rho_pore, where the pores hold the fluid and, since
#               the gas phases of 2026-10-28, the pore gas bubble: rho_pore = (1 - s_g) rho_fluid +
#               s_g rho_gas, with gas saturation s_g = min(Vgas / Vpore, 1) (Vgas_Vpore column) and
#               rho_gas the bubble's mass (sum of its g_ moles x molar mass) over its volume. s_g is
#               capped at 1: a bubble larger than the pores (pore_gas_fate = "retain") is not physical,
#               and the excess is not counted. Such rows are flagged (gas_exceeds_pores) and warned of.

MINERAL_VOLUMES_CSV <- "data/mineral_molar_volumes.csv"
# IOM density (g/cm3) from its atomic H/C: Delarue et al. (2026), EPSL 691, 120227, Eq. 10,
# rho = 1472 - 345 (H/C) + 94 / ((H/C) + 0.122) kg/m3, a fit to kerogens and coals (within 100 kg/m3
# for H/C below 1; the bulk IOM's H/C is 0.753, 1.320 g/cm3). Compression is left out: their
# K0 = 10.2 GPa gives 0.25% at the core's 26 MPa. See docs/MODEL.md (Rock density).
iom_density_delarue2026 <- function(hc) (1472 - 345 * hc + 94 / (hc + 0.122)) / 1000
# The bulk IOM's C and H (mol per kg of IOM): iom_bulk_formula() in iom_module.R
# (tests/bulk_density_check.R checks they agree).
IOM_BULK_C_H <- c(C = 60.611, H = 45.639)
# C and H (mol per mol) of a KINETICS -formula, pseudo-elements resolved as the IOM channels use
# them: Mtg is CH4 ("Mtg 1 H 0.435" = C 1 H 4.435), Sg is S ("Sg 1 H 1" = S 1 H 1, iom_module.R).
iom_formula_CH <- function(formula) {
  tk <- strsplit(trimws(formula), "\\s+")[[1]]
  el <- tk[c(TRUE, FALSE)]; n <- as.numeric(tk[c(FALSE, TRUE)])
  bad <- setdiff(el, c("C", "H", "O", "N", "S", "Mtg", "Sg"))
  if (length(bad)) stop("iom_formula_CH: unknown element in an IOM -formula: ", paste(bad, collapse = ", "))
  c(C = sum(n[el %in% c("C", "Mtg")]), H = sum(n[el == "H"]) + 4 * sum(n[el == "Mtg"]))
}

# The volumes table: phase, formula, molar mass (g/mol), Vm (cm3/mol).
read_mineral_volumes <- function(file = MINERAL_VOLUMES_CSV) {
  v <- utils::read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  v[!is.na(v$Vm_cm3_mol), c("phase", "formula", "molar_mass_g_mol", "Vm_cm3_mol")]
}

# Element gram formula weights from a database's SOLUTION_MASTER_SPECIES (element rows only),
# pseudo-elements (Mtg, Hdg, ...) included.
element_gfw <- function(database) {
  l <- sub("#.*", "", readLines(database, warn = FALSE, encoding = "UTF-8"))
  i <- which(trimws(l) == "SOLUTION_MASTER_SPECIES"); j <- which(trimws(l) == "SOLUTION_SPECIES")[1]
  t <- strsplit(trimws(l[(i + 1):(j - 1)]), "\\s+")
  t <- t[vapply(t, function(x) length(x) >= 5 && !grepl("\\(", x[1]), logical(1))]
  stats::setNames(as.numeric(vapply(t, `[`, "", 5)), vapply(t, `[`, "", 1))
}
# Mass (g/mol) of a KINETICS -formula ("C 1 H 1 O 1.3115", "Mtg 1 H -2.2").
kinetics_formula_mass <- function(formula, gfw) {
  t <- strsplit(trimws(formula), "\\s+")[[1]]
  el <- t[c(TRUE, FALSE)]; n <- as.numeric(t[c(FALSE, TRUE)])
  if (any(!el %in% names(gfw))) stop("no gram formula weight for: ", paste(setdiff(el, names(gfw)), collapse = ", "))
  sum(gfw[el] * n)
}

# Solid density from phase moles (named numeric, mol) and a mass of IOM (g). Returns
# list(rho_solid, mass_g, volume_cm3, breakdown: phase, mol, mass_g, volume_cm3, mass_frac, vol_frac).
# A phase with moles but no molar volume is an error (so nothing is silently dropped).
solid_density <- function(moles, iom_g = 0, iom_density = iom_density_delarue2026(IOM_BULK_C_H[["H"]] / IOM_BULK_C_H[["C"]]),
                          volumes = read_mineral_volumes()) {
  moles <- moles[!is.na(moles) & moles > 0]
  miss <- setdiff(names(moles), volumes$phase)
  if (length(miss)) stop("no molar volume in ", MINERAL_VOLUMES_CSV, " for: ", paste(miss, collapse = ", "))
  k <- match(names(moles), volumes$phase)
  b <- data.frame(phase = names(moles), mol = unname(moles),
                  mass_g = unname(moles) * volumes$molar_mass_g_mol[k],
                  volume_cm3 = unname(moles) * volumes$Vm_cm3_mol[k], stringsAsFactors = FALSE)
  if (iom_g > 0) b <- rbind(b, data.frame(phase = "IOM", mol = NA_real_, mass_g = iom_g, volume_cm3 = iom_g / iom_density))
  m <- sum(b$mass_g); v <- sum(b$volume_cm3)
  b$mass_frac <- b$mass_g / m; b$vol_frac <- b$volume_cm3 / v
  list(rho_solid = m / v, mass_g = m, volume_cm3 = v, breakdown = b[order(-b$vol_frac), ])
}

# The solids a generated .pqi starts with: KINETICS 1 -m0 (minerals; IOM channels with their
# -formula) and EQUILIBRIUM_PHASES 1 initial amounts, and the initial IOM mass from the
# "# Organic content: X wt.%" header (per kg of rock; 0 if absent).
pqi_initial_solids <- function(pqi) {
  l <- readLines(pqi, warn = FALSE)
  blk <- function(kw) { i <- grep(paste0("^", kw, "$"), l)[1]; if (is.na(i)) return(character(0))
    j <- i + 1; while (j <= length(l) && !grepl("^[A-Z_]+( |$)", l[j])) j <- j + 1; l[(i + 1):(j - 1)] }
  kin <- blk("KINETICS 1"); ep <- blk("EQUILIBRIUM_PHASES 1")
  nm <- character(0); m0 <- numeric(0); fo <- character(0)
  for (s in kin) {
    if (grepl("^\\s+[A-Za-z]", s)) { nm <- c(nm, trimws(s)); m0 <- c(m0, NA); fo <- c(fo, NA) }  # an entry (any indent)
    else if (grepl("^\\s+-m0\\s", s)) m0[length(m0)] <- as.numeric(sub("^\\s+-m0\\s+", "", s))
    else if (grepl("^\\s+-formula\\s", s)) fo[length(fo)] <- trimws(sub("^\\s+-formula\\s+", "", s))
  }
  epl <- grep("^\\s+[A-Za-z]", sub("#.*", "", ep), value = TRUE)
  ept <- strsplit(trimws(epl), "\\s+")
  epm <- stats::setNames(vapply(ept, function(t) as.numeric(t[3]), 0), vapply(ept, `[`, "", 1))
  org <- suppressWarnings(as.numeric(sub(".*: ([0-9.]+) wt.*", "\\1", grep("^# Organic content: ", l, value = TRUE)[1])))
  iom <- grepl("^IOM_", nm)
  # the SELECTED_OUTPUT -totals elements (the option line and its continuation lines)
  ti <- grep("^\\s+-totals\\s", l)[1]; tl <- character(0)
  if (!is.na(ti)) { j <- ti; repeat { tl <- c(tl, sub("#.*", "", l[j])); j <- j + 1
    if (j > length(l) || grepl("^\\s*-", l[j]) || !grepl("^\\s", l[j])) break } }
  totals <- setdiff(unlist(strsplit(trimws(sub("^\\s+-totals", "", tl)), "\\s+")), "")
  list(totals = totals, minerals = c(stats::setNames(m0[!iom], nm[!iom]), epm[epm > 0]),
       iom = data.frame(name = nm[iom], m0 = m0[iom], formula = fo[iom], stringsAsFactors = FALSE),
       iom_mass0_g = if (is.na(org)) 0 else org / 100 * 1000,
       database = sub("^# Database: ", "", grep("^# Database: ", l, value = TRUE)[1]))
}

# Rock and bulk density through a run. csv: a SELECTED_OUTPUT file; its .pqi (same name) must be
# beside it, for the IOM formulas and initial amounts. porosity: the model's porosity
# (Simulator$params$porosity, 0.32 by default). Rows: t = 0 (from the .pqi) and every porewater
# row that holds the rock (solution 1 with nonzero solids: each REACT step). fluid_density: NULL
# uses the run's Rho_soln column (PHREEQC's RHO; under the LLNL aqueous model this is ~1 g/cm3 for
# water at any T and P, see README), else a number in g/cm3. iom_density: "delarue2026" (default)
# evaluates Delarue et al. (2026) Eq. 10 at the remaining IOM's H/C on each row (the bulk IOM less
# what the channels have released); a number fixes it (g/cm3). origin_pqi: for a restart, the first
# run's .pqi, which holds the IOM's t = 0 inventory (a restart's own .pqi starts its IOM channels at
# what was left); the t = 0 row is then the restart's start. Returns list(summary, breakdown).
bulk_density <- function(csv, porosity, pqi = sub("\\.csv$", ".pqi", csv), iom_density = "delarue2026",
                         fluid_density = NULL, volumes = read_mineral_volumes(), origin_pqi = pqi) {
  if (missing(porosity) || !is.numeric(porosity) || porosity < 0 || porosity >= 1)
    stop("porosity: give the model's porosity (Simulator$params$porosity)")
  init <- pqi_initial_solids(pqi)
  n_iom_start <- NULL                                       # IOM channel amounts at this run's start
  if (!identical(normalizePath(origin_pqi, mustWork = FALSE), normalizePath(pqi, mustWork = FALSE))) {
    org <- pqi_initial_solids(origin_pqi)
    n_iom_start <- init$iom$m0[match(org$iom$name, init$iom$name)]
    if (anyNA(n_iom_start)) stop("bulk_density: the restart's .pqi lacks IOM channels of origin_pqi")
    init$iom <- org$iom; init$iom_mass0_g <- org$iom_mass0_g
  }
  gfw <- element_gfw(init$database)
  iom_M <- vapply(init$iom$formula, kinetics_formula_mass, 0, gfw = gfw)
  if (identical(iom_density, "delarue2026")) {
    iom_CH <- if (nrow(init$iom)) t(vapply(init$iom$formula, iom_formula_CH, c(C = 0, H = 0))) else matrix(0, 0, 2, dimnames = list(NULL, c("C", "H")))
    CH0 <- IOM_BULK_C_H * init$iom_mass0_g / 1000           # mol in the rock at t = 0
    iom_hc <- function(n_iom) { ch <- CH0 - colSums((init$iom$m0 - n_iom) * iom_CH); ch[["H"]] / ch[["C"]] }
    iom_rho <- function(hc) iom_density_delarue2026(hc)
  } else if (is.numeric(iom_density) && length(iom_density) == 1 && iom_density > 0) {
    iom_hc <- function(n_iom) NA_real_; iom_rho <- function(hc) iom_density
  } else stop('iom_density: "delarue2026" or a density in g/cm3')
  x <- utils::read.delim(csv, check.names = FALSE, strip.white = TRUE); names(x) <- trimws(names(x))
  x$state <- trimws(x$state)
  # Equilibrium-phase columns, by index. A phase named like an element (Ni) shares its name with
  # that element's -totals column, which comes first: use such a phase only if the name appears
  # twice, and take the last. Fe metal is in its own Fe_metal column.
  ep_names <- setdiff(intersect(names(x), volumes$phase), "Fe")
  ep_names <- ep_names[!ep_names %in% init$totals | vapply(ep_names, function(p) sum(names(x) == p) > 1, logical(1))]
  ep_idx <- vapply(ep_names, function(p) max(which(names(x) == p)), 0L)
  kin_cols <- grep("^k_", names(x), value = TRUE)
  kin_min <- intersect(sub("^k_", "", kin_cols), volumes$phase)
  rows <- list(); bd <- list()
  add <- function(label, time, moles, iom_g, hc, rho_f, s_g = 0, rho_g = 0, vg_vp = 0) {
    rho_iom <- iom_rho(hc)
    s <- solid_density(moles, iom_g, rho_iom, volumes)
    rho_p <- (1 - s_g) * rho_f + s_g * rho_g
    rows[[length(rows) + 1]] <<- data.frame(row = label, Time_Years = time, rho_solid = s$rho_solid,
                                             rho_fluid = rho_f, gas_saturation = s_g, rho_gas = rho_g,
                                             Vgas_Vpore = vg_vp, gas_exceeds_pores = vg_vp > 1,
                                             rho_bulk = (1 - porosity) * s$rho_solid + porosity * rho_p,
                                             solid_mass_g = s$mass_g, solid_volume_cm3 = s$volume_cm3, iom_g = iom_g,
                                             iom_HC = hc, iom_density = rho_iom)
    bd[[length(bd) + 1]] <<- cbind(row = label, s$breakdown)
  }
  rho_col <- if ("Rho_soln" %in% names(x)) "Rho_soln" else NULL
  if (is.null(n_iom_start)) {
    add("t0", 0, init$minerals, init$iom_mass0_g, iom_hc(init$iom$m0), if (!is.null(fluid_density)) fluid_density else 1.0)  # pure water
  } else {                                                  # a restart: its start, IOM measured from the origin
    add("t0", NA_real_, init$minerals, max(0, init$iom_mass0_g - sum((init$iom$m0 - n_iom_start) * iom_M)), iom_hc(n_iom_start),
        if (!is.null(fluid_density)) fluid_density else 1.0)
  }
  for (r in which(x$state == "react" & trimws(x$soln) == "1")) {
    mol <- c(stats::setNames(as.numeric(x[r, paste0("k_", kin_min)]), kin_min),
             stats::setNames(as.numeric(unlist(x[r, ep_idx])), ep_names),
             if ("Fe_metal" %in% names(x)) c(Fe = as.numeric(x$Fe_metal[r])))
    if (sum(mol[mol > 0], na.rm = TRUE) == 0) next          # a MIX row: no solids
    n_iom <- if (nrow(init$iom)) as.numeric(x[r, paste0("k_", init$iom$name)]) else numeric(0)
    iom_g <- max(0, init$iom_mass0_g - sum((init$iom$m0 - n_iom) * iom_M))
    rho_f <- if (!is.null(fluid_density)) fluid_density else if (!is.null(rho_col)) as.numeric(x[[rho_col]][r]) else 1.0
    # the pore gas bubble (gas phase in use on this REACT row)
    vg <- if ("Vgas_L" %in% names(x)) as.numeric(x$Vgas_L[r]) else 0
    vg_vp <- if (vg > 0 && "Vgas_Vpore" %in% names(x)) as.numeric(x$Vgas_Vpore[r]) else 0
    s_g <- min(vg_vp, 1)                                   # gas beyond the pore volume is not counted
    mg <- sum(vapply(grep("^g_", names(x), value = TRUE), function(k) {
      fe <- formula_elements(sub("\\(g\\)$", "", sub("^g_", "", k)))
      as.numeric(x[r, k]) * sum(gfw[names(fe)] * fe) }, 0))  # g
    rho_g <- if (vg > 0) mg / (vg * 1000) else 0           # g / cm3
    add(as.character(r), if ("Time_Years" %in% names(x)) as.numeric(x$Time_Years[r]) else NA_real_, mol, iom_g, iom_hc(n_iom), rho_f, s_g, rho_g, vg_vp)
  }
  s <- do.call(rbind, rows)
  if (any(s$gas_exceeds_pores))
    warning(sprintf("bulk_density: the pore gas exceeds the pore volume on %d of %d rows (Vgas/Vpore up to %.2f); gas saturation is capped at 1 and the excess is not counted (see gas_exceeds_pores)",
                    sum(s$gas_exceeds_pores), nrow(s), max(s$Vgas_Vpore)))
  # End-state check against the present-day core density (Bagheri et al. 2025 Table 5): a check on the
  # final state, not an input.
  s$core_density_ok <- s$rho_bulk >= BAGHERI_CORE_DENSITY_G_CM3[1] & s$rho_bulk <= BAGHERI_CORE_DENSITY_G_CM3[2]
  if (!s$core_density_ok[nrow(s)])
    warning(sprintf("bulk_density: the final rho_bulk %.3f g/cm3 is outside Bagheri et al. (2025)'s core density %.2f-%.2f g/cm3",
                    s$rho_bulk[nrow(s)], BAGHERI_CORE_DENSITY_G_CM3[1], BAGHERI_CORE_DENSITY_G_CM3[2]))
  list(summary = s, breakdown = do.call(rbind, bd))
}

# Present-day core bulk density, Bagheri et al. (2025) Table 5 (both dissipation scenarios), g/cm3.
BAGHERI_CORE_DENSITY_G_CM3 <- c(2.160, 2.400)

# Porosity that gives a core bulk density rho_core from a solid of density rho_solid and a pore fluid of
# rho_fluid (all g/cm3): rho_core = (1 - phi) rho_solid + phi rho_fluid, so phi = (rho_solid - rho_core) /
# (rho_solid - rho_fluid). Vectorised. At t = 0 with 15 wt% IOM (rho_solid 2.920) Bagheri's 2.16-2.40 gives
# phi 0.271-0.396; less IOM means a denser solid and a higher phi for the same core density.
phi_from_core_density <- function(rho_solid, rho_core, rho_fluid = 1.0) (rho_solid - rho_core) / (rho_solid - rho_fluid)

# ---- Element totals of the whole system ---------------------------------------
#
# system_totals() adds up each element over a row's solution (molality x water mass, plus the water
# itself for H and O), equilibrium phases, kinetic reactants (minerals, Schreibersite, IOM channels by
# their -formula) and, where the row has a gas phase, its components. The restart and gas-phase
# mass-balance checks compare these.

# Element counts of a chemical formula as PHREEQC writes it: "Mg3Si2O5(OH)4", "CaSO4:2H2O",
# "Na.96Al.96Si2.04O6:H2O", "Mtg:6H2O"; or a KINETICS -formula "C 1 H 1 O 1.3115".
formula_elements <- function(f) {
  f <- trimws(f)
  if (grepl("\\s", f)) {                                   # "El n El n ..."
    t <- strsplit(f, "\\s+")[[1]]
    return(tapply(as.numeric(t[c(FALSE, TRUE)]), t[c(TRUE, FALSE)], sum))
  }
  add <- function(a, b, k = 1) { for (e in names(b)) a[e] <- (if (is.na(a[e])) 0 else a[e]) + k * b[[e]]; a }
  group <- function(s) {
    out <- numeric(0); i <- 1
    while (i <= nchar(s)) {
      ch <- substr(s, i, i)
      if (ch == "(") {
        depth <- 1; j <- i
        while (depth > 0) { j <- j + 1; c2 <- substr(s, j, j); depth <- depth + (c2 == "(") - (c2 == ")") }
        m <- regmatches(substring(s, j + 1), regexpr("^[0-9.]*", substring(s, j + 1)))
        out <- add(out, group(substr(s, i + 1, j - 1)), if (nzchar(m)) as.numeric(m) else 1); i <- j + 1 + nchar(m)
      } else {
        m <- regmatches(substring(s, i), regexpr("^[A-Z][a-z]*[0-9.]*", substring(s, i)))
        if (!length(m)) stop("cannot parse formula ", f)
        el <- sub("[0-9.]+$", "", m); n <- sub("^[A-Za-z]+", "", m)
        out <- add(out, stats::setNames(if (nzchar(n)) as.numeric(n) else 1, el)); i <- i + nchar(m)
      }
    }
    out
  }
  tot <- numeric(0)
  for (part in strsplit(f, ":", fixed = TRUE)[[1]]) {
    k <- regmatches(part, regexpr("^[0-9.]+", part)); body <- sub("^[0-9.]+", "", part)
    tot <- add(tot, group(body), if (length(k) && nzchar(k)) as.numeric(k) else 1)
  }
  tot
}

# Element totals (mol) of row i of a selected-output data frame x (read with check.names = FALSE).
# elements: the -totals elements (pqi_initial_solids(pqi)$totals); iom: pqi_initial_solids(pqi)$iom.
# H and O: the H_total_mol / O_total_mol columns (water included) where the run punches them; else
# the -totals H and O (which leave out water) plus water as mass_H2O / 18.01528 g/mol (approximate).
system_totals <- function(x, i, elements, iom = NULL, volumes = utils::read.csv(MINERAL_VOLUMES_CSV, stringsAsFactors = FALSE)) {
  nm <- names(x); w <- as.numeric(x$mass_H2O[i])
  tot <- stats::setNames(vapply(elements, function(e) as.numeric(x[i, match(e, nm)]) * w, 0), elements)
  if (all(c("H_total_mol", "O_total_mol") %in% nm)) tot[c("H", "O")] <- c(as.numeric(x$H_total_mol[i]), as.numeric(x$O_total_mol[i]))
  else tot[c("H", "O")] <- tot[c("H", "O")] + c(2, 1) * w / 0.01801528
  addf <- function(n, f) { if (is.na(n) || n == 0) return(); fe <- formula_elements(f)
    for (e in names(fe)) tot[e] <<- (if (is.na(tot[e])) 0 else tot[e]) + n * fe[[e]] }
  form <- stats::setNames(volumes$formula, volumes$phase)
  # equilibrium phases: the last column of each name (Ni shares its name with the Ni total)
  for (p in setdiff(intersect(unique(nm), volumes$phase), c("Fe", elements)))
    addf(as.numeric(x[i, max(which(nm == p))]), form[[p]])
  for (p in intersect(elements, volumes$phase)) if (sum(nm == p) > 1) addf(as.numeric(x[i, max(which(nm == p))]), form[[p]])
  # gases listed as equilibrium phases (an infinite reservoir at fixed partial pressure)
  for (g in nm[grepl("\\(g\\)$", nm) & !grepl("^(g_|d_)", nm)]) addf(as.numeric(x[i, max(which(nm == g))]), sub("\\(g\\)$", "", g))
  if ("Fe_metal" %in% nm) addf(as.numeric(x$Fe_metal[i]), "Fe")
  for (k in grep("^k_", nm, value = TRUE)) {
    p <- sub("^k_", "", k)
    f <- if (!is.null(iom) && p %in% iom$name) iom$formula[match(p, iom$name)] else if (p %in% names(form)) form[[p]] else NA
    if (is.na(f)) stop("system_totals: no formula for kinetic reactant ", p)
    addf(as.numeric(x[i, k]), f)
  }
  for (g in grep("^g_", nm, value = TRUE)) addf(as.numeric(x[i, g]), sub("\\(g\\)$", "", sub("^g_", "", g)))
  tot
}

# The formal oxidation state of Fe in a mineral formula, by charge balance with fixed valences for the other
# elements (O -2, H +1, Cl -1, F -1, Si +4, Al +3, Mg/Ca/Mn/Ni +2, Na/K +1, Cr +3; S -2 without O, -1 in
# pyrite/marcasite, +6 with O; P +5 and C +4 with O). Schreibersite is entered zero-valent (Fe3P, as the
# generators treat it). Fe-free formulas return NA.
fe_valence <- function(formula, phase = "") {
  if (phase == "Schreibersite") return(0)
  f <- unlist(formula_elements(formula)); has <- function(e) e %in% names(f) && f[[e]] > 0
  if (!has("Fe")) return(NA_real_)
  v <- c(O = -2, H = 1, Cl = -1, F = -1, Si = 4, Al = 3, Mg = 2, Ca = 2, Mn = 2, Ni = 2, Na = 1, K = 1, Cr = 3, Ba = 2, Sr = 2)
  hasO <- has("O")
  if (has("S")) v["S"] <- if (hasO) 6 else if (phase %in% c("Pyrite", "Marcasite")) -1 else -2
  if (has("P")) v["P"] <- if (hasO) 5 else -3
  if (has("C")) v["C"] <- if (hasO) 4 else NA
  other <- setdiff(names(f), "Fe")
  if (any(is.na(v[other]))) stop("fe_valence: no valence for ", paste(other[is.na(v[other])], collapse = ", "), " in ", formula)
  -sum(f[other] * v[other]) / f[["Fe"]]
}

# The run's H2 budget per orbit (item 6, PI 2026-10-29), on the run's normalised basis (mol per kg of the
# reacting box's rock). Two counts of the H2 generated since t = 0:
#  - lower bounds, as before: H2_from_metal (Fe and Ni metal to Fe(II) and Ni(II), 1 H2 each) and
#    H2_from_schreibersite (zero-valent Fe3P to 3 Fe(II) and phosphate, 5.5 H2 each);
#  - by electron balance (item 4 of 10-31, PI 2026-11-01), exact for Fe: H2_from_Fe = (1/2) the change of
#    sum(n z_Fe) over every Fe-bearing phase, with z_Fe from each formula (fe_valence(): magnetite 8/3,
#    Fe(III) phases 3, troilite and other Fe(II) phases 2, metal and schreibersite 0) in the porewater's
#    phases (REACT row: equilibrium phases and kinetic reactants) and the ocean's (MIX 4 row), and the dissolved
#    Fe of the porewater and ocean counted as Fe(II) (the CSV has no Fe(3) column; dissolved_Fe_mol bounds the
#    error, at most 1/2 H2 per mol); H2_from_Ni = Ni metal to Ni(II); H2_from_P = (5/2) per schreibersite P to
#    phosphate; H2_generated their sum;
# and where the H2 is at the end of each orbit: dissolved in the ocean (MIX 2 row) and in its gas (MIX 4 row),
# dissolved in the porewater (MIX 1 row) and in the pore gas (the REACT row under "retain"; none under
# "to_ocean", moved to the ocean; <base>_poregas.csv under "retain_to_pore_volume"); H2_present their sum, and
# H2_other_redox = H2_generated - H2_present, the H2 the other redox couples (carbon, sulfur, nitrogen) took up
# (positive) or gave (negative).
# csv: a tidal_kinetic run's CSV (its .pqi beside it). Returns one row per orbit.
h2_budget <- function(csv, pqi = sub("\\.csv$", ".pqi", csv)) {
  init <- pqi_initial_solids(pqi)
  fate <- sub("^# pore_gas_fate: ", "", grep("^# pore_gas_fate: ", readLines(pqi, n = 80), value = TRUE)[1])
  x <- utils::read.delim(csv, check.names = FALSE, strip.white = TRUE); names(x) <- trimws(names(x))
  r <- x[trimws(x$state) == "react", , drop = FALSE]; n <- nrow(r) %/% 4
  pg <- if (identical(fate, "retain_to_pore_volume")) { y <- utils::read.delim(sub("\\.csv$", "_poregas.csv", csv), check.names = FALSE,
                                                                              strip.white = TRUE); names(y) <- trimws(names(y)); y } else NULL
  last <- function(d, i, nm) { k <- which(names(d) == nm); if (!length(k)) 0 else as.numeric(d[i, max(k)]) }
  g <- function(m, k) if (is.na(m) || m == 0) 0 else m
  fe0 <- g(init$minerals["Fe"]); ni0 <- g(init$minerals["Ni"]); sch0 <- g(init$minerals["Schreibersite"])
  # Fe-bearing phases and their Fe valence; the t = 0 sum of n z_Fe over the initial solids
  v <- read_mineral_volumes(); v <- v[grepl("Fe", v$formula), ]
  zFe <- stats::setNames(mapply(function(f, ph) tryCatch(fe_valence(f, ph), error = function(e) NA_real_), v$formula, v$phase), v$phase)
  unknown <- names(zFe)[is.na(zFe)]; zFe <- zFe[!is.na(zFe)]   # Fe phases fe_valence() cannot assign: must be absent
  nFe <- stats::setNames(vapply(v$formula, function(f) { e <- unlist(formula_elements(f)); if ("Fe" %in% names(e)) e[["Fe"]] else 0 }, 0), v$phase)
  nz0 <- sum(vapply(names(init$minerals), function(ph) if (ph %in% names(zFe)) init$minerals[[ph]] * nFe[[ph]] * zFe[[ph]] else 0, 0))
  nz <- function(row, kinetic) {   # sum n z_Fe over a row's Fe phases (equilibrium-phase columns; k_ columns too)
    s <- 0
    for (ph in intersect(unknown, names(row))) { m <- as.numeric(row[1, max(which(names(row) == ph))])
      if (!is.na(m) && m != 0) stop("h2_budget: no Fe valence for ", ph, ", present at ", m, " mol") }
    for (ph in intersect(names(zFe), names(row))) if (ph != "Fe") { m <- as.numeric(row[1, max(which(names(row) == ph))]); if (!is.na(m)) s <- s + m * nFe[[ph]] * zFe[[ph]] }
    if (kinetic) for (k in grep("^k_", names(row), value = TRUE)) { ph <- sub("^k_", "", k)
      if (ph %in% names(zFe)) s <- s + as.numeric(row[[k]]) * nFe[[ph]] * zFe[[ph]] }
    s
  }
  do.call(rbind, lapply(seq_len(n), function(o) {
    i <- 4 * (o - 1)
    fe <- if ("Fe_metal" %in% names(r)) as.numeric(r$Fe_metal[i + 1]) else 0
    ni <- last(r, i + 1, "Ni"); ni <- if (ni0 > 0) ni else 0
    sch <- if ("k_Schreibersite" %in% names(r)) as.numeric(r$k_Schreibersite[i + 1]) else 0
    pore_gas <- switch(fate, retain = last(r, i + 1, "g_H2(g)"), retain_to_pore_volume = last(pg, o, "g_H2(g)"), 0)
    fe_d <- as.numeric(r$Fe[i + 3]) * as.numeric(r$mass_H2O[i + 3]) + as.numeric(r$Fe[i + 4]) * as.numeric(r$mass_H2O[i + 4])
    nz_end <- nz(r[i + 1, , drop = FALSE], TRUE) + nz(r[i + 2, , drop = FALSE], FALSE) + 2 * fe_d
    b <- data.frame(orbit = o, Time_Years = as.numeric(r$Time_Years[i + 1]),
               H2_from_metal = (fe0 - fe) + (ni0 - ni), H2_from_schreibersite = 5.5 * (sch0 - sch),
               H2_from_Fe = (nz_end - nz0) / 2, H2_from_Ni = ni0 - ni, H2_from_P = 2.5 * (sch0 - sch),
               ocean_dissolved = as.numeric(r$m_H2[i + 4]) * as.numeric(r$mass_H2O[i + 4]), ocean_gas = last(r, i + 2, "g_H2(g)"),
               porewater_dissolved = as.numeric(r$m_H2[i + 3]) * as.numeric(r$mass_H2O[i + 3]), pore_gas = pore_gas,
               dissolved_Fe_mol = fe_d)
    b$H2_generated <- b$H2_from_Fe + b$H2_from_Ni + b$H2_from_P
    b$H2_present <- b$ocean_dissolved + b$ocean_gas + b$porewater_dissolved + b$pore_gas
    b$H2_other_redox <- b$H2_generated - b$H2_present
    b
  }))
}

# Scale a normalised budget (mol per kg of the box's rock) to mol in the reacting box, the whole permeable
# layer and the whole core: x rock mass, V (rho_core - phi rho_water) with the Simulator's geometry. sim: the
# run's Simulator.
scale_to_core <- function(budget, sim, cols = setdiff(names(budget), c("orbit", "Time_Years"))) {
  cfg <- sim$config; R <- cfg$core_radius; phi <- sim$params$porosity
  rock <- function(L) 4 / 3 * pi * (R^3 - max(0, R - L)^3) * (cfg$core_density - phi * CONSTANTS$FLUID$density)
  m <- c(box = rock(sim$reactive_thickness()), layer = rock(sim$params$layer_thickness), core = rock(R))
  out <- budget[c("orbit", "Time_Years")]
  for (s in names(m)) for (cc in cols) out[[paste0(cc, "_", s, "_mol")]] <- budget[[cc]] * m[[s]]
  attr(out, "rock_kg") <- m
  out
}

# Column-by-column comparison of two runs' rows (the comparison tests; PI 2026-11-01). A value passes if
# within `rel` relative or `abs_floor` absolute (molal). Log-scale columns (pe, pH, la_*, si_*) are compared only
# where the quantity they are the log of is poised: the relevant species above `gate` (1e-10 molal) in both
# runs; there they pass within `log_abs` (1e-4) log units, and below the gate the molalities' absolute floor
# covers them (the log of a trace or unpoised quantity carries solver noise: la_SO4-2 differed by 1.2e-4 log
# units at <= 4e-9 molal). The relevant species (log_column_gate()):
#   la_X   m_X;
#   pH     the larger of m_H+ and m_OH- (water poises pH; always compared in practice);
#   si_X   the smallest dissolved total (molal) of the elements of X's formula other than H and O;
#   pe     a poising redox pair with both members above the gate: H2 (against water), Fe+2 + FeOH+ / Fe+3,
#          SO4-2 / HS- + H2S, HCO3- + CO3-2 + CO2 / Mtg, N2 / NH4+ + NH3 (any one suffices).
# Returns, per column, the largest relative difference among the values that fail (0 if none); a log-scale
# column that fails is reported by its difference in log units.
is_log_column <- function(name) name %in% c("pe", "pH") | startsWith(name, "la_") | startsWith(name, "si_")
log_column_gate <- function(x, col, gate = 1e-10, formulas = NULL) {
  m <- function(s) { k <- paste0("m_", s); v <- rep(0, nrow(x)); for (kk in k[k %in% names(x)]) v <- v + as.numeric(x[[kk]]); v }
  if (startsWith(col, "la_")) { s <- sub("^la_", "", col); if (!paste0("m_", s) %in% names(x)) stop("log_column_gate: no m_", s); return(m(s) > gate) }
  if (col == "pH") return(pmax(m("H+"), m("OH-")) > gate)
  if (col == "pe") {
    pairs <- list(cbind(m("H2"), Inf), cbind(m(c("Fe+2", "FeOH+")), m("Fe+3")), cbind(m("SO4-2"), m(c("HS-", "H2S"))),
                  cbind(m(c("HCO3-", "CO3-2", "CO2")), m("Mtg")), cbind(m("N2"), m(c("NH4+", "NH3"))))
    return(Reduce(`|`, lapply(pairs, function(pp) pmin(pp[, 1], pp[, 2]) > gate)))
  }
  if (startsWith(col, "si_")) {
    ph <- sub("^si_", "", col); if (is.null(formulas)) formulas <- { v <- read_mineral_volumes(); stats::setNames(v$formula, v$phase) }
    if (!ph %in% names(formulas)) stop("log_column_gate: no formula for ", ph)
    el <- setdiff(names(unlist(formula_elements(formulas[[ph]]))), c("H", "O"))
    if (!all(el %in% names(x))) stop("log_column_gate: no total column for ", paste(setdiff(el, names(x)), collapse = ", "), " (", ph, ")")
    return(Reduce(pmin, lapply(el, function(e) as.numeric(x[[e]]))) > gate)
  }
  stop("log_column_gate: not a log-scale column: ", col)
}
column_differences <- function(A, B, cols, abs_floor = 1e-12, log_abs = 1e-4, gate = 1e-10) {
  v <- read_mineral_volumes(); fm <- stats::setNames(v$formula, v$phase)
  vapply(cols, function(c) { a <- as.numeric(A[[c]]); b <- as.numeric(B[[c]])
    if (is_log_column(c)) { k <- abs(a - b) > log_abs & log_column_gate(A, c, gate, fm) & log_column_gate(B, c, gate, fm)
      return(if (any(k, na.rm = TRUE)) max(abs(a - b)[k], na.rm = TRUE) else 0) }
    k <- abs(a - b) > abs_floor
    if (!any(k, na.rm = TRUE)) 0 else max(abs(a - b)[k] / pmax(abs(a), abs(b))[k], na.rm = TRUE) }, 0)
}

# Approach to steady state (item 7, PI 2026-10-29), per orbit, on the run's normalised basis: the net
# elemental flux into the ocean (the change of the ocean side, its water, minerals and gas, over the orbit:
# the ocean exchanges only with the porewater), and the net change of the secondary minerals in the
# porewater (REACT row) and in the ocean (MIX 4 row), in moles of formula units. A steady state shows as
# fluxes that stop changing, not necessarily zero (docs/MODEL.md).
steady_state_diagnostics <- function(csv, pqi = sub("\\.csv$", ".pqi", csv), elements = c("C", "Mtg", "Na", "Mg", "Si", "Fe", "Ca", "S", "Cl", "N")) {
  init <- pqi_initial_solids(pqi)
  fate <- sub("^# pore_gas_fate: ", "", grep("^# pore_gas_fate: ", readLines(pqi, n = 80), value = TRUE)[1])
  x <- utils::read.delim(csv, check.names = FALSE, strip.white = TRUE); names(x) <- trimws(names(x))
  v <- read_mineral_volumes()
  r <- x[trimws(x$state) == "react", , drop = FALSE]; n <- nrow(r) %/% 4
  ini <- x[trimws(x$state) == "i_soln", , drop = FALSE]
  el <- intersect(c(elements, "H", "O"), names(x))
  ocean <- function(o) whole_system_parts(r[seq_len(4 * o), , drop = FALSE], init, fate)$ocean[el]
  oc0 <- if (nrow(ini) >= 2) stats::setNames(vapply(el, function(e) as.numeric(ini[2, e]) * as.numeric(ini$mass_H2O[2]), 0), el) else NULL
  if (!is.null(oc0)) oc0[c("H", "O")] <- c(as.numeric(ini$H_total_mol[2]), as.numeric(ini$O_total_mol[2]))
  primaries <- c(names(init$minerals), "Fe", "Ni")
  sec <- function(i) { p <- setdiff(intersect(names(r), v$phase), c(primaries, init$totals))
    sum(vapply(p, function(q) as.numeric(r[i, max(which(names(r) == q))]), 0), na.rm = TRUE) }
  prev <- oc0; p1 <- 0; p2 <- 0
  do.call(rbind, lapply(seq_len(n), function(o) {
    i <- 4 * (o - 1); oc <- ocean(o)
    flux <- if (is.null(prev)) rep(NA_real_, length(el)) else oc - prev
    s1 <- sec(i + 1); s2 <- sec(i + 2)
    row <- data.frame(orbit = o, Time_Years = as.numeric(r$Time_Years[i + 1]), sec_minerals_pore_change = s1 - p1, sec_minerals_ocean_change = s2 - p2)
    for (k in seq_along(el)) row[[paste0("flux_to_ocean_", el[k])]] <- flux[k]
    prev <<- oc; p1 <<- s1; p2 <<- s2
    row
  }))
}
