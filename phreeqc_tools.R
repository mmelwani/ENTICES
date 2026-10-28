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
#               remaining IOM at IOM_DENSITY_G_CM3
#   rho_bulk  = (1 - porosity) rho_solid + porosity rho_pore, where the pores hold the fluid and, since
#               the gas phases of 2026-10-28, the pore gas bubble: rho_pore = (1 - s_g) rho_fluid +
#               s_g rho_gas, with gas saturation s_g = min(Vgas / Vpore, 1) (Vgas_Vpore column) and
#               rho_gas the bubble's mass (sum of its g_ moles x molar mass) over its volume. s_g is
#               capped at 1: a bubble larger than the pores (pore_gas_fate = "retain") is not physical.

MINERAL_VOLUMES_CSV <- "data/mineral_molar_volumes.csv"
# Density of the remaining IOM (g/cm3). UNSOURCED, see README section 7 (Rock density):
# a measured density for chondritic IOM or kerogen is needed. 1.4 g/cm3 is a placeholder.
IOM_DENSITY_G_CM3 <- 1.4

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
solid_density <- function(moles, iom_g = 0, iom_density = IOM_DENSITY_G_CM3, volumes = read_mineral_volumes()) {
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
# water at any T and P, see README), else a number in g/cm3. Returns list(summary, breakdown).
bulk_density <- function(csv, porosity, pqi = sub("\\.csv$", ".pqi", csv), iom_density = IOM_DENSITY_G_CM3,
                         fluid_density = NULL, volumes = read_mineral_volumes()) {
  if (missing(porosity) || !is.numeric(porosity) || porosity < 0 || porosity >= 1)
    stop("porosity: give the model's porosity (Simulator$params$porosity)")
  init <- pqi_initial_solids(pqi)
  gfw <- element_gfw(init$database)
  iom_M <- vapply(init$iom$formula, kinetics_formula_mass, 0, gfw = gfw)
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
  add <- function(label, time, moles, iom_g, rho_f, s_g = 0, rho_g = 0) {
    s <- solid_density(moles, iom_g, iom_density, volumes)
    rho_p <- (1 - s_g) * rho_f + s_g * rho_g
    rows[[length(rows) + 1]] <<- data.frame(row = label, Time_Years = time, rho_solid = s$rho_solid,
                                             rho_fluid = rho_f, gas_saturation = s_g, rho_gas = rho_g,
                                             rho_bulk = (1 - porosity) * s$rho_solid + porosity * rho_p,
                                             solid_mass_g = s$mass_g, solid_volume_cm3 = s$volume_cm3, iom_g = iom_g)
    bd[[length(bd) + 1]] <<- cbind(row = label, s$breakdown)
  }
  rho_col <- if ("Rho_soln" %in% names(x)) "Rho_soln" else NULL
  add("t0", 0, init$minerals, init$iom_mass0_g, if (!is.null(fluid_density)) fluid_density else 1.0)  # pure water
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
    s_g <- if (vg > 0 && "Vgas_Vpore" %in% names(x)) min(as.numeric(x$Vgas_Vpore[r]), 1) else 0
    mg <- sum(vapply(grep("^g_", names(x), value = TRUE), function(k) {
      fe <- formula_elements(sub("\\(g\\)$", "", sub("^g_", "", k)))
      as.numeric(x[r, k]) * sum(gfw[names(fe)] * fe) }, 0))  # g
    rho_g <- if (vg > 0) mg / (vg * 1000) else 0           # g / cm3
    add(as.character(r), if ("Time_Years" %in% names(x)) as.numeric(x$Time_Years[r]) else NA_real_, mol, iom_g, rho_f, s_g, rho_g)
  }
  list(summary = do.call(rbind, rows), breakdown = do.call(rbind, bd))
}

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
