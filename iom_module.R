# =============================================================================
# iom_module.R  --  Standalone insoluble organic matter (IOM) decomposition
#                    module for ENTICES.
#
# PURPOSE OF THIS SEPARATION
#   The IOM treatment (channel definitions, activation energies, calibration
#   against Miller et al. 2025) is under active revision and is independent,
#   in principle, of the hydrotidal circulation physics (Darcy flow, tidal
#   exchange fractions, mineral kinetics). Keeping it in its own file means:
#     - editing IOM channels/Ea/composition cannot accidentally touch the
#       circulation code, and vice versa;
#     - this file can be unit-tested and validated against Miller's lab data
#       in complete isolation from ENTICES;
#     - swapping in a different organic-matter treatment later (e.g. a
#       different meteorite analogue, or decoupled-redox product species)
#       is a one-file change.
#
# INTERFACE CONTRACT (this is the entire surface ENTICES should depend on)
#   iom_default_config()       -> config list/data.frame (channel definitions)
#   iom_rates_block(config)    -> character string: PHREEQC RATES entries
#   iom_kinetics_block(cfg, iom_mass_kg) -> character string: KINETICS entries
#   iom_punch_headings(config) -> character string: USER_PUNCH heading fragment
#   iom_punch_lines(config, start_line) -> character string: USER_PUNCH PUNCH lines
#   iom_validate_config(config) -> TRUE or stops with a descriptive error
#
#   ENTICES calls these five functions and splices their return values into
#   its own template at the appropriate points (see INTEGRATION POINTS below).
#   ENTICES should not read or modify config$... fields directly beyond
#   passing the config object through; if a caller needs a different channel
#   set, it should build a whole replacement config via iom_default_config()
#   -> modify -> iom_validate_config(), not patch individual strings.
#
# STATUS: channel structure and calibration described in
#   PATCH_C_final_and_PATCH_E_redox.md. NOT YET VALIDATED against Miller's
#   own residue-trajectory data (Fig. 9) -- see iom_selftest.R companion
#   script, which runs this module standalone (no ENTICES, no circulation
#   physics) against the isothermal conditions Miller used.
# =============================================================================

# -----------------------------------------------------------------------------
# Configuration: the calibrated channel set
# -----------------------------------------------------------------------------
# Bulk bounded IOM composition (per kg IOM), used only for the balance check
# in iom_validate_config(), not otherwise referenced by ENTICES.
#   NOTE: N subscript here (2.142) comes from the gram composition (30 g N/kg
#   IOM; grams sum to exactly 1000 g/mol). An earlier code comment quoted
#   N_3.284 (46 g N) -- UNRESOLVED, confirm against source composition before
#   trusting nitrogen-channel absolute amounts.
iom_bulk_formula <- function() {
  c(C = 60.611, H = 45.639, O = 9.939, N = 2.142, S = 1.154)
}

#' Default IOM channel configuration.
#'
#' Each row is one decomposition channel: a first-order Arrhenius kinetic
#' reactant whose `-formula` is what the SOLID LOSES (not the product
#' molecule -- see PATCH_C_final_and_PATCH_E_redox.md sec 1 for why: Miller's
#' own O mass balance does not close without water contributing oxygen to the
#' measured CO2, so giving PHREEQC the solid's release stoichiometry lets it
#' draw the remainder from water and speciate at the ambient redox state).
#'
#' m0_per_kg is mol of channel-formula-units released per kg of IOM at full
#' exhaustion. The UNRELEASED remainder (bulk minus all channels) is treated
#' as inert char and is NOT part of this config -- it never appears in
#' PHREEQC because it never reacts. iom_validate_config() checks that no
#' element is over-drawn.
#'
#' Calibration: Miller et al. (2025) GCA 390:38-56, Tables 2/5/9, Murchison
#' IOM. Ea values fit directly against the ABSOLUTE m0_per_kg used here and
#' the ABSOLUTE measured 350C yield (48h, A = 1e13 /s), i.e. solving
#'   released(350C) = m0_per_kg * (1 - exp(-k*t)),  k = A*exp(-Ea/RT)
#' for Ea given released(350C) = Miller's measured 350C value and m0_per_kg
#' = the exhaustion amount used here (NOT Miller's raw 500C yield -- those
#' differ because the channel formulas carry H and O beyond the carbon Miller
#' measured, e.g. IOM_CO2's m0 of 4.169 vs Miller's measured CO2 yield of
#' 3.43 mol C/kg -- the extra amount is the H and O released alongside that
#' carbon). CORRECTED 2026-09: an earlier version fit Ea to the RATIO of
#' Miller's two raw measured yields (3.04/3.43) and then applied that Ea to
#' the larger, formula-scaled m0, which silently changed what fraction was
#' released -- verified by iom_selftest.R Part B against the absolute
#' 350C/500C values, not just the ratio. CO2: 214 -> 216.2 kJ/mol; CH4:
#' 228 -> 229.4 kJ/mol. See PATCH_C_final_and_PATCH_E_redox.md for the
#' original (now-superseded) derivation.
#' UPDATED 2026-09-17: IOM_N's Ea has now been through the same absolute-
#' value check. Miller's Table 9 gives NH4+ release for Murchison at 3 kbar
#' as both a direct mol/kg figure (umol N/mg sample == mol/kg) and as a %
#' of starting IOM nitrogen; converting the % via Table 1's Murchison N
#' content (2.48 wt%, i.e. 1.7705 mol N/kg IOM) reproduces the same absolute
#' values (0.568 mol/kg @350C, 0.723 mol/kg @500C) as the direct column,
#' which is a good cross-check on both tables. Refitting Ea_N against the
#' absolute 350C value (m0_per_kg = 0.875, unchanged) gives 217.3 kJ/mol
#' (was 215 kJ/mol placeholder). Predicted 500C release comes out at ratio
#' 1.21 vs Miller's measured value -- the same overshoot seen for CO2 (1.22)
#' and CH4 (1.22), for the identical reason (m0 is a fixed exhaustion pool,
#' not Miller's raw yield basis); this consistency across three
#' independently-fit channels is a check on the method, not a coincidence.
#'   IOM_S remains fully UNCALIBRATED -- Miller does not measure H2S yield
#'   at all. Its Ea is still just pinned to IOM_N's (now 217.3 kJ/mol) as a
#'   placeholder with no independent support; treat it accordingly.
#'
#' CAVEATS carried over verbatim from the original derivation:
#'   - two-point fit only; no 250C Murchison data exist, and single Ea per
#'     channel is known to under-predict low-temperature (~250C) CO2 release
#'     seen in Miller's syn-IOM samples. Low-Ea sub-pools are the planned fix
#'     for cold-Enceladus runs and are NOT yet in this config.
#'   - S channel is UNCALIBRATED (Miller does not measure H2S yield).
#' @return data.frame, one row per channel.
iom_default_config <- function() {
  data.frame(
    name      = c("IOM_CO2", "IOM_CH4", "IOM_N", "IOM_S"),
    m0_per_kg = c(4.169, 1.908, 0.875, 0.462),
    formula   = c("C 1 H 1 O 1.3115", "C 1 H 4.435", "N 1 H 2", "S 1 H 1"),
    Ea_J      = c(216.2e3, 229.4e3, 217.3e3, 217.3e3),
    logA      = c(13, 13, 13, 13),
    stringsAsFactors = FALSE
  )
}

#' Parse a "-formula" string ("C 1 H 1 O 1.3115") into a named numeric vector.
iom_parse_formula <- function(f) {
  toks <- strsplit(trimws(f), "\\s+")[[1]]
  els  <- toks[c(TRUE, FALSE)]
  amts <- as.numeric(toks[c(FALSE, TRUE)])
  stats::setNames(amts, els)
}

#' Validate a channel config: structural checks + element balance sanity
#' (no channel, individually or combined, may draw more of any element than
#' the bulk IOM formula contains). Called by ENTICES before generation, and
#' by the self-test script.
#' @param config data.frame as returned by iom_default_config()
#' @return TRUE (invisibly) on success; stops with a descriptive message on
#'   any failure.
iom_validate_config <- function(config) {
  req_cols <- c("name", "m0_per_kg", "formula", "Ea_J", "logA")
  missing_cols <- setdiff(req_cols, names(config))
  if (length(missing_cols) > 0)
    stop("iom_validate_config: config missing column(s): ",
         paste(missing_cols, collapse = ", "))
  if (nrow(config) == 0) stop("iom_validate_config: config has no channels")
  if (any(duplicated(config$name)))
    stop("iom_validate_config: duplicate channel name(s): ",
         paste(config$name[duplicated(config$name)], collapse = ", "))
  if (any(config$m0_per_kg < 0)) stop("iom_validate_config: negative m0_per_kg")
  if (any(config$Ea_J <= 0)) stop("iom_validate_config: non-positive Ea_J")

  bulk <- iom_bulk_formula()
  released <- stats::setNames(rep(0, length(bulk)), names(bulk))
  for (i in seq_len(nrow(config))) {
    f <- tryCatch(iom_parse_formula(config$formula[i]), error = function(e)
      stop("iom_validate_config: could not parse formula for channel '",
           config$name[i], "': ", config$formula[i]))
    unknown_els <- setdiff(names(f), names(bulk))
    if (length(unknown_els) > 0)
      stop("iom_validate_config: channel '", config$name[i],
           "' formula has element(s) not in bulk IOM formula: ",
           paste(unknown_els, collapse = ", "))
    for (el in names(f)) {
      released[el] <- released[el] + config$m0_per_kg[i] * f[el]
    }
  }
  over <- released > bulk + 1e-9
  if (any(over))
    stop("iom_validate_config: channel(s) over-draw element(s) beyond bulk IOM: ",
         paste(sprintf("%s (released %.4f > bulk %.4f)",
                       names(bulk)[over], released[over], bulk[over]),
               collapse = "; "),
         ". Residue would be negative -- check m0_per_kg / formulas.")
  invisible(TRUE)
}

#' Implied inert residue composition (bulk minus all channels at full
#' exhaustion). Diagnostic only -- not used in PHREEQC generation, since the
#' residue never reacts and therefore never appears in a RATES/KINETICS block.
#' @return named numeric vector, same elements as iom_bulk_formula().
iom_residue_formula <- function(config) {
  iom_validate_config(config)
  bulk <- iom_bulk_formula()
  released <- stats::setNames(rep(0, length(bulk)), names(bulk))
  for (i in seq_len(nrow(config))) {
    f <- iom_parse_formula(config$formula[i])
    for (el in names(f)) released[el] <- released[el] + config$m0_per_kg[i] * f[el]
  }
  bulk - released
}

# -----------------------------------------------------------------------------
# PHREEQC text generation (the actual interface ENTICES uses)
# -----------------------------------------------------------------------------

#' PHREEQC RATES entries for every channel: first-order Arrhenius decay,
#' rate = 10^logA * exp(-Ea/(R*TK)) * M, where M is the remaining moles of
#' that channel's KINETICS reactant (supplied by PHREEQC as `M` inside a
#' RATES Basic program).
#' @param config data.frame as returned by iom_default_config()
#' @return character(1): the RATES entries, ready to splice into a PHREEQC
#'   RATES block (after any mineral RATES entries; order does not matter to
#'   PHREEQC).
iom_rates_block <- function(config) {
  iom_validate_config(config)
  paste0(apply(config, 1, function(p) {
    sprintf(paste0(
      "%s\n -start\n",
      " 1  REM Arrhenius 1st-order IOM decomposition channel\n",
      " 2  REM formula = what the SOLID loses; see PATCH_C_final doc\n",
      " 10 k = 10^(%s) * exp(-%s / (8.314 * TK))\n",
      " 20 rate = k * M\n",
      " 30 if (M <= 0) then rate = 0\n",
      " 40 moles = rate * TIME\n",
      " 50 SAVE moles\n -end\n"),
      p[["name"]], p[["logA"]], p[["Ea_J"]])
  }), collapse = "")
}

#' PHREEQC KINETICS entries for every channel, scaled to the IOM mass present.
#' @param config data.frame as returned by iom_default_config()
#' @param iom_mass_kg numeric(1): kg of IOM in the reacting system (i.e.
#'   (organic_wt_percent/100) * total_rock, however ENTICES's circulation
#'   physics defines total_rock -- this module does not know or care what
#'   that normalization is, it only multiplies by whatever kg value it is given)
#' @return character(1): the KINETICS entries, ready to splice inside a
#'   PHREEQC KINETICS block (after any mineral entries; -steps/-cvode
#'   settings on the enclosing block apply to these too).
iom_kinetics_block <- function(config, iom_mass_kg) {
  iom_validate_config(config)
  if (!is.numeric(iom_mass_kg) || length(iom_mass_kg) != 1 || iom_mass_kg < 0)
    stop("iom_kinetics_block: iom_mass_kg must be a single non-negative number")
  paste0(apply(config, 1, function(p) {
    m0 <- as.numeric(p[["m0_per_kg"]]) * iom_mass_kg
    sprintf("   %s\n      -formula %s\n      -m0 %.6e\n      -step_divide 1\n",
            p[["name"]], trimws(p[["formula"]]), m0)
  }), collapse = "")
}

#' USER_PUNCH heading fragment for the IOM channels (one "<name>_mol" per
#' channel), to be appended to whatever headings ENTICES already punches.
#' @param config data.frame as returned by iom_default_config()
#' @return character(1), e.g. " IOM_CO2_mol IOM_CH4_mol IOM_N_mol IOM_S_mol"
#'   (leading space so it concatenates cleanly onto an existing heading list).
iom_punch_headings <- function(config) {
  iom_validate_config(config)
  paste0(" ", paste(sprintf("%s_mol", config$name), collapse = " "))
}

#' USER_PUNCH PUNCH-statement lines for the IOM channels, i.e.
#' "   210 PUNCH KIN(\"IOM_CO2\")\n" etc. Line numbers are auto-assigned
#' starting at start_line in steps of 10, so callers should pass a start_line
#' safely above whatever line numbers their own PUNCH statements already use.
#' @param config data.frame as returned by iom_default_config()
#' @param start_line integer(1): first Basic line number to use (default 500,
#'   chosen to sit well above typical mineral PUNCH blocks; ENTICES should
#'   pass an explicit value if its own numbering could reach that high).
#' @return character(1): the PUNCH lines, ready to splice into an existing
#'   USER_PUNCH -start ... -end block.
iom_punch_lines <- function(config, start_line = 500) {
  iom_validate_config(config)
  if (!is.numeric(start_line) || length(start_line) != 1 || start_line %% 1 != 0)
    stop("iom_punch_lines: start_line must be a single integer")
  paste0(vapply(seq_len(nrow(config)), function(j) {
    sprintf("   %d PUNCH KIN(\"%s\")\n", start_line + 10 * (j - 1), config$name[j])
  }, character(1)), collapse = "")
}

# -----------------------------------------------------------------------------
# INTEGRATION POINTS (for whoever wires this into ENTICES)
# -----------------------------------------------------------------------------
# In the ENTICES template generator, replacing the current inline IOM block:
#
#   cfg <- iom_default_config()               # or a caller-supplied variant
#   iom_mass_kg <- (self$organic_wt_percent / 100) * self$phreeqc_params$total_rock
#
#   iom_rates    <- if (self$organic_wt_percent > 0) iom_rates_block(cfg)    else ""
#   iom_kinetics <- if (self$organic_wt_percent > 0) iom_kinetics_block(cfg, iom_mass_kg) else ""
#   iom_headings <- if (self$organic_wt_percent > 0) iom_punch_headings(cfg) else ""
#   iom_punch    <- if (self$organic_wt_percent > 0) iom_punch_lines(cfg, start_line = 200) else ""
#
#   Splice iom_rates at the end of the RATES section (after mineral entries).
#   Splice iom_kinetics inside KINETICS 1, after mineral entries, before
#     -steps / -cvode settings.
#   Append iom_headings to the existing USER_PUNCH -headings line.
#   Splice iom_punch inside USER_PUNCH -start ... -end, after the mineral
#     PUNCH statements and before -end. Choose start_line above ENTICES's own
#     highest PUNCH line number (currently ~110 in the reviewed template, so
#     a default of 200 leaves headroom for a few dozen mineral additions
#     before colliding -- confirm against the current file before wiring in).
#
# ENTICES's restart function should read the config with iom_default_config()
# (or store/reload whatever config a run used) and use the same three
# generator calls, feeding -m0 from the recovered "<name>_mol" CSV columns
# instead of iom_mass_kg -- i.e. build a config with m0_per_kg replaced by
# per-channel absolute recovered moles and call iom_kinetics_block(recovered_cfg,
# iom_mass_kg = 1) (mass_kg = 1 makes m0_per_kg act as the absolute value).
