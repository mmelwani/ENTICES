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
# STATUS (2026-09-18): channel structure is now distributed-Ea sub-pools for
#   CO2/CH4/CHn (single-Ea for N/S), calibrated per PATCH_F_vitrimat2018.md
#   and PATCH_F_addendum_N_and_CHn.md -- see iom_default_config()'s docstring
#   for the full derivation and history. Validated so far against Miller's
#   absolute Murchison gas yields (350C, 500C, both now exact fit points) and
#   the 350C residue trajectory (independent of the fit; see PROGRESS.md
#   Part 4). NOT YET validated against ENTICES/PHREEQC itself, and Option 3
#   (coupled/decoupled redox, PROGRESS.md Part 3) is not implemented -- this
#   config's -formula values are coupled-mode only. See iom_selftest.R, which
#   runs this module standalone (no ENTICES, no PHREEQC, no circulation
#   physics) against the isothermal conditions Miller used.
# =============================================================================

# -----------------------------------------------------------------------------
# Configuration: the calibrated channel set
# -----------------------------------------------------------------------------
# Bulk bounded IOM composition (per kg IOM), used only for the balance check
# in iom_validate_config(), not otherwise referenced by ENTICES.
#   NOTE (RESOLVED 2026-09-18, PATCH_F_addendum_N_and_CHn.md sec 1): N = 2.142
#   is correct. Cody et al. (2024), via Miller Table 1, gives Murchison N/C =
#   0.036 +/- (2.48 +/- 0.04 wt% N), which on this normalized 1000 g/mol basis
#   implies N = 2.182 -- within 1.8% of 2.142, inside Cody's own uncertainty.
#   The old rival figure (N_3.284, "46 g N/kg") was traced to a units mixup:
#   46 g is this formula's HYDROGEN gram value (H 45.639 mol x 1.0079 g/mol =
#   46.0 g), misread as nitrogen. 2.142 stays; 2.182 would be a defensible
#   sub-2%-change tidy-up for exact Cody-consistency, not a fix, and has not
#   been applied.
iom_bulk_formula <- function() {
  c(C = 60.611, H = 45.639, O = 9.939, N = 2.142, S = 1.154)
}

#' Default IOM channel configuration.
#'
#' Each row is one decomposition channel: a first-order Arrhenius kinetic
#' reactant whose `-formula` is what the SOLID LOSES (not the product
#' molecule -- see PATCH_C_final_and_PATCH_E_redox.md sec 1 / PROGRESS.md
#' Part 1 for why: Miller's own O mass balance does not close without water
#' contributing oxygen to the measured CO2, so giving PHREEQC the solid's
#' release stoichiometry lets it draw the remainder from water and speciate
#' at the ambient redox state).
#'
#' m0_per_kg is mol of channel-formula-units released per kg of IOM at full
#' exhaustion. The UNRELEASED remainder (bulk minus all channels) is treated
#' as inert char and is NOT part of this config -- it never appears in
#' PHREEQC because it never reacts. iom_validate_config() checks that no
#' element is over-drawn.
#'
#' HISTORY (kept because each stage found a real, previously-undetected
#' error -- see PROGRESS.md and the PATCH_F docs for full derivations):
#'   1. Original single-Ea fit (PATCH_C): CO2/CH4 Ea fit to the RATIO of
#'      Miller's two measured yields (3.04/3.43), applied to a differently-
#'      scaled m0 -- an 18% error, caught by iom_selftest.R Part B.
#'   2. CORRECTED (2026-09): refit CO2/CH4 to ABSOLUTE 350C yields instead of
#'      the ratio, at A=1e13. CO2: 214 -> 216.2 kJ/mol; CH4: 228 -> 229.4.
#'      IOM_N refit the same way: 215 -> 217.3 kJ/mol. All still single-Ea,
#'      so 500C was a PREDICTION, and overshot Miller's measured yield by
#'      ~22% for all three channels -- attributed to m0 being a fixed
#'      exhaustion pool rather than Miller's raw yield basis.
#'   3. SUPERSEDED (2026-09-18, PATCH_F_vitrimat2018.md): single-Ea was wrong
#'      by 3-5 ORDERS OF MAGNITUDE at Enceladus temperatures (0-100C), not a
#'      refinement issue. Reason: our only calibration data are 48h at 350C/
#'      500C, and at those temperatures a real activation-energy distribution's
#'      low-Ea tail is already fully converted -- so single-Ea data contain
#'      almost no information about exactly the tail that controls cold-case
#'      release. Also A was wrong: 1e13 is the 1989 Vitrimat value; Burnham
#'      (2019) Table 1 gives A = 2e15 /s ("Vitrimat 2018"), explicitly the
#'      revision meant for hydrous pyrolysis (which is what Miller ran).
#'
#' CURRENT SCHEME (2026-09-18): CO2, CH4 and CHn are each split into several
#' sub-pool rows sharing one -formula but each with its own Ea_J, at a
#' logA = 15.301 (A = 2e15 /s) common to every row in this config (N and S
#' included). Sub-pool WEIGHTS (fraction of the channel's total m0 in each
#' Ea bin) are Burnham (2019) Table 1's Vitrimat-2018 distributions --
#' fitted to vitrinite, not IOM, but validated below. Each channel's total
#' m0 and a uniform kcal/mol SHIFT applied to every bin are fit to Miller's
#' two absolute Murchison measured yields (350C, 500C) by 2-point nonlinear
#' solve (shift given the ratio of the two yields; m0 then follows). Unlike
#' the old single-Ea fit, THIS makes both 350C and 500C exact fit points, so
#' the old ~22% 500C overshoot for CO2/CH4 is resolved by construction, not
#' just explained.
#'   IOM_CO2: shift +3.587 kcal/mol, m0 4.1687 mol/kg (ENTICES basis).
#'   IOM_CH4: shift +3.736 kcal/mol, m0 2.2760 mol/kg (ENTICES basis).
#'   Both independently re-derived (not copied from the patch doc) in a
#'   2026-09-18 VS Code session; agreed with PATCH_F_vitrimat2018.md's
#'   hand-computed values to within 0.2%.
#'
#' VALIDATION (the actual reason to trust the distribution shape at all):
#' the UNSHIFTED Vitrimat-2018 CO2 distribution predicts the ratio between
#' Miller's HC113 250C and 500C yields to 0.991 (see iom_selftest.R Part B),
#' with zero fitting to our data. That is the strongest evidence this module
#' has for the shape (not the shift/m0, which ARE fit to Murchison).
#'
#' IOM_CHn is new: Burnham/Vitrimat's "oil" (C2+) channel, previously omitted
#' entirely (silently assuming zero C2+ production). Miller's GC does not
#' report C2+ species, so CHn's m0 CANNOT be calibrated against Murchison
#' data at all -- unlike CO2/CH4/N, it is a straight transfer from Burnham's
#' own vitrinite stoichiometry (c(oil) = 2% of total C, formula C 1 H 1.8),
#' using the CO2 channel's fitted shift for lack of any independent CHn
#' constraint. Treat its epistemic status as strictly weaker than CO2/CH4/N.
#'
#' IOM_N stays single-Ea (no distribution data exists for it): m0 unchanged
#' at 0.875 mol/kg (ENTICES basis). Ea = 242.8 kJ/mol.
#'
#' CORRECTED 2026-09-18 (caught while re-deriving this for the new A, not
#' inherited from any patch doc -- both the ORIGINAL 217.3 kJ/mol fit at
#' A=1e13 and my own first attempt at refitting it to 244.8 kJ/mol at the
#' new A shared the same latent bug): m0_per_kg is on the ENTICES-normalized
#' bulk basis (bulk N = 2.142 mol/kg), which is NOT the same basis as
#' Miller's raw Murchison measurement (real Murchison N = 1.7705 mol/kg,
#' Table 1) -- they differ by a basis scale factor of 2.142/1.7705 = 1.2098,
#' analogous to (and independently derived from) the carbon-basis scale
#' factor of 1.2154 used for CO2/CH4/CHn above. Both fits computed
#' "target_fraction = Miller's raw 350C value / m0_per_kg" directly, silently
#' dividing a Murchison-basis absolute value by an ENTICES-basis pool size --
#' dimensionally inconsistent, and it forced Ea to secretly absorb a ~21%
#' basis-conversion factor that has nothing to do with the actual activation
#' energy. This means the WIDELY-REPEATED "~20-22% overshoot at 500C, because
#' m0 is a fixed exhaustion pool not on Miller's raw yield basis" explanation
#' recorded in this file's history above and in PROGRESS.md/HANDOVER docs is
#' the right SIGN but the wrong MECHANISM for at least the IOM_N case: it
#' isn't an intentional, benign consequence of pool sizing, it is this basis
#' bug, caught only by explicitly checking the arithmetic before reusing it
#' at the new A rather than just re-deriving the old number at a new A.
#'   Correct method: since m0=0.875 is calibrated so that 500C/48h
#'   approximates full exhaustion of the reactive N pool (via Murchison's own
#'   raw 500C measurement, scaled up), the physically meaningful fit target
#'   is the DIMENSIONLESS ratio of Miller's two raw measurements,
#'   0.5676/0.7229 = 0.7852, matched to (1-exp(-k*t)) -- not an absolute
#'   value divided by a differently-scaled m0. That gives Ea_N = 242.8
#'   kJ/mol (not 244.8). Checked: predicted/n_scale reproduces Miller's raw
#'   350C and 500C values to within 0.06%.
#'   IOM_CO2/IOM_CH4/IOM_CHn above do NOT share this bug -- their shift/m0
#'   fit already used the ratio of Miller's two raw measurements to solve
#'   for shift first (basis-independent by construction), then scaled the
#'   resulting Murchison-basis m0 up to ENTICES basis as a separate,
#'   dimensionally-clean step. Re-verify this claim rather than trust it
#'   before extending the pattern to any new channel.
#'   IOM_S remains fully UNCALIBRATED -- Miller does not measure H2S yield at
#'   all. Its Ea is just pinned to IOM_N's (242.8 kJ/mol) as a placeholder
#'   with zero independent support; treat it accordingly in any
#'   sulfur-dependent output.
#'
#' STILL OPEN (see PROGRESS.md Part 4/5 and the PATCH_F docs for detail):
#'   - Murchison's fitted shift (+3.6) and HC113's own gas-yield fit (~0.0)
#'     disagree -- most likely real syn-IOM-vs-meteoritic material
#'     difference (corroborated by Miller's own Table 10 needing HC113 +8
#'     kcal/mol vs BS89), but NOT reconciled with Miller's own Table 10
#'     value for HC113, which is a residue fit at geologic (10C/Myr) heating
#'     rather than our 48h isothermal gas-yield fit. Do not average the two;
#'     Murchison is used for production here, HC113 is a separate upper
#'     bound, not merged in.
#'   - CHn's m0 (2% of C, from coal) and its borrowed shift are unconstrained
#'     by any IOM-specific measurement.
#'   - Option 3 (coupled/decoupled redox at low temperature): design settled
#'     in OPTION3_decision.md, formulas TESTED IN PHREEQC 2026-09-20
#'     (PHREEQC_TEST_RESULTS_2026-09-20.md). Only IOM_CH4 and IOM_S swap
#'     formula; IOM_N deliberately stays coupled; IOM_CO2/IOM_CHn decouple
#'     via Hdg at the solution level. See DECOUPLED_FORMULAS below for the
#'     measured numbers and the Mtg-vs-Sg asymmetry that makes their two
#'     formulas follow opposite rules.
#'     KNOWN LIMITATION, unchanged: iom_validate_config() cannot validate a
#'     decoupled config -- its element-balance check tokenizes -formula
#'     strings expecting real periodic-table symbols, so "Mtg"/"Sg" are
#'     rejected as unknown elements (iom_selftest.R Part C asserts this). It
#'     is left as an honest failure rather than quietly taught to recognise
#'     pseudo-species, because the balance it would then report would be
#'     meaningless: Mtg's hydrogens are outside the H mass balance by
#'     construction, so "elements released" is not comparable between modes
#'     as a config-string property. The meaningful check is the PHREEQC one,
#'     now done -- see the test results doc.
#' @param redox_mode named character vector, e.g. c(C="coupled", N="coupled",
#'   S="coupled") (the default -- fully backward compatible, identical output
#'   to calling with no argument). Each element independently "coupled" or
#'   "decoupled":
#'     S="decoupled"  -> IOM_S switches to "Sg 1 H 1".
#'     C="decoupled"  -> IOM_CH4 switches to "Mtg 1 H 0.435"; IOM_CO2 and
#'                       IOM_CHn keep their coupled formulas, because their
#'                       decoupling is achieved by the CALLER entering the H2
#'                       inventory as Hdg in the PHREEQC solution, which no
#'                       -formula can express. attr(cfg,
#'                       "carbon_decoupling_note") records this reminder.
#'     N="decoupled"  -> refused with an error. Not a gap: NH3/NH4+ is
#'                       already the stable N species at Enceladus pore
#'                       conditions, and the only decoupled N species in the
#'                       database (Ntg = N2) is the wrong end of the ladder.
#' @return data.frame, one row per sub-pool/channel. Carries formula (the
#'   ACTIVE one, selected by redox_mode -- this is what every other function
#'   in this file reads) plus formula_coupled/formula_decoupled for
#'   reference/switching.
iom_default_config <- function(redox_mode = c(C = "coupled", N = "coupled", S = "coupled")) {
  R_GAS <- 8.314
  CAL_TO_J <- 4184
  logA <- log10(2e15)  # 15.301; Burnham (2019) Table 1, "Vitrimat 2018", A=2e15/s

  # Vitrimat-2018 sub-pool shapes (Burnham 2019 Table 1), Ea bins in kcal/mol
  # BEFORE the Murchison shift below is added.
  co2_bins <- c(44, 46, 48, 50, 52, 54, 56); co2_w <- c(10, 15, 15, 15, 15, 15, 15) / 100
  ch4_bins <- c(52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76)
  ch4_w    <- c(2, 5, 8, 10, 12, 15, 12, 10, 8, 6, 5, 4, 3) / 100
  chn_bins <- c(50, 52, 54, 56, 58, 60, 62); chn_w <- c(5, 15, 30, 20, 15, 10, 5) / 100
  stopifnot(abs(sum(co2_w) - 1) < 1e-9, abs(sum(ch4_w) - 1) < 1e-9, abs(sum(chn_w) - 1) < 1e-9)

  # Fitted shift + total m0 (ENTICES basis), 2026-09-18 -- see the docstring
  # above for how, and iom_selftest.R Part B for the independent re-derivation.
  co2_shift <- 3.587; co2_m0_tot <- 4.1687
  ch4_shift <- 3.736; ch4_m0_tot <- 2.2760
  chn_shift <- co2_shift  # borrowed; no independent CHn constraint
  chn_m0_tot <- 0.02 * iom_bulk_formula()[["C"]]  # Burnham's c(oil) = 2% of C

  # DECOUPLED-mode formulas -- TESTED IN PHREEQC 2026-09-20. Full results and
  # measured numbers in PHREEQC_TEST_RESULTS_2026-09-20.md; design rationale
  # in OPTION3_decision.md. Only two channels get a decoupled formula:
  #
  #   IOM_CH4 -> "Mtg 1 H 0.435"  VERIFIED CORRECT. Reproduces the coupled
  #     formula's aqueous state to every printed digit: m(H2) 7.175e-04,
  #     pH 11.5000, pe -11.425, all identical to "C 1 H 4.435". Writing
  #     "Mtg 1 H 4.435" instead would double-count methane's four hydrogens
  #     and give 3.8x too much H2.
  #
  #   IOM_S -> "Sg 1 H 1"  CORRECTED from OPTION3_decision.md's "Sg 1 H -1",
  #     which was wrong by 2 mol H per mol S. Measured H contributed per mole:
  #     coupled "S 1 H 1" = +1.000; "Sg 1 H 1" = +1.000 (match); "Sg 1 H -1"
  #     = -1.000 (the solid would ABSORB an H instead of releasing one);
  #     "Sg 1" = 0.000. The bad variant drives pe to +15.4 and destroys the
  #     H2 reservoir.
  #
  # WHY THE TWO DIFFER (this is the non-obvious part, and what the design doc
  # got wrong): the database defines their master species differently --
  #     Mtg   Mtg    0   Mtg     16.032   <- master species has NO H
  #     Sg    H2Sg   1   H2Sg    34.08    <- master species IS H2Sg
  # so Mtg's four hydrogens sit OUTSIDE the H mass balance (supply only the
  # surplus H), while H2Sg's two sit INSIDE it and are drawn from solution
  # (supply the same H the coupled formula does). Same-looking species,
  # opposite treatment. Do not "make them consistent."
  #
  # NOT decoupled, deliberately:
  #   IOM_N   -- stays coupled. At Enceladus pore conditions NH3/NH4+ is
  #     already the stable N species, so coupled equilibrium moves released
  #     N(-2) to N(-3): one electron, the direction kinetics would take it
  #     anyway. There is no artefact to fix, and the only decoupled N species
  #     available (Ntg = N2) is the WRONG end of the ladder. See
  #     OPTION3_decision.md sec 1.
  #   IOM_CO2 / IOM_CHn -- keep their coupled formulas. Carbon decoupling is
  #     achieved at the SOLUTION/DATABASE level by entering H2 as Hdg, not by
  #     a carbon species (there is no decoupled oxidised-carbon species, and
  #     adding one would silently break carbonate-mineral equilibria).
  #     Tested and confirmed: with Hdg, carbon reduction drops ~26,000x
  #     (10.8% -> 0.0004% of C converted to methane) while pe stays sensibly
  #     set by the mineral assemblage (-8.26) and calcite still precipitates
  #     normally. Because this is a solution-level switch it is NOT expressed
  #     in any -formula here -- ENTICES must enter its H2 inventory as Hdg.
  #     CAVEAT, measured: decoupling H2 is not redox-neutral for the rest of
  #     the system -- magnetite precipitation ~2.9x and pyrrhotite dissolution
  #     ~3.6x higher than the coupled case. Report the mineral assemblage
  #     alongside any coupled/decoupled carbon comparison.
  DECOUPLED_FORMULAS <- list(
    IOM_CH4 = "Mtg 1 H 0.435",
    IOM_S   = "Sg 1 H 1"
  )

  subpool_rows <- function(prefix, bins, w, shift, m0_tot, formula_coupled, family) {
    data.frame(
      name            = sprintf("%s_%d", prefix, bins),
      m0_per_kg       = w * m0_tot,
      formula_coupled = formula_coupled,
      formula_decoupled = if (is.null(DECOUPLED_FORMULAS[[prefix]])) NA_character_ else DECOUPLED_FORMULAS[[prefix]],
      family          = family,
      Ea_J            = (bins + shift) * CAL_TO_J,
      logA            = logA,
      stringsAsFactors = FALSE
    )
  }

  cfg <- rbind(
    subpool_rows("IOM_CO2", co2_bins, co2_w, co2_shift, co2_m0_tot, "C 1 H 1 O 1.3115", "C"),
    subpool_rows("IOM_CH4", ch4_bins, ch4_w, ch4_shift, ch4_m0_tot, "C 1 H 4.435", "C"),
    subpool_rows("IOM_CHn", chn_bins, chn_w, chn_shift, chn_m0_tot, "C 1 H 1.8", "C"),
    data.frame(name = "IOM_N", m0_per_kg = 0.875, formula_coupled = "N 1 H 2",
               formula_decoupled = NA_character_, family = "N",
               Ea_J = 242.8e3, logA = logA, stringsAsFactors = FALSE),
    data.frame(name = "IOM_S", m0_per_kg = 0.462, formula_coupled = "S 1 H 1",
               formula_decoupled = DECOUPLED_FORMULAS[["IOM_S"]], family = "S",
               Ea_J = 242.8e3, logA = logA, stringsAsFactors = FALSE)
  )

  # --- select the ACTIVE formula per redox_mode (default: all coupled, i.e.
  # identical behavior to before redox_mode existed) ---
  allowed <- c("coupled", "decoupled")
  for (el in names(redox_mode)) {
    if (!redox_mode[[el]] %in% allowed)
      stop("iom_default_config: redox_mode['", el, "'] must be 'coupled' or 'decoupled'")
  }
  # Nitrogen decoupling is a decision made AGAINST, not a gap: the only
  # decoupled N species in the database is Ntg (N2), which is the wrong end
  # of the redox ladder for the reduced N this channel releases. Coupled
  # equilibrium already puts it where kinetics would. Refuse rather than
  # silently doing nothing.
  if (identical(redox_mode[["N"]], "decoupled"))
    stop("iom_default_config: redox_mode N='decoupled' is deliberately not ",
         "supported. NH3/NH4+ is already the stable N species at Enceladus ",
         "pore conditions, so coupled equilibrium gives the right answer; the ",
         "only decoupled N species available (Ntg = N2) would oxidise the ",
         "released nitrogen away from it. See OPTION3_decision.md sec 1.")
  # A row keeps its coupled formula unless a decoupled one is defined for it.
  # That fallback is intentional, not a gap:
  #   IOM_CO2 / IOM_CHn -- decoupled via Hdg at the solution level (the
  #     caller must enter H2 as Hdg); no -formula change is involved.
  #   IOM_N             -- deliberately stays coupled (rejected above).
  # Only IOM_CH4 (Mtg) and IOM_S (Sg) actually swap formula.
  use_decoupled <- redox_mode[cfg$family] == "decoupled" & !is.na(cfg$formula_decoupled)
  cfg$formula <- ifelse(use_decoupled, cfg$formula_decoupled, cfg$formula_coupled)
  attr(cfg, "redox_mode") <- redox_mode
  attr(cfg, "carbon_decoupling_note") <-
    if (identical(redox_mode[["C"]], "decoupled"))
      paste("Carbon decoupling is a SOLUTION-level switch: the caller must enter",
            "the H2 inventory as Hdg. No -formula in this config encodes it.",
            "IOM_CH4 has switched to Mtg; IOM_CO2/IOM_CHn keep coupled formulas.")
    else NA_character_
  cfg
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
#     CAUTION (2026-09-18): the distributed-Ea config now has 29 rows, not 4
#     -- iom_punch_lines() at step 10 needs ~290 lines of Basic line-number
#     space (e.g. 200..3080 rather than 200..230). Re-check this headroom
#     claim against ENTICES's current template before wiring in; it may no
#     longer hold if ENTICES's own PUNCH numbering has grown since it was
#     last checked.
#
# ENTICES's restart function should read the config with iom_default_config()
# (or store/reload whatever config a run used) and use the same three
# generator calls, feeding -m0 from the recovered "<name>_mol" CSV columns
# instead of iom_mass_kg -- i.e. build a config with m0_per_kg replaced by
# per-channel absolute recovered moles and call iom_kinetics_block(recovered_cfg,
# iom_mass_kg = 1) (mass_kg = 1 makes m0_per_kg act as the absolute value).
