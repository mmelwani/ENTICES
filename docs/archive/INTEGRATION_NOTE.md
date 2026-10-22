# IOM module split — what this is and what it found

## What you asked for

Separate the IOM (organic decomposition) treatment from the rest of ENTICES,
so hydrotidal-circulation edits can't accidentally touch organics and
organics edits can't accidentally touch circulation.

## What's here

- **iom_module.R** — the standalone module. Interface is deliberately
  minimal, per your choice: five functions.
    - `iom_default_config()` — the channel definitions (data.frame)
    - `iom_rates_block(config)` — PHREEQC RATES text
    - `iom_kinetics_block(config, iom_mass_kg)` — PHREEQC KINETICS text
    - `iom_punch_headings(config)` / `iom_punch_lines(config, start_line)`
      — USER_PUNCH text
    - `iom_validate_config(config)` — structural + element-balance checks,
      called internally by all of the above so a bad config fails loudly
      wherever it's used, not silently downstream in a PHREEQC error.

  ENTICES should treat this as a black box: build or receive a config,
  pass it to these four generator functions, splice the returned strings
  into its own template. It should not reach into `config$...` fields
  directly for anything except passing the object through. The
  "INTEGRATION POINTS" comment block at the bottom of the file spells out
  exactly where each string needs to be spliced into the current template,
  including a specific caution about USER_PUNCH line numbers colliding
  with the mineral PUNCH statements already there.

- **iom_selftest.R** — runs the module with no ENTICES, no PHREEQC, no
  circulation physics. Two parts: (A) fast structural checks (config
  validity, balance arithmetic, string generation — catches things like a
  channel that draws more of an element than the bulk IOM contains); (B) a
  from-scratch reimplementation of the Arrhenius rate law in R, integrated
  at Miller's own experimental conditions, to check whether the calibrated
  channels reproduce his measured yields.

## What Part B found before you even ran it

I hand-traced Part B's arithmetic to make sure the script itself was
correct before handing it to you, and that trace surfaced a real
calibration bug in the numbers that were previously written up as final
(in `PATCH_C_final_and_PATCH_E_redox.md`):

The Ea values for IOM_CO2 and IOM_CH4 had been fit to match the *ratio*
Miller measured between 350C and 500C (3.04/3.43 for CO2), then that same
Ea was applied to `m0_per_kg` values that are scaled to bulk IOM elemental
content, not to Miller's raw carbon-only yields. Those two m0 bases differ
because the channel formulas carry extra H and O beyond the carbon Miller
measured (e.g. IOM_CO2 is `C 1 H 1 O 1.3115`, not just `C 1`). Applying a
ratio-fitted Ea to a differently-scaled m0 silently changes what fraction
gets released. Concretely: predicted CO2 at 350C came out to 3.60 mol/kg
against the 3.04 the fit was supposed to reproduce — an 18% error, in the
one number we'd already told you was calibrated.

**Fixed:** refit both Ea values directly against the absolute calibrated
m0 (4.169, 1.908) and the absolute measured 350C yield, rather than a
ratio. CO2: 214 -> 216.2 kJ/mol. CH4: 228 -> 229.4 kJ/mol. Confirmed by
hand (Python) that this reproduces Miller's 350C values almost exactly
(3.034 vs 3.04, 0.1847 vs 0.184) before writing it into the module. This
is now the version in `iom_default_config()`, with the correction noted
inline and dated.

Two things this doesn't fix, flagged in the code:
- IOM_N and IOM_S Ea were not independently re-checked this way — Miller
  reports NH3 as a percentage of IOM nitrogen, not directly as mol/kg, so
  the same absolute-value check needs a slightly different derivation.
  Treat those two as less trustworthy than IOM_CO2/IOM_CH4 for now.
- The residue-trajectory comparison (Part B's second half, against
  Miller's Table 5 H/C and O/C) still has placeholder `NA` values — I
  didn't have exact Table 5 Murchison 350C/500C figures pinned down
  precisely enough to hardcode with confidence in this pass. The script
  will run and print predicted values regardless; you or Lucas should
  fill in the measured `H_C`/`O_C` in `miller_residue` from the paper (or
  I can go pull them properly next).

## How to run it

```
cd wherever/you/put/both/files
Rscript iom_selftest.R
```

No packages beyond base R. Should take well under a second.

## What's NOT done

- The module hasn't been wired into `ENTICES_LF_patched.R` (or whatever
  Lucas's current file is) — you chose to design the interface first and
  apply it once his current file is confirmed. The "INTEGRATION POINTS"
  section in `iom_module.R` is what that wiring should follow.
- Low-Ea sub-pools for cold-temperature CO2 release (the gap Miller's own
  syn-IOM 250C data expose) are still not in the config.
- The residue-trajectory validation is scaffolded but not filled in.
- IOM_N/IOM_S calibration still needs the same absolute-value fix applied
  to IOM_CO2/IOM_CH4.

## Suggested next step

Run `iom_selftest.R` yourself once you confirm base R has no surprises in
your environment, then let's fill in the real Table 5 numbers and pull the
IOM_N/IOM_S fix the same way — that closes out the validation this module
was specifically built to make possible.
