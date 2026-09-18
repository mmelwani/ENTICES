# Handover: IOM module — first execution, IOM_N recalibration, benchmark hardening

Date: 2026-09-17
Branch: `iom-module` (created off `main`; not yet merged)
Prior context: `INTEGRATION_NOTE.md` (module split rationale) and a longer
handover doc pasted into the originating conversation (not saved as a file
in this repo — ask the user if you need its full text; this document
summarizes what it asked for and what was actually done in response).

## What the prior handover asked for

1. Run `iom_selftest.R` — this would be the module's first real execution
   (it had been authored by an instance of Claude without R access, so
   nothing in `iom_module.R` had ever actually been run before this).
2. Confirm Part B's CO2/CH4 comparison at 350C/500C matches the documented
   expectation (350C should match Miller closely since Ea was fit to it;
   500C legitimately overshoots Miller's raw yield for a documented reason
   — do not "fix" that without checking first).
3. Get Miller et al. (2025) Table 5 Murchison residue H/C, O/C at 350C/500C
   to fill in placeholder `NA`s and get a real residue comparison.
4. Only after that: recalibrate IOM_N/IOM_S the same way IOM_CO2/IOM_CH4
   had already been fixed (absolute measured value in, Ea out — not a
   ratio fit), checked against a second temperature.

## What was actually done, in order

### 1. First execution found a real bug — in the test script, not the module

Running `Rscript iom_selftest.R` failed on Part A's last check with
`Error in if (ok) "PASS" else "FAIL" : missing value where TRUE/FALSE
needed`. Root cause: `iom_selftest.R`'s "KINETICS block m0 scales linearly"
check used `regmatches(kb, gregexpr("-m0 [0-9.eE+-]+", kb))`, which matches
the literal `-m0 ` prefix *and* the number, so `as.numeric()` on the
captured string (`"-m0 4.169000e+00"`) returned `NA`.

**Fixed** by stripping the prefix before coercion (`sub("-m0 ", "", ...)`
in `iom_selftest.R`). Verified by hand in a scratch script that
`iom_kinetics_block()` itself was already correct — `-m0` does scale
exactly linearly with `iom_mass_kg` (4.169 → 8.338 at 2×). **The module was
not broken; the test harness was.** Worth independently confirming this
diagnosis rather than taking it on faith — re-run
`iom_kinetics_block(iom_default_config(), iom_mass_kg = X)` for a couple of
X values and check the `-m0` values by eye if you want a second check.

### 2. Part B (gas yields) matched the documented expectation exactly

350C: CO2 ratio 1.00, CH4 ratio 1.00 (fit points). 500C: CO2 ratio 1.22,
CH4 ratio 1.22 (expected overshoot, per the documented m0-basis
explanation: `m0_per_kg` includes H/O beyond the carbon Miller measured, so
full exhaustion of the pool exceeds his raw carbon-only yield). This
matched the prior handover's prediction, so no changes were needed here.

### 3. Residue trajectory — real Table 5 data filled in (partially)

The user supplied Miller Table 5 as a pasted table. The Murchison-specific
rows are `350-3-M` (3 kbar) and `500-3-M` (3 kbar). Read directly:
- 350C: H/C = 0.63 ± 0.08, O/C = 0.096 ± 0.011
- 500C: **H/C and O/C are blank in Table 5 itself** — that Murchison row
  only reports N/C (0.015) and C/N wt%, not H/C or O/C. This is a genuine
  gap in Miller's own data, not a transcription omission — do not go
  looking for a 500C Murchison H/C/O/C value elsewhere in the paper without
  first checking whether Miller reports one at all.

`miller_residue` in `iom_selftest.R` now has the real 350C values; 500C
stays `NA` with an inline comment explaining why.

Result: predicted residue at 350C is H/C=0.703 vs measured 0.63 (+11%),
O/C=0.104 vs measured 0.096 (+8%). Both a bit high, same direction. This is
the **first genuinely independent validation point** (residue was not part
of the two-point Ea fit), and the small systematic overshoot is consistent
with the model's known gap: a single Ea per channel can't reproduce the
low-Ea CO2/CH4 release Miller's own data implies at lower temperatures (see
`iom_module.R`'s caveats and the prior handover's item 3). **This has not
been investigated further or "fixed"** — it's flagged, not resolved.

### 4. IOM_N recalibrated against absolute Table 9 values; IOM_S left alone

User supplied Miller Table 1 (starting composition) and Table 9 (NH4+
yields) as pasted tables.

- Table 1: Murchison N = 2.48 ± 0.14 wt% → 1.7705 mol N / kg IOM total.
- Table 9, Murchison 3 kbar rows: `350-3-M` gives both a direct
  µmol N/mg sample figure (0.57, which is numerically `== mol/kg`) and "%
  of IOM nitrogen in NH4+" (32.06%); `500-3-M` gives 0.72 µmol/mg and
  40.83%. **These two independent readings of the same underlying release
  agree to within ~0.5%** when cross-checked via Table 1's N content
  (0.5676 mol/kg vs 0.57 direct; 0.7229 vs 0.72 direct) — this cross-check
  is itself worth re-verifying if anything about it looks off, since it's
  the basis for trusting the transcription.

Refit `Ea_N` the same way IOM_CO2/IOM_CH4 were fixed: solve
`released(350C) = m0_per_kg * (1 - exp(-k*t))` for Ea given the *absolute*
350C value (0.5676 mol/kg) and the existing `m0_per_kg = 0.875` (not
refit — only Ea was refit here). Result: **215 kJ/mol → 217.3 kJ/mol**.
500C prediction: ratio 1.21 vs Miller's measured value — the same
overshoot magnitude as CO2 (1.22) and CH4 (1.22), for the same m0-basis
reason. This cross-channel consistency is a reason for confidence, not
proof of correctness — if you want to sanity-check it independently,
redo the two-point fit by hand (Python or R, doesn't matter) rather than
trusting the number in the code comment.

**IOM_S was deliberately left untouched.** Miller does not measure H2S
yield anywhere in the paper (confirmed by the prior handover, not
re-derived here). `Ea_S` remains pinned to `Ea_N`'s value (now 217.3
kJ/mol, was 215) purely as a placeholder with zero independent support.
Do not treat `Ea_S` as calibrated just because a number sits in the config.

### 5. Part B's diagnostic printout converted to asserted pass/fail checks

Previously Part B only printed numbers for a human to eyeball. Per the
user's explicit request, added `pass()`-based assertions:

- **350C** (fit points): predicted must match Miller's absolute measured
  value within **1%**, for CO2, CH4, and N. (Tolerance sized to cover the
  ~0.3–0.4% numerical error of the 2000-step explicit-Euler integration in
  `run_isothermal()` — verified by hand in a scratch script that the raw
  ratios are 0.9981/1.0037/1.0033, not exactly 1.0000, before picking 1%
  rather than something tighter.)
- **500C**: does **not** assert against Miller's raw yield (that offset is
  deliberate — see above). Instead asserts the actual invariant the
  Ea/logA combination is supposed to produce: each channel reaches
  **≥99.9% of its own `m0_per_kg`** within the 48h experiment window. This
  was checked to be a real physical result, not a numerical-instability
  artifact: at 500C, `k*dt` for the Euler step is ≈2 (i.e. the true rate
  constant is so large that the ODE's analytic solution also predicts
  ~100% depletion well within one timestep) — confirmed by hand before
  relying on it. If you change `n_steps` in `run_isothermal()` down from
  2000, re-verify this claim; a coarser timestep could change the clamping
  behavior at 500C in ways that aren't obviously wrong but should be
  checked.
- **Residue** (350C only, since 500C measured values don't exist): a
  **loose 25% sanity bound**, not a tight target — deliberately loose,
  since the model is known to be missing the low-Ea tail. This catches
  gross errors (sign flips, unit mixups) without failing on the known,
  already-flagged ~10% overshoot.

All 16 checks (10 Part A + 6 Part B assertions, including 2 residue sanity
checks) currently pass; `Rscript iom_selftest.R` exits 0.

## Files changed on this branch (vs. `main`, which has neither file)

- `iom_module.R` — new file. Only change from the version described in the
  prior handover: `Ea_J` for `IOM_N` and `IOM_S` updated from `215e3` to
  `217.3e3`, with an inline dated comment explaining the refit (see
  `iom_default_config()`'s doc comment, "UPDATED 2026-09-17" block).
- `iom_selftest.R` — new file. Changes from the version described in the
  prior handover:
  - fixed the `-m0` regex/coercion bug (section 1 above)
  - `miller_residue` filled in with real Murchison 350C values (section 3)
  - `miller_measured`'s `NH3_pct_of_N` given precise values (32.06/40.83
    instead of rounded 32/41) and `murchison_N0` computed from Table 1
    instead of left as an approximate placeholder
  - added an `N_pred`/`N_meas`/`ratio` column to the printed comparison
    table
  - converted Part B from diagnostic-only to asserted checks (section 5)
- `INTEGRATION_NOTE.md` — carried over unchanged from the prior handover;
  not touched this session. **It has not been updated to reflect the
  IOM_N recalibration or the new asserted checks** — if this branch merges,
  that file's "what's NOT done" list should be revisited (IOM_N item is now
  partially resolved).

Note: `Primordial_CI_chondrite_mineralogy.Rmd` and its `.nb.html` also show
as modified/untracked in this repo — that is unrelated pre-existing work by
the user, not touched by this session, and was deliberately excluded from
the commit on this branch.

## Git state

Created branch `iom-module` off `main` (which had never had these files
committed at all — no history to worry about). One commit so far:
`5d5797b — Add standalone IOM decomposition module, split from ENTICES`,
containing `iom_module.R`, `iom_selftest.R`, `INTEGRATION_NOTE.md` as they
stood *before* the benchmark-hardening described in section 5. **The Part
B assertion changes (section 5) are, as of this document, uncommitted** —
check `git status` / `git diff` on this branch before assuming they're
saved anywhere.

## What is still NOT done (unchanged from or newly identified beyond the prior handover)

1. Low-Ea sub-pools for cold (~250C) conditions — not started.
2. Wiring into ENTICES's actual template generator — not started. See
   "INTEGRATION POINTS" comment block at the bottom of `iom_module.R`.
3. The 350C residue overshoot (+11% H/C, +8% O/C) is flagged, not
   explained beyond "probably the missing low-Ea tail" — that's an
   inference, not a verified diagnosis. Don't repeat it as settled fact.
4. 500C residue H/C, O/C remain permanently unverifiable against Murchison
   specifically unless a different Miller table or a different source has
   the number — check before assuming it's obtainable.
5. `IOM_S` is exactly as uncalibrated as before this session — nothing
   changed there except its Ea number moved because it's pinned to
   `Ea_N`, not because anyone calibrated H2S.
6. This branch has not been merged to `main`, and the user has not yet
   confirmed they want the Part B assertion changes committed. Don't
   assume permission to merge or push anything from this document alone —
   confirm with the user first.

## Suggested first step for whoever picks this up

Run `Rscript iom_selftest.R` yourself (R is at
`C:\Program Files\R\R-4.5.2\bin\Rscript.exe` on this machine, not on
PATH — the plain `Rscript` command will fail with "not recognized").
Confirm the same 16/16 pass result described here before trusting anything
else in this document. Then decide, with the user, whether to commit the
Part B changes and/or merge this branch.
