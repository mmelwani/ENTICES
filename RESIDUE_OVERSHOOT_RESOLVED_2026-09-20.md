# The residual 350 C residue overshoot is a starting-composition artefact, not a model error

Date: 2026-09-20
Closes: `PROGRESS.md` Part 5 item 6 ("chase the residual +7.8% residue
overshoot, or accept it"); revises the diagnosis in Part 4's validation state.
Method: R only, no new data. Reproduce with `residue_ratio.R` (below).

---

## Result

The +7.7% overshoot in predicted 350 C residue H/C is **inherited almost
entirely from the model's starting composition**, not produced by the
decomposition model. Judged on the dimensionless `RX/C` ratios — which
`PROGRESS.md` and the Miller-table handovers already say are the correct
cross-sample comparison — the distributed-Ea model reproduces Miller's
Murchison 350 C residue to **+0.1%**.

| | model | Miller (350-3-M) | error |
|---|---|---|---|
| **starting** H/C | 0.7530 | 0.70 | +7.6% |
| **final** H/C (350 C) | 0.6785 | 0.63 | **+7.7%** |
| **R(H/C)** = final/initial | 0.9011 | 0.90 ± 0.12 | **+0.1%** |
| **R(O/C)** = final/initial | 0.5544 | 0.53 ± 0.06 | +4.6% (inside stated uncertainty) |

The decisive line:

    predicted/measured FINAL    H/C = 1.0770
    model/measured   STARTING   H/C = 1.0757
    difference contributed by the model = +0.13%

So the decomposition treatment adds **0.13%** of discrepancy. The other 7.6%
was already present before any IOM reacted.

O/C tells the same story in the other direction: the model starts 8.9% too
O-**poor** (0.164 vs 0.180) and ends 5.3% too O-poor.

## What this means for the missing-low-Ea-tail hypothesis

`PROGRESS.md` Part 4 currently records the residue re-check as "partially
supported, not confirmed", with three candidate causes for the residual: m0
basis, unreacted channel mass counted as solid, sample heterogeneity.

**It is the first of those, and it is now settled.** The residual is a
composition-basis offset between `iom_bulk_formula()` and the Murchison
material Miller measured. It is not evidence for or against a missing low-Ea
tail, because it is present at zero reaction extent.

The honest statement is therefore stronger than "partially resolved": on the
correct (dimensionless) comparison, the 350 C residue check **passes**, and it
remains the one genuinely independent validation point since residue was never
used in any fit.

Two things this does **not** license:
- It does not retroactively confirm the low-Ea tail. The tail is justified by
  the 250 C HC113 result (0.991) and by the distributed-vs-single-Ea gap at
  Enceladus temperatures, not by this.
- It does not mean the absolute numbers are fine to quote. Any absolute
  residue H/C or O/C comparison against Murchison will carry this ~8% offset.
  Quote `RX/C`, or reconcile the composition (below).

## Where the starting composition comes from, and the choice it implies

`iom_bulk_formula()` returns `C 60.611, H 45.639, O 9.939, N 2.142, S 1.154`.
Traced to `Primordial_CI_chondrite_mineralogy.Rmd` (~line 1316), which states:

    # Target formula: C60.611H45.639O9.941N3.284S1.154

and decomposes it into KerogenC128 (C128H68O7) + pyridine (C5H5N) + S + O —
the same Kerogen/Pyridine phases that appear in ENTICES's
`EQUILIBRIUM_PHASES` list. So this is a **CI-chondrite-derived** composition,
while the decomposition kinetics are calibrated against **Murchison**, a CM2.
CI and CM IOM genuinely differ in H/C and O/C, so the offset may be entirely
deliberate.

**That is a scientific choice for the PI, not a bug to fix, and it needs an
explicit decision:**

- **Keep CI composition + Murchison kinetics.** Defensible if the core is
  modelled as CI. Then all residue validation must use `RX/C`, and that should
  be stated wherever residue agreement is claimed.
- **Switch the bulk formula to Murchison** (H/C 0.70, O/C 0.18) for internal
  consistency with the calibration. This would change the element budget,
  the channel m0 values on the ENTICES basis, and the absolute residue
  predictions.

Either is defensible. What is not defensible is quoting absolute residue
agreement while carrying a CI composition against CM data.

## A live inconsistency found while tracing this

The source Rmd's target formula still carries **N3.284**, the value
`PATCH_F_addendum_N_and_CHn.md` identified as a transcription error (the
hydrogen gram figure, 46 g, read as nitrogen) and closed in favour of
**N 2.142**. `iom_module.R` uses 2.142; `Primordial_CI_chondrite_mineralogy.Rmd`
still uses 3.284.

So the addendum's own open question — "confirm against whatever source the
ENTICES bulk composition came from (possibly
`Primordial_CI_chondrite_mineralogy.Rmd`)" — is now answered: **that source
does use 3.284**, and it disagrees with the module. One of the two needs
updating. Note the Rmd solves for its Kerogen/Pyridine/S/O coefficients
*from* that target, so changing N there changes the pyridine coefficient and
hence the mineral endmember breakdown — i.e. this is not a one-line edit, and
it is the user's call whether it matters. Flagged, not touched.

---

## Reproduce

`residue_ratio.R` — sources `iom_module.R`, integrates the channels
analytically (closed form, not stepped: for an isothermal hold `k` is constant
so `m(t) = m0*exp(-k*t)` is exact), rebuilds the residue as inert + unreacted
channel mass, and differences against Miller Table 5 row `350-3-M`
(H/C 0.63 ± 0.08, RH/C 0.90 ± 0.12, O/C 0.096 ± 0.011, RO/C 0.53 ± 0.06).
Nothing was tuned.
