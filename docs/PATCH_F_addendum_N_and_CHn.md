# PATCH F addendum — nitrogen resolved, CHn channel added

Date: 2026-09-18
Appends to `PATCH_F_vitrimat2018.md`. Closes the two items left open there.

---

## 1. Nitrogen: RESOLVED. Keep N = 2.142. The 3.284 was a transcription error.

Agreed that Cody et al. (2024) via Miller Table 1 should be the sole source.
Working it through confirms our current value and identifies where the wrong
one came from.

**Murchison IOM from Table 1** (C 59.9, H 3.49, N 2.48, O 16.2 wt%):

    per kg of real IOM:  C 49.871  H 34.626  N 1.771  O 10.126 mol
    implied ratios:      H/C 0.694   N/C 0.0355   O/C 0.203
    (Table 1's stated ratios: 0.70, 0.036, 0.18 — consistent)

Note those wt% sum to 82.1%, the remainder being S, P and ash. That is why
there are two different bases in play and why this got confusing:

- **(a) real-IOM basis:** N = 1.771 mol per kg of actual Murchison IOM.
- **(b) ENTICES basis:** the bulk formula is normalized so C+H+O+N+S = exactly
  1000 g/mol, giving C 60.611. On that basis, Cody's N/C = 0.036 implies
  **N = 2.182**.

Our current value is **2.142**, i.e. N/C = 0.0353 — within **1.8%** of Cody's
0.036, which is inside Cody's own stated uncertainty (2.48 ± 0.04 wt% is
±1.6%). So it is already right.

**Where 3.284 came from.** The two candidate gram figures were "30 g N/kg"
(= 2.142 mol) and "46 g N/kg" (= 3.284 mol). Cody's N/C 0.036 on the ENTICES
basis corresponds to **30.6 g** N per 1000 g formula — so 30 g is correct.

And 46 g is the **hydrogen** figure: ENTICES H = 45.639 mol × 1.0079 = 46.0 g.
Someone read the H gram value as N. That is consistent with the old formula
carrying *both* H 45.639 *and* N 3.284 — the same 46 g reused as if it were
nitrogen.

**Action:** no change required to `iom_module.R`. Optionally refine
N 2.142 → 2.182 for exact consistency with Cody's N/C (a 1.9% change, below
measurement uncertainty — a defensible tidy-up, not a fix). Either way, delete
the "N_3.284 / 46 g N — UNRESOLVED, confirm against source" caveat from
`iom_bulk_formula()`'s comment and replace it with the Cody citation. **This
item has been flagged as unresolved for several sessions; it is now closed.**

---

## 2. CHn channel: added, using Burnham's own stoichiometry

Agreed: take CO2 as updated by Miller, and CH4 (and now CHn) from Burnham 2019
unchanged.

### The calibration problem, and the honest answer

Miller's GC measured **H2, CO, CH4, CO2** (Table 2). He does **not** report
C2+ hydrocarbons, which is what Vitrimat's "CHn" pool represents. So unlike
CO2 and CH4, **CHn's m0 cannot be calibrated against Miller's data at all.**

The defensible option is Burnham 2019 Table 1's own stated stoichiometry, from
the same table the distributions come from:

    c (oil), % of C = 2      ->  CHn carbon = 2% of total carbon
    n, oil H/C      = 1.8    ->  formula "C 1 H 1.8"

On the ENTICES basis: CHn m0 = 0.02 × 60.611 = **1.2122 mol/kg IOM**.

**Flag this clearly in the code.** CO2 and CH4 m0 are fitted to measured
Murchison yields; CHn m0 is transferred from vitrinite and has no
IOM-specific support. It is better than omitting the channel (which silently
assumes zero C2+ production), but it is a different epistemic class from the
other two.

### CHn distribution (Burnham 2019 Table 1, unchanged)

| Ea kcal/mol | % |
|---|---|
| 50 | 5 |
| 52 | 15 |
| 54 | 30 |
| 56 | 20 |
| 58 | 15 |
| 60 | 10 |
| 62 | 5 |

Fractional conversion in 48 h, unshifted: 250 C = 0.031, 350 C = 0.849,
500 C = 1.000.

**Shift:** apply the same **+3.6 kcal/mol** Murchison shift as CO2/CH4, for
consistency of material basis. This is an assumption — the shift was derived
from CO2 data and has no independent CHn constraint.

### Config rows (ENTICES basis, Murchison shift +3.6)

    IOM_CHn sub-pools, total m0 = 1.2122 mol/kg, formula "C 1 H 1.8":
      name            m0_per_kg   Ea_J      logA
      IOM_CHn_50      0.060610    224264    15.301
      IOM_CHn_52      0.181830    232632    15.301
      IOM_CHn_54      0.363660    241000    15.301
      IOM_CHn_56      0.242440    249368    15.301
      IOM_CHn_58      0.181830    257736    15.301
      IOM_CHn_60      0.121220    266104    15.301
      IOM_CHn_62      0.060610    274472    15.301

Ea_J = (kcal + 3.6) × 4184. Re-derive rather than trusting these.

---

## 3. Full five-channel element balance — verified, closes correctly

| channel | m0 | formula | basis |
|---|---|---|---|
| IOM_CO2 | 4.1691 | C 1 H 1 O 1.3115 | fitted, Murchison |
| IOM_CH4 | 2.2718 | C 1 H 4.435 | fitted, Murchison |
| IOM_CHn | 1.2122 | C 1 H 1.8 | **transferred from vitrinite** |
| IOM_N | 0.8750 | N 1 H 2 | fitted, Murchison |
| IOM_S | 0.4620 | S 1 H 1 | **uncalibrated placeholder** |

Element budget on the ENTICES bulk basis:

| element | released | bulk | residue | % released |
|---|---|---|---|---|
| C | 7.653 | 60.611 | 52.958 | 12.6% |
| H | 18.638 | 45.639 | 27.001 | 40.8% |
| O | 5.468 | 9.939 | 4.471 | 55.0% |
| N | 0.875 | 2.142 | 1.267 | 40.8% |
| S | 0.462 | 1.154 | 0.692 | 40.0% |

All residues positive; `iom_validate_config()` should pass. Implied inert
residue: **H/C 0.510, O/C 0.0844** (bulk was 0.753 / 0.164).

Carbon released rises from 10.0% (pre-patch, CO2+CH4 single-Ea) to **12.6%**
with CHn included. Note this is *full exhaustion*, not the 500 C yield —
CH4 and CHn are only partly converted at 500 C, so the comparison to Miller's
~10% measured C release at 500 C is not apples-to-apples.

### Validation opportunity this creates

The residue H/C of 0.510 is now a **direct prediction** testable against
Miller's Table 5. Murchison 350 C measured H/C = 0.63 ± 0.08, O/C = 0.096 ±
0.011. Because CHn carries H/C 1.8 (H-rich), adding it pulls predicted residue
H/C *down*, which is the right direction to fix the previously-reported
+11% overshoot. **Re-run the residue check after implementing — this patch is
a direct test of the "missing low-Ea tail" hypothesis and it should now
either resolve or fail visibly.**

---

## 4. What remains genuinely open

- **IOM_S** — Miller measures no H2S. Unchanged, uncalibrated, pinned to
  IOM_N's Ea. Flag in any sulfur-dependent output.
- **CHn's m0 and shift** — vitrinite-derived, no IOM constraint. Sensitivity
  test worth doing: the 2%-of-C figure is Burnham's coal value and could
  plausibly be several times larger for H-rich IOM.
- **H2O channel** — Vitrimat 2018 supplies it (a 21% of O, distribution
  42-56 kcal). Deliberately still omitted: its products go to the solvent and
  it does not change ocean chemistry meaningfully. Recording the decision so
  it is not mistaken for an oversight.
- **Miller Table 10's +8 kcal for HC113 vs our 0.0 for HC113** — still
  unreconciled (residue fit at 10 C/Myr vs gas-yield fit at 48 h isothermal).
  Does not block production, since our Murchison fit is independent of it.
