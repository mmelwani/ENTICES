# PATCH F — Vitrimat 2018 parameters and distributed Ea

Date: 2026-09-18
Source: Burnham (2019), Org. Geochem. 131, 50-59, **Table 1 "Parameters for
Vitrimat 2018 (vitrinite)"**, read directly from the project PDF (a normal
PDF; `pdftotext -layout` works).
Supersedes: the Ea values in `iom_default_config()` and the BS89-shape fitting
approach in `PROGRESS.md` Part 4 / `HANDOVER_distributed_Ea_and_option3.md`
Task 2.

This is a specification, not a file. The `iom_module.R` copy available to this
session predates the VS Code session's IOM_N recalibration (it still has
`Ea_J = 215e3`), so shipping a full replacement file would silently revert
that work. Apply these changes to the branch version.

---

## 1. Confirmed: A = 2e15 s^-1, and our assumed 1e13 was the OLD value

Burnham (2019) p.52, and Table 1's `A, s^-1` row: **A = 2 x 10^15 s^-1** for
all four product channels in Vitrimat 2018.

The paper also states explicitly why we had 1e13 (p.51-52): *"both Vitrimat
and Easy%Ro used A = 1 x 10^13 s^-1, which is approximately equal to the
transition..."* — that is the **1989** value. The abstract describes the 2018
revision as using "substantially higher frequency factors," and notes the
larger A improves "agreement with hydrous and confined pyrolysis without
changing the results for geological maturation." Miller's experiments *are*
hydrous pyrolysis, so 2e15 is the right choice for our application.

`logA = log10(2e15) = 15.3010`.

### Effect on our existing single-Ea fits

A and Ea trade off, so refitting the same 350 C data at the new A shifts every
Ea by **+27.4 kJ/mol** (not the "+21" I estimated in an earlier comment — that
figure was a hand-written guess left in a print statement, the computation is
+27.4):

| channel | Ea at A=1e13 | Ea at A=2e15 |
|---|---|---|
| IOM_CO2 | 216.2 | 243.6 |
| IOM_CH4 | 229.4 | 256.9 |
| IOM_N | 217.3 | 242.8 [CORRECTED — was 244.8, basis bug; see PROGRESS.md Part 4] |

**Nothing we previously validated breaks** — the 350 C fits are preserved by
construction. But low-temperature extrapolation changes by ~1-2 orders of
magnitude, in the *opposite* direction from the distribution effect (higher A
with matched 350 C behaviour means *less* low-T release for a single Ea).

**If you keep single-Ea for any channel, its Ea must be updated.** Leaving
Ea at the A=1e13 value while setting logA=15.301 would be wrong twice over.

---

## 2. The Vitrimat 2018 distributions (Burnham 2019 Table 1)

All four channels, percent of precursors by Ea in kcal/mol. **This includes
the CH4 distribution** that was item 5 on the "need to source" list — it is
here, and it supersedes BS89 for CH4 as well as CO2:

| Ea kcal/mol | H2O | CO2 | CHn | CH4 |
|---|---|---|---|---|
| 42 | 10 | | | |
| 44 | 10 | 10 | | |
| 46 | 15 | 15 | | |
| 48 | 15 | 15 | | |
| 50 | 15 | 15 | 5 | |
| 52 | 15 | 15 | 15 | 2 |
| 54 | 10 | 15 | 30 | 5 |
| 56 | 10 | 15 | 20 | 8 |
| 58 | | | 15 | 10 |
| 60 | | | 10 | 12 |
| 62 | | | 5 | 15 |
| 64 | | | | 12 |
| 66 | | | | 10 |
| 68 | | | | 8 |
| 70 | | | | 6 |
| 72 | | | | 5 |
| 74 | | | | 4 |
| 76 | | | | 3 |

All four columns sum to 100%. Other Table 1 parameters, for reference:
initial H/C 0.90, final H/C 0.20, initial O/C 0.35, a(water) 21% of O,
b(CO2) 74% of O, c(oil) 2% of C, n(oil H/C) 1.8.

Note this is **vitrinite**, not IOM. That matters for interpreting section 3.

---

## 3. THE KEY RESULT: Vitrimat 2018's CO2 distribution independently
reproduces Miller's 250 C data to 1%

This is the strongest validation the organic module has had, and it was not
fitted to anything of ours.

Vitrimat 2018 CO2 fractional conversion in 48 h, as published, no adjustment:

    250 C: 0.4581     350 C: 0.9999     500 C: 1.0000

Miller's HC113 measured CO2 (Table 2, mol/kg): 250 C mean of four runs
= 1.158; 500 C = 2.55. Measured ratio **0.4539**.

    Vitrimat 2018 predicted 250/500 ratio : 0.4581
    Miller HC113 measured                 : 0.4539
    agreement                             : 0.991

**A distribution Burnham fitted to vitrinite coal predicts Miller's IOM CO2
release split to within 1%.** For comparison, our single-Ea model gets this
ratio wrong by a factor of **1016** (at A=1e13) or **2766** (at A=2e15).

This is why distributed Ea is not a refinement. It is the difference between a
model that reproduces an independent measurement and one that is wrong by
three orders of magnitude.

---

## 4. A real tension — report it, do not average it away

Fitting the Vitrimat 2018 CO2 *shape* with a uniform shift + m0:

| sample | fit | result |
|---|---|---|
| Murchison (350+500 C) | shift **+3.6** kcal/mol, m0 3.431 mol/kg | 350 C 3.038 (meas 3.04); 500 C 3.431 (meas 3.43); RMS 0.002 |
| HC113 (250+500 C) | shift **0.0** kcal/mol, m0 2.528 mol/kg | 250 C 1.166 (meas 1.158); 500 C 2.546 (meas 2.55) |

Each sample fits its own data almost perfectly. **They disagree on the shift.**

Why, mechanically: at shift 0 both 350 C and 500 C are essentially fully
converted (0.9999 / 1.0000), so the model *cannot* reproduce Murchison's
observed 350/500 difference (3.04 vs 3.43, ratio 0.886) — it predicts 1.000.
Forcing that difference requires pushing Ea up until 350 C is only partly
converted, which then destroys the 250 C agreement.

**This is most likely real material difference, not model error.** HC113 is
synthetic (H/C 1.13, O/C 0.27, Table 1); Murchison is meteoritic (H/C 0.70,
O/C 0.18). Miller's own Table 10 independently needed HC113 shifted **+8
kcal/mol** relative to BS89 — so the two materials genuinely differ in Ea, and
in the same direction we find here (syn-IOM more refractory than a
vitrinite-derived baseline... though note Miller's +8 for HC113 vs our 0.0 for
HC113 is itself unreconciled; his is a residue fit at 10 C/Myr, ours a gas-yield
fit at 48 h isothermal).

**Do not fit one shift to both samples.** Carry them as separate
parameterizations and report the spread.

### What the spread costs at Enceladus temperatures

CO2 fraction released over 4.5 Gyr:

| T | Murchison (+3.6) | HC113 (0.0) | spread |
|---|---|---|---|
| 0 C | 0.0000 | 0.0002 | 759x |
| 25 C | 0.0004 | 0.0876 | 229x |
| 50 C | 0.0960 | 0.3501 | 3.6x |
| 75 C | 0.3623 | 0.6269 | 1.7x |
| 100 C | 0.6389 | 0.9042 | 1.4x |

Below ~50 C the material choice dominates everything. Above ~75 C it barely
matters. **Use Murchison for production** (our bulk composition is
Murchison-like) and run HC113 as the upper bound.

### And what the single-Ea model was getting wrong

Against the Murchison distribution fit, over 4.5 Gyr:

    0 C   : single-Ea low by ~3e7 x
    25 C  : low by ~6e6 x
    50 C  : low by ~8e5 x
    100 C : low by 29 x

---

## 5. Concrete config change

Replace the single `IOM_CO2` and `IOM_CH4` rows with sub-pools. Set
`logA = 15.301` on **every** row (including IOM_N and IOM_S). Sub-pool `m0`
values below are on the **ENTICES bulk basis** (Murchison fit x 1.215, where
1.215 = 60.611/49.87 = ENTICES bulk C per kg / Murchison C per kg).

`Ea_J` values are `(kcal + 3.6) * 4184` for the Murchison fit.

    IOM_CO2 sub-pools (total m0 = 4.1691 mol/kg IOM):
      name            m0_per_kg   formula             Ea_J      logA
      IOM_CO2_44      0.416910    C 1 H 1 O 1.3115    199160    15.301
      IOM_CO2_46      0.625365    C 1 H 1 O 1.3115    207528    15.301
      IOM_CO2_48      0.625365    C 1 H 1 O 1.3115    215896    15.301
      IOM_CO2_50      0.625365    C 1 H 1 O 1.3115    224264    15.301
      IOM_CO2_52      0.625365    C 1 H 1 O 1.3115    232632    15.301
      IOM_CO2_54      0.625365    C 1 H 1 O 1.3115    241000    15.301
      IOM_CO2_56      0.625365    C 1 H 1 O 1.3115    249368    15.301

    IOM_CH4 sub-pools (total m0 = 2.2718 mol/kg IOM; shift +3.7):
      IOM_CH4_52      0.045436    C 1 H 4.435         233472    15.301
      IOM_CH4_54      0.113590    C 1 H 4.435         241840    15.301
      IOM_CH4_56      0.181744    C 1 H 4.435         250208    15.301
      IOM_CH4_58      0.227180    C 1 H 4.435         258576    15.301
      IOM_CH4_60      0.272616    C 1 H 4.435         266944    15.301
      IOM_CH4_62      0.340770    C 1 H 4.435         275312    15.301
      IOM_CH4_64      0.272616    C 1 H 4.435         283680    15.301
      IOM_CH4_66      0.227180    C 1 H 4.435         292048    15.301
      IOM_CH4_68      0.181744    C 1 H 4.435         300416    15.301
      IOM_CH4_70      0.136308    C 1 H 4.435         308784    15.301
      IOM_CH4_72      0.113590    C 1 H 4.435         317152    15.301
      IOM_CH4_74      0.090872    C 1 H 4.435         325520    15.301
      IOM_CH4_76      0.068154    C 1 H 4.435         333888    15.301

**RE-DERIVE THESE.** They were computed in Python in a session with no R. The
recipe: m0_total = (least-squares fit of Vitrimat-2018-shape conversion
fractions to Miller's two Murchison yields) x 1.215; per-sub-pool m0 =
m0_total x (Table 1 percent / 100); Ea_J = (Ea_kcal + shift) x 4184. If your
numbers differ, yours are probably right — say so.

Note the CH4 m0 (2.27) now exceeds the old single-Ea value (1.908) because
fitting both temperature points with a distribution that is only 27% converted
at 350 C requires a larger pool. Check this against
`iom_validate_config()`'s element-balance assertion — combined C release is now
4.169 + 2.272 = 6.44 mol/kg, up from 6.08. **That is still well under bulk C
(60.611), so the residue stays positive, but re-run the assertion.**

IOM_N and IOM_S stay single-pool (no distribution data for either), but **both
need `logA = 15.301` and Ea refit to 242.8 [CORRECTED — was 244.8, basis bug; see PROGRESS.md Part 4] kJ/mol** (IOM_N; IOM_S remains
pinned to IOM_N and uncalibrated).

---

## 6. Validation now available that was not before

1. **250 C CO2 ratio** (HC113): independent, passes at 0.991 with the
   unshifted distribution. Add as a self-test assertion.
2. **Residue trajectories**: Table 5 is now fully in hand, including syn-IOM
   at 250 C and HC113 at 350 C (the run I previously reported as missing does
   exist: `350-10-Mar2023`, 10 kbar, H/C 0.66, RH/C 0.58, O/C 0.047,
   RO/C 0.18). Full table is in the project as `Miller 2025 - Table 5.txt`.
   The `RX/C` columns are the ones to compare against — dimensionless.
3. **Murchison 350 C residue** remains the main independent check
   (H/C 0.63 +/- 0.08, O/C 0.096 +/- 0.011). Worth re-running now that the
   low-Ea tail exists: the previous +11%/+8% overshoot was attributed to the
   missing tail, and this patch is the direct test of that hypothesis.
4. **Murchison 500 C** H/C and O/C are confirmed genuinely blank in Table 5
   (that row has only N/C 0.015, RN/C 0.42, C 55.9 wt%, N 0.99 wt%). Do not
   look for them again.

---

## 7. Still open

- **Miller's +8 kcal/mol for HC113 (Table 10) vs our 0.0 for HC113.** Both are
  HC113, both derived from Miller's own data, and they disagree. His is a
  residue-trajectory fit at 10 C/Myr geologic heating; ours is a gas-yield fit
  at 48 h isothermal. Worth understanding before the paper, but it does not
  block production: our Murchison fit does not depend on it.
- **Nitrogen composition discrepancy**, now three-way and still unresolved:
  our bulk formula uses N 2.142 per formula unit; an old comment said N 3.284;
  Miller's Table 1 gives Murchison N = 2.48 +/- 0.04 wt%, which is 1.77 mol/kg
  — i.e. N/C 0.036, so on ENTICES's C 60.611 basis, **N 2.18**. That is close
  to 2.142 and nowhere near 3.284, so 2.142 is probably right and 3.284 was
  the error — but confirm against whatever source the ENTICES bulk composition
  came from (possibly `Primordial_CI_chondrite_mineralogy.Rmd`).
- **H2O and CHn channels** are now available from Vitrimat 2018 and are not in
  our config at all. H2O release is chemically uninteresting for us (it goes
  into the solvent) but CHn is real organic carbon we are currently ignoring.
  Worth deciding deliberately rather than by omission.