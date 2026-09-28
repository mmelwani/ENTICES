# Handover: distributed-Ea implementation + Option 3 (for Claude Code / VS Code)

Date: 2026-09-18 (REVISED same day — see "SUPERSEDED" notice below)
Companion docs: `PATCH_F_vitrimat2018.md` (READ THIS FIRST — it supersedes
Task 2's fitting approach and every Ea value quoted here), `PROGRESS.md` (full rationale — read Parts 1, 2 and 4 before
starting), `iom_module.R`, `iom_selftest.R`, `HANDOVER_for_vscode_claude.md`
(earlier handover; still accurate on norms and environment).

Read `PROGRESS.md` Part 1 first. It explains why channel `-formula` values are
elemental ("what the solid loses") and not product molecules. That is the most
likely thing to get "helpfully" broken.

---

## SUPERSEDED — read `PATCH_F_vitrimat2018.md` before Task 2

The user supplied Burnham (2019) and Miller Tables 1 and 5 after this document
was written. Three things changed:

1. **A = 2e15 s^-1, not 1e13.** Confirmed from Burnham 2019 Table 1. Our 1e13
   was the *1989* Vitrimat value, explicitly superseded for hydrous pyrolysis.
   `logA = 15.301`. Every Ea in this document is therefore wrong by
   **+27.4 kJ/mol** and must be refit. The 350 C fits still hold (A/Ea trade
   off), so nothing previously validated breaks.
2. **Use the Vitrimat 2018 distributions, not BS89 shape-fitting.** Burnham
   2019 Table 1 gives full Ea distributions for H2O, CO2, CHn **and CH4**
   (item 5 on the "need to source" list below — now sourced). Its CO2
   distribution reproduces Miller's HC113 250 C/500 C split to **0.991**
   with no adjustment, versus a factor-of-1000 error for single-Ea. Task 2's
   "fit BS89 shape with a shift" recipe is replaced by the explicit sub-pool
   config in `PATCH_F_vitrimat2018.md` section 5.
3. **Task 3's tension is resolved but replaced by a narrower one.** It was
   partly the wrong A. What remains: Murchison wants shift +3.6 kcal/mol,
   HC113 wants 0.0, and each fits its own data to <1%. That is real material
   difference (syn-IOM vs meteoritic), corroborated by Miller's own Table 10
   needing +8 for HC113. **Do not fit one shift to both.** Use Murchison for
   production, HC113 as an upper bound; below 50 C the two differ by
   200-750x, which is the honest uncertainty on cold-case organic release.

---

## Why you're being asked to do this

Two findings changed the priority order since the last session.

**1. The single-Ea approximation is wrong by ~1000x at the temperatures that
matter.** Our Ea values were fitted to 48 h experiments at 350 and 500 C. At
those temperatures the low-Ea end of the real activation-energy distribution
is *fully converted*, so our calibration data contain almost no information
about it. But it is precisely that low-Ea tail that controls organic release
at Enceladus temperatures (0-100 C over Myr). Measured directly against
Miller's own 250 C data: the model predicts 0.045% conversion where the
experiment shows 45%. **A factor of ~1000.** The good 350/500 C agreement is
not evidence of correctness at low T — structurally it cannot be.

**2. The coupled-redox treatment imposes a large spurious electron demand.**
At Enceladus conditions PHREEQC reduces all released carbon to CH4, drawing
~24 mol electrons (~12 mol H2-equivalent) per kg IOM from dissolved H2 and
from oxidising reduced minerals. At low temperature CO2/CH4 (and
sulfate/sulfide, and ammonia/N2/nitrate) are kinetically inhibited and should
*not* equilibrate. See `PROGRESS.md` Parts 2 and 3.

Item 1 changes how much carbon exists to speciate, so it comes first —
otherwise Option 3 work gets re-run.

---

## TASK ORDER

### Task 0 — Environment and baseline (do first, ~10 min)

R is at `C:\Program Files\R\R-4.5.2\bin\Rscript.exe`, not on PATH.
Confirm `iom_selftest.R` still passes 16/16 before changing anything. Also
check `git status` on branch `iom-module` — per the last handover, the Part B
assertion changes were **uncommitted**. Establish a clean baseline and
confirm with the user before committing or merging anything.

### Task 1 — Redox budget diagnostic (cheap, may invalidate existing results)

Take an existing ENTICES run with `organic_wt_percent > 0` and the same case
with organics off. Compare, at matched timesteps:
- dissolved H2
- pe
- Fe-bearing mineral assemblage (magnetite, Fe(III) phases), and pyrrhotite

Do **not** just check H2. The electron demand may be satisfied by mineral
oxidation instead, in which case H2 looks fine and the mineralogy is wrong.
Report what you find; don't adjust anything yet.

### Task 2 — Distributed Ea for IOM_CO2 (the main implementation job)

Replace the single `IOM_CO2` channel with N sub-pools sharing one `-formula`
but each with its own `Ea_J` and a fraction of `m0_per_kg`.

Design constraints:
- Keep the existing interface (`iom_rates_block`, `iom_kinetics_block`, etc.)
  unchanged in signature. Sub-pools are a config-level change: more rows.
- Name them predictably (`IOM_CO2_1` ... `IOM_CO2_7`) so `iom_punch_headings`
  and `iom_punch_lines` keep working and the restart logic can recover them.
- `iom_validate_config()`'s element-balance check must still pass: the
  sub-pools' summed `m0_per_kg` must equal the parent channel's.
- Add a self-test assertion: sub-pool total release must equal the
  single-channel release in the high-temperature limit (where both exhaust).

Use the **BS89 distribution shape** as a prior and fit only a uniform shift
and `m0`. The shape (weights, 2 kcal/mol spacing) is:

    offset from lowest bin:  0    2    4    6    8   10   12   kcal/mol
    weight %:                5   15   25   25   15   10    5

Fits already computed (verify these independently — see "Verify, don't
trust" below):

| dataset | fit | 250 C | 350 C | 500 C |
|---|---|---|---|---|
| Murchison, 350+500 C | shift +0.6, m0 3.40 | — | 3.043 (meas 3.04) | 3.400 (meas 3.43) |
| HC113, 250+500 C | shift −2.9, m0 2.50 | 1.164 (meas 1.158) | — | 2.500 (meas 2.55) |

**Use the Murchison fit for production** (our bulk IOM composition is
Murchison-like). Scale m0 to the ENTICES basis by x1.215 (= 60.611/49.87,
ENTICES bulk C / Murchison C per kg): 3.40 -> ~4.13, close to the current
4.169.

### Task 3 — Reconcile a tension before trusting Task 2

The two fits above disagree in direction, and this needs explaining, not
averaging:
- Our HC113 gas-yield fit wants Ea **lower** than Murchison's (−2.9 vs +0.6).
- Miller's Table 10 has HC113 needing Ea **+8 kcal/mol higher** than BS89.

These are not directly comparable — Miller fitted *residue trajectories* at a
*geologic* heating rate (10 C/Myr), we fitted *gas yields* at *48 h
isothermal*. But the sign disagreement means at least one of the following is
true, and you should work out which: (a) the two observables genuinely imply
different Ea; (b) A ≠ 1e13 for Vitrimat, shifting everything; (c) our m0
basis for HC113 is wrong because 500 C is not actually its exhaustion point;
(d) pressure matters (HC113 250 C data are 1 and 3 kbar, the 500 C point is
10 kbar). **(c) and (d) are the most likely and the cheapest to check.**

Do not proceed to production runs until this is understood or explicitly
parked with the user's agreement.

### Task 4 — Option 3 (coupled/decoupled redox)

Full design is in `PROGRESS.md` Part 3. Four steps, in order, each verified
before the next:
1. Database work **in isolation**: add decoupled pseudo-species (Mtg-style
   methane, Amm-style ammonium, and sulfide) to a copy of
   `Core11_idealgas_mod_v2.dat`. Test with a hand-written PHREEQC input: put
   CO2 and CH4 in a strongly reducing solution, confirm they do **not**
   interconvert with the new database and **do** with the old one. No ENTICES,
   no module.
2. Module support: add a decoupled-mode formula column and a per-element
   `redox_mode` (decouple C, N, S independently — they have different
   inhibition temperatures).
3. Element-conservation assertion: total C, H, O, N, S released must be
   identical in both modes; only partitioning differs. **Write this before
   the database work.**
4. Wire into ENTICES and run dual cases.

Note the semantic inversion: in coupled mode `-formula` is what the solid
loses and PHREEQC picks products; in decoupled mode product identity becomes
an *input*, because nothing downstream re-partitions it. Draft decoupled
formulas are in `PROGRESS.md` Part 3 and are **unverified** — PHREEQC's
acceptance of `C(4)` notation in a KINETICS `-formula`, and the
charge/electron bookkeeping when mixing a decoupled pseudo-element with free
H, both need testing.

### Task 5 — Later
CH4 channel distributed the same way (BS89 puts CH4 precursors higher, to
~74 kcal/mol — but Table 10 only gives the CO2 distribution, see below).
Then wire the module into ENTICES's template generator
(`INTEGRATION POINTS` block at the bottom of `iom_module.R`).

---

## DATA APPENDIX — everything extracted so far

The Miller PDF at `Miller2025_light.pdf` has a **ZIP header, not a PDF one**
(page-by-page archive). Standard PDF tools fail; `unzip` works and yields
`1.txt` ... `19.txt`. Don't rediscover this.

Units note: Table 2 is in `mmoles/mg sample`, which equals **mol/kg** exactly.
Table 9 is in `µmoles/mg`, which also equals mol/kg. Convenient, and already
verified two independent ways (direct µmol/mg vs "% of IOM N" agreed to 0.5%).

### Table 2 — gas yields (mol/kg sample). 48 h duration, confirmed p.3.

| Experiment | T C | P kbar | Sample | H2 | CH4 | CO2 |
|---|---|---|---|---|---|---|
| 250#3 | 250 | 1 | HC095 | 4.59e-7 | n.d. | 4.57e-5 |
| 250-1-1-Mar2023 | 250 | 1 | HC113 | 1.95e-7 | 4.04e-6 | 1.09e-3 |
| 250-1-2-Mar2023 | 250 | 1 | HC113 | 9.56e-7 | 5.27e-6 | 1.19e-3 |
| 250#2 | 250 | 3 | HC096 | 1.27e-6 | 8.83e-7 | 1.85e-4 |
| 250-3-1-Mar2023 | 250 | 3 | HC113 | 8.74e-7 | 3.47e-6 | 1.11e-3 |
| 250-3-2-Mar2023 | 250 | 3 | HC113 | 1.73e-6 | 4.76e-6 | 1.24e-3 |
| 350_1_1 | 350 | 1 | HC083 | 4.24e-6 | 1.77e-4 | 4.01e-3 |
| 350-3-M | 350 | 3 | **Murchison** | 6.00e-6 | 1.84e-4 | 3.04e-3 |
| 500-3-M | 500 | 3 | **Murchison** | 3.00e-5 | 1.57e-3 | 3.43e-3 |
| 500#4 | 500 | 10 | HC096 | 2.44e-5 | 3.18e-3 | 1.01e-3 |
| 500-10 | 500 | 10 | HC113 | 3.04e-5 | 3.40e-3 | 2.55e-3 |

Careful: these are `mmoles/mg`, so e.g. Murchison 350 C CO2 = 3.04e-3
mmol/mg = **3.04 mol/kg**. The exponents above are as printed in the table.

### Table 1 — starting compositions (partial)
- **Murchison IOM**: C 59.9 wt%, H/C 0.70, N/C 0.036, O/C 0.18, N 2.48±0.14 wt%
  -> per kg: C 49.87, H 34.91, O 8.98, N 1.77 mol
- **HC113, HC083, HC095, HC096**: NOT EXTRACTED — needed, see below.

### Table 5 — residue compositions (partial)
- **Murchison 350-3-M**: H/C 0.63±0.08, O/C 0.096±0.011
- **Murchison 500-3-M**: H/C and O/C are **BLANK in the table** (only N/C
  0.015 reported). This is a gap in the published data, not a transcription
  miss. Do not hunt for it.
- **HC113 250-1-1-Mar2023**: H/C 0.75, RH/C 0.67, N/C 0.033, RN/C 0.33,
  O/C 0.145, RO/C 0.54, C 74.3 wt%, H 4.67 wt%, N 2.83 wt%, O 14.4 wt%
- **HC113 250-1-2-Mar2023**: H/C 0.78, RH/C 0.69, N/C 0.033
- **HC095 250#3**: H/C 0.73, RH/C 0.76, N/C 0.045, RN/C 0.79, O/C 0.218,
  C 67.4 wt%, H 4.08, N 3.54, O 19.6
- Others not extracted.

The R-prefixed columns (RH/C, RO/C, RN/C) are residue-ratio / starting-ratio.
**Use these for cross-sample comparison** — dimensionless, so they don't
confound material differences with model error.

### Table 9 — NH3/NH4+ yields
- Murchison 350-3-M: 0.57 µmol/mg (= 0.57 mol/kg), 32.06% of IOM N
- Murchison 500-3-M: 0.72 µmol/mg (= 0.72 mol/kg), 40.83% of IOM N

### Table 10 — CO2 precursor Ea distributions
| Ea kcal/mol | BS89 % | HC113 (Miller) % |
|---|---|---|
| 42 | 5 | |
| 44 | 15 | |
| 46 | 25 | |
| 48 | 25 | |
| 50 | 15 | 5 |
| 52 | 10 | 15 |
| 54 | 5 | 25 |
| 56 | | 25 |
| 58 | | 15 |
| 60 | | 10 |
| 62 | | **5 (inferred)** |

The 62 entry was cut off in the source paste. Inferred as 5 because both
columns then sum to 100 and the distributions have identical shape
(5,15,25,25,15,10,5) — HC113 is BS89 shifted +8 kcal/mol. Near-certain but
**confirm against the paper**.

Context that matters: these are for Vitrimat run at **10 C/Myr geologic
heating**, fitted to **residue trajectories**, for **HC113 synthetic IOM**.
They are NOT interchangeable with our 48 h isothermal Murchison gas-yield
fits. Tested: HC113's distribution under-predicts Murchison's 350 C yield by
3.6x. Fig. 9b separately shows all distributions "shifted 22 kcal/mol higher"
to demonstrate what would match the lab data — a different claim from Table
10's +8; both are real and answer different questions.

### Current module calibration (single-Ea, to be replaced for CO2)
| channel | m0_per_kg | Ea kJ/mol | basis |
|---|---|---|---|
| IOM_CO2 | 4.169 | 216.2 | absolute 350 C Murchison yield |
| IOM_CH4 | 1.908 | 229.4 | absolute 350 C Murchison yield |
| IOM_N | 0.875 | 217.3 | absolute 350 C Murchison NH4+ |
| IOM_S | 0.462 | 217.3 | **UNCALIBRATED** — pinned to Ea_N, no data exists |

Constants: **A = 2e15 /s** (Burnham 2019 Table 1, CONFIRMED — supersedes the
1e13 assumed throughout the rest of this document), `logA = 15.301`,
R = 8.314, duration 48 h = 172800 s.

---

## WHAT THE USER CAN SOURCE (ask — don't guess these)

Mohit has the paper and the literature and has offered to supply numbers. The
following are either missing or assumed, and **guessing them will silently
corrupt the calibration** — this exact failure mode has already produced one
18% error and one 1000x error in this project.

**RESOLVED since this document was written** (all now in
`PATCH_F_vitrimat2018.md`):

1. ~~Vitrimat pre-exponential factor A~~ — **A = 2e15 s^-1**, Burnham 2019
   Table 1. Our 1e13 was the 1989 value.
2. ~~Table 1 starting compositions for syn-IOMs~~ — supplied, in the project as
   `Miller 2025 - Table 1.txt`. HC113: H/C 1.13, N/C 0.100, O/C 0.27, C 55.1
   wt%. HC083: H/C 0.83, N/C 0.073, O/C 0.38, C 54.5 wt%. HC095: H/C 0.95,
   N/C 0.057, C 58.1 wt%. HC096: H/C 0.96, N/C 0.057, C 67.3 wt%.
   Murchison: H/C 0.70, N/C 0.036, O/C 0.18, C 59.9 wt%, N 2.48 wt%.
3. ~~Whether HC113 has a 350 C run~~ — **yes**: `350-10-Mar2023`, 10 kbar,
   residue H/C 0.66, RH/C 0.58, O/C 0.047, RO/C 0.18. (Note it is 10 kbar, not
   3, so pressure is a confound when comparing to the 250 C 1/3 kbar runs.)
4. ~~Table 5 residues for syn-IOMs~~ — full table supplied as
   `Miller 2025 - Table 5.txt`. Use the `RX/C` columns.
5. ~~BS89 CH4-precursor distribution~~ — better than asked: **Vitrimat 2018**
   gives H2O, CO2, CHn and CH4 distributions (Burnham 2019 Table 1,
   reproduced in `PATCH_F_vitrimat2018.md` section 2).
6. ~~Table 10 "62 kcal = 5%" inference~~ — user confirmed correct.
7. **Murchison 500 C residue H/C and O/C: confirmed genuinely absent** from
   Table 5 (that row has only N/C, RN/C, C wt%, N wt%). Stop looking.

**Still open — ask the user:**

A. **Miller's Table 10 gives HC113 as +8 kcal/mol vs BS89; our gas-yield fit
   gives HC113 shift 0.0.** Both from Miller's own data on the same sample.
   His is a residue-trajectory fit at 10 C/Myr geologic heating; ours is a
   gas-yield fit at 48 h isothermal. Does not block production (our Murchison
   fit is independent of it) but should be understood before publication.
B. **Nitrogen composition, now three-way.** Our bulk formula uses N 2.142 per
   formula unit; an old comment said N 3.284; Miller Table 1 implies N/C 0.036
   = **N 2.18** on ENTICES's C 60.611 basis. So 2.142 is probably right and
   3.284 was an error — but the ENTICES bulk composition's actual source
   (possibly `Primordial_CI_chondrite_mineralogy.Rmd`) should confirm it.
   Outstanding for several sessions now.
C. **Do we want the H2O and CHn channels?** Vitrimat 2018 supplies both. H2O
   is chemically uninteresting here (goes to solvent) but **CHn is real organic
   carbon we currently ignore entirely** — an explicit decision, not an
   omission by default.
D. **Pressure.** HC113's 250 C runs are 1 and 3 kbar; its 500 C run is 10 kbar;
   its 350 C residue run is 10 kbar. Miller calls pressure "the weaker
   variable" but that is qualitative. If pressure suppresses CO2 yield, some of
   the apparent temperature trend is a pressure artefact.

## Verify, don't trust

Every number in this document was computed either by hand in Python or
transcribed from a text extraction of the PDF, in an environment with **no R
and no PHREEQC**. Nothing here has been executed as R code.

The project's track record on this is the reason for the warning: two
separate calibration errors have been found *only* by recomputing a quantity
independently rather than trusting the derivation — an 18% Ea error (ratio
fit applied to a differently-scaled m0) and the ~1000x low-temperature error
described above. Both looked fine in documentation.

So: re-derive the fits in Task 2 yourself rather than hardcoding my shift and
m0 values. If your numbers differ from the table above, **your numbers are
probably right** — say so, show the arithmetic, and flag the discrepancy
rather than matching mine.

Report discrepancies plainly. Mohit has asked explicitly for pushback over
agreement, and every real error in this project has surfaced that way.
