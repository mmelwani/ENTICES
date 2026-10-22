# Handover for review — work since PROGRESS.md's 2026-09-20 revision

Date: 2026-09-30
From: the VS Code (Claude Code) session working on the `iom-module` and
`lucas-integration` branches.
For: the Claude review/planning session that maintains `PROGRESS.md`,
`OPTION3_decision.md` and `NEXT_STEPS.md`.

**Everything below has had only self-review.** No second reader has checked
any of it. §7 lists the claims I'd most want independently re-derived, and §6
lists two places where I already caught and corrected my own overstatement —
so please treat the rest with the same skepticism.

This document is self-contained: it restates the numbers rather than
pointing at files you may not have. File names are given so you can ask for
any of them.

---

## 0. What in PROGRESS.md is now stale

| PROGRESS.md currently says | Actual state |
|---|---|
| Residue re-check "partially resolved", +11% → +7.8% (Part 4, Part 5 "Done") | **Resolved** — it's a starting-composition artefact; passes to +0.1% on `RX/C`. See §2. |
| Next #1: "Redox budget diagnostic… the only item still untouched" | **Done.** See §3. |
| Next #3: "Re-run Test B sims 3/4 with ENTICES's real secondary-phase list" | **Done** — but it can't size the mineral perturbation. See §5. |
| Next #4: "Wire the module into ENTICES's template generator" | **Done**, as a PR awaiting Lucas's review. See §4. |
| Next #6: "Chase the residual +7.8% residue overshoot" | **Done** (same as row 1). |
| Part 3 caveat: decoupling H2 → "Magnetite 2.9×, pyrrhotite 3.6×" | Still the **only** measurement of that effect, and it comes from a 2-mineral toy assemblage. The real-assemblage re-run couldn't measure it (§5). Treat as unconfirmed. |
| "Wire the module in" / "Lucas's code" generally | Lucas's actual code turned out to be 3 calibration generations stale and not sourcing the module at all. See §4. |

Nothing in `OPTION3_decision.md`'s per-channel decisions is contradicted.
Two **new** gaps of the same kind as the `Amm` finding have appeared (§5.3,
§5.4) and should probably be added to it.

---

## 1. Git / GitHub state

Everything was local-only until 2026-09-26 — which is why Lucas had never
seen any of it. Now pushed to `github.com/mmelwani/ENTICES`:

- **PR #1** `iom-module` → `main`: the IOM module, self-test, all calibration
  work through 2026-09-20, and the docs. 27/27 self-test checks passing.
- **PR #2** `lucas-integration` → `iom-module` (base is `iom-module` because it
  depends on #1; retarget to `main` after #1 merges): Lucas's script added
  unmodified in one commit, the wiring change in a separate commit, plus the
  Test B3 re-run. Needs Lucas's review before merging — it changes his
  production calibration.

Lucas has since reported (2026-09-30): *"the generator code works for me and
the pqi input file looks good"* — so the wiring is confirmed from his side.
He cannot yet **run** the generated input because neither database is
complete (§5.1).

Lucas's own fork (`lucasmfifer-code/ENTICES`, branch `experimental`) is an
**orphan history** — root commit with no parent, containing only his two
files. It can't be merged; its script is byte-identical to what he emailed,
so nothing was lost by working from the emailed copy.

---

## 2. The 350 C residue overshoot is a starting-composition artefact

`RESIDUE_OVERSHOOT_RESOLVED_2026-09-20.md`, reproduced by `residue_ratio.R`.

The distributed-Ea + CHn model predicts 350 C residue H/C 0.6785 against
Miller's measured 0.63 (+7.7%). But the model's **starting** H/C is already
7.6% above Murchison's:

    model/measured STARTING  H/C = 0.7530 / 0.70 = 1.0757
    predicted/measured FINAL H/C = 0.6785 / 0.63 = 1.0770
    contributed by the decomposition model       = +0.13%

On the dimensionless ratios (Miller Table 5, row `350-3-M`):

| | model | measured | error |
|---|---|---|---|
| RH/C | 0.9011 | 0.90 ± 0.12 | +0.1% |
| RO/C | 0.5544 | 0.53 ± 0.06 | +4.6% (inside uncertainty) |

So on the comparison your documents already specify as correct, the check
**passes**. It is not evidence for or against the low-Ea tail — the offset
exists at zero reaction extent.

**Where the offset comes from:** `iom_bulk_formula()` (C 60.611, H 45.639,
O 9.939, N 2.142, S 1.154) traces to `Primordial_CI_chondrite_mineralogy.Rmd`
~line 1316, which states target `C60.611H45.639O9.941N3.284S1.154` and
decomposes it into KerogenC128 + pyridine + S + O. So the bulk composition is
**CI-derived**, while the kinetics are calibrated to **Murchison (CM2)**.
Possibly deliberate — but absolute residue agreement can't be quoted while it
stands. **Decision for the PI:** keep CI and always report `RX/C`, or switch
the bulk formula to Murchison.

**Related, and it closes a question `PATCH_F_addendum` left open:** that Rmd
still carries **N3.284**, the value the addendum traced to a misread
hydrogen gram figure. So the source disagrees with the module (N 2.142).
Not changed — the Rmd solves its phase coefficients from that target, so
editing N there changes the mineral endmember breakdown.

---

## 3. Redox budget diagnostic — done without Lucas's runs

`REDOX_BUDGET_DIAGNOSTIC_2026-09-20.md`; inputs and `compare.R` in
`phreeqc_inputs/redox_budget/`.

`phreeqc_inputs/` already held a real ENTICES-generated input at 15 wt%
organics (`enceladus_k1.00e-12_d600_t7.50e+04_circ1.pqi`). Its exact
organics-off twin was made by deleting only the 9-line `REACTION 1` block
(verified by diff). Both run clean in PHREEQC 3.8.6.

| | organics OFF | organics ON |
|---|---|---|
| pH | 8.617 | 5.792 |
| pe | −7.334 | −3.101 |
| m(H2) | 2.36e-06 | 3.62e-09 (**653× lower**) |
| m(Fe2+) | 1.44e-06 | 5.66e-05 (**39× higher**) |
| kinetic primary minerals (incl. Magnetite, Pyrrhotite) | — | **no measurable change** (largest: Fayalite −0.018%) |
| secondary | Greenalite, Cronstedtite-7A, Boehmite present | all three **absent**; Pyrite 56× higher; Beidellite-Fe appears |

**Confirmed:** H2 depletion and a +4.2 pe shift — the Part 2b electron-sink
signature.

**Not confirmed in the form Part 2b poses it:** "electron demand satisfied by
Fe(II) mineral oxidation" predicts Fe2+ falling; it rises 39×, and the
kinetic Fe minerals don't move. My interpretation — offered as an
interpretation, not a demonstrated mechanism — is that the organics also
drop pH by ~2.8 units, and the vanishing phases (Greenalite, Cronstedtite,
Boehmite) are ones that destabilise at lower pH, releasing Fe that would
otherwise be locked in clays. So it's a **compound perturbation** (pH shock +
H2 depletion), not a clean redox pair. Separating the two needs a
pH-buffered run, which I didn't do.

**Caveats that limit what this shows:**
- It uses the **old inline** organic treatment (elemental `REACTION` block with
  aromatic/thiophene pseudo-species), not the calibrated module. It answers
  "is existing output contaminated" (yes), not "how does the new module
  behave".
- It's a single `-steps 1` step of ~75 kyr, and the kinetic primary minerals
  barely react in *either* run — so "no primary-mineral difference" may just
  mean too little reaction for a difference to show.
- The old inline block injects its own H2 (0.0106 mol), so the IOM isn't a pure
  electron sink in that treatment.

---

## 4. Lucas's code, and wiring the module in

`FROM_LUCAS_20260924_review.md`.

Lucas sent `ENTICES_LF_20260917.R`, an example run
(`enceladus_k1.00e-07_d194200_T150_org15_orb12.pqi` + `.out`) and his database
`coreclath_CH4coupled.dat`.

**His script did not source `iom_module.R`.** It had its own hardcoded channel
table, repeated in three places (main path, bulk restart, tidal restart), with
the original PATCH_C values:

    m0_per_kg = c(4.169, 1.908, 0.875, 0.462)
    Ea_J      = c(214e3, 228e3, 215e3, 215e3)
    logA      = c(13, 13, 13, 13)

i.e. pre-2026-09-17: wrong A, single-Ea CO2/CH4, and the CH4 value that
carried the original 18% ratio-fit bug. His own comment flagged it as
provisional ("Sub-pools to follow").

**Real-output evidence of the Part 2 artefact.** In his example run, porewater
at 150 C right after the IOM channels react:

    pe  = −11.692
    CH4 = 4.327e-07 mol/kg
    CO2 = 2.551e-16 mol/kg      → CH4/CO2 ≈ 1.7×10⁹

That's the coupled-carbon-to-methane collapse, seen in a production run
rather than predicted.

**The wiring change** (PR #2, commit `d59b630`): added
`source("iom_module.R")` and replaced all three hardcoded tables with
`iom_default_config()`. Nothing else in his script was touched. I left his
hand-rolled RATES/KINETICS/USER_PUNCH/restart loops alone because they
already iterate generically over `self$iom_pools`' rows; his PUNCH lines use
`200 + 10*j`, and the template's own mineral PUNCH lines stop at 110, so 29
channels (lines 210-490) don't collide.

**Verified by execution, not reading:** instantiated his R6 classes and called
`generate_phreeqc_input_tidal_kinetic()` with his example parameters (k=1e-7,
d=194200, T=150 C, 15 wt%, 12 orbits, ungrouped). The generated `.pqi` has 29
IOM channels in KINETICS, 29 RATES entries all at `logA = 15.30103` (zero at
`logA = 13`), and 29 `<channel>_mol` USER_PUNCH headings. One pre-existing
warning (a missing salinity-profile text file) is unrelated.

Worth noting: his restart design (punch `IOM_*_mol`, rebuild KINETICS from
recovered moles) independently matches what `iom_module.R`'s INTEGRATION
POINTS block recommends.

---

## 5. Databases, and the real-assemblage Test B re-run

`PHREEQC_TEST_B3_real_assemblage_2026-09-27.md`; input
`test_B3_real_assemblage.pqi`.

### 5.1 Full diff: `Core11_idealgas_mod_v4.dat` (Mohit) vs `coreclath_CH4coupled.dat` (Lucas)

A full section-by-section diff, not a keyword search:

- **PHASES:** 402 shared. Mohit's only: `Fluorapatite`, `Fluorite`. Lucas's
  only: `CH4_hydrate`, `CO2_hydrate`, `H2S_hydrate`.
- **SOLUTION_MASTER_SPECIES:** differ only by `F`.
- **SOLUTION_SPECIES:** Mohit's adds the F⁻ complexes (CaF⁺, HF, HF2⁻, MgF⁺,
  MnF⁺, NaF, HPO3F⁻, H2PO3F, PO3F⁻², SiF6⁻²). **Plus one difference nobody had
  flagged:** the `CO2 = CO2` identity species (sets CO2(aq)'s `-Vm`, PS01) is
  **commented out in Mohit's, active in Lucas's**. `log_k = 0` either way, so
  only the pressure/volume correction differs — small but real at ENTICES
  pressures.
- **RATES:** functionally identical. (Two rate-law lines differ only in
  comments; Fayalite's acid-term `^1` vs `^1.0` is the same fix made
  independently on both sides.)
- Both carry the full `Hdg`/`Oxg`/`Mtg`/`Ntg`/`Sg` set.

**Consequence Lucas has now hit:** neither database is runnable for the
target configuration — Mohit's lacks the clathrates, Lucas's lacks the
fluorine species. He is merging them. **Which merged file becomes canonical,
and whether the `CO2 = CO2` definition should be active, needs a decision.**

### 5.2 Test B re-run on the real assemblage

The original Test B "realistic" cases used Magnetite + Pyrrhotite only.
Re-run with Lucas's database, the **real 78-phase ocean secondary list**
(`EQUILIBRIUM_PHASES 2`, extracted from his actual output — the porewater
list `EQUILIBRIUM_PHASES 1` has 81 and differs by −Fe, −Ni, −Lawrencite,
−Quartz, +Magnetite, +3 clathrates), and his real ocean pressure (66.7 atm;
the original tests set no pressure).

Note `Pyrrhotite` is in **neither** real secondary list — in production it's
kinetic (via `Troilite`), not something that precipitates.

| | coupled, no minerals | Hdg, no minerals | coupled, real 78-phase | Hdg, real 78-phase |
|---|---|---|---|---|
| pe | −9.204 | +3.000 | −9.136 | −8.283 |
| % of C → CH4 | 12.5% | 0% | 10.75% | 0.0000607% |

Suppression with the real assemblage ~177,000× (toy: ~26,000×); pe shift
+0.85 (toy: +0.84). **The core `Hdg` result holds.**

**What it does NOT show — please weigh this carefully:**
- Magnetite never saturates (SI −6.6 to −9.7). Pyrite is the only phase that
  does, and it precipitates **exactly 1.000e-05 mol in both the coupled and
  decoupled case** — all the Fe (only 1e-5 mol was supplied). The test is
  **Fe-limited**, so the mineral transfers are identical by construction.
  Therefore:
  - the pe difference is not caused by any mineral difference;
  - **it cannot measure the "~3× Fe/S perturbation" PROGRESS.md quotes.** That
    figure remains a toy-assemblage result, unconfirmed.
- The suppression factor is sensitive to the arbitrary small inventories
  chosen (1e-3 C, 1e-5 Fe, 1e-4 S); report it as "several orders of
  magnitude," not as a precise number.
- Calcite check (Hdg + Ca): SI 1.33 → 0.00, **8.055e-04** mol precipitated.
  Carbon still sees carbonates under `Hdg`.

### 5.3 New gap: clathrates can't see decoupled methane

Found by reading, not running. Lucas's
`CH4_hydrate: CH4:6H2O = CH4 + 6 H2O` uses **coupled** `CH4`. Methane held as
`Mtg` shares no mass balance with `CH4`, so if carbon decoupling is ever
switched on, clathrate formation won't see it. B3 couldn't test whether this
matters: every clathrate stayed far undersaturated (SI −2.8 to −107) at 1e-3
mol/kgw C. Whether production ocean concentrations approach clathrate
saturation is a question about the full model.

### 5.4 New gap, raised by Lucas today: sulfide minerals can't see decoupled sulfide

Lucas asked whether `HS-` in mineral PHASES definitions should be swapped to
`HSg-` for a decoupled-sulfur database. Checked:

    Troilite    FeS + H+ = Fe+2 + HS-
    Pyrrhotite  FeS + H+ = Fe+2 + HS-
    Pyrite      FeS2 + H2O = 0.25 H+ + 0.25 SO4-2 + Fe+2 + 1.75 HS-

25 PHASES reference `HS-` (Pyrite, Pyrrhotite, Troilite, Sphalerite,
Chalcopyrite, Bornite, Covellite, …, plus generic entries like Th2S3, EuS,
CrS that are irrelevant to Enceladus).

**My answer to him (not yet sent as a decision, pending review):** don't swap
in place. `Troilite` is the kinetic primary sulfide in production and
releases **coupled** `HS-`; rewriting Pyrite/Pyrrhotite to `HSg-` would stop
them responding to the rock's own sulfide. It's the same class of gap as
`Amm` and §5.3: a decoupled pseudo-element is invisible to anything written
in its coupled counterpart. As things stand, `IOM_S`'s decoupled sulfide has
**no mineral sink** and will accumulate as `H2Sg`/`HSg-`.

Whether that matters is a design question: FeS/pyrite precipitation from H2S
is fast even at low T and isn't the S(−2)↔S(+6) redox step decoupling exists
to block — so zero sink is unrealistic. But `IOM_S` is already uncalibrated
and small. If a sink is wanted, the fix is **parallel** phases (e.g.
`Pyrite_Sg`, `Pyrrhotite_Sg`) written in `HSg-`, for the 2-3 minerals that
matter, not editing the existing 25.

**One more thing worth your check:** `Sg`'s master species is `H2Sg`, so a
parallel `Pyrite_Sg` would be a reaction in `HSg-`/`H2Sg`. PHREEQC phases are
single-reaction, so `Pyrite` and `Pyrite_Sg` would be separate phases with
separate SIs — two pyrite pools that can't exchange. Whether that's
acceptable physically (it's the same compromise the `Mtg`/`CH4` split already
makes) is worth a sentence in `OPTION3_decision.md`.

---

## 6. Corrections I made to my own work during this period

Listed so you know where my first drafts over-reached:

1. **B3 "Pyrite buffers pe".** My first write-up (and the PR #2 description,
   and the commit message `4966349`) said Pyrite "does the actual
   pe-buffering" and that the toy test "got the right answer for the wrong
   mechanism." That's not supported: pyrite takes all the Fe in both cases,
   identically. The B3 doc has been corrected (2026-09-30); **the commit
   message and PR #2 description still contain the overstatement.**
2. **B3 calcite number.** I quoted 8.16e-04 (the B2 value at no pressure); the
   66.7 atm value is 8.055e-04. Corrected in the doc.

---

## 7. What I'd most like independently re-derived

In order of how much depends on each:

1. **`Ea_N = 242.8` kJ/mol (the basis-bug fix).** Claim: fitting
   `m0 = 0.875` (ENTICES basis) against Miller's raw 0.5676 mol/kg (Murchison
   basis) mixed two bases differing by 2.142/1.7705 = 1.2098. The correct
   target is the ratio 0.5676/0.7229 = 0.7852 → Ea 242.8, not 244.8.
   Consequence claimed: the "~20-22% overshoot at 500 C" long attributed to
   exhaustion-pool sizing was this basis mismatch. Note the 500 C ratio of
   1.00 after correction is **by construction** (m0 was calibrated with 500 C
   as the exhaustion proxy), not independent validation.
2. **The Mtg/Sg asymmetry** (`Mtg 1 H 0.435` right, `Sg 1 H 1` right, opposite
   rules because `Mtg`'s master species has no H and `Sg`'s is `H2Sg`).
   Measured in `test_A4_formulas.pqi`; the criterion differs between the two
   channels (aqueous-state match for Mtg, ΔH match for Sg) — check that I
   chose the right criterion for each.
3. **The residue argument (§2).** Specifically, whether differencing starting
   vs final ratios is legitimate given unreacted channel mass is counted as
   solid.
4. **The redox-budget interpretation (§3).** I attribute the Fe response to pH
   without having run a pH-buffered control.
5. **The wiring's "minimal diff" claim (§4)** — that Lucas's hand-rolled
   generators are genuinely row-generic everywhere, including the restart
   paths I only read, not executed.

---

## 8. Decisions pending (collected)

1. CI vs Murchison bulk composition — and with it, the Rmd's N3.284.
2. `Hdg` mineral-perturbation: accept and report, or restrict to
   organics-active cases? (Its magnitude is still only a toy-test number.)
3. `Hdg` globally vs below a temperature threshold.
4. Whether to actually switch Option 3 on in `ENTICES_LF_20260917.R` (the
   wiring keeps default all-coupled).
5. Canonical merged database, and whether `CO2 = CO2` is active.
6. Decoupled-sulfide mineral sink (§5.4): none, or parallel `_Sg` phases?
7. Decoupled-methane clathrate (§5.3): same question for `CH4_hydrate`.

---

## 9. Files, if you want to ask for any

New docs: `FROM_LUCAS_20260924_review.md`,
`PHREEQC_TEST_B3_real_assemblage_2026-09-27.md`,
`RESIDUE_OVERSHOOT_RESOLVED_2026-09-20.md`,
`REDOX_BUDGET_DIAGNOSTIC_2026-09-20.md`.
Scripts/inputs: `residue_ratio.R`, `test_B3_real_assemblage.pqi`,
`phreeqc_inputs/redox_budget/{org_on,org_off}.pqi`, `compare.R`.
Code: `iom_module.R`, `iom_selftest.R`, `ENTICES_LF_20260917.R` (wired).
Databases: `Core11_idealgas_mod_v4.dat`, `coreclath_CH4coupled.dat`.
