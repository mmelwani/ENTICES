# Tidal-kinetic convergence diagnostic: for review

Date: 2026-10-06
From: the VS Code (Claude Code) session working on the `rock-composition` branch.
For: the Claude review/planning session.

**Everything below has had only self-review.** §6 separates what was measured
from what I inferred. §7 lists what I'd most want re-derived.

The document is self-contained: numbers are restated, not referenced.

---

## 0. Context

### State of the branches

- PRs #1 (`iom-module`) and #2 (`lucas-integration`) are merged to `main`.
- Lucas pushed four commits to `rock-composition` and opened PR #3, which isn't merged yet.
  - The commits are `ca8368f`, `c84775b`, `236b42b` and `9a7ced5`; the branch tip is `9a7ced5`.
  - They contain:
    - a tephroite rate law (Casey 1993, via Pincus 2026);
    - `coreclath_CH4uncoupled.dat`;
    - `mineral_dissolution_constants.r`;
    - CHANGELOG updates.
- Before those, Lucas's commit `beef340` was verified in the 2026-10-02 session:
  - **What it changed:** bulk_kinetic now punches the IOM columns, and `depletion_order` was rewritten.
  - **Bulk run then bulk restart:** both run in PHREEQC with 0 errors.
  - **IOM carry-over:** the 29 IOM pools carry over exactly.
  - **Decay:** the fastest pool decays at 150 °C.

### Lucas's report

Lucas reports very slow convergence, or none, in `generate_phreeqc_input_tidal_kinetic`. His CHANGELOG attributes it to the K-feldspar/albite base-weathering mechanism: "pH runs away to ~11.5 within the first step, triggering millions of Newton-Raphson iterations."

He proposes two guards:
- **(a)** a minimum orbit step size;
- **(b)** automatically moving some primaries to equilibrium at t = 0, based on T and pH.

---

## 1. Test setup

- **Code and tools:** PR #3 tip, PHREEQC 3.8.6, `coreclath_CH4coupled.dat`.
- **Generator:** `generate_phreeqc_input_tidal_kinetic`, called with:

  | Argument | Value |
  |---|---|
  | `porewater_temp` | 160 |
  | `n_orbital_cycles` | 2 |
  | `mode` | `"ungrouped"` |
  | `substeps` | 1 |
  | `organic_wt_percent` | 15 |
  | `CH4_redox_override` | `"coupled"` |
  | `k` | 1e-7 |
  | `layer_thickness` | 194200 |

- **Porewater (SOLUTION 1):** 0.1574 kg of pure water at pH 7, pe 4, 160 °C and 70 atm.
- **Rock:** the scaled rock (×0.85) is in contact with that porewater.
  - KINETICS 1 has 38 reactants:
    - 8 silicate/sulfide primaries;
    - schreibersite;
    - 29 IOM sub-pools.
  - EQUILIBRIUM_PHASES 1 holds Fe 1.998 mol, Ni 0.2388 mol and Lawrencite 0.00956 mol, all `dissolve_only`.
  - It also holds about 75 secondary phases at `0 0`, including 11 gases.
- **Step size:** one orbit is a single KINETICS step of 118,800 s, solved with `-cvode true` and `-bad_step_max 5000`.

---

## 2. Result 1: reproduced; no single ingredient is the cause

The baseline did not finish the **first** orbit's kinetic step within 10 minutes. By then it had done about 16,000 equilibrium solves. Each solve took about 40 iterations, with `Overflow` in revise_guesses and infeasible-solution messages.

Then I removed one ingredient at a time:

| Variant (one change from baseline) | Finished orbit 1 within cap? |
|---|---|
| No IOM (`organic_wt_percent = 0`) | No (10 min) |
| Fe and Ni set to 0 mol | No (10 min) |
| Schreibersite m0 = 0 | No (10 min) |
| Albite and K-feldspar `kbase = 0` | No (10 min) |
| `-cvode false` (Runge–Kutta) | No (10 min; about 184,000 solves) |
| All 11 gas phases removed from EP1 | No (15 min) |
| Magnetite added to EP1 | No (15 min) |
| Magnetite added and gases removed | No (15 min) |

So the slowness doesn't come from IOM, metal, schreibersite, the base mechanisms or the integrator choice alone.

---

## 3. Result 2: where the time goes (first orbit, sub-stepped)

To see inside the first orbit, I split it into log-spaced sub-steps that sum to 118,800 s. The increments were 1, 1, 3, 5, 10, 30, 50, 100, 300, 500, 1000, … s, with `INCREMENTAL_REACTIONS true` and the logfile off. Each sub-step writes one output row.

Two runs, each capped at 15 minutes:

**Baseline**

| Cumulative t (s) | pH | pe | H₂O (kg) | H₂(g) (mol) | Albite (mol) | K-fsp (mol) | Forsterite (mol) |
|---|---|---|---|---|---|---|---|
| 1 | 4.54 | −4.48 | 0.1034 | 2.997 | 0.2705 | 0.01744 | 1.3658 |
| 200 | 6.09 | −6.03 | 0.1029 | 2.991 | 0.2700 | 0.01744 | 1.3585 |
| 500 | 8.80 | −8.74 | 0.1022 | 2.961 | 0.2560 | 0.01743 | 1.3548 |
| 1000 | 11.48 | −11.42 | 0.0909 | 2.381 | 0.0125 | 0.01569 | 1.3502 |
| (cap reached) | | | | | | | |

**Magnetite added to EP1, gases removed**

| Cumulative t (s) | pH | pe | H₂O (kg) | H₂(aq) (molal) | Albite (mol) | K-fsp (mol) | Forsterite (mol) |
|---|---|---|---|---|---|---|---|
| 1 | 4.91 | −6.35 | 0.1544 | 1.08 | 0.2705 | 0.01744 | 1.3659 |
| 200 | 6.06 | −7.50 | 0.1540 | 1.08 | 0.2700 | 0.01744 | 1.3591 |
| 500 | 8.92 | −10.36 | 0.1530 | 1.09 | 0.2592 | 0.01743 | 1.3541 |
| 1000 | 11.18 | −12.97 | 0.1362 | 4.68 | 0.0374 | 0.01645 | 1.3496 |
| 2000 | 11.16 | −12.96 | 0.1350 | 4.83 | 0.0001 | 0.01381 | 1.3406 |
| 5000 | 11.18 | −12.98 | 0.1328 | 4.91 | 0 | 0.00440 | 1.2709 |
| (cap reached) | | | | | | | |

### Reading of the tables

- **About 500 s of model time took about 10 minutes of wall-clock time.** The 500–1000 s interval alone took that long, in both runs.
- **Then pH rises and albite collapses.**
  - Forsterite dissolution raises pH from about 6 to 9.
  - Albite's base mechanism then drives pH to about 11.2–11.5.
  - About 95% of the albite dissolves within that same 500 s.
- **K-feldspar follows**, about 75% gone by 5,000 s. Forsterite then speeds up.
- **In the baseline, 0.58 mol of H₂(g) is consumed in the 500–1000 s interval.** I infer, without checking the phase output, that this is H₂ reducing the t = 0 hematite into Fe(II) silicates once silica becomes available.
- **In the gas-free variant, 0.50 mol of Fe metal corrodes** in the same interval, and dissolved H₂ climbs to about 4.9 molal.

**Conclusion.** Lucas's identification of the pH runaway driven by albite and K-feldspar is correct. But it isn't a numerical artefact of the base mechanism. It's a real transient: these phases dissolve in less than about 1.5 h at 160 °C, against a 33 h step. Each equilibrium solve is expensive throughout, and the hard work falls in that window.

### Hand estimate of timescales (not from PHREEQC)

The specific rate is (kacid + kneut + kbase)·SSA·mw, giving s⁻¹ per mol present.

- Albite at pH 9, 160 °C: kbase ≈ 10⁻⁶·⁵⁸ mol m⁻² s⁻¹, so the specific rate is about 3.4e-4 s⁻¹ and τ is about 3,000 s.
- K-feldspar at pH 11.5, 160 °C: the specific rate is about 3.5e-4 s⁻¹ and τ is about 3,000 s.

These are consistent with the sub-stepped results, and short compared with an orbit.

---

## 4. Result 3: two setup problems that are not convergence problems

### 4.1 Magnetite is absent from the porewater's EQUILIBRIUM_PHASES 1

It is present only in EP2 (ocean).

Lucas's CHANGELOG (2026-09-17/18) says Magnetite was "demoted to an ordinary secondary… added to whichever secondary mineral list each function was missing it from". In `tidal_kinetic`'s EP1, it isn't there.

### 4.2 Gas phases at `SI = 0` give each gas an unlimited reservoir at 1 atm partial pressure

The system is at 70 atm. To isolate the effect, I equilibrated SOLUTION 1 with EP1 only, with no kinetics:

| EP1 as used | Fe metal reacted | Products | H₂O left (kg) | pH | pe |
|---|---|---|---|---|---|
| As generated (gases at 0, no Magnetite) | 1.998 of 1.998 mol | Hematite 0.999 mol; H₂(g) 2.997 mol | 0.1034 | 4.52 | −4.46 |
| Plus Magnetite | 1.998 of 1.998 mol | Magnetite 0.666 mol; H₂(g) 2.664 mol | 0.1094 | 4.35 | −4.29 |
| Plus Magnetite, gases removed | 0.126 mol (6%) | Magnetite 0.042 mol; H₂(aq) 1.08 molal | 0.1544 | 4.90 | −6.34 |

These are the consequences for every tidal_kinetic run on this generator:

- **At t = 0, before any kinetics, all the metal corrodes.** It goes to Fe(III) (hematite) and consumes about a third of the porewater. That's because the 1 atm H₂ sink never lets H₂ build up.
- **The low pH at t = 0 comes from lawrencite.** Its Cl is released as HCl-equivalent once the Fe goes into hematite. The fluid has no other buffering.
- **Removing the gases is not the fix.** Dissolved H₂ then climbs to about 4.9 molal (§3), well above any plausible solubility at 70 atm.
  - The likely fix is gas phases at SI = log₁₀(P_total) ≈ 1.85, or a fixed-pressure GAS_PHASE.
  - That is a modelling decision. I haven't tested either.

**Possible tie-in.** This may also bear on the earlier finding that production H₂ comes from coupled metal corrosion. If the gas sink is at 1 atm, the H₂ budget in production runs is set by this artefact. I haven't checked how many of Lucas's production outputs are affected.

---

## 5. Assessment of Lucas's proposed guards

- **(a) Minimum orbit step size: wrong direction.**
  - The first orbit needs finer sub-steps, not coarser ones.
  - The log-spaced stepping in §3 is a cheap way to get there. It makes progress visible and lets the run be cut at known times.
  - But it didn't make the run finish faster.
- **(b) Primaries to equilibrium at t = 0: defensible, with a better criterion.**
  - The data support it: albite and K-feldspar are gone in under about 1.5 h at 160 °C.
  - I suggest switching a primary to `dissolve_only` equilibrium when its specific dissolution timescale at the step's T is much shorter than the step length (say τ < step/100). That would replace a fixed T/pH table.
  - **Unresolved problem: pH is an output, so the switch criterion can't use it directly.** The base mechanism dominates only after pH has run away. One option is to evaluate τ at the maximum plausible pH; another is a short probe run.
- **Fix §4 first.** The t = 0 state (Fe(III), pH 4.5, a third of the water gone) is the state the solver starts from. Changing it may change the transient.

---

## 6. Measured vs inferred

| Measured (PHREEQC output) | Inferred / not checked |
|---|---|
| No variant finished orbit 1 within its cap (§2) | That per-solve cost, not CVODE step count, is the main cost driver. Suggested by about 40 iterations plus `Overflow` per solve, but I didn't profile |
| All the trajectory values in §3 | That the 0.58 mol H₂ loss is hematite reduction. I didn't inspect phase output at that step |
| All the t = 0 values in §4 | That fixing §4 would change the convergence picture. Untested |
| Magnetite absent from EP1 (grep of the generated .pqi) | Whether other generators have the same EP1 omission. Not checked |
| | The τ estimates in §3 are hand calculations from the PK04 constants in the RATES blocks |
| | `coreclath_CH4uncoupled.dat` (new in PR #3) is untested |

---

## 7. What I'd most like re-derived or challenged

1. **The 1 atm interpretation of gas phases at `SI = 0` in a 70 atm EQUILIBRIUM_PHASES block, and its consequence for the t = 0 metal.** Is SI = log₁₀(70) the right replacement? Or should the porewater be closed to gas, with H₂ handled some other way? This is the finding with the widest reach.
2. **Whether letting both feldspars dissolve within an hour at 160 °C is physically reasonable** with the PK04 base-mechanism extrapolation and SSA = 5 m²/g. Or is the SSA × mw × M scaling (Lucas's 2026-09-10 change) inflating it?
3. **The water-to-rock setup.** 0.157 kg of unbuffered pure water against about 0.85 kg-equivalent of reactive rock, at pH 7 and pe 4. Is that a reasonable initial porewater, or does it guarantee a violent t = 0?
4. **The criterion proposed in §5(b)**, and how to break the dependence on pH.

## 8. Smaller items from PR #3, for completeness

These are the tephroite rate-law issues. None of them is linked to convergence: Lucas found tephroite's ΔMoles was about 5e-14, and my no-change ablations didn't involve it.

- **No supersaturation guard.** The rate law has no `IF SR("Tephroite") > 1 THEN rate = 0`, so it can precipitate. Every other primary has the guard.
- **Double-counted rate at pH 7.** `kneut` uses the same pre-exponential, with ACT(H⁺) fixed at 1e-7, and is *added* to `kacid`. So the rate at pH 7 is twice the intended value, rather than flat at it.
- **Only in one generator.** The rate law is only in `tidal_kinetic`. The other three generators still carry the k = 0 placeholder, as Lucas's CHANGELOG notes.

## 9. Files

- The diagnostic inputs and outputs are in the session scratchpad (`conv/`), not in the repo. Available on request:
  - the baseline .pqi;
  - each variant's single-line diff;
  - the sub-stepped CSVs.
- The repo state is PR #3's tip, `9a7ced5`. This document is the only new file, and it isn't committed.
