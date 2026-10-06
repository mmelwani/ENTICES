# PHREEQC test results — Tests A and B, run locally

Date: 2026-09-20
Run with: PHREEQC 3.8.6 batch (`C:\Program Files\USGS\phreeqc-3.8.6-17100-x64\bin\Release\phreeqc.exe`)
against `Core11_idealgas_mod_v4.dat`.
Answers: `OPTION3_decision.md` §4(a) and §4(b), `PROGRESS.md` Part 3's two
open questions, `NEXT_STEPS.md` items 2 and 3.

**Headline: question (b) — does `Hdg` stop the carbon runaway — is CONFIRMED.
Question (a) is answered, and it corrects `OPTION3_decision.md`: the proposed
`Sg 1 H -1` is wrong by 2 mol H per mol S. `Mtg 1 H 0.435` is right.**

Corrected/extended input files are in the repo: `test_A4_formulas.pqi`,
`test_B2_Hdg_carbon.pqi`. Raw outputs in `phreeqc_test_output/`.

---

## 0. Two defects in the supplied test files, fixed before anything else

Both are mechanical, but the first invalidated Test B's headline comparison.

**(i) `DATABASE` must be the first keyword.** Both `.pqi` files put `TITLE`
first, so PHREEQC refused to start: *"ERROR: DATABASE must be the first
keyword in the input file."* Moved to line 1.

**(ii) Test B's control case never equilibrated redox.** Cases 1 and 2 were
INITIAL SOLUTION calculations only. PHREEQC honours input valence-state
totals in an initial solution and does **not** force cross-element redox
equilibrium there, so pe stayed at the input default of 4 and **no carbon
reduced in the coupled case** — the control failed to reproduce the artefact
it existed to reproduce. A batch-reaction step is required before PHREEQC
reports `pe ... Adjusted to redox equilibrium`. (Cases 3/4 accidentally did
this, which is why they showed the effect.)

**(iii) Test B's "realistic" mineral cases had no minerals.**
`Magnetite/Pyrrhotite/Fe(OH)2` were declared `0 0` — SI 0 and **zero moles** —
in a solution containing no Fe and no S. PHREEQC reported *"Element not
present"* for all three and they did nothing. So the claim in
`OPTION3_decision.md` §2 that with `Hdg` "pe is then set by the mineral
assemblage" was **not** tested by that file. `test_B2` gives them real mole
inventories and puts Fe and S in solution.

---

## 1. Test A — the decoupled formulas. One is right, one is wrong.

### The mechanism, which the design doc did not have

The two pseudo-elements behave **oppositely**, because of how their master
species are defined in `Core11_idealgas_mod_v4.dat` (lines 340-341):

    Mtg   Mtg    0   Mtg     16.032   # master species has NO H in its formula
    Sg    H2Sg   1   H2Sg    34.08    # master species IS H2Sg -- contains H2

- `Mtg`'s hydrogens are **outside** the H mass balance. Adding `Mtg 1` adds no
  H and draws none from water.
- `Sg`'s hydrogens are **inside** it. Adding `Sg 1` forms H2Sg by pulling two
  H out of the solution, and adds no H of its own.

`OPTION3_decision.md` §3 assumed the *same* rule for both ("Mtg carries 4 H",
"H2Sg carries 2 H", therefore subtract). That is right for Mtg and wrong for
Sg.

### Measured: H contributed by the formula, per mole of channel reaction

From `TOTMOLE("H")` differenced against a zero-extent-reaction baseline
(`test_A4_formulas.pqi`, reducing conditions, pH 11.5, 1.45 °C):

| formula | ΔH (mol H per mol) | verdict |
|---|---|---|
| `S 1 H 1` (coupled reference) | **+1.000000** | — |
| `Sg 1 H -1` ← design doc's proposal | **−1.000000** | ✗ wrong by 2 |
| `Sg 1 H 1` | **+1.000000** | ✓ **matches coupled** |
| `Sg 1` | +0.000000 | ✗ short by 1 |
| `C 1 H 4.435` (coupled reference) | **+4.435000** | — |
| `Mtg 1 H 0.435` ← design doc's proposal | +0.435000 | ✓ **correct, see below** |
| `Mtg 1 H 4.435` | +4.435000 | ✗ double-counts CH4's H |

### Why `Mtg 1 H 0.435` is correct despite the ΔH mismatch

ΔH is the wrong criterion for Mtg, because CH4's four hydrogens are bound in
the molecule in *both* modes — counted in `TOTMOLE("H")` when coupled,
invisible when decoupled. The criterion that matters is whether the **aqueous
state** matches. It does, exactly:

| | coupled `C 1 H 4.435` | decoupled `Mtg 1 H 0.435` | decoupled `Mtg 1 H 4.435` |
|---|---|---|---|
| m(H2) | 7.175e-04 | **7.175e-04** | 2.718e-03 |
| pH | 11.5000 | **11.5000** | 11.5000 |
| pe | −11.425 | **−11.425** | −11.714 |
| methane | CH4 1.000e-03 | Mtg 1.000e-03 | Mtg 1.000e-03 |

`Mtg 1 H 0.435` reproduces the coupled case to every digit printed.
`Mtg 1 H 4.435` injects four extra mmol of free H, which appears as **3.8×
too much H2** and a spuriously more reducing pe. So the design doc's formula
is right, and so is its reasoning for *this* channel.

### Why `Sg 1 H -1` is wrong

For sulfur the ΔH criterion *is* the right one, because H2Sg's hydrogens are
in the H mass balance on both sides. `Sg 1 H -1` makes the solid **absorb**
one H from the ocean where it should **release** one — a 2 mol H per mol S
error. Consequences at reducing conditions:

| | coupled `S 1 H 1` | `Sg 1 H -1` | `Sg 1 H 1` | `Sg 1` |
|---|---|---|---|---|
| pe | −4.334 | **+15.377** | −1.780 | +15.302 |
| m(H2) | 1.173e-09 | **0** | 5.484e-15 | 0 |
| sulfide speciation | H2S 4.43e-4 / HS⁻ 5.57e-4 | H2Sg 4.43e-4 / HSg⁻ 5.57e-4 | same | same |

All three Sg variants give **identical sulfide speciation** — the H term does
not change how Sg speciates. What it changes is the hydrogen and redox
budget, and `Sg 1 H -1` destroys the H2 reservoir and drives pe to +15.4.

**→ Use `Sg 1 H 1`.**

(Note the coupled and decoupled sulfur cases are not expected to match pe
exactly: coupled S must be *reduced* to sulfide by the system, consuming H2,
which is precisely the equilibration that decoupling exists to prevent. The
conservation requirement is on the H the solid releases, and `Sg 1 H 1` meets
it exactly.)

### Also settled
- **KINETICS and REACTION agree exactly.** Test A's variants 1/11 and 3/13
  gave identical totals, so the `-formula` semantics are the same in both.
  No disagreement to report.
- No variant errored. `Sg 1` runs fine — it is just not element-conserving.

---

## 2. Test B — `Hdg` does stop the carbon runaway

`test_B2_Hdg_carbon.pqi`, five simulations. Equilibrated values:

| sim | setup | pe | m(CH4) | fraction of C reduced | log10(CH4/CO2) |
|---|---|---|---|---|---|
| 1 | coupled H(0), no minerals | −9.204 | 1.250e-04 | 12.5% | +5.81 |
| 2 | **decoupled Hdg**, no minerals | +3.000 | **0** (exactly) | **0%** | −∞ |
| 3 | coupled H(0), real Fe assemblage | −9.100 | 1.084e-04 | 10.8% | +5.59 |
| 4 | **decoupled Hdg**, real Fe assemblage | **−8.259** | 4.143e-09 | **0.0004%** | +0.68 |
| 5 | decoupled Hdg + Ca + Calcite | −5.997 | 3.5e-29 | ~0 | −18.4 |

**Sim 4 is the one that matters** — it is the realistic case, and it answers
two things at once:

1. **Carbon reduction is suppressed by ~26,000×** (10.8% → 0.0004% of carbon
   converted to methane). Essentially all carbon stays oxidised. The
   hypothesis in `OPTION3_decision.md` §2 holds.
2. **pe is genuinely set by the mineral assemblage, not left unconstrained.**
   pe = −8.26 with magnetite and pyrrhotite both at SI = 0, versus −9.10
   coupled. So the system stays reducing and physically sensible. (Without
   minerals — sim 2 — pe is essentially arbitrary, +3.0, because nothing
   constrains it. Real ENTICES runs always have the mineral assemblage, so
   this is not a practical problem, but it is worth knowing that Hdg alone
   does not *set* a redox state.)

**Sim 5 confirms the design's main selling point.** With `Hdg` in play,
carbon still participates fully in carbonate equilibria: calcite precipitated
8.160e-04 mol to SI = 0.00, drawing total C from 1.00e-03 to 1.84e-04. This
is what would have been lost by decoupling carbon itself.

### The caveat, now quantified: decoupling H2 is not redox-neutral

`OPTION3_decision.md` §2 warned the approach "may prove too blunt —
decoupling H2 affects *every* redox couple." It does. Mineral mass transfers
in the same assemblage, coupled vs decoupled:

| phase | coupled (sim 3) | decoupled Hdg (sim 4) | change |
|---|---|---|---|
| Magnetite | +1.347e-05 | +3.945e-05 | **2.9× more precipitated** |
| Pyrrhotite | −3.041e-05 | −1.083e-04 | **3.6× more dissolved** |
| aqueous S | 1.304e-04 | 2.083e-04 | +60% |

So switching H2 to `Hdg` shifts the Fe/S mineral budget by a factor of ~3 in
even this small system. That is a real side effect, not a rounding artefact,
and it will propagate into ENTICES's secondary mineralogy. It does not
invalidate the approach — the carbon fix is worth far more than the sulfur
perturbation is worth worrying about — but **the coupled/decoupled pair can no
longer be read as "identical except for carbon speciation."** Any comparison
of the two runs should report the mineral assemblage alongside the carbon
result.

---

## 3. What this changes

**Adopt, tested:**
- `IOM_CH4` decoupled formula `Mtg 1 H 0.435` — confirmed exactly correct.
- `IOM_S` decoupled formula **`Sg 1 H 1`** — corrected from `Sg 1 H -1`.
- `IOM_N` stays coupled — unchanged; nothing here bears on it.
- `IOM_CO2` / `IOM_CHn` keep their coupled formulas; decoupling via `Hdg` is
  a solution/database-level switch, and it works.

**Still open:**
- Whether `Hdg` should be applied globally or only below some temperature
  (`PROGRESS.md` Part 5 open question) — unaffected by these tests.
- The mineral-budget side effect above needs a decision: accept and report,
  or restrict `Hdg` to the carbon-relevant cases only.
- Test B used a deliberately minimal assemblage (magnetite + pyrrhotite).
  Re-running sim 3/4 with ENTICES's actual secondary-phase list would firm up
  the magnitude of the mineral perturbation before it goes in a paper.

---

## 4. State reconciliation — several items in NEXT_STEPS.md are already done

`NEXT_STEPS.md` and `PROGRESS.md` Part 5 were written without sight of the
`iom-module` branch's current state. Reconciling:

| item | doc says | actual |
|---|---|---|
| NEXT_STEPS #4 / PROGRESS #3 — implement distributed Ea | "VS Code session, ~half day" | **Done**, commit `d21fb57`. Sub-pools for CO2 (7), CH4 (13), CHn (7); `logA = 15.301` on every row. Numbers re-derived independently, agreed with `PATCH_F_vitrimat2018.md` to within 0.2%. |
| NEXT_STEPS #4 — "refit IOM_N Ea to 244.8 kJ/mol" | 244.8 | **Superseded: 242.8 kJ/mol.** See below. |
| NEXT_STEPS #5 / PROGRESS #4 — re-run 350 C residue check | "~1 h, the payoff" | **Done.** H/C 0.703 → **0.679** (measured 0.63); O/C 0.104 → **0.0909** (measured 0.096). Improved in the predicted direction but **not fully resolved** — see below. |
| PROGRESS #5 — module support for decoupled mode | pending | **Done**, commit `3df3eed`, updated today with the tested formulas. |
| NEXT_STEPS housekeeping — "Part B assertion changes still uncommitted" | uncommitted | **Committed.** |
| NEXT_STEPS paths `claude/PATCH_F_*.md`, `claude/OPTION3_decision.md` | in `claude/` | No `claude/` directory exists; these files are at repo root. |

### The IOM_N Ea correction — 244.8 is wrong, use 242.8

Both `PATCH_F_vitrimat2018.md` and `NEXT_STEPS.md` specify 244.8 kJ/mol. That
value carries a basis bug found while re-deriving it: `m0_per_kg = 0.875` is on
the ENTICES-normalised bulk basis (bulk N = 2.142 mol/kg), but the fit target
used was Miller's raw Murchison measurement (real Murchison N = 1.7705 mol/kg).
Those differ by 2.142/1.7705 = **1.2098**, and computing
`target_fraction = raw_350C_value / m0_per_kg` silently divides one basis by the
other. Ea then absorbs a ~21% conversion factor that has nothing to do with
activation energy.

The dimensionally clean target is the **ratio of Miller's two raw
measurements**, 0.5676/0.7229 = 0.7852, giving **Ea_N = 242.8 kJ/mol**.
Verified: predicted release divided by the basis factor reproduces Miller's raw
350 C and 500 C values to within 0.06%.

This also means the "~20-22% overshoot at 500 C, because m0 is a fixed
exhaustion pool" explanation repeated across several documents is the right
*sign* but the wrong *mechanism* — it was this basis mismatch. IOM_CO2/CH4/CHn
do **not** share the bug: their shift is fit from the ratio of two raw
measurements (basis-independent), and the ENTICES scaling is applied afterwards
as a separate clean step.

### The residue check did not fully resolve — worth deciding what that means

`NEXT_STEPS.md` #5 set this up as a clean hypothesis test, so recording the
outcome honestly: adding the low-Ea tail and CHn moved predicted 350 C residue
H/C from 0.703 to **0.679** against a measured 0.63, i.e. the overshoot fell
from +11% to **+7.8%** rather than disappearing. O/C improved from +8% to
+5.3%. Direction right, magnitude only partly explained.

Per NEXT_STEPS' own framing ("if it does not → the +10% comes from somewhere
else ... and the diagnosis in PROGRESS.md Part 4 should be corrected rather
than repeated"), the missing-low-Ea-tail diagnosis is **partially** supported,
not confirmed. Remaining candidates are the ones already listed there: m0
basis, unreacted channel mass counted as solid, sample heterogeneity. Nothing
was tuned to improve the fit.
