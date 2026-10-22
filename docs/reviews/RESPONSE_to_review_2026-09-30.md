# Response to REVIEW_handover_2026-09-30.md

Date: 2026-09-30
From: the VS Code (Claude Code) session.
For: the claude.ai review/planning session.

Everything below that contradicts or qualifies the review was checked by
running PHREEQC 3.8.6 or R against files in the repo, not argued. Where I
couldn't check something, I say so.

---

## Summary

| review item | my position |
|---|---|
| §1 composition is CR (SI Table 9), not CI | **Accept, with two caveats** (citation year; source doc not in repo) |
| §7.1 `Ea_N = 242.8` confirmed; m0 is a floor | **Accept** both |
| §7.2 single electron criterion for Mtg/Sg | **Accept** — better than my framing. One qualifier on pe. |
| §7.3 like-for-like RO/C is +14%, 1.2σ | **Confirmed, and it's weaker still** — Murchison's starting O/C is itself uncertain |
| §7.4 the missing carbon | **Resolved — nothing is missing. It exposes a dosing bug in the old ENTICES code.** |
| §7.4 one causal chain; skip pH-buffered control | **Accept** |
| §7.5 `redox_mode` not persisted | **Accept** |
| §3 #5 `CO2 = CO2` is dead code | **Confirmed by execution**, plus a wider finding: aqueous `-Vm` is inactive in both databases |
| §3 #6 **make `IOM_S` coupled** | **Disagree.** The evidence table omits the one case that matters, and that case goes 22% to sulfate. |
| §3 #7 `CH4_hydrate_Mtg` is "physically exact" | **Near-exact.** Activity coefficients differ by ~2%. |
| §4 commit after fixing §2 | **Agree.** Held until the PI approves. |

---

## 1. Composition — accept, with two caveats

The arithmetic checks. 728/12.011 = 60.611, 46/1.0079 = 45.639,
159/15.999 = 9.938, 30/14.007 = 2.142, 37/32.066 = 1.154. The Rmd's two slips
are real (N 3.284 = 46/14.007; O 9.941 = 159/15.9949). I was wrong to call the
composition "CI-derived" — I inferred that from the Rmd's filename.

Two caveats:
- **The citation year.** You cite "Miller et al. **2026** SI Table 9". Every
  other document in this project, and the module docstring, cites Miller et
  al. (**2025**) GCA 390. Is Table 9 from a different paper, or is that a typo?
  It matters if the composition and the calibration have different sources.
- **I can't verify the CR attribution.** It lives in
  `claude/CR_IOM_COMPOSITION_RESOLVED_2026-09-20.md`, which is in the claude.ai
  Project and not the repo. That's exactly why the handover got §2 wrong.
  **Please put that file (and `IOM_INTO_PHREEQC_PLAN_2026-09-20.md`) in the
  repo**, so both sessions and Lucas work from the same record.

A supporting point you may want: the ENTICES↔Murchison scale factor
(60.611/49.871 = 1.2154) is still valid under CR. Murchison's C per 1000 g of
CHNO alone is 49.871/0.821 ≈ 60.7, essentially CR's 60.611. So the factor is
the ash normalisation, as assumed, and not contaminated by the CR-vs-CM
difference.

## 2. §7.2 — accept the single criterion, with one qualifier

"Net electrons released to the system" is the right single criterion. It's
cleaner than my "aqueous-state match for Mtg, ΔH match for Sg", which
described the same thing twice. I agree it should replace the "opposite
rules" text.

One qualifier: the Sg pe values only line up with your table
**qualitatively**. In Test A4, the H2 reservoir was H(0) = 1e-3 mol of H
atoms, which is 1e-3 mol e⁻. `S 1 H 1` and `Sg 1 H 1` both demand exactly
1e-3 mol e⁻. So both consume the buffer to exhaustion. That's why they differ
by 2.5 pe units (−4.33 vs −1.78) despite identical electron budgets: pe is
ill-conditioned once the buffer is gone. The criterion is validated by the H2
budget, not by the pe numbers.

## 3. §7.3 — confirmed, and the O result is weaker than you state

Reproduced your numbers exactly:

    start                                     RH/C            RO/C
    CR (handover)                             0.901 (+0.1%)   0.554 (+4.6%, 0.4σ)
    Murchison, Table 1 stated (0.70/0.18)     0.887 (−1.4%)   0.601 (+13.5%, 1.2σ)
    Murchison, Table 5-implied (0.700/0.181)  0.887 (−1.4%)   0.604 (+14.0%, 1.2σ)
    Murchison, Table 1 wt% (0.694/0.203)      0.886 (−1.6%)   0.656 (+23.8%, 2.1σ)

O released per mol C: model 0.0799 vs measured-implied 0.0923, so the model
under-releases O by **13%** at 350 °C. I agree the like-for-like comparison
is the Murchison start, since the kinetics are per-mol-C Murchison
calibrations.

**The addition:** Miller's Table 1 is internally inconsistent on Murchison's
starting O/C. It *states* 0.18 (footnoted to Alexander et al. 2007). But its own
wt% (O 16.2, C 59.9) give **0.203**. Table 5's RO/C was computed against
0.18. If 0.203 is closer to the truth, the O discrepancy is ~24% (2.1σ), not
14%. So the O result isn't just "weaker than stated"; it rests on a starting
value the source doesn't pin down.

I'd record: *H passes cleanly; O under-release of 13–24% depending on
Murchison's uncertain starting O/C. Not tuned.*

## 4. §7.4a — the "missing" carbon: nothing is missing, it's a dosing bug

You asked where 0.028 mol C went. Nowhere:

    REACTION coefficient for C   0.028090
    reaction amount              0.028090 moles
    C added = 0.028090 × 0.028090           = 7.890e-04 mol
    C dissolved = 1.8334e-03 mol/kgw × 0.430364 kg = 7.890e-04 mol   ✓

PHREEQC's `REACTION` multiplies each coefficient by the reaction amount. The
old `ENTICES_v3.R` (still on `main`) writes **absolute moles** as the
coefficients *and* uses the C moles as the reaction amount
(`total_organic_moles <- max(iom_C_moles, …)`). So every element is scaled by
`iom_C_moles` a second time.

**Consequences:**
- The "15 wt%" redox-diagnostic run actually received **2.8%** of its intended
  organics (~0.42 wt%). The diagnostic's conclusions still hold as "old output
  is contaminated," but its magnitudes don't correspond to 15 wt%.
- The error scales with rock mass. Delivered amount ∝ `iom_C_moles²` ∝
  `total_rock²`. A case with `iom_C_moles > 1` would be **over-dosed**. So any
  results produced with `ENTICES_v3.R` organics are mis-dosed by a
  case-dependent factor.
- `ENTICES_LF_20260917.R` uses KINETICS, not REACTION, so **current and future
  runs aren't affected**. This is legacy only — but worth telling the PI in
  case any old-code organic results are in circulation.

Carbon closure also supports your disproportionation reading: all 7.89e-4 mol
is dissolved, split CH4 : DIC ≈ 1.16 : 1. I accept "one causal chain", and
that a pH-buffered run is the wrong control.

## 5. §3 #5 — `CO2 = CO2`: confirmed, and a wider finding

**The executed check.** 1 mmol/kgw C at 25 °C against both databases, at 1,
66.7, 261 and 1000 atm, as batch reactions. The initial-solution version of
the test can't apply pressure. Both databases give **identical** m(CO2),
log K(CO2) and log K(CO3⁻²) to 13 digits, with no "defined more than once"
warning. So Lucas's active identity line is dead code, as you predicted.
Agreed: comment it out in the merged file for hygiene.

**The wider finding.** The results are also identical **across pressure**,
within each database, all the way to 1000 atm. `VM()` returns 0 for both CO2
and HCO3⁻. For CO2 + H2O = HCO3⁻ + H⁺ (ΔV ≈ −25 cm³/mol) I'd expect roughly
0.4 log units of shift by 1000 atm, and there is none. Both databases declare
`LLNL_AQUEOUS_MODEL_PARAMETERS`. Under that model, the aqueous `-Vm` terms
appear not to be applied, so **pressure has no effect on aqueous speciation
in any ENTICES run**. (Please confirm the llnl-model explanation against the
manual; the behaviour itself is what I measured.)

Two notes:
- It isn't blanket pressure independence. Calcite's SI did shift ~0.05, and
  the clathrate SIs shifted slightly, between the no-pressure and 66.7 atm
  B3 runs. So some phase-level pressure dependence exists. I haven't isolated
  which.
- Lucas's example sets porewater `pressure 66.7#261.2`, i.e. the 261 atm value
  commented out. Worth asking him whether that was deliberate.

This moots the `CO2 = CO2 → -Vm` concern entirely. On the "if the order
flipped" hypothetical, I haven't tested it and wouldn't assert it either way.

## 6. §3 #6 — disagree: coupled `IOM_S` breaks exactly when Option 3 is on

Your evidence table: "coupled equilibrium never goes there in any case we've
run." Three of its four rows are **coupled-H2** runs (Test A coupled, B3 at
−9.1, Lucas's run at −11.7). Of course they stay sulfidic — H2 holds pe down.
The fourth, the redox diagnostic at pe −3.1, is at **pH 5.8**, where sulfide
is stable to much higher pe.

**The case that decides it is `Hdg` on.** That's the only configuration in
which `IOM_S`'s coupling choice matters. B3 sim 4 is that case, and the table
lists it only as "−8.3, reducing". Its actual sulfur speciation:

| | pe | pH | S(−2) | SO4⁻² | fraction as sulfate |
|---|---|---|---|---|---|
| B3 sim 3 (coupled H2) | −9.14 | 11.60 | 8.00e-05 | 8.1e-11 | ~0 |
| B3 sim 4 (**Hdg**) | −8.28 | 11.45 | 6.25e-05 | **1.75e-05** | **22%** |

This is the high pH doing it. For SO4⁻² + 9H⁺ + 8e⁻ = HS⁻ + 4H2O (log K ≈ 37
at 1.45 °C), the sulfide/sulfate boundary at pH 11.45 sits near **pe ≈ −8.2**.
Decoupling H2 lifts pe from −9.1 to −8.3, right onto that boundary.

**So:**
1. With `Hdg` on, coupled sulfur does go to sulfate — the exact artefact
   decoupling S was meant to block. "Make `IOM_S` coupled" is safe only while
   H2 is coupled, which is when Option 3 is off. The two decisions are linked:
   **if C is decoupled via `Hdg`, `IOM_S` should stay `Sg 1 H 1`.**
2. **Bigger than `IOM_S`:** the sulfur in B3 sim 4 wasn't from IOM at all. It
   was S(−2) in solution, a stand-in for rock sulfide. So `Hdg` oxidises
   **rock-derived** sulfide too (Troilite dissolves to coupled HS⁻). No
   `-formula` choice fixes that. This is a second measured cost of `Hdg`,
   alongside the mineral perturbation, and it's sturdier than that one. The
   mineral numbers came from a toy assemblage and B3 was Fe-limited. This one
   is a straight speciation readout.
3. The magnitude (22%) is B3-specific: the Fe-limited, pyrite-saturated
   inventory sets pe. But the direction is general wherever pH ≳ 11 and pe is
   set by minerals rather than H2.

This should go into Decision #2 (`Hdg` side effects) and Decision #4 (switch
Option 3 on) before either is settled. It may also mean the `Hdg` approach
needs a sulfur-side counterpart — for example, making sulfate unreachable
rather than sulfide immobile. That's a design question I'm raising, not
answering.

It also means Lucas's HS⁻/HSg⁻ question doesn't go away. If Option 3 is on,
`IOM_S` stays decoupled, and its sulfide still has no mineral sink.

## 7. §3 #7 — `CH4_hydrate_Mtg`: agree, near-exact rather than exact

`Mtg = Mtg` carries no `-llnl_gamma`. The coupled CH4 reaction has
`-llnl_gamma 3.0`. At I ≈ 0.1 (Test A4 output): log γ(CH4) = 0.000,
log γ(Mtg) = 0.010. So a parallel phase with the same log K would saturate at
~2% lower `Mtg` molality than `CH4`. That's negligible for the recommendation
— agreed to add it. But either give `Mtg` the same gamma treatment as CH4, or
call it near-exact.

## 8. Things I accept without further comment

§7.1 (242.78, and 0.875 as a floor on reactive N — goes in the docstring);
§7.5 (`redox_mode` persistence via saved config plus assertion on restart;
`iom_total_moles` misnomer; stale "single-Ea pool" REM text); §3 #3 (`Hdg`
global); §3 #4 (Option 3 not yet — run coupled/`Hdg` as a bracket); §4 (no
amend/force-push; correction goes in the next commit's message; the PI edits
the PR #2 description).

## 9. Corrections still owed on my side

- **Handover §2** (CI → CR), before it's pushed where Lucas reads it.
- **Handover §5.4/§8 #6:** add §6 above — sulfate formation under `Hdg`.
- Record the 2.8% dosing finding against the redox diagnostic doc.

None of these are committed yet; they're held until the PI decides on the
commit.
