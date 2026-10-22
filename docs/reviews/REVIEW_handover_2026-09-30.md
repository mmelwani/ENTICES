# Review of HANDOVER_for_review_2026-09-30.md

Date: 2026-09-30
Reviewer: the claude.ai planning session (read the local repo on branch
`lucas-integration` @ `4966349` plus uncommitted working files; re-derived
numbers in Python; could **not** run PHREEQC).
Reads with: `claude/CR_IOM_COMPOSITION_RESOLVED_2026-09-20.md`,
`claude/IOM_INTO_PHREEQC_PLAN_2026-09-20.md` (both in the claude.ai Project, not
the repo — which is why the handover missed the first one).

---

## Verdict

The work is sound, and several of its self-corrections were right. Two things
change the handover's conclusions:

1. **§2 is wrong about where the composition comes from.** Decision #1 (CI vs
   Murchison) has already been settled.
2. **§7.3, the residue check, passes less cleanly than claimed.** Like-for-like,
   RO/C is +14% (1.2σ), not +4.6%.

I'd also change one recommendation: **make `IOM_S` coupled.** That removes the
sulfide-mineral gap (§5.4) entirely rather than patching it.

---

## 1. §2 — the bulk composition is Miller et al. 2026 SI Table 9 (CR IOM), not CI

Settled on 2026-09-20 with the PI: `iom_bulk_formula()` is SI Table 9
(C 728, H 46, O 159, N 30, S 37 mg/g) converted to mol/kg, and it matches to
1 part in 10⁴ on all five elements. Table 9 is the mean IOM of Alexander et al.
2007's seven CR chondrites, chosen because CR IOM is the least altered group.

The handover traces the formula to the Rmd, but the Rmd is a *derivative* of
Table 9, not the source. It carries two conversion slips:

    Rmd N 3.284 = 46 / 14.007     <- the H gram value divided by N's atomic mass
    Rmd O 9.941 = 159 / 15.9949   <- 16O isotopic mass, not O's atomic weight 15.999

Lucas's own code comment (`ENTICES_LF_20260917.R` ~line 1131) quotes the Table 9
gram values directly.

**Consequences**
- **Decision #1 goes away.** The composition stays as it is.
- The +7.6% starting H/C offset against Murchison is the real **CR-vs-CM
  material difference**. It is not a derivation artefact.
- The Rmd's N3.284 is a known slip in the Rmd only. Fixing it is the PI's call,
  since the Rmd derives its mineral endmembers from that target.
- **Edit §2 before the handover is pushed where Lucas will read it.**

## 2. The §7 re-derivations

### 7.1 `Ea_N = 242.8` — confirmed

With A = 2×10¹⁵ s⁻¹ and 48 h isothermal, solving f(350)/f(500) = 0.5676/0.7229
= 0.7852 gives **242.78 kJ/mol**. At that Ea:
- predicted raw 350 °C yield is 0.5671 (−0.09%);
- predicted 500 °C yield is 0.7232 (+0.05%).

The buggy target, 0.5676/0.875 = 0.6487, gives **244.78** exactly, so the
diagnosis of the old value is also confirmed.

Two points to add:
- **f(500) = 1.00000**, so the 500 °C match is by construction (as the handover
  says). In effect this is a **one-point fit at 350 °C**.
- m0 = 0.7229 × 1.2098 = 0.8746 assumes N release is complete by 500 °C/48 h.
  That makes 0.875 mol/kg, 41% of IOM N, a **floor** on the reactive N pool.
  This doesn't matter at Enceladus temperatures, but it belongs in the docstring.

### 7.2 Mtg/Sg asymmetry — confirmed, and there is only one criterion

The handover worries that Mtg was judged on aqueous-state match and Sg on ΔH.
Both are the same test: **electrons released to the system**, assuming C ends
as CH4 (−4), S as sulfide (−2) and H as H(+1).

    formula                 net e- to system
    C 1 H 4.435 (coupled)        +0.435
    Mtg 1 H 0.435                +0.435   ✓ match
    Mtg 1 H 4.435                +4.435   ✗ (=> Test A's 3.8x excess H2)
    S 1 H 1 (coupled)            −1.000
    Sg 1 H 1                     −1.000   ✓ match
    Sg 1 H −1                    −3.000   ✗ (=> pe +15.4, H2 destroyed)
    Sg 1                         −2.000   ✗ (=> pe +15.3)

The rules only look opposite because Mtg enters as neutral CH4-equivalent
(0 e⁻), while Sg enters as H2Sg (S −2, which needs 2 e⁻). Test A's pe results
line up with this table.

**Recommendation:** replace the "opposite rules" text in `OPTION3_decision.md`
and the module docstring with this single criterion. It is easier to apply to
the next pseudo-species.

### 7.3 Residue ratio — passes, but O is weaker than stated

RX/C is **not** invariant to the starting composition when a large fraction of
an element is released. At 350 °C/48 h the model releases **44–49% of the O**.
The handover's +0.1% / +4.6% compares a CR start with Murchison residue data.
Like-for-like, using a Murchison start on Table 5's own basis
(H/C 0.63/0.90 = 0.700, O/C 0.096/0.53 = 0.181) with the same releases:

    start                       RH/C (meas 0.90±0.12)     RO/C (meas 0.53±0.06)
    CR / Table 9 (handover)     0.901  (+0.1%)            0.554  (+4.6%, 0.4σ)
    Murchison (like-for-like)   0.888  (−1.4%)            0.604  (+14%, 1.2σ)

**Correct statement:** consistent within uncertainty. H is excellent. The model
under-releases O by roughly 15% relative to what was observed. "Resolved, it's
a starting-composition artefact" is true for H only.

Don't tune anything. Record the O result as the thing to watch. It sits with
the known point that water supplies part of the CO2 oxygen.

### 7.4 Redox budget diagnostic — agree it's contaminated; the pH-buffered control is the wrong next step

From `org_on.csv`: m(CH4) 9.84e-04, and m(CO2) + m(HCO3⁻) = 8.49e-04. The
dissolved organic carbon splits **about 1:1 between CH4 and CO2**. That is what
elemental C(0) injected by the old `REACTION` block would do if it
disproportionated (2C + 2H2O → CH4 + CO2). The pH drop from 8.6 to 5.8 then
follows from carbonic acid.

So "pH shock + H2 depletion" is probably **one causal chain**, not two
independent perturbations. A pH-buffered run would suppress the mechanism
rather than isolate it.

What's worth doing instead:
- **(a)** A carbon mass balance on `org_on`. The block adds 0.02809 mol C, but
  only ~8×10⁻⁴ mol is dissolved (in 0.43 kg water). Where did the rest go?
- **(b)** Re-run the diagnostic with the **new module**, coupled and `Hdg`, once
  a database runs. That is the real Step 1. The current result answers only
  "old output is contaminated" (yes).

Sulfur, which matters for §3 below: **sulfate is negligible in both runs.**
SO4²⁻ is 6.7×10⁻¹⁶ against 3.4×10⁻⁸ for H2S + HS⁻, even in the most oxidised
case (pe −3.1, pH 5.8).

### 7.5 Wiring "minimal diff" — confirmed by reading

- All three paths (main, bulk restart, tidal restart) iterate generically over
  `self$iom_pools`.
- Both restarts take `-m0` from the `<name>_mol` CSV columns, with a guard that
  stops if the columns are missing.
- The tidal-restart KINETICS splice is newline-safe (checked at lines 3414–3421).
- `iom_total_moles` is misnamed (it holds kg), but the units are correct.
- Cosmetic: the RATES `REM` still says "single-Ea pool".

**One latent gap.** `redox_mode` is never persisted. A restart in a fresh R
object calls `iom_default_config()`, which is all-coupled. A decoupled run
would therefore silently restart coupled, and the column names are identical in
both modes, so the guard can't catch it. This is harmless while Option 3 is off
(Decision #4), but it needs fixing before it's turned on:
- save the config next to the CSV (for example `saveRDS`);
- assert on restart that the saved config matches.

---

## 3. Recommendations on the §8 decisions

| # | decision | recommendation |
|---|---|---|
| 1 | CI vs Murchison composition | **Dissolved** — composition is SI Table 9 (CR). See §1. |
| 2 | Accept `Hdg` mineral perturbation? | Unchanged: accept and report. Note its size is still only a toy-test figure. B3 was Fe-limited, so it couldn't measure it. |
| 3 | `Hdg` global vs thresholded | Unchanged: **global**. A threshold adds a free parameter and a discontinuity inside the tidal cycle. |
| 4 | Switch Option 3 on? | **Not yet.** First: a runnable database, #5–#7 settled, and `redox_mode` persisted. Then run coupled and `Hdg` as a **bracket**, not one or the other. |
| 5 | Canonical DB; `CO2 = CO2` | Both files define CO2(aq) again, fully, later in SOLUTION_SPECIES (`HCO3- + H+ = CO2 + H2O`, with its own `-Vm`). If PHREEQC keeps the last definition, Lucas's active identity line is **dead code** and the two databases agree on CO2(aq). **Needs a 2-minute executed check** (below). Either way, **comment it out in the merged file.** If the order ever flipped, CO2(aq) would become a log_k = 0 identity species, which is badly wrong. |
| 6 | Decoupled-sulfide mineral sink | **Make `IOM_S` coupled, for the same reason `IOM_N` is.** See below. |
| 7 | Decoupled-methane clathrate | Add **one** parallel phase, `CH4_hydrate_Mtg` (`Mtg:6H2O = Mtg + 6 H2O`, same log_k and analytic terms). Clathration isn't a redox step, so this is physically exact. Under `Hdg`, coupled CH4 is ~0 (~10⁻⁴ %), so nearly all methane is Mtg and the split-saturation error is negligible. |

### Why coupled `IOM_S` (Decision #6)

Decoupling exists to block an equilibrium that low-T kinetics wouldn't reach.
For sulfur, that means sulfide → sulfate. But **coupled equilibrium never goes
there in any case we've run:**

| case | pe | outcome |
|---|---|---|
| Test A, coupled | −4.3 | all sulfide |
| B3 | −8.3 / −9.1 | reducing |
| Lucas's run | −11.7 | reducing |
| redox diagnostic, organics on (most oxidised case) | −3.1 | SO4²⁻ 7 orders of magnitude below sulfide |

So decoupled S gains nothing, and it costs the mineral sink:
- IOM sulfide can't reach Troilite, Pyrite or Pyrrhotite, and piles up as HSg⁻;
- parallel `_Sg` phases would split pyrite saturation between two pools that
  can't exchange. That's a poor split here, because rock sulfide (Troilite)
  dominates.

With `IOM_S` coupled, Lucas's HS⁻/HSg⁻ question doesn't come up: nothing in
production would carry `Sg`.

- **Condition to watch:** if the **ocean** box ever goes oxidising (for example
  from surface-delivered oxidants), coupled S could form sulfate that kinetics
  would block. Keep `Sg 1 H 1` available in the config for that case.
- **Revise `OPTION3_decision.md` §1** to decouple C only, via `Hdg` + Mtg.

---

## 4. Commit / push

**Yes — commit the B3 fix and the handover to `lucas-integration` and push, but
correct handover §2 first.** Lucas reads PR #2, and §2 as written would send
him after a CI-vs-Murchison decision that doesn't exist.

For the pyrite overstatement already in `4966349` and in the PR description:
- **Don't amend or force-push** a branch Lucas has pulled.
- Put the correction in the new commit's message.
- The PI edits the PR description on GitHub.

## 5. Tasks for the Claude Code session

1. **CO2(aq) check (executed).** Using PHREEQC 3.8.6:
   - Run an identical 1 mmol/kgw C solution at pH 7, 25 °C, at 1 atm and at
     66.7 atm, against each database.
   - Compare m(CO2), log K and `-Vm` for CO2 in the species printout.
   - Look for any "defined more than once" warning.
   - Pass: identical results.
2. Correct handover §2 (per §1 of this review), then commit and push as in §4.
3. Carbon mass balance on `org_on` (§7.4a). Record where the 0.028 mol C went.
4. Persist `redox_mode` / config and assert it on restart (§7.5).
5. If the PI accepts §3's #6 and #7:
   - default `IOM_S` to coupled even when `C = "decoupled"`;
   - add `CH4_hydrate_Mtg` to the merged database;
   - update `OPTION3_decision.md` and the module docstring (including §7.2's
     single criterion).
