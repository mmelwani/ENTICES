# Next steps — what to do, who does it, and why in this order

> **SUPERSEDED 2026-09-20 (same day).** Items 2-5 below are **DONE** — the two
> PHREEQC tests were run (`PHREEQC_TEST_RESULTS_2026-09-20.md`), distributed Ea
> and decoupled-mode support are committed, and the residue check has been
> re-run. Only item 1 (the redox budget diagnostic) remains untouched.
> **Use `PROGRESS.md` Part 5 for the current list.** This file is kept because
> its per-task rationale and pre-registered pass/fail framing are still the
> record of what was asked and why — notably §5's "if it does not resolve, the
> diagnosis should be corrected rather than repeated," which is what actually
> happened (+11% → +7.8%, only partially resolved).
>
> Also corrected since writing: `Ea_N` is **242.8** kJ/mol, not 244.8 (a
> composition-basis bug), and the `IOM_S` decoupled formula is **`Sg 1 H 1`**,
> not `Sg 1 H -1` (wrong by 2 mol H per mol S).

Date: 2026-09-20
Reads with: `PROGRESS.md` Part 5 (the canonical list), `claude/OPTION3_decision.md`
Companion files produced with this note: `test_A_Sg_hydrogen.pqi`,
`test_B_Hdg_carbon.pqi`

---

## The short version

Three things are ready to run right now, and they are cheap. Everything else
is blocked behind them or is bookkeeping.

| # | task | who | effort | blocks |
|---|---|---|---|---|
| 1 | Redox budget diagnostic on existing output | Lucas (has the runs) | ~1 h | nothing — but may invalidate results already in hand |
| 2 | Test B: does `Hdg` stop the carbon runaway | anyone with PHREEQC | ~30 min | the entire Option 3 carbon design |
| 3 | Test A: does `Sg 1` demand its own H | anyone with PHREEQC | ~30 min | `IOM_S` decoupled formula |
| 4 | Implement distributed Ea | VS Code session | ~half day | the residue re-check, and all production runs |
| 5 | Re-run Murchison 350 C residue check | VS Code session | ~1 h | nothing; it is the payoff of #4 |

Do 1–3 first. They are all small, and two of them can change the design.

---

## 1. Redox budget diagnostic — do this first because it is nearly free

**What.** Take an existing ENTICES run with `organic_wt_percent > 0` and the
matching run with organics off. At matched timesteps, compare:

- dissolved H2
- pe
- the Fe-bearing mineral assemblage (magnetite, Fe(OH)2, Fe(III) phases)
- pyrrhotite

**Why it is first.** It needs no new code, no new database, and no new runs —
the output may already exist. And it tests whether the ~24 mol e⁻/kg IOM
electron demand (`PROGRESS.md` Part 2b) has been quietly corrupting results
already produced. If it has, that is worth knowing before anything else is
built on top.

**The trap to avoid.** Do not check H2 alone. The electron demand can be
satisfied by oxidising Fe(II) minerals instead, in which case H2 looks
untouched and the mineralogy is wrong. Check both, plus pe.

**What counts as a finding.** Any systematic difference in pe or in Fe-phase
abundance between the organics-on and organics-off runs, beyond what the
added mass alone explains. Report the numbers; don't adjust anything yet.

---

## 2 & 3. The two Option 3 PHREEQC tests — inputs are written and attached

I have written both as runnable `.pqi` files rather than describing them,
since that removes the main friction. Each needs one edit before running:
**replace `/path/to/Core11_idealgas_mod_v4.dat` with the real path.** Both
require that database (the decoupled species `Hdg`, `Sg` are not in stock
`phreeqc.dat`).

### Test B first, not Test A

`test_B_Hdg_carbon.pqi` decides the whole carbon-decoupling design, so it has
the highest leverage per minute. Five cases:

1. coupled H2, oxidised carbon added — should reproduce the artefact
2. decoupled `Hdg`, otherwise identical — should retain oxidised carbon
3. and 4. the same pair **with a reducing Fe(II) mineral assemblage present**,
   because minerals are the other electron source named in Part 2b
5. a calcite check, confirming that carbon still participates in carbonate
   saturation when `Hdg` is used (this is the design's main advantage over
   adding a decoupled carbon species)

**The single number that decides it:** `log10(CH4/CO2)` in cases 2 and 4
versus cases 1 and 3. Reference point from this project's own earlier output:
at pe −11.5, C(−4)/C(+4) was 3.1e-16 / 7.9e-33, about 17 orders of magnitude
favouring methane. Case 1 should reproduce roughly that; case 2 should not.

**If it partly works, that is still informative.** Watch pe in case 2. If pe
stays very low and carbon still reduces, the hypothesis is wrong in a useful
way — the mineral assemblage, not H2, is setting pe — and the fallback is the
carbonate database edit with its calcite caveat accepted. Cases 3 and 4 exist
precisely to catch that.

### Then Test A

`test_A_Sg_hydrogen.pqi` tries three variants (`Sg 1 H -1`, `Sg 1 H 1`,
`Sg 1`) via KINETICS, plus two via REACTION as a cross-check. The deciding
diagnostic is whether H2Sg forms at the expected ~1e-3 in all three, or only
where H is supplied explicitly.

**A failure is a result.** If a variant errors out, record the exact message —
a charge-balance or element-not-found error on `Sg 1` would itself settle the
question. Note also that if KINETICS and REACTION disagree, that disagreement
is the finding and should be reported rather than reconciled by picking one.

---

## 4. Implement distributed Ea — the largest single correctness gain

Spec: `claude/PATCH_F_vitrimat2018.md` §5, with the CHn rows in
`claude/PATCH_F_addendum_N_and_CHn.md` §2.

Three parts:
- replace single `IOM_CO2` / `IOM_CH4` with sub-pools on the Vitrimat 2018
  distributions (7 and 13 rows respectively), plus `IOM_CHn` (7 rows)
- set `logA = 15.301` on **every** row, including IOM_N and IOM_S
- refit `IOM_N` Ea to 242.8 [CORRECTED — was 244.8, basis bug; see PROGRESS.md Part 4] kJ/mol

**Re-derive the numbers, do not copy mine.** They were computed in Python in a
session with no R. If your values differ from the tables, yours are probably
right — say so rather than matching. This has already caught two errors in this
project (an 18% Ea error and a ~1000x low-temperature error), both of which
looked fine in documentation.

**Sanity checks that should pass afterwards:**
- `iom_validate_config()` still passes; all residues positive
- combined C release is now 7.653 mol/kg (12.6% of IOM carbon), up from 6.08
- sub-pool m0 sums equal the parent channel totals
- a new assertion: sub-pool total release equals single-channel release in the
  high-temperature limit where both exhaust

---

## 5. Re-run the Murchison 350 C residue check — the payoff

This is the one genuinely independent validation point (residue was never used
in any fit). The previous single-Ea model overshot: H/C 0.703 vs measured 0.63
(+11%), O/C 0.104 vs 0.096 (+8%), both high and in the same direction. That was
attributed to the missing low-Ea tail.

**Step 4 is the direct test of that hypothesis.** Adding the low-Ea tail should
pull predicted release up and residue H/C down. CHn helps too, since it is
H-rich (H/C 1.8). So:

- if the overshoot largely resolves → the diagnosis was right, and the
  distributed model is validated against data it was not fitted to
- if it does not → the +10% comes from somewhere else (m0 basis, unreacted
  channel mass counted as solid, sample heterogeneity) and the diagnosis in
  `PROGRESS.md` Part 4 should be corrected rather than repeated

Either outcome is publishable-grade information. Do not tune anything to make
it come out right.

---

## What I would NOT do yet

- **Wire the module into ENTICES** (`INTEGRATION POINTS`). Steps 2–4 can all
  change the config shape; wiring now means rewiring later.
- **Add any new database species.** Test B may make that unnecessary. If it
  fails, the carbonate edit is the fallback — but decide after the test.
- **Syn-IOM 250 C residue comparison.** Worth doing, but it carries the
  syn-IOM-vs-Murchison confound and is a refinement, not a blocker.
- **Chase IOM_S calibration.** Miller measures no H2S. Nothing to do until a
  source appears; flag it in sulfur-dependent output instead.

---

## One piece of housekeeping worth not forgetting

`iom_module.R` and `INTEGRATION_NOTE.md` have drifted. The copy available to
the chat sessions predates the VS Code IOM_N recalibration, which is why the
patch documents are written as *instructions* rather than replacement files.
Whoever merges the `iom-module` branch should reconcile `INTEGRATION_NOTE.md`'s
"what's NOT done" list against `PROGRESS.md` Part 5 — the IOM_N item there is
now stale, and several others have moved.

Also still uncommitted, per the 2026-09-17 handover: the Part B assertion
changes on the `iom-module` branch. Check `git status` before assuming they are
saved anywhere.

---

## Open questions that need you, not a test

1. **Decoupling: global or temperature-dependent?** Current recommendation is
   a global per-element switch, letting the coupled/decoupled bracket carry the
   uncertainty rather than introducing an arbitrary temperature threshold. Worth
   a decision before step 5 of `PROGRESS.md` Part 5.
2. **CHn sensitivity.** Its m0 comes from Burnham's vitrinite value (2% of C)
   with no IOM constraint. For H-rich IOM it could plausibly be several times
   larger. Worth one sensitivity run once the distributed model is in.
3. **Miller Table 10's +8 kcal/mol for HC113 vs our 0.0.** Unreconciled — his
   is a residue fit at 10 C/Myr, ours a gas-yield fit at 48 h isothermal. Does
   not block anything (the Murchison fit is independent) but should be
   understood before publication.