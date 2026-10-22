# Getting the organics into PHREEQC accurately — what is left, in order

Date: 2026-09-20
Supersedes: `NEXT_STEPS.md` (items 2–5 done), `PROGRESS.md` Part 5 (partly done)
Assumes: `iom_bulk_formula()` = Miller et al. 2026 SI Table 9, correct as given
(provenance of its S = 37 mg/g still to be chased — see
`claude/CR_IOM_COMPOSITION_RESOLVED_2026-09-20.md` §3; it does not block anything
below)

---

## Where things actually stand

The module is in better shape than the planning docs suggest. Verified from
`iom_module.R` and `PHREEQC_TEST_RESULTS_2026-09-20.md`:

- **Composition:** correct, traced to source. ✓
- **Kinetic scheme:** distributed Ea, 29 rows, `logA = 15.301`. ✓
- **Ea values:** CO2/CH4 fit to both Miller yields exactly; `Ea_N` = 242.8
  (basis bug fixed); `IOM_S` uncalibrated and flagged. ✓/flagged
- **PHREEQC `-formula` semantics:** tested in PHREEQC. `Mtg 1 H 0.435` and
  `Sg 1 H 1` confirmed; `Hdg` suppresses carbon reduction 26,000×. ✓
- **Text generators:** five functions, 16 self-test assertions passing. ✓
- **Wired into ENTICES:** **no.** Not one line.
- **Run in ENTICES:** **no.** The module has never produced a PHREEQC input
  that ENTICES executed.

So the remaining work is not more calibration. It is **integration and
whole-system verification**, and one unmeasured physical effect that could
invalidate results either way.

My own element budget, computed from the current config (independent of the
module's own validator):

    element   bulk    released   residue   % released
    C        60.611    7.6569    52.9541     12.6%
    H        45.639   18.6568    26.9822     40.9%
    O         9.939    5.4673     4.4717     55.0%
    N         2.142    0.8750     1.2670     40.8%
    S         1.154    0.4620     0.6920     40.0%

All residues positive; C release 12.6%, matching the module's stated value.

---

## Step 1 — The redox budget diagnostic. Do this first.

This has been item #1 on every plan since 2026-09-20 and has never been done.
It is still first, and the reason is sharper than before.

**The number that makes it urgent.** The released carbon has a weighted mean
oxidation state of **−0.72**, not +4:

    channel   formula            released C oxidation state   m0 (mol/kg)
    IOM_CO2   C 1 H 1 O 1.3115        +1.62                    4.169
    IOM_CH4   C 1 H 4.435             −4.435                   2.276
    IOM_CHn   C 1 H 1.8               −1.80                    1.212
    ------------------------------------------------- weighted mean −0.72

In **coupled** mode, PHREEQC will equilibrate that carbon with the ambient
redox state. If the pore fluid is oxidising enough to push it toward CO2, the
system must supply **(4 − (−0.72)) × 7.657 = 36.1 mol e⁻ per kg IOM**, i.e.
~18 mol H2-equivalent per kg IOM. That is enormous. It is roughly 50% larger
than the ~24 mol e⁻/kg figure in `PROGRESS.md` Part 2b, because that figure
used the CO2 channel's +1.6 alone and did not weight in the two reduced
channels.

**Those electrons come from somewhere, and PHREEQC will find them silently.**
Per the PI's own correction: not from IOM-derived H2 specifically, but from the
whole geochemical system — oxidising Fe(II) in the mineral assemblage,
consuming serpentinization H2, shifting sulfide to sulfate. Test B already
measured this happening in a *two-mineral* toy system: switching H2 to `Hdg`
moved magnetite 2.9× and pyrrhotite 3.6×. In a real ENTICES run with a full
secondary assemblage and 36 mol e⁻/kg on the table, the perturbation could be
large.

**What to do.** Take an existing ENTICES run pair — `organic_wt_percent > 0`
and the matching zero-organics run — and at matched timesteps compare:

1. dissolved H2
2. pe
3. **the Fe-bearing mineral assemblage** (magnetite, Fe(OH)2, any Fe(III) phase)
4. pyrrhotite and total aqueous S
5. the C(−4)/C(+4) split

**The trap.** Do not check H2 alone. If the electron demand is met by oxidising
Fe(II) minerals, H2 looks untouched and the *mineralogy* is quietly wrong. That
is the failure mode that would survive into a paper.

**Why first:** needs no new code, no new database, possibly no new runs — the
output may already exist. And it tests whether results already in hand are
contaminated. Everything below is building on top of that answer.

**What counts as a finding:** any systematic difference in pe or Fe-phase
abundance beyond what the added mass alone explains. Report numbers; change
nothing yet.

---

## Step 2 — Decide coupled vs decoupled, and global vs temperature-thresholded

Two decisions, both the PI's, both blocking Step 3. Step 1's numbers should
inform them, which is why they come second.

**(a) Is `Hdg` on in production?** Test B settled that it works. It did not
settle whether to use it. The trade:

- *Coupled:* CO2–CH4 equilibrate at low T, which is thermodynamically correct
  but kinetically wrong — the PI's original complaint, that released carbon
  spectates as CH4 when it should stay CO2.
- *Decoupled (`Hdg`):* carbon stays oxidised (12.5% → 0.0004% reduced), calcite
  equilibria still work, pe still set by minerals (−8.26). But the Fe/S budget
  moves ~3×.

My recommendation: **run both as a bracket, not one as truth.** The pair
straddles the real answer (kinetic inhibition is neither complete nor absent
over Gyr), and the honest statement is a range. That makes the Fe/S
perturbation a reported model uncertainty rather than a hidden choice.

**(b) Global or temperature-thresholded `Hdg`?** `PROGRESS.md` Part 5's open
question. Recommendation unchanged: **global**. A temperature threshold
introduces a free parameter with no measurement behind it, and a discontinuity
in the middle of a tidal cycle. The coupled/decoupled bracket already carries
the uncertainty that a threshold would be pretending to resolve.

**(c) Accept or restrict the `Hdg` mineral perturbation?** Recommendation:
**accept and report**, with Step 1's numbers quantifying it in the real
assemblage. Restricting `Hdg` to "carbon-relevant cases" is not well-defined —
it is one solution-level switch, not a per-reaction one.

---

## Step 3 — Wire the module into ENTICES

Only after Step 2, because the config shape depends on those answers.

The integration points are documented at the foot of `iom_module.R`. Four
splices: RATES, KINETICS, USER_PUNCH headings, USER_PUNCH lines. Plus the
restart path.

**Three specific hazards, in decreasing order of how quietly they fail:**

1. **`USER_PUNCH` line-number collision.** The config is **29 rows**, not the 4
   it had when the integration note was written. `iom_punch_lines()` steps by
   10, so starting at 200 occupies **200–480**, not 200–230. The module's own
   note says "~290 lines of space (e.g. 200..3080)" — that is wrong arithmetic
   (29 rows × 10 = 290 *of span*, so 200..480). Either way: check ENTICES's
   highest existing PUNCH line number and leave real headroom. A collision
   silently overwrites a mineral PUNCH statement and corrupts a CSV column.

2. **29 KINETICS reactants may not integrate.** This is the one I would expect
   to bite. The previous config had 4; CVODE now has 29 additional stiff ODEs
   whose rate constants span ~11 orders of magnitude at low T (Ea from 47.6 to
   79.7 kcal/mol). Watch for non-convergence, and be ready to set `-cvode true`
   with a tighter tolerance or to prune negligible sub-pools. **Do not prune by
   eye** — prune only bins whose contribution over the full run is demonstrably
   below the CSV's precision, and record which were dropped.

3. **`iom_validate_config()` cannot validate a decoupled config.** Known and
   deliberate: it rejects `Mtg`/`Sg` as unknown elements. So in decoupled mode
   ENTICES gets **no element-balance check at all**. Do not teach the validator
   to accept pseudo-elements (the balance it reported would be meaningless).
   Instead: validate the *coupled* config as a proxy before switching formulas,
   and assert that the two configs differ only in the `formula` column.

4. **Restart path.** `-m0` must come from the recovered `<name>_mol` CSV
   columns, not recomputed from `iom_mass_kg` — otherwise a restarted run
   resets the IOM pools to full and silently re-releases everything already
   released. The module documents the `iom_mass_kg = 1` trick for this. This is
   a correctness bug waiting to happen in any multi-segment run.

---

## Step 4 — Whole-system verification: three tests with pre-registered outcomes

The module is unit-tested. ENTICES+module has never been tested. These are the
checks that catch integration errors, and each has a defined pass/fail *before
running* so the result cannot be rationalised afterwards.

**Test 1 — Mass conservation end to end.**
Sum the punched `<name>_mol` columns at the end of a run. Each must equal
`m0_per_kg × iom_mass_kg` minus what remains, and total released C must not
exceed 7.657 × `iom_mass_kg`. *Pass:* agreement to CSV precision. *Fail:* any
element over-drawn → a `-formula` or `-m0` splice is wrong.

**Test 2 — The zero-organics limit.**
Run `organic_wt_percent = 0`. *Pass:* byte-identical PHREEQC input to
pre-integration ENTICES, and identical output. This is the cheapest possible
test that the splices did not disturb anything, and it should be run first.

**Test 3 — The isothermal Miller reproduction, through ENTICES.**
Force ENTICES to 350 °C isothermal for 48 h with Murchison-equivalent IOM mass,
and check it reproduces `iom_selftest.R`'s yields. *Pass:* within a few percent.
*Fail:* the discrepancy is pure integration error, since the module alone gets
this right — which makes it the single most diagnostic test available.

I would add one more, which nothing in the current plan covers:

**Test 4 — Does the IOM change the ocean at all?**
The actual science question. Run the two-box model to steady state with and
without organics and compare ocean-side CO2, CH4, NH3/NH4⁺, H2S and pH. *This
has no pre-registered pass criterion* — it is the result, not a check. But it
should be run early, at low resolution, because if the answer is "no measurable
change at any plausible organic fraction," that reframes the project's central
question and is worth knowing before polishing the integration.

---

## Step 5 — Sensitivities, once the above works

In priority order, by how much they could move the answer:

1. **`IOM_CHn` m0.** Currently 2% of C, Burnham's *coal* value, with no IOM
   constraint at all. For H-rich CR IOM it could plausibly be several times
   larger. It is 1.212 of 7.657 mol/kg released C — 16% — so a 3× change moves
   total C release by ~30%. Highest-leverage unconstrained number in the config.
2. **`IOM_S` rate.** Uncalibrated (Miller measures no H2S). Bracket it; do not
   present sulfur output as constrained.
3. **Vitrimat distribution *width*.** Open item from the composition work: the
   Ea *centres* are on chondritic data, the *widths* are Burnham's vitrinite
   widths. Test by widening/narrowing the CO2 distribution at fixed mean and
   seeing how much low-T release moves. This is the one that could be orders of
   magnitude, given that single-Ea vs distributed was 3–5 orders.
4. **Organic mass fraction** itself.

---

## What I would explicitly NOT do

- **Do not add any new database species.** Test B made that unnecessary.
- **Do not re-fit anything to improve the residue agreement.** The +7.8% H/C
  overshoot is now explained as the real CR-vs-CM material difference
  (`CR_IOM_COMPOSITION_RESOLVED` §7), not an error. Tuning it away would be
  fitting to a known material mismatch.
- **Do not chase the HC113 +8 kcal/mol vs our 0.0 discrepancy** before
  publication-prep. Residue fit at 10 °C/Myr vs gas-yield fit at 48 h
  isothermal; it blocks nothing.
- **Do not renormalise `iom_bulk_formula()`.** It now matches a published table
  exactly; changing it breaks that agreement and every calibration on it.

---

## Ordering summary

    1. Redox budget diagnostic on existing output        (no new code; may invalidate)
    2. Decide: Hdg on/off, global/thresholded, accept perturbation   (PI)
    3. Wire module into ENTICES                          (4 splices + restart)
    4. Tests 2, 1, 3, 4 in that order                    (cheapest/most diagnostic first)
    5. Sensitivities: CHn m0, IOM_S rate, Ea width, organic fraction

Steps 1 and 2 are the only ones that can change the design. Everything after is
construction and verification.
