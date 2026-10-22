# Review of CONVERGENCE_DIAGNOSTIC_for_review_2026-10-06.md

Date: 2026-10-06
From: the claude.ai review session.
For: the PI, the VS Code (Claude Code) session and Lucas.
Read:
- the diagnostic;
- the PHREEQC v3 manual (GAS_PHASE, Example 7);
- `ENTICES_LF_20260917.R` at `rock-composition` @ `9a7ced5` (local checkout);
- Lucas's production `.pqi`;
- RECONSTRUCTION_STATUS.md.

I could not run PHREEQC. Arithmetic was re-derived by hand; nothing here is
new PHREEQC output.

---

## Verdict

The diagnostic is careful and its main conclusions hold.

- **§7.1, the gas-phase finding, is confirmed, and it's wider than stated.**
  Exsolved gas is also locked out of pore-to-ocean exchange. That bears
  directly on ocean H₂/CH₄ — the Cassini comparison.
- **§7.2 (feldspar rate) and §7.3 (water/rock) are physically reasonable.** The
  fast first-orbit transient is real chemistry. The model just resolves it
  inside a single 33 h step.
- **§7.4: I'd replace the per-mineral equilibrium switch with a one-time
  spin-up run.** That removes the pH-dependence problem rather than working
  around it, and it doubles as the reconstruction's still-open forward check
  against King et al.

---

## 1. §7.1 — gas phases at SI = 0: confirmed, and two consequences beyond t = 0

### What the manual says

Example 7 (manual p. 321):

> "Another way to model gas-liquid reactions in PHREEQC is to maintain a fixed
> partial pressure by using the EQUILIBRIUM_PHASES data block…
> **Conceptually, an infinite gas reservoir is assumed**… If the gas reservoir
> is finite and the total pressure is constant, as in gas bubbles… in
> groundwater, then a **fixed-pressure gas phase** should be used."

And for GAS_PHASE `-fixed_pressure`: "a gas bubble forms as soon as **the sum
of the partial pressures** of the component gases exceeds the specified
pressure."

So `H2(g) 0 0` in EQUILIBRIUM_PHASES is an infinite H₂ reservoir/sink at
1 atm, whatever the solution's pressure. The diagnostic's reading is right.

### The replacement

**Use a `GAS_PHASE -fixed_pressure` at the porewater pressure, not
`SI = log10(70)` per gas.** SI = 1.85 fixes the first problem but keeps the
second. Each gas could independently reach 70 atm, so the partial pressures
could sum to ~770 atm, and the reservoir is still infinite. The fixed-pressure
gas phase is PHREEQC's built-in model of exactly this situation. A bubble
forms only when Σpᵢ exceeds P, and its amount is finite and tracked. The
database is ideal-gas, so Henry's-law solubility applies (no Peng–Robinson
needed).

Two things to settle with it:

- **Which pressure.** Lucas's porewater has `pressure 66.7#261.2`, with the
  261 atm value commented out. The bubble criterion now depends on P directly.
  H₂ solubility at 160 °C is ~1×10⁻³ mol kg⁻¹ atm⁻¹, so a bubble forms at
  ~0.07 molal at 70 atm and ~0.26 molal at 261 atm. Use the pressure that
  matches the porewater's depth, and keep the porewater SOLUTION and the gas
  phase consistent.
- **Where the bubble goes** (see below).

### Wider consequence 1 — the gas is locked in the porewater

Every tidal path (lines 1855–1907, 2408–2455, 4413–4437) does:

    REACT step:  USE equilibrium_phases 1 ... SAVE equilibrium_phases 1
    MIX steps:   USE equilibrium_phases none

So whatever exsolves into EP1's gas "phases" is **saved with the porewater
assemblage and never exchanged with the ocean.**
- It re-dissolves only into later porewater, which holds porewater pH₂ at
  1 atm.
- The ocean sees H₂, CH₄ and CO₂ only through that 1 atm buffer.

In the T150 baseline, ~3 mol of H₂ sat in EP1 after t = 0. **Ocean H₂ and CH₄
— the quantities compared with Cassini — are therefore set by this artefact,
in every tidal_kinetic run.** This is the "possible tie-in" the diagnostic
flagged, and the code confirms it.

Physically, a free gas in a porous core is buoyant and migrates upward. Once
there is a fixed-pressure GAS_PHASE, the design needs a rule for it. Three
options:

- (a) Transfer the bubble's moles to the ocean each orbit, in R, as a REACTION
  on solution 2.
- (b) Vent it (count it as lost to the plume).
- (c) Leave it in place as a documented end-member.

MIX doesn't move gas phases, so any of these is an explicit choice.

### Wider consequence 2 — the ocean side has the same construction

Lucas's EP2 (the ocean list) carries the same gases at `0 0` and is also
`SAVE`d. The original two-box design had ocean-side exsolution standing for
**plume loss**. But:
- 1 atm isn't the vent condition (the plume source is near the triple point,
  ~0.006 atm);
- a saved, re-dissolvable gas phase isn't a loss.

If plume loss is the intent, the exsolved gas has to be **removed** each
orbit, not saved. Worth a decision alongside the porewater one.

### And the t = 0 corrosion

This one is partly physical. With magnetite available and a 70 atm bubble,
the metal should go to magnetite instead of hematite. It will still corrode
substantially at t = 0, because the Fe/magnetite equilibrium H₂ fugacity at
160 °C is far above 70 atm. That follows from treating metal as
`dissolve_only` equilibrium, which the reconstruction itself records as
unsourced (§8 item 11). Native Fe probably does corrode fast in 160 °C water,
so the fixes don't stop t = 0 corrosion. What they change is its **products
and side effects**:
- no Fe(III);
- far less water consumed;
- no lawrencite-HCl acid pulse, which arose only because Fe²⁺ was oxidised
  into hematite, leaving Cl⁻ with H⁺.

---

## 2. §7.2 — feldspar dissolving within an hour is plausible chemistry

### Arithmetic

Re-derived (PK04 albite base branch: log k₂₅ = −15.60, Ea 71 kJ/mol,
n = −0.572):
- at 160 °C the Arrhenius factor is +3.88 log units;
- at pH 9 the pH term is +5.15;
- so log k = **−6.57**, matching the diagnostic's −6.58.

With SSA 5 m²/g × 262.2 g/mol = 1,311 m²/mol, the specific rate is
**3.5×10⁻⁴ s⁻¹ and τ ≈ 2,900 s.** Confirmed.

### Is the SSA × mw × M scaling inflating it?

No. It is dimensionally just "surface area per mol × mol present". If
anything is high, it's the **value** 5 m²/g:
- that corresponds to ~0.46 µm spheres (d = 6/(ρ·SSA), ρ = 2.6 g/cm³);
- plausible for chondrite matrix, high for chondrule-hosted feldspar;
- lab-to-field rate discrepancies are often 1–3 orders of magnitude.

### Does it matter?

Even 1,000× slower gives τ ≈ 3×10⁶ s (about a month). That is still
instantaneous on any Enceladus timescale. So the transient is real chemistry.
The problem is only that it's resolved inside the first 33 h step. **This is a
numerics problem, not a rate-law problem.**

## 3. §7.3 — the water/rock setup is reasonable

- 0.157 kg of water per 0.85 kg of rock is W/R ≈ 0.18 by mass. That implies
  ~32% porosity (ρ_rock 2.6). It's within the range used for Enceladus' porous
  core.
- Starting from pure water at pH 7, pe 4 is arbitrary but harmless: the rock
  sets pH and pe within seconds.
- A violent t = 0 is what fresh anhydrous rock plus water at 160 °C should
  do.
- The one thing that was *not* physical, the acid pulse, is a gas-phase
  artefact (§1).

## 4. §7.4 — replace the per-mineral switch with a one-time spin-up

The diagnostic's suggestion (switch a primary to equilibrium when τ is much
smaller than the step) is sound, but it can't use pH, because pH is an output.
A cleaner design removes that problem:

**Spin-up.**
- Before the tidal cycles, run one batch alteration of the fresh rock at the
  porewater T, with the log-spaced sub-steps the diagnostic already built.
- Run it until the fast primaries are exhausted or rates fall below some
  fraction per orbit.
- `SAVE` the solution, equilibrium phases, kinetics and gas phase, and start
  orbit 1 from that state.

Why this over a switch:

- **Nothing to predict.** The spin-up computes pH, so there's no circular
  criterion and no fixed T/pH table.
- **Every phase stays kinetic.** Depleted primaries keep their tiny remaining
  M, so the restart machinery and IOM carry-over don't change.
- **It's the reconstruction's open forward check.** The spin-up's end state is
  exactly what RECONSTRUCTION_STATUS §8 item 8 asks for: hydrate the Scenario R
  assemblage and test it against King et al.
  - 81–84 vol% phyllosilicate;
  - 6–10 vol% magnetite;
  - 4–7 vol% sulfide.

  That check is still listed as fully open, with no input written. The spin-up
  gives it for free, once per porewater T.
- **It isolates the cost.** One slow run per T, then the tidal loop starts
  from a quiet state.

**Caveat: spin-up may still be slow.** 500 s of model time took ~10 min, even
sub-stepped. The ~40 iterations per solve and the `Overflow`/infeasible
messages point to the gas phases and the ~75 `0 0` phases making each Newton
solve hard. The 1 atm H₂ reservoir holding ~3 mol is an obvious suspect. So I
agree with the diagnostic's order: **fix §1 first, then time a spin-up.** If
it's still slow after the gas fix, the next suspect is the size of the
`0 0` list. Prune phases that never come within, say, SI −5 in the
sub-stepped CSVs.

**On guard (a), a minimum orbit step size:** I agree it goes the wrong way.

## 5. Tephroite rate law (PR #3)

Agreed on all three points. One refinement on `kneut`:
- The reconstruction (§4) says there is **no reliable neutral branch**, and
  that Casey's pH > 7 data shouldn't be trusted.
- So the summed form is wrong: it doubles the rate at pH 7, then holds a flat
  floor equal to the pH-7 acid rate at higher pH.
- If a floor at pH ≥ 7 is intended (so that the rate doesn't vanish at high
  pH), write it as the acid branch evaluated at `max(ACT("H+"), 1e-7)`, not
  as a sum.
- Either way, flag it as a modelling choice, since no source supports it.

Add the `SR > 1` guard, as for every other primary.

## 6. Decisions for the PI

1. **Porewater gas:** GAS_PHASE `-fixed_pressure`, at which pressure (66.7 or
   261 atm)?
2. **Fate of the porewater bubble:** transfer to the ocean, vent, or leave in
   place?
3. **Ocean-side gas:** keep exsolution as plume *loss*? If so, remove the gas
   each orbit instead of saving it, and pick the pressure.
4. **Spin-up instead of per-mineral switching.** If yes, it also serves as the
   King et al. forward check.
5. **Tephroite at pH ≥ 7:** acid branch only, or an explicit floor.

## 7. For Claude Code / Lucas, in order

1. Add Magnetite to EP1 in every generator. Check the other three for the
   same omission, which the diagnostic didn't check.
2. Implement decision 1 (fixed-pressure gas phase) in the tidal_kinetic
   porewater step. Re-run the t = 0 table and the sub-stepped first orbit,
   and report the wall-clock time.
3. If it's still slow, prune the `0 0` list by observed SI (§4 caveat).
4. Prototype the spin-up at 160 °C and compare its end state with King et al.
5. Tephroite fixes (§5).
6. Then audit which existing production outputs carry the 1 atm H₂ artefact.
   Any ocean H₂/CH₄ already quoted is suspect.
