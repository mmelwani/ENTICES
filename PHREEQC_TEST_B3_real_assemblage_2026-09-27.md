# Test B, re-run against the real database and the real secondary mineral list

Date: 2026-09-27
Supersedes the "realistic" cases (sims 3-4) in `PHREEQC_TEST_RESULTS_2026-09-20.md`
Test 2, which used a 2-mineral toy assemblage (Magnetite + Pyrrhotite) and
`Core11_idealgas_mod_v4.dat`. Prompted by the question: which database and
mineral list am I actually testing against, and are there other differences
between Mohit's database and Lucas's beyond the known fluorapatite/clathrate
split?

Files: `test_B3_real_assemblage.pqi`, output in
`phreeqc_test_output/test_B3_real_assemblage.out`.

---

## 0. Database diff, done properly this time

Ran a full diff of `Core11_idealgas_mod_v4.dat` (Mohit's) against
`coreclath_CH4coupled.dat` (Lucas's, confirmed 2026-09-24 as what he actually
used for his shared example run), not just a keyword search. Both are
8200+-line files derived from the same lineage (RATES sections are
functionally identical; SOLUTION_MASTER_SPECIES differ only by the `F`
element).

**Known differences, confirmed:**
- Mohit's has `Fluorapatite` and `Fluorite` (PHASES), plus the supporting F⁻
  aqueous species (CaF⁺, HF, HF2⁻, MgF⁺, MnF⁺, NaF, HPO3F⁻, H2PO3F, PO3F⁻²,
  SiF6⁻²) — all dated additions from 2025.
- Lucas's has three clathrate phases Mohit's lacks: `CH4_hydrate`,
  `CO2_hydrate`, `H2S_hydrate` (marked "#Added by LF").
- Both databases carry the full `Hdg`/`Oxg`/`Mtg`/`Ntg`/`Sg` decoupled-species
  set identically — confirmed again here, not assumed.

**One difference neither of you had flagged:** the `CO2 = CO2` self-identity
species definition (sets CO2(aq)'s molar volume via the PS01 reference) is
**commented out** in Mohit's copy and **active** in Lucas's. Since `log_k=0`
either way this doesn't change equilibrium constants, only the molar-volume
(`-Vm`) correction applied to CO2 at pressure — a small effect, but real, and
worth reconciling since ENTICES runs at tens to hundreds of atm. Everything
else in the 104-line raw diff was either encoding-artifact noise (a `°`
character rendering differently) or two rate-law comments that turned out to
apply to numerically identical formulas.

**A structural point worth flagging now, not discovered by running anything
— just by reading the clathrate definition:** `CH4_hydrate` is written as
`CH4:6H2O = CH4 + 6 H2O`, using **coupled** `CH4`. If Option 3's decoupled
carbon treatment (`Mtg`) is ever switched on in production, methane held as
`Mtg` would not be seen by this clathrate reaction at all, since `Mtg` and
`CH4` don't share a mass balance. Whether that matters depends on how close
the system actually gets to clathrate saturation (see below) — recording it
here so it doesn't get rediscovered later as a surprise.

## 1. What "the actual secondary mineral list" is

There are two, not one — extracted directly from Lucas's real output
(`phreeqc_inputs/from_lucas_20260924/enceladus_..._orb12.pqi.out`), not
reconstructed from the template source:

- **Porewater** (`EQUILIBRIUM_PHASES 1`, 150 °C): 81 phases, includes `Fe`,
  `Ni`, `Lawrencite` (metal/chloride phases) and `Quartz`; no clathrates.
- **Ocean** (`EQUILIBRIUM_PHASES 2`, 1.45 °C): 78 phases — same list minus
  those four, plus `Magnetite` and the three clathrates.

Test B's solutions are at 1.45 °C, matching the ocean box, so this re-run
uses the **ocean list** (78 phases) and Lucas's real pressure for that box
(66.7 atm, added to every SOLUTION block here — the original test_B/B2 files
omitted pressure entirely, which matters for clathrate saturation
specifically).

Note `Pyrrhotite` — used as one of only two minerals in the original toy
test — **is not in either real secondary list**. In production it's a
kinetic primary phase (fed via `Troilite`'s dissolution), not something that
precipitates from solution. The real secondary reduced-sulfur phase is
`Pyrite`, which turns out to matter (below).

## 2. Result: the core finding holds, and holds more strongly

| | sim 1 (coupled, no minerals) | sim 2 (Hdg, no minerals) | sim 3 (coupled, real 78-phase assemblage) | sim 4 (Hdg, real 78-phase assemblage) |
|---|---|---|---|---|
| pe | −9.204 | +3.000 | −9.136 | **−8.283** |
| m(CH4) | 1.250e-04 | 0 | 1.0749e-04 | **6.072e-10** |
| % of C reduced | 12.5% | 0% | 10.75% | **0.0000607%** |

**Suppression factor with the real assemblage: ~177,000×** (vs. ~26,000× in
the toy-assemblage test). **pe shift: +0.85 units** (vs. +0.84 before) —
essentially identical to the toy test's magnitude, which is itself a useful
cross-check: the earlier conclusion wasn't an artefact of an oversimplified
mineral list.

**What's actually doing the buffering changed, and that's informative.**
Magnetite never saturates in this run (`SI_Magnetite` stays around −6.6 to
−9.7 throughout) — with only 1e-5 mol Fe in solution and no Fe3+/Fe2+ mix
established, there isn't enough iron in this dilute system to reach
magnetite saturation. **Pyrite does saturate** (`SI_Pyrite = 0.000` at the
end of both sim 3 and sim 4), consuming the S(-2) and buffering the system
instead. So the specific mineral responsible for the pe shift depends on
what's actually in solution and doesn't have to be Magnetite — the toy
test's choice of Magnetite+Pyrrhotite got the right answer for the wrong
specific mechanism. Worth knowing before generalizing "Magnetite buffers pe"
as a rule.

## 3. Sim 5 (calcite) — unchanged, still confirms

Recomputed with the corrected pressure (66.7 atm rather than uncontrolled):
Calcite goes from SI = 1.33 (initial, supersaturated) to SI = 0.00 (final),
precipitating 8.16e-04 mol and drawing total C from 1e-3 to 1.9e-4 mol/kgw.
Carbon still fully participates in carbonate equilibria with `Hdg` in play —
same conclusion as the original test, now confirmed at the real pressure.

## 4. What this does NOT settle

**The clathrate/Mtg gap remains untested, not resolved.** `CH4_hydrate`,
`CO2_hydrate` stayed deeply undersaturated in every simulation here (SI from
−2.8 to −107) — this dilute test (1e-3 mol/kgw total C) never got close to
clathrate saturation regardless of coupled/decoupled mode, so it can't
confirm or refute whether the coupled-CH4-only clathrate reaction actually
matters in practice. That depends on whether real ENTICES ocean
concentrations of methane get anywhere near clathrate saturation at 66.7 atm
and 1.45 °C — a question about the full model's actual state, not something
a small isolated test like this can answer. Flagged, not resolved.

**Still open from before, unaffected by this re-run:** the `Hdg`
mineral-perturbation decision (accept the Fe/S budget shift, or restrict
`Hdg` to organics-active cases), and whether decoupling is global or
temperature-gated.

## 5. Reproduce

`test_B3_real_assemblage.pqi`, run against `coreclath_CH4coupled.dat`
(now in the repo root). The 78-phase list and the pressure value are both
pulled directly from the real production output rather than retyped from
the template, so if ENTICES's phase list changes, this file's list will
need re-extracting the same way rather than hand-edited.
