# Redox budget diagnostic — organics-on vs organics-off, matched ENTICES case

Date: 2026-09-20
Closes: `PROGRESS.md` Part 5 item 1 (the last item from the original
`NEXT_STEPS.md` list). Run with PHREEQC 3.8.6 batch against
`Core11_idealgas_mod_v4.dat`.

**This was not actually blocked on Lucas's runs.** `phreeqc_inputs/` already
contains a real ENTICES-generated input at 15 wt% organics
(`enceladus_k1.00e-12_d600_t7.50e+04_circ1.pqi`). Its organics-off twin was
built by removing only the 9-line `REACTION 1` block (verified: the diff
between the two files is exactly that block; nothing else differs — same
SOLUTION, same KINETICS primary minerals, same EQUILIBRIUM_PHASES list).

**Caveat that matters for interpretation:** this file uses the OLD inline
organic treatment (elemental `REACTION` with aromatic/thiophene pseudo-species
— `[(aro)-O-(aro)]`, `[(6)(CB)(CB)S]` — added in one step, not the calibrated
distributed-Ea `iom_module.R` channels, which are not wired into ENTICES yet).
So this diagnoses **current production ENTICES's** behaviour, which is what
was asked for, but it is not a test of the new module specifically.

---

## Headline: H2 depletion is confirmed and large. Fe response is NOT what the simple hypothesis predicted.

| quantity | organics OFF | organics ON | ratio (on/off) |
|---|---|---|---|
| pH | 8.617 | 5.792 | — |
| pe | −7.334 | −3.101 | — (+4.2 units, more oxidising) |
| m(H2) | 2.364e-06 | 3.617e-09 | **0.00153** (653× depletion) |
| m(Fe2+) | 1.436e-06 | 5.662e-05 | **39.4×** (more, not less) |
| m(Fe3+) | 6.2e-27 | 4.5e-21 | negligible both |
| aqueous S | 3.555e-07 | 3.367e-08 | 0.095 (90% down) |

**H2 depletion (653×) and pe becoming more oxidising (+4.2 units) are exactly
what `PROGRESS.md` Part 2b predicts**: released organic carbon draws
electrons from the H2 reservoir to reduce toward CH4, leaving the system less
reducing overall even though the carbon that formed is itself highly reduced.
This part of the hypothesis holds.

**But the "electron demand shows up as mineral over-oxidation" half does
not hold in the simple form PROGRESS.md posed it.** Fe2+ *increases* 39×
with organics on — the opposite of what "Fe(II) minerals oxidised to supply
electrons" would produce (that would show Fe2+ falling and Fe3+ rising).
Kinetic primary minerals (Forsterite, Fayalite, Enstatite, Ferrosilite,
**Pyrrhotite**, Anorthite, Albite, **Magnetite**) show **no measurable
difference between the two runs at all** (largest change: Fayalite,
−0.018%). Whatever electron budget effect exists, it is not visible in the
kinetic mineral dissolution rates over this single reaction step.

## What is actually driving the Fe/S mineral shift: pH, not redox

| secondary phase | off | on | change |
|---|---|---|---|
| Pyrite | 1.387e-07 | 7.722e-06 | **+7.58e-06** (56× more) |
| Greenalite (Fe-serpentine) | 6.032e-06 | **0** | disappears entirely |
| Cronstedtite-7A (Fe-serpentine) | 3.246e-06 | **0** | disappears entirely |
| Boehmite (AlOOH) | 3.444e-07 | **0** | disappears entirely |
| Beidellite-Fe (Fe-smectite) | 0 | 1.833e-07 | appears |

The pH drop (8.62 → 5.79, driven by the organic REACTION's own acidity/redox
inputs) is large enough to explain this pattern without invoking a
redox-electron-donor role for Fe: Greenalite, Cronstedtite and Boehmite are
all phases whose stability drops sharply below neutral pH, and their
disappearance frees Fe that would otherwise be locked in secondary silicates
— a plausible, simpler explanation for the Fe2+ increase than "iron minerals
oxidised to donate electrons." Pyrite's large increase (tying up more S into
a sulfide mineral) is consistent with S availability rather than a
redox-electron story either.

**This is worth being precise about, because it changes the actionable
conclusion.** The organics-on run in *current* ENTICES is not simply "H2 goes
down, Fe(II) goes up to compensate" as a clean redox pair — it is a
**compound perturbation**: the organic addition pushes pH and redox at the
same time, and the Fe/mineral response looks driven primarily by the pH
shock, with the H2/pe shift as a separate, genuinely redox effect running in
parallel. Disentangling them would need a controlled run (e.g. buffered pH)
that is out of scope here, but the two should not be conflated into one
"electron demand" narrative.

## What this does NOT settle

- This is one matched pair, one timestep (`-steps 1`, a single ~75,000 year
  increment), from the OLD inline treatment. It shows the artefact **exists
  and is measurable** in production ENTICES output, which was the ask — it
  does not establish its magnitude across parameter space, nor say anything
  about the new distributed-Ea module (not wired in).
- Kinetic primary minerals showing essentially zero reaction progress in
  *either* run over one step is itself worth flagging to whoever owns
  ENTICES's kinetics: it suggests this particular case is far from
  equilibration at this timestep, which is a separate question from the
  organics comparison but affects how much weight to put on "no primary
  mineral difference" as a finding (there may be too little reaction happening
  in either case, in this one step, for a difference to show up yet).

## Reproduce

Scripts in `phreeqc_inputs/redox_budget/`: `org_on.pqi` (unmodified copy of
the repo's organics-on input), `org_off.pqi` (REACTION block removed,
diffed to confirm nothing else changed), `compare.R` (final-row comparison
across redox indicators, kinetic minerals, and all secondary phases with a
non-trivial difference).
