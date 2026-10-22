# Primordial CI Chondrite Mineralogy — Status
Last updated: 2026-09-15 (supersedes 2026-09-14)

## Purpose
Reconstruct the pre-alteration mineral assemblage of CI chondrite dust as
input to PHREEQC water-rock interaction models. Deliverable is a phase
list with (a) mass/volume modes, (b) kinetic rate laws for reactants,
(c) thermodynamic entries only for phases that may precipitate.

Secondary phases and carbon are OUT OF SCOPE here: secondaries are the
PHREEQC output, and carbon is constrained separately via IOM.

## Reading rules
- **[V]** Verified: traced to a source, or arithmetic shown and checkable.
- **[D]** Decision: a modelling choice, with who/why.
- **[U]** Unverified: needs checking before use.
- **[X]** Known wrong / discarded. Kept so it is not reintroduced.

---

## 1. HEADLINE RESULT

**[D] 2026-09-13, final: Lodders et al. (2025) Table 4 is the primary
composition** (A-quality-rated Mg and Si, real 1-sigma uncertainties).
Alexander (2019) Table 1/6 is kept only as a bracketing comparison, not
recomputed further. Scenario C remains a comparison, not a fit — see §5.

**[D] 2026-09-14: P and Schreibersite added to the fit.** P is now in
`elements_R` and Schreibersite (Fe3P) in `phases_R`. See §4 for the
finalized reaction chemistry and §8 for what's still unsourced (Vm, rate
constant). This does not reopen §2's Fo/En question: Schreibersite
carries no Mg or Si.

**[V]** Scenario R modes, NNLS point-estimate fit on the primary
(Lodders) target `b_CI_Lodders` (`recompute_O()` then `drop_TiCo()` on
`b_bulk_CI_LoddersT4`, elements_R, weight="relative"), wt% of the
anhydrous reservoir:

| Phase | wt% | vol% (anhyd) | vol% (bulk basis) |
|---|---|---|---|
| Enstatite MgSiO3 | 22.84 | 26.53 | 11.93 |
| Forsterite Mg2SiO4 | 22.61 | 26.23 | 11.80 |
| Troilite FeS | 20.95 | 16.94 | 7.62 |
| Fe metal | 13.19 | 6.24 | 2.81 |
| Albite NaAlSi3O8 | 8.34 | 11.89 | 5.35 |
| Diopside CaMgSi2O6 | 5.94 | 6.76 | 3.04 |
| Schreibersite Fe3P | 0.93 | **0.47** | 0.21 |
| Anorthite CaAl2Si2O8 | 1.73 | 2.34 | 1.05 |
| Ni metal | 1.65 | 0.69 | 0.31 |
| Magnesiochromite | 0.71 | 0.60 | 0.27 |
| K-Feldspar | 0.57 | 0.83 | 0.37 |
| Tephroite Mn2SiO4 | 0.52 | 0.46 | 0.21 |
| Rhodonite | 0.00 | — | — |

**[V] RHO_CI_GRAIN sourced 2026-09-15** (§4, §8 item 4; was
`RHO_CI_BULK`, a physically backwards name — see §8 item 4 for why):
bulk-basis vol% shifted up uniformly by **+1.65%** across every phase
(2.42→2.46 g/cm3 = ×1.016529), verified column-by-column against the
previous 2.42-based table — max deviation from the uniform ratio is in
the 4th decimal place (rounding only), not a disproportionate shift in
any one phase. Superseded the pre-2026-09-15 bulk-basis column above.

**[V] Vm sourced 2026-09-15** (§4, §8 item 1a): Schreibersite's vol% is
no longer excluded. Every other phase's vol% shifted down by the same
uniform 0.9953 ratio (simple renormalisation to include a 12th
volume-contributing phase, verified column-by-column, not a
disproportionate shift in any one phase) — this table's vol% columns
supersede the 2026-09-14 version by that small, explicable amount.

Closure **1000.00 of 1000 mg/g** (was 998.54; the ~1.46 mg/g P shortfall
this table used to carry is gone now that P/Schreibersite are in the
fit). Max |resid_pct| across all 13 fit elements: **4.0×10⁻¹⁴%** —
floating-point noise, exact to machine precision. **This is a stronger
version of the same structural-overdetermination pattern flagged in §4:
the system remained exactly solvable after adding a 13th element and a
12th phase because Schreibersite is the only phase carrying P, so P
alone fixes its abundance, then Troilite (only S host) and Fe metal
absorb the rest of Fe — a near-triangular structure, not evidence the
extended phase list is uniquely correct.**

`to_volume()` now correctly excludes Schreibersite from the vol% columns
(warning: "No Vm for: Schreibersite") rather than propagating NA into
every other phase's vol% — see §6 for the bug this replaced.

**[V] Old (no P) vs new (P/Schreibersite added) side by side, wt%, both
on the Lodders composition.** Every phase except Iron is
**byte-identical** — schreibersite carries only Fe and P, so nothing
else in the fit can move. Only Iron drops, by exactly the Fe mass
Schreibersite now claims.

| Phase | Old (no P) | New (with P) | Δ |
|---|---|---|---|
| Forsterite | 22.609 | 22.609 | 0.000 |
| Enstatite | 22.845 | 22.845 | 0.000 |
| Diopside | 5.945 | 5.945 | 0.000 |
| Anorthite | 1.729 | 1.729 | 0.000 |
| Albite | 8.344 | 8.344 | 0.000 |
| K-Feldspar | 0.571 | 0.571 | 0.000 |
| Magnesiochromite | 0.713 | 0.713 | 0.000 |
| Tephroite | 0.525 | 0.525 | 0.000 |
| Troilite | 20.948 | 20.948 | 0.000 |
| Iron | 13.978 | 13.189 | **−0.789** |
| Nickel | 1.649 | 1.649 | 0.000 |
| Schreibersite | 0 (absent) | **0.935** | +0.935 |

**[D]** Drop Rhodonite from phases_R: returns exactly zero, all Mn goes
to tephroite.
**[D]** Fe and Ni(+ now Schreibersite) are fitted as separate native
phases because the fit cannot know Fe/Ni alloy or that some Fe leaves as
phosphide. Report Fe-Ni alloy as **14.84 wt%**, Ni/(Fe+Ni) = **11.11
wt%** (was 15.63 wt%, 10.55 wt% before P was added — the alloy fraction
drops and Ni-richness rises because Schreibersite draws only on Fe).
That straddles the kamacite/taenite boundary (~7.5 wt% Ni), so it is
likely two-phase. Cosmetic for PHREEQC, where both are
EQUILIBRIUM_PHASES dissolve_only.

### Derived results worth stating
**[V]** Feldspar is **An15Ab80Or5** — sodic oligoclase (molar, on the
Lodders fit: An 15.5%, Ab 79.4%, Or 5.1% — unchanged by the P addition,
since Ca/Na/K/Al are untouched by Schreibersite). Earlier drafts had
An80Ab20, wrong by a factor of five in An. PK04's Oligoclase rate law
covers An 10-30 directly.
**[V]** Troilite on the bulk basis is **7.50 vol%**, unchanged by the P
addition (S is untouched by Schreibersite) vs King et al. 4-7 vol%
sulfide. The gap is the right size for non-sulfide S (Alexander et al.
2021: IOM carries 2-7% of bulk S; plus sulfate and S0). Not a problem.
**[V]** Metal is **2.36 mmol Fe/g** (was 2.50 before P/Schreibersite);
full oxidation to magnetite (3Fe + 4H2O → Fe3O4 + 4H2) yields **~3.1
mmol H2 per gram of rock** (was ~3.3). Some of the Fe budget now leaves
as Schreibersite instead of native metal, so this number legitimately
drops — it is not a regression. Directly relevant to habitability, and a
sanity check — if PHREEQC does not produce it, something is wrong.
**[D]** Note in the discussion: 14.8 wt% Fe-Ni metal + 0.9 wt%
schreibersite, 21.0 wt% troilite, enstatite ~ forsterite, reduced
throughout — mineralogically this resembles an enstatite chondrite
despite CI bulk chemistry. State it rather than letting a reader
discover it.

**This is now the FINAL primordial composition for PHREEQC input.** The
table above, together with the phase roles/kinetics in §7, is closed: no
further NNLS re-runs, no further sigma-sourcing (§2, §8), pending only
the forward alteration check against King et al. (§8 item 8) and the two
still-open Schreibersite sub-items (Vm, rate constant — §8) that do not
block the mass-balance table above.

---

## 2. THE OLIVINE/PYROXENE SPLIT — handle with care

**[X] DO NOT CLAIM "enstatite exceeds forsterite" on a mass or volume
basis.** It was claimed in an earlier draft and withdrawn.

**[V]** Closed-form norm (mmol/g), anorthite terms cancel:
Si_rem = Si - 3Na - 3K - 2Ca - 0.5Mn
Mg_rem = Mg - Ca + 0.5(Al - Na - K) - 0.5Cr
Fo = Mg_rem - Si_rem
En = 2*Si_rem - Mg_rem

At b_CI: Si_rem = 3.8952, Mg_rem = 5.4631, Fo = 1.5680, En = 2.3272
mmol/g = 220.6 and 233.6 mg/g. Matches NNLS to 0.2 mg/g.

**[V]** Crossover points for r = Mg_rem/Si_rem:
- molar: r = 1.5. Actual r = 1.4025, margin **6.6%**.
- **mass/volume: r = 1.4164. Actual r = 1.4025, margin 0.99%.**

The mass crossover is the one that matters, because Fo and En are
differences (Fo ∝ r−1, En ∝ 2−r) and multiplying by 140.69 and 100.39
brings two ~225 mg/g quantities into near-cancellation.

**[V] The 1% margin is INSIDE the inter-compilation spread.** Lodders
(2003) molar Mg/Si = 1.0398 vs Alexander Table 1 = 1.0303. On Lodders'
composition r = 1.4162, within 0.01% of the crossover — forsterite and
enstatite come out essentially equal by mass.

**[V]** Monte Carlo, n = 2000, correlated errors (sigma_shared = 0.04),
O recomputed from perturbed cations. Alexander-based run:

| Phase | 2.5% | 50% | 97.5% |
|---|---|---|---|
| Forsterite | 15.08 | 21.98 | 28.89 |
| Enstatite | 15.97 | 23.41 | 30.82 |
| Diopside | 5.37 | 5.94 | 6.53 |
| Anorthite | 1.21 | 1.68 | 2.13 |
| Albite | 7.78 | 8.33 | 8.89 |
| Troilite | 19.57 | 21.58 | 23.90 |
| Iron | 11.34 | 13.47 | 15.58 |

**[V]** Diagnostics all pass: cor(Fo, En) = **−0.976**; closed form
matches NNLS to <0.2 mg/g; max residual 0.036%. The MC is correct.

**[V]** Lodders-based run (2026-09-12, after fixing the sigma_rel/O bugs,
§6), on `b_bulk_CI_LoddersT4` with its own Table 4 1-sigma column
(`sigma_abs_Lodders / b_bulk_CI_LoddersT4`, same `sigma_shared = 0.04`):

| Phase | 2.5% | 50% | 97.5% |
|---|---|---|---|
| Forsterite | 20.685 | 22.530 | 24.813 |
| Enstatite | 19.348 | 22.946 | 25.935 |
| Diopside | 4.984 | 5.923 | 6.916 |
| Anorthite | 0.766 | 1.751 | 2.684 |
| Albite | 6.935 | 8.301 | 9.951 |
| Troilite | 18.104 | 20.932 | 24.370 |
| Iron | 11.454 | 13.993 | 16.096 |

Fo/(Fo+En) = 0.4442 / **0.4955** / 0.5620 — a dead heat at the median.
Diagnostics: cor(Fo,En) = −0.972; closed form − NNLS = +0.082 mg/g (Fo),
−0.173 mg/g (En), both <0.2 mg/g; max |resid_pct| = 0.0355%. All match
the diagnostics above to rounding — the MC machinery is unchanged and
still correct.

### What to report
- **Fo + En ≈ 45.4 wt% of the reservoir is tightly constrained.**
- **The split between them, Fo/(Fo+En), sits at essentially 0.50** — both
  Alexander- and Lodders-based point estimates (0.4854, 0.4974) and both
  Monte Carlo medians (0.4955 on Lodders) agree on this — with a wide
  uncertainty band (2.5–97.5%: roughly 0.44–0.56). CI bulk chemistry sits
  almost exactly on the olivine/pyroxene mass balance crossover: on the
  primary (Lodders) composition, Mg_rem/Si_rem = **1.4139** vs 1.4164
  crossover, **0.18% away**; on the Alexander bracket, 1.4025, 0.99%
  away. **[V] Corrected 2026-09-14:** an earlier version of this bullet
  gave "1.4025 actual" with no composition label, which read as the
  (now-primary) Lodders value but was actually Alexander's — found while
  fixing an unrelated diagnostic-chunk bug in the montecarlo .Rmd (§6).
- **Do not report a P(En>Fo) or P(Fo>En) figure anywhere, in this file or
  downstream.** It is a probability statement about a coin flip and adds
  no information beyond "the split is ~50/50, uncertain." `p_En_gt_Fo`
  has been removed from `mc_fit()`/`mc_fit2()`'s returned summary (see
  the montecarlo .Rmd) so it cannot be picked up and quoted later by
  accident.
- Olivine and pyroxene are **comparable, ~22 wt% each** on either
  composition. This is the real result.
- Show both Alexander and Lodders values; they bracket the crossover.

**Deprioritized, 2026-09-13.** Sourcing sigma(Mg/Si) more precisely
(formerly §8 item 1) and investigating why the point probability landed
where it did on one particular run (formerly §8 item 11) are both
retired: the Fo/En split is reported as ~50/50 with wide uncertainty
regardless of the exact sigma, and further precision here does not
change the phase assemblage used downstream (§1).

---

## 3. THE BASIS FACTOR — verified, and load-bearing

**[V]** Table 6 CI ÷ Table 1 CI, across 14 elements:

Mg 1.4675, Al 1.4643, Si 1.4673, Ca 1.4709, Ti 1.4773, Cr 1.4695,
Mn 1.4635, K 1.4652, Na 1.4677, P 1.4721, S 1.4673, Fe 1.4676,
Ni 1.4717, Co 1.4644

Mean **1.4683 ± 0.0038** (0.26% RSD; scatter is Table 1 rounding — Al is
given to 2 s.f.). So `ANHYDROUS_TO_BULK = 0.6811`.

**[V] CONSEQUENCE: Table 6's CI column is Table 1's CI column times a
single scalar. It carries NO independent information.** Monte Carlo must
perturb Table 1 and rescale, never perturb Table 6 as a separate
measurement.

**[V]** Any comparison of model modes against King et al. (2015) MUST
apply this factor. Without it every phase is overstated by ~47%.

---

## 4. Other verified facts

### Oxygen balance
**[V]** All lithophiles oxidised, Fe/Ni/Co/P reduced: O = 18.317 mmol/g =
293.06 mg/g vs 293.2 tabulated (0.05%).
**[V] INTERPRETATION IS CIRCULAR.** Table 6's caption states the
reservoir is reduced and anhydrous by construction. This is a consistency
check on the phase list, NOT evidence about the real precursor's redox
state. Do not present it as such.
**[V]** Also: the Scenario R fit is square (12 elements, 11 active
phases). An exact non-negative solution was guaranteed to exist, so the
~0.02% residuals and the weighting-insensitivity (22.051 vs 22.091) are
**not validation**. Do not report them as such.

### Aluminium budget
**[V]** Albite hosting all Na needs 83.0 mg/g albite carrying 8.55 mg/g
Al = **69% of the Al budget**; orthoclase hosting all K needs 5.70 mg/g
carrying 0.55 mg/g = 4.5%. Alkalis take ~74% of Al. K is cheap; Na is the
binding constraint. Albite consumes 3 Si per Na = 17% of the Si budget.

### Charnoz et al. (2026) Nature 652, 925-
**[V]** Elements tracked: H He O Mg Si Fe Al Na K Ni Ca Cr S C. **No P,
no Cl** — cannot constrain either.
**[V]** SI Sect. 3.4 reaction list forms **K2Si2O5** ("K. silicate").
**[V]** Microcline IS in their thermodynamic data but does NOT appear in
model output — K-feldspar was available to form and did not. Positive
evidence, not absence. K therefore goes to the amorphous phase.
**[V]** Anorthite DOES appear in type C at both fast and slow nebular
reaction rates (user-verified in SI).
**[X]** An earlier draft claimed "no feldspar in Scenario C" by
generalising one main-text sentence about plagioclase. Wrong.
**[U]** Does albite appear in type C? Main text singles it out as blocked
by Si/Al availability. If absent, Na joins K in the amorphous phase.
**[V]** Charnoz et al. state primary amorphous silicates cannot be
modelled for want of thermodynamic data — a gap they acknowledge.

### Schreibersite reaction — finalized 2026-09-14

**[X→V] Correction to the original justification.** An earlier version
of this file argued Fe3P is "Fe(+I) + P(-III), not Fe(0) + P(0)," and
used that to derive 11 electrons. Pasek & Lauretta (2005), Astrobiology
5, 515-535 (now in Project Knowledge) use the Fe(0)/P(0) formalism
instead. Both give 11 electrons total — the formal split between Fe and
P is a bookkeeping choice, not a chemical fact, since Fe3P is neutral
either way. The electron count was right; the justification for it was
an unnecessary and, it turns out, non-matching detour. Drop it.

**[V] Independent confirmation of the electron count.** Pasek & Lauretta
propose (their Reaction 7):

    Fe3P + 7 H2O = Fe3O4 + H3PO3 + 5.5 H2

Balanced (Fe 3, P 1, O 7, H 14 both sides — checked here, not just taken
on trust). **5.5 H2 per Fe3P** matches the coefficient derived
independently in this project before their paper was available — a real
independent check, not a coincidence of rounding.

**[D] Final reaction form: O2-consuming, matching Core11's native-Fe
entry style.**

    Fe3P + 2.75 O2 + 4 H+ = 3 Fe+2 + HPO4-2 + 1.5 H2O

Balanced on Fe (3), P (1), O (5.5), H (4), charge (+4 both sides —
checked here). O2 convention matches Core11's native Fe entry; HPO4-2
matches its Fluorapatite entry. Equivalent under PHREEQC's O(0)/H(0)
redox coupling to the H2-explicit form:

    Fe3P + 4 H2O + 4 H+ = 3 Fe+2 + HPO4-2 + 5.5 H2

(also checked: Fe 3, P 1, O 4, H 12, charge +4 both sides) — same mass
transfer, no physical O2 or H2 reservoir implied by either; both are
electron-accounting devices, same as O2 in the database's own native-Fe
entry two lines away. Chose O2-consuming for consistency with that
existing entry, not because it's more "real."

**[D] O2/H2 species choice for the eventual PHREEQC input (not yet
written — flag for whoever writes it):** use the redox-COUPLED O2/H+
species, NOT Core11's redox-uncoupled `Oxg`/`Hdg` tracer gases. Using
`Hdg` would decouple schreibersite's H2 from the Fe-metal corrosion H2
already being tracked for the ~3.1 mmol H2/g habitability number (§1) —
those need to stay on the same books.

**[V] P(V)-only endpoint is justified, not just a database limitation.**
Core11 has no P(III) species (master species is `P HPO4-2`; confirmed by
direct grep — no phosphite, no hypophosphite anywhere in the file).
Pasek & Lauretta's Fig. 4: phosphonate oxidation is rapid in the
presence of Fe metal, negligible with magnetite alone or in controls.
This assemblage carries ~14.8 wt% Fe-Ni metal (§1), so the reduced P
intermediates the database can't represent would likely be catalytically
oxidized to orthophosphate in this system anyway. The P(V) simplification
is therefore reasonably justified here, not merely forced by the
database — worth stating both reasons in the methods.

**[U] 2026-09-15 correction: the paragraph above is not fully
reconciled and its "reasonably justified" wording should not be trusted
without the note below.** Astrobio2005.pdf (the full paper) was read
directly this session, not just excerpts. Pasek & Lauretta test 5
solution types; in Mg/Ca-CHLORIDE (saline) solutions — with or without
added organics — **phosphonate is 100% of detected P**; orthophosphate,
hypophosphate, and pyrophosphate all fall below detection, and the solid
residue also contains whitlockite, (Ca,Fe,Mg)3(PO4)2, as a secondary
precipitate. This is the closest analog in the paper to THIS project's
own assumed fluid chemistry (§8 item 3: Cl in the initial fluid, ~520
mol ppm, following Neveu et al. 2015) — and it points toward a
**persistently reduced corrosion product (phosphonate), not P(V)**,
under conditions resembling this model's own. This has NOT been
reconciled with the Fig. 4 finding cited above (Fe-metal catalyzes
phosphonate oxidation): the two may not be in tension (different
experimental conditions; unclear whether the Fe-metal catalysis
experiment was itself run in saline solution) — but the paragraph above
currently cites Fig. 4 alone without checking it against this more
directly fluid-relevant result from the same paper. **[D] Practical
note, not yet actioned:** if phosphonate is the more fluid-realistic
corrosion product for this model's assumed Cl-bearing fluid, the P(V)
approximation UNDERSTATES persistence of reduced P rather than
overstating it — the opposite direction from what "database can only do
P(V)" implies to a reader. Whether this matters depends on how much
P-speciation detail the eventual PHREEQC run needs — a modelling-scope
decision, not something more literature search alone resolves.

**[X] SUPERSEDED, same day — the paragraph originally here claimed no
Fe3P sample mass exists anywhere in the paper. That was wrong.** It was
based on the H2-quantification passage alone (250 mL sealed container,
DI water, Ar-flushed, RGA pumping 35 bar gas/day, H2 via ideal gas law
from a 1 mbar partial-pressure rise at 293 K over 1 day → 4 +/- 2e-4 mol
H2/day), which indeed states no mass. It did not check Table 2, which
does. Kept here, struck, so the error and its correction are both on
the record rather than the wrong version silently vanishing.

**[V] Table 2 sample mass correctly matched to the H2-release
measurement.** A further correction, to an assumption an earlier session
made when it DID look at Table 2: the H2-quantification passage says
the experiment was DI water flushed with **Ar**, not air. Table 2 has
two DI-water rows (1.0473 g, Air; 1.2197 g, Argon) — **1.2197 g, pH 7,
20C, Argon is the correct match**, not the Air row assumed previously.

Bulk fractional dissolution rate, derived from the matched row (checked
independently, not just re-stated):
- moles Fe3P = 1.2197 g / 198.509 g/mol = 6.144e-3 mol
- moles Fe3P consumed/day = (4e-4 mol H2/day) / 5.5 = 7.273e-5 mol/day
- fractional rate = 7.273e-5 / 6.144e-3 = **1.184%/day**
- propagating the stated H2 uncertainty (+/- 2e-4, ~50% relative) gives
  a wide range, **0.59% to 1.78%/day** — carry this through any
  downstream use.

**[U] Still cannot be converted to mol m-2 s-1.** Table 2 gives sample
mass, fluid, atmosphere, pH, and temperature for all 11 rows — **no
grind size, mesh fraction, or BET surface area for any row** (checked
directly). A surface-area-normalized rate requires an assumed grain
size not present anywhere in this source. Do not estimate one and
present it as sourced.

**[U] The 1.18%/day figure is a single-day linear extrapolation from one
gas-pressure measurement, not a measured rate over an extended period**
— Table 2 doesn't show whether other rows/durations exist for repeat
sampling of the same condition over multiple days. If the paper
elsewhere reports H2 evolution over more than one day for this
condition, that would give a real time series; worth checking Fig. 3
(referenced in the text as the H2 evolution plot) directly if a more
complete rate curve is wanted — not done this session.

**[D] Recommendation, not yet actioned:** if a PHREEQC rate law needs
Schreibersite dissolution at some numeric rate rather than remaining
excluded, use the bulk fractional rate above directly against the
Schreibersite mass in the model (not a surface-area-normalized form) as
a simple first-order approximation, explicitly flagged wherever used as:
(a) derived from one lab measurement at pH 7, 20C, in dilute Ar-purged
water — not the model's actual fluid composition, temperature, or W:R;
(b) a single-day linear extrapolation; (c) not mechanistically tied to
grain size, so it will not respond correctly to any surface-area scaling
elsewhere in the model. A bounded, honestly-caveated placeholder, not a
first-principles rate law.

Item (1b), §8: **not fully open, not fully resolved** — a usable bulk
(mass-normalized) rate with a wide, honest uncertainty band exists; a
mechanistic (surface-area-normalized) rate does not and cannot be
derived from this source. Note also: the paper reports corrosion as
substantially slower under Ar than air elsewhere; the matched H2-release
row is the Ar condition, so this project's ~1.18%/day figure is
presumably on the slower (Ar) branch rather than the faster Air-row one
— **not independently confirmed by comparing the two Table 2 rows'
rates directly**, since only the Ar row has a matched H2 measurement to
derive a rate from at all.

**[X] Schreibersite Vm/density — was not sourced.** A value of
7.100 g/cm3 appeared in an earlier discarded (v0) version of this
project with no citation. Never reused (§6, "Entire v0 approach").
`phase_lib` in the main .Rmd carried Vm = NA for Schreibersite on
purpose from 2026-09-14 — `to_volume()` excluded it correctly (bug fixed
that day, see §6) rather than propagating NA into every other phase's
vol%.

**[V] Vm sourced 2026-09-15**, see §8 item 1a for the full derivation.
`to_volume()` no longer excludes Schreibersite; verified all other
phases' vol% are unchanged *in the sense that matters* — every other
phase's vol% dropped by the exact same uniform ratio (0.995316 =
1 − Schreibersite's new 0.468% share), i.e. simple renormalisation to
include a 12th volume-contributing phase, not NA corruption or a
disproportionate shift in any one phase. This re-confirms the
2026-09-14 NA-propagation bug fix is still working correctly — see §6
for the full before/after check.

### Tephroite kinetics — 2026-09-15, replaces the Fayalite proxy

**[X→V] The Fayalite proxy for Tephroite (§7) is corrected.** Direction
was flagged as wrong two sessions ago (Pincus et al. 2026); this session
sources and cross-checks a real, independent replacement rather than
leaving a qualitative correction in place.

**[V] Primary source: Pincus et al. (2026) ACS Earth Space Chem. 10,
1174-1184, Table S1**, citing Casey et al. (1993) GCA 57, 785-793 as the
underlying data: log kH+ (measured, 25 C) = −5.4 mol m-2 s-1, log A
(Arrhenius pre-exponential) = 7.6 s-1, mean Mn-O bond distance = 2.21 Å.
SI methods state PK04-style entries in this table are referenced to
pH 2 (log kH+ = log AH+ − (Ea/RT)log(e) + n·log{H+}). Cross-check on
Fayalite in the same table (log kH+ = −6.8, sourced to PK04): evaluating
Core11's own Fayalite acid branch at pH 2 (not pH 4) gives
10^(−4.80)·[10^-2]^1.0 = 10^−6.8 — exact match. This resolves a
provisional `[U]` flagged two sessions ago (a naive pH-4 evaluation gave
10^−8.8 and looked like a >1000x discrepancy vs. Pincus; that was a
reference-pH mismatch, not a real inconsistency).

**[V] Independently re-derived and cross-checked against Casey (1993)'s
own raw data** (`Casey1993_Appendix_tephroite_dissolution_rates.csv`,
user-supplied clean data — NOT the OCR-corrupted PDF appendix table used
in earlier sessions, which should not be used for numeric work).
Regression of log(rate) vs pH at 25 C, synthetic tephroite only,
bin-averaged replicates, pH 1.78-5.2: slope (reaction order, n) = −0.47
(matches Casey's own stated 0.47 ± 0.06 at 25 C exactly); intercept
(pH=0, mol/cm2/s) = −8.60 → mol/m2/s = −4.60; evaluated at pH 2
(mol/m2/s) = −5.54. −5.54 vs. Pincus's −5.4: agreement to 0.14 log
units, from two independent derivations off the same underlying
dataset. A real cross-check, not a restatement.

**[U] Temperature dependence NOT independently confirmed.** Casey report
the acid-branch order n increasing with T: 0.47±0.06 (25C), ~0.60 (35C),
0.64±0.06 (45C) — from their own careful matched-replicate-set
regression (5 pH sets: 2.5, 3.5, 4.2, 5.15, 6.2). A quick re-fit this
session using simple pH-bin averages across all rows at each
temperature (not Casey's matched-set method) gave slopes clustering near
−0.42 to −0.47 at all three temperatures — did NOT reproduce the stated
T-trend. Casey's paper explicitly warns the slope is "sensitive to the
number and pH of the samples included"; their careful method is trusted
over this session's rough one. **Do not treat the flat-slope re-fit as
evidence the T-dependence is smaller than Casey reports** — more likely
an artifact of a less careful method.

**[V] High-pH unreliability confirmed IN THE RAW DATA, not just
asserted.** At pH 8.0-8.1, the two experiments in the CSV show pH-stat
vs. batch-method disagreement up to 1.8x on the same run (3.4 vs 5.5
×10⁻¹³; 6.24 vs 10 ×10⁻¹³ mol/cm2/s) — a concrete instance of Casey's
stated "uncertainties as high as 80% at near-neutral/basic conditions."
Do not fit or trust a kneut/kbase branch from pH>7 data in this dataset.

**[V] Run-to-run scatter at fixed nominal conditions is real and
extent-correlated**, not just noise: three replicates at pH 2.5-2.51,
45C give 10.9, 7.4, 12.8 (×10⁻¹⁰ mol/cm2/s) at 26%, 17%, 27% mass
consumed respectively — rate tracks reaction progress, consistent with
Casey's leached-layer mechanism discussion. State this as an uncertainty
on the rate constant itself (~1.7x run-to-run at "identical" conditions),
not just on the regression.

**[D] Recommended RATES entry, acid branch only** (no reliable neutral
branch — see high-pH note above): sourced n=0.47, pH-2-referenced log
kH+ = −5.4 (Pincus/Casey). **Ea is UNRESOLVED** — do not assign a
numeric Ea from either this session's rough T-refit or by borrowing
Fayalite's Ea (94.4 kJ/mol) as an unlabeled placeholder. Casey's own
concrete quantification of what a temperature-independent approximation
simplifies away: a 10 K rise multiplies the rate by ~2.2x at pH 2,
~1.5x at pH 3-4, only ~1.1x at pH>5.

**[V] Critical scope caveat, applies to the entire dataset:** all Casey
(1993) experiments were run rigorously anoxic by design — N2-degassed,
O2 scrubbed to <6e-6 mol/L via Fe filings/hot Cu, CO2 stripped with
NaOH. The paper states Mn2+ oxidation to Mn3+/Mn4+ "severely retards"
dissolution. This rate law is valid ONLY under reducing conditions
consistent with Scenario R's reduced-precursor assumption (§1) — it does
NOT transfer to any oxidizing scenario without correction.

**[V] Ancillary, not currently used:** BET surface area of the
experimental powder (25-75 μm sieved, synthetic and natural) = 0.11
m2/g. Available for a future bulk-consumption-rate cross-check.

**[V] Implemented and hand-checked 2026-09-15, `Core11_idealgas_mod_v4.dat`
RATES block.** One correction made against the source material above
while transcribing it into PHREEQC's `kacid = 10^(log k0) * ACT("H+")^n`
form: **Pincus's tabulated −5.4 is the rate already evaluated at pH 2**
(their own stated SI convention), not the pre-exponential `log k0` that
slot needs — using −5.4 directly as `log k0` and then also multiplying
by `ACT("H+")^0.47` double-counts the pH dependence. The same relation
the source material itself uses to validate the Fayalite entry
(Core11's pre-exponential −4.80, at pH 2 gives −4.80 + 1.0×(−2) = −6.8,
matching Pincus's tabulated fayalite value) inverts to give the correct
Tephroite pre-exponential: log k0 = log kH+ + n×2 = −5.4 + 0.47×2 =
**−4.46**. Hand-checked: `10^(-4.46) * ACT("H+")^0.47` at pH 2, 25 C
(TK=298.15, so the Ea/T term = exp(0) = 1 regardless of the Ea=0
placeholder) = 10^(−4.46)×(10⁻²)^0.47 = 10^(−4.46−0.94) = **10^(−5.4)**,
exact — confirms the transcription, not just the arithmetic. The RATES
block's own comments (lines 7-12) document this derivation in place so
a future reader working from the .dat file alone sees it, not only this
file. `Ea = 0` is an explicit, heavily-commented **null placeholder**
(zero temperature correction) — not a sourced value, not Fayalite's Ea
borrowed unlabeled. No PHASES change needed: Tephroite's thermodynamic
entry (Vm = 47.52, Webmineral) already existed and is untouched.

### Magnesiochromite kinetics — literature check, 2026-09-15: no usable
### rate law found, and the absence is itself informative

**[V] Checked three papers, none provide a usable PHREEQC rate law,
for three distinct reasons — this outcome supports rather than
contradicts the existing "(none)" proxy (§7):**

Razavi et al. (2026) J. Sustain. Metall. 12, 3009-3026: industrial
mechanical-activation/acid-leaching study, 100-400C, concentrated
H2SO4 (50-300 kg/t). No pH-mechanism dependence reported (acid DOSE
varied, not pH). Reported Ea=18.78 kJ/mol is for DIFFUSION THROUGH A
PRODUCT LAYER (shrinking-core model), explicitly not surface-reaction
controlled — incompatible with the ACT("H+")^n surface-mechanism form
used throughout this database. Rate constants in min-1 (shrinking-core
fractional conversion), not mol m-2 s-1. Not usable at any level.

Ao et al. (2024) J. Environ. Manage. 349, 119475: chromite oxidative
dissolution via Mn(II/IV)-catalyzed Cr(III)->Cr(VI), pH 5/7/9, 25C,
anoxic/oxic/anoxic-oxic, 24-2350 h batch. WRONG REDOX REGIME for
Scenario R (reduced/anoxic) — this paper's entire subject is oxidative
Cr mobilization requiring Mn(III/IV) oxide formation, which itself
needs oxidizing conditions. No rate constants (A/Ea/n) reported at any
condition — concentration-vs-time curves only. Not usable.

Schindler et al. (2018) Am. Mineral. 103, 271-283: natural chromitite
(chromite + clinochlore) dissolution, acidic/neutral/alkaline, 25C,
4-week batch. KEY FINDING: chromite itself is REFRACTORY within the
rock and does not control Cr release — dissolution is dominated by the
associated clinochlore (Mg-Al phyllosilicate) matrix, and chromite is
released as intact NANOPARTICLES during clinochlore alteration rather
than dissolving congruently. No chromite-specific rate law is reported
because the paper's own conclusion is that this isn't the
rate-controlling phase. Directly and independently supports treating
Magnesiochromite as refractory/very-slow relative to the rest of the
Scenario R assemblage, consistent with the existing "(none)" proxy —
this is NOT the same situation as Tephroite (above), where a proxy was
actively wrong in direction and a real replacement existed to find.

**[D] No change to §7's Magnesiochromite kinetics.** "(none)" stands,
now with a literature basis for why: three independent studies across
industrial, oxidative-environmental, and natural-weathering contexts
all fail to produce a usable pure-chromite aqueous rate law, and one
(Schindler) explains why — it isn't the reactive phase. At 0.71 wt% of
the assemblage (§1), treating it as effectively inert on the timescale
of the model (or requiring an arbitrarily slow bound if a nonzero rate
is needed for completeness) is a defensible simplification with
citable support, not an unaddressed gap.

### Native Fe, Ni kinetics — checked 2026-09-15, no source found

**[U] No independently sourced rate law exists for aqueous Fe/Ni metal
corrosion in this project.** Searched `Core11_idealgas_mod_v4.dat`
directly (no native-metal RATES block exists) and PK04 (Palandri &
Kharaka 2004 does not cover metallic corrosion at all — it's a
silicate/carbonate/sulfide compilation). No paper currently in this
project addresses aqueous Fe/Ni metal corrosion kinetics.

Current model treatment (EQUILIBRIUM_PHASES, dissolve_only, §7) is
mechanistically plausible but NOT independently sourced: metallic
corrosion is generally electrochemical (coupled anodic Fe → Fe2+ +
oxidation / cathodic H+ or H2O reduction) and typically fast relative
to silicate network-breaking reactions at comparable surface area — a
standard qualitative argument in corrosion science, not a number
derived from any source in this project.

**Two caveats not yet checked:**
1. Whether `dissolve_only` on the metal composes correctly with the
   model's existing EQUILIBRIUM_PHASES secondary Fe minerals (magnetite,
   sulfides) rather than trapping Fe with nowhere to go — needs a
   PHREEQC input review, out of scope for this reconstruction-only
   thread.
2. Passivation (protective oxide/hydroxide surface layers) can slow
   metal corrosion substantially in some conditions and is NOT accounted
   for by the "always fast" assumption — plausibly unimportant in a
   strongly reducing, O2-free system (nothing to build a passivating
   oxide from), but this is reasoning, not a checked fact for this
   specific composition and fluid.

If a genuine kinetic rate is ever needed rather than equilibrium
`dissolve_only`, corrosion-science literature on native Fe-Ni meteorite
weathering (terrestrial or Antarctic finds) or steel/iron corrosion rate
studies would be the place to look — neither has been searched yet. Not
done this session.

### Glass / amorphous kinetics
**[V]** Gislason & Oelkers (2003): A_A = 10^-5.6 mol Si cm-2(geo) s-1,
E_A = 25.5 kJ/mol; roughness 92 = 23000/250 cm2/g, for THEIR ground
powder only.
**[V]** Mazer & Walther (1994): Ea = 22.7 ± 3.7 kcal/mol = 94.99 kJ/mol.
**[V]** Wolff-Boenisch et al. (2006) GCA 70, 858-870, pH 4, 25 C,
log(r/n_Si), mol Si m-2 s-1:
  geo glass −6.91 − 6.77(Si:O) r2=0.71; geo mineral −5.15 − 15.43(Si:O)
  r2=0.47; BET glass −7.12 − 10.04(Si:O) r2=0.53; BET mineral −7.91 −
  10.26(Si:O) r2=0.28.
**[V]** Glasses and minerals of **Si:O ≤ 0.28 dissolve at similar rates**
(95% CIs overlap). The ≥1.6 order-of-magnitude glass enhancement applies
only to Si-rich (rhyolitic, Si:O ~0.40) compositions.
**[V]** Forsterite, anorthite and wollastonite are **excluded** from the
mineral regression — they do not dissolve by Si-O breaking. Do not
extrapolate that line to forsterite Si:O.
**[V]** Gislason & Oelkers: basaltic glass rates at 5<pH<9, 0-30 C are
NOT significantly larger than enstatite, forsterite or intermediate
feldspar rates.

### Molar volumes (all sourced)
Core11 `-Vm`: forsterite 43.79, enstatite 31.276, diopside 66.09,
anorthite 100.79, albite 100.25, K-feldspar 108.87, magnetite 44.524,
magnesiochromite 43.564, lizardite 107.31, greenalite 115, troilite
19.07, Fe 7.092, Ni 6.588, spinel 39.71, akermanite 92.81, rutile 18.82.
Webmineral: tephroite 47.52, rhodonite 35.87.
Derived from Handbook of Mineralogy density: gehlenite 274.197/3.038 =
**90.26**; fluorapatite 504.30/3.175 = **158.8**; schreibersite (Fe3P
endmember) 198.509/7.44 = **26.68** (Fasiska & Zweil 1967, via Handbook
of Mineralogy — §8 item 1a).
**[U]** Fayalite 46.39 — still not confirmed against Core11.

---

## 5. Scenario C — use published modes, do not fit

**[D]** Decided 2026-09-11: take type C phase proportions **directly from
Charnoz's figures**. Do not run the NNLS fit.

**[D]** Stated position for the paper: doubt that Charnoz's model outputs
represent the primordial CC composition — closed system, solar
composition, no amorphous silicates, and a phase list that cannot
reproduce non-solar bulk chemistry. Say this explicitly; readers will ask
why the two scenarios are not weighted equally.

**[X]** Why the fit failed: 15 phases, 9 elements, underdetermined by six,
and exactly six phases returned zero — fayalite, magnetite, iron,
lizardite, spinel, gehlenite. Four of those are phases type C is DEFINED
by having. Greenalite at 30 wt% absorbed all non-troilite Fe. Root cause:
O is not a fit element, so nothing penalises an over-oxidised solution.

**[V]** The one useful output — the O deficit decomposes cleanly:
Fe0 → Fe2+ in greenalite +38.7; structural OH (2.91 wt% H2O) +25.8;
Na,K oxides not required −2.7; **total 61.8 mg/g** (reported 61.67).
Interpretable, but for an assemblage that is not Charnoz's.

---

## 6. Errors found (do not reintroduce)

**[X]** *Core11 fayalite kinetics.* Acid mechanism was missing its
ACT("H+")^n. Corrected form (user-verified against PK04 Table 20):

kacid = 10^(-4.80) exp(-94.4e3/8.314(1/TK-1/298.15)) * ACT("H+")^1.0
kneut = 10^(-12.80)exp(-94.4e3/8.314(1/TK-1/298.15))

**[U]** PK04 Table 20 reportedly gives NO neutral mechanism for fayalite,
yet kneut exists with an Ea identical to the acid branch. Trace it.
Lower priority now: fayalite is absent from Scenario R.

**[X]** *Glass script* (`Amorphous_glass_silicate_dissolution.Rmd`):
Branch B is geometric, Branch A is anchored on Mazer with no roughness
conversion — if Mazer is BET the branches differ by ~92x. The handover
claims normalisations were unified; the code does not do it. Anchor
`k_linear_ref <- 10^(-13.48)` at 85 C could not be located in Mazer &
Walther. Units are mol cm-2 s-1 where PK04/Core11 use mol m-2 s-1 (1e4).
Citation "Fajian et al." appears as both 2026 and 2006, and may be a
given name rather than a surname.

**[X]** *Assistant errors in this project:*
- "GEMS are ~100 nm" — **fabricated, no source**. The conclusion built on
  it ("grain size dominates by 3 orders of magnitude") is void.
- "Amorphous silicates dissolve orders of magnitude faster" —
  contradicted by Gislason & Oelkers and Wolff-Boenisch at the relevant
  compositions.
- "No feldspar in Charnoz type C" — over-generalised from one sentence.
- "No microcline in Charnoz's phase list" — wrong, it is in SI Table 2.
- P(En > Fo) ≈ 0.94 predicted — computed the margin in MOLES while the
  code compared WT%. The MC was right; the prediction was wrong.
- Predicted the correlation fix would move P to ~0.94 while having
  written the independence caveat into the code comments myself.
- Troilite Vm 18.20 — wrong; 19.07 correct.
- CIPW ordering rules imported from igneous petrology without
  justification for nebular condensates.
- O-balance presented as evidence when it is circular.
- **Entire v0 approach**: density-ratio "mass conservation" that cancels
  algebraically, invented hydration expansion factors (1.35, 1.30),
  invented magnetite partitioning (50/30/20, then 35/20/30/15), invented
  `mag_factor 0.8` and `anorthite_correction 0.8`, 70:30
  serpentine:saponite falsely attributed to Tomeoka & Buseck, circular P
  mass balance. **Do not mine v0 for numbers.**

### 2026-09-12 session: bugs found and fixed in the live .Rmd files

All six items below were reported by the user as present in the live
files; each was independently re-checked in the actual code before being
touched, per the "run it, don't recall it" rule. Two (#3, #5) turned out
to already be fixed and were left alone. Full diff context: [Primordial_CI_chondrite_mineralogy.Rmd](Primordial_CI_chondrite_mineralogy.Rmd),
[Primordial_CI_chondrite_mineralogy_montecarlo.Rmd](Primordial_CI_chondrite_mineralogy_montecarlo.Rmd).

**[X→V] #1 sigma_rel unit bug — CONFIRMED and FIXED**, in the montecarlo
file. Lodders et al. (2025) Table 4 ABSOLUTE 1-sigma values (mg/g, e.g.
Si = 4.00) were assigned directly into `sigma_rel`, which `mc_fit()` /
`mc_fit2()` use as a lognormal `sd` — a FRACTION. Reproduced the exact
failure before fixing: **P(En > Fo) = 0.0659** (handoff's "garbage
0.066" — matched almost to the last digit) with `NNLS hit maxit = 500`
warnings and Forsterite/Enstatite collapsing to 0 in most draws. Fix:
kept the absolute values in `sigma_abs_Lodders`, computed
`sigma_rel_Lodders <- sigma_abs_Lodders / b_bulk_CI_LoddersT4`. Resulting
fractions (all matched the handoff's expected ranges): Si 3.75%, Mg
3.14%, Fe 3.24%, Al 5.67%, Ca 5.95%, S 8.88%, Na 10.32%, O 1.72%.

**[X→V] #2 Oxygen input bug — CONFIRMED and FIXED**, in the montecarlo
file. `b_bulk_CI_LoddersT4["O"] <- 465.7` is Lodders' TOTAL bulk CI
oxygen (hydrous: water, organics, sulfate, ferric phases) — far too much
for the reduced-anhydrous reservoir. `recompute_O()` (deriving O from the
perturbed cations, same formula as the main .Rmd's "Check 1") already
existed but was not consistently reached: a stray, never-executed
fragment (`b_bulk_CI["O"] <- NA; tgt <- drop_TiCo(recompute_O(bp))`)
referencing an undefined `bp` sat between the two composition blocks.
Run as originally written, top to bottom, **the file crashes outright**:
`Error in recompute_O(bp) : could not find function "recompute_O"`
(recompute_O is defined further down, in the next chunk). Fix: deleted
the stray fragment; both `mc_fit()` and `mc_fit2()` now call
`recompute_O()` on the perturbed cations before `drop_TiCo()`, never
reading O from any table.

**[V] #3 dead check_inputs() — CHECKED, ALREADY FIXED, no change made.**
`check_inputs()` is a top-level function (not nested inside
`fit_scenario()`) and its call is already the first line of
`fit_scenario()` (main .Rmd). Grepped for a second definition; none
found. The handoff's description does not match the file's current
state — noted rather than acted on.

**[X→V] #4 no-op `* 1.4683` scaling — CONFIRMED and FIXED**, in the
montecarlo file (`drop_TiCo(bp * 1.4683)`, two call sites). Verified
algebraically (any constant scalar on `drop_TiCo`'s input cancels in its
closing `b*1000/sum(b)` renormalisation) AND numerically
(`max(abs(drop_TiCo(x) - drop_TiCo(x*1.4683)))` = 5.7e-14 on a random
perturbed vector — floating-point noise, not a real difference). Removed
both instances; `tgt <- drop_TiCo(bp)` now, no change to any reported
number.

**[V] #5 stale duplicate phases_R/elements_R/phases_C/elements_C —
CHECKED, NOT FOUND in the live .Rmd.** Grepped the whole project; each of
the four is defined exactly once in `Primordial_CI_chondrite_mineralogy.Rmd`.
Older duplicate-looking copies exist only in RStudio's own cache
(`.Rproj.user/…`), which is not part of the executed script. No change
made.

**[X→V] #6 leftover Ti/Rutile (CIPW chunk) — CONFIRMED and FIXED.** The
first code chunk of the main .Rmd was an entire superseded
single-scenario CIPW-ordering calculation (hand-typed densities,
`Rutile`, `Orthoclase`) predating the NNLS two-scenario approach. Besides
reintroducing `Ti`/`Rutile` (explicitly dropped elsewhere via
`drop_TiCo()`), it used CIPW allocation ordering — already flagged above
in this same section as an assistant error not to reintroduce — and
defined its own global `AW` (16 elements, no H) that would silently
shadow the second chunk's `AW` (17 elements, with H) if chunks were run
out of order or individually. Removed in full; replaced with a short
markdown note recording why (this file, in place). No downstream
variable depended on it; the NNLS fit results are unchanged (verified by
re-running — see below).

**New, not on the original list:** the stray `recompute_O(bp)` /
undefined-`bp` fragment under #2 is itself a distinct bug (a crash, not a
unit error) — flagged here since it wasn't named in the handoff and would
otherwise go unrecorded. Also: after the #1/#2 fixes, the Scenario R MC
re-run on the Lodders composition gave P(En>Fo) = 0.554, inside the band
the handoff said not to expect — **superseded 2026-09-13: this
investigation is deprioritized (§2), and P(En>Fo) is no longer reported
or computed at all** (`p_En_gt_Fo` removed from the MC functions'
returned summary). This paragraph is kept only as the historical record
of what the bug fix produced, not as an open question.

**Full re-run after all fixes, `Primordial_CI_chondrite_mineralogy.Rmd`
then `_montecarlo.Rmd`, both end to end:** Scenario R wt%, closure, and
residuals are **byte-for-byte unchanged** from the values already in §1
(Forsterite 22.051, Enstatite 23.381, closure 998.5/1000, max
|resid_pct| 0.0355%) — expected, since bug 6 only removed dead code that
nothing downstream read. Scenario C and the O-deficit decomposition
(61.67 mg/g) are likewise unchanged. The Monte Carlo results are new (§2
above) since #1/#2 previously prevented a valid run on the Lodders
composition.

### 2026-09-14 session: adding P/Schreibersite, and what that uncovered

The requested change (add "P" to `elements_R`, "Schreibersite" to
`phases_R`, a Schreibersite row to `phase_lib`) was genuinely mechanical
in itself — see §1 for the result, byte-identical to the pre-P fit on
every phase except Iron, exactly as predicted. Three things found while
doing it were not mechanical and are logged here.

**[X→V] The "final" Lodders table in §1 was never actually produced by
any tracked file.** STATUS.md documented a Lodders-based headline result
as of 2026-09-13, but the main .Rmd's §12 "Run" still called
`fit_scenario(phases_R, b_CI, elements_R)` — `b_CI` is the Alexander
target. The 2026-09-13 table had in fact been computed by an untracked
ad-hoc script, not by this project's own files. Running the file as it
stood would have silently regressed the headline result back to
Alexander and called it final. Fixed: `b_bulk_CI_LoddersT4` and
`recompute_O()` are now defined in the main .Rmd itself (previously only
in the montecarlo .Rmd — moved, not duplicated, to keep one source of
truth), `b_CI_Lodders <- drop_TiCo(recompute_O(b_bulk_CI_LoddersT4))` is
computed there, and §12 now fits `b_CI_Lodders` as `fit_R` (primary),
with a separate `fit_R_Alexander` computed alongside for the bracket
comparison. This is not a new modelling decision — Lodders-primary was
already decided 2026-09-13 — it is making the code match a decision that
had only been made in prose.

**[X→V] `to_volume()` bug: a phase with a known `phase_lib` entry but
Vm = NA silently corrupted every other phase's vol%, not just its own.**
The function's only exclusion check was for phases missing from
`phase_lib` entirely; a phase present with `Vm = NA` (exactly
Schreibersite's situation) passed that check, then `rho = M/NA`,
`cm3 = .../rho`, and `sum(cm3)` all became `NA`, making `vol_pct = NA`
for **every** phase in the table via `100 * cm3 / sum(cm3)`. This is the
"silently produces a number" failure mode the handoff explicitly warned
against, except worse — it would have looked like a crash-adjacent NA
column, not a wrong number, so probably would have been caught, but it
was still a real bug in shipped code. Fixed: `to_volume()` now checks
Vm-per-phase separately from library membership, excludes NA-Vm phases
from both the numerator and the `sum(cm3)` denominator, and warns
naming them (verified: `to_volume(fit_R)` now warns "No Vm for:
Schreibersite" and returns a clean 11-row table with no NA).

**[X→V] Montecarlo diagnostic chunk was validating `fit_R` against the
wrong composition.** `_montecarlo.Rmd`'s diagnostic chunk (chunk 3)
compared `norm_closed(b_CI)` (Alexander closed-form) against `fit_R`'s
NNLS output. This was correct when `fit_R` targeted `b_CI` (before this
session); once `fit_R` was repointed at `b_CI_Lodders`, the same code
would have compared two different compositions and reported a spurious
~0.08-0.17 mg/g "mismatch" that is really just the Alexander-vs-Lodders
gap, not a fit problem. Fixed: the check now uses `norm_closed(b_CI_Lodders)`
against `fit_R`, plus a separate Alexander-bracket line kept for
reference. Re-run result: closed form − NNLS now agrees to **0.003–0.005
mg/g** (tighter than the 0.082/−0.173 mg/g this file previously reported
for this exact check, because that number was measuring the wrong thing).
Also surfaced: the true Lodders-composition crossover distance is
Mg_rem/Si_rem = 1.4139 (0.18% from the 1.4164 crossover), not 1.4025
(that is the Alexander value) — §2's "What to report" bullet has been
corrected, since it previously gave "1.4025 actual" with no composition
label attached.

None of this reopens §2 (Fo/En is still reported as ~50/50 with no
probability figure) or changes any decision from the 2026-09-13 session
— it makes the code match decisions already recorded in prose, and fixes
one latent bug that P/Schreibersite happened to expose.

### 2026-09-15 session: two verification checks, one expectation corrected

**[V] Check (a): the "wrong composition" diagnostic, re-derived fresh,
side by side.** The 2026-09-14 notes above (last paragraph) asserted
that the OLD diagnostic mismatch was ~0.08/−0.17 mg/g and that fixing it
tightens the number to ~0.003/0.005 mg/g. Both halves of that claim were
checked again today, from scratch:

    norm_closed(b_CI_Lodders) - fit_R$mg_g[Fo,En]  =  0.0032, 0.0046  mg/g
    norm_closed(b_CI)         - fit_R$mg_g[Fo,En]  =  -5.49, +5.19    mg/g

The first line confirms the 2026-09-14 fix (0.003–0.005 mg/g, as
claimed). **The second line does NOT reproduce "~0.08/−0.17" as this
file's own prose predicted — it is roughly 30x larger.** This is not a
new bug: `fit_R` today targets `b_CI_Lodders` *and* includes
Schreibersite (both added since the 0.08/−0.17 number was first
measured), so `norm_closed(b_CI) - fit_R$mg_g` no longer isolates "the
old bug's residual" — it now ALSO contains the full Alexander-vs-Lodders
Fo/En gap (Alexander NNLS − Lodders NNLS = −5.58, +5.36 mg/g, matching
§1's own "old vs new" table to rounding). −5.58 + ~0.08 loosely
recomposes to ≈−5.49, i.e. the old 0.08/−0.17 signal is still in there,
just swamped by a ~30x larger, already-documented, unrelated effect.
**Correcting this file's own prediction, not the code**: a handoff that
asks to "reproduce ~0.08/−0.17" against today's `fit_R` is describing a
comparison that no longer isolates what it once did. Re-derive fresh
before trusting any such expectation stated in prose (per this file's
own ground rules) rather than assuming an old figure still applies once
`fit_R`'s target has changed again.

**[V] Check (b): which phases are structurally forced vs. jointly fit.**
Worked through `build_A(phases_R, elements_R)` column by column
(element → unique vs. shared hosts). Every phase in Scenario R is
**forced** by a unique-element elimination chain, except one pair:

| Phase | Forced by (unique host) | Or: jointly fit |
|---|---|---|
| Magnesiochromite | Cr (only host) | — |
| Tephroite | Mn (only host) | — |
| Albite | Na (only host) | — |
| K-Feldspar | K (only host) | — |
| Troilite | S (only host) | — |
| Schreibersite | P (only host) | — |
| Nickel | Ni (only host) | — |
| Anorthite | Al, net of Albite+K-Feldspar's share (forced once those two are fixed by Na, K) | — |
| Diopside | Ca, net of Anorthite's share (forced once Anorthite is fixed) | — |
| Iron | Fe, net of Troilite+Schreibersite's share (forced once those two are fixed) | — |
| **Forsterite** | — | **yes — Mg_rem & Si_rem, 2 elements/2 phases** |
| **Enstatite** | — | **yes — same 2x2 system as Forsterite** |

This is a stronger statement than §1's existing "near-triangular"
caveat: it is not just P→Schreibersite and S→Troilite that are forced —
**the entire 12-phase system reduces to a triangular elimination chain
with exactly one genuine 2-variable joint fit at the end (Fo/En)**, which
is exactly the pair §2 already treats with the closed-form derivation
and the correlated-error Monte Carlo. Every other phase's "fit" result
is arithmetic, not evidence of a multi-element mass-balance success —
consistent with, and now fully explaining, why the residuals are
machine-precision zero (§1) rather than merely small.

---

## 7. Phase roles and database status

| Phase | Role | Thermo | Kinetics | Proxy |
|---|---|---|---|---|
| Forsterite, Enstatite, Diopside, Anorthite, Albite, K-Feldspar | reactant | n/a | PK04 | direct; albite → Oligoclase for An15 |
| Fayalite | reactant (C) | yes | fixed, see §6 | — |
| Troilite | reactant | yes | none | pyrrhotite PK04; **SR("Troilite")**, not "Pyrrhotite". Vacancy-bearing pyrrhotite dissolves faster → over-predicts early S |
| Greenalite | reactant (C) + product | yes | none | lizardite PK04; Fe-serpentine likely faster → biased slow |
| Native Fe, Ni | reactant | yes | none — **[U] checked 2026-09-15, §4: no source found** (Core11 and PK04 both lack native-metal corrosion kinetics; not searched further) | EQUILIBRIUM_PHASES dissolve_only; corrosion not rate-limiting is a standard qualitative corrosion-science argument, not a sourced number. **[U]** verify dissolve_only syntax; **[U]** composition with secondary Fe minerals unchecked; **[U]** passivation not ruled out |
| Tephroite | reactant | yes (pre-existing) | own rate law (2026-09-15, §4): acid branch only, log k0=−4.46 (pre-exp; Pincus/Casey's pH-2-referenced value is −5.4, see §4 for the distinction), n=0.47, Ea unresolved (null placeholder), anoxic conditions only | Pincus et al. (2026)/Casey et al. (1993), no longer fayalite |
| Magnesiochromite | reactant | ? | **(none) — checked, not just absent** (2026-09-15, §4): 3 papers reviewed (Razavi 2026, Ao 2024, Schindler 2018), none give a usable ACT("H+")^n rate law; Schindler's own finding is that chromite is refractory and not rate-controlling in a natural matrix | treat as effectively inert at 0.71 wt% of the assemblage, or bound arbitrarily slow if a nonzero rate is needed for completeness |
| Schreibersite | reactant | none | none, rate TBD (§4) | RATES-only, no EQUILIBRIUM_PHASES entry, no SR() term — placeholder log_k never evaluated |
| Magnetite, Lizardite, Pyrrhotite, saponites, nontronites | product | yes | n/a | — |
| GEMS (amorphous) | reactant | **none by design** | forsterite anchor × f_enh | never in EQUILIBRIUM_PHASES; no SR() term |
| Fe-Ni nanoinclusions in GEMS | reactant | — | separate phase | keep out of GEMS stoichiometry |

Placeholder log_k in RATES-only PHASES entries is never evaluated
**provided** the rate law omits SR(). Document this in each entry.

---

## 8. Open items

1. **Schreibersite reaction chemistry — RESOLVED 2026-09-14** (was item
   1). See §4. Two sub-items remain open, tracked separately:
   - **(1a) Vm — RESOLVED 2026-09-15.** Handbook of Mineralogy, Fe3P
     endmember (matching phase_lib's formula exactly), D(calc.) = 7.44
     g/cm3 (Fasiska & Zweil 1967, Trans. AIME 239, 924-925). Vm =
     198.509/7.44 = 26.68 cm3/mol, computed in the script from
     `formula_mass("Fe3P")`, not hardcoded. The Handbook's other three
     analyses (D = 7.0-7.3 g/cm3) are Ni-bearing natural schreibersite
     (18-43 wt% Ni) and were NOT used, since mixing a Ni-substituted
     density with the Ni-free Fe3P formula in phase_lib would be a basis
     mismatch. **[U]**: if a Ni-bearing (Fe,Ni)3P composition is ever
     adopted instead (the Fe-Ni metal alloy already carries ~11 wt% Ni,
     §1), this density and the phase_lib formula both need revisiting
     together.
   - **(1b) Rate constant — PARTIALLY resolved 2026-09-15, corrected
     same day.** An earlier pass this session wrongly claimed no Fe3P
     sample mass exists anywhere in the paper (it checked only the
     H2-quantification passage, not Table 2 — see §4 for the struck
     record). Table 2 DOES give sample mass (1.2197 g, Ar-flushed DI
     water, pH 7, 20C — matching the H2-release measurement) for all 11
     rows, but NO grind size, mesh fraction, or BET surface area for any
     row. Net result: a **bulk (mass-normalized) fractional dissolution
     rate is derivable — 1.18%/day, range 0.59-1.78%/day propagating the
     stated H2 uncertainty** — but a mechanistic (surface-area-
     normalized, mol m-2 s-1) rate is still not extractable and still
     should not be estimated (§4 has the full derivation and the
     recommended, heavily-caveated way to use the bulk rate if one is
     needed). Not fully open (there is now a usable number), not fully
     resolved (it isn't the form most PHREEQC rate laws in this database
     use).
   - **(1c) NEW 2026-09-15: P(V)-only justification not fully
     reconciled.** §4's "P(V) is reasonably justified" claim cites only
     Pasek & Lauretta's Fig. 4 (Fe-metal catalyzes phosphonate
     oxidation). The same paper's saline (Mg/Ca-chloride) solution
     experiments — the closer analog to this project's own assumed
     Cl-bearing initial fluid (item 3 below) — show phosphonate as 100%
     of detected P, i.e. a persistently REDUCED product. Not yet checked
     against each other (see §4 for the detail: possibly different,
     non-conflicting conditions, but unverified). Do not state P(V) as
     "reasonably justified" without resolving this.
2. **GEMS**: composition from `amorphous_si_o()` once run; confirm Si:O
   ≤ 0.28 (forsterite 0.25, enstatite 0.333) so f_enh defaults to 1.
   **Grain size has NO source** — needs the TEM literature. f_enh from
   Wolff-Boenisch Table 1: use the MEASURED forsterite row against the
   glass regression, never the mineral regression.
3. **Chlorine**: absent from Alexander Table 6 and Charnoz's element set.
   Chlorapatite cannot host it (all 1.45 mg/g P as chlorapatite yields
   only 0.55 mg/g Cl and needs zero schreibersite). **[D]** Put Cl in the
   initial fluid, following Neveu et al. (2015), 520 mol ppm from Wasson
   & Kallemeyn (1988).
4. **RHO_CI_GRAIN — RESOLVED 2026-09-15** (was "RHO_CI_BULK = 2.42
   unsourced placeholder"). The old name was physically backwards, not
   just unsourced: King et al. (2015)'s XRD modal percentages are volume
   fractions of the SOLID material only — porosity is invisible to XRD —
   so the correct conversion constant is **grain density**, not the
   porosity-inclusive bulk/porous density a name like "BULK" suggested.
   Renamed `RHO_CI_BULK` → `RHO_CI_GRAIN` everywhere (grepped the whole
   project; the string no longer appears in any live file), value set to
   **2.46 g/cm3** — Friedrich (2024), ideal gas pycnometry (the method
   the paper itself identifies as most reliable per Britt & Consolmagno
   2003; CCl4/water immersion methods underestimate grain density and
   were not used): Ivuna 2.55, Orgueil 2.44, plus Flynn et al. (1999)
   2.46 on an independent Orgueil sample; sample-mass-weighted mean 2.46
   ± 0.05. `RHO_CI_GRAIN_sigma <- 0.05` recorded alongside for future MC
   propagation, not yet wired into anything. Friedrich's BULK/porous
   density (1.61-2.10 g/cm3, ~15-34% macroscopic porosity plus additional
   unresolved microporosity) is explicitly NOT what this constant should
   be — using it would inflate every bulk-basis vol% to account for pore
   space King et al.'s XRD never measured.

   **[U] The CMWG Report on Material Density and Porosity cross-check
   was NOT performed.** The handoff asked for this document to be
   checked "already in Project Knowledge" — that is a claude.ai
   knowledge-base feature belonging to a different (remote) session; it
   is not accessible from this environment, and no local copy of the
   report exists in this project's files (checked). 2.46 g/cm3 is used
   on Friedrich (2024) alone. If the CMWG report gives a different or
   additional CI grain density, it has not been reconciled against this
   value — flagging for the user/other session to check, per the
   handoff's own instruction not to resolve this silently.
5. **Surface areas** for every kinetic phase — that parameter set does
   not exist yet, not just for GEMS.
6. **Albite in Charnoz type C?** (§4)
7. **Fayalite Vm** vs Core11.
8. **Forward check against King et al.** — hydrate the Scenario R
   assemblage and test against 81-84 vol% phyllosilicate, 6-10 magnetite,
   4-7 sulfide. This IS the PHREEQC model and the only genuinely
   independent test of the reconstruction, and the actual next step now
   that §1's phase assemblage is closed. **[U] 2026-09-15: attempted and
   deferred** — this machine has only "Phreeqc Interactive" (GUI-only,
   no console/batch mode, no IPhreeqc COM, no R `phreeqc` package), so
   PHREEQC cannot be run headlessly here. No input file was drafted; the
   task is fully open, not partially done.
9. **Phase-list sensitivity** (variants with/without anorthite,
   wollastonite for diopside, alkalis to amorphous, Mn to olivine) —
   cheap, and the split in §2 is phase-list dependent.
10. **RESOLVED, 2026-09-13: composition source decided.** Lodders et al.
    (2025) Table 4 is primary; Alexander (2019) is the bracketing
    comparison only. `b_bulk_CI_LoddersT4` / `sigma_rel_Lodders` are the
    active choice in the montecarlo file; `b_bulk_CI_AlexanderT1` /
    `sigma_rel_assumed` are kept solely for the comparison in §1/§2.
11. **Native Fe, Ni kinetics — checked 2026-09-15, still open** (§4,
    §7). No source found in Core11 or PK04 for aqueous Fe/Ni metal
    corrosion; current `dissolve_only` treatment rests on a qualitative
    corrosion-science argument (electrochemical corrosion typically fast
    relative to silicate dissolution), not a sourced rate. Three
    sub-items, none checked yet: (a) `dissolve_only` syntax for the
    actual PHREEQC build; (b) whether it composes correctly with
    secondary Fe minerals (magnetite, sulfides) rather than trapping Fe;
    (c) passivation not ruled out for this specific composition/fluid.
    If a real rate is ever needed, look at native Fe-Ni meteorite
    weathering or steel/iron corrosion-rate literature — neither
    searched yet.

Sourcing sigma(Mg/Si) further and investigating the Fo/En point
probability are deprioritized, not open — see the closing note in §2.

## 9. Process notes
- The Scenario R fit is square: near-zero residuals are structural, not
  validation. Say so in the methods before a reviewer does.
- Cross-check numeric results against a second LLM.
- Ask "did you check this, or is it recalled?" — every fabricated value
  in this project was presented with the same confidence as verified ones.
- **Recompute every prose number from actual code output.** In v0 the
  prose said metal 8-12 vol% while the code returned 1.3.
- When a prediction and a simulation disagree, suspect the prediction:
  in §2 the MC was right twice and the analytic prediction wrong twice.
- **2026-09-15: a careless PowerShell text-replace corrupted this file's
  encoding** (every non-ASCII character — §, —, →, ×, Δ — was mangled
  into multi-byte garbage via an unintended double-encoding round-trip).
  Caught immediately because the file is read after every edit, not
  assumed correct. Reconstructed from this session's own edit history
  rather than left corrupted or silently "fixed" with a lossy algorithmic
  reversal that was tried first and did not work cleanly. Lesson: prefer
  UTF-8-safe tools (this project's own Edit/Write tooling) over shell
  text-replacement one-liners for any file containing non-ASCII text.
