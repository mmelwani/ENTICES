# PROGRESS — IOM module, redox treatment, and the state of ENTICES

Living document: the current state. History lives in git and in `docs/reviews/`.

**Last updated: 2026-10-26.** This revision folds in the review exchanges of
2026-09-30 → 2026-10-26 (all in `docs/reviews/`) and the merged branches
`fixes`, `consistency`, `pairing` and `eq-parity` (merged as `255981f`). It **overturns**
three conclusions of the 2026-09-20 version:

- `Hdg` carbon decoupling is **not** working in production. It works in the
  tests but acts on nothing in a real run (Part 3).
- The bulk IOM composition is **CR-chondrite IOM**, Miller et al. 2026 SI
  Table 9. It is not CI- or kerogen-derived (Part 4).
- The 350 C residue check **passes for H, but O is under-released by
  13–30%**. It is not "a starting-composition artefact" (Part 4).

Scope has widened: the rock, the databases, the gas phases and the initial
fluid now matter as much as the IOM. Parts 5–9 are new. **Part 10 is the
decision list**: start there if you only read one section.

Reading guide:

| Part | Topic | Status |
|---|---|---|
| 0 | What changed since 2026-09-20 | new |
| 1 | What the IOM channels release | unchanged, still the most important concept |
| 2 | The low-temperature speciation problem | unchanged; diagnostic status updated |
| 3 | Option 3 (redox decoupling) | **rewritten** |
| 4 | IOM calibration | composition and validation updated |
| 5 | Rock and primary minerals | new |
| 6 | Databases and numerics | new |
| 7 | Gas phases and H₂ | new |
| 8 | Initial fluid, seed, attribution | new |
| 9 | Code state | new |
| 10 | **Decisions and next actions** | replaces the old Part 5 |
| 11 | Process lessons | new |

Companion documents (paths on `main`):
- **`docs/PATCH_F_vitrimat2018.md`:** A = 2e15, distributed Ea, the 250 C
  validation.
- **`docs/PATCH_F_addendum_N_and_CHn.md`:** nitrogen and the CHn channel.
- **`docs/OPTION3_decision.md`:** the per-channel decoupling design. Its `Hdg`
  premise is superseded by Part 3 here.
- **`docs/CR_IOM_COMPOSITION_RESOLVED_2026-09-20.md`:** the composition
  provenance.
- **`docs/PHREEQC_TEST_RESULTS_2026-09-20.md`** and
  **`docs/PHREEQC_TEST_B3_real_assemblage_2026-09-27.md`:** Tests A, B2
  and B3.
- **`reconstruction/RECONSTRUCTION_STATUS.md`:** the primordial rock.
- **`CHANGELOG.md`:** the code change log (Lucas's and ours).

---

# PART 0 — WHAT CHANGED SINCE 2026-09-20

| Finding | Effect | Part |
|---|---|---|
| Production H₂ is generated in situ by Fe/Ni metal corrosion, in coupled form. `Hdg` only decouples H₂ *entered* as an inventory, and no production run enters any | Option 3's carbon fix is a no-op as specified | 3 |
| `Hdg` raises pe to the mineral buffer, and at pH ~11.5 that oxidises sulfide to sulfate (22% in B3) | A second cost of `Hdg`, including for rock-derived sulfide | 3 |
| `coreclath_CH4uncoupled.dat` (used below 150 °C) has no coupled CH₄ and no C(−4). But C₂H₆, C₂H₄, acetate etc. remain | A *partial* "remove reduced carbon" design is already in production; its leak paths are untested | 3, 10 |
| Composition = Miller 2026 SI Table 9 (CR IOM) exactly | Composition settled; CR-composition / CM-kinetics split documented | 4 |
| Residue, like-for-like with Murchison: H passes, O under-released 13–30% | Validation status corrected | 4 |
| Rock = reconstruction (Lodders 2025, Scenario R), now generated from `reconstruction/rock_mol_per_kg.csv`. Schreibersite added (kinetic), magnesiochromite inert, lawrencite kept with its Fe taken from metal | Rock settled | 5 |
| The 10-06 convergence hang was tephroite precipitating (no `SR > 1` guard). It wasn't feldspars or gases | Fixed; full runs take minutes | 6 |
| EQUILIBRIUM_PHASES gases at `0 0` are a 1 atm infinite reservoir, saved with the porewater and never exchanged with the ocean | **Every ocean H₂ number so far is set by this artefact** | 7 |
| Aqueous `-Vm` is inactive under the LLNL model: aqueous speciation has no pressure dependence; minerals and gases do | Pressure matters mainly for gases and mineral SIs | 6 |
| The legacy `legacy/ENTICES_v3.R` `REACTION` dosing was off by a factor of `iom_C_moles` (35.6× under-dose in the redox-diagnostic case) | All v3 organic results unreliable; the redox diagnostic is legacy | 6 |
| Pure-water start kept on purpose (attribution). Trace C seed 1e-14 mol/kgw for 25 °C numerics. Paired runs with a same-rock blank give the organics' contribution | Attribution method in place | 8 |
| PHREEQC version matters for minor species (Mn 130× between 3.5.0 and 3.8.6) | Standardised on 3.8.6; stamped in every `.pqi` | 6 |

**Corrections to statements in this file's 09-20 version:**
- "`Hdg` carbon decoupling — CONFIRMED". That holds for the mechanism, but
  not for production (Part 3).
- "Calcite… 8.16e-04 mol". The 66.7 atm value is 8.055e-04.
- "Magnetite 2.9×, pyrrhotite 3.6×" under `Hdg`. This is a toy-assemblage
  result only. B3 was Fe-limited and couldn't measure it, and pyrrhotite isn't
  in production lists. Unconfirmed.
- "Residue partially resolved, +7.8%". See Part 4.
- Part 5's 09-20 to-do list: all done or superseded (Part 10).

---
# PART 1 — WHAT THE IOM CHANNELS RELEASE (read this before touching anything)

This section exists because the answer is counter-intuitive and has been
re-asked more than once. It is the single most important thing to understand
about the module.

## The channels release ELEMENTS, and (in coupled mode) the solver decides the products.

Every channel's `-formula` — for example `IOM_CO2`'s `C 1 H 1 O 1.3115` —
adds **elemental totals** to the PHREEQC solution. For a normal (coupled)
element, PHREEQC then decides what aqueous species those elements become, by
thermodynamic equilibrium at whatever pH and pe the solution happens to have
at that instant.

The one exception, which Part 3 exploits deliberately: if a formula names a
**decoupled pseudo-element** (`Mtg`, `Sg`, `Ntg`, `Hdg`, `Oxg`), that material
has its own mass balance and does *not* re-equilibrate. So "elements, never
molecules" is the right description of coupled mode, but it is not a hard
limit of PHREEQC — it is a property of which elements you use.

The channel *names* (`IOM_CO2`, `IOM_CH4`, `IOM_N`, `IOM_S`, `IOM_CHn`) are
labels for **which elements, in what ratio, released at what rate**. They do
**not** force the named molecule to appear in solution. `IOM_CO2` is not "the
channel that makes CO2"; it is "the channel calibrated against the population
of bonds whose breakdown Miller measured as CO2 yield."

## Why it is built this way

**CORRECTED 2026-09-18.** An earlier version of this section gave two reasons,
the first of which was wrong. It claimed PHREEQC "mechanically cannot do
otherwise" and that "there is no syntax by which it adds a molecule." Checked
against the PHREEQC v3 manual (KINETICS, `-formula`), both parts of that are
false:

- **Phase and molecular names ARE accepted.** The manual: "*formula — Chemical
  formula **or the name of a phase** to be added by the kinetic reaction... A
  phase name may be entered independent of case.*" Its own examples use
  `-formula FeS2 1.0 FeAs2 0.001` and `-formula CH2O(NH3)0.1 0.5`. So writing
  `-formula CO2` is perfectly legal syntax.
- **And persistence IS achievable** — via decoupled pseudo-elements (`Mtg`,
  `Sg`, `Ntg`, `Hdg`, `Oxg`), which have their own mass balances and do not
  re-equilibrate. That is the entire basis of Part 3. The old wording
  contradicted this document's own Option 3 design.

What is true, and what the manual does support, is narrower:

**1. Naming a molecule buys you its stoichiometry, not its identity.** The
manual is explicit that when a phase name is used, "*the formula for that phase
is then used for the **stoichiometry** of the reaction*," and that a formula
"*may be considered as adding or removing **native elements** from the system in
the given stoichiometry.*" So `-formula CO2` adds 1 C + 2 O as elements, and
those elements then speciate by equilibrium exactly as `C 1 O 2` would. Writing
the molecule's name does not make the molecule persist — only a *decoupled
element* does that.

**2. The load-bearing reason is chemical, not mechanical.** Even setting syntax
aside, the product molecule is the wrong stoichiometry to hand PHREEQC.
Checking Miller et al.'s own mass balance at 350C: the residue loses
4.50 mol O/kg, but the measured CO2 requires 6.08 mol O/kg. Roughly a quarter
of the oxygen in the product CO2 does not come from the kerogen at all — it
comes from the water, because these are hydrous pyrolysis experiments. So
`-formula CO2` would assert that the solid gives up 2 O per C, which the
residue data show it does not. The physically correct quantity is **what the
solid gives up** (`C 1 H 1 O 1.3115`), letting PHREEQC draw the remaining
oxygen from water.

This is why the formulas look the way they do. The fractional oxygen is the
*solid's* measured oxygen loss per carbon; the balance needed to make CO2 — if
CO2 is what forms — comes from water via PHREEQC's own H/O bookkeeping.

**Do not "simplify" these formulas to product molecules.** It is the most
likely well-intentioned change someone would make. It is legal syntax, it will
run without error, and it will be quantitatively wrong on the oxygen budget.
That combination — legal, silent, wrong — is what makes it worth this much
space.

## The consequence, and why it is now a problem

Because speciation is decided by equilibrium, at Enceladus pore-water
conditions (reducing, H2-rich from mineral reaction, high pH) **essentially
all released organic carbon speciates as CH4 — regardless of which channel
released it.** Carbon released by the carboxyl-derived `IOM_CO2` channel
arrives as elements and is immediately converted to methane by the solver,
because that is what is thermodynamically favoured there.

That is correct PHREEQC behaviour and it is what the model was asked to do.
**It is not correct physics at low temperature.** See Part 2.

---
# PART 2 — THE LOW-TEMPERATURE SPECIATION PROBLEM

## Statement of the problem

The CO2/CH4 redox couple is famously sluggish. Abiotic interconversion between
oxidised and reduced carbon requires high temperature, effective mineral
catalysis, or biological mediation. Below roughly 150-200C, CO2 and CH4 do not
equilibrate with each other on laboratory timescales and often not on
geological ones either — they persist in metastable coexistence. This is
standard in the hydrothermal organic geochemistry literature (McCollom,
Seewald, Shock and co-workers on metastable equilibria). *Caveat: these
citations are from memory and were never verified from a source; confirm
before they go in a paper. The underlying physical result is not in doubt; the
specific attributions should be checked.*

Nitrogen and sulfur couples are likewise kinetically inhibited abiotically:
NH4+/N2/NO3- interconversion, and HS-/SO4-2 equilibration, are both very slow
at low temperature without microbial mediation.

**However — see `docs/OPTION3_decision.md` — this does not mean all four channels
need decoupling.** What matters per element is whether coupled equilibrium
puts the released material where kinetics would have left it. For nitrogen it
does; for carbon it emphatically does not. That asymmetry is the basis of the
revised Part 3.

## Two consequences, one of which is quantitatively severe

**(a) Predicted carbon speciation is wrong.** The model predicts an ocean with
essentially no oxidised carbon. Cassini observes both CO2 and CH4 in the
plume. A fully-coupled model at low pe *cannot* reproduce coexisting CO2 and
CH4 — not because the chemistry is wrong, but because the model forbids the
metastability that reality exhibits.

**(b) The coupled model imposes a large spurious electron demand on the whole
system.** Forcing released carbon down to CH4 requires electrons, and PHREEQC
takes them from wherever the system's electron reservoir is — dissolved H2,
*and* oxidation of reduced minerals (Fe(II) -> Fe(III), magnetite formation,
pyrrhotite oxidation). It is **not** specifically "H2 from the IOM being
destroyed"; the IOM is the sink, the rest of the geochemical system is the
source.

Magnitude, computed from the actual released stoichiometry (not from assuming
the carbon arrives as CO2):

    channel    C oxidation state as released    e- to reach CH4
    IOM_CO2    +1.623                           +5.623 per mol  -> +23.44 mol e-/kg IOM
    IOM_CH4    -4.435 (already below CH4)        -0.435 per mol  ->  -0.83 mol e-/kg IOM
    IOM_N      N at -2 -> NH4+ at -3                             ->  +0.88 mol e-/kg IOM
    IOM_S      S at -1 -> HS- at -2                              ->  +0.46 mol e-/kg IOM
    -------------------------------------------------------------------------------
    NET                                                    ~24.0 mol e- per kg IOM
                                                  = ~12.0 mol H2-equivalent per kg IOM
       = 0.60 mol H2-eq per kg rock at 5 wt% organics
       = 1.80 mol H2-eq per kg rock at 15 wt% organics

(An earlier version of this note said 16.7 mol H2/kg IOM. That assumed the
released carbon was already at +4, i.e. CO2. It is not — the solid releases
carbon at about +1.6, because only ~1.3 oxygens come off per carbon.)

Note the target: ~24 mol e⁻/kg IOM is what it takes to drive all released
carbon to **CH₄**. Driving all of it to **CO₂** instead needs
(4 − (−0.72)) × 7.657 ≈ **36 mol e⁻/kg IOM**, because the three carbon channels
release carbon at a weighted mean oxidation state of −0.72. Which one applies
depends on where equilibrium puts the carbon.

**Status of the redox-budget diagnostic (2026-10-26).** It was run on 09-20,
but on legacy code: `legacy/ENTICES_v3.R`, with its 35.6× dosing bug and the
superseded rock (Part 6). So it shows only that old output was contaminated.
**It has not been redone on current code.** It is covered by test B4 case 1
(Part 10, N2). Still check the whole budget (H₂, pe, Fe phases), not H₂ alone.
And until Part 7 is fixed, H₂ is capped by the 1 atm gas reservoir, so H₂ alone
can't reveal the demand.

## Therefore: Option 3 moves from "interesting bracket" to "required"

Previously scoped as a nice-to-have bracket on the ocean's degree of redox
equilibrium. It is now the only way to model the low-temperature cases
correctly. The bracket interpretation still stands — but the decoupled run is
no longer the exotic end-member; at low temperature **it is the physically
appropriate one**, and the coupled run is the artefact.

**Still true, but the mechanism chosen in 09-18 (`Hdg`) doesn't act in
production.** See Part 3.

---
# PART 3 — OPTION 3: REDOX DECOUPLING (REWRITTEN 2026-10-26)

## What still stands from the 2026-09-18/20 design

- **The decoupled pseudo-elements.** The databases contain five, each with its
  own mass balance:

      Hdg (H2)   Oxg (O2)   Mtg (CH4)   Sg (H2S, species H2Sg)   Ntg (N2)

  There is no decoupled ammonium and no decoupled oxidised carbon. `C(-4)`,
  `N(-3)` and `S(-2)` are redox *states* of the shared elements, not
  decoupled species.
- **N stays coupled.** At Enceladus pore conditions NH₃/NH₄⁺ is the stable N
  species, so coupled equilibrium puts the released N(−2) where kinetics would.
  `Ntg` would oxidise it to N₂, which is the wrong direction.
- **The decoupled formulas are correct:** `IOM_CH4` → `Mtg 1 H 0.435`,
  `IOM_S` → `Sg 1 H 1`. Both were tested in PHREEQC (Test A4).

## One criterion for every decoupled formula: electrons delivered to the system

Count electrons with C ending as CH₄ (−4), S as sulfide (−2) and H as
H(+1). `Mtg` enters as neutral CH₄-equivalent (0 e⁻); `Sg` enters as H₂Sg
(S −2, which needs 2 e⁻).

    formula                 net e- to system   verdict
    C 1 H 4.435 (coupled)        +0.435
    Mtg 1 H 0.435                +0.435        correct
    Mtg 1 H 4.435                +4.435        wrong (Test A4: 3.8x excess H2)
    S 1 H 1 (coupled)            −1.000
    Sg 1 H 1                     −1.000        correct
    Sg 1 H −1                    −3.000        wrong (pe +15.4, H2 destroyed)
    Sg 1                         −2.000        wrong (pe +15.3)

The 09-20 "opposite rules for Mtg and Sg" was this one criterion described
twice. Use the electron count for any future pseudo-species. It is validated
by the H₂ budget, not by pe: in Test A4 the H₂ buffer was exhausted, which
makes pe ill-conditioned.

## Why `Hdg` doesn't work in production

**The design (09-18).** Enter the H₂ inventory as `Hdg`, so that H₂ can't
reduce released carbon. pe is then set by the mineral assemblage. **The tests
(B2, B3)** supplied H₂ as an *initial-solution inventory* with toy element
amounts. They confirmed the mechanism:

| Test | Assemblage | Suppression of C → CH₄ |
|---|---|---|
| B2 | Magnetite + pyrrhotite | ~26,000× |
| B3 | Lucas's real 78-phase ocean list, 66.7 atm | ~1.8×10⁵× |

The B3 factor depends on the arbitrary inventories, so report it as "several
orders of magnitude". Calcite still precipitates under `Hdg` (8.055e-04 mol
at 66.7 atm).

**Production never enters H₂.** Lucas's inputs start from pure water. All H₂
comes from corroding native metal (~2.0 mol Fe and 0.24 mol Ni per unit rock
at 15 wt% organics), written in coupled form:

    Fe + 2 H+ + 0.5 O2 = Fe+2 + H2O

So every electron from corrosion enters the coupled pool. Keeping that H₂
coupled is deliberate: `reconstruction/RECONSTRUCTION_STATUS.md` §4 [D] requires it, so that
metal and schreibersite H₂ stay on the same books as the ~3.1 mmol H₂/g
habitability figure. **As specified, `Hdg` therefore acts on nothing in a
production run.**

The artefact it was meant to fix is real in production output. In Lucas's
150 °C example run, pe is −11.69 and CH₄/CO₂ ≈ 1.7×10⁹.

**Rewriting the metal to release `Hdg` looks worse.** Untested, but from what
is measured:
- **It contradicts the reconstruction decision** above.
- **It oxidises sulfide.** With H₂ hidden, pe rises to the mineral buffer. At
  pH ~11.45 the sulfate/sulfide boundary is near pe −8.2, which is where B3
  sim 4 sat: 22% of the S(−2) went to sulfate, and that S was the stand-in
  for rock sulfide (Troilite releases coupled HS⁻).
- **It probably fails to hold carbon anyway** *(inferred)*. Secondary
  Fe²⁺ → magnetite still releases coupled H₂, and at low T even buffer-level
  H₂ reduces CO₂.

**The Fe/S mineral perturbation under `Hdg` is unmeasured.** The 2.9×/3.6×
figures came from a 2-mineral toy. B3 was Fe-limited: pyrite took all
1.000e-05 mol Fe in both modes, so it can't size the effect.

## Design (d): remove coupled reduced carbon, keep H₂ coupled

**The idea.** CO₂ can't be reduced if the coupled C element has no reduced
species. Carbonates are unaffected (C(+4) stays coupled), H₂ stays coupled,
sulfide stays stable, and methane exists only as `Mtg`.

**It's partly built already.** `coreclath_CH4uncoupled.dat`, the default below
150 °C, already has no coupled CH₄, no `C(-4)` master and no `CH4(g)`, and its
methane clathrate is `Mtg_hydrate`. But the other reduced-carbon species
remain (Claude Code's 10-23 section diff):
- the masters `C(-3)`, `C(-2)` and `C(+2)`;
- C₂H₆, C₂H₄, CO, formate, acetate (48 lines including complexes), CH₃OH and
  others.

So coupled CO₂ may simply be reducing to ethane or acetate instead of methane
below 150 °C. **This is untested and is the cheapest high-value test
available** (Part 10).

**The full scope** (09-30 count, on the coupled file):
- ~73 aqueous reactions;
- 4 secondary master species;
- 8 phases in the production lists: CH₄(g), C₂H₄(g), C₂H₆(g), C₃H₈(g),
  CO(g), Glycine, Citric_Acid and CH4_hydrate.

**Costs and choices:**
- **Glycine, citric acid and the light hydrocarbons** look like deliberate
  organic products, so removing them is a modelling choice.
- **Formate** equilibrates with CO₂ + H₂ far faster than methanation, so it is
  a candidate to re-admit as a sensitivity case.
- **`IOM_CHn` (`C 1 H 1.8`)** needs a home: `Mtg 1 H -2.2`, drop the channel,
  or a decoupled C₂⁺ pseudo-species.
- **Required check: a leak check.** After any edit, every C valence state other
  than C(+4) must total 0.

## What this means for `IOM_S`

- **If carbon is handled by (d):** H₂ stays coupled, pe stays low and sulfide
  stays sulfide. Then `IOM_S` can be **coupled**, which gives IOM sulfide its
  mineral sink (pyrite/troilite).
- **If `Hdg` is used:** keep `IOM_S` as `Sg 1 H 1`. The consistent version then
  needs a fully S-decoupled database, with all reduced S (including kinetic
  Troilite) written in `HSg⁻`. Pyrite (S −1, written with SO₄²⁻) needs a
  rebuilt log K.

## Owed code fix

`iom_module.R` still instructs the caller to "enter H2 as Hdg" (docstring and
`carbon_decoupling_note`, lines ~211–348). That instruction is a no-op in
production. Rewrite it once decision D3 (Part 10) is made.

---
# PART 4 — IOM CALIBRATION STATE (REVISED 2026-10-26)

**Superseded: the single-Ea values and A = 1e13 that this section previously
carried.** Full derivation in `docs/PATCH_F_vitrimat2018.md`.

## The pre-exponential factor was wrong

Burnham (2019) Table 1 gives **A = 2 × 10^15 s^-1** for Vitrimat 2018.
Our assumed 1e13 was the **1989** value, which Burnham explicitly supersedes;
the paper notes the larger A improves agreement with hydrous pyrolysis, which
is exactly Miller's experiment type. `logA = 15.301`.

Refitting the same 350 C data at the correct A raises every Ea by
**+27.4 kJ/mol**. Nothing previously validated breaks (A and Ea trade off),
but any channel left at an old Ea with the new logA would be wrong twice.

## Use the Vitrimat 2018 distributions, not single Ea

Burnham 2019 Table 1 supplies full Ea distributions for H2O, CO2, CHn **and
CH4** (kcal/mol, percent of precursors):

| Ea | H2O | CO2 | CHn | CH4 |
|---|---|---|---|---|
| 42 | 10 | | | |
| 44 | 10 | 10 | | |
| 46 | 15 | 15 | | |
| 48 | 15 | 15 | | |
| 50 | 15 | 15 | 5 | |
| 52 | 15 | 15 | 15 | 2 |
| 54 | 10 | 15 | 30 | 5 |
| 56 | 10 | 15 | 20 | 8 |
| 58 | | | 15 | 10 |
| 60 | | | 10 | 12 |
| 62 | | | 5 | 15 |
| 64-76 | | | | 12,10,8,6,5,4,3 |

### THE KEY VALIDATION RESULT

Vitrimat 2018's CO2 distribution, as published with no adjustment, predicts
**45.81%** conversion in 48 h at 250 C. Miller's HC113 measured **45.39%**
(mean of four runs, relative to its 500 C point). **Agreement 0.991.**

A distribution Burnham fitted to *vitrinite coal* reproduces Miller's *IOM*
CO2 release split to within 1%. Our single-Ea model gets that same ratio wrong
by a factor of **1016**. Two unrelated labs, two unrelated materials — this is
the strongest independent validation the organic module has had.

### A real tension: the shift differs by sample

| sample | fit | result |
|---|---|---|
| Murchison (350+500 C) | shift **+3.6** kcal/mol, m0 3.431 | 350 C 3.038 (meas 3.04); 500 C 3.431 (meas 3.43); RMS 0.002 |
| HC113 (250+500 C) | shift **0.0**, m0 2.528 | 250 C 1.166 (meas 1.158); 500 C 2.546 (meas 2.55) |

Each fits its own data to under 1%, and they disagree. Mechanically: at shift 0
both 350 and 500 C are essentially fully converted, so Murchison's observed
350/500 difference cannot be reproduced. This is most likely **real material
difference** (HC113 synthetic, H/C 1.13, O/C 0.27; Murchison meteoritic, H/C
0.70, O/C 0.18), corroborated by Miller's Table 10 independently needing +8
kcal/mol for HC113. **Do not fit one shift to both.** Use Murchison for
production, HC113 as upper bound.

Cost of that spread at Enceladus temperatures (CO2 fraction over 4.5 Gyr):

    T        Murchison(+3.6)   HC113(0.0)   spread
    0 C       0.0000            0.0002        759x
    25 C      0.0004            0.0876        229x
    50 C      0.0960            0.3501        3.6x
    75 C      0.3623            0.6269        1.7x
    100 C     0.6389            0.9042        1.4x

Below ~50 C the material choice dominates everything; above ~75 C it barely
matters.

### And what single-Ea was getting wrong

Against the Murchison distribution fit, over 4.5 Gyr: single-Ea is low by
~3e7x at 0 C, ~6e6x at 25 C, ~8e5x at 50 C, 29x at 100 C.

## Current channel set (five channels)

| channel | m0_per_kg | formula | basis |
|---|---|---|---|
| IOM_CO2 | 4.1691 | C 1 H 1 O 1.3115 | fitted, Murchison, distributed Ea |
| IOM_CH4 | 2.2718 | C 1 H 4.435 | fitted, Murchison, distributed Ea |
| IOM_CHn | 1.2122 | C 1 H 1.8 | **transferred from vitrinite** (Burnham c(oil)=2% of C, n=1.8) |
| IOM_N | 0.8750 | N 1 H 2 | fitted, Murchison, single Ea **242.8** kJ/mol (was 244.8 — basis bug, see below) |
| IOM_S | 0.4620 | S 1 H 1 | **UNCALIBRATED** — pinned to Ea_N, no data exists |

Element budget on the ENTICES basis — verified, all residues positive:

    element  released    bulk   residue  %released
    C           7.653  60.611    52.958      12.6%
    H          18.638  45.639    27.001      40.8%
    O           5.468   9.939     4.471      55.0%
    N           0.875   2.142     1.267      40.8%
    S           0.462   1.154     0.692      40.0%

Implied inert residue H/C 0.510, O/C 0.0844 (bulk 0.753 / 0.164).

**Epistemic ordering, worth keeping straight:** CO2/CH4/N are fitted to
measured Murchison yields; CHn is transferred from vitrinite with no
IOM-specific support; IOM_S has none at all. Flag the last two in any output
that depends on them.

### Ea_N corrected 242.8 (was 244.8) — a basis bug, found 2026-09-20

The 244.8 value in `docs/PATCH_F_vitrimat2018.md` and `docs/archive/NEXT_STEPS.md`
mixed two composition bases. `m0_per_kg = 0.875` is on the ENTICES-normalised
basis (bulk N 2.142 mol/kg), but the fit target used was Miller's **raw**
Murchison measurement (real Murchison N 1.7705 mol/kg). Those differ by
2.142/1.7705 = 1.2098, so computing `target = raw_350C / m0_per_kg` silently
divided one basis by the other and Ea absorbed a ~21% conversion factor that has
nothing to do with activation energy.

The dimensionally clean target is the **ratio of Miller's two raw
measurements**, 0.5676/0.7229 = 0.7852, giving **Ea_N = 242.8 kJ/mol**.
Independently verified: 242.8 reproduces Miller's raw 350 C value exactly
(0.5676 predicted vs 0.5676 measured), whereas 244.8 predicts 0.4689 — **17.4%
low**.

Two related points, checked:
- **CO2/CH4/CHn do NOT share this bug.** Their shifts are fitted to the ratio of
  two raw measurements (basis-independent), with the ENTICES ×1.215 scaling
  applied afterwards as a separate step.
- **Ea_N independently re-derived 2026-09-30:** 242.78 kJ/mol. The old 244.8
  is exactly what the mixed-basis target gives. With A = 2e15, f(500 C, 48 h)
  = 1.0000, so this is effectively a **one-point fit at 350 C**. m0 = 0.875
  (41% of IOM N, matching Miller's 40.8% NH₄⁺ at 500 C) assumes complete
  release by 500 C, so it is a **floor** on the reactive N pool.
- **`m0_N = 0.875` is fine and should not be "fixed."** The basis-converted 500 C
  exhaustion value is 0.7229 × 1.2098 = 0.8746, which agrees with 0.875 to
  0.06%. It arrived by a different route (the ~40%-of-N release figure) but is
  self-consistent with how Ea_N is now fitted.

This also means the "~20-22% overshoot at 500 C because m0 is a fixed exhaustion
pool" explanation, repeated across several documents, had the right **sign** but
the wrong **mechanism** — it was this basis mismatch.

## Bulk composition — RESOLVED 2026-09-20: CR-chondrite IOM (Miller et al. 2026)

**The source.** `iom_bulk_formula()` (C 60.611, H 45.639, O 9.939, N 2.142,
S 1.154 mol/kg) is **Miller et al. 2026 SI Table 9** (728 C, 46 H, 159 O,
30 N, 37 S mg/g). It matches to 1 part in 10⁴, and closure is exactly
1000 g. Table 9 is the mean IOM of Alexander et al. 2007's seven CR
chondrites, chosen because CR IOM is the least altered.

**Two earlier labels were wrong.** It was called "kerogen-like" (09-20,
withdrawn) and "CI-derived" (09-30 handover, withdrawn).
- `reconstruction/`'s v0-era Rmd target (N 3.284, O 9.941) is Table 9
  converted with two slips: 46/14.007 and 159/15.9949. It is not a source.

**Two papers.** Composition comes from Miller 2026 (Titan building blocks).
Kinetics come from Miller 2025 GCA (Murchison/HC113 pyrolysis). So composition
is CR and kinetics are CM.
- That gap is smaller than the inter-lab spread on Murchison itself: H/C 0.699
  (Miller) vs 0.588 (Alexander 2007).
- The ENTICES↔Murchison C scale factor (1.2154) is the ash normalisation,
  unaffected by CR vs CM.

**Open (D10).** Table 9's S/C = 0.019 can't be reproduced from Alexander
2007's CR rows (mean 0.036, estimated by difference) or Alexander 2021
(0.021–0.090, measured directly).

## Nitrogen — RESOLVED (was open for several sessions)

Cody et al. (2024) via Miller Table 1 is the sole source: Murchison IOM
N = 2.48 ± 0.04 wt%, N/C = 0.036. On the ENTICES normalised basis (C 60.611)
that is **N = 2.182**; our value of **2.142** is within 1.8%, inside Cody's
own uncertainty. So 2.142 is correct.

The old 3.284 figure came from **misreading the hydrogen gram value as
nitrogen**: ENTICES H = 45.639 mol × 1.0079 = 46.0 g, and 46/14.007 = 3.284.
Not a competing measurement — a transcription error. Optionally refine 2.142 →
2.182 for exact consistency; either way delete the "UNRESOLVED" caveat from
`iom_bulk_formula()`.

## Validation state

- **250 C CO2 ratio (HC113): 0.991** — independent, passes. Add as a self-test
  assertion.
- **Murchison 350 C residue — REVISED 2026-09-30.** Against Murchison's
  measured residue, the model's absolute overshoot (H/C 0.679 vs 0.63) is
  mostly the CR-vs-CM starting difference. But `RX/C` ratios are **not**
  invariant to the start when a large fraction of an element is released.
  At 350 C/48 h the model releases 44–49% of the O.

  So the right test is **like-for-like**, from a Murchison start on Table 5's
  basis (H/C 0.700, O/C 0.181):

  | Quantity | Model | Measured | Error |
  |---|---|---|---|
  | RH/C | 0.888 | 0.90 ± 0.12 | −1.4% |
  | Residue O/C | 0.109 | 0.096 | +14% |

  - **H passes.**
  - **O is under-released by 13%.** If Murchison's start is really 0.203
    (Miller Table 1's own wt% give that; the table *states* 0.18), the residue
    O/C is 0.133 vs 0.096 and the under-release is **30%**.
  - **Record:** H passes; O under-released by 13–30%, depending on Murchison's
    uncertain starting O/C. Plausibly related to water supplying CO₂ oxygen in
    the experiments. Not tuned.
- **Murchison 500 C** H/C and O/C are confirmed genuinely blank in Miller's
  Table 5 (that row reports only N/C 0.015, RN/C 0.42, C 55.9 wt%, N 0.99
  wt%). Stop looking for them.
- Full Table 5 including syn-IOM 250 C and HC113 350 C (`350-10-Mar2023`,
  10 kbar, H/C 0.66, RH/C 0.58, O/C 0.047, RO/C 0.18) is in the Claude Project
  (not in this repo) as `Miller 2025 - Table 5.txt`. Compare using the dimensionless `RX/C` columns.

---
# PART 5 — ROCK AND PRIMARY MINERALS (NEW)

## Composition

**Source.** `reconstruction/` (moved into the repo 10-21), Lodders et al.
(2025) Table 4, Scenario R. The reconstruction Rmd writes
`reconstruction/rock_mol_per_kg.csv`, from which `ROCK_MOL_PER_KG` is
generated. `tests/rock_table_check.R` checks the two agree (1e-4; worst
1.1e-5).

**Phases** (mol per unit `total_rock`, i.e. per kg anhydrous rock, before the
organic scaling):

| Phase | mol per unit rock | Treatment |
|---|---|---|
| Enstatite | 2.2757 | kinetic |
| Forsterite | 1.6070 | kinetic |
| Troilite | 2.3830 | kinetic |
| Albite | 0.3182 | kinetic |
| Diopside | 0.2745 | kinetic |
| Anorthite | 0.06214 | kinetic |
| K-feldspar | 0.02052 | kinetic |
| Tephroite | 0.02599 | kinetic |
| Fe metal | 2.3504 | `dissolve_only` equilibrium; the reconstruction's 2.3616 minus lawrencite's Fe |
| Ni metal | 0.2809 | `dissolve_only` equilibrium |
| Lawrencite | 0.01125 | `dissolve_only` equilibrium |
| Schreibersite | 0.04709 | kinetic, `-formula Fe 3 P 1`, placeholder bulk rate (reconstruction §4/§8 1b) |
| Magnesiochromite | 0.0371 (CSV only) | inert: stays in the mass basis, so the other phases aren't renormalised; not passed to PHREEQC |

Total Fe matches the reconstruction exactly (4.886 mmol/g).

**Intentional choices (PI):**
- Cl enters as lawrencite in the rock, not in the fluid.
- `organic_wt_percent = 15` is deliberate and will be swept, to test primordial
  organic inventories beyond carbonaceous chondrites. CI is about 3 wt% of
  bulk, or about 4 wt% on the anhydrous basis; use that as the low end of the
  sweep.

## Metal corrosion

Fe and Ni metal are `dissolve_only` equilibrium phases, so they corrode
**completely at t = 0** of the first orbit. This was measured on 10-21 with
Magnetite in the porewater list:
- the products are 0.666 mol magnetite and 2.66 mol H₂;
- roughly 30% of the porewater is consumed;
- pH is about 4.5 at t = 0, because lawrencite's Fe goes into magnetite. It's
  buffered within seconds.

That makes the H₂ history a single pulse at the start. No metal corrosion rate
has been sourced (reconstruction §8 item 11). See decision D6.

## Rate laws

**PK04 for the silicates, with the `SR > 1` no-precipitation guard on every
primary** (schreibersite has no log K).

**Tephroite** follows Pincus et al. 2026 / Casey et al. 1993:
- **Acid branch:** `10^(-4.46)·ACT(H+)^0.47`, evaluated at
  `max(ACT("H+"), 1e-7)`. That gives a floor equal to the pH-7 rate above pH 7,
  with no summed neutral term.
- **Ea is pH-dependent:** `Ea = −60.11·ln(pH) + 106.17` kJ/mol, floored at 0.
  - It reproduces Casey's quoted temperature sensitivity: 2.3×, 1.5× and 1.1×
    per 10 K at pH 2, 3.5 and 5.
  - It implies **no temperature dependence above pH 5.85**, i.e. at all
    Enceladus conditions.
  - It is low stakes (0.5 wt%), but it's an extrapolation from 25–45 °C data.

**Surface areas are unsourced** (reconstruction §8 item 5): 0.1 m²/g for
forsterite, enstatite, diopside and tephroite; 5 m²/g for the feldspars,
anorthite and troilite. That 50× split controls the relative rates. See
decision D7.

## `rec_eqb_minerals` (optional, default off)

**What it does.** It promotes a primary to `dissolve_only` equilibrium if
2 × t½ < 1 orbit. The reference pH is 9 for base-branch minerals (albite,
K-feldspar) and 11 for acid-branch ones. At 160 °C it promotes albite
(t½ ≈ 2,000 s) but not K-feldspar (2.2×10⁵ s at pH 9).

**What it isn't.** It isn't needed for convergence. Lucas's TRUE/FALSE pair
agrees within ~5% on bulk chemistry. Minor phases differ by about as much as
the PHREEQC version does. The speed-up was 1.7× on Lucas's machine and
slower on ours.

---

# PART 6 — DATABASES AND NUMERICS (NEW)

## Databases (`databases/`)

**`coreclath_CH4coupled.dat`** (≥150 °C by default) and
**`coreclath_CH4uncoupled.dat`** (<150 °C). Chosen by `resolve_ch4_mode()`
in every generator, including the equilibrium ones since `eq-parity`.

**Shared by both:**
- the clathrates;
- fluorine, added 10-23 and identical to the old Core11;
- the five decoupled pseudo-elements, restored in the uncoupled file after
  `2e13e34` deleted them;
- ideal-gas treatment for all gases. The decoupled gases' critical constants
  are commented out, because PHREEQC can't mix ideal and Peng–Robinson gases
  in one GAS_PHASE.

**What differs.** The uncoupled file has no `C(-4)`, coupled CH₄ or `CH4(g)`,
and has `Mtg_hydrate` (Part 3).

**Other database changes:**
- `Core11_idealgas_mod_v4.dat` is retired to `legacy/`; tests A4/B2
  reproduce `docs/PHREEQC_TEST_RESULTS_2026-09-20.md` on the coupled file.
- Dead `CO2 = CO2` commented out.
- Duplicate CH₄(g)/C₂H₆(g) definitions removed. Both changes were checked as
  byte-identical in output.
- `tests/db_load_check.R`: every database must load with 0 errors. It would
  have caught `2e13e34`.

## Pressure

**Aqueous pressure dependence is off.** Aqueous `-Vm` is inactive under the
LLNL activity model these databases use (measured 09-30: speciation identical
from 1 to 1000 atm). Mineral log K and gases *are* pressure-corrected, so
SI shifts of ~0.05–0.1 at production pressures are partly an aqueous/mineral
inconsistency. For carbonate at 261 atm the missing aqueous term is
~0.11 log units.

**Generators use different porewater pressures:**
- `bulk_equil` porewater mode reacts at the profile's 261.2 atm;
- `tidal_kinetic`, `bulk_kinetic` and the restarts use the ocean's 66.7 atm.

Pure water at 300 °C fails at 66.7 atm, which is below its saturation
pressure. See decision D2.

## Convergence

**Cause.** The 10-06 hang was **tephroite precipitating**, because its rate
law lacked the `SR > 1` guard. Its SI reaches 0 at pH 8.8, where the run
stalled. Feldspars and gases were not the cause.

**After the fix (2 orbits):**

| Run | Wall time |
|---|---|
| `tidal_kinetic`, 100 °C | 660 s |
| `tidal_kinetic`, 160 °C | 302 s |
| `bulk_kinetic`, 25–160 °C | 95–290 s |
| equilibrium generators | under 2 s |

**At 25 °C,** `tidal_kinetic` needs the trace C seed (Part 8).

## PHREEQC version

**Standard: 3.8.6.** Every generated `.pqi` is stamped with the target version
and the database, and `run_phreeqc()` warns on a mismatch. Rerunning
Lucas's 3.5.0 outputs on 3.8.6 moved dissolved Mn by 130× and sulfide by
~11×. Bulk chemistry was stable.

## Legacy

**`legacy/ENTICES_v3.R`** is unreachable from `main`. Its `REACTION` block
wrote absolute element moles as coefficients *and* used the C moles as the
reaction amount. Delivered moles were therefore the intended amount ×
`iom_C_moles`:
- ×0.02809 in the redox-diagnostic case;
- over-dosed whenever `iom_C_moles` > 1.

All organic results from v3 are mis-dosed.

**`docs/REDOX_BUDGET_DIAGNOSTIC_2026-09-20.md`** ran on v3, with the dosing bug
and the superseded rock (fayalite/magnetite, no troilite, no metal). Its only
surviving conclusion is that old output was contaminated.
- Carbon closed exactly at 7.89e-4 mol.
- CH₄ : DIC ≈ 1.16 : 1, consistent with elemental C(0) disproportionating.
- The pH drop was carbonic acid, so it's one causal chain, and a pH-buffered
  control would have been the wrong test.

---

# PART 7 — GAS PHASES AND H₂ (NEW)

## The artefact

Every generator lists 11 gases at `0 0` in EQUILIBRIUM_PHASES. Per the
manual (Example 7) that is an **infinite reservoir at 1 atm partial pressure**,
whatever the solution pressure.

**The gas never reaches the ocean:**
- the REACT step saves EP1 with the porewater;
- the MIX steps use `equilibrium_phases none`.

**What that does to H₂:**
- H₂ exsolved in the porewater stays there, about 3 mol per unit rock after
  t = 0.
- It caps porewater dissolved H₂ at about 1×10⁻³ molal (the 1 atm
  solubility).
- **Ocean H₂ is set by that cap.** Measured 10-26 on two blanks with different
  rock amounts: the ratio of their ocean H₂ gains per orbit (0.839, 0.802)
  equals the ratio of porewater m(H₂) × exchanged water mass. So more rock
  gave *less* ocean H₂: its extra H₂ went into the trapped gas, and its
  porewater had less water.

**The ocean side (EP2)** has the same construction and is also saved. So
ocean exsolution isn't the "plume loss" the two-box design intended: the gas
can re-dissolve, and 1 atm isn't the vent condition (~0.006 atm).

**Consequence.** **No ocean H₂ number from any run so far is meaningful**,
including paired differences. CH₄ and CO₂ are also listed at `0 0`, so DIC and
CH₄ are affected whenever their fugacities reach 1 atm. That hasn't been
checked per run.

## The fix, tested but not adopted

A fixed-pressure GAS_PHASE at the porewater pressure: a bubble forms only when
Σpᵢ exceeds P, and its amount is finite. Tested 10-21 at 160 °C, 70 atm,
bubbles left in place:

| Quantity | Result |
|---|---|
| Dissolved H₂ | 7.3×10⁻² molal (×70) |
| Metal corrosion at t = 0 | still 100%, to magnetite |
| Bubble | 1.35 L, about 12× the porewater volume |
| Full 2-orbit run time | 6× slower (267 s vs 44 s) |

It is technically unblocked: all gases are now ideal. It still needs decision
D1 on pressure, bubble fate and ocean-side loss.

---

# PART 8 — INITIAL FLUID, SEED, ATTRIBUTION (NEW)

**Pure water stays the default, on purpose.** Then all ocean carbon can be
attributed to the organics, which is the test of whether IOM degradation can
produce Cassini's carbonate.

**A second scenario:** an ocean seeded with cometary ices (CO₂ etc.) before
circulation and maturation get going. The `initial_ocean` argument is
implemented; its values haven't been chosen (D8).

**Trace C seed (`trace_C_molal`, default NULL).** Pure water at 25 °C fails on
~10⁻²¹ mol of IOM carbon in the ocean step. A seed of **1e-14 mol/kgw**
(recommended in the README) makes `tidal_kinetic` and its restart complete.
- The seed is stamped in the header and carried down restart chains.
- At 25 °C over 2 orbits, the seed is all the carbon there is. That's the
  expected physics: the IOM releases ~7×10⁻²¹ mol.

**Attribution: `pair_organic_run()`.** It runs the same case with the organics
absent. The blank uses `mineral_basis_wt_percent` so that the **rock is
identical**. Without that, the blank has 1/0.85 the rock and 18% more metal;
that wrongly gave +2.7e-9 mol/kgw of "organic" H₂ at 160 °C. The seed cancels
in the difference.

160 °C, 2 orbits, organics' contribution to the ocean:

| Ocean quantity | Contribution (mol/kgw) |
|---|---|
| DIC | 6.7e-12 |
| CH₄ | 5.0e-9 |
| NH₄⁺ | 3.4e-14 |
| H₂ | −1.6e-11 (≈ 0; meaningless until Part 7 is fixed) |

**Attributing with a CO₂-bearing initial ocean** also uses the paired
difference. Simple bookkeeping fails once carbon enters carbonates or gas.

---

# PART 9 — CODE STATE (NEW)

**Branch state.** `main` @ `255981f` (PRs #1–#9 merged: IOM module, Lucas
integration, rock composition, fluorine, cleanup, fixes, consistency,
pairing, eq-parity). `eq-parity` (clathrates in `tidal_equil` porewater, and
`bulk_equil` CH₄ mode) merged as `255981f`.

**Layout:**

| Path | Contents |
|---|---|
| top level | `ENTICES_LF_20260917.R` (main program), `iom_module.R`, `mineral_dissolution_constants.r`, `phreeqc_tools.R`, `organic_contribution.R`, `parameter_sweep.R`, `ENTICES_runs.Rmd` (a sandbox run in chunks; depends on local plot scripts) |
| `databases/` | the two `.dat` files |
| `reconstruction/` | the rock module |
| `tests/` | the tests |
| `examples/` | example runs |
| `docs/` (`reviews/`, `archive/`) | design docs, review exchanges, superseded docs |
| `legacy/` | retired files |
| `scripts/` | standalone scripts |

**Tests:** `tests/iom_selftest.R`, `tests/rock_table_check.R` and
`tests/db_load_check.R`, plus PHREEQC test inputs A4/B2/B3. All pass.

**Generators:** `tidal_kinetic`, `tidal_equil`, `bulk_equil`, `bulk_kinetic`,
`bulk_kinetic_restart` and `tidal_kinetic_restart`. They are now consistent in:
- CH₄ mode and database selection;
- clathrate gating;
- Quartz/temperature gating;
- the Magnetite listing;
- one `Fe` column per CSV (metal is `Fe_metal`).

The equilibrium generators are labelled **inorganic limiting cases** (no IOM,
no schreibersite).

**Verification (`docs/reviews/STATE_OF_MAIN_2026-10-24.md`, plus later reruns):**
- All generators run with 0 errors at 25/100/160 °C, plus restarts. 25 °C
  `tidal_kinetic` needs the seed.
- Restart read-back is exact for all 38 kinetic reactants.
- Element mass balance closes to 1e-8–1e-16 over a 160 °C orbit.

**Still open in code:**
- the `iom_module.R` `Hdg` docstring (Part 3);
- renaming `ENTICES_LF_20260917.R` → `ENTICES.R` (Lucas's agreement);
- database consolidation into one source (Lucas's call).

**Production runs made before `fixes` must be regenerated.** The tephroite
rate, rock table and `rec_eqb` default have all changed since.

---

# PART 10 — DECISIONS AND NEXT ACTIONS

## Decisions (PI), in order of how much depends on them

| # | Decision | Evidence | Recommendation |
|---|---|---|---|
| D1 | **Gas handling**: fixed-pressure GAS_PHASE; what happens to the porewater bubble (stay, move to ocean, vent); whether ocean-side gas is removed as plume loss | Part 7 | **Gates every H₂ result.** Run the bracket: bubble stays vs bubble moves to the ocean. Remove ocean-side gas each orbit if it is meant as plume loss |
| D2 | **Porewater pressure**: profile (261 atm) or 66.7 atm, uniformly across generators | Part 6; 300 °C fails at 66.7 atm | Profile pressure everywhere; it sets the bubble threshold once D1 is in |
| D3 | **Carbon decoupling mechanism**: `Hdg` (no-op as specified), design (d), or the existing partial (d) in the uncoupled database | Part 3 | Run the leak test (N1) first; it may show the partial (d) needs only a few more species removed |
| D4 | `IOM_S` coupled or `Sg 1 H 1` | Part 3 | Follows D3: coupled if (d) |
| D5 | `IOM_CHn` routing and formate, if (d) | Part 3 | Decide after N1 |
| D6 | **Metal corrosion**: keep `dissolve_only` (H₂ pulse at t = 0) or add a rate | Part 5 | Needed before any H₂-vs-time claim; the rate is unsourced |
| D7 | **Surface areas** (0.1 vs 5 m²/g, unsourced) | Part 5 | Largest unsourced kinetic input; source them or run a sensitivity |
| D8 | **Cometary-ice initial ocean**: composition (CO₂/H₂O ratio, NH₃, from which comets) | Part 8 | — |
| D9 | **Spin-up** before tidal cycling; doubles as the King et al. forward check (reconstruction §8 item 8) | 10-06 review §4 | Worth doing as a modelling option; not needed for convergence |
| D10 | Chase where SI Table 9's S = 37 mg/g comes from. S/C 0.019 is below the CR mean of Alexander 2007 (0.036, by difference) and below Alexander 2021 (direct, 0.021–0.090) | `docs/CR_IOM_COMPOSITION_RESOLVED_2026-09-20.md` §3 | Stacks on the uncalibrated `IOM_S` rate |
| D11 | `rec_eqb_minerals` on or off | Part 5 | Off by default; a modelling choice |

## Next actions (code and tests), in order

- **N1. Leak test of `coreclath_CH4uncoupled.dat`.** At 25 and 100 °C, with
  coupled H₂ from metal corrosion, does coupled C reduce to C₂H₆, C₂H₄,
  acetate, formate or CO?
  - Report every C valence total other than C(+4).
  - Cheap, and it decides how far D3 is from done.
- **N2. B4: Option 3 on a real production input.**
  - (1) As-is.
  - (2) `IOM_CH4` → Mtg and `IOM_S` → Sg, with no metal rewrite. Expect
    carbon still to reduce wherever coupled reduced C exists.
  - (3) Design (d) with the leak check.
- **N3.** Rewrite the `iom_module.R` `Hdg` instructions per D3.
- **N4.** Implement D1/D2, then rerun the paired 160 °C case and the 25 °C
  seeded case.
- **N5.** Regenerate production runs on current `main` (PHREEQC 3.8.6, stamped).
- **N6.** Sensitivities: CHn m0 (Burnham's coal value of 2% of C), the
  Vitrimat Ea *width* (centres are fit to Murchison; widths are vitrinite's),
  surface areas (D7).
- **N7.** Spin-up and the King et al. comparison (D9).

## Still open (science), unchanged from 09-20

- **`IOM_S` rate:** uncalibrated; Miller measures no H₂S.
- **CHn:** m0 and shift transferred from vitrinite.
- **HC113:** Miller Table 10's +8 kcal/mol for HC113 vs our 0.0 is
  unreconciled. It doesn't block production.
- **H₂O channel:** omitted deliberately.
- **Kinetic shape for chondritic IOM.** Alexander 2007 says chondritic IOM
  follows a different maturation path from kerogens, with N in amines rather
  than pyridines. The Ea *centres* are fit to Murchison; the distribution
  *shape* is untested.

---

# PART 11 — PROCESS LESSONS (NEW)

1. **Findings that lived only in chat or the Project got lost:**
   - the dosing bug was never written down;
   - the CR composition provenance was missed by the code session;
   - RECONSTRUCTION_STATUS was invisible to it.

   Durable findings go in the repo. `docs/reviews/` now holds the exchanges.
2. **Executed checks beat reasoning, every time. The misdiagnoses that only
   execution caught:**
   - the convergence misdiagnosis (both sessions);
   - "only CH4_hydrate in the phase lists" (a regex that silently didn't
     match);
   - "kerogen-like" and then "CI-derived" composition.

   Test the ablation that removes each ingredient, including the one that
   breaks a pattern (the only primary without a guard).
3. **One source per number.** The rock table was hand-copied from a 2-decimal
   spreadsheet three times; it is now generated and tested. The same applies
   to databases: three diverging copies became two, plus a load test.
4. **Stamp provenance.** PHREEQC version, database, seed and initial ocean go
   in every `.pqi`, because a version change alone moved Mn 130×.
5. **Compare like with like.** The residue check mixed a CR start with
   Murchison data. The organic pairing mixed different rock masses. Both
   produced confident, wrong numbers.
