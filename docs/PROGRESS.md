# PROGRESS — IOM module, redox treatment, and the state of ENTICES

Living document: the current state. History lives in git and in `docs/reviews/`.

**Last updated: 2026-10-28** (`gas-ice-docs`: profile pressure, initial pe,
mineral volumes and rock density, the mineral `-Vm` strip, restart closure,
fixed-pressure gas phases, the cometary-ice option; Parts 6–10).
The 2026-10-27 revision folded in the review exchanges of
2026-09-30 → 2026-10-27 (all in `docs/reviews/`) and the merged branches
`fixes`, `consistency`, `pairing`, `eq-parity`, `docs-progress`,
`ch4-mode-iom`, `carbon-d` and `carbon-default` (merged as `c33a77f`).

**The carbon treatment is now decided** (Part 3). Coupled carbon is C(+4)
only, by default and at every temperature (`carbon_mode = "c4_only"`). Methane
exists only as decoupled `Mtg`, and H₂ stays coupled.

This revision **overturns** three conclusions of the 2026-09-20 version:

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

**Added 2026-10-26/27:**

| Finding | Effect | Part |
|---|---|---|
| Below 150 °C, `IOM_CH4` still entered as coupled `C 1 H 4.435`, though the uncoupled database has no coupled CH₄. So the IOM methane channel produced **no methane**; its carbon went to ethane, methanol or carbonate | Fixed in `ch4-mode-iom`: IOM formulas follow the database. All earlier sub-150 °C runs carry this error | 3 |
| PHREEQC's valence totals file formate, acetate, methanol, propane and others under C(4) | Leak checks must re-bin by actual oxidation state (`carbon_by_oxidation_state()`) | 6 |
| At realistic C (1e-3 mol/kgw) the uncoupled database leaks badly: 25 °C porewater carbon ≈100% ethane; ~1% reduced at 100 °C. C₂ species scale with [C]², so trace-C tests hide them | The partial design (d) is not a carbon treatment | 3 |
| `coreclath_C4only.dat`: every coupled non-C(+4) carbon entry removed. Leak check passes (100.000000% C(+4) in all blocks) | **D3 decided:** `c4_only` is the default at all temperatures | 3 |
| The 25 °C porewater becomes a ~2.5 molal NaOH brine (pH 14, I 2.4) as albite dissolves with no Na sink. It flushes over ~3×10⁴ orbits (~100 yr) | A transient; outside the activity model while it lasts | 5 |
| PHREEQC's *printed* SI tables are offset by 0.143 × the H₂O in each phase's reaction (25 °C brine). CSV SIs are correct | Don't read SIs from `.out` tables | 6 |

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

## Design (d): remove coupled reduced carbon, keep H₂ coupled — ADOPTED 2026-10-27

**The idea.** CO₂ can't be reduced if the coupled C element has no reduced
species. Carbonates are unaffected (C(+4) stays coupled), H₂ stays coupled,
sulfide stays stable, and methane exists only as decoupled `Mtg`.

### The partial version leaked (measured 10-26)

`coreclath_CH4uncoupled.dat`, the old default below 150 °C, had no coupled CH₄
but kept every other reduced carbon species. Leak test on generated
`tidal_kinetic` inputs with production metal corrosion, re-binned by *actual*
oxidation state:

| Total C | Case | Where the coupled carbon went |
|---|---|---|
| 1e-3 mol/kgw | 25 °C porewater | **≈100% ethane** (both orbits) |
| 1e-3 mol/kgw | 100 °C porewater | carbonate 99.0%; ethane 0.65%, acetate 0.16%, propane 0.15%, formate 0.07% |
| trace C | 25 °C | looked clean (99.95% carbonate in orbit 2) |
| either | ocean | ≥99.99% carbonate |

The trace-C result was an artefact. Ethane holds two carbons
(2 HCO₃⁻ → C₂H₆), and [C₂H₆]/[carbonate]² was identical at both carbon levels.
So **C₂ species grow as [C]²**, and trace-C tests hide them. That's the same
artefact as methane through a different species.

**Two traps this exposed:**
- **PHREEQC's valence totals can't be used for a leak check.** Formate,
  acetate, methanol, propane and others are filed under C(4). Use
  `carbon_by_oxidation_state()` (Part 6).
- **The IOM methane channel made no methane below 150 °C** (fixed in
  `ch4-mode-iom`, PR #11). The database was chosen by temperature, but the IOM
  config stayed coupled. `IOM_CH4` entered `C 1 H 4.435`, with no coupled CH₄
  to become, so ~30% of released IOM carbon was forced into ethane, methanol
  or carbonate. Now the IOM formulas follow the database (`resolve_iom_pools()`),
  and Mtg is punched and carried across restarts.

### The full version: `coreclath_C4only.dat` (PR #12)

**What's removed:** every coupled carbon entry whose carbon isn't +4:
- 70 species, including 45 acetate complexes, all organics and the
  cyanides;
- 16 phases: Glycine, Alanine, Citric_Acid, the C₂–C₃ and CO gases,
  graphite, 5 kerogens, pyridine, TiC, UC;
- the masters C(−3), C(−2) and C(+2).

SCN⁻ and urea are kept (carbon +4).

**Leak check passes:** coupled carbon is **100.000000% C(+4) in all 10
blocks** at 25 and 100 °C with 1e-3 C, and C + Mtg closes to ≤1.4e-5.

**Against the uncoupled database:**

| Quantity | Change under `c4_only` |
|---|---|
| Porewater | almost unchanged (H₂ is held by the trapped H₂(g)) |
| Ocean H₂ | 33–41× higher: reducing carbonate had been consuming it |
| 25 °C porewater carbon | goes briefly to dolomite instead of ethane |

**IOM under `c4_only`:**
- `IOM_CH4` → `Mtg 1 H 0.435`;
- `IOM_CHn` → `Mtg 1 H -2.2`, which imposes a 2.2 e⁻ reduction by hand
  (C −1.8 → −4), documented in `iom_module.R`;
- `IOM_N` and `IOM_S` are coupled.

**Formate sensitivity** (re-admitted, 100 °C, 1e-3 C): 0.07% of porewater C.
It still cuts ocean H₂ ~5×, through H₂ + HCO₃⁻ → HCOO⁻.

### Decision (D3–D5, 2026-10-27)

- **`carbon_mode = "c4_only"` is the default at every temperature.**
- **Coupled carbon** (`"coupled"`, or the old `"temperature"` rule with
  `ch4_cutoff_C`) is reported as an **upper bound** on abiotic methanation. At
  150 °C it gave CH₄/CO₂ ≈ 10⁹ in Lucas's run. Methanation is kinetically slow
  at 150–200 °C unless catalysed, and the rock's Fe–Ni metal is a possible
  catalyst.
- **Formate** is a sensitivity case.

## What this means for `IOM_S`

**Decided: coupled.** Under `c4_only`, H₂ stays coupled, pe stays low and
sulfide stays sulfide, so IOM sulfide keeps its mineral sink
(pyrite/troilite).
- The `Sg 1 H 1` formula stays in `iom_module.R` for any future
  `Hdg`/oxidising scenario.
- The S-decoupled-database route (all reduced S in `HSg⁻`, a rebuilt pyrite
  log K) is **not needed**.

## Owed code fix — done 10-27

`iom_module.R` no longer tells the caller to "enter H2 as Hdg". The docstring
and comments describe `c4_only`, keeping the `Hdg` tests as history, and
`carbon_decoupling_note` now says what each carbon mode does (N3, `docs-1027`).

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
  - Origin (checked 10-29 against Casey, Hochella & Westrich 1993, GCA 57,
    785–793): a least-squares fit of Ea against ln(pH) to Fig. 4's labels,
    54, 26, 22 and 4 kJ/mol at pH 2.5, 3.5, 4.2 and 5.15, with 0 taken at
    pH 6.2 (Fig. 4 says "< 4"). That reproduces −60.11 and 106.17 exactly.
  - Not from the appendix rate table
    (`Casey1993_Appendix_tephroite_dissolution_rates.csv`): an Arrhenius fit
    to its 25/35/45 °C replicate sets gives 58, 33, 47, 10 and −5 kJ/mol
    at those pH values, and −66.8·ln(pH) + 122.8.
  - The fit is not forced through (6.2, 0): it reaches 0 at pH 5.85.
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

## Na sink and the early porewater brine

**Na zeolites.** Natrolite and Mesolite were added to every generator's
secondary lists (`carbon-default`). At 100 °C they replace Analcime as the Na
sink (tidal_kinetic, 2 orbits: Natrolite 0.096 mol, Mesolite 0.058 mol).

**At 25 °C no Na phase is saturated** while albite dissolves; corrected SIs
are Saponite-Fe-Na −1.64 and Analcime −1.74. The Al–Si assemblage
(chamosite, greenalite, lizardite, tobermorite) holds Al and Si too low
*(inferred)*. So Na⁺ and OH⁻ accumulate. 10 orbits at 25 °C, `c4_only`:

| Orbit | Porewater I (molal) | pH | Porewater Na (molal) | Ocean Na (molal) | Ocean pH |
|---|---|---|---|---|---|
| 2 | 1.47 | 13.9 | 1.56 | 1.6e-5 | 10.0 |
| 3 | 2.20 | 14.1 | 2.52 (albite exhausted at orbit 4) | 4.0e-5 | 10.4 |
| 10 | 2.40 | 14.1 | 2.55 | 2.1e-4 | 11.1 |

**It's a transient.** Only 3.3×10⁻⁵ of the porewater is exchanged per orbit,
so it flushes over ~3×10⁴ orbits (~100 yr). Grouped mode compounds the
exchange (x_step = 1 − (1 − x_tide)^G), so long runs flush it. While it
lasts, ionic strength exceeds the B-dot model's ~1 molal range: treat that
window's activities and pH as qualitative.

**N9, closed 10-27.** The ocean's pH 11.1 with ~2×10⁻⁴ molal Na is what that
much NaOH gives at the ocean's 1.45 °C, where pKw = 14.87. The ~pH 10.3 once
quoted here used the 25 °C value, 14.0. About 2×10⁻⁴ molal OH⁻ gives
pH ≈ 14.87 − 3.75 ≈ 11.1. The ocean is charge-balanced (cations and anions
agree to ~1e-7 eq), and its alkalinity is ~91% OH⁻, then silicate and
phosphate. A pure-water ocean is unbuffered, so its pH is hypersensitive:
don't read it as an Enceladus pH.

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
- **Peng–Robinson for every gas** (since 10-28, `scripts/gas_eos.py`): each
  gas has `-T_c`, `-P_c`, `-Omega`; decoupled gases take their analogue's
  (Mtg←CH₄, Hdg←H₂, H2Sg←H₂S, Ntg←N₂, Oxg←O₂, as in phreeqc.dat 3.8.6), and
  H₂'s Ω is −0.225 (phreeqc.dat), not the +0.225 the entry had. PHREEQC uses
  PR only if every component of a gas phase has all three.
- **No mineral molar volumes** (since 10-28, `scripts/strip_mineral_vm.py`,
  reversible): see Pressure below. The volumes are kept in
  `data/mineral_molar_volumes.csv`.

**What differs.** The uncoupled file has no `C(-4)`, coupled CH₄ or `CH4(g)`,
and has `Mtg_hydrate` (Part 3).

**`coreclath_C4only.dat`** (the default since `carbon-default`) is the
uncoupled file with every coupled carbon entry that isn't C(+4) commented out
(`#c4only#`, listed in its header). That's 70 species, 16 phases, and the
masters C(−3), C(−2) and C(+2).
- **Kept:** SCN⁻ and urea, whose carbon is +4.
- **Regenerated by** `scripts/make_c4only.py`; `tests/c4only_db_check.R`
  checks the committed file matches it.
- **Generated inputs** drop the removed phase names automatically
  (`restrict_to_database()`).

**Carbon modes** (`carbon_mode` on all six generators):

| Mode | Database | Use |
|---|---|---|
| `"c4_only"` | `coreclath_C4only.dat` | **default at every temperature** |
| `"temperature"` | the old rule: coupled ≥ `ch4_cutoff_C` (150 °C), uncoupled below | reproduces old defaults byte for byte |
| `"coupled"` | `coreclath_CH4coupled.dat` | an **upper bound on abiotic methanation** (Fe–Ni metal could catalyse it at higher T) |
| `"ch4_uncoupled"` | `coreclath_CH4uncoupled.dat` | leaks (Part 3); diagnostic only |

**Restarts** read the parent's carbon mode from its `.pqi` header. They refuse
a mismatch unless `carbon_mode` is given explicitly, and then they warn,
because the parent's solutions would be re-speciated.

**Carbon accounting.** `carbon_by_oxidation_state()` (`phreeqc_tools.R`)
re-bins carbon by its actual oxidation state, counting `Mtg` as C(−4),
decoupled. PHREEQC's own valence totals file formate, acetate, methanol,
propane, amino acids and cyanides under C(4), so **never use them for a leak
check**. Test: `tests/carbon_rebin_check.R`.

**Printed SI tables.** In `.out` files they are offset by 0.143 × the H₂O in
each phase's reaction (measured on the 25 °C brine; present phases print
above 0). The CSV SIs are correct, and those are what ENTICES reads. The cause
hasn't been found; it may be worth a reproducible example to USGS.

**Other database changes:**
- `Core11_idealgas_mod_v4.dat` is retired to `legacy/`; tests A4/B2
  reproduce `docs/PHREEQC_TEST_RESULTS_2026-09-20.md` on the coupled file.
- Dead `CO2 = CO2` commented out.
- Duplicate CH₄(g)/C₂H₆(g) definitions removed. Both changes were checked as
  byte-identical in output.
- `tests/db_load_check.R`: every database must load with 0 errors. It would
  have caught `2e13e34`.

## Pressure (D2, implemented 10-27)

**From the profile, everywhere.** Every generator's porewater is at the
profile's mean pressure over the permeable layer (261.2 atm for the Enceladus
profile, now in `profiles/`), the ocean at the ocean floor (66.7 atm). No
fallback: a missing profile is an error. `REACTION_PRESSURE` pins each step's
pressure (MIX averages pressure by mass). Both are stamped in the header, and
restarts refuse a parent at different pressures (`allow_pressure_change`).

**Pressure acts only through gases.** Under the LLNL aqueous model PHREEQC
skips aqueous molar volumes, water density and dielectric pressure terms
(prep.cpp, utilities.cpp). With aqueous volumes zero, a mineral's reaction
volume was minus its whole molar volume, so every SI moved by
−Vm·ΔP/(ln10·RT) (66.7 → 261.2 atm, 160 °C: forsterite −0.10, albite −0.23,
antigorite −4.16), a physically meaningless shift. The mineral `-Vm` lines
were stripped 10-28 (database change for Lucas to confirm). Aqueous speciation
and mineral equilibria are now at the database's reference pressure (1 bar to
100 °C, the saturation curve above): with real dissolution volumes of tens of
cm³/mol, leaving out 260 bar costs ~0.1–0.2 log units. After the strip the
160 °C production case gives byte-identical output at 66.7 and 261.2 atm
(until the gas phase).

**300 °C.** Pure water at 300 °C failed at *every* pressure, not only below
saturation: the cause was the initial `pe 4` at the run temperature (~10³
molal O₂). Initial fluids are now defined at 25 °C with pe 4 and heated by
REACTION_TEMPERATURE (`tests/initial_fluid_check.R`; a full 300 °C orbit runs).
Under the old setting the initial porewater carried O(0) 1.1e-22 molal at
100 °C, 9.0e-12 at 160 °C and **4.9 at 250 °C** (0.77 mol O in 0.157 kg):
**old runs at ≥ 250 °C need a caveat**; below that it is negligible against the
~2.6 mol H₂ of the metal.

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

**At 25 °C,** `tidal_kinetic` needs the trace C seed (Part 8). So does a
cometary-ice ocean with a carbon-free porewater (Part 8).

**Gas phases (10-28).** `-high_precision` (needed for restart precision) also
sets PHREEQC's mole-balance criterion to 1e-12; with a gas phase the metal-pulse
step then ran > 15 min. A KNOBS after SELECTED_OUTPUT restores 1e-8 (the last
setting read wins, read.cpp; the manual says "the smaller of"). Gas-phase runs
are ~3.5× slower than before: Lucas's 12-orbit case at 150–160 °C takes
320–515 s (116 s with gases as equilibrium phases). `KNOBS -logfile` writes
every solver iteration: ~300 MB per 12-orbit run, and 11.9 GB from one stalled
run. **Off by default since 10-29** (PI); `logfile = TRUE` in any generator
turns it on.

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

## The artefact (before 10-28)

Every generator listed 11 gases at `0 0` in EQUILIBRIUM_PHASES. Per the
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

## As built (D1, 10-28)

**Fixed-pressure GAS_PHASEs, starting empty** (Parkhurst & Appelo 2013,
GAS_PHASE): porewater GAS_PHASE 1 at the porewater pressure, ocean GAS_PHASE 2
at the ocean pressure (never emptied: plume loss is unknown and may be
negligible); one GAS_PHASE 1 in the bulk generators. A bubble forms only when
Σpᵢ > P. Components are the gases the EP lists carried, filtered by carbon
mode (c4_only: Mtg(g) CO₂(g) H₂(g) H₂S(g) N₂(g) NH₃(g) SO₂(g); Ntg(g) with a
cometary ocean); clathrates stay equilibrium phases. All Peng–Robinson (Part 6).
H₂ at 160 °C, 261 atm: PR dissolves 8.3% more than ideal at saturation, and
its bubble's molar volume is 8.7% larger (`tests/gas_eos_check.R`).

**`pore_gas_fate`** (tidal generators):
- `"retain"`: the bubble stays with the porewater; MIX moves solution only. **When `Vgas_Vpore` > 1 this is an upper
  bound on retention, not a physical state.**
- `"to_ocean"` (**default**, PI 2026-10-29): after each porewater step the bubble's moles go to the ocean
  step, and the porewater starts the next orbit with an empty phase. PHREEQC
  has no keyword that adds one gas phase to another, so the porewater step's
  USER_PUNCH stores the component moles (`PUT`) and `KINETICS 3` adds them to
  the ocean solution at a constant rate over 1 s (exact: a constant rate
  integrates exactly); they then partition with GAS_PHASE 2 at the ocean
  pressure. Only one gas reservoir is ever in contact with a solution.

Why `"to_ocean"` is the default (PI, 2026-10-29): under `"retain"` the pores hold about 2.2 pore
volumes of gas, which is physically impossible, so `"retain"` is an upper bound on retention.
Carbon attribution is the same in both fates (organic − blank ocean C agrees to 3 significant
figures at 150 and 160 °C). `"to_ocean"` without plume loss saturates the ocean in H₂, so it is an
upper bound on ocean H₂. A third fate, `"retain_to_pore_volume"` (keep gas up to the pore volume,
send the excess to the ocean), is planned once the metal corrosion rate (D6) is decided, since the
bubble is metal H₂.

**`ocean_gas_loss_fraction_per_orbit`** was not added: removing a fraction of a
gas phase at run time would need the same PUT/GET route plus a negative
transfer, which is not cheap.

**Diagnostics** every row: `Vgas_L`, `P_gas_atm`, `Vgas_Vpore` (pore volume =
porewater mass at 1 kg/L), component moles (`g_` columns).

**Checks (10-28):** whole-system element totals per orbit (porewater, pore
gas, ocean, ocean gas, solids, IOM) close to ≤ 9e-12 relative in both fates
(`tests/gas_mass_balance_check.R`); with no bubble (25 °C, metal removed) the
two fates agree to the solver's tolerance; restarts carry both gas phases.

**The bracket.** Lucas's production case (k 1e-7, full-core layer, 12 orbits,
sub 1, 15 wt% IOM, seed 1e-14), with its same-rock blank, `c4_only`:

| T, fate | Pore bubble, orbit 1 → 12 (L) | Vgas/Vpore, orbit 12 | Ocean gas (L) | Ocean H₂ (molal) | Ocean pH | Mtg (molal) | ΣS(−2) | SO₄²⁻ | Organic − blank ocean C | Run / blank (s) |
|---|---|---|---|---|---|---|---|---|---|---|
| 150, retain | 0.325 → 0.348 | 2.21 | 0 | 1.41e-5 | 11.3 | 2.1e-14 | 6.1e-10 | ~1e-28 | 3.48e-8 | 412 / 412 |
| 150, to_ocean | 0.325 → 4.8e-4 | 0.0031 | 0.82 | 6.52e-2 | 11.3 | 4.7e-10 | 6.3e-10 | 0 | 3.48e-8 | 437 / 371 |
| 160, retain | 0.337 → 0.356 | 2.26 | 0 | 1.47e-5 | 11.3 | 9.9e-14 | 9.0e-10 | ~1e-28 | 1.29e-7 | 515 / 411 |
| 160, to_ocean | 0.337 → 3.3e-4 | 0.0021 | 0.82 | 6.52e-2 | 11.3 | 1.8e-9 | 9.4e-10 | 0 | 1.29e-7 | 421 / 320 |

Ocean values at orbit 12. DIC from the organics is the same in both fates;
ocean H₂ differs by 4 orders of magnitude, and organic Mtg by ~2×10⁴ (under
`to_ocean` it is also in the ocean gas: 6.2e-9 mol at 150 °C, 2.4e-8 at
160 °C). Before the gas phase (gases as equilibrium phases) the same case had
ocean H₂ 7.7e-8 molal at orbit 12.

Against the 10-21 test (160 °C, 70 atm, ideal): porewater H₂ is now 0.19–0.20
molal (saturation at 261 atm, PR), not 7.3e-2; the orbit-1 bubble 0.32–0.34 L,
2.1× the pore volume, not 1.35 L and 12×.

**The bubble is the metal.** Until the metal corrosion rate is decided (D6),
the bubble is the t = 0 pulse: Fe and Ni are `dissolve_only` equilibrium
phases, so they corrode completely in step 1 (~2.3 mol H₂ in the bubble). The
organic run and its blank have the same bubble to 5 digits.

---

# PART 8 — INITIAL FLUID, SEED, ATTRIBUTION (NEW)

**Pure water stays the default, on purpose.** Then all ocean carbon can be
attributed to the organics, which is the test of whether IOM degradation can
produce Cassini's carbonate.

**A second scenario (D8, implemented 10-28 as an option):**
`initial_ocean = cometary_ice_ocean(v = 0.02, include_NH3 = FALSE,
include_O2 = TRUE)`. Miller et al. (2026) Table 5: comet 67P ice (Rubin et al.
2019) with ices more volatile than NH₃ at 273 K (Fray & Schmitt 2009) depleted
by v_NH3 = 0.02, which comes from **Titan's** atmospheric N budget (Miller et
al. 2026, eqs. 21–24), not from Enceladus. The ices react into the pure-water
ocean at ocean T and P (REACTION, not SOLUTION totals: the mixture is out of
redox equilibrium). Result, `c4_only`:
- coupled C 0.0866 molal, all C(+4) (CO + H₂O → CO₂ + H₂ is forced: there is
  no C(+2)); Mtg 0.0102; pH 1.85 (sulfuric acid from oxidised sulfide);
  S(+6) 0.0094, S(−2) 0.0028 molal; H₂ ≈ 0; no bubble at 66.7 atm. The hand
  estimate (0.0086 / 0.0036) left out that `Mtg 2 H -2` (ethane) draws 2 e⁻
  per C₂H₆ from the coupled system;
- without O₂: no sulfate, H₂ 0.031 molal, pH 3.80; with NH₃: 7.44e-3 molal
  N(−3).

Caveats: at 1.45 °C both the CO shift and the O₂ oxidation are equilibrium
assumptions (water–gas shift is slow in cold water without a catalyst); CO is
40% of the carbon. **NH₃ is pending the PI**: Table 5 lists none and its Total
N excludes it, though the text applies v_NH3 to NH₃. v = 1 (undepleted 67P) is
an unrealistic upper bound (~2.6 molal CO₂ + 1.7 molal CO).

**Attribution with a cometary ocean.** `pair_organic_run()` passes the same
`initial_ocean` (and `pore_gas_fate`) to the run and the blank. The ocean
starts with ~0.087 molal DIC, so the organic contribution is organic − blank
and can't be read from the ocean's carbon. Comparison with Glein & Waite
(2020), GRL 47, e2019GL085885 (README §7, Carbon): they adopt Postberg et
al. (2009)'s envelopes, 0.2 molal Cl with 0.1 molal TDC (HCO₃⁻ + CO₃²⁻) or
0.05 molal Cl with 0.01 molal TDC, and prefer pH 8.5–9, log a(CO₂) −4.6 to
−3.2. The cometary option has 0.087 molal DIC, but at pH 1.85 it is 0.0866
molal CO₂(aq) and 2e-6 molal HCO₃⁻. The pure-water default's organic
contribution is 1.29e-7 molal after 12 orbits at 160 °C; the long runs (N8)
are what can say whether it approaches the observed range. `organic_contribution()` also reports the ocean gas phase (mol).

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

**Charge balance.** `charge_balance` (`"pH"`, an element listed in
`initial_ocean`, or NULL = unchanged) balances the ocean SOLUTION 2. Use it
for any cometary composition: an unbalanced 1e-3 C ocean is off by ~4e-4 eq.
At the production seed (1e-14) the imbalance is negligible, and the porewater
seed is left unbalanced.

**The pure-water ocean is unbuffered,** so its pH swings far on tiny inputs
(Part 5). That's acceptable for attribution, but it isn't an Enceladus pH.

---

# PART 9 — CODE STATE (NEW)

**Branch state.** `gas-ice-docs` (this revision) is open against `main` @
`f9a0eb8` (PR #14, `docs-1027`, merged). Before that, `main` @ `c33a77f` (PRs #1–#13 merged: IOM module, Lucas
integration, rock composition, fluorine, cleanup, fixes, consistency,
pairing, eq-parity, docs-progress, ch4-mode-iom, carbon-d, carbon-default).
`carbon-default` (`c4_only` default, restart mode check, `make_c4only.py`, the
Na zeolites) merged as `c33a77f`.

**Tests:** the five original self-tests (`iom_selftest.R`,
`rock_table_check.R`, `db_load_check.R`, `c4only_db_check.R`,
`carbon_rebin_check.R`) and, from `gas-ice-docs`: `pressure_check.R`,
`initial_fluid_check.R`, `bulk_density_check.R`, `restart_closure_check.R`,
`gas_eos_check.R`, `gas_mass_balance_check.R`, `cometary_ocean_check.R`
(see `tests/README.md`).

**Restarts (fixed 10-28).** Both restart paths dropped K (K-feldspar), Mn
(tephroite), Cl (lawrencite) and Ni at every restart (no generator punched
them) and carried everything else at 4–5 significant figures. Every generator
now punches all 36 database elements at high precision, and restarts carry
them at 15 digits; a 2-orbit run and 1 + 1 restarted agree on every element's
system total (`restart_closure_check.R`). **Ask Lucas which production runs
used `tidal_kinetic_restart` or `bulk_kinetic_restart`**: those lost K, Mn, Cl
and Ni at each restart.

**Rock density.** `bulk_density()` (`phreeqc_tools.R`) with
`data/mineral_molar_volumes.csv` (README, Interpreting outputs).

**Layout:**

| Path | Contents |
|---|---|
| top level | `ENTICES_LF_20260917.R` (main program), `iom_module.R`, `mineral_dissolution_constants.r`, `phreeqc_tools.R`, `organic_contribution.R`, `parameter_sweep.R`, `ENTICES_runs.Rmd` (a sandbox run in chunks; depends on local plot scripts) |
| `databases/` | the three `.dat` files |
| `reconstruction/` | the rock module |
| `tests/` | the tests |
| `examples/` | example runs |
| `docs/` (`reviews/`, `archive/`) | design docs, review exchanges, superseded docs |
| `legacy/` | retired files |
| `scripts/` | standalone scripts |

**Tests:** the five self-tests above, plus PHREEQC test inputs A4/B2/B3. All
pass.

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
- renaming `ENTICES_LF_20260917.R` → `ENTICES.R` (Lucas's agreement);
- database consolidation into one source (Lucas's call).

**Production runs made before `fixes` must be regenerated.** The tephroite
rate, rock table and `rec_eqb` default have all changed since.

---

# PART 10 — DECISIONS AND NEXT ACTIONS

## Decisions (PI), in order of how much depends on them

| # | Decision | Evidence | Recommendation |
|---|---|---|---|
| D1 | ~~Gas handling~~ | Part 7 | **Decided and implemented 10-28:** fixed-pressure GAS_PHASEs, ocean gas kept; `pore_gas_fate` bracket run (Part 7). **Default `"to_ocean"`** (PI, 10-29) |
| D2 | ~~Porewater pressure~~ | Part 6 | **Decided and implemented 10-27:** the profile everywhere |
| D3 | ~~Carbon decoupling mechanism~~ | Part 3 | **Decided 10-27:** `c4_only` at every temperature; coupled runs are the methanation upper bound |
| D4 | ~~`IOM_S` coupled or `Sg 1 H 1`~~ | Part 3 | **Decided:** coupled (H₂ stays coupled, so sulfide is stable) |
| D5 | ~~`IOM_CHn` routing and formate~~ | Part 3 | **Decided:** CHn as `Mtg 1 H -2.2` (a 2.2 e⁻ reduction imposed by hand, documented); formate as a sensitivity case only |
| D6 | **Metal corrosion**: keep `dissolve_only` (H₂ pulse at t = 0) or add a rate | Part 5 | Needed before any H₂-vs-time claim; the rate is unsourced |
| D7 | **Surface areas** (0.1 vs 5 m²/g, unsourced) | Part 5 | Largest unsourced kinetic input; source them or run a sensitivity |
| D8 | **Cometary-ice initial ocean** | Part 8 | **Implemented 10-28 as an option** (Miller et al. 2026 Table 5); **NH₃ pending the PI** |
| D9 | **Spin-up** before tidal cycling; doubles as the King et al. forward check (reconstruction §8 item 8) | 10-06 review §4 | Worth doing as a modelling option; not needed for convergence |
| D10 | Where SI Table 9's S = 37 mg/g comes from (S/C 0.019, below Alexander 2007's CR mean 0.036 and Alexander 2021's 0.021–0.090) | `docs/CR_IOM_COMPOSITION_RESOLVED_2026-09-20.md` §3; `docs/CR_IOM_TABLE9_SULFUR_ADDENDUM_2026-10-27.md` (the Project's §8 addendum) | **The PI is checking it:** a zero-filled average of the seven CR rows reproduces Table 9's sulfur (see the addendum) |
| D11 | `rec_eqb_minerals` on or off | Part 5 | Off by default; a modelling choice |

## Next actions (code and tests), in order

- ~~N1. Leak test~~ **Done 10-26** (Part 3).
- ~~N2. B4~~ **Superseded:** `Hdg` is dropped, and the leak check on
  production-style inputs replaces B4.
- ~~N3. Rewrite the `iom_module.R` `Hdg` instructions and the
  `carbon_decoupling_note` for `c4_only`.~~ **Done 10-27** (`docs-1027`).
- ~~N4. Implement D1/D2~~ **Done 10-28** (`gas-ice-docs`; bracket in Part 7).
- **N5.** Regenerate production runs on current `main` (PHREEQC 3.8.6,
  `c4_only`, stamped). **Every earlier run is affected** by at least one of:
  - the IOM_CH4 bug below 150 °C;
  - the coupled-carbon artefact above it;
  - the gas artefact.
- **N6.** Sensitivities:
  - CHn m0 (Burnham's coal value of 2% of C);
  - the Vitrimat Ea *width*;
  - surface areas (D7);
  - formate re-admitted (cuts ocean H₂ ~5× at 100 °C via H₂ + HCO₃⁻ →
    HCOO⁻);
  - `carbon_mode = "coupled"` as the methanation bound.
- **N7.** Spin-up and the King et al. comparison (D9).
- **N8.** Confirm in a long grouped run that the 25 °C porewater brine flushes
  as expected (Part 5).
- ~~N9. pH 11 with 2e-4 molal Na~~ **Closed 10-27** (Part 5: pKw at 1.45 °C).
- **N10.** D6 (metal corrosion rate): until it is decided, every gas-phase
  result is dominated by the t = 0 metal pulse (Part 7).
- **N12 (after D6).** A third `pore_gas_fate`, `"retain_to_pore_volume"`: keep
  gas up to the available pore volume and send the excess to the ocean. Build it
  once D6 is fixed, because the bubble is metal H₂.
- **N11.** Fluorite and Fluorapatite have no molar volume (not production
  phases). ~~The IOM density for `bulk_density()` is UNSOURCED~~: **sourced 10-29**,
  Delarue et al. (2026) Eq. 10 at the residual IOM's H/C (1.320 g/cm³ at t = 0;
  README §7).

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
