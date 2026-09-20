# PROGRESS — IOM module and redox treatment

Living document. Supersedes `PROGRESS_twobox_superseded.md` for everything
concerning organic matter.

**Last updated: 2026-09-18**, consolidating four findings that had been living
in separate patch documents (Vitrimat 2018 / A-value; the nitrogen resolution;
the CHn channel; the Option 3 database correction). Those documents remain as
the detailed derivations — this file is the current state.

Companion documents, in reading order:
- `claude/PATCH_F_vitrimat2018.md` — A = 2e15, distributed Ea, the 250 C
  validation result, sub-pool config rows
- `claude/PATCH_F_addendum_N_and_CHn.md` — nitrogen closed out, CHn channel
- `OPTION3_decision.md` — per-channel decoupling decision (supersedes Part 3's
  original sketch)
- `HANDOVER_option3_database_findings_2026-09-18.md` — the database reading
  that prompted the Option 3 revision

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

**However — see `OPTION3_decision.md` — this does not mean all four channels
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

**Action item, still the cheapest test available:** compare a run with
`organic_wt_percent > 0` against the same case with organics off, and look at
the *whole redox budget*, not just H2: dissolved H2, pe, and the Fe-bearing
mineral assemblage. The electron demand may show up as mineral over-oxidation
rather than — or as well as — H2 depletion, so checking H2 alone could miss it.

## Therefore: Option 3 moves from "interesting bracket" to "required"

Previously scoped as a nice-to-have bracket on the ocean's degree of redox
equilibrium. It is now the only way to model the low-temperature cases
correctly. The bracket interpretation still stands — but the decoupled run is
no longer the exotic end-member; at low temperature **it is the physically
appropriate one**, and the coupled run is the artefact.

---

# PART 3 — OPTION 3 DESIGN (REVISED 2026-09-18)

**The original four-species sketch in this section was wrong and has been
replaced.** Full reasoning in `OPTION3_decision.md`; summary here.

## What the real database has

`Core11_idealgas_mod_v4.dat` (note: **v4**, earlier drafts of this file said
v2) contains five genuinely decoupled pseudo-elements, under a section
labelled "Redox-uncoupled gases from phreeqc.dat":

    Hdg (H2)   Oxg (O2)   Mtg (CH4)   Sg (H2S, species H2Sg)   Ntg (N2)

There is **no `Amm`** (decoupled ammonium) and **no decoupled oxidised-carbon
species**. `Sg` has its own acid-base ladder added by MMD in 2024, so the
decoupled-sulfide mechanism is already-exercised prior work.

Critically: `C(-4)`, `N(-3)`, `S(-2)` are **redox states of the shared
elements** C, N, S — not decoupled species. In PHREEQC's
`SOLUTION_MASTER_SPECIES` the element name defines the electron-balance group,
so writing `-formula C(4) 1 ...` sets only the *starting* valence and the
solver still moves it down the shared ladder. That is precisely the coupled
behaviour Option 3 exists to prevent.

## Decision, per channel

The four channels are **not** in the same situation, so a uniform
"decouple everything" design was wrong from the start.

| channel | decision | decoupled formula | why |
|---|---|---|---|
| `IOM_CH4` | decouple | `Mtg 1 H 0.435` | `Mtg` exists and genuinely decouples. Correct as originally drafted. |
| `IOM_S` | decouple | `Sg 1 H -1` (pending test) | Name corrected from `S(-2)`. Sulfide/sulfate equilibration genuinely inhibited at low T. |
| `IOM_N` | **leave coupled** | *(none)* | Coupled equilibrium already gives the right answer — see below. |
| `IOM_CO2` | decouple **via `Hdg`** | *(no new species)* | Coupled gives exactly the wrong answer, but the fix is not a carbon species. |
| `IOM_CHn` | follows `IOM_CO2` | *(same mechanism)* | Same carbon problem, same solution. |

### Nitrogen stays coupled — this simplifies the design

`IOM_N` releases `N 1 H 2`, i.e. N at oxidation state **−2**. NH3/NH4+ is
N(−3); N2 is N(0). At Enceladus pore conditions (pe ≈ −11.5, pH ≈ 11.5,
H2-rich) **NH3/NH4+ is the thermodynamically stable nitrogen species** — N2
and NOx are favoured only under oxidising conditions. So coupled equilibrium
moves N from −2 to −3: one electron, in the same direction kinetics would take
it, to the same species. **There is no artefact to fix.**

This also disposes of `Ntg`: pinning released N as N2 would *oxidise* it from
−2 to 0, moving it away from the correct answer. `Ntg` is worse than leaving N
coupled. **Do not add a decoupled ammonium species** — database work for no
change in result.

*Caveat: this depends on the pore fluid staying strongly reducing. If an
oxidising case is ever modelled, revisit.*

### Carbon is the opposite case, and `Hdg` may be the answer

From this project's own output (`k1.00e-07_d194200_orb12`) at pe = −11.5:
**C(−4)/C(+4) = 3.1e-16 / 7.9e-33** — seventeen orders of magnitude, i.e. all
carbon to methane. `IOM_CO2` releases carbon at +1.6, so coupled equilibrium
drives a 5.6-electron reduction per mole. That is the Part 2b artefact, and
for carbon it is real and large.

Adding a decoupled oxidised-carbon element is unattractive: carbonate has
three protonation states (a substantial edit, not a one-liner like `Sg`), and
**decoupled carbon would no longer see calcite, dolomite or any carbonate
mineral** — carbonate precipitation would silently stop working for exactly
the carbon we care about. Routing `IOM_CO2` into `Mtg` is also wrong: it
imposes the same reduction by hand.

**Proposed alternative, using a species already present:** the mechanism
driving carbon to CH4 is the large H2 reservoir setting pe. `Hdg` (decoupled
H2) already exists. If H2 is decoupled it cannot act as the electron donor,
and pe is then set by the mineral assemblage instead. No new database species;
targets the actual mechanism; and it directly addresses Part 2b's electron
demand, which was being met by H2 plus mineral oxidation.

**This is a hypothesis, not a result.** It must be tested (see below). It may
prove too blunt — decoupling H2 affects *every* redox couple, including
mineral reactions we may want coupled. Fallback is the carbonate database edit
with the calcite caveat accepted and documented.

## Element conservation — the invariant, verified for both changed channels

    IOM_CH4  coupled  C 1 H 4.435     H total = 4.435
             decoupled Mtg 1 H 0.435  Mtg carries 4 H + 0.435 free = 4.435  OK
    IOM_S    coupled  S 1 H 1         H total = 1
             decoupled Sg 1 H -1      H2Sg carries 2 H − 1 free  = 1        OK

Negative coefficients are legal — confirmed against the PHREEQC manual text in
this project (KINETICS `-formula`: "a negative stoichiometric coefficient and
a positive value for reaction progress gives a negative mole transfer, which
removes reactants from the aqueous solution"), and its own `Fe_di_ox` example
uses −1.0.

Only `IOM_CH4` and `IOM_S` need a second formula column, so the config change
is smaller than originally assumed. But `redox_mode` must still be
**per-element**, because H2 decoupling is a database-level switch while
Mtg/Sg are formula-level.

## Two things that need a PHREEQC run, not more argument

**(a) Does `Sg 1` implicitly demand H2Sg's two hydrogens?** Determines whether
the right form is `Sg 1 H -1`, `Sg 1 H 1`, or `Sg 1` alone. Test: add a known
amount via each variant and check total H and S against hand arithmetic.

**(b) Does decoupling H2 (`Hdg`) actually stop the carbon runaway?** Test in
isolation, no ENTICES: reducing solution with known H2 inventory, add oxidised
carbon, run once with `H2` and once with `Hdg`, compare C(+4)/C(−4) and pe.
This is the decisive cheap test and should precede any module changes.

---

# PART 4 — CALIBRATION STATE (REVISED 2026-09-18)

**Superseded: the single-Ea values and A = 1e13 that this section previously
carried.** Full derivation in `claude/PATCH_F_vitrimat2018.md`.

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
| IOM_N | 0.8750 | N 1 H 2 | fitted, Murchison, single Ea 244.8 kJ/mol |
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
- **Murchison 350 C residue**: previous single-Ea prediction overshot at
  H/C 0.703 vs measured 0.63 (+11%) and O/C 0.104 vs 0.096 (+8%). Both high,
  attributed to the missing low-Ea tail. **This patch is the direct test of
  that hypothesis** — and because CHn is H-rich (H/C 1.8), adding it pulls
  predicted residue H/C *down*, the right direction. Re-run after implementing;
  it should now resolve or fail visibly.
- **Murchison 500 C** H/C and O/C are confirmed genuinely blank in Miller's
  Table 5 (that row reports only N/C 0.015, RN/C 0.42, C 55.9 wt%, N 0.99
  wt%). Stop looking for them.
- Full Table 5 including syn-IOM 250 C and HC113 350 C (`350-10-Mar2023`,
  10 kbar, H/C 0.66, RH/C 0.58, O/C 0.047, RO/C 0.18) is in the project as
  `Miller 2025 - Table 5.txt`. Compare using the dimensionless `RX/C` columns.

---

# PART 5 — STATUS AND NEXT ACTIONS

## Done
- IOM module split from ENTICES; first real execution completed.
- Test-harness regex bug found and fixed (module itself was correct).
- Ea corrected from a ratio fit to an absolute fit (18% error at 350 C).
- IOM_N recalibrated the same way.
- Part B converted to 16 asserted checks; passing.
- **A corrected to 2e15 from Burnham 2019** (+27.4 kJ/mol on every Ea).
- **Vitrimat 2018 distributions adopted**, validated independently at 250 C
  (0.991).
- **CHn channel added**; five-channel element balance verified.
- **Nitrogen discrepancy resolved** and its origin identified.
- **Option 3 design corrected** per-channel; N stays coupled, `Hdg` proposed
  for carbon.

## Next, in recommended order
1. **Redox budget diagnostic** (Part 2b). Cheapest test; may invalidate
   existing organic-enabled H2 and mineral predictions. Check H2, pe, *and*
   Fe phases.
2. **The two Option 3 PHREEQC tests** (Part 3): does `Sg 1` demand its own H,
   and does `Hdg` stop the carbon runaway. Both in isolation, no ENTICES.
   Test (b) decides the whole carbon-decoupling design, so it is high leverage
   for very little work.
3. **Implement distributed Ea** (Part 4 / `PATCH_F_vitrimat2018.md` §5):
   sub-pools for CO2, CH4, CHn; update logA to 15.301 everywhere; refit IOM_N
   Ea to 244.8. Re-derive the numbers rather than copying them.
4. **Re-run the Murchison 350 C residue check** — the direct test of the
   missing-low-Ea-tail hypothesis.
5. **Module support for decoupled mode** (second formula column for CH4/S
   only; per-element `redox_mode`), plus the element-conservation assertion.
6. Wire the module into ENTICES's template generator (`INTEGRATION POINTS` in
   `iom_module.R`).
7. Syn-IOM 250 C residue comparison via R-ratios, carrying the
   material-difference confound explicitly.

## Still open
- **IOM_S** uncalibrated; Miller measures no H2S. Will stay that way unless a
  source appears.
- **CHn m0 and shift** are vitrinite-derived with no IOM constraint. The
  2%-of-C figure is Burnham's coal value and could plausibly be several times
  larger for H-rich IOM — worth a sensitivity test.
- **H2O channel** deliberately omitted (products go to the solvent; no
  meaningful effect on ocean chemistry). Recorded so it is not mistaken for an
  oversight.
- **Miller Table 10's +8 kcal/mol for HC113 vs our 0.0 for HC113** —
  unreconciled. His is a residue fit at 10 C/Myr; ours a gas-yield fit at 48 h
  isothermal. Does not block production; our Murchison fit is independent.
- **Decoupling: global or temperature-dependent?** Global is cleaner and gives
  an honest bracket; a threshold is more physical but introduces an arbitrary
  discontinuity. Current recommendation: global switch, per-element
  configurable, let the bracket carry the uncertainty.
- **`iom_module.R` / `INTEGRATION_NOTE.md` drift.** The module copy available
  to the chat session predates the VS Code IOM_N recalibration, so patch
  specifications are written as instructions rather than replacement files.
  Whoever merges the `iom-module` branch should reconcile
  `INTEGRATION_NOTE.md`'s "what's NOT done" list against this file.