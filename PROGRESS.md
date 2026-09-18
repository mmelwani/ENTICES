# PROGRESS — IOM module and redox treatment

Living document. Supersedes `PROGRESS_twobox_superseded.md` for everything
concerning organic matter. `INTEGRATION_NOTE.md` in the `iom-module` branch
predates the IOM_N recalibration and the Part B assertion work and should be
reconciled against this file when that branch merges.

Last updated: 2026-09-17, after the first real execution of the module
(VS Code session) and the low-temperature speciation problem surfacing.

---

# PART 1 — WHAT THE IOM CHANNELS RELEASE (read this before touching anything)

This section exists because the answer is counter-intuitive and has been
re-asked more than once. It is the single most important thing to understand
about the module.

## The channels release ELEMENTS. Never molecules.

Every channel's `-formula` — for example `IOM_CO2`'s `C 1 H 1 O 1.3115` —
adds **bare elemental totals** to the PHREEQC solution. PHREEQC then decides
what aqueous species those elements become, by thermodynamic equilibrium at
whatever pH and pe the solution happens to have at that instant.

The channel *names* (`IOM_CO2`, `IOM_CH4`, `IOM_N`, `IOM_S`) are labels for
**which elements, in what ratio, released at what rate**. They do **not**
force the named molecule to appear in solution. `IOM_CO2` is not "the channel
that makes CO2"; it is "the channel calibrated against the population of
bonds whose breakdown Miller measured as CO2 yield."

## Why it is built this way

Two independent reasons, both of which have to hold:

**1. PHREEQC mechanically cannot do otherwise.** `-formula` in a KINETICS
block is an elemental-composition addition to solution totals. There is no
syntax by which it adds a molecule that persists as that molecule. Whatever
you write, the products are decided afterwards by the equilibrium solver.

**2. Even if it could, the product molecule would be the wrong thing to
add.** Checking Miller et al.'s own mass balance at 350C: the residue loses
4.50 mol O/kg, but the measured CO2 requires 6.08 mol O/kg. Roughly a quarter
of the oxygen in the product CO2 does not come from the kerogen at all — it
comes from the water, because these are hydrous pyrolysis experiments. So the
physically correct quantity to hand PHREEQC is **what the solid gives up**,
and to let the water contribution and the redox speciation emerge from
PHREEQC's own thermodynamics.

This is why formulas look like `C 1 H 1 O 1.3115` and not like `CO2`. The
fractional oxygen is the *solid's* oxygen loss per carbon, measured; the
balance needed to make CO2 (if CO2 is what forms) is drawn from water by
PHREEQC.

**Do not "simplify" these formulas to product molecules.** It is the most
likely well-intentioned change someone would make, and it would be wrong on
both counts above.

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

The CO2/CH4 redox couple is famously sluggish. Abiotic interconversion
between oxidised and reduced carbon requires high temperature, effective
mineral catalysis, or biological mediation. Below roughly 150-200C, CO2 and
CH4 do not equilibrate with each other on laboratory timescales and often not
on geological ones either — they persist in metastable coexistence. This is
standard in the hydrothermal organic geochemistry literature (McCollom,
Seewald, Shock and co-workers on metastable equilibria in hydrothermal
systems). *Caveat: citations here are from memory and cannot be verified from
the environment this was written in — confirm them before they go in a paper.
The underlying physical result is not in doubt; the specific attributions
should be checked.*

The same is true of the other redox-active elements the IOM releases:
- **N:** NH4+ / N2 / NO3- interconversion is strongly kinetically inhibited
  abiotically. (Stock `phreeqc.dat` ships a decoupled `Amm` species for
  exactly this reason.)
- **S:** HS- / SO4-2 equilibration is very slow at low temperature without
  microbial mediation.

So **all four** of our channels release elements into redox couples that are
kinetically frozen at Enceladus temperatures. Full redox coupling is not a
minor idealisation here; it is the wrong limit for the low-temperature cases,
which are a substantial part of the parameter space being modelled.

## Two consequences, one of which is quantitatively severe

**(a) Predicted carbon speciation is wrong.** The model predicts an ocean
with essentially no oxidised carbon. Cassini observes both CO2 and CH4 in the
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
carbon at about +1.6, because only ~1.3 oxygens come off per carbon. The
corrected figure is ~12 mol H2-equivalent.)

**Action item:** compare a run with `organic_wt_percent > 0` against the same
case with organics off, and look at the *whole redox budget*, not just H2:
dissolved H2, pe, and the Fe-bearing mineral assemblage (magnetite, Fe(III)
phases). The electron demand may show up as mineral over-oxidation rather
than — or as well as — H2 depletion, so checking H2 alone could miss it.
Whichever way it manifests, at low temperature it is an artefact.

## Therefore: Option 3 moves from "interesting bracket" to "required"

Previously scoped as a nice-to-have bracket on the ocean's degree of redox
equilibrium. It is now the only way to model the low-temperature cases
correctly. The bracket interpretation still stands and is still valuable —
but the decoupled run is no longer the exotic end-member; at low temperature
**it is the physically appropriate one**, and the coupled run is the artefact.

---

# PART 3 — OPTION 3 DESIGN (coupled / decoupled redox)

## The mechanism

PHREEQC decouples a redox couple when the valence states are defined as
**separate master species with distinct element names**, so they do not share
an electron balance. Stock `phreeqc.dat` does this for `Mtg` (methane) and
`Amm` (ammonium). We need the equivalent in the project database
(`Core11_idealgas_mod_v2.dat`), at minimum for methane, ammonium, and
sulfide.

## The design consequence that must be understood before implementing

Decoupling **changes the semantics of the channel formulas.**

- **Coupled mode:** formula = what the solid loses, as elements. PHREEQC
  decides products. (This is Part 1, and it is correct for that mode.)
- **Decoupled mode:** the model must now *specify* which decoupled species
  the released carbon / nitrogen / sulfur enters as, because nothing
  downstream will re-partition it. Product identity becomes an input, not an
  output.

So the config needs a **second formula column** — a decoupled-mode formula
per channel — rather than reusing one formula for both modes. Sketch:

    channel     coupled formula          decoupled formula (DRAFT, UNVERIFIED)
    IOM_CO2     C 1 H 1 O 1.3115         C(4) 1 H 1 O 1.3115
    IOM_CH4     C 1 H 4.435              Mtg 1 H 0.435
    IOM_N       N 1 H 2                  Amm 1
    IOM_S       S 1 H 1                  S(-2) 1 H 1

Note `IOM_CH4`'s decoupled form: `Mtg` already carries CH4's four hydrogens,
so only the surplus 0.435 H is released separately. **These drafts are
unverified.** PHREEQC's acceptance of valence-state notation like `C(4)` in a
KINETICS `-formula`, and the exact charge/electron bookkeeping when mixing a
decoupled pseudo-element with free H, both need testing before use.

**UPDATE 2026-09-18, checked against the real database, not stock
phreeqc.dat:** `HANDOVER_option3_database_findings_2026-09-18.md` found that
two of these four are wrong, not just unverified. `Core11_idealgas_mod_v4.dat`
has no `Amm` species at all (only `Ntg` = decoupled N2, the wrong end of the
ladder for IOM_N's reduced-nitrogen release), and `S(-2)` is a coupled-ladder
valence tag, not a decoupled species (the real one is `Sg` = H2S). `Mtg` is
confirmed correct and already in the database (added by the user in 2024).
Read that handover doc before touching this table or `iom_module.R`'s
`DECOUPLED_FORMULAS` -- it lays out two real design options for IOM_N/CO2
rather than a fix, since fixing it needs a decision this document can't make
on its own.

## The invariant that makes this testable

**Total elements released from the solid must be identical in both modes.**
Only the redox partitioning may differ. This is a strong, cheap check:

> Run the same case in coupled and decoupled mode. Sum C, H, O, N, S across
> solution + gas + precipitated phases. The totals must match to rounding.
> If they do not, the decoupled formulas are unbalanced.

Write that check before the database work, not after.

## Implementation order (each step verifiable before the next)

1. **Database first, in isolation.** Add the decoupled species to a copy of
   `Core11_idealgas_mod_v2.dat`. Test with a trivial hand-written PHREEQC
   input: put CO2 and CH4 into a strongly reducing solution and confirm they
   do **not** interconvert; confirm the coupled database *does* convert them.
   No ENTICES, no module. This isolates database errors from everything else.
2. **Module support.** Add the decoupled formula column and a `redox_mode`
   argument to `iom_default_config()`; make it **per-element configurable**
   (e.g. decouple C and N but not S) rather than a single global flag, since
   the three couples have different inhibition temperatures and we will want
   to test them separately.
3. **Element-conservation check** (the invariant above) as a new assertion in
   `iom_selftest.R`.
4. **Then** wire into ENTICES and run the dual cases.

## Interpretation, once it runs

The coupled/decoupled pair brackets how far the ocean can be from redox
equilibrium. Three outcomes, all publishable:
- Cassini ratios fall inside the bracket — the position within it constrains
  the effective degree of equilibration.
- Ratios sit at the decoupled end — kinetic delivery dominates; the ocean
  does not equilibrate on the circulation timescale.
- Ratios fall outside entirely — something else is needed (different source,
  plume fractionation, ice-shell processing).

The A–B gap is also, directly, the redox disequilibrium free energy available
to metabolism. Worth reserving output columns for that now even though the
calculation comes later.

---

# PART 4 — VALIDATION STATE AND THE RESIDUE DATA QUESTION

## Where calibration currently stands

| channel | m0_per_kg | Ea (kJ/mol) | basis |
|---|---|---|---|
| IOM_CO2 | 4.169 | 216.2 | absolute 350C yield, Murchison |
| IOM_CH4 | 1.908 | 229.4 | absolute 350C yield, Murchison |
| IOM_N | 0.875 | 217.3 | absolute 350C NH4+ yield, Murchison |
| IOM_S | 0.462 | 217.3 | **UNCALIBRATED** — pinned to Ea_N, no data |

All fits are **two-point, single-Ea, A fixed at 1e13**. Miller does not
measure H2S at all, so `IOM_S` has zero independent support; the number
exists only so the channel runs. Flag it in any sulfur-dependent output.

## The one genuinely independent validation point so far

Gas yields at 350C are circular — Ea was fit to them. The **residue
composition** was not used in the fit, so it is the only real test:

- 350C Murchison: predicted H/C 0.703 vs measured 0.63 (+11%);
  predicted O/C 0.104 vs measured 0.096 (+8%).

Both high, same direction. The model under-releases, leaving a residue too
rich in H and O. That is the expected signature of a **missing low-Ea tail** —
material that should have decomposed by 350C but hasn't, because a single Ea
tuned to the 350/500C endpoints cannot also capture early, low-temperature
release.

**This is an inference, not a verified diagnosis.** It is consistent, but a
+10% offset could also come from the m0 basis, from the residue's unreacted
channel mass being counted as solid, or from sample heterogeneity. Do not
repeat it as settled.

## Table 10 and the activation-energy distribution — TESTED, with a warning

Table 10 (transcribed; the 62 kcal entry was cut off in the paste and is
inferred as 5 — reasoning below):

    Ea (kcal/mol)   BS89 %    HC113 (Miller) %
    42               5
    44              15
    46              25
    48              25
    50              15         5
    52              10        15
    54               5        25
    56                        25
    58                        15
    60                        10
    62                         5   <- inferred

Both columns sum to 100 with that inference, and the two distributions have
**identical shape** (5,15,25,25,15,10,5) — Miller's HC113 tuning is simply the
BS89 distribution **shifted +8 kcal/mol**. That makes the inferred value
near-certain, but confirm against the paper.

### WARNING: Table 10 cannot be dropped into our framework as-is

Miller's Vitrimat curves are run at a **geologic heating rate of 10 C/Myr**
(Fig. 9 caption: "temperature increases from 20 C to greater than 600 C at a
rate of 10 C/Ma"), and the HC113 column is tuned to **HC113 synthetic IOM**,
not Murchison. Our Ea values are fitted to **48 h isothermal** experiments on
**Murchison**. Different materials, different time regimes.

Tested directly. Using our m0 = 4.169 with each distribution, at 48 h:

    model                        350C      500C    (measured: 3.04 / 3.43)
    BS89 distribution            3.840     4.169
    HC113 distribution (Miller)  0.836     4.168
    our single Ea (216.2 kJ)     3.034     4.169

Miller's HC113 distribution **under-predicts Murchison's 350 C yield by 3.6x**.
It is calibrated for a more refractory synthetic material on Myr timescales.
Do not adopt it directly.

Related: Fig. 9b shows the same distributions "all shifted 22 kcal/mol higher"
to demonstrate what would be needed to match the lab data — so the earlier
note about a +22 kcal/mol (+92 kJ/mol) shift was correct and refers to Fig. 9b,
not to Table 10's +8. Both are real; they answer different questions.

### What does work: fit the BS89 *shape*, not its absolute values

Keep BS89's distribution shape as a prior and fit only a uniform shift plus
m0 to our own Murchison data. Still a two-parameter fit — same cost as the
current single-Ea approach — but distributed:

    best fit: shift = +0.6 kcal/mol, m0 = 3.40 mol/kg (Murchison basis)
      350C: 3.043 (measured 3.04)
      500C: 3.400 (measured 3.43)      RMS error 0.03
    Ea range: 42.6-54.6 kcal/mol = 178-228 kJ/mol

This fits **both** measured points, whereas the single-Ea version matched
350 C and overshot 500 C by 22%.

### THE HEADLINE RESULT: this choice dominates the Enceladus prediction

Extrapolating each model to Enceladus-relevant conditions (1 Myr exposure):

    T        distribution     single-Ea      ratio
    25 C      ~0               ~0            2.5e5 x
    50 C       9e-4 mol/kg     1.2e-8        7.8e4 x
    100 C      0.768 mol/kg    5.8e-4        1.3e3 x

**Single-Ea versus distributed-Ea changes low-temperature organic release by
three to five orders of magnitude.** This is not a refinement — it is the
dominant uncertainty in the entire organic contribution to ocean chemistry,
and it sits precisely in the temperature range the Enceladus cases occupy.

### And the epistemic problem that follows

Per-bin conversion in the 48 h experiments, at the fitted shift:

    Ea kcal    kJ/mol    wt%    f@350C    f@500C
    42.6       178.2      5     1.0000    1.0000
    44.6       186.6     15     1.0000    1.0000
    46.6       195.0     25     1.0000    1.0000
    48.6       203.3     25     1.0000    1.0000
    50.6       211.7     15     0.9547    1.0000
    52.6       220.1     10     0.4596    1.0000
    54.6       228.4      5     0.1152    1.0000

The **low-Ea bins — 70% of the mass — are fully converted at both 350 C and
500 C.** Our calibration data therefore contain essentially **no information
about the shape of the low-Ea end of the distribution**. Only the top two bins
are partially converted at 350 C, so the data constrain the high-Ea end a
little and the low-Ea end not at all.

Yet the low-Ea tail is exactly what determines Enceladus-temperature release.
**We have been fitting the part of the distribution that does not matter for
the application, and the part that does matter is unconstrained by our data.**

This also means the good 350/500 C agreement is **not** evidence that the
model will behave correctly at Enceladus temperatures. It cannot be.

### Which retroactively justifies going after the 250 C data

This inverts the earlier recommendation. Adopting Miller's distribution is out
(wrong material, wrong time regime). Fitting our own distribution shape
requires data at a temperature where the low-Ea bins are only *partially*
converted — and 250 C is the only such data in the paper. So the **syn-IOM
250 C results become the single most valuable remaining calibration target**,
not a secondary check.

Cautions when using them:
1. **Use the R-ratios (RH/C, RO/C, RN/C)**, which are dimensionless and
   comparable across samples with different starting compositions.
2. **Syn-IOM is not Murchison** — Miller needed +8 kcal/mol to fit HC113, so
   the synthetic material is demonstrably more refractory than the meteoritic
   material. A 250 C syn-IOM constraint cannot be transferred to Murchison
   without accounting for that offset, and the offset is itself known only at
   geologic heating rates. State this confound explicitly.

### Residue data availability (unchanged)

Murchison 500 C H/C and O/C are **blank in Miller's Table 5**; that row
reports only N/C. A gap in the published data, not a transcription miss.


---

# PART 5 — STATUS AND NEXT ACTIONS

## Done
- IOM module split from ENTICES; first real execution completed.
- Test-harness regex bug found and fixed (module itself was correct).
- IOM_CO2/IOM_CH4 Ea corrected from a ratio fit to an absolute fit (an 18%
  error at 350C, found only by building the standalone test).
- IOM_N recalibrated the same way (215 -> 217.3 kJ/mol).
- Part B converted from diagnostic printout to 16 asserted checks; passing.
- Murchison 350C residue comparison established as first independent
  validation point.

## Next, in recommended order
1. **Check the redox budget in existing output** (Part 2b). Cheapest test,
   and it may invalidate existing organic-enabled H2 and mineral-assemblage
   predictions. Look at H2, pe, *and* Fe phases, not H2 alone.
2. **Move IOM_CO2 to a distributed Ea** (Part 4). The single-Ea version is
   wrong by 3-5 orders of magnitude at Enceladus temperatures, which is a
   larger error than anything else currently on this list. Implement as
   sub-pools using the BS89 shape with a fitted shift (+0.6 kcal/mol,
   m0 3.40 Murchison basis / ~4.13 ENTICES basis). Note this is a
   *provisional* fit: the shape is a prior, not constrained by our data.
3. **Extract the syn-IOM 250 C data** (Table 2 yields and Table 5 residue
   R-ratios). This is now the highest-value calibration target, since it is
   the only data that can constrain the low-Ea tail. Carry the
   syn-IOM-vs-Murchison confound explicitly.
4. **Build Option 3**, in the four-step order in Part 3 — database in
   isolation first. Still required for low-T speciation realism, but note
   that item 2 changes *how much* organic carbon there is to speciate, so
   doing 2 first avoids re-running everything.
5. Apply the same distributed treatment to IOM_CH4 once CO2 is working
   (BS89 gives CH4 precursors higher Ea, up to ~74 kcal/mol).
6. Wire the module into ENTICES's template generator (`INTEGRATION POINTS`
   block in `iom_module.R`).
7. IOM_S remains uncalibrated and will stay that way unless a source for H2S
   yield appears.

## Open questions for the group
- Should decoupling be a global mode, or temperature-dependent (decouple below
  some T)? Global is cleaner and gives an honest bracket; a threshold is more
  physical but introduces an arbitrary discontinuity. Current recommendation:
  global mode switch, per-element configurable, and let the bracket carry the
  uncertainty.
- Unresolved since early on: the nitrogen composition discrepancy (30 g N/kg
  -> N_2.142 vs an older comment citing N_3.284). Still not confirmed against
  the source composition.
