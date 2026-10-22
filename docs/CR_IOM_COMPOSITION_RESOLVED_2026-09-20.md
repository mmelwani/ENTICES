# The IOM composition is correct: it is Miller et al. (2026) SI Table 9, exactly

Date: 2026-09-20
Prompted by: PI — "our paper (Miller et al. 2026) used organics from CR
chondrites instead, for Titan's building blocks, since they appear to be less
altered", then "Miller+ (2026)'s IOM composition (SI Table 9) is H = 46 mg/g,
O = 159 mg/g, C = 728 mg/g, N = 30 mg/g and S = 37 mg/g."

**Supersedes `claude/IOM_COMPOSITION_PROBLEM_2026-09-20.md` §1, §2, §3, §6, §7,
§8** and the "switch to Murchison" recommendation in `IOM_SULFUR_SOURCED` §3.
Also supersedes an earlier version of *this* file, which recommended changing
sulfur — that recommendation is withdrawn; see §3.
Method: arithmetic. Reproducible from the numbers below.

---

## Short answer

`iom_bulk_formula()` is **Miller et al. (2026) SI Table 9**, converted to mol/kg.
Not approximately — to 1 part in 10⁴ on all five elements:

    el    SI Table 9    /atomic wt    iom_bulk_formula()   ratio
    C      728 mg/g     12.011          60.611  mol/kg     1.0000
    H       46 mg/g      1.008          45.635  (45.639)    1.0001
    O      159 mg/g     15.999           9.938  ( 9.939)    1.0001
    N       30 mg/g     14.007           2.142               1.0001
    S       37 mg/g     32.060           1.154               0.9999
           ----------
           1000 mg/g

**Nothing in the composition needs changing.** There was no error in the code.
The provenance is now traced rather than inferred, and it is the PI's own
published composition for CR-chondrite IOM.

## 1. What this retracts

Two of my earlier conclusions were wrong, and they were wrong in opposite
directions. Recording both, because the failure mode matters more than the
numbers.

**(a) `IOM_COMPOSITION_PROBLEM` §1: "this is a terrestrial kerogen composition."**
False. The reasoning was: H/C ≈ 0.75 with O/C ≈ 0.16 sits on the Type III
vitrinite path, therefore it is vitrinite. But CR-chondrite IOM sits in nearly
the same place on a Van Krevelen plot, so the test had no discriminating power.
I compared against a *remembered* CI (0.70/0.22) and a half-remembered CR
(0.78/0.16); when the formula landed close to the latter I reached for the
kerogen explanation rather than the obvious one. The real CI rows (Alexander
2007: Orgueil 0.673/0.181, Ivuna 0.723/0.154) also contradict the remembered
values, so the "CI is more oxidised, therefore the offset has the wrong sign"
argument fails at its premise as well.

**(b) The first version of this file: "sulfur is low by 1.8×, change it to
S/C = 0.036."** Also false — S = 37 mg/g is the published value. That
recommendation came from comparing against my own unweighted mean of Alexander
2007's CR rows, and treating my average as the authority over the paper. It was
not. §3 below keeps the one genuine question that comparison surfaced.

**The common error in both:** using recalled literature values, and then my own
recomputation from a source, as the reference against which to judge a number
whose provenance I had not actually traced. The provenance was one question away
throughout.

## 2. The closure is the source's, not ours

`IOM_COMPOSITION_PROBLEM` §3 treated "C+H+O+N+S = exactly 1000 g" as a modelling
choice ENTICES made, and identified it as the origin of the error class that
produced the 17% `Ea_N` bug. Half of that is wrong: the closure is how SI Table 9
reports the composition. It arrives closed.

What remains true is the narrower point: converting *from* this closed basis *to*
Murchison's measured 82.1 wt% basis (the ×1.215 factor) is ours, and that
conversion is what caused the `Ea_N` error. So the hazard is real but it lives in
the conversion step, not in the composition. Worth stating that way in PROGRESS,
because "our normalisation is a source of error" invites someone to renormalise
the formula, which would now break agreement with the paper.

## 3. The one open question: where SI Table 9's sulfur comes from

C, H, O and N reproduce from Alexander et al. (2007) Table 2 straightforwardly.
Unweighted mean of the seven CR rows:

    meteorite    C wt%   H/C    N/C    O/C    S/C     (ratios x100 atomic)
    MET00426      69.0   78.3   3.39   14.4   n.m.
    QUE99177      67.8   80.3   3.65   14.6   n.m.
    EET92042      66.5   75.7   3.83   16.6   3.2
    GRA95229      63.0   68.9   4.39   22.0   3.1
    Al Rais       61.4   80.6   3.21   16.2   7.4
    GRO95577      75.5   78.6   3.16   11.4   0.7
    LEW85332      56.7   71.9   3.19   20.4   n.m.
    -------------------------------------------------
    mean          65.70  0.7633 0.0355 0.1651 0.0360   (S: n = 4 of 7)

    Table 9 implies   H/C 0.7529   N/C 0.0353   O/C 0.1640   S/C 0.0190

H/C, N/C and O/C agree to within 1.3% — well inside the group's own scatter
(H/C sd ±0.044, O/C ±0.036, N/C ±0.0045). Good.

**Sulfur does not reproduce.** Table 9's S/C = 0.0190 sits below both the
four-CR mean (0.0360) and the mean excluding Al Rais (0.0233, the most
aqueously altered CR and a 0.074 outlier that dominates the average). Three
possible explanations, in rough order of likelihood:

1. A different source for S. Alexander 2021 measures IOM S directly and finds no
   systematic CM/CR/C2 difference, group-wide S/C = 0.021–0.090, mean ≈ 0.041 —
   higher still, so this does not explain 0.019 either.
2. A subset, weighting, or exclusion I have not reconstructed.
3. A different normalisation of S specifically — e.g. correcting for the sulfide
   contamination both papers warn about.

**Why it is worth resolving.** Alexander 2007 §2.3.3 states its S values "*were
not determined directly, but were estimated by difference after summation of the
H, C, N, and ash weight fractions*" — a method requiring the absence of
oxidisable phases, which the same paragraph says fails often enough to produce
>100 wt% totals in many residues. So the CR S/C numbers I compared against are
the weakest column in that table, and a lower published value may well be better
grounded than my average of them. But `IOM_S` is also the channel with **no
kinetic calibration at all** (Miller measures no H2S), so its m0 is now the one
number in the formula I cannot trace to a source. Composition and rate
uncertainty stack on the same channel.

**This is a question for the PI, not a defect.** No change recommended.

## 4. The CR-composition / CM-kinetics mismatch

Still real: composition is CR, kinetics are fitted to Miller's Murchison (CM2)
hydrous-pyrolysis experiments. But it is smaller than I implied earlier, and
smaller than the analytical spread on Murchison itself:

    ratio    CR mean    Murchison (Miller T1)    Murchison (Alexander 2007)
    H/C      0.763      0.699                    0.588
    O/C      0.165      0.203                    0.183
    N/C      0.0355     0.0355                   0.0327

CR vs Miller's Murchison: +9% H/C, −19% O/C. But the *two* Murchison analyses
differ from each other by 19% on H/C — different isolation procedures and labs.
The inter-lab spread on one meteorite exceeds the inter-group difference, which
makes it hard to argue the CR/CM mismatch is the dominant error term.

It is a caveat to state, not a blocker. And CR is the right choice on the
science: Alexander 2007 finds CR IOM the least thermally processed of any group
(highest H/C, most anomalous D and ¹⁵N, "*shows little or no evidence of having
been heated*").

## 5. What is genuinely unresolved: the kinetic distribution shape

Alexander 2007 §3.3:

> "Terrestrial kerogens evolve along a path that is quite distinct from the
> distribution of chondritic IOM compositions (Fig. 1). Unlike terrestrial
> kerogens, the IOM H/C and N/C ratios are correlated and the H/C ratio
> decreases before the O/C ratio."

and attributes the N behaviour to functionality: IOM N in **amines** (labile),
coal N in **pyridines**. Our whole kinetic scheme is Vitrimat 2018, fitted to
vitrinite. Two consequences, which must be kept separate:

**(a) The Ea *distribution shape* is untested against chondritic data.** This
stands. `IOM_N`'s Ea = 242.8 kJ/mol was refitted to Miller's *Murchison*
measurement rather than inherited from Burnham, so the *centre* is on chondritic
data — but it is a single Ea riding a vitrinite-shaped `logA = 15.301`, and the
*width* of the CO2/CH4/CHn distributions is Burnham's vitrinite width throughout.
Nothing has tested that width.

**(b) Fig. 1's ordering argument does not apply to us.** Alexander's "H/C falls
before O/C" describes the *dry metamorphic* sequence across meteorite parent
bodies. Miller's experiments are 48 h *hydrous* pyrolysis, where water is a
reagent — this project already found ~a quarter of product-CO2 oxygen comes from
water, not from IOM. And Miller's 350 °C Murchison residue gives RH/C 0.90
against RO/C 0.53: O/C falling *faster*, the coal-like ordering, opposite to
Alexander's dry trend. Since Enceladus pore water is liquid water at modest T for
Gyr, Miller's hydrous experiments are the right analogue and Fig. 1 is not the
correct test. I nearly mis-transferred this; recording it so the next reader does
not.

So (b) defuses the framing, (a) survives as a real open item. It should be
flagged as a question about distribution *width*, distinct from the Ea *values*.

## 6. Actions

1. **Change nothing in `iom_bulk_formula()`.** It is Miller 2026 SI Table 9.
2. **Correct the provenance statements** in `IOM_COMPOSITION_PROBLEM` (§1 "looks
   like terrestrial kerogen", §3 "the normalisation is a modelling choice", §6–§7
   "switch to Murchison") and in `RESIDUE_OVERSHOOT_RESOLVED` (§"Where the
   starting composition comes from", which describes it as CI-like). Both should
   read: CR-chondrite IOM, Miller et al. 2026 SI Table 9, from Alexander et al.
   2007's seven CRs.
3. **Ask the PI where Table 9's S = 37 mg/g comes from** (§3), and record the
   answer next to `IOM_S`'s uncalibrated-rate flag.
4. **State the CR/CM mismatch** with §4's numbers showing it is smaller than the
   inter-lab spread on Murchison.
5. **Open item: Vitrimat distribution width** (§5a), explicitly separated from
   the Ea values, which are on chondritic data.
6. **Leave `Ea_N` = 242.8 alone.** Fitted to Murchison, not inherited from coal,
   so the amine-vs-pyridine argument does not touch it.

## 7. What the residue offsets now mean

`RESIDUE_OVERSHOOT_RESOLVED` found starting H/C +7.6% and O/C −8.9% against
Murchison's measured values, and treated these as possible composition errors.
They are not errors: they are the real CR-vs-CM difference (§4), and they are
expected. The document's numerical result stands; its interpretation should
change from "the starting composition may be wrong" to "the starting material is
CR and the residue data are CM, so an offset of this size and sign is the
correct prediction." That makes the R-ratio workaround the right approach rather
than a patch, since ratios cancel the material difference.
