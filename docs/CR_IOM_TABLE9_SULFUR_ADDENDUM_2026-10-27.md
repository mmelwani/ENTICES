# Where SI Table 9's IOM sulfur comes from: addendum (2026-10-27)

Companion to `docs/CR_IOM_COMPOSITION_RESOLVED_2026-09-20.md` (§3), and copied
from §8 of the claude.ai Project copy of that file.

Prompted by: PI, "I don't even remember what Table 9 is, let alone where its
sulfur comes from."

**What SI Table 9 is.** Miller et al. (2026) Table 1's caption says the
elemental abundances of the accretionary components are given in
Supplementary Table 9. Those components are IOM, anhydrous cometary dust,
hydrated asteroidal dust and ices. The IOM row is C 728, H 46, O 159, N 30,
S 37 mg/g. Main text §2.1.2 says it is "the average CHONS elemental
abundances for IOM from the seven CR chondrites reported by Alexander et al.
(2007)".

**The test.** Average the seven CR rows of Alexander 2007 Table 2 as wt%. Each
meteorite's H, N, O and S wt% comes from its C wt% and its atomic ratio. Then
renormalise CHONS to 1000 mg/g. Only S is missing for three meteorites
(MET00426, QUE99177 and LEW85332 are n.m.), so there are two ways to handle
those three:

    case                              C      H     N     O      S    (mg/g)
    SI Table 9 (published)           728    46    30    159    37
    missing S counted as 0           727    46.7  30.0  157.6  38.6
    missing S skipped (mean of 4)    707    45.4  29.2  153.2  65.7

Counting the three missing values as **zero** reproduces Table 9:
- C, H, N and O agree within 1%;
- S agrees within 4%;
- skipping the missing values overshoots S by 1.8×.

The raw CHONS sum is 90.4 wt%. The renormalisation removes the remaining
~9.6%, which is ash and unmeasured matter. This is arithmetic on the published
table, not confirmation; only the spreadsheet behind Table 9 can confirm it.

**What it implies if confirmed.**
- Table 9's S is about **1.7× low** for the CR group. Averaging only measured
  values gives S ≈ 66 mg/g.
- Alexander et al. (2021/22) measured IOM S directly. Their four CRs give
  S/C 0.021–0.043, mean 0.035, which works out to about 65 mg/g on the same
  closure.
- So two independent routes give about 65 mg/g, not 37.

**Effect on ENTICES.**
- **Small for the S budget.** IOM S is 1.154 mol/kg IOM. At 15 wt% organics
  that is about 0.17 mol S per kg rock, against troilite's 2.38. Doubling it
  adds about 7% to total S.
- **Larger for the IOM channel itself.** `IOM_S`'s m0 doubles, and that channel
  still has no kinetic calibration.
- **A re-closure moves C, H, O and N by about 3%**, which shifts `IOM_CO2`,
  `IOM_CH4` and the rest.

**Decision for the PI:**
- keep Table 9 as published, for consistency with Miller 2026; or
- correct S, as a footnoted ENTICES deviation from the paper.

Either way, record the answer next to `IOM_S`.
