# data/

## `mineral_molar_volumes.csv`

A molar volume for every non-gas PHASES entry of the three coreclath databases, plus Schreibersite
(which the generators add as a KINETICS `-formula`, not a PHASE). `bulk_density()`
(`phreeqc_tools.R`) uses it. Regenerate it with

```
python scripts/export_mineral_volumes.py --ds636 <tc-ds636.txt> --hp636 <hp636ver.dat>
```

Columns: `phase`, `formula` (as written in the database), `molar_mass_g_mol` (from the Core11
element weights), `Vm_cm3_mol`, `source`, `priority_level`, `databases` (which coreclath files define
the phase), `production` (`yes` if a generator lists it in EQUILIBRIUM_PHASES or KINETICS).

All volumes are 25 °C, 1 bar values. No thermal or compressibility correction is applied.

### Source priority (PI decision, 2026-10-28)

The first source that has a phase is used.

1. `legacy/Core11_idealgas_mod_v4.dat`, `-Vm`: 381 of the 388 rows, and every production phase but
   Schreibersite and the hydrates. The four hydrates (CH4/Mtg/CO2/H2S), which only the coreclath databases have,
   take their coreclath value at this level: CH4 and Mtg `Vm 136.112` (written without the hyphen),
   CO2 `-vm 137` (Aya et al. 1997, in the comment), H2S `-vm 133.99` (Carroll 2009, in the comment).
2. Green et al. (2025) `tc-ds636.txt`, from the bundle `tc-thermoinput-igG25-H18corr-2024-12-21`
   (PI's copy). THERMOCALC units: V in kJ/kbar, which is J/bar, so ×10 gives cm³/mol.
3. Perple_X `hp636ver.dat` (TC-DS636, 3 May 2023), the raw file from
   https://github.com/jadconnolly/Perple_X/blob/main/datafiles/hp636ver.dat. V0 in J/bar, ×10.
4. Values the PI supplied: **Schreibersite**, 0.60221 × 366.10 / 8 = 27.56 cm³/mol (ρ = 7.20 g/cm³),
   from mindat.org (min-3582) unit-cell data, tetragonal I-4, a = 9.05 Å, c = 4.47 Å, V = 366.10 Å³,
   Z = 8. Fe₃P end-member assumed: mindat gives (Fe,Ni)₃P, ENTICES uses Fe₃P.
   [citation needed: the primary crystal-structure reference mindat lists; mindat refused automated
   access, so it was not checked]

Levels 2 and 3 are only used through an explicit name map, and only where the stoichiometry is
proportional to the PHREEQC formula (the ratio scales the volume). In the current table no phase
needed levels 2 or 3; they serve as a cross-check.

**Unfilled:** Fluorite and Fluorapatite. Core11 defines both without a `-Vm`, and neither ds636 nor
hp636 has them. Neither is in a production list. For the PI to source.

### Cross-checks

Level 2/3 values that differ from Core11 by more than 2% (Core11 kept):

| Phase | Core11 | ds636 / hp636 | Difference |
|---|---|---|---|
| Celadonite | 157.10 | 139.57 | −11.2% |
| Tridymite | 26.12 | 28.00 | +7.2% |
| Merwinite | 104.40 | 98.47 | −5.7% |
| Troilite | 19.07 | 18.19 | −4.6% |
| Greenalite | 115.00 | 119.80 | +4.2% |
| Tephroite | 47.52 | 48.99 | +3.1% |
| Rhodonite | 35.87 | 34.94 | −2.6% |
| Clinochlore-14A | 207.11 | 211.40 | +2.1% |

Troilite, Greenalite, Tephroite and Clinochlore-14A are production phases.

Clathrates against Sloan & Koh (2007), *Clathrate Hydrates of Natural Gases*, 3rd ed., Table 2.5a
(p. 78), structure I lattice parameter a (the 0 °C value unless marked), 46 H₂O per cell, so for the
formula X·6H₂O, Vm = 0.60221 a³ / (46/6):

| Hydrate | a (Å) | Sloan & Koh Vm | Database Vm | Difference |
|---|---|---|---|---|
| CH4 (and Mtg) | 11.981 (Ogienko et al. 2006, cited by S&K) | 135.09 | 136.112 | +0.8% |
| CO2 | 12.07 | 138.12 | 137 | −0.8% |
| H2S | 12.02 | 136.41 | 133.99 | −1.8% |

All within 2%, so the database values are kept. The database formulas fix the hydration number at
6; real sI hydrates run about 5.75–6.4 (the CO2_hydrate comment says so).

### Provenance of Core11's volumes

Core11 (Neveu 2017, extended by Diab 2021 and Melwani Daswani 2025) names molar-volume sources in its
header for some entries: WebMineral (38 solids), Marion et al. 2005, 2008, 2009 (Arcanite, Gypsum,
Niter, Thenardite, Melanterite, Alum-K), McColm 2013 (CaUO4), Hovis et al. 2004 (NH4-muscovite) and
WebElements (K, U). Many `-Vm` lines carry a comment naming `thermo.com.V8.R6+.tdat` or
`Webmineral.com`, copied into the `source` column. For the rest: [provenance needed].
