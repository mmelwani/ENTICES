# reconstruction/

Reconstruction of the primordial (pre-alteration) CI chondrite mineralogy that
supplies ENTICES's rock composition. It is kept as a separate module: ENTICES
depends only on its output, `rock_mol_per_kg.csv`.

| File | What it is |
|---|---|
| `Primordial_CI_chondrite_mineralogy.Rmd` | The reconstruction (2026-09-16 version). A non-negative least-squares fit over a fixed phase list. The headline result is Scenario R against Lodders et al. (2025) Table 4. Its last chunk writes `rock_mol_per_kg.csv`. |
| `Primordial_CI_chondrite_mineralogy_montecarlo.Rmd` | Monte Carlo uncertainty propagation. Requires objects from the main Rmd. |
| `RECONSTRUCTION_STATUS.md` | Status, verified results, errors found, and the open items (§8). |
| `rock_mol_per_kg.csv` | Machine-readable output, one row per phase, with columns `phase, formula, mw_g_mol, wt_pct, mol_per_kg, role, scenario`. |
| `notes/` | Tephroite kinetics notes and the hand-off to Lucas. |
| `archive/v0_uncommitted_edit_2026-09.patch` | An uncommitted edit (density citations) to the discredited v0 Rmd, saved before it was removed. **Don't use v0 for numbers**; the v0 file itself is in git history. |

## Regenerating the CSV

Knit `Primordial_CI_chondrite_mineralogy.Rmd` from this folder. Then run
`Rscript tests/rock_table_check.R` from the repo root to check ENTICES's
`ROCK_MOL_PER_KG` against it.

The check allows 1e-4 relative or 1e-6 absolute. It applies two
ENTICES-side steps before comparing:
- ENTICES adds lawrencite (FeCl₂) and takes its Fe out of the metal.
- Magnesiochromite is inert: it stays in the mass basis but isn't reacted.

## Still open

These are tracked in `RECONSTRUCTION_STATUS.md` §8:
- the King et al. forward check (§8 item 8);
- surface areas for every kinetic phase (§8 item 5);
- the schreibersite rate (§8 item 1b).

The tephroite activation energy is now fitted by Lucas, as
`Ea = −60.11·ln(pH) + 106.17` kJ/mol floored at 0; see `CHANGELOG.md`. §4
still says "UNRESOLVED" and needs updating.
