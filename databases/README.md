# PHREEQC databases

| File | Used by | Notes |
|---|---|---|
| `coreclath_CH4coupled.dat` | `carbon_mode = "coupled"` (older spelling `CH4_redox_override = "coupled"`), or `carbon_mode = "temperature"` with porewater T ≥ `ch4_cutoff_C` (150 °C) | Lucas's production database: Core11 plus clathrate phases. All carbon redox-coupled: an upper bound on abiotic methanation. |
| `coreclath_CH4uncoupled.dat` | `carbon_mode = "ch4_uncoupled"` (older spelling `CH4_redox_override = "decoupled"`), or `carbon_mode = "temperature"` with porewater T < `ch4_cutoff_C`; species are `Mtg`/`Mtg(g)` | Decoupled-methane variant. |
| `coreclath_C4only.dat` | **Default** (`carbon_mode = NULL` or `"c4_only"`), at every temperature | `coreclath_CH4uncoupled.dat` with every coupled carbon species, phase and master species that is not C(+4) commented out (`#c4only#`); the header lists each removal. Methane stays decoupled as `Mtg`. Generators drop the removed names from their inputs (`restrict_to_database()`). |

The generator writes the database it expects into each `.pqi` header
(`# Recommended database: ...`), chosen by `resolve_ch4_mode()`.

`Core11_idealgas_mod_v4.dat` (Mohit's base database) was retired to `legacy/` on
2026-10-23: `coreclath_CH4coupled.dat` contains every Core11 entry plus the
clathrate phases, and the tests now use it.

Consolidating the two coreclath files into one canonical source plus generated
variants is pending, and is Lucas's call.
