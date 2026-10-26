# PHREEQC databases

| File | Used by | Notes |
|---|---|---|
| `coreclath_CH4coupled.dat` | Every generator when CH₄ is redox-coupled: porewater T ≥ `ch4_cutoff_C` (150 °C), or `carbon_mode = "coupled"` (older spelling `CH4_redox_override = "coupled"`) | Lucas's production database: Core11 plus clathrate phases. |
| `coreclath_CH4uncoupled.dat` | Every generator when CH₄ is decoupled (porewater T < `ch4_cutoff_C`, 150 °C by default, or `carbon_mode = "ch4_uncoupled"`, older spelling `CH4_redox_override = "decoupled"`); species are `Mtg`/`Mtg(g)` | Decoupled-methane variant. |
| `coreclath_C4only.dat` | Every generator with `carbon_mode = "c4_only"` | `coreclath_CH4uncoupled.dat` with every coupled carbon species, phase and master species that is not C(+4) commented out (`#c4only#`); the header lists each removal. Methane stays decoupled as `Mtg`. Generators drop the removed names from their inputs (`restrict_to_database()`). |

The generator writes the database it expects into each `.pqi` header
(`# Recommended database: ...`), chosen by `resolve_ch4_mode()`.

`Core11_idealgas_mod_v4.dat` (Mohit's base database) was retired to `legacy/` on
2026-10-23: `coreclath_CH4coupled.dat` contains every Core11 entry plus the
clathrate phases, and the tests now use it.

Consolidating the two coreclath files into one canonical source plus generated
variants is pending, and is Lucas's call.
