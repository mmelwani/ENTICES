# PHREEQC databases

| File | Used by | Notes |
|---|---|---|
| `coreclath_CH4coupled.dat` | Every generator when CH₄ is redox-coupled: porewater T ≥ 150 °C, or `CH4_redox_override = "coupled"` | Lucas's production database: Core11 plus clathrate phases. |
| `coreclath_CH4uncoupled.dat` | Every generator when CH₄ is decoupled (porewater T < 150 °C by default, or `CH4_redox_override = "decoupled"`); species are `Mtg`/`Mtg(g)` | Decoupled-methane variant. |

The generator writes the database it expects into each `.pqi` header
(`# Recommended database: ...`), chosen by `resolve_ch4_mode()`.

`Core11_idealgas_mod_v4.dat` (Mohit's base database) was retired to `legacy/` on
2026-10-23: `coreclath_CH4coupled.dat` contains every Core11 entry plus the
clathrate phases, and the tests now use it.

Consolidating the two coreclath files into one canonical source plus generated
variants is pending, and is Lucas's call.
