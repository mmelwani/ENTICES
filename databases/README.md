# PHREEQC databases

| File | Used by | Notes |
|---|---|---|
| `coreclath_CH4coupled.dat` | Every generator when CH₄ is redox-coupled: porewater T ≥ 150 °C, or `CH4_redox_override = "coupled"` | Lucas's production database: Core11 plus clathrate phases. |
| `coreclath_CH4uncoupled.dat` | Every generator when CH₄ is decoupled (porewater T < 150 °C by default, or `CH4_redox_override = "decoupled"`); species are `Mtg`/`Mtg(g)` | Decoupled-methane variant. |
| `Core11_idealgas_mod_v4.dat` | The PHREEQC tests in `tests/` (except B3) | Mohit's base database; has no clathrates. |

The generator writes the database it expects into each `.pqi` header
(`# Recommended database: ...`), chosen by `resolve_ch4_mode()`.

Consolidating these into one canonical source plus generated variants is
pending, and is Lucas's call.
