# ENTICES_LF_edited.R — Changelog

## stitch_plots.R

- Added `run_mode` toggle ("bulk" / "tidal_kinetic") to control phase splitting behaviour
- Fixed k_* patching bug: save `.is_react_pre` before patching so MIX rows are not misclassified as REACT
- Added `plot_pe` block: pe vs time for porewater/ocean/bulk phases
- Added `plot_Eh` block: Eh (V) computed from pe using phase-appropriate temperature (T_pw / T_o)
- Added `chem_filter_min` boolean list: drop species from chemistry plots whose max value never exceeds 1e-9
- Added react phase as third source in tidal chemistry plots (solid/dashed/dotted) so H2(g) appears
- Added `plot_mix4_secondary_min` block: equilibrium minerals from MIX4 (soln=4) rows
- Added `plot_H2_rate` block: rate of change of H2(aq) and H2(g) in mol/year (linear scale)
- Added IOM kinetic phases (IOM_labile, IOM_mid, IOM_refract) to primary moles and primary rates plots
- Converted zero mineral moles to NA in secondary minerals and primary moles plots (breaks lines at depletion rather than connecting through zero)
- Temperature reading now conditional on `run_mode`: tidal reads T_pw and T_o; bulk reads only T_o

## ENTICES_LF_edited.R

- Added IOM kinetic phases to `generate_phreeqc_input_tidal_kinetic_restart`: pools carry over between restarts, reading remaining moles from k_IOM_* columns with fallback to initial fractions
- Fixed `generate_phreeqc_input` SELECTED_OUTPUT `-kinetics` line to include IOM pool names
- Added porewater temperature tag (e.g. T100) to output filenames in both main generator and restart generator
- Added IOM_labile, IOM_mid, IOM_refract as valid `depletion_mineral` options in restart generator
- Fixed Enstatite neutral-weathering activation energy in porewater restart: corrected from -79.0e3 to -80.0e3 J/mol to match all other generators
- Merged `generate_phreeqc_input_porewater_kinetic_restart` and `generate_phreeqc_input_hydrosphere_kinetic_restart` into single `generate_phreeqc_input_bulk_kinetic_restart(mode = c("porewater","hydrosphere"))`; originals saved to restart_functions_backup.R
- Merged `generate_phreeqc_input_hydrosphere_kinetic` and `generate_phreeqc_input_porewater_kinetic` into single `generate_phreeqc_input_bulk_kinetic(mode = c("porewater","hydrosphere"))`; mode controls water mass source (water_rock_ratio vs combined_mass_norm) and mineral moles scaling (1 kg vs total_rock kg); originals saved to kinetic_functions_backup.R
- Added `dissolve_only` keyword to all primary minerals in EQUILIBRIUM_PHASES blocks except Magnetite, across all four generator functions that use primary EP entries (`equil`, `hydrosphere_equil`, `bulk_kinetic_restart`, `tidal_kinetic_restart`)
- Added 28 new equilibrium minerals to all 4 primary generator functions (`generate_phreeqc_input`, `bulk_kinetic`, `bulk_kinetic_restart`, `tidal_kinetic_restart`); total secondary mineral count raised from 49 to 77; minerals added: Aragonite, C2H4(g), C2H6(g), C3H8(g), Calcite, CH4(g), Chalcedony, Chamosite, Chrysotile, Citric_Acid, Clinochlore-14A, Clinochlore-7A, Clinozoisite, CO(g), CO2(g), Daphnite-14A, Daphnite-7A, Dawsonite, Diaspore, Dolomite, Fe2(SO4)3, Glycine, Gyrolite, Huntite, Hydromagnesite, Lizardite, N2(g), NH3(g)
- Added Magnetite to EQUILIBRIUM_PHASES 2 (ocean/seafloor phase) in `generate_phreeqc_input`; SELECTED_OUTPUT `-equilibrium_phases` updated; EP1 (porewater) unchanged
- Added Magnetite to EQUILIBRIUM_PHASES 2 (ocean) in `generate_phreeqc_input_tidal_kinetic_restart`; EP2 now also includes any `equil_minerals` primary minerals (with `dissolve_only`, except Magnetite) sourced from the ocean CSV state; Magnetite inserted alphabetically between Lizardite and Melanterite; SELECTED_OUTPUT updated; EP1 and `bulk_kinetic_restart` unchanged

### Patch C — Miller-calibrated IOM channels (2026-08-20)

Replaces the three composition pools (`IOM_labile`, `IOM_mid`, `IOM_refract`) in
`generate_phreeqc_input` with four channels calibrated against Miller et al. (2025)
GCA 390:38-56 (Murchison IOM, Tables 2, 5, 9).

**Affected function:** `generate_phreeqc_input` (~lines 942–976)

| Channel  | m0 / kg IOM | -formula (what the solid loses) | Ea (kJ/mol) |
|----------|-------------|----------------------------------|-------------|
| IOM_CO2  | 4.169       | C 1 H 1 O 1.3115                | 214         |
| IOM_CH4  | 1.908       | C 1 H 4.435                     | 228         |
| IOM_N    | 0.875       | N 1 H 2                         | 215         |
| IOM_S    | 0.462       | S 1 H 1                         | 215 (uncal) |

Key design decisions:
- `-formula` is **what the solid loses**, not the product molecule. PHREEQC draws
  remaining oxygen from water and speciates at the prevailing pe.
- Only ~10% of IOM carbon and ~40% of IOM N is ever released (inert char remainder
  is deliberately absent). The old 3-pool model would dissolve 100% over hot/long runs.
- Channels stored in `self$iom_pools`; set it before calling to override defaults.
- Arrhenius parameters: A = 10^13 s^-1 per channel; Ea fitted two-point (350 vs
  500 °C, 48 h). S channel is uncalibrated (Miller does not report H2S yields).
- Single-Ea caveat: cannot reproduce the low-Ea CO2 tail Miller sees at 250 °C
  (sub-pools to follow).

**Supersedes** the entries in the ENTICES_LF_edited.R section above that reference
`IOM_labile`, `IOM_mid`, `IOM_refract`, and `k_IOM_*` columns.

---

### Patch D1 — USER_PUNCH IOM column emission (2026-08-20)

**Affected function:** `generate_phreeqc_input` cycle-block builder (~lines 1493–1508,
1531–1544)

Each per-cycle `USER_PUNCH` block now emits four additional columns:

    IOM_CO2_mol  IOM_CH4_mol  IOM_N_mol  IOM_S_mol

written via `PUNCH KIN("IOM_CO2")` etc. at BASIC line numbers 210, 220, 230, 240.
Column set is generated from `self$iom_pools`, so a non-default channel configuration
stays consistent automatically. Columns are present only when `organic_wt_percent > 0`.

These columns are the input that Patch D2 (below) reads on restart.

---

### Patch D2 — IOM channel recovery on restart (2026-08-20)

**Affected functions:**
- `generate_phreeqc_input_bulk_kinetic_restart` (~lines 2748–2813, 2932–2985, 3091–3096)
- `generate_phreeqc_input_tidal_kinetic_restart` (~lines 3409–3476, 3484–3503, 3569–3571, 3696–3704)

Both restart functions now:

1. **Recover IOM channel moles** from the `IOM_*_mol` columns written by Patch D1.
   - Bulk restart reads from `last_row` (the final CSV row).
   - Tidal restart reads from `row_react` (the porewater row of the final cycle).
2. **Rebuild RATES** — same Arrhenius rate laws as the primary generator; regenerated
   from `self$iom_pools`, not read from the CSV.
3. **Rebuild KINETICS** — each channel's `-m0` is set to the recovered remaining moles
   (not the original initial moles).
4. **Re-emit USER_PUNCH columns** — IOM headings and `PUNCH KIN(...)` lines are
   appended to every restart cycle block, so the next restart continues to have
   what it needs.
5. **Guard** — if `organic_wt_percent > 0` but any `IOM_*_mol` column is absent from
   the CSV, `stop()` with a message naming the missing columns and pointing to
   Patch D1. Silently restarting with zero organics is the failure mode most likely
   to survive undetected into a published run.
6. **`depletion_order`** in tidal restart updated from old pool names to
   `IOM_CO2`, `IOM_CH4`, `IOM_N`, `IOM_S` (~line 3569).

`self$iom_pools` is re-seeded to the default 4-channel table in both restart
functions if it is `NULL` at call time (e.g. if the object was reconstructed
without re-running the primary generator).

---

### Surface-area/molecular-weight rate scaling (2026-09-10)

Ported the updated Palandri & Kharaka rate law from `mineral_dissolution_constants.r`
into the `RATES` blocks of all four primary generator functions
(`generate_phreeqc_input`, `generate_phreeqc_input_bulk_kinetic`,
`generate_phreeqc_input_bulk_kinetic_restart`, `generate_phreeqc_input_tidal_kinetic_restart`),
for all 8 primary minerals (Forsterite, Fayalite, Enstatite, Ferrosilite,
Pyrrhotite, Anorthite, Albite, Magnetite).

- Old: `k = kacid + kneut [+ kbase]` — mol m⁻² s⁻¹, the bare PK04 intrinsic rate constant.
- New: `k = (kacid + kneut [+ kbase]) * SSA * mw * M` — mol s⁻¹, scaled by specific
  surface area, molecular weight, and current remaining moles.

`M` is PHREEQC's built-in "current moles of this kinetic reactant" — the BASIC-code
equivalent of the R script's `mol` parameter, so the port is exact rather than
approximate.

| Mineral | SSA (m²/g) | mw (g/mol) |
|---|---:|---:|
| Forsterite | 0.1 | 140.6715 |
| Fayalite | 0.1 | 203.7555 |
| Enstatite | 0.1 | 100.3725 |
| Ferrosilite | 0.1 | 131.9145 |
| Pyrrhotite | 5 | 87.913 |
| Anorthite | 5 | 278.164 |
| Albite | 5 | 262.1798 |
| Magnetite | 5 | 231.517 |

Only the `k = ...` line changed in each `RATES` block; `rate = k * (1 - SR(...))`,
`moles = rate * TIME`, and `SAVE moles` are untouched, since `M` is folded into `k`
(mirroring how `k_primary_mineral()` in `mineral_dissolution_constants.r` already
folds `mol` into its returned `k`, so `dissolution_per_step()` needs no separate
moles multiplication).

`generate_phreeqc_input_bulk_kinetic_restart` and `generate_phreeqc_input_tidal_kinetic_restart`
had byte-identical `.rate_laws` list blocks for these 8 minerals, so both were
patched together per mineral rather than duplicating the edit.

**Verification:** generated a `.pqi` from all 4 functions and inspected the `RATES`
text directly — 32/32 expected new `SSA`/`mw` lines present (8 minerals × 4
functions), 0 old-style `k = kacid + kneut` lines remaining anywhere in the file.

### CH4/Mtg redox-coupling split by porewater temperature (2026-09-10)

Added a `CH4_redox_override` parameter to all four primary generator functions
(`generate_phreeqc_input`, `generate_phreeqc_input_bulk_kinetic`,
`generate_phreeqc_input_bulk_kinetic_restart`, `generate_phreeqc_input_tidal_kinetic_restart`)
and two new R6 helper methods on `PhreeqcIntegrator`:

- `resolve_ch4_mode(porewater_temp_C, CH4_redox_override)` — decides coupled vs.
  decoupled and returns the species names/recommended database to use.
- `apply_ch4_mode(text, ch4_mode)` — rewrites literal `CH4`/`CH4(g)` tokens to
  `Mtg`/`Mtg(g)` in already-rendered `.pqi` text when decoupled.

**Behavior:** the porewater temperature resolved by each generator (`porewater_temp_C`)
now decides, once per generated file, whether CH4 is treated as redox-active:

| `CH4_redox_override` | Behavior |
|---|---|
| `NULL` (default) | `porewater_temp_C >= 150`: CH4 coupled (as before). `< 150`: decoupled — Mtg used instead. |
| `"coupled"` | CH4 always redox-active, regardless of temperature. |
| `"decoupled"` | Mtg always used instead of CH4, regardless of temperature. |

When decoupled, every `CH4`/`CH4(g)` token in the generated file (EQUILIBRIUM_PHASES
entries, SELECTED_OUTPUT `-molalities`/`-activities`/`-equilibrium_phases`/`-gases`
lists) is rewritten to `Mtg`/`Mtg(g)`. The substitution is applied once to the fully
rendered file text (word-boundary regex, gas form matched first), so it covers both
the static hardcoded mineral lists and the vector-built `secondary_minerals` lists
without needing per-location template changes. `IOM_CH4` (the kinetic organic-matter
channel name) is unaffected — the boundary regex excludes it — but note that channel's
`-formula C 1 H 4.435` still releases plain element `C`, unchanged by this feature; under
a decoupled database that may no longer route to methane the way it does when coupled.

Each generated file's header now carries a comment noting which database to run it
with, e.g.:

    # Recommended database: coreclath_CH4uncoupled.dat (Mtg decoupled -> aqueous/gas species are Mtg/Mtg(g); porewater T = 100.0 C)

using the two new database files `coreclath_CH4coupled.dat` / `coreclath_CH4uncoupled.dat`.

**Verification:** generated `.pqi` files from all 4 functions at T=100 (decoupled),
T=200 (coupled), and with `CH4_redox_override` forcing each mode against the opposite
temperature — confirmed 0 CH4 tokens / Mtg tokens present in decoupled output and vice
versa in coupled output, correct database noted in each header, `IOM_CH4` left intact
in all cases, and an invalid `CH4_redox_override` value raises an error.

### Clathrate stability header note (2026-09-11)

Added two new `PhreeqcIntegrator` helper methods and wired their output into the
`.pqi` header of all four primary generator functions (`generate_phreeqc_input`,
`generate_phreeqc_input_bulk_kinetic`, `generate_phreeqc_input_bulk_kinetic_restart`,
`generate_phreeqc_input_tidal_kinetic_restart`):

- `resolve_clathrate_stability(porewater_temp_C, porewater_pressure_atm)` — for CH4,
  CO2, and H2S, computes the clathrate dissociation pressure at `porewater_temp_C` and
  a `P_stable` flag (`porewater_pressure_atm >= dissociation_pressure_atm`). CO2 and
  H2S also get a `threshold_temperature_C` (the quadruple-point-like break between
  their two-branch dissociation-pressure formulas) and a `T_stable` flag
  (`porewater_temp_C <= threshold_temperature_C`); CH4 has neither, since its
  dissociation-pressure formula is a single branch with no such threshold.
- `clathrate_stability_note(porewater_temp_C, porewater_pressure_atm)` — formats
  `resolve_clathrate_stability()`'s results into `.pqi` header comment lines: one
  pressure-based line per species, plus a temperature-based line **only for species
  that are already `P_stable`** (CH4 never gets one), and a final summary line of
  which species are considered stable (`P_stable`, and `T_stable` where applicable)
  vs. not. Example output:

      # Porewater pressure 261 atm >= dissociation pressure 43.3 atm, so CH4 clathrates may be stable.
      # Porewater pressure 261 atm >= dissociation pressure 24.2 atm, so CO2 clathrates may be stable.
      # Porewater temperature 5 C < threshold temperature 9.85 C, so CO2 clathrates remain stable.
      # Porewater pressure 261 atm >= dissociation pressure 1.7 atm, so H2S clathrates may be stable.
      # Porewater temperature 5 C < threshold temperature 28.5 C, so H2S clathrates remain stable.
      # Clathrates considered stable: CH4, CO2, H2S. Clathrates considered unstable: none.

  Note the species name in this note is subject to the same CH4→Mtg rewrite as
  everything else in the file (see previous entry): in decoupled-mode output it reads
  "Mtg clathrates", not "CH4 clathrates".

`generate_phreeqc_input_bulk_kinetic` and `generate_phreeqc_input_bulk_kinetic_restart`
previously computed `pore_temp_C` from the profile file but never its pressure
counterpart; added `pore_pressure_atm <- mean_layer_pressure_MPa(...) / 0.101325`
(70 atm fallback, mirroring the existing `pore_temp_C`/70 °C pattern) for future use.
The clathrate note in these two (single combined solution) functions is fed
`ocean_pressure_atm` rather than the new `pore_pressure_atm`, matching the pressure
these functions already use as their actual reaction pressure; `generate_phreeqc_input`
and `generate_phreeqc_input_tidal_kinetic_restart` (which have distinct porewater vs.
ocean solutions) feed it `pore_pressure_atm`.

**Known gap (not addressed here):** the `IOM_CH4` kinetic channel's
`-formula C 1 H 4.435` still releases plain element `C`; under a decoupled database
that may not route to `Mtg` the way it speciates through the CH4 redox network when
coupled. Flagged in the previous entry; user is circling back to it separately.

### Clathrate phases added to EQUILIBRIUM_PHASES when stable (2026-09-11)

Clathrate stability is no longer just reported — a clathrate phase (`CH4_hydrate` /
`Mtg_hydrate` matching `ch4_mode`, `CO2_hydrate`, `H2S_hydrate`, per
`coreclath_CH4uncoupled.dat`'s `PHASES` block) is now added to the EQUILIBRIUM_PHASES
block and the matching SELECTED_OUTPUT `-equilibrium_phases` list, for each species
that's stable (`P_stable`, and `T_stable` where applicable) at that block's own
conditions. New `resolve_ch4_mode()` field `hydrate_phase` supplies the correct
CH4- vs Mtg-hydrate name; new method `stable_clathrate_phases(temp_C, pressure_atm,
ch4_mode)` returns the stable phase names, appended into the same
`extra_pw_minerals`/`extra_ocean_minerals`/`secondary_minerals` vectors already used
for temperature-gated minerals, so no separate template wiring was needed.
`generate_phreeqc_input` and `generate_phreeqc_input_tidal_kinetic_restart` check
ocean and porewater independently (their own temperature/pressure); the two
bulk-hydrosphere functions check once against `ocean_pressure_atm`, consistent with
the pressure source already used for their stability note.

### Threshold for mineral depletion changed from 0 -> 1e-9 (2026-09-14)

Just applies to the tidal kinetic restart function. Due to update to the mineral 
dissolution rates (now incorporating SSA, specific surface areas), kinetic dissolution
rates can decrease over time as the mineral is depleted. Thus, minerals tend to 
asymptote towards 0 moles in the rock phase rather than reach 0 directly. This threshold
allows for the STOP flags to correctly catch when a given minerals is at low enough
abundance that we can shift it to the equilibrium phase for the next run. 

### Fixed fractional rate denominator in `last_dissolution_rate` (2026-09-14)

`mineral_dissolution_constants.r`'s `last_dissolution_rate()` computed `rate_frac_s`
(fractional dissolution rate, s⁻¹) as `rate_mol_s / k_curr`, where `k_curr` is the
`k_<mineral>` column's value at the last REACT row — i.e. the moles remaining *after*
that step's `dk_<mineral>` change was already applied. Since `dk_<mineral> = k_curr -
k_prior`, the correct denominator (moles available at the *start* of that step) is
`k_prior = k_curr - dk`, not `k_curr`. Fixed to compute `k_prior` and divide by that.

Verified with a synthetic two-row CSV (`k` 90→80, `dk = -10`): fractional rate now
correctly comes out to `10/90`, not the old `10/80`.

Added `latest_mineral_moles(csv_path)`: returns `k_<mineral>` (remaining moles) at
the same last REACT row `last_dissolution_rate()` reads `dk_<mineral>` from — i.e.
`k_curr`, not the `k_prior` used as the fractional-rate denominator above.

### Fayalite kacid: added missing H+ term (2026-09-15)

Fayalite's acid-mechanism rate constant was missing its `ACT("H+")^1.0` (H+ activity)
term — present in PK04 Table 20 but absent from Table 23 and from Core11, so it never
made it into any of these. Added `* ACT("H+")^1.0` (`* H^1.0` in the R version) to
`kacid` in both `coreclath_CH4coupled.dat` and `coreclath_CH4uncoupled.dat`, all 4
`ENTICES_LF_*.R` RATES blocks (`generate_phreeqc_input`, `generate_phreeqc_input_bulk_kinetic`,
`generate_phreeqc_input_bulk_kinetic_restart`, `generate_phreeqc_input_tidal_kinetic_restart`),
and `mineral_dissolution_constants.r`'s `k_primary_mineral()`.

### Primary mineral assemblage replaced with enstatite-chondrite-like composition (2026-09-17/18)

Replaced the primary mineral set in all 4 generator functions (`generate_phreeqc_input`,
`generate_phreeqc_input_bulk_kinetic`, `generate_phreeqc_input_bulk_kinetic_restart`,
`generate_phreeqc_input_tidal_kinetic_restart`) with a new composition table (initial
moles are mol-per-unit-normalised-rock coefficients, same convention as the old
constants they replace):

| Mineral | Type | mol/kg-rock | mw (g/mol) | SSA (m²/g) | Rate law source |
|---|---|---:|---:|---:|---|
| Enstatite | kinetic | 2.291858 | 100.3725 | 0.1 | unchanged (existing PK04) |
| Forsterite | kinetic | 1.618863 | 140.6715 | 0.1 | unchanged (existing PK04) |
| Troilite | kinetic | 2.400644 | 87.913 | 5 | ex-Pyrrhotite's rate law (hexagonal pyrrhotite, PK04) — no separate PK04 troilite entry |
| Albite | kinetic | 0.319418 | 262.1798 | 5 | unchanged (existing PK04) |
| Diopside | kinetic | 0.276319 | 216.55 | 0.1 | `coreclath_CH4coupled.dat`'s own RATES block (PK04) |
| Anorthite | kinetic | 0.062643 | 278.164 | 5 | unchanged (existing PK04) |
| K-Feldspar | kinetic | 0.020630 | 278.33 | 5 | `coreclath_CH4coupled.dat`'s own RATES block (PK04) |
| Tephroite | kinetic | 0.025938 | 201.96 | 0.1 | **none yet** — `kacid=kneut=0` placeholder (TODO in RATES block); user to supply Mn-olivine constants |
| Fe | equilibrium | 2.380130 | 55.845 | — | n/a (EQUILIBRIUM_PHASES, `dissolve_only`) |
| Ni | equilibrium | 0.283206 | 58.693 | — | n/a (EQUILIBRIUM_PHASES, `dissolve_only`) |
| Lawrencite | equilibrium | 0.011250 | 126.75 | — | n/a (EQUILIBRIUM_PHASES, `dissolve_only`) |

Removed entirely as primaries: Fayalite, Ferrosilite, Pyrrhotite, Magnetite. **Magnetite**
is not gone — it's demoted to an ordinary secondary/precipitate phase (`Magnetite 0 0`,
no initial moles), added to whichever secondary mineral list each function was missing
it from. **Schreibersite** (0.047194 mol/kg-rock in the source table) is deliberately
not implemented yet — it's meant to be incorporated as instantaneous dissolution
products rather than a phase, which is being handled as a separate follow-up.

Fe/Ni/Lawrencite are always-equilibrium and never kinetic-switchable (i.e. not part of
`equil_minerals`/`all_primary` in the two restart functions); in the two functions with
separate porewater/ocean solutions (`generate_phreeqc_input`,
`generate_phreeqc_input_tidal_kinetic_restart`) they're porewater-only (EQUILIBRIUM_PHASES
1), not added to the ocean side (EQUILIBRIUM_PHASES 2).

**Bug found and fixed while wiring this in:** native `Fe` as an EQUILIBRIUM_PHASES name
collides with the pre-existing `-totals Fe` (total dissolved iron) SELECTED_OUTPUT
column — PHREEQC writes both under the literal header `Fe`. In
`generate_phreeqc_input_bulk_kinetic_restart` (`read.delim` without `check.names = FALSE`)
this becomes `Fe`/`Fe.1`; in `generate_phreeqc_input_tidal_kinetic_restart`
(`check.names = FALSE`) both stay literally `Fe`. Both restart functions' recovery
helpers (`get_equil_moles`, `get_ep_moles`) previously took the *first* match, i.e. the
dissolved-iron total, not the mineral's remaining moles. Fixed to take the *last*
match, since `-equilibrium_phases` always comes after `-totals` in these templates.
Verified against a synthetic duplicate-column CSV.

**Known gaps, deliberately left as-is:**
- `generate_phreeqc_input`'s early-exit STOP logic still hardcodes the now-removed
  `Fayalite`/`Pyrrhotite` — left in place as a marker at the user's request, to be
  revisited once a depletion-trigger mineral is chosen for the new composition. Since
  neither name exists in KINETICS any more, `KIN()` returns 0 for both, which makes
  the whole early-exit block permanently inert (never triggers) rather than erroring.
- The two restart functions' `depletion_order` fallback (used only when the caller
  doesn't pass `depletion_mineral` explicitly) was updated to the new mineral names,
  ordered by ascending initial abundance — a provisional, unverified default. Not
  needed when `depletion_mineral` is passed explicitly (the normal usage pattern).

### mineral_dissolution_constants.r: added the 4 new kinetic minerals (2026-09-18)

Added `Troilite`, `Diopside`, `K-Feldspar`, `Tephroite` branches to `k_primary_mineral()`
(same rate laws/mw/SSA as the ENTICES RATES ports above — Troilite reuses Pyrrhotite's,
Tephroite is the `k=0` placeholder), and to the `primary_minerals` lists in
`last_dissolution_rate()`/`latest_mineral_moles()` and the standalone plotting section's
`minerals` vector. The old 8-mineral branches (Fayalite, Ferrosilite, Pyrrhotite,
Magnetite included) are kept as-is, not removed — both assemblages are available as
options in this file; only the ENTICES generator functions were switched over to the
new set exclusively.

Verified: `k_primary_mineral("Troilite", ...)` == `k_primary_mineral("Pyrrhotite", ...)`
(shared rate law), `Tephroite` returns 0, `Diopside`/`K-Feldspar` return sensible nonzero
values, and old minerals are unaffected; `last_dissolution_rate()`/`latest_mineral_moles()`
correctly read `k_Troilite`/`dk_Troilite` etc. from a synthetic CSV.

### generate_phreeqc_input_equil / generate_phreeqc_input_hydrosphere_equil: same mineral update (2026-09-18)

Extended the primary mineral assemblage replacement (see above) to the two
equilibrium-only generator functions. Considered merging them the way
`generate_phreeqc_input_bulk_kinetic` merged its porewater/hydrosphere variants, but
they differ more fundamentally than that precedent (1 vs. 2 solutions, single-step vs.
multi-cycle tidal MIX exchange, 1 vs. 2 REACTION_TEMPERATURE blocks) rather than just
water-mass source — kept them separate.

All 11 minerals are equilibrium-only in both (no kinetics exist in either function),
each `dissolve_only`; `generate_phreeqc_input_equil`'s EQUILIBRIUM_PHASES 1
(porewater) carries all 11, EQUILIBRIUM_PHASES 2 (ocean) carries none, matching how
Fe/Ni/Lawrencite are scoped in the kinetic functions.

Also modernized both functions' secondary-mineral lists, which had drifted out of
sync with the other 4 functions (missing most of the 28 minerals from the "Added 28
new equilibrium minerals" changelog entry; `generate_phreeqc_input_equil` additionally
had orphaned entries — KerogenC128–515, Na2CO3, Pyridine, Siderite, Thermonatrite,
Ice, Monohydrocalcite — not present anywhere else). `generate_phreeqc_input_equil`
now matches `generate_phreeqc_input`'s exact list (partially commented, matching that
sibling); `generate_phreeqc_input_hydrosphere_equil` now matches
`generate_phreeqc_input_bulk_kinetic`'s (fully active). Magnetite demoted to ordinary
secondary in both, consistent with the other 4 functions.

**Verified:** generated `.pqi` from both functions — EQUILIBRIUM_PHASES/SELECTED_OUTPUT/
USER_PUNCH all correctly reflect the new 8 kinetic-named-but-equilibrium + 3
always-equilibrium primaries, modernized secondary lists present, Magnetite correctly
demoted, no leftover references to the old mineral names.

### Renamed the two equilibrium-only generators for disambiguation (2026-09-18)

`generate_phreeqc_input_equil` → `generate_phreeqc_input_tidal_equil` (two solutions,
tidal MIX-cycling time series); `generate_phreeqc_input_hydrosphere_equil` →
`generate_phreeqc_input_bulk_equil` (one solution, single-step) — matching the
`_tidal_*`/`_bulk_*` naming already used by the kinetic generators. Internal helper
`create_phreeqc_template_equil` renamed to `create_phreeqc_template_tidal_equil` to
match. Generated `.pqi`/`.csv` filename patterns are unchanged (still `..._equil...`
/ `..._hydrosphere_equil...`) — only the R function names changed. Updated the one
stale comment reference in `sweep_plots.R`. Verified both renamed functions still
generate valid output.

### Renamed the base kinetic generator for naming consistency (2026-09-18)

`generate_phreeqc_input` → `generate_phreeqc_input_tidal_kinetic`, aligning it with
`generate_phreeqc_input_tidal_kinetic_restart` (its restart counterpart). Internal
helper `create_phreeqc_template` renamed to `create_phreeqc_template_tidal_kinetic`
to match. Updated the two external call sites that reference this function directly:
`check1_regression.R` and `parameter_sweep.R`. Generated `.pqi` filenames are
unchanged — only the R function names changed. Verified the renamed function still
generates valid output.

### generate_phreeqc_input_bulk_kinetic: wired in IOM pools (2026-09-30)

`generate_phreeqc_input_bulk_kinetic` had no IOM support at all — no RATES/KINETICS
entries, no `IOM_*_mol` USER_PUNCH columns — while its own restart counterpart,
`generate_phreeqc_input_bulk_kinetic_restart`, hard-requires those columns whenever
`organic_wt_percent > 0` (PATCH D2's guard). Any run through `bulk_kinetic` with
organics enabled left `bulk_kinetic_restart` with nothing to recover, failing loudly
with "IOM columns missing from prev_csv" rather than silently dropping organics.

Mirrored `generate_phreeqc_input_tidal_kinetic`'s existing wiring exactly: IOM
RATES/KINETICS blocks generated generically from `self$iom_pools` (lazily
initialized via `iom_default_config()` if unset), IOM channel names appended to the
`SELECTED_OUTPUT -kinetics` line, and PATCH D1-style `IOM_*_mol` USER_PUNCH columns
so the restart function can read them back.

**Verified:** generated a `.pqi` with `organic_wt_percent = 1` — now contains all 29
IOM channels' RATES/KINETICS entries and `IOM_*_mol` PUNCH columns (previously
zero). Fed a synthetic CSV carrying those columns into
`generate_phreeqc_input_bulk_kinetic_restart`, which previously threw on this exact
input; it now succeeds.

### generate_phreeqc_input_tidal_kinetic_restart: updated depletion_order (2026-10-01)

Replaced the provisional `depletion_order` fallback (abundance-ordered, unverified,
and referencing the old 4-channel `IOM_CO2`/`IOM_CH4` names that no longer exist in
`iom_pools`) with an order based on each phase's actual dissolution/degradation rate
at 150 C, pH 7, SR=0: primary minerals fastest-to-slowest, Troilite forced last
(its rate depends on ACT(Fe3+), which swings it from fastest to slowest primary
across a plausible activity range — not trusted as an early pick), then the 29 IOM
channels fastest-to-slowest. Tephroite omitted for now (still a kacid=kneut=0
placeholder rate; to be inserted once a real rate is supplied) — it remains a valid
kinetic phase, just not auto-selectable as `depletion_mineral`/`proxy_mineral` until
then; explicit `depletion_mineral = "Tephroite"` still works.
