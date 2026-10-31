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

### Tephroite: real rate law added, superseding the kacid=kneut=0 placeholder (2026-10-06)

Replaced Tephroite's placeholder rate (`kacid=kneut=0`, i.e. never dissolves) with
a real pH- and temperature-dependent rate law: Pincus et al. (2026) ACS Earth
Space Chem 10, 1174-1184, Table S1, citing Casey et al. (1993) GCA 57, 785-793
(synthetic tephroite, anoxic conditions). Acid mechanism (`H^0.47`) with a
pH-dependent activation energy (`Ea = -60.11*ln(-log10(ACT(H+))) + 106.17`,
clamped at 0 — Ea decreases with increasing pH per Casey's Fig 4 data, assumed
0 at pH~6.2). No base branch: Casey's pH>7 data are internally inconsistent
between methods by up to 1.8x and not fittable, so the neutral-mechanism rate
is held flat at its pH=7 value for all pH>7, consistent with how the other
primaries without a base mechanism behave. Anoxic-only: Mn2+ oxidation
"severely retards" dissolution per the source, so this rate law should not be
reused under oxidizing conditions without revisiting.

**Scope:** landed in `generate_phreeqc_input_tidal_kinetic`'s RATES block only
(the function this was tested against). `generate_phreeqc_input_bulk_kinetic`
and both restart functions' `.rate_laws` lists still carry the old
`kacid=kneut=0` placeholder — unlike Schreibersite's rate (below), Tephroite's
RATES block is independently copy-pasted per function rather than sourced from
one shared constant, so each needs updating separately. Also ported 1:1 into
`mineral_dissolution_constants.r`'s `k_primary_mineral()` (verified against the
source's own hand-check: `kacid` at pH 2, 25 C reproduces `10^(-5.4)` exactly).

**Verified:** generated a `.pqi` from `generate_phreeqc_input_tidal_kinetic` and
confirmed the RATES block renders correctly (no stray braces/unescaped quotes
breaking the template); cross-checked `mineral_dissolution_constants.r`'s port
against the same hand-check. Running the generated `.pqi` through PHREEQC at
160 C surfaced a severe convergence problem, unresolved at this time.

### Schreibersite: implemented as a kinetic phase, superseding the equilibrium-products plan (2026-09-30)

Schreibersite (Fe3P) is now a kinetic primary with a placeholder bulk
first-order rate (`k = 1.370e-7 s^-1`, `-formula Fe 3 P 1`) — a flat rate with
no temperature, pH, or surface-area dependence, sourced from a single 1-day lab
H2-evolution measurement at pH 7, 20 C in Ar-purged DI water
(RECONSTRUCTION_STATUS §4/§8 item 1b). This **supersedes** the plan described
earlier in this changelog ("Schreibersite... deliberately not implemented
yet... meant to be incorporated as instantaneous dissolution products rather
than a phase") — that approach was dropped in favor of giving it a real (if
provisional) kinetic rate alongside the other primaries. Unlike Tephroite,
Schreibersite's rate law and KINETICS entry are defined once as shared
constants (`SCHREIBERSITE_RATE`/`SCHREIBERSITE_RATE_BODY`/
`schreibersite_kinetics_entry()`) and referenced identically by all four
kinetic-capable generator functions (`generate_phreeqc_input_tidal_kinetic`,
`generate_phreeqc_input_bulk_kinetic`, and both restart functions) — no
per-function propagation gap. Treat any schreibersite-dependent output as
provisional until a real rate law replaces this placeholder.

### Magnetite: added to EQUILIBRIUM_PHASES wherever it was missing (2026-10-07)

Audited Magnetite (the demoted-ex-primary) across all 6 generator functions and
found two gaps: it was entirely absent from `generate_phreeqc_input_bulk_kinetic_restart`,
and missing from EQUILIBRIUM_PHASES 1 (porewater) in all three tidal functions
(`tidal_kinetic`, `tidal_equil`, `tidal_kinetic_restart`), which previously only
had it on the ocean side. All 6 functions now include it consistently in every
applicable EQUILIBRIUM_PHASES block; verified no duplicate entries resulted.

> **Correction (2026-10-23):** "verified no duplicate entries" holds for the
> EQUILIBRIUM_PHASES blocks, but not for SELECTED_OUTPUT: in
> `generate_phreeqc_input_tidal_kinetic_restart`, Magnetite is in
> `secondary_minerals` and is also `append`ed after Lizardite in `.all_eq_so`,
> so it is listed twice in `-equilibrium_phases`. Fixed in the 2026-10-23
> fixes branch.

### Tephroite: no-precipitation guard added — resolves the convergence problem previously blamed on K-Feldspar/Albite (2026-10-19)

The Tephroite rate law's final line was `rate = k * (1 - SR("Tephroite"))`, with
no guard against `SR > 1`. Once the solution becomes supersaturated with
Tephroite (observed directly in earlier diagnostics: SI crept from -0.09 up
through 0.00 over the course of a run), `(1 - SR)` goes negative, so PHREEQC's
RATES Basic script asks it to `SAVE` a negative `moles` value — i.e. run the
dissolution-calibrated rate expression in reverse as an ad hoc precipitation
rate, which it was never fit for. Added the same guard already used by every
other mineral's rate law: `IF SR("Tephroite") > 1 THEN rate = 0 ELSE rate = k *
(1 - SR("Tephroite"))`.

**This fixes the earlier convergence problem.** With this guard in place, the same class of run converges — pointing to the missing Tephroite guard (not K-Feldspar/Albite's `kbase` terms, previously speculated about) as the actual cause. The extended investigation that followed (small water-rock-ratio restart chains failing even after Albite was moved to an equilibrium phase, traced through several restart generations) was consequently chasing a downstream symptom of this same missing guard, not an independent K-Feldspar/Albite-specific problem.

### generate_phreeqc_input_bulk_kinetic: rec_eqb_minerals auto-promotion (2026-10-19)

Added `rec_eqb_minerals = TRUE` and `rec_eqb_margin = 2`. At the given
`porewater_temp`, computes each of the 8 switchable primaries' kinetic
half-life at a conservatively alkaline reference pH of 11 (via
`k_primary_mineral()` in `mineral_dissolution_constants.r`, now `source()`'d
in so this reuses the exact same rate laws as the RATES blocks rather than a
re-derived copy). Any primary whose half-life, times `rec_eqb_margin`, is
still shorter than one orbital period is moved from `KINETICS` to
`EQUILIBRIUM_PHASES` (`dissolve_only`) instead — treating it as kinetic in that
regime is both numerically fragile (CVODE resolving a near-instantaneous
reaction inside one macro step) and physically no different from equilibrium.
Schreibersite is excluded from the check (no Fe3P phase/log K in either
database, so it can never be an equilibrium primary). The generated `.pqi`
header records which minerals were promoted and why. `rec_eqb_minerals =
FALSE` reproduces the previous unconditional-kinetic behavior exactly.

**Verified:** at 160 C, `rec_eqb_margin=2` promotes Albite and K-Feldspar only
(their half-lives are orders of magnitude under the threshold even doubled);
`rec_eqb_margin=1` additionally promotes Forsterite and Anorthite (borderline
cases just inside 1x but not 2x the orbital period); `rec_eqb_minerals=FALSE`
reproduces the original all-kinetic output exactly. No duplicate USER_PUNCH
line numbers in any case.

> **Correction (2026-10-23):** pH 11 is conservative (it under-promotes) only
> for acid-branch minerals, whose rate is slowest there. For albite and
> K-feldspar, whose base branch dominates, pH 11 gives the *fastest* rate, so
> it over-promotes: at pH 9 albite is ~14x and K-feldspar ~44x slower than
> assumed, which the x2 margin does not cover. See the 2026-10-23 entry that
> changes the reference pH.

### generate_phreeqc_input_bulk_equil: added mode/water_rock_ratio and water_temp (2026-10-19)

This function previously had no equilibrium equivalent to
`generate_phreeqc_input_bulk_kinetic`'s `mode = c("porewater",
"hydrosphere")` split, and no temperature control at all (always hardcoded to
`274.6 - 273.15` ocean-floor temperature). Added both, mirroring
`bulk_kinetic`'s existing pattern: `mode`/`water_rock_ratio` for the
reservoir choice, and `water_temp = "ocean"` (accepting `"ocean"`,
`"profile"`, or a numeric degrees-C value) for temperature, with the same
`T<value>` filename tag `bulk_kinetic` already uses.

**Verified:** omitting `mode` reproduces the exact prior output
(`water = 44.1234526469118`, filename unchanged in structure other than the
new `Tocean` tag); `mode = "porewater"` with `water_rock_ratio` produces the
expected porewater-only solution; `water_temp` accepts `"ocean"`/`"profile"`/
numeric correctly and rejects invalid strings with the same error style as
`bulk_kinetic`'s `porewater_temp`.

### Tephroite rate law and ZERO_MOLES_FLOOR propagated into both restart functions (2026-10-19)

Two fixes to `generate_phreeqc_input_bulk_kinetic_restart` and
`generate_phreeqc_input_tidal_kinetic_restart`'s `.rate_laws` lists:

- **Tephroite's real rate law** (Pincus et al. 2026, see the 2026-10-06 entry
  above) replaces the `kacid=kneut=0` placeholder that entry explicitly left
  in both restart functions. All four kinetic-capable generators now carry
  the same, current Tephroite rate law.
- **`ZERO_MOLES_FLOOR <- 1e-15`** added to both functions' `get_k_moles()`.
  A kinetic mineral dissolved down to floating-point dust (e.g. `1.06e-21`
  mol) is chemically depleted, but writing that literal value into a
  `dissolve_only` EQUILIBRIUM_PHASES entry still makes PHREEQC try to satisfy
  that phase's saturation condition with essentially no mass to work with —
  observed directly under a low-`water_rock_ratio` restart chain (Albite at
  `SI = -4.95`, undersaturated, but `m0 = 1.06e-21`): PHREEQC couldn't
  reconcile "dissolve toward SR=1" with "only this much mass exists," and
  failed to converge on the very first re-equilibration, independent of time
  step size. Clamping to exactly 0 removes the phase from the system
  entirely and converged cleanly in the same scenario. Added to
  `bulk_kinetic_restart` where the failure was observed; added defensively
  (not yet observed as a live problem there) to `tidal_kinetic_restart`,
  which uses the identical `dissolve_only` pattern. Any genuine remaining
  mineral mass in this system is many orders of magnitude above this
  threshold, so it only ever catches floating-point dust.

**Verified:** generated a restart `.pqi` from `bulk_kinetic_restart` with
Tephroite's real rate law rendering correctly; reproduced the exact
`k_Albite = 1.0607e-21` dust value that caused the original failure and
confirmed it now writes as `0.000000e+00` in `EQUILIBRIUM_PHASES`.

### mineral_dissolution_constants.r: last_dissolution_rate()/latest_mineral_moles() accept a minerals override (2026-10-19)

Both functions now take an optional `minerals` argument — pass a custom
vector of mineral names to restrict to those, or omit it (`NULL` default) to
get the original behavior (the full current primary list: Forsterite,
Enstatite, Anorthite, Albite, Troilite, Diopside, K-Feldspar, Tephroite,
IOM_CO2, IOM_CH4, IOM_N, IOM_S). Verified against a synthetic CSV: the
default call still returns `NA` for minerals missing from the data (as
before), and `minerals = c("Forsterite", "Albite")` correctly restricts both
functions' output to just that pair.

### coreclath_CH4coupled.dat and coreclath_CH4uncoupled.dat: fluorine added (2026-10-23)

Both coreclath databases were missing fluorine entirely relative to
`Core11_idealgas_mod_v4.dat`: no `F` master species, no `F-` solution
species, none of the 9 F-bearing aqueous complexes (CaF+, HF, HF2-, MgF+,
MnF+, NaF, SiF6-2, HPO3F-, H2PO3F, PO3F-2), and neither of the two
F-bearing phases (Fluorite, Fluorapatite). Ported all of the above from
`Core11_idealgas_mod_v4.dat` into both coreclath databases, at the same
alphabetical positions Core11 uses.

**Verified:** confirmed via a full phase-name and solution-species diff
that, aside from this fluorine gap and the two databases' intentional
CH4/clathrate coupling differences (`CH4(g)`/`CH4_hydrate` vs. the
decoupled `Mtg_hydrate` scheme), `Core11_idealgas_mod_v4.dat` and both
coreclath databases now carry identical phase and species lists, and that
their `RATES` blocks were already functionally identical (prior diffs were
comment text and `°` character encoding only).

> **Correction (2026-10-23, later):**
> - "identical phase and species lists" was not true of
>   `coreclath_CH4uncoupled.dat`: it lacked the five decoupled-species identity
>   reactions (`Hdg = Hdg`, `Oxg = Oxg`, `Mtg = Mtg`, `Ntg = Ntg`,
>   `H2Sg = H2Sg`), removed in 2e13e34, so it failed to load (11 errors).
>   Restored in the fixes branch; `tests/` now has a load test for every
>   database.
> - There are 10 F-bearing aqueous complexes (the list above names 10), not 9.
> - The coupled file also carried a `CO2 = CO2` identity not in Core11 (dead
>   code; now commented out).

### Fixes branch: databases, tephroite, pressure, rock table, rec_eqb, provenance (2026-10-23)

One commit per item; each commit message has its pass check.

- **coreclath_CH4uncoupled.dat loads again.** The Hdg/Oxg/Mtg/Ntg/H2Sg identity
  reactions removed in 2e13e34 are restored. `tests/db_load_check.R` now requires
  every `databases/*.dat` to load with 0 errors.
- **Databases:** earlier duplicate `CH4(g)`/`C2H6(g)` definitions deleted (PHREEQC
  used the last; output unchanged); dead `CO2 = CO2` commented out in the coupled
  file; RATES comments' mis-encoded degree signs fixed and Fayalite written `^1.0`
  everywhere; critical constants on the five decoupled gases commented out (all
  gases ideal, so they can share a fixed-pressure GAS_PHASE) and the `Ntg(g)`
  `T_c` typo fixed.
- **Core11 retired** to `legacy/`: `coreclath_CH4coupled.dat` contains all of it.
  Tests A4/B2 now use the coupled file and reproduce PHREEQC_TEST_RESULTS.
- **Tephroite:** single acid term at `max(ACT("H+"), 1e-7)` in all four kinetic
  generators and `k_primary_mineral()`; the summed `kneut` (2x at pH 7) is gone.
- **bulk_equil** porewater mode reacts at `pore_pressure_atm`.
- **tidal_kinetic_restart** lists Magnetite once in SELECTED_OUTPUT.
- **ROCK_MOL_PER_KG** regenerated from `reconstruction/rock_mol_per_kg.csv`
  (changes up to 3.1e-4 relative); `tests/rock_table_check.R` passes at 1e-4.
- **rec_eqb_minerals:** default FALSE; reference pH 9 for the base-branch
  feldspars, pH 11 for acid-branch minerals. At 160 C, K-feldspar is no longer
  promoted (its half-life at pH 9 is 2.2e5 s); albite still is.
- **Provenance:** every generated `.pqi` carries `# Target PHREEQC version: 3.8.6`
  and `# Database: <file>`; `run_phreeqc()` (new `phreeqc_tools.R`) warns when the
  executable's version differs.
- **ENTICES_runs.Rmd** has a sandbox header listing the 15 local scripts it needs;
  **legacy/ENTICES_v3.R** stops if sourced.

### Consistency branch: Fe_metal column, equilibrium generators, initial fluid composition (2026-10-25)

One commit per item; each commit message has its pass check.

- **Fe_metal:** Fe metal is punched as `Fe_metal` / `d_Fe_metal` instead of being
  listed in `-equilibrium_phases`, where its column was a second "Fe" next to the
  `-totals` Fe (dissolved iron). Every CSV now has one Fe column. The phase is
  still `Fe` in the databases and EQUILIBRIUM_PHASES. Both restart readers use
  `ep_column()`, which still reads CSVs written before this change.
- **tidal_equil:** new `porewater_temp` ("ocean" default, "profile" or degrees C)
  and `CH4_redox_override`; CH4 mode and database come from `resolve_ch4_mode()`,
  as in the kinetic generators. The default (1.45 C) is therefore now decoupled
  and stamped for `coreclath_CH4uncoupled.dat`. Filenames gain the T tag.
  Secondary minerals are temperature-gated per block (separate commit).
- **bulk_equil:** `temp_gated_minerals(water_temp_C)` appended, as in bulk_kinetic.
- **Equilibrium generators** say in their `.pqi` header and source comment that
  they are inorganic limiting cases (no IOM, no Schreibersite).
- **trace_C_molal** (tidal_kinetic, bulk_kinetic; default NULL): C seed in every
  initial SOLUTION, stamped in the header. 1e-12 lets tidal_kinetic at 25 C
  finish (1e-16 also does); NULL still fails as before.
- **initial_ocean** (tidal_kinetic, tidal_equil; default "pure_water"): named list
  of element totals in mol/kgw for the ocean SOLUTION 2 only, stamped in the
  header. No values chosen.
- **pair_organic_run()** (new `organic_contribution.R`): a run and its
  `organic_wt_percent = 0` blank, and the organics' contribution to ocean DIC,
  CH4, NH4+ and H2 by difference. The blank has 1/0.85 the minerals at 15 wt%.

### Pairing branch: same-rock blank, restart headers, equilibrium clathrates, seed (2026-10-26)

One commit per item; each commit message has its pass check.

- **mineral_basis_wt_percent** (`calculate_phreeqc_parameters()`, default NULL =
  unchanged): minerals scaled by (100 - basis)/100 instead of
  (100 - organic_wt_percent)/100; stamped in the header when set.
  `pair_organic_run()` now uses it, so its blank has the run's mineral inventory
  with the organics simply absent.
- **Restart headers** stamp the parent's `trace_C_molal` and `initial_ocean`
  (`# parent ...:`), read from the parent `.pqi`; values pass down restart chains.
- **Equilibrium clathrates:** tidal_equil's ocean block and bulk_equil's single
  solution get the clathrate phases stable at their conditions, gated by
  `stable_clathrate_phases()` as in the kinetic generators. tidal_equil's
  porewater block gets none.
- **README** recommends `trace_C_molal = 1e-14` (tidal_kinetic and its restart
  both run at 25 C).

### ch4-mode-iom branch: IOM carbon formulas follow the database, carbon re-binning (2026-10-26)

One commit per item; each commit message has its pass check.

- **IOM config tied to the CH4 mode:** `resolve_iom_pools(ch4_mode)` in
  tidal_kinetic, bulk_kinetic and both restarts. With
  `coreclath_CH4uncoupled.dat` (below 150 C) IOM_CH4 is `Mtg 1 H 0.435`
  (redox_mode C decoupled, N and S coupled) and `Mtg` is punched in `-totals`;
  otherwise all coupled and inputs are byte-identical to before. Restarts carry
  and punch the parent's Mtg. A caller-set `iom_pools` is kept.
- **`carbon_by_oxidation_state()`** in `phreeqc_tools.R`: carbon from a PHREEQC
  output or selected-output file re-binned by the actual C oxidation state
  (PHREEQC files formate, acetate, methanol etc. under C(4)). Test:
  `tests/carbon_rebin_check.R` on the N1 leak test output.

### carbon-d branch: C(+4)-only database, carbon_mode, ocean charge balance (2026-10-26)

One commit per item; each commit message has its pass check.

- **`databases/coreclath_C4only.dat`:** `coreclath_CH4uncoupled.dat` with every
  coupled carbon species, phase and master species that is not C(+4) commented
  out (70 species, 16 phases, C(-3)/C(-2)/C(+2); listed in its header). Mtg and
  every C(+4) entry kept; SCN- and urea kept (their C is +4).
  `database_names()` and `restrict_to_database()` in `phreeqc_tools.R`. Test:
  `tests/c4only_db_check.R`.
- **`carbon_mode`** = NULL | `"coupled"` | `"ch4_uncoupled"` | `"c4_only"` and
  **`ch4_cutoff_C`** = 150 on all six generators. NULL is the old temperature
  rule (inputs byte-identical); `CH4_redox_override` still works. `c4_only`
  selects the C4only database, drops the names it lacks from the input, and
  sets IOM_CHn to `Mtg 1 H -2.2` (a 2.2 e- reduction imposed by hand; see
  `iom_module.R`).
- **`charge_balance`** = NULL | `"pH"` | an `initial_ocean` element, for the
  ocean SOLUTION 2 of tidal_kinetic and tidal_equil. NULL (default) is
  unchanged.

### carbon-default branch: c4_only default, restart mode check, Na zeolites (2026-10-26)

One commit per item; each commit message has its pass check.

- **`carbon_mode = "c4_only"` is the default at every temperature.** Coupled
  runs (`"coupled"`) are an upper bound on abiotic methanation: equilibrium
  reduces carbonate with H2 there. `"coupled"` and `"ch4_uncoupled"` stay
  selectable; the old temperature rule is `carbon_mode = "temperature"` (with
  `ch4_cutoff_C`), and `CH4_redox_override` still works. Inputs generated with
  the default before this change are reproduced by `carbon_mode = "temperature"`.
- **Restarts check the parent's carbon mode** (`parent_carbon_mode()`: the
  database on the parent `.pqi`'s `# Database:` line, else its
  `# Recommended database:` line). A mismatch is an error unless the mode is
  given explicitly, then a warning; an unreadable parent mode is a warning.
  Header line `# parent carbon_mode:`.
- **`scripts/make_c4only.py`** regenerates `coreclath_C4only.dat` from
  `coreclath_CH4uncoupled.dat`; `tests/c4only_db_check.R` checks the committed
  file matches it. The database header now names the script (otherwise
  unchanged).
- **Natrolite and Mesolite** (Na zeolites) added to the secondary minerals of
  every generator (porewater and ocean EP blocks, SELECTED_OUTPUT lists, both
  restarts), ungated like Analcime. At 100 °C both were supersaturated in the
  orbit-2 porewater (corrected SI +1.08, +1.30) without being listed.

### docs-1027 branch: PROGRESS 2026-10-27, iom_module.R carbon notes (2026-10-27)

- `docs/PROGRESS.md` replaced with the 2026-10-27 revision (cited paths
  checked; merge state corrected).
- **N3:** `iom_module.R` no longer instructs entering H2 as `Hdg`. The
  docstring and comments describe `c4_only` (the default), with the `Hdg`
  tests kept as history. `attr(cfg, "carbon_decoupling_note")` describes
  `c4_only` and `ch4_uncoupled` (NA when all coupled). `iom_selftest.R`
  checks the `c4_only` formulas and notes.

### gas-ice-docs branch: profile pressure (2026-10-27)

- **Pressures come from the profile in all six generators.** The porewater
  SOLUTION (and a bulk `porewater`-mode solution) is at the profile's mean
  pressure over the permeable layer, 261.2 atm for the Enceladus profile; the
  ocean SOLUTION (and a bulk `hydrosphere`-mode solution) is at the ocean-floor
  pressure, 66.7 atm. Before, five generators wrote 66.7 atm into the
  porewater (`pressure 66.7#261.2` in the tidal ones); only `bulk_equil` used
  the layer pressure.
- **No fallback pressure.** `profile_conditions()` stops if `profile_file`
  does not exist (the generators used to fall back to 70 atm). The Enceladus
  profile is now in `profiles/`, and is the default `profile_file`.
- **`REACTION_PRESSURE 1/2`** in the tidal generators, used by each step
  beside `REACTION_TEMPERATURE` (MIX averages pressure by mass, so the
  porewater would otherwise drift towards the ocean pressure).
- **Header stamps** `# Porewater pressure:` and `# Ocean pressure:`, beside the
  version and database stamps. **Restarts** compare them with the parent's
  (`check_restart_pressures()`): a mismatch is an error, or a warning with
  `allow_pressure_change = TRUE`; an unstamped parent is a warning. Header line
  `# parent pressures:`.
- `tests/pressure_check.R` checks all of the above (inputs only).
- **Initial fluids at 25 C.** The initial porewater and ocean SOLUTIONs are
  defined at 25 C with pe 4 (no oxidant or reductant) and brought to the run
  temperature by REACTION_TEMPERATURE (added to `bulk_kinetic` and
  `bulk_equil`, which had none). Defined at the run temperature, pe 4 put
  O(0) into the porewater: 1.1e-22 molal at 100 C, 9.0e-12 at 160 C, 4.9 at
  250 C (0.77 mol O in the 0.157 kg porewater), and did not converge at
  300 C at any pressure. Restarts are unchanged (they start from the parent's
  reacted state). `tests/initial_fluid_check.R` runs each generator's initial
  porewater at 25, 160 and 300 C.
- `profiles/README.md`: provenance of the Enceladus profile (PlanetProfile).
- **Mineral molar volumes and rock density.** `scripts/export_mineral_volumes.py` writes
  `data/mineral_molar_volumes.csv`: a volume for every non-gas phase of the coreclath databases plus
  Schreibersite, by the PI's source priority (Core11, then Green et al. 2025 ds636, then Perple_X
  hp636, then PI-supplied values). Every production phase is covered; Fluorite and Fluorapatite are
  unfilled. `bulk_density()` (`phreeqc_tools.R`) gives rho_solid, rho_bulk and a per-phase breakdown
  through a run; the IOM density (1.4 g/cm3) is UNSOURCED. Every generator now punches `Rho_soln`
  (PHREEQC's RHO); `FE_METAL_PUNCH_*` are renamed `SHARED_PUNCH_*`. `tests/bulk_density_check.R`.
- **Mineral molar volumes stripped from the databases** (database change for Lucas to confirm).
  `scripts/strip_mineral_vm.py` comments out (`#vm# `) the `-Vm`/`Vm` line of every non-gas PHASES
  entry in the coupled and uncoupled databases (384 each, 1 unhyphenated); `make_c4only.py` applies
  the same pass (flag `strip_mineral_vm`, `--no-strip-mineral-vm`). With aqueous volumes inactive
  under the LLNL model, the term shifted each mineral's log K by its whole Vm·ΔP/(ln10·RT)
  (antigorite −4.16 at 261 atm). Mineral equilibria are now at the database reference pressure
  (~0.1–0.2 log unit omission at 260 bar). Check: the 160 °C production case at 66.7 vs 261.2 atm
  gives byte-identical CSVs, 265 warnings each (was 350 / 641), 116 s (was 481 / 672 s).
- **Restarts carry every element, at full precision** (mass-balance fix). Both restart paths
  (`tidal_kinetic_restart`, `bulk_kinetic_restart`) built their SOLUTIONs from a fixed list,
  Al C Ca Fe Mg N Na P S Si (+ Mtg), from a CSV with 4-5 significant figures: K (K-feldspar), Mn
  (tephroite), Cl (lawrencite) and Ni were lost at every restart, since no generator punched them
  (`#K #Mn #Cl` were commented out of `-totals`; the equilibrium generators did not punch Mtg either).
  Now every generator punches all 36 database elements (`totals_block()`, pseudo-elements included)
  with `-high_precision true`; restarts take every element the parent punched
  (`restart_elements()`, warning about any ENTICES element a pre-2026-10-28 parent lacks) and write
  pH, pe, water, element totals, kinetic and equilibrium-phase amounts with 15 digits. MIX fractions
  are written with 15 digits too (at 7 they did not sum to 1, losing up to ~5e-8 of the porewater per
  cycle). `ep_column()` takes the last of a duplicated name (Ni is now both a total and a phase).
  Note: `-high_precision` also tightens PHREEQC's mole-balance convergence from 1e-8 to 1e-12 of each
  element's total (KNOBS, Parkhurst & Appelo 2013, p. 112).
  `system_totals()` and `formula_elements()` (phreeqc_tools.R) sum each element over solution,
  solids, IOM and gas. `tests/restart_closure_check.R`: 2 orbits vs 1 + restart 1, tidal and bulk,
  160 C, 15 wt% IOM: every element total agrees within 1e-9 relative or the solver tolerance (max
  4.8e-13 for the rock-dominated totals).

### gas-ice-docs branch, Stage 2: fixed-pressure gas phases (2026-10-28)

- **Gases leave EQUILIBRIUM_PHASES.** Every `X(g) 0 0` equilibrium phase (an infinite reservoir at
  1 atm, which capped dissolved gases at their 1-atm solubility and kept what exsolved) is gone from
  all six generators and both restarts. Each reservoir has a fixed-pressure GAS_PHASE instead,
  starting empty: porewater GAS_PHASE 1 at the porewater pressure, ocean GAS_PHASE 2 at the ocean
  pressure (tidal generators; never emptied), and one GAS_PHASE 1 at the solution's pressure in the
  bulk generators. Components (`gas_components()`): the gases the EP lists carried, filtered by
  carbon mode (`c4_only`: Mtg(g) CO2(g) H2(g) H2S(g) N2(g) NH3(g) SO2(g)). Clathrates stay
  equilibrium phases.
- **`pore_gas_fate`** (`tidal_kinetic`, `tidal_equil`, `tidal_kinetic_restart`): `"retain"`
  (default) keeps the pore bubble with the porewater; `"to_ocean"` moves its moles into the ocean
  step each orbit (USER_PUNCH `PUT` + `KINETICS 3` at a constant rate over 1 s; `gas_transfer_block()`)
  and starts the porewater's next orbit with an empty gas phase. Stamped in the header; a restart
  warns if it differs from its parent's.
- **Peng-Robinson for every gas** (`scripts/gas_eos.py`; `--ideal` reverses it): ordinary gases get
  their own commented constants back, the decoupled gases their analogues' (Mtg←CH4, Hdg←H2,
  H2Sg←H2S, Ntg←N2, Oxg←O2, as in phreeqc.dat 3.8.6); H2's Omega is −0.225 (phreeqc.dat), not the
  entry's +0.225. H2 at 160 °C, 261.2 atm: PR dissolves 8.3% more at saturation and its bubble's
  molar volume is 8.7% larger than ideal (`tests/gas_eos_check.R`).
- **Diagnostics punched every row:** `Vgas_L`, `P_gas_atm` (0 when there is no bubble),
  `Vgas_Vpore` (gas volume over the pore volume, from the porewater mass at 1 kg/L), and each
  component's moles (`g_` columns, `-gases`). `H_total_mol`, `O_total_mol` (total H and O with water;
  the `-totals` H and O leave water out) for mass balance.
- **Convergence:** `-high_precision` (item 4) also sets PHREEQC's mole-balance criterion to 1e-12;
  with a gas phase the metal-pulse step then ran for > 15 min. A KNOBS after SELECTED_OUTPUT puts
  back the default 1e-8 (the last setting read wins: PHREEQC read.cpp; the manual says "the smaller
  of"). Gas-phase runs are still ~10x slower than before (2 orbits at 160 °C: ~140-280 s).
- **Restarts** rebuild both gas phases from the parent's last cycle (`GAS_PHASE_MODIFY`, moles to
  15 digits; `gas_restore_block()`).
- `organic_contribution()` also reports the ocean gas phase's components (mol).
- `bulk_density()` includes the pore gas bubble: rho_pore = (1 - s_g) rho_fluid + s_g rho_gas, gas
  saturation s_g = min(Vgas/Vpore, 1).
- `ocean_gas_loss_fraction_per_orbit` was not added: removing a fraction of a gas phase's moles at
  run time has no simple PHREEQC route (see PROGRESS).
- `tests/gas_mass_balance_check.R`: whole-system element totals per orbit, both fates (≤ 1e-11);
  the no-bubble identity (25 °C, no metal).

### gas-ice-docs branch, Stage 3: cometary-ice initial ocean (option) (2026-10-28)

- **`initial_ocean = cometary_ice_ocean(v = 0.02, include_NH3 = FALSE, include_O2 = TRUE)`**
  (tidal_kinetic, tidal_equil). Miller et al. (2026) Table 5: comet 67P ice (Rubin et al. 2019) with
  ices more volatile than NH3 depleted by v = 0.02 (from Titan's N budget, not an Enceladus
  constraint). The default ocean stays pure water. The ices react into the pure-water ocean at
  ocean T and P in a simulation of their own before the first orbit (REACTION, not SOLUTION totals:
  the mixture is out of redox equilibrium), with the ocean's gas phase, which is saved;
  `<base>_initial_ocean.csv` records the result. Methane and ethane enter as Mtg (`Mtg 2 H -2` for
  ethane) unless carbon is coupled; N2 as Ntg (and Ntg(g) joins the gas components); CO as CO2
  minus O (REACTION reads "CO" as the cobalt phase Co).
- Result (c4_only, v = 0.02): coupled C 0.0866 molal, all C(+4); Mtg 0.0102; pH 1.85; H2 ~5e-14;
  S(+6) 0.0094, S(-2) 0.0028 molal; no bubble at 66.7 atm. The hand estimate (0.0086 / 0.0036)
  left out the ethane term (2 e- per C2H6 drawn from the coupled system). Without O2: no sulfate,
  H2 0.031 molal, pH 3.80. With NH3: 7.44e-3 molal N(-3), pH 2.08. `tests/cometary_ocean_check.R`.
- `include_NH3` defaults to FALSE (Table 5 as printed); pending the PI.
- A cometary ocean with a carbon-free porewater (no `trace_C_molal`, no IOM) left the first porewater
  step running for > 20 min without converging (100 C, 1 orbit); with `trace_C_molal = 1e-14` it
  completes. `tidal_kinetic` warns in that case. The cause is not understood.

### gas-ice-docs branch, Stage 5: documentation (2026-10-28)

- `README.md` rewritten for regenerating production: what ENTICES is, requirements, step-by-step
  running (self-tests first, outputs and their columns, pairs, sweeps, restarts, runtimes, common
  failures), worked examples, a parameter table (default, units, source, status), model choices
  with sources, interpreting outputs (rock density, gas, carbon), open decisions, repository map,
  references (checked against the PDFs; unverifiable ones marked `[citation needed]`).
- `ENTICES_runs.Rmd` is now the worked examples (25 C smoke test, Lucas's production case, a 100 C
  pair, the cometary ocean, the pore-gas bracket). Lucas's notebook moved to
  `examples/ENTICES_runs_lucas_2026-10.Rmd`.
- Roxygen-style headers on the six generators (purpose, every argument with units, default and
  justification, return value, example). Surface areas, Fe/Ni instantaneous corrosion and the tidal
  stress amplitude are tagged UNSOURCED in the code; the ocean temperature carries its source.
- `profiles/README.md`: the profile's provenance and the verified Styczinski et al. (2023) citation.
- `tidal_equil` gains `trace_C_molal` (default NULL): with a cometary ocean and no seed every case
  failed ("Cl has not converged"); with 1e-14 they run, except `to_ocean` + `c4_only`, which crashes
  PHREEQC (segmentation fault in the ocean step).

### gas-ice-docs branch, Stage 6: tests and the generator matrix (2026-10-28)

- All 14 self-tests pass (`tests/README.md`), including the long `restart_closure_check.R` and
  `gas_mass_balance_check.R`, on the final code.
- `docs/MATRIX_2026-10-28.md`: 102 runs (every generator x carbon mode x pore_gas_fate x initial
  ocean x 25/100/160 C, 1 orbit, plus restarts): 87 converge; the failures are the unseeded cometary
  `tidal_equil` cases (rerun with the seed: converge) and `tidal_equil` + `to_ocean` + cometary +
  `c4_only`, which crashes PHREEQC 3.8.6.

### gas-ice-docs branch: decisions and sources, 2026-10-29

- **`pore_gas_fate` defaults to `"to_ocean"`** (PI): `"retain"` puts ~2.2 pore volumes of gas in the
  pores (an upper bound on retention); carbon attribution is the same in both (to 3 significant
  figures); `"to_ocean"` without plume loss is an upper bound on ocean H2. `"retain_to_pore_volume"`
  is planned after D6.
- **IOM density from Delarue et al. (2026) Eq. 10** (PI): `bulk_density(iom_density = "delarue2026")`,
  the default, evaluates rho = 1472 - 345 (H/C) + 94 / ((H/C) + 0.122) kg/m3 at the residual IOM's
  H/C on each row (bulk IOM less the channels' releases; 1.320 g/cm3 at t = 0); a number fixes it.
  Replaces the UNSOURCED 1.4 g/cm3. Lucas's 160 C case: rho_solid 2.977 -> 2.920 g/cm3 at t = 0
  (-1.9%), 2.806 -> 2.759 at orbit 12; rho_bulk 2.344 -> 2.305 at t = 0. Compression left out
  (K0 = 10.2 GPa: 0.25% at 26 MPa). Summary gains `iom_HC` and `iom_density`.
- **`bulk_density()` flags Vgas > Vpore** (PI): new `Vgas_Vpore` and `gas_exceeds_pores` summary
  columns and a warning; gas beyond the pore volume is not counted (s_g capped at 1, as before).
- **Core density sources** (PI): README section 7 compares rho_bulk with Bagheri et al. (2025) Table 5
  (core 2160-2400 kg/m3, 183-213 km): 2.305 g/cm3 at t = 0 and 2.255 at orbit 12 (`"to_ocean"`) are
  inside, 1.880 (`"retain"`) below. Iess et al. (2014)'s ~2.4 g/cm3 is an example structure for
  C/MR2 = 0.335, not a measurement. The config's radius, mass, core radius and core density are the
  PlanetProfile profile's; its [citation needed] now says so (neither paper gives those values).
- **Citations checked against the PDFs** (PI): Casey, Hochella & Westrich (1993) GCA 57 (tephroite;
  not the Nature 366 chain-silicate paper); Pincus, Peters, Yi & Kermani (2026); Glein & Waite (2020);
  White & Brantley (2003). The tephroite Ea fit (-60.11 ln(pH) + 106.17) is from Casey Fig. 4's
  labels (54, 26, 22, 4 kJ/mol; 0 at pH 6.2), not the appendix rate table; the REM now says it
  reaches 0 at pH 5.85. README section 7 compares the oceans' carbon with Glein & Waite (2020).
- **KNOBS `-logfile` off by default** (PI): every generator takes `logfile = FALSE`; `TRUE` writes
  `phreeqc.log` (258-346 MB per 12-orbit production run; 11.9 GB from one stalled run).
  tests/pressure_check.R checks all six generators.
- **Flow driving documented as a decision (D12) and a sweep design** (PI; not run, N6): README section 6
  cites Liao et al. (2020), Choblet et al. (2017), Zandanel et al. (2021) and Fisher et al. (2024)
  as range constraints with their limitations (checked against the PDFs), tabulates the flow
  equivalents (mean Darcy flux, exchange per orbit, displacement, circulation time) over k and
  dP, and flags production k = 1e-7 m2 (top of Liao's range; 243 yr ocean circulation against
  Choblet's 25-250 Myr). Open for the PI: oscillatory flow as net exchange, and the porosity
  factor in the Simulator's flow rate.
- **Matrix follow-up** (PI): 87 of the 90 distinct cases converge counting the seeded `tidal_equil`
  reruns (the 102 rows include 12 reruns); the 3 failures are the PHREEQC 3.8.6 segfault, now
  reduced to a 20-line reproducer in `tests/known_failures/` (crashes with the stock llnl.dat). The
  unseeded cases' 0.019125 mol Cl is the porewater's lawrencite (2 x 0.0095625). Correction: the
  full cometary ocean with a carbon-free porewater is slow, not stuck (1,459 s unseeded vs 394 s
  seeded; the earlier "> 20 min" was a timeout); a CO2-only ocean does stall unseeded (no progress
  in 3,000 s) and completes in 267 s seeded. `docs/MATRIX_2026-10-28.md`.

### exchange-physics branch, 2026-10-29

- **Darcy flow no longer multiplied by porosity** (PI): `flow_velocity` is the Darcy flux, per unit
  total area, and `flow_area` is the core's whole surface, so the flow is q * A. The extra factor
  made the exchange 1/porosity (3.125x) too small. Lucas's production case (k 1e-7, full core):
  exchange per orbit 3.342e-5 -> 1.044e-4 of the porewater; circulation time 243 -> 78 yr; auto
  grouping factor 242 -> 77. Every tidal run's exchange changes; earlier results need rerunning.
- **Exchange per half-orbit = min(poroelastic, Darcy, pore volume)** (PI): replaces the Darcy-only
  exchange (`Simulator$calculate_exchange()`, `sim$exchange`). The poroelastic term is
  alpha * eps_v * V_layer with eps0 from Liao et al. (2020) (Fig. 2 caption form, reproduces their Fig.
  2a; the 3/5 in their text does not) and a degree-2 dilatation factor (0.115 for the full core); core
  rigidity and viscosity are Bagheri et al. (2025) placeholders (17 GPa, 1e12 Pa s); biot_alpha 0.2.
  Per orbit the halves compound, so x <= 1. At production k the poroelastic limit binds: exchange
  1.25e-6 of the pore volume per half-orbit, circulation 3,250 yr (78 yr Darcy-only); Darcy binds
  below k ~ 2.4e-9. The pore-volume limit cannot bind at real porosities. Header line and metrics
  report the binding limit. `tests/exchange_check.R`.
- **Known failure (open):** at 25 C the ocean step fails to converge when the porewater fraction mixed
  in is ~1e-6 (3.3e-6 fails, 3.3e-7 and 3.3e-5 run; the pre-branch code fails the same way at k
  1e-8). The new production exchange (2.5e-6 per orbit) falls in this window: the 25 C smoke run
  now gives 18 errors; 100 and 160 C (1 orbit) run with 0 errors.
- **Darcy term over the poroelastic skin depth** (PI): the pressure gradient acts over L_eff =
  min(layer_thickness, delta), delta = sqrt(c/omega) with Liao et al. (2020) Table A1's c (K_f 2.2 GPa,
  K_s 10 GPa; new config fields). delta: 28 m (k 1e-14), 0.9 km (1e-11), 88 km (1e-7). The Darcy volume
  now grows as sqrt(k); Darcy/poro crossover 2.4e-9 -> 1.2e-11 m2. Circulation: 1.1e5 yr (1e-14), 3.5e4
  (1e-13), 1.1e4 (1e-12), 3,530 (1e-11); production unchanged (poro binds, 3,250 yr). Cross-check:
  S dP delta A = Darcy / 2.02 at every k. Note alpha 0.2 < phi 0.32 (not physical for a porous solid).
- **25 C ocean-step convergence fixed** (PI): KNOBS -diagonal_scale true around every ocean step (MIX 4)
  of tidal_kinetic, tidal_equil and the tidal restart, off after it. Without it the 25 C ocean step
  failed at scattered porewater fractions (3e-8, 2e-6, 5e-6 of 12 tested from 1e-8 to 1e-3, and the
  production 2.5e-6); with it all converge. A 160 C orbit is unchanged (1e-8 relative); set globally
  it moved the porewater m_H2S by 0.9%, so it is ocean-step only. The 25 C 1-orbit and 2-orbit runs
  give 0 errors. Trace seeds were not tried (not needed).
- **The reacting box is the skin** (PI): `Simulator$new(reactive_depth = "skin")` (default) or "layer" (old).
  `sim$reactive_thickness()` = min(L, delta) sets the porewater T and P (profile means over it), the rock
  mass and pore volume, and the exchange fraction (relative to the skin's pore volume). Poro binds at
  every k: x_half 2.1e-6 (k 1e-14 to 1e-11, skin 28-883 m, 1.4 C, 66.7 atm; circulation 4.5e6-1.4e5 yr)
  and 1.4e-6 at production (88 km, 36.7 C, 177.3 atm; 3,410 yr). Tidal pumping alone reaches only cold
  near-seafloor rock at low k. Header "Reacting box" line; pressure stamp "(profile mean over the
  reacting box)". pressure_check and initial_fluid_check use reactive_depth = "layer" for their
  full-layer pressures; exchange_check gains skin-box checks.
- **Porosity and the Biot coefficient** (PI): `Simulator$new(phi0 = 0.32, phi_c = 0.5, alpha = NULL)`
  (`porosity` still accepted). alpha = min(1, phi/phi_c) from the critical-porosity model (Nur et al. 1998,
  The Leading Edge 17, 357-362; cited via Dvorkin & Nur 2000); `alpha =` overrides, below phi is an error.
  The config's free `biot_alpha` (0.2, below phi) is removed. The box's solid mass is V (rho_core -
  phi rho_water). `phi_from_core_density()`; `bulk_density()` flags (`core_density_ok`) and warns when
  the final rho_bulk is outside Bagheri's 2.16-2.40 g/cm3. At phi0 0.32 (alpha 0.64): skin 23 m-74 km;
  production circulation 1,117 yr (skin box); x_half 6.66e-6 at low k for any phi0 (= eps_v/phi_c).
  Bagheri gives no frame/grain moduli (one fixed core bulk modulus), so alpha is not computed from them.
- **Porosity that evolves with the chemistry** (PI): `run_tidal_evolution()` (`porosity_evolution.R`),
  rigid frame, phi = 1 - V_solid/V_total, updated between chunks (k via `k_of_phi`, alpha, storage,
  delta, box, exchange); chunks rerun at half length if phi changes > 1%; phi_min 0.01 and `clogged`.
  Option, default FALSE. `bulk_density(origin_pqi =)` carries the IOM inventory across restarts.
  `tidal_kinetic_restart(stop_on_depletion = TRUE)` (the USER_PRINT STOP once depletion_mineral is
  exhausted, which PHREEQC reports as 2 errors) can now be turned off. 160 C, 24 orbits: phi 0.32 ->
  0.219 (15 wt% IOM) / 0.165 (none), k 1e-7 -> 2.5e-8 / 9.1e-9 (Kozeny-Carman), delta 74 -> 45 / 32 km,
  circulation 1,117 -> 2,014 / 3,334 yr; poro binds throughout; the core closes, no clogging. Water:
  PHREEQC's porewater runs 37% (IOM) / 20% (none) short of the solids' pore volume after orbit 1.
  `tests/porosity_evolution_check.R`.
- **The tidal pore pressure is derived from the strain** (PI, 10-30): dP = alpha eps_v M, M = 1/S
  (Simulator$tidal_pressure(); tidal_strain()). config `max_tidal_stress` defaults to NULL; a number
  overrides it. Placeholders (alpha 0.64): pole 6.2e4 Pa; box mean 1.2e4 (thin skin), 8.6e3 (k 1e-7
  skin), 7.2e3 Pa (layer). At Liao et al.'s Fig. B3 case the formula gives 1.1e3 Pa against their
  ~3e4 (their ocean-loading term is not in it). Storage volume = strain volume in the skin, Darcy =
  2.2x storage: the strain binds at every k in the skin box. Header and `sim$exchange` report dP.
- **Strain factor follows the box** (PI, 10-30): eps_v/eps0 = (mean |P2| / 2 = 0.192) x the radial mean of (r/R)^2
  over the reacting box: 0.192 for a thin skin (0.138 for the 74 km skin at k 1e-7), 0.115 for the full layer.
  Already how tidal_strain() evaluates it; documented and tested (exchange_check 3d).
- **phi_c = 0.4 by default** (PI, 10-30): Dvorkin & Nur (2000)'s table (36-40% for sorted grain packs and
  sandstones) after Nur et al. (1998); phi >= phi_c gives alpha = 1 (a suspension, not a frame). At phi0
  0.32, alpha 0.80: skin box x_half 8.3e-6 at low k, production 6.1e-6 (69 km skin, 31 C, 156 atm),
  circulation 912 yr. README tables recomputed with dP from the strain.
- **Refill of the pores from the ocean** (PI, 10-30): `tidal_kinetic_restart(refill_kg =)` mixes ocean water
  into the porewater (negative: the reverse) before the first orbit, without output rows; run_tidal_evolution()
  refills to m_pw(0) phi(t)/phi0 between chunks and reports it next to the tidal exchange. After orbit 1 it
  moves 0.033-0.044 kg (normalised) against 1.9e-6 kg of tidal exchange per orbit. restart_closure_check:
  whole-system totals conserved for +0.02 / -0.01 kg refills. `porosity_evolution` defaults to TRUE in
  run_tidal_evolution() (FALSE in the quick self-test).
- **The 10-29 porosity test overrode the temperature** (porewater_temp = 160 on the skin box, whose own T is
  31 C at k 1e-7). Rerun at the skin's temperature (new defaults, refill on): phi 0.32 -> 0.279 (15 wt% IOM)
  / 0.259 (none) in 24 orbits; at 160 C -> 0.222 / 0.169. The core closes, no clogging, poro binds.
- **Ocean tidal loading at the seafloor** (PI, 10-31): Simulator$tidal_load_pressure(), rho_w g_c |u_ocean - u_core|
  with u_r = h2 U2/g at the pole (Bagheri et al. 2025 Eqs. 1-2; surface |h2| 0.0459, lag 2.95 deg, Table 3
  dissipative core; core h2 = (5/3) k2* of the uniform Maxwell core). Pole 213 Pa (108-459 over eta), box
  mean 25-41 Pa; added to the strain dP as an amplitude (upper bound) and its stored volume to the poro
  limit. 0.3% of the strain term at every k and eta. `ocean_load = FALSE` drops it. It does not explain
  Liao's ~3e4 Pa (README corrected). Flag: the core h2 behind eps0 (0.29) is 6x the whole-body h2 (0.046).
- **Water drawn in by the refill** (PI, 10-31): `refill_summary()` (porosity_evolution.R) per chunk and
  cumulatively, scaled to the whole core by the box's rock mass (new `layer_mass_kg` column), as a fraction
  of the profile's ocean (`profile_ocean_mass_kg()`, 2.13e19 kg), and the ocean's concentration factor.
  24-orbit runs: 1.27e18-1.82e18 kg, 6.0-8.5% of the ocean, almost all in orbit 1; concentration factor
  1.06-1.09 (model's own ocean 1.08-1.12).
- **strain_scale** (PI, 10-31): config factor on eps0 (default 1) for the N6 sweep; the core h2 bracket is
  strain_scale 0.046/0.29 = 0.16 to 1, documented as the largest uncertainty in the exchange. At 0.16 the
  thin-skin exchange is 1.3e-6 per half-orbit (from 8.3e-6). exchange_check 3f.
