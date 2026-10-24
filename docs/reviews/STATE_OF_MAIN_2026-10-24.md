# State of main: for review

Date: 2026-10-24
From: the VS Code (Claude Code) session.
For: the claude.ai review session and the PI. Lucas should see §1 and §4.

**Code state described:** `main` @ `bd6c6dc` (PRs #4 and #5 merged) plus the `fixes` branch (19 commits, pushed, PR not yet merged). All Stage 3–5 measurements were made on `fixes` @ `0979d8d`.

**Tools:** PHREEQC 3.8.6 and R 4.5.2 on Windows. Every number below is **measured** (PHREEQC/R output, a diff or a grep) unless it is marked *inferred*.

---

## 1. What changed

### 1.1 PR #4 (Lucas, merged)

PR #4 added fluorine to both coreclath databases:
- the `F` master species;
- `F-` plus 10 complexes;
- Fluorite and Fluorapatite.

These are identical to Core11's entries. PR #4 did **not** restore the uncoupled database's decoupled species, so that file still failed to load (11 errors).

### 1.2 PR #5 (cleanup, merged)

PR #5 reorganised the repo:

| Path | Contents |
|---|---|
| `databases/`, `tests/`, `examples/` | databases, tests and example runs |
| `docs/` (`reviews/`, `archive/`) | design docs, review exchanges, superseded docs |
| `legacy/`, `scripts/` | retired files and standalone scripts |
| `reconstruction/` | the reconstruction module; its Rmd now writes `rock_mol_per_kg.csv` |

It also added `tests/rock_table_check.R`, a README and per-folder READMEs.

### 1.3 `fixes` branch (one commit per item; each commit message records its pass check)

| Item | Change | Pass check (measured) |
|---|---|---|
| Uncoupled database | `Hdg`/`Oxg`/`Mtg`/`Ntg`/`H2Sg` identity reactions restored (removed in `2e13e34`) | Loads with 0 errors; 100 °C decoupled `tidal_kinetic` runs 2 orbits, 0 errors (635 s) |
| Duplicate gases | Earlier `CH4(g)` / `C2H6(g)` definitions deleted (PHREEQC uses the last) | Byte-identical output before/after with C(4), C(−4), C(−3) at 25 and 150 °C |
| `CO2 = CO2` | Commented out in the coupled file (dead) | Byte-identical carbonate speciation |
| Encoding | RATES comments' mis-encoded `°` fixed; Fayalite `^1.0` everywhere | All 58 RATES identical across files |
| Core11 | Retired to `legacy/`; tests A4/B2 now on `coreclath_CH4coupled.dat` | Every number in `PHREEQC_TEST_RESULTS` reproduced; only difference is 3 extra hydrate SI lines |
| CHANGELOG | Correction notes under the 10-07, 10-19 and 10-23 entries | Originals untouched |
| Load test | `tests/db_load_check.R` + `phreeqc_tools.R` | Passes on repo; fails (11 errors) on the `2e13e34` file |
| Decoupled gases | `-T_c/-P_c/-Omega` commented out on all five; `Ntg(g)` typo fixed | 100 °C fixed-pressure GAS_PHASE with H₂(g)+Mtg(g), and with all five decoupled gases, runs (before: "cannot mix ideal and Peng_Robinson") |
| Tephroite | One acid term at `max(ACT(H+), 1e-7)` in all four kinetic generators and `k_primary_mineral()` | k/M (s⁻¹), hand = PHREEQC = R: 25 °C pH 5 3.128e-6, pH 7/11 3.591e-7; 160 °C pH 5 1.023e-5, pH 7/9/11 3.591e-7 (old pH 7 value 7.18e-7) |
| `bulk_equil` | Porewater mode reacts at `pore_pressure_atm` | With the profile file: porewater 261.2 atm, hydrosphere 66.7 atm; both run |
| `tidal_kinetic_restart` | Magnetite listed once in SELECTED_OUTPUT | Once per block; restart runs; one CSV column |
| Rock table | `ROCK_MOL_PER_KG` regenerated from the reconstruction CSV (lawrencite Fe still subtracted) | `rock_table_check` passes at 1e-4 (worst 1.1e-5); changes ≤ 3.1e-4 relative |
| `rec_eqb_minerals` | Default `FALSE`; reference pH 9 for albite and K-feldspar, pH 11 for the rest | At 160 °C K-feldspar is no longer promoted (t½ 2.2e5 s at pH 9; was 5.1e3 s at pH 11); albite still is (1974 s) |
| `ENTICES_runs.Rmd` | Sandbox header listing the 15 local scripts it sources | Matches grep |
| Provenance | Every `.pqi` stamped `# Target PHREEQC version: 3.8.6` and `# Database: <file>`; `run_phreeqc()` warns on a mismatch | All six generators stamped; a 3.5.0-stamped copy warns |
| `ENTICES_v3.R` | Unreachable: stale comment removed; `source()` stops unless explicitly opted in | 0 `source()` calls and 0 code mentions outside `legacy/` |

Housekeeping, also on the branch:
- **Stray outputs removed:** a run CSV and a `.sel` file that my own 10-22 test runs had committed in `3ae539a`. Root-level run outputs are now git-ignored.
- **Comment fix:** a doc comment I had misplaced, now back above `resolve_ch4_mode()`.

---

## 2. What was verified

### 2.1 Self-tests (measured)

| Check | Result |
|---|---|
| `tests/iom_selftest.R` | 0 failures |
| `tests/rock_table_check.R` | 0 failures |
| `tests/db_load_check.R` | 0 failures |
| Test inputs A4, B2, B3 (from the repo root) | 0 errors each |

### 2.2 Verification matrix

**Setup:** 2 orbits, 15 wt% IOM, `k = 1e-7`, full-core layer, with the PlanetProfile file. CH₄ is coupled at 160 °C and decoupled at 25 and 100 °C (the defaults).

**Warnings.** Every run gives about 17–100 "Element X is contained in <phase> (which has 0.0 mass)" warnings (PHREEQC prints at most 100). They are benign. Other warning types are noted in the table.

**Parents:**

| Generator | 25 °C | 100 °C | 160 °C |
|---|---|---|---|
| `tidal_kinetic` | **FAILED** in orbit 1, ocean MIX 4 step (35 s): "C has not converged. Total: 6.84e-21" | 0 errors, 660 s | 0 errors, 302 s; 41 "numerical method failed with this set of parameters" and 22 "maximum iterations exceeded" warnings, all recovered |
| `bulk_kinetic` (porewater, W/R 1) | 0 errors, 204 s | 0 errors, 290 s | 0 errors, 95 s |
| `bulk_equil` (porewater, W/R 1) | 0 errors, <1 s | 0 errors, <1 s | 0 errors, <1 s |
| `tidal_equil` | No temperature argument: one run at the ocean temperature, 0 errors, 1 s (14 recovered numerical warnings) | — | — |

**Restarts** (from the parents above):

| Restart | 25 °C | 100 °C | 160 °C |
|---|---|---|---|
| `tidal_kinetic_restart` | Not generated: the parent never completed an orbit. The generator refuses, which is correct | 0 errors, 93 s | 0 errors, 119 s |
| `bulk_kinetic_restart` | 0 errors, 114 s | 0 errors, 35 s | 0 errors, 37 s |

**Isolating the 25 °C failure (measured):**
- The same input forced to coupled CH₄ fails identically.
- The same input with `organic_wt_percent = 0` completes (10 rows, 0 errors, 38 s).

**Inferred cause:** at 25 °C the IOM pools release only ~10⁻²⁰ mol C per orbit, and PHREEQC cannot converge the ocean-side equilibrium on that trace carbon total with carbonate phases in EP2. This is a robustness problem in how IOM behaves at low T, not a database problem.

### 2.3 Mass conservation

**Setup:** the `tidal_kinetic` 160 °C input, first-orbit porewater REACT step (kinetics + EP1, closed before any MIX).

The test compares two quantities for each element:
- the change in the whole system, `SYS(E)` final minus initial, where initial = the porewater plus the EP1 metal and lawrencite;
- the amount released by the kinetic reactants, Σ (m0 − KIN) × stoichiometry, using IOM `-formula`, schreibersite `Fe 3 P 1` and the mineral formulas.

Values are printed at 15 significant figures.

| Element | Change in system | Released by kinetics | Relative residual |
|---|---|---|---|
| C | 1.683315e-05 | 1.683315e-05 | 6e-13 |
| H | 1.683624e-05 | 1.683624e-05 | 1e-10 |
| O | 5.630187 | 5.630187 | 3e-16 |
| N | 1.633178e-10 | 1.633178e-10 | 2e-8 |
| S | 1.607162e-06 | 1.607162e-06 | 4e-10 |
| Fe | 1.940294e-03 | 1.940294e-03 | 9e-13 |
| P | 6.462289e-04 | 6.462289e-04 | 7e-15 |
| Mg / Si | 1.526460 / 1.694486 | same | ~1e-16 |

- **Fe** starts at 1.997823 mol metal + 0.0095625 mol lawrencite in the system. Its change equals the troilite + schreibersite release.
- **IOM release this orbit:** 1.68e-5 mol C.
- **Scope:** MIX steps conserve by construction in PHREEQC, and weren't separately checked.

### 2.4 Restart read-back (measured)

- **`tidal_kinetic_restart`:** all 38 kinetic `-m0` values (9 minerals + 29 IOM pools) equal the parent's REACT-row `k_*` / `IOM_*_mol` values exactly.
- **`bulk_kinetic_restart`:** the same, except one 8e-9 mol K-feldspar remainder that differs by 1.2e-5 relative.
- **Ni:** carries over exactly.
- **Promoted primary:** K-feldspar promoted in a tidal restart comes back as an EP1 phase with exactly the parent's 6.4525e-4 mol and leaves KINETICS.
- **Molar masses:** the PUNCH multipliers and the restart divisors are the same constants in every generator (schreibersite 198.509).

---

## 3. Generator consistency

All six were generated under the same conditions: 160 °C, coupled, 15 wt% IOM, with the profile file (porewater 261.2 atm, ocean 66.7 atm).

| Aspect | tidal_kinetic | tidal_equil | bulk_equil | bulk_kinetic | bulk_kin_restart | tidal_kin_restart | Intended? |
|---|---|---|---|---|---|---|---|
| Primaries | 8 kinetic + schreibersite; Fe/Ni/lawrencite EP1 | all 8 EP1; no schreibersite | as tidal_equil | as tidal_kinetic | as tidal_kinetic | as tidal_kinetic | Yes; schreibersite has no log K |
| Rate laws | identical bodies, all 9 kinetic phases, all 4 kinetic generators; SSA 0.1 or 5; SR guard on all but schreibersite | — | — | identical | identical | identical | Yes. SSA values are unsourced (RECONSTRUCTION_STATUS §8 item 5) |
| IOM | 29 pools | none | none | 29 | 29 | 29 | **Unclear**: equilibrium generators ignore organics; undocumented |
| CH₄ mode / database | `resolve_ch4_mode()` | always coupled | always coupled | `resolve_ch4_mode()` | same | same | **No** for <150 °C equilibrium runs |
| Porewater T | `porewater_temp` | **ocean T (1.45 °C), no argument** | `water_temp` | `porewater_temp` | same | same | **No / undocumented** |
| EP1 secondaries | 79 incl. Quartz (temperature-gated ≥150 °C) | 78 | **78, no Quartz at 160 °C** (no gating) | 79 | 79 | 79 | bulk_equil gating: **no** |
| EP2 (ocean) | 81: no Quartz, + 3 hydrates | 78 | — | — | — | 81 | Yes (stability-gated) |
| Gases | 11 gases at `0 0` in every EP block; no GAS_PHASE | same | same | same | carried | carried | **Open decision** |
| SOLUTION pressure | porewater **66.7** (ocean value), ocean 66.7 | 66.7 / 66.7 | **261.2** (porewater) | **66.7** | 66.7 | 66.7 | **Open decision** |
| SELECTED_OUTPUT duplicates | none | none | none | none | none | none | — |
| PUNCH numbers | 10–110 (+85), IOM 210–490 | 10–110 (+85) | 10–110 (+85) | 10–110 (+95), 210–490 | 10–120, 210–490 | 10–120, 210–490 | No collisions; headings = PUNCH count in all 9 blocks |

Also: every kinetic CSV has two columns headed `Fe`, one for dissolved-Fe totals and one for the Fe metal phase. Lucas's readers take the last match, so it works, but it's fragile.

---

## 4. Open decisions (PI), with evidence

1. **Porewater pressure, and which generators use it.**
   - Only `bulk_equil` porewater mode reacts at the profile's 261.2 atm.
   - `tidal_kinetic`, `bulk_kinetic` and both restarts react the porewater at the ocean's 66.7 atm.
2. **Gas handling** (fixed-pressure GAS_PHASE, bubble fate, ocean-side plume loss). Unchanged: 1 atm `0 0` gas phases everywhere.
   - Now unblocked technically: all gases in both databases are ideal, and a 70 atm GAS_PHASE with H₂(g)/Mtg(g) runs.
   - Measured on 10-21: at 160 °C, a 70 atm gas phase raised dissolved H₂ 70-fold (to 0.073 molal). The full 2-orbit file ran 6× slower (267 s versus 44 s).
3. **IOM at low temperature.** `tidal_kinetic` at 25 °C fails on trace IOM carbon (6.8e-21 mol) in the ocean step; it runs with organics off.
   - Options: a release floor below which pools are dropped; skipping IOM below a temperature; or accepting no 25 °C tidal runs.
   - This touches IOM handling, so I didn't change it.
4. **Equilibrium generators:**
   - should they carry organics, follow `resolve_ch4_mode()` below 150 °C, and take a porewater temperature (`tidal_equil`)?
   - should `bulk_equil` apply the Quartz temperature gating?
5. **Metal corrosion kinetics and the spin-up.** Fe metal is still a `dissolve_only` equilibrium phase, fully corroded within the first orbit. Not addressed.
6. **`rec_eqb_minerals`** now defaults to off. Whether to turn it on is a modelling choice; it isn't needed for convergence.
7. **Surface areas** (0.1 versus 5 m²/g) are unsourced and decide promotion outcomes and rates (STATUS §8 item 5).
8. **Database consolidation** into one source plus generated variants. Lucas's call.

## 5. Not done

- No modelling decision was implemented. That covers gas phases, corrosion kinetics, spin-up, IOM calibration and Option 3.
- No production runs were changed or re-run. Results produced before the tephroite, rock-table and `rec_eqb` changes differ slightly from what main now generates.
- `PROGRESS.md` still needs the review-session update before `docs/reviews/` is archived.
