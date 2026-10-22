# Lucas's merged updates, the repo cleanup, and Lucas's rec_eqb timing pair: for review

Date: 2026-10-22
From: the VS Code (Claude Code) session.
For: the claude.ai review session. Lucas should see §1 and §3.

**Everything below has had only self-review.** Each section separates what was
**measured** (PHREEQC or R output, grep or diff of files) from what I
**inferred**.

Repo state: `main` @ `3ce883f` (PR #3 merged). The cleanup is on branch
`cleanup` (3 commits on top of `3ce883f`, pushed, PR to open).

---

## 1. Lucas's updates merged since `2d55fbb`

Two commits arrived via PR #3:
- `2e13e34`: `mineral_dissolution_constants.r`, `coreclath_CH4uncoupled.dat`;
- the merge `3ce883f`, which also brought in `ENTICES_runs.Rmd`, `parameter_sweep.R`, `CHANGELOG.md` and `ENTICES_LF_20260917.R`.

### 1.1 What we asked for

| Item | Status | Evidence |
|---|---|---|
| `ENTICES_runs.Rmd` no longer sources `ENTICES_v3.R` | **Done** | All 9 `source()` calls of the main file use `ENTICES_LF_20260917.R` |
| `parameter_sweep.R` missing | **Done** | Now committed |
| CHANGELOG: tephroite guard as the cause of the hang | **Done** | Entry dated 2026-10-19; the 10-06 entry's "K-Feldspar/Albite" blame is withdrawn |
| CHANGELOG: schreibersite entry | **Done** | Entry dated 2026-09-30, superseding the instant-products plan |
| `depletion_order` comment's missing CHANGELOG reference | **Resolved** | The comment's "CHANGELOG 2026-10-01" now matches an entry of that date. That entry describes the ordering but has no rate-ranking table |
| Database: decoupled gases made ideal | **Not done** | All five decoupled gases keep active `-T_c/-P_c/-Omega` in all three `.dat` files |
| Database: `Ntg(g)` `T_c` typo | **Partly done** | Fixed in `coreclath_CH4uncoupled.dat` only; still `T_c` (no dash) in `coreclath_CH4coupled.dat` and `Core11_idealgas_mod_v4.dat` |

### 1.2 What we flagged that remains (measured, by reading the code on `main`)

- **`rec_eqb_minerals` still evaluates every primary at pH 11** (`k_primary_mineral(m, H = 1e-11, ...)`). For albite and K-feldspar, pH 11 gives the *fastest* rate, so the ×2 margin isn't conservative for them.
- **The tephroite `kneut` is still summed with `kacid`**, as `k = (kacid + kneut)` with `kneut` fixed at pH 7. It is now in all four kinetic generators, because Lucas ported the rate law into both restarts.
- **`bulk_equil` porewater mode still uses ocean pressure**: SOLUTION 1 has `pressure {round(ocean_pressure_atm, 1)}`, while `pore_pressure_atm` is computed but unused there.
- **`tidal_kinetic_restart` lists Magnetite twice in SELECTED_OUTPUT.** Magnetite is in `secondary_minerals` and is still `append`ed after Lizardite in `.all_eq_so` (line ~4405).

### 1.3 New problems

**(a) `coreclath_CH4uncoupled.dat` no longer loads. This is the serious one.**
- **Measured:** a bare `SOLUTION 1 / END` run with this database gives 11 errors: "Elements in species have not been tabulated" and "Reaction for species has not been defined" for `Hdg`, `Mtg` and `H2Sg`.
- **Measured:** the same test on the `2d55fbb` version of the file gives 0 errors.
- **Cause, measured by diff:** `2e13e34` deleted the "5 Redox-uncoupled gases" SOLUTION_SPECIES block (`Hdg = Hdg`, `Oxg = Oxg`, `Mtg = Mtg`, `Ntg = Ntg`, `H2Sg = H2Sg`, with their `-dw/-Vm`). The commit message says it removed "repeated species definitions". In the previous file each identity reaction appeared **once**, so they weren't repeats. These are the master species' defining reactions.
- **Inferred:** every decoupled-CH₄ run fails at load. That is the default below 150 °C (`resolve_ch4_mode()`), or any run with `CH4_redox_override = "decoupled"`.
- **Fix:** restore that block.

**(b) `ENTICES_runs.Rmd` sources 15 helper scripts that aren't in the repo.** Examples are `minerals_plot.R`, `chemistry_plot.R`, `check1_regression.R` and `stitch_plots.R`. Those chunks will fail for anyone but Lucas. Line 698 also sources `mineral_dissolution_constants.R` (capital R), while the file is `.r`. That works on Windows and default macOS, but fails on case-sensitive file systems.

### 1.4 Runs from `main` (measured)

- **`iom_selftest.R`:** 0 failures.
- **`tidal_kinetic`** (160 °C, 15 wt% organics, 2 orbits, coupled): generated and ran in 42 s with 0 errors, 11 CSV rows.
- **`bulk_kinetic`** (porewater mode, W/R 1, 160 °C, 2 orbits): ran in 27 s with 0 errors. Its header shows `rec_eqb_minerals` promoted albite and K-feldspar.
- **Side effects of sourcing:** none. The plotting code in `mineral_dissolution_constants.r`, which `bulk_kinetic` sources, is wrapped in a function and doesn't run on `source()`.

---

## 2. Cleanup (branch `cleanup`)

### 2.1 Commits

| Commit | What |
|---|---|
| `78730a8` | The pending B3 write-up correction: pyrite claim withdrawn, calcite 8.055e-04 |
| `29a60ff` | Reorganisation (layout below), path fixes, `.gitignore` patterns, README files |
| `c31760a` | `reconstruction/` module, CSV hand-off and `tests/rock_table_check.R` |

### 2.2 Layout

```
ENTICES_LF_20260917.R  ENTICES_runs.Rmd  iom_module.R  mineral_dissolution_constants.r
parameter_sweep.R  CHANGELOG.md  README.md  LICENSE
databases/       3 .dat + README (which generator uses which)
reconstruction/  Rmd (2026-09-16), montecarlo Rmd, STATUS, rock_mol_per_kg.csv,
                 notes/ (tephroite), archive/v0_uncommitted_edit_2026-09.patch
tests/           iom_selftest.R, rock_table_check.R, test_*.pqi, output/*.out, README
examples/        from_lucas_20260924/, redox_budget/
docs/            current design docs; reviews/ (this exchange); archive/
legacy/          ENTICES_v3.R + README (dosing bug)
scripts/         Circulation_timescale_bounds.R, residue_ratio.R
```

### 2.3 Removed

All of these remain in git history:
- the duplicate `HANDOVER_IOM_validation_2026-09-17.txt`;
- the PHREEQC manual PDF, now linked from the README;
- the v0 `Primordial_CI_chondrite_mineralogy.Rmd`. The PI's uncommitted edit to it is saved as a patch, and I checked that the patch applies cleanly.

**Not imported** from `~/Documents/Enceladus_hydrotidal`:
- the xlsx (PI's decision);
- `.RData`, `.Rhistory`, `.Rproj*`, the `.nb.html`;
- `From Lucas 20260727/`;
- `Superseded/`;
- the duplicate Core11.

### 2.4 Deleted (PI-approved)

- the local `From Lucas 20260924/` folder (all four files verified byte-identical to blobs in history);
- branches `iom-module`, `lucas-integration` and `rock-composition`, locally and on GitHub, after confirming all three are merged into `main`.

### 2.5 References

Every moved filename was grepped.

**Fixed:**
- **The generator's recommended-database string** now reads `databases/...`. This is the only edit to `ENTICES_LF_20260917.R`.
- **Test `DATABASE` lines** are now relative to the repo root.
- **Path mentions** in `tests/iom_selftest.R`, `examples/redox_budget/compare.R` and the living docs in `docs/`.

**Left alone:**
- **The dated records in `docs/archive/` and `docs/reviews/`** keep the old flat paths; `docs/README.md` says so.
- **`test_A_Sg_hydrogen.pqi` and `test_B_Hdg_carbon.pqi`** (the superseded originals) keep their `/path/to/` placeholders.
  - **Measured:** with a real database path, `test_A_Sg_hydrogen.pqi` hung PHREEQC at 12 GB of memory.
  - `tests/README.md` marks both as records only, not to be run.

### 2.6 Post-move checks (measured)

| Check | Result |
|---|---|
| `Rscript tests/iom_selftest.R` | 0 failures |
| `test_A4`, `test_B2`, `test_B3` .pqi from the repo root | 0 errors each |
| `scripts/residue_ratio.R`, `examples/redox_budget/compare.R` | run, same outputs as before |
| Generated `tidal_kinetic` (160 °C, 2 orbits) | 47 s, 0 errors; header names `databases/coreclath_CH4coupled.dat` |
| `Rscript tests/rock_table_check.R` | **4 failures, left failing on purpose** (§2.7) |

### 2.7 The rock-table check fails for 4 phases (measured)

`reconstruction/rock_mol_per_kg.csv` comes from running the Rmd's own code. Its Scenario R fit reproduces the 12 wt% values, summing to 100%.

The check compares `ROCK_MOL_PER_KG` against the CSV. ENTICES's Fe metal is compared after adding its lawrencite back; magnesiochromite is inert and must be absent from the table. The tolerance is 1e-4 relative or 1e-6 absolute, as the PI set. Seven active phases plus the magnesiochromite and no-extra-entries checks pass. Four phases fail:

| Phase | ENTICES | CSV | rel |
|---|---|---|---|
| Tephroite | 0.025995 | 0.025987 | 3.0e-4 |
| K-feldspar | 0.020515 | 0.020521 | 2.9e-4 |
| Anorthite | 0.062148 | 0.062135 | 2.1e-4 |
| Schreibersite | 0.047101 | 0.047093 | 1.6e-4 |

**Inferred cause:** `ROCK_MOL_PER_KG` was built from STATUS §1's wt% values rounded to 3 decimals, not from the fit. For example, tephroite was entered as 0.525 against the fit's 0.52483, which alone is 3.2e-4 relative. The 6-decimal mol/kg rounding (≤ 2.4e-5) isn't the issue. Small phases lose the most to wt% rounding.

**Decision owed (PI/Lucas):** either regenerate `ROCK_MOL_PER_KG` from the CSV (changes ≤ 3e-4 relative), or loosen the tolerance to about 5e-4. I recommend regenerating, since that is the point of the checked hand-off.

---

## 3. Lucas's runs in `From Lucas 20261022`

The folder holds two `.pqi`/`.out` pairs and no CSVs. They are a `rec_eqb_minerals` TRUE/FALSE timing pair. Neither is committed.

### 3.1 What they are (measured, from headers, diff and outputs)

**Generator:** `generate_phreeqc_input_bulk_kinetic`.

**Parameters:**

| Parameter | Value |
|---|---|
| `mode` | `"porewater"` (no ocean) |
| `water_rock_ratio` | 0.1574 |
| Rock | 1 kg |
| Temperature | 160 °C |
| `pressure` | 66.7 atm |
| Organics | 0 wt% (`org0`) |
| Orbits | 3 |
| Steps | `-steps 3.564e5 in 3 steps` (one per orbit) |
| Integrator | CVODE |

- **Database:** `coreclath_CH4coupled.dat`, run on **PHREEQC 3.5.0** (Lucas's Mac); ours is 3.8.6.
- **The only difference between the two inputs:** albite and K-feldspar.
  - TRUE: `dissolve_only` EQUILIBRIUM_PHASES, with `EQUI()` punches.
  - FALSE: KINETICS with RATES blocks.
- **Errors:** 0 in both.
- **Warnings:** 100 in each, all of the "Element X is contained in <phase> (which has 0.0 mass)" kind (PHREEQC prints at most 100). They're benign.
- **Wall time on Lucas's machine:** TRUE 161 s, FALSE 279 s (×1.7).

### 3.2 Known issues

**Measured:**
- **Gases at `0 0`:** yes, all 11 gas phases, including `H2(g) 0 0`. H₂(g) holds 3.08 mol after orbit 1 and 3.15 mol after orbit 3, in both runs. This mode has no ocean, so "locked out of the ocean" doesn't apply. The 1 atm infinite reservoir does: dissolved H₂ stays at about 8–9×10⁻⁴ molal throughout, instead of about 0.07 molal at 70 atm (§2.5 of the 10-21 response).
- **Pressure:** 66.7 atm, the hand-entered value, not 261.
- **Fe metal:** fully gone by the end of orbit 1 in both runs (2.350 mol, no Fe left in the assemblage). Magnetite is 0.72 mol.
- **Water loss:** porewater falls from 0.1574 kg to 0.072 kg (−54%) by the end of orbit 1, and to 0.059 kg by orbit 3.

**Inferred:**
- **The metal corrodes at t = 0.** These outputs only print end-of-orbit states, so the timing isn't measured here. It is inferred from the 10-21 time-zero test, where the same EP1 corroded all the metal at t = 0.
- **The magnetite in EP1 and the tephroite guard are both present (measured from the inputs), so this run doesn't carry the 10-06 hang.**

### 3.3 Does `rec_eqb_minerals` change the chemistry?

**Measured:** end-of-orbit values from Lucas's outputs, TRUE versus FALSE.

| Quantity (end of orbit 3) | TRUE | FALSE | rel diff |
|---|---|---|---|
| pH | 10.94 | 10.97 | 0.03 units |
| pe | −10.88 | −10.91 | |
| H₂O (kg) | 0.05886 | 0.05874 | 0.2% |
| H₂ (aq, molal) | 8.32e-4 | 8.27e-4 | 0.6% |
| H₂(g) (mol) | 3.152 | 3.152 | 0 |
| Magnetite (mol) | 0.7895 | 0.7895 | 0 |
| Lizardite (mol) | 0.9085 | 0.9253 | 1.8% |
| Analcime (mol) | 0.2870 | 0.2846 | 0.8% |
| Na (molal) | 0.725 | 0.765 | 5.2% |
| K (molal) | 0.3485 | 0.3493 | 0.2% |
| Si (molal) | 0.0239 | 0.0251 | 4.8% |
| Al (molal) | 0.0304 | 0.0321 | 5.3% |
| Clinochlore-14A (mol) | 0.0458 | 0.0387 | 16% |
| Clinozoisite (mol) | 0.0302 | 0.0356 | 15% |
| Diopside remaining (mol) | 0.2745 | 0.2636 | 4.0% |
| Forsterite remaining (mol) | 0.1474 | 0.1456 | 1.2% |
| S(−2) (molal) | 2.4e-6 | 2.6e-5 | 91% (11×) |
| Mn (molal) | 2.5e-4 | 9.1e-4 | 73% (3.7×) |

These runs have no carbon (`org0`), so CH₄ and DIC don't exist; I compared H₂, SiO₂, Na and K instead.

**Largest relative differences:**
- **Traces:** dissolved sulfide (11×) and Mn (3.7×).
- **Minor secondaries:** clinochlore and clinozoisite (about 15%). In orbit 1, tobermorite (0.0013 vs 0.0033) and chamosite also differ.
- **Major-element chemistry:** within about 5%. That covers pH, pe, H₂, the gas, magnetite, lizardite, Na, K, Cl and the primaries.
- **K-feldspar:** 0.01155 mol (TRUE, as an equilibrium phase) versus 0.01125 (FALSE, kinetic) after orbit 1, about 3%. Both are gone by orbit 3.

**Measured: the same inputs re-run on PHREEQC 3.8.6 with our copy of the database.**
- **Timing reverses:** TRUE 75 s, FALSE 42 s.
- **The same minor quantities move as much between PHREEQC versions as between TRUE and FALSE:**

| Quantity | Lucas 3.5.0 TRUE | Ours 3.8.6 TRUE |
|---|---|---|
| Mn (end of orbit 3) | 2.5e-4 | 3.2e-2 (130×) |
| Clinochlore-14A | 0.0458 | 0.0307 |
| Saponite-Fe-Na (orbit 1) | 0 | 0.0196 |
| Diopside remaining | 0.2745 | 0.2522 |

  Major quantities agree across versions: pH within 0.1, H₂(g) and magnetite identical by orbit 3, lizardite within 4%.

**Inferred:**
- **Diopside, tephroite and K-feldspar each dissolve partly in orbit 1 and then stop.** Each is guarded by `SR > 1 → rate 0`. So how much of each dissolves depends on the transient path inside the first macro step, which is resolved by CVODE inside one 118,800 s step. The trace and minor-phase differences (Mn from tephroite, Ca–Al phases fed by diopside and anorthite, sulfide) track that path.
  - **Measured:** diopside remaining is 0.2745 / 0.2636 / 0.2522 across three runs. Tephroite is 0.02599 in three runs and 0.02504 in our TRUE run, which matches its 130× Mn.
- **So `rec_eqb_minerals` doesn't change the bulk chemistry**, to within about 5%. Its effect on minor phases and traces isn't distinguishable from solver and version path-dependence, so neither run's minor-phase values should be quoted to better than a factor of a few at this step size.
- **The speed-up isn't robust either.** It is 1.7× on Lucas's 3.5.0 and reverses to 0.56× on our 3.8.6. Promotion now looks like a modelling choice, not a convergence fix. The actual hang was the tephroite guard, as both CHANGELOG and the 10-21 test agree.
- **Next test to settle the path-dependence:** finer `-steps` (say 30 per orbit) on both versions. If the minor phases converge as the step shrinks, they are step-size artefacts.

---

## 4. Decisions owed

1. **Restore the uncoupled database's SOLUTION_SPECIES block (Lucas).** This blocks every run below 150 °C.
2. **`ROCK_MOL_PER_KG`:** regenerate from the CSV, or loosen the test tolerance (PI/Lucas).
3. **Database PR edits** (all decoupled gases ideal, `Ntg` typo in the other two files): still owed. The 10-21 review asked for them.
4. **The pending items in §1.2.**
5. **Unchanged from 10-06/10-21:** porewater pressure, bubble fate, ocean-side gas, spin-up and the tephroite floor.
