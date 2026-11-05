# tests/

Run everything from the repo root.

| File | How to run | Write-up |
|---|---|---|
| `iom_selftest.R` | `Rscript tests/iom_selftest.R` | IOM module calibration checks. |
| `db_load_check.R` | `Rscript tests/db_load_check.R` | Every `databases/*.dat` loads in PHREEQC with 0 errors. |
| `rock_table_check.R` | `Rscript tests/rock_table_check.R` | `ROCK_MOL_PER_KG` vs `reconstruction/rock_mol_per_kg.csv`. |
| `carbon_rebin_check.R` | `Rscript tests/carbon_rebin_check.R` | `carbon_by_oxidation_state()` (`phreeqc_tools.R`) on the N1 leak test, `N1_leak_T100.pqi` and `output/N1_leak_T100.{out,csv}`. |
| `c4only_db_check.R` | `Rscript tests/c4only_db_check.R` | `databases/coreclath_C4only.dat` (header, C(+4)-only speciation) and `restrict_to_database()` (`phreeqc_tools.R`) on `N1_leak_T100.pqi`; `scripts/make_c4only.py` regenerates the database exactly (needs Python, else skipped). |
| `pressure_check.R` | `Rscript tests/pressure_check.R` | Profile pressures in every generator's input, the header stamps, the missing-profile error and the restart pressure check (generates inputs only). |
| `initial_fluid_check.R` | `Rscript tests/initial_fluid_check.R` | Each generator's initial porewater brought from 25 C to 25, 160 and 300 C: converges, at the profile pressure, with negligible O(0) and H(0). |
| `bulk_density_check.R` | `Rscript tests/bulk_density_check.R` | `data/mineral_molar_volumes.csv` covers every production phase; `solid_density()` closes mass and volume on a hand-built case; the Scenario R rock matches the reconstruction's phase library. |
| `restart_closure_check.R` | `Rscript tests/restart_closure_check.R` | A 2-orbit run vs a 1-orbit run plus a 1-orbit restart (tidal and bulk, 160 C, 15 wt% IOM): every element's system total agrees within 1e-9 relative or the solver's mole-balance tolerance. Runs PHREEQC, about 6 minutes. |
| `gas_eos_check.R` | `Rscript tests/gas_eos_check.R` | Every gas a generator puts in a GAS_PHASE has Peng-Robinson constants in all three databases; H2 at 160 C, 261.2 atm, PR against ideal. |
| `gas_mass_balance_check.R` | `Rscript tests/gas_mass_balance_check.R` | Whole-system element totals per orbit for both `pore_gas_fate` options (160 C, 15 wt% IOM, <= 1e-9); the no-bubble identity (25 C, no metal). Runs PHREEQC, about 15 minutes. |
| `cometary_ocean_check.R` | `Rscript tests/cometary_ocean_check.R` | The cometary-ice initial ocean: C(+4)-only leak check, element closure, the expected numbers, the include_O2 and include_NH3 toggles. |
| `exchange_check.R`, `kozeny_carman_check.R`, `porosity_evolution_check.R` | `Rscript tests/<name>` | The exchange physics, Kozeny-Carman k and the porosity-evolution driver (minutes each). |
| `grouped_validation_check.R` | `Rscript tests/grouped_validation_check.R` | Grouped time steps against one orbit per step to orbit 1000 (ocean C(+4), pH, H2, Na within 1%). ~1.5 h with the phosphate phases (the 1000 one-orbit steps take ~70 min). |
| `pipeline_check.R` | `Rscript tests/pipeline_check.R` | `run_case.R` on two short cases (organic and blank), `summary.csv` columns and units, the RUNINFO and `.pqi` commit stamps, all six figures and an overlay. ~10 min. |
| `resume_check.R` | `Rscript tests/resume_check.R` | A short case killed mid-chunk and rerun: resumes after the last complete chunk, duplicates none, ends with the uninterrupted run's `summary.csv` (1e-4 relative or 1e-6 absolute). Then the staged runs: raising `t_end` extends the finished run (organic and blank in lockstep, earlier chunks unchanged), lowering it keeps the run whole, and a run ended at a steady state stays done. ~20 min. |
| `paths_check.R` | `Rscript tests/paths_check.R` | Paths handed to PHREEQC: relative input, output and workdir; a missing output directory stops at once; the version check before chunk 1; `run_case.R` with relative paths. ~30 min on a loaded machine. |
| `liveness_check.R` | `Rscript tests/liveness_check.R` | The run watchdog: gives up only when files are stale and the process is gone, on consecutive checks over a minimum span; an unanswered process check counts as alive. `variant_status()`: only runs that ended at `t_end` reopen. ~1 min. |
| `run_all.R` | `Rscript tests/run_all.R [--quick]` | **All of the above** in turn, each in its own process; logs in `output/run_all/`; ends `ALL n TESTS PASSED`. |
| `test_A4_formulas.pqi` | `phreeqc tests/test_A4_formulas.pqi out.txt databases/coreclath_CH4coupled.dat` | `docs/PHREEQC_TEST_RESULTS_2026-09-20.md` |
| `test_B2_Hdg_carbon.pqi` | same pattern, coupled coreclath | `docs/PHREEQC_TEST_RESULTS_2026-09-20.md` |
| `test_B3_real_assemblage.pqi` | same pattern, `databases/coreclath_CH4coupled.dat` | `docs/PHREEQC_TEST_B3_real_assemblage_2026-09-27.md` |
| `output/*.out` | | The outputs those write-ups cite. |

`test_A_Sg_hydrogen.pqi` and `test_B_Hdg_carbon.pqi` are the **original,
superseded** inputs, kept as records only. They contain the defects found on
2026-09-20: their `DATABASE` keyword isn't first in the file, and their control
case is never equilibrated. A4 and B2 replace them. **Don't run them.** With a
real database path, `test_A_Sg_hydrogen.pqi` was seen to hang PHREEQC.
