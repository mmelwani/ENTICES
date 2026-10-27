# tests/

Run everything from the repo root.

| File | How to run | Write-up |
|---|---|---|
| `iom_selftest.R` | `Rscript tests/iom_selftest.R` | IOM module calibration checks. |
| `db_load_check.R` | `Rscript tests/db_load_check.R` | Every `databases/*.dat` loads in PHREEQC with 0 errors. |
| `rock_table_check.R` | `Rscript tests/rock_table_check.R` | `ROCK_MOL_PER_KG` vs `reconstruction/rock_mol_per_kg.csv`. |
| `carbon_rebin_check.R` | `Rscript tests/carbon_rebin_check.R` | `carbon_by_oxidation_state()` (`phreeqc_tools.R`) on the N1 leak test, `N1_leak_T100.pqi` and `output/N1_leak_T100.{out,csv}`. |
| `c4only_db_check.R` | `Rscript tests/c4only_db_check.R` | `databases/coreclath_C4only.dat` (header, C(+4)-only speciation) and `restrict_to_database()` (`phreeqc_tools.R`) on `N1_leak_T100.pqi`; `scripts/make_c4only.py` regenerates the database exactly (needs Python, else skipped). |
| `test_A4_formulas.pqi` | `phreeqc tests/test_A4_formulas.pqi out.txt databases/coreclath_CH4coupled.dat` | `docs/PHREEQC_TEST_RESULTS_2026-09-20.md` |
| `test_B2_Hdg_carbon.pqi` | same pattern, coupled coreclath | `docs/PHREEQC_TEST_RESULTS_2026-09-20.md` |
| `test_B3_real_assemblage.pqi` | same pattern, `databases/coreclath_CH4coupled.dat` | `docs/PHREEQC_TEST_B3_real_assemblage_2026-09-27.md` |
| `output/*.out` | | The outputs those write-ups cite. |

`test_A_Sg_hydrogen.pqi` and `test_B_Hdg_carbon.pqi` are the **original,
superseded** inputs, kept as records only. They contain the defects found on
2026-09-20: their `DATABASE` keyword isn't first in the file, and their control
case is never equilibrated. A4 and B2 replace them. **Don't run them.** With a
real database path, `test_A_Sg_hydrogen.pqi` was seen to hang PHREEQC.
