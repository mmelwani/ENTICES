# tests/

Run everything from the repo root.

| File | How to run | Write-up |
|---|---|---|
| `iom_selftest.R` | `Rscript tests/iom_selftest.R` | IOM module calibration checks. |
| `rock_table_check.R` | `Rscript tests/rock_table_check.R` | `ROCK_MOL_PER_KG` vs `reconstruction/rock_mol_per_kg.csv`. |
| `test_A4_formulas.pqi` | `phreeqc tests/test_A4_formulas.pqi out.txt databases/coreclath_CH4coupled.dat` | `docs/PHREEQC_TEST_RESULTS_2026-09-20.md` |
| `test_B2_Hdg_carbon.pqi` | same pattern, coupled coreclath | `docs/PHREEQC_TEST_RESULTS_2026-09-20.md` |
| `test_B3_real_assemblage.pqi` | same pattern, `databases/coreclath_CH4coupled.dat` | `docs/PHREEQC_TEST_B3_real_assemblage_2026-09-27.md` |
| `output/*.out` | | The outputs those write-ups cite. |

`test_A_Sg_hydrogen.pqi` and `test_B_Hdg_carbon.pqi` are the **original,
superseded** inputs, kept as records only. They contain the defects found on
2026-09-20: their `DATABASE` keyword isn't first in the file, and their control
case is never equilibrated. A4 and B2 replace them. **Don't run them.** With a
real database path, `test_A_Sg_hydrogen.pqi` was seen to hang PHREEQC.
