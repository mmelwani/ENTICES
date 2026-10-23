# ENTICES

Enceladus Tidally-Induced Circulation and Exchange Simulator.

ENTICES models water–rock reaction in Enceladus' porous core. An R model of
tidally driven porewater circulation sets the flow, temperature and pressure;
it then writes PHREEQC input files that react a reconstructed primordial rock
(kinetic and equilibrium primaries, insoluble organic matter) with porewater,
exchanging fluid with the ocean each orbit. You run those inputs with PHREEQC.

## Layout

| Path | What it is |
|---|---|
| `ENTICES_LF_20260917.R` | **Entry point.** Circulation model and all PHREEQC input generators (`generate_phreeqc_input_*`). |
| `ENTICES_runs.Rmd` | Lucas's run notebook: example calls of the generators and restarts. |
| `iom_module.R` | Insoluble organic matter (IOM) decomposition module: calibrated sub-pools, sourced by the main file. |
| `mineral_dissolution_constants.r` | Mineral rate laws in R (used by `rec_eqb_minerals`) and rate plots. |
| `phreeqc_tools.R` | Locates and runs PHREEQC from R; warns if it isn't the expected version (3.8.6). |
| `parameter_sweep.R` | Parameter sweep helper used by `ENTICES_runs.Rmd`. |
| `CHANGELOG.md` | Change log for the main file. |
| `databases/` | PHREEQC thermodynamic databases. See `databases/README.md` for which to use when. |
| `reconstruction/` | Primordial CI chondrite mineralogy reconstruction: the source of the rock composition. |
| `tests/` | Self-tests and the PHREEQC test inputs, with their outputs in `tests/output/`. |
| `examples/` | Example production run (`from_lucas_20260924/`) and the redox budget diagnostic inputs. |
| `scripts/` | Standalone analysis scripts. |
| `docs/` | Design decisions, test write-ups and progress (`docs/PROGRESS.md`). `docs/reviews/` holds the review exchanges and `docs/archive/` superseded documents. |
| `legacy/` | `ENTICES_v3.R`: an old version kept for reference only. **Do not use** (see `legacy/README.md`). |

## Running

Run everything from the repo root.

1. Generate inputs from R, for example:

   ```r
   source("ENTICES_LF_20260917.R")
   enc <- EnceladusConfig$new()
   sim <- Simulator$new(k = 1e-7, layer_thickness = 194200, config = enc)
   sim$run_simulation()
   ph <- PhreeqcIntegrator$new(sim)
   ph$calculate_phreeqc_parameters(organic_wt_percent = 15)
   ph$generate_phreeqc_input_tidal_kinetic(n_orbital_cycles = 12, porewater_temp = 150,
                                           mode = "ungrouped", output_dir = "phreeqc_inputs")
   ```

   Generated files go to `phreeqc_inputs/`, which is git-ignored.

2. Run PHREEQC 3.8 on the generated file. Each file's header names the
   database it needs, `databases/coreclath_CH4coupled.dat` or
   `databases/coreclath_CH4uncoupled.dat`:

   ```
   phreeqc phreeqc_inputs/<file>.pqi phreeqc_inputs/<file>.pqi.out databases/coreclath_CH4coupled.dat
   ```

## Self-tests

```
Rscript tests/iom_selftest.R        # IOM module calibration checks
Rscript tests/rock_table_check.R    # ROCK_MOL_PER_KG vs reconstruction/rock_mol_per_kg.csv
Rscript tests/db_load_check.R       # every databases/*.dat loads in PHREEQC with 0 errors
```

Each prints `SUMMARY: 0 failure(s)` when passing. `db_load_check.R` needs PHREEQC; set
`PHREEQC_EXE` if it isn't at the default Windows install path or on the `PATH`
(see `phreeqc_tools.R`).

## PHREEQC manual

Parkhurst & Appelo (2013), *Description of input and examples for PHREEQC
version 3*, USGS Techniques and Methods 6-A43:
https://pubs.usgs.gov/tm/06/a43/
