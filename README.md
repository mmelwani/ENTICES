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
| `profiles/` | The PlanetProfile interior profile that sets temperatures and pressures (`profiles/README.md`). |
| `data/` | `mineral_molar_volumes.csv`: molar volumes for rock density (`data/README.md`). |
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
   database it needs, `databases/coreclath_CH4coupled.dat`,
   `databases/coreclath_CH4uncoupled.dat` or `databases/coreclath_C4only.dat`
   (chosen by `carbon_mode`; see `databases/README.md`):

   ```
   phreeqc phreeqc_inputs/<file>.pqi phreeqc_inputs/<file>.pqi.out databases/coreclath_CH4coupled.dat
   ```

## Carbon mode

`carbon_mode` picks how carbon redox is treated, and with it the database
(`databases/README.md`):

- `"c4_only"` (**default**, every temperature): coupled carbon is C(+4) only and
  reduced carbon enters as decoupled methane (`Mtg`) from the IOM.
- `"coupled"`: all carbon redox-coupled, so equilibrium reduces carbonate with
  H2. **Coupled runs are an upper bound on abiotic methanation**, not the
  expected case.
- `"ch4_uncoupled"`: methane decoupled (`Mtg`), other reduced carbon (ethane,
  acetate, formate, ...) coupled.
- `"temperature"`: the default before 2026-10-26, `"ch4_uncoupled"` below
  `ch4_cutoff_C` (150 °C) and `"coupled"` at or above it.

Restarts read the parent's mode from the database named in its `.pqi` header and
refuse to run in a different one unless `carbon_mode` (or `CH4_redox_override`)
is given, in which case they warn. The restart header records it as
`# parent carbon_mode:`.

## Initial fluid

The initial fluid is pure water by default, so all ocean carbon comes from the
organics. Two arguments change that:

- `trace_C_molal` (tidal_kinetic, bulk_kinetic): a C seed in every initial
  SOLUTION. Pure water with no C at all fails at 25 °C ("C has not converged")
  once the organics release their first trace of C. **Recommended production
  value: `trace_C_molal = 1e-14`**, the smallest value tested with both
  tidal_kinetic and its restart at 25 °C (2 orbits each, 0 errors). Keep it
  small: at 25 °C the organics release far less C than any seed, so ocean C is
  the seed.
- `initial_ocean` (tidal_kinetic, tidal_equil): `"pure_water"` (default) or a
  named list of element totals in mol/kgw for the ocean SOLUTION 2 only, e.g.
  `list(C = ..., N = ...)`.
- `charge_balance` (tidal_kinetic, tidal_equil): charge-balance the ocean
  SOLUTION 2 on `"pH"` or on one of the `initial_ocean` elements (PHREEQC
  adjusts its total). Default `NULL`: no charge balance, as before. A
  `trace_C_molal` seed enters as unbalanced HCO3-; `"pH"` balances it in the
  ocean (pH 4.8 at 1e-3), not in the porewater SOLUTION 1.

Both are stamped in the `.pqi` header, and restarts stamp their parent's values
as `# parent trace_C_molal:` / `# parent initial_ocean:`.

`pair_organic_run()` (`organic_contribution.R`) runs a model and an otherwise
identical blank with the organics absent (same mineral inventory), and returns
the organics' contribution to ocean DIC, CH4, NH4+ and H2 by difference.

## Interpreting outputs

### Rock density

PHREEQC reports no solid volumes, so rock density is computed in R from a run's selected output:

```r
source("phreeqc_tools.R")
d <- bulk_density("phreeqc_inputs/<file>.csv", porosity = sim$params$porosity)
d$summary     # per REACT step (and t = 0 from the .pqi): rho_solid, rho_fluid, rho_bulk
d$breakdown   # per phase: mass and volume, and their fractions
```

- `rho_solid` = Σ nᵢMᵢ / Σ nᵢVᵢ over every mineral with nonzero moles (equilibrium phases, kinetic
  primaries, schreibersite, Fe/Ni metal, clathrates) and the remaining IOM. Molar volumes come from
  `data/mineral_molar_volumes.csv` (sources in `data/README.md`), at 25 °C and 1 bar. Leaving out
  thermal expansion and compression costs about 1% for silicates at 160 °C and 26 MPa, less than the
  porosity and IOM-density uncertainty.
- Remaining IOM mass = initial IOM − Σ (m0ᵢ − nᵢ) × M(formulaᵢ) over the IOM channels (each
  channel's `-formula` is what the solid loses). Its density, `IOM_DENSITY_G_CM3` = 1.4 g/cm³, is
  **UNSOURCED** [citation needed: a measured density of chondritic IOM or kerogen].
- `rho_bulk` = (1 − φ) ρ_solid + φ ρ_fluid, with φ the model's porosity (`Simulator$params$porosity`,
  0.32) and ρ_fluid PHREEQC's solution density (`Rho_soln` column). Under the LLNL aqueous model
  PHREEQC does not compute water's density from T and P, so `Rho_soln` is about 1 g/cm³ at any
  temperature; give `fluid_density =` to override.

For Lucas's 160 °C production case (12 orbits, 15 wt% IOM, φ = 0.32, ρ_fluid = 1): ρ_solid
2.985 g/cm³ and ρ_bulk 2.349 g/cm³ at t = 0; 2.795 and 2.221 after 12 orbits, mostly from lizardite
(26% of the solid volume by orbit 12) taking up water. The profile itself gives a silicate-layer
density of 2353 kg/m³ (`rhoSil_kgm3` in its header).

**Comparing with Enceladus.** A published core density is a bulk density (porosity and pore fluid
included), so compare it with `rho_bulk`, not `rho_solid`. The usual reference is Iess et al. (2014)
[citation needed: the core density and its uncertainty, to be taken from the paper; the PDF was not
available here].

## Self-tests

```
Rscript tests/iom_selftest.R        # IOM module calibration checks
Rscript tests/rock_table_check.R    # ROCK_MOL_PER_KG vs reconstruction/rock_mol_per_kg.csv
Rscript tests/db_load_check.R       # every databases/*.dat loads in PHREEQC with 0 errors
Rscript tests/carbon_rebin_check.R  # carbon by actual oxidation state, on the N1 leak test output
Rscript tests/c4only_db_check.R     # coreclath_C4only.dat: C(+4)-only carbon; restrict_to_database()
Rscript tests/pressure_check.R      # profile pressures in every generator; restart pressure check
Rscript tests/initial_fluid_check.R # initial porewater from 25 C to 25/160/300 C: converges, no O2
Rscript tests/bulk_density_check.R  # molar volumes, bulk_density(), t = 0 vs the reconstruction
Rscript tests/restart_closure_check.R # restarts conserve every element (runs PHREEQC, ~6 min)
Rscript tests/gas_eos_check.R        # Peng-Robinson constants; H2 PR vs ideal
Rscript tests/gas_mass_balance_check.R # gas phases conserve mass, both pore_gas_fate options (~15 min)
```

Each prints `SUMMARY: 0 failure(s)` when passing. `db_load_check.R` and `c4only_db_check.R` need PHREEQC; set
`PHREEQC_EXE` if it isn't at the default Windows install path or on the `PATH`
(see `phreeqc_tools.R`).

## PHREEQC manual

Parkhurst & Appelo (2013), *Description of input and examples for PHREEQC
version 3*, USGS Techniques and Methods 6-A43:
https://pubs.usgs.gov/tm/06/a43/
