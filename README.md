# ENTICES

Enceladus Tidally-Induced Circulation and Exchange Simulator.

## 1. What ENTICES is

ENTICES models water–rock reaction in Enceladus' porous core. An R model of tidally driven
porewater circulation (a two-layer gravity field, tidal stress over each 1.375-day orbit, Darcy
flow through the porous core) sets how much porewater exchanges with the ocean per orbit. It
writes PHREEQC inputs that react a reconstructed primordial rock (kinetic silicates and sulfide,
Fe–Ni metal, schreibersite and insoluble organic matter (IOM)) with porewater at the core's
temperature and pressure, and exchange that porewater with the ocean every orbit. PHREEQC runs
the inputs, and ENTICES' R tools read the results.

The science question is whether IOM degradation in the core adds carbon to the ocean (the
carbonate Cassini saw in the plume) on a timescale that constrains how long hydrothermal
circulation has lasted. ENTICES therefore starts, by default, from pure water and attributes ocean
carbon to the organics by difference against an otherwise identical run without them
(`pair_organic_run()`).

Documentation: this README (requirements, how to run, worked examples, repository map, references),
[PRODUCTION_GUIDE.md](PRODUCTION_GUIDE.md) (production runs), [docs/MODEL.md](docs/MODEL.md) (the model's
physics and chemistry, parameters and outputs) and [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md) (test and
restart internals, known failures).

## 2. Requirements

- R (tested with 4.5.2) with the packages the code loads: `tidyverse`, `conflicted`, `dplyr`,
  `patchwork`, `glue`, `R6` (`ENTICES.R`). The rate plots in `mineral_dissolution_constants.r` and
  the worked examples also need `ggplot2`.
- PHREEQC 3.8.6. `phreeqc_tools.R` looks for it in the `PHREEQC_EXE` environment variable, then
  at `C:/Program Files/USGS/phreeqc-3.8.6-17100-x64/bin/Release/phreeqc.exe`, then on the `PATH`.
  Every generated `.pqi` is stamped `# Target PHREEQC version: 3.8.6`, and `run_phreeqc()` warns
  if the executable reports another version, because minor phases and trace species differ between
  versions. Earlier outputs made with PHREEQC 3.5.0 on a Mac and rerun on 3.8.6 in Windows and Linux
  moved dissolved Mn 130× and sulfide ~11×, while bulk chemistry was stable. Treat the warning as
  "not reproducible as stamped".
- Python 3, only for the database scripts (`scripts/make_c4only.py`, `strip_mineral_vm.py`,
  `gas_eos.py`, `export_mineral_volumes.py`).
- Databases (`databases/`, see `databases/README.md`), chosen by `carbon_mode` and stamped in
  each `.pqi` header:

  | `carbon_mode` | Database |
  |---|---|
  | `"c4_only"` (default) | `coreclath_C4only.dat` |
  | `"coupled"` | `coreclath_CH4coupled.dat` |
  | `"ch4_uncoupled"` | `coreclath_CH4uncoupled.dat` |
  | `"temperature"` | uncoupled below `ch4_cutoff_C` (150 °C), coupled above |

## 3. How to run (regenerating production)

Run everything from the repo root.

Production runs are described in `PRODUCTION_GUIDE.md`. A **case** is a small R file in `cases/`.
`cases/reference.R` lists every field with its units, its range and the section of this README that
justifies it, and `cases/sweep/` contains the one-at-a-time sweep, in which each case inherits the
reference. `Rscript run_case.R cases/<name>.R` validates the case and runs it with
`run_tidal_evolution()` (grouped chunks, [docs/MODEL.md §1.5](docs/MODEL.md#15-long-run-stepping)) together with its same-rock **blank**. The
blank has the same minerals without the IOM and follows the organic run's φ, k and time steps chunk by
chunk, so the two share the exchange and the output times. The run writes `runs/<name>/organic/`,
`blank/`, `summary.csv` (one row per chunk end), `summary_units.csv` and `RUNINFO.txt` (commit, PHREEQC
version, database, machine, wall time, and what ended the run).

A run ends at the first of these end conditions:
- a steady state (`ss_tol` 0.01 on the ocean's C(+4), pH and Na over one circulation time);
- the backstop `t_end_yr` = 1e8 yr;
- clogging (φ ≤ `phi_min`, with the clog time as a `summary.csv` column);
- φ reaching φc.

Runs are resumable. A checkpoint (`state.rds`) is written after every chunk, and the same command
continues from the last complete chunk without repeating one (`tests/resume_check.R` kills a run
mid-chunk to test this). `Rscript run_sweep.R [--parallel N]` runs and resumes every case in
`cases/sweep/`. `Rscript make_figures.R` draws the six paper figures from any set of runs, with no
observational data, and its `overlay` argument takes the user's bands and points. `Rscript
tests/run_all.R` runs the whole test suite. Every generated `.pqi` header and RUNINFO records the
code's commit (`# Code: ENTICES commit <sha>`, with "(dirty: uncommitted changes)" if tracked files
were modified).

1. Pull `main`.
2. Set paths. PHREEQC as above (`Sys.setenv(PHREEQC_EXE = "...")` if needed).
3. Run the self-tests first. Each prints `SUMMARY: 0 failure(s)` when passing (`tests/README.md`
   describes each, and [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md#comparing-two-runs-column-by-column) gives the
   rule for comparing two runs):

   ```
   Rscript tests/iom_selftest.R           # IOM module calibration checks
   Rscript tests/rock_table_check.R       # ROCK_MOL_PER_KG vs reconstruction/rock_mol_per_kg.csv
   Rscript tests/db_load_check.R          # every databases/*.dat loads in PHREEQC with 0 errors
   Rscript tests/carbon_rebin_check.R     # carbon by actual oxidation state, on the N1 leak test output
   Rscript tests/c4only_db_check.R        # coreclath_C4only.dat: C(+4)-only carbon; restrict_to_database()
   Rscript tests/pressure_check.R         # profile pressures in every generator; restart pressure check
   Rscript tests/initial_fluid_check.R    # initial porewater from 25 C to 25/160/300 C: converges, no O2
   Rscript tests/bulk_density_check.R     # molar volumes, bulk_density(), t = 0 vs the reconstruction
   Rscript tests/gas_eos_check.R          # Peng-Robinson constants; H2 PR vs ideal
   Rscript tests/cometary_ocean_check.R   # the cometary-ice initial ocean
   Rscript tests/restart_closure_check.R  # restarts conserve every element (runs PHREEQC, ~6 min)
   Rscript tests/gas_mass_balance_check.R # gas phases conserve mass, both pore_gas_fate options (~15 min)
   ```

4. Choose the profile. `ph$profile_file` (default
   `profiles/EnceladusProfile_Seawater_10.0ppt_Tb272.4578K_PorousRock.txt`, a PlanetProfile model,
   `profiles/README.md`) sets the porewater pressure (layer mean, 261.2 atm), the ocean pressure
   (ocean floor, 66.7 atm) and, with `porewater_temp = "profile"`, the porewater temperature. There
   is no default pressure, so a missing profile is an error.
5. Set parameters ([docs/MODEL.md §2](docs/MODEL.md#2-parameters)): the circulation (`k`, `layer_thickness`), the organics
   (`calculate_phreeqc_parameters(organic_wt_percent = 15)`), the temperature and the run length.
6. Generate the `.pqi`:

   ```r
   source("ENTICES.R")
   sim <- Simulator$new(k = 1e-7, layer_thickness = 194200, config = EnceladusConfig$new())
   sim$run_simulation()
   ph <- PhreeqcIntegrator$new(sim)
   ph$calculate_phreeqc_parameters(organic_wt_percent = 15)
   f <- ph$generate_phreeqc_input_tidal_kinetic(output_dir = "phreeqc_inputs", mode = "ungrouped",
                                                n_orbital_cycles = 12, porewater_temp = 160, substeps = 1,
                                                trace_C_molal = 1e-14)
   ```

7. Run with the database stamped in its header (`# Database:`):

   ```r
   db <- sub("^# Database: ", "", grep("^# Database: ", readLines(f, n = 20), value = TRUE))
   r  <- run_phreeqc(f, db, output = sub("\\.pqi$", ".out", f), workdir = dirname(f))
   ```

   or `phreeqc <file>.pqi <file>.out <repo>/databases/coreclath_C4only.dat` from the input's
   directory.
8. Outputs (beside the `.pqi`): `<base>.csv`, the SELECTED_OUTPUT (tab-separated), `<base>.out`,
   and `<base>_initial_ocean.csv` for a cometary ocean. The CSV has these columns:
   - `state`, `soln`, `step`: `i_soln` rows are the initial solutions, followed by four `react` rows per
     orbit of a tidal run. These are REACT (soln 1: porewater + rock), MIX 4 (soln 4: ocean + exchanged
     porewater + seafloor phases + ocean gas), MIX 1 (soln 1: resampled porewater) and MIX 2 (soln 2:
     the ocean at the end of the orbit);
   - `pH`, `pe`, `Alk`, `mu`, `mass_H2O` (kg);
   - element totals (molal) for all 36 database elements, and `H_total_mol`, `O_total_mol` (moles,
     water included, whereas the `H` and `O` totals leave water out);
   - `m_` / `la_` species molalities and log activities;
   - equilibrium-phase moles by name and `d_` their change, `k_` kinetic reactants (moles) and `dk_`
     their change, `<mineral>_remain_g`, `IOM_*_mol`, and `Fe_metal` (the `Fe` column is dissolved Fe,
     and `Ni` appears twice, dissolved then the metal);
   - `si_` primary-mineral SIs. Use these, not the `.out` tables ([docs/DEVELOPMENT.md](docs/DEVELOPMENT.md#known-failures));
   - gas: `pressure`, `total mol`, `volume`, `g_<gas>` moles of the gas phase in use, and `Vgas_L`,
     `P_gas_atm`, `Vgas_Vpore` (gas volume / pore volume);
   - `Rho_soln` (PHREEQC's solution density).
9. Paired runs give the organics' contribution by difference against a same-rock blank:

   ```r
   source("organic_contribution.R")
   pr <- pair_organic_run(ph, "tidal_kinetic", mode = "ungrouped", n_orbital_cycles = 12, porewater_temp = 160,
                          trace_C_molal = 1e-14, output_dir = "phreeqc_inputs/pair160")
   pr$contribution  # DIC, CH4, NH4+, H2, total C (mol/kgw), ocean gas (mol): run, blank, difference
   ```

10. Sweeps: `parameter_sweep.R` generates inputs over a grid of `k` and layer thickness.
11. Restarts continue a run from its CSV (`prev_csv`), with the parent's `.pqi` beside it:
    `generate_phreeqc_input_tidal_kinetic_restart()`, `generate_phreeqc_input_bulk_kinetic_restart()`.
    They pass on every element at 15 digits, the equilibrium phases, kinetic reactants, IOM channels
    and both gas phases. They refuse a parent run in another carbon mode or at other pressures unless
    told (`carbon_mode`, `allow_pressure_change`). Their internals: [docs/DEVELOPMENT.md, Restart
    internals](docs/DEVELOPMENT.md#restart-internals).

12. Runtimes (this machine, PHREEQC 3.8.6): the 160 °C example case (§4, 12 orbits, 15 wt% IOM)
    takes 320–515 s at 150–160 °C, 2 orbits at 160 °C take 140–290 s, and the equilibrium generators
    take seconds. The first orbit (the metal pulse into the gas phase) is the slowest.
13. Common failures: [docs/DEVELOPMENT.md, Known failures](docs/DEVELOPMENT.md#known-failures).

## 4. Worked examples

`ENTICES_runs.Rmd` has each example as code, with narrative and expected results. **The 160 °C
example case** (example 2: k = 1e-7 m², full-core layer, 12 orbits, 15 wt% IOM) is the case the rest
of the documentation quotes. All expected numbers below include the phosphate phases.

1. A 25 °C smoke test (minutes): `tidal_kinetic`, 2 orbits, 15 wt% IOM, `trace_C_molal = 1e-14`.
   At orbit 2 (`to_ocean`), expect a porewater pH of 10.27, Na 0.115 molal, ionic strength 0.20 and
   4.3e-4 mol hydroxyapatite (pH 11.04 and ionic strength 0.26 without phosphate phases). Expect an
   ocean pH of 7.41, H₂ 0.065 molal (saturated), DIC 9.2e-15 (the seed) and 0.90 L of ocean gas, and
   a pore bubble of 0.44 L (2.8 pore volumes) after orbit 1, with none in orbit 2.
2. The 160 °C example case: k = 1e-7 m², full-core layer (69 km reacting skin), 12 orbits, 160 °C,
   15 wt% IOM, with the layered strain at K_d (exchange 4.5e-8 per half-orbit in the 69 km skin). At
   orbit 12, under the default `pore_gas_fate = "to_ocean"`, expect an ocean pH of 8.52, ocean H₂
   0.065 molal (saturated), ocean DIC 2.9e-10 molal (the seed is 1e-14) and 0.84 L of ocean gas. The
   pore bubble is 0.56 L (3.6 pore volumes) after the first orbit and 5e-4 L by orbit 12. ρ_solid
   goes from 2.92 to 2.76 g/cm³ and ρ_bulk from 2.31 to 2.24, and the pores have 0.0024 mol
   hydroxyapatite. Without the phosphate phases, the strain sensitivities were these: at K = 100 GPa
   the ocean reached pH 8.15 and DIC 1.5e-10, and under `strain_model = "liao"` pH 10.63 and DIC
   3.4e-8.
3. A 100 °C organic vs blank pair with `pair_organic_run()` (2 orbits). Expect the organics to add
   1.25e-15 molal DIC to the ocean (the blank's is the seed, 9.6e-15), 9.1e-15 molal CH₄ and 2.8e-16
   molal NH₄⁺. The ocean gas has 2.16 mol H₂ in both, from the metal. **An organic contribution of
   ~1e-15 molal is what this exchange rate (~5e-8 of the pore volume per half-orbit) and these skin
   temperatures give over a few orbits. It is not a measure of what the organics can supply.** The
   long runs (`docs/PROGRESS.md`) answer the organics question, and this example does not.
4. The cometary-ice ocean: `initial_ocean = cometary_ice_ocean()` (NH₃ included by default). The
   ocean starts at pH 2.08 with 0.087 molal DIC and 7.4e-3 molal NH₄⁺.
5. The pore-gas bracket: the 160 °C example case with `pore_gas_fate = "retain"` and `"to_ocean"`,
   each with its blank (layered strain at K_d). At orbit 12, ocean H₂ is 1.8e-8 (retain) against
   6.5e-2 molal (to_ocean, saturated, with 2.4 mol H₂ in the ocean gas). The organic DIC is 2.5e-10
   against 2.9e-10 molal, because under to_ocean the pore bubble's CO₂ (1.6e-11 mol) also goes to the
   ocean. Organic methane (Mtg) is 1.3e-16 against 1.9e-9 molal, plus 2.6e-8 mol in the ocean gas
   under to_ocean. The ocean pH is 8.52 under both.

## 5. Repository map

| Path | What it is |
|---|---|
| `ENTICES.R` | Entry point. Circulation model and all PHREEQC input generators (`generate_phreeqc_input_*`). |
| `ENTICES_runs.Rmd` | The worked examples (§4). |
| `iom_module.R` | IOM decomposition channels: calibrated sub-pools. |
| `mineral_dissolution_constants.r` | Mineral rate laws in R (used by `rec_eqb_minerals`) and rate plots. |
| `phreeqc_tools.R` | Runs PHREEQC, carbon re-binning, database name tools, `bulk_density()`, `system_totals()`. |
| `organic_contribution.R` | `pair_organic_run()`, `organic_contribution()`. |
| `parameter_sweep.R` | Parameter sweep helper. |
| `CHANGELOG.md` | Change log. |
| `PRODUCTION_GUIDE.md` | For production runs: setup, running cases and the sweep, figures, outputs, failures. |
| `run_case.R`, `run_sweep.R`, `make_figures.R` | Production: one case with its blank, the sweep, the figures (§3). |
| `production_kit.R` | Case files, validation, the run/resume machinery, `summary.csv` (used by the scripts above). |
| `porosity_evolution.R` | `run_tidal_evolution()`: chunked runs with evolving porosity, grouped steps, checkpoints. |
| `cases/` | `reference.R` and `sweep/`: the production cases. |
| `runs/`, `figures/` | Outputs (not in git). |
| `databases/` | PHREEQC databases (`databases/README.md`). |
| `profiles/` | The PlanetProfile interior profile (`profiles/README.md`). |
| `data/` | `mineral_molar_volumes.csv` (`data/README.md`). |
| `reconstruction/` | The primordial CI-chondrite rock reconstruction. |
| `tests/` | Self-tests and PHREEQC test inputs (`tests/README.md`). |
| `examples/` | An example production run and the run notebook by L. Fifer (`ENTICES_runs_lucas_2026-10.Rmd`), and the redox budget diagnostic inputs. |
| `scripts/` | Database scripts and standalone analyses. |
| `docs/` | The model (`docs/MODEL.md`), development notes (`docs/DEVELOPMENT.md`), design decisions, test write-ups, progress (`docs/PROGRESS.md`), `reviews/`, `archive/`. |
| `legacy/` | Retired files (`legacy/README.md`). **Do not use for runs.** |

## 6. References

- Alexander, C.M.O'D., Fogel, M., Yabuta, H., Cody, G.D. (2007). The origin and evolution of
  chondrites recorded in the elemental and isotopic compositions of their macromolecular organic
  matter. *Geochim. Cosmochim. Acta* 71, 4380–4403. doi:10.1016/j.gca.2007.06.052
- Alexander, C.M.O'D., et al. (2022). Sulfur abundances and isotopic compositions in bulk
  carbonaceous chondrites and insoluble organic matter. *Meteorit. Planet. Sci.* 57, 334–351.
  doi:10.1111/maps.13746
- Bagheri, A., Simons, M., Park, R.S., Berne, A., Hemingway, D., Melwani Daswani, M., Vance, S.D.
  (2025). Exploring the interior structure and mode of tidal heating in Enceladus. *Planet. Sci. J.*
  6, 245. doi:10.3847/PSJ/ae0cab
- Blanc, Ph., Lassin, A., Piantone, P., Azaroual, M., Jacquemet, N., Fabbri, A., Gaucher, E.C. (2012).
  Thermoddem: A geochemical database focused on low temperature water/rock interactions and waste
  materials. *Appl. Geochem.* 27, 2107–2116. doi:10.1016/j.apgeochem.2012.06.002
- Burnham, A.K. (2019). Kinetic models of vitrinite, kerogen, and bitumen reflectance. *Org.
  Geochem.* 131, 50–59. doi:10.1016/j.orggeochem.2019.03.007
- Casey, W.H., Hochella, M.F., Jr., Westrich, H.R. (1993). The surface chemistry of manganiferous
  silicate minerals as inferred from experiments on tephroite (Mn₂SiO₄). *Geochim. Cosmochim.
  Acta* 57, 785–793. (Not Casey et al. 1993, *Nature* 366, on chain-silicate leaching.)
- Delarue, C., Reynard, B., Sotin, C., Fellah, C., Cardon, H., Montagnac, G., Confortini, G.,
  Ferreiro Mählmann, R. (2026). Metamorphic and density evolution of carbonaceous matter: a
  potential cometary composition of Titan. *Earth Planet. Sci. Lett.* 691, 120227 (Eq. 10: IOM
  density from H/C).
- Choblet, G., Tobie, G., Sotin, C., Běhounková, M., Čadek, O., Postberg, F., Souček, O. (2017).
  Powering prolonged hydrothermal activity inside Enceladus. *Nature Astron.* 1, 841–847.
  doi:10.1038/s41550-017-0289-8
- Fisher, A.T., Dickerson, K.L., Blackman, D.K., Randolph-Flagg, N.G., German, C.R., Sotin, C.
  (2024). Sustaining hydrothermal circulation with gravity relevant to ocean worlds. *J. Geophys.
  Res. Planets* 129, e2023JE008202. doi:10.1029/2023JE008202
- Fray, N., Schmitt, B. (2009). Sublimation of ices of astrophysical interest: a bibliographic
  review. *Planet. Space Sci.* 57, 2053–2080.
- Glein, C.R., Waite, J.H. (2020). The carbonate geochemistry of Enceladus' ocean. *Geophys. Res.
  Lett.* 47, e2019GL085885. doi:10.1029/2019GL085885
- Green, E.C.R., Holland, T.J.B., Powell, R., Weller, O.M., Riel, N. (2025). *J. Petrol.*,
  doi:10.1093/petrology/egae079 (its `tc-ds636` dataset cross-checks the molar volumes).
- Iess, L., Stevenson, D.J., Parisi, M., Hemingway, D., Jacobson, R.A., Lunine, J.I., Nimmo, F.,
  Armstrong, J.W., Asmar, S.W., Ducci, M., Tortora, P. (2014). The gravity field and interior
  structure of Enceladus. *Science* 344(6179), 78–80. doi:10.1126/science.1250551
- King, A.J., et al. (2015). Modal mineralogy of CI and CI-like carbonaceous chondrites by X-ray
  diffraction. *Geochim. Cosmochim. Acta* 165, 148–160.
- Liao, Y., Nimmo, F., Neufeld, J.A. (2020). Heat production and tidally driven fluid flow in the
  permeable core of Enceladus. *J. Geophys. Res. Planets* 125, e2019JE006209.
  doi:10.1029/2019JE006209
- Lodders, K., et al. (2025). Solar system elemental abundances from the solar photosphere and
  CI-chondrites. *Space Sci. Rev.* 221, 23. doi:10.1007/s11214-025-01146-w
- Miller, K.E., Foustoukos, D.I., Cody, G.D., Alexander, C.M.O'D. (2025). Experimental heating of
  complex organic matter at Titan's interior conditions supports contributions to atmospheric N₂
  and CH₄. *Geochim. Cosmochim. Acta* 390, 38–56.
- Miller, K.E., Melwani Daswani, M., Sotin, C., et al. (2026). Titan's refractory core evolution:
  Implications for organics in its subsurface ocean. *Icarus* 449, 116961.
  doi:10.1016/j.icarus.2026.116961
- Palandri, J.L., Kharaka, Y.K. (2004). A compilation of rate parameters of water–mineral
  interaction kinetics for application to geochemical modeling. USGS Open-File Report 2004-1068.
- Parkhurst, D.L., Appelo, C.A.J. (2013). Description of input and examples for PHREEQC version 3.
  USGS Techniques and Methods 6-A43. https://pubs.usgs.gov/tm/06/a43/
- Pasek, M.A., Lauretta, D.S. (2005). Aqueous corrosion of phosphide minerals from iron meteorites:
  a highly reactive source of prebiotic phosphorus on the surface of the early Earth.
  *Astrobiology* 5, 515–535.
- Peng, D.-Y., Robinson, D.B. (1976). A new two-constant equation of state. *Ind. Eng. Chem.
  Fundam.* 15, 59–64.
- Pincus, L.N., Peters, C.A., Yi, C., Kermani, H.M. (2026). Manganese and cobalt in olivine-group
  silicates offer potential for natural H₂ generation. *ACS Earth Space Chem.* 10, 1174–1184.
  doi:10.1021/acsearthspacechem.5c00301 (Table S1: tephroite log k(H⁺) = −5.4, from Casey et al.
  1993)
- Rubin, M., et al. (2019). Elemental and molecular abundances in comet 67P/Churyumov-Gerasimenko.
  *Mon. Not. R. Astron. Soc.* 489, 594–607. doi:10.1093/mnras/stz2086
- Sloan, E.D., Koh, C.A. (2007). *Clathrate Hydrates of Natural Gases*, 3rd ed. CRC Press. Table
  2.5a, p. 78.
- Styczinski, M.J., Vance, S.D., Melwani Daswani, M. (2023). PlanetProfile: Self-consistent interior
  structure modeling for ocean worlds and rocky dwarf planets in Python. *Earth Space Sci.* 10,
  e2022EA002748. doi:10.1029/2022EA002748
- White, A.F., Brantley, S.L. (2003). The effect of time on the weathering of silicate minerals:
  why do weathering rates differ in the laboratory and field? *Chem. Geol.* 202, 479–506.
  doi:10.1016/j.chemgeo.2003.03.001
- Zandanel, A., Truche, L., Hellmann, R., Myagkiy, A., Choblet, G., Tobie, G. (2021). Short
  lifespans of serpentinization in the rocky core of Enceladus: implications for hydrogen
  production. *Icarus* 364, 114461. doi:10.1016/j.icarus.2021.114461
