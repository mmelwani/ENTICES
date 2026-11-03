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
`run_tidal_evolution()` (grouped chunks, §6 *Long runs*) together with its same-rock **blank**. The
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
   describes each):

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

   `column_differences()` (`phreeqc_tools.R`) compares two runs column by column. A value passes
   within 1e-6 relative or 1e-12 molal absolute. Log-scale columns (pe, pH, la_, si_) are compared only
   where the species behind them is above 1e-10 molal in both runs, and there they pass within 1e-4 log
   units. The species behind each column are these:
   - la_X: m_X;
   - pH: the larger of m_H⁺ and m_OH⁻;
   - si_X: the scarcest dissolved element of X other than H and O;
   - pe: any poising redox pair with both members above 1e-10 (H₂, Fe²⁺/Fe³⁺, SO₄²⁻/HS⁻ + H₂S,
     carbonate/Mtg or N₂/NH₄⁺ + NH₃).

   Below that gate the molality's absolute floor covers them, because the log of a trace or unpoised
   quantity has solver noise (la_SO₄²⁻ differed by 1.2e-4 log units at ≤ 4e-9 molal).

4. Choose the profile. `ph$profile_file` (default
   `profiles/EnceladusProfile_Seawater_10.0ppt_Tb272.4578K_PorousRock.txt`, a PlanetProfile model,
   `profiles/README.md`) sets the porewater pressure (layer mean, 261.2 atm), the ocean pressure
   (ocean floor, 66.7 atm) and, with `porewater_temp = "profile"`, the porewater temperature. There
   is no default pressure, so a missing profile is an error.
5. Set parameters (§5): the circulation (`k`, `layer_thickness`), the organics
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
   - `si_` primary-mineral SIs. Use these, not the `.out` tables (below);
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
    told (`carbon_mode`, `allow_pressure_change`).

    The solutions are passed on verbatim. Every `tidal_kinetic` input and its restart end by DUMPing
    solutions 1 and 2 as SOLUTION_RAW to `<base>_dump.txt` beside the CSV, and the restart reads them
    back verbatim, so total H, total O and the charge balance are the parent's. A new temperature
    enters only through REACTION_TEMPERATURE. For a parent without a dump, the restart warns and
    redefines the solutions from pH, pe and totals at the new temperature, which changes total H by
    4–6e-6 at every restart. In the 24-orbit skin-temperature evolution, the per-chunk drift in H is
    ≤ 2.1e-9 with the dump and 4.0e-6 without it, and a restart at 130 °C from a 160 °C parent keeps
    the whole-system totals to 5.4e-10 (`tests/restart_closure_check.R`).

    The whole-system element totals of a restart agree with the continuous run to ~1e-9 (closure),
    but its solution composition agrees only to ~1e-5 relative (porewater Na 6e-6 after a 1 + 1-orbit
    restart at 25 °C, under every `pore_gas_fate`). This is the kinetic integrator's accuracy, and no
    digit is lost. CVODE restarts its step history at a restart, and its own tolerance moves the
    continuous run by ~4e-4 in Na between `-tol` 1e-8 and 1e-12, whereas the restart difference does
    not converge with it (6e-6, 1e-6, 1.3e-5 at 1e-8, 1e-10, 1e-12). The solutions and the pore gas
    are passed on verbatim.
12. Runtimes (this machine, PHREEQC 3.8.6): the 160 °C example case (§4, 12 orbits, 15 wt% IOM)
    takes 320–515 s at 150–160 °C, 2 orbits at 160 °C take 140–290 s, and the equilibrium generators
    take seconds. The first orbit (the metal pulse into the gas phase) is the slowest.
13. Common failures:
    - *No convergence, or a step running for tens of minutes:* check that the input has the KNOBS
      `-convergence_tolerance 1e-8` after SELECTED_OUTPUT. At 25 °C, or with a cometary ocean, give
      `trace_C_molal = 1e-14`.
    - *The version warning:* rerun on PHREEQC 3.8.6.
    - *25 °C ocean step:* without extra scaling, the 25 °C ocean step (MIX 4) fails to converge
      ("Numerical method failed on all combinations", with Al, C, Ca, Cl, Fe, K, Mg, Mn, Na, Ni, P, S
      and Si "not converged" at ~1e-15–1e-8 mol) for scattered porewater fractions. Of fractions from
      1e-8 to 1e-3, 3e-8, 2e-6 and 5e-6 failed while their neighbours ran. An earlier version of the
      exchange code failed at 3.3e-6 too, and `retain` fails the same way, so the cause is not the gas
      transfer. Every ocean step therefore runs with KNOBS `-diagonal_scale true`, switched off after it
      (`OCEAN_STEP_KNOBS_ON/OFF`). All fractions tested converge, and a converging 160 °C orbit is
      unchanged (1e-8 relative). Set for the whole run instead, the option shifted the porewater step's
      m_H₂S by 0.9%. Trace seeds were not needed.
    - *A cometary ocean with a carbon-free porewater:* `tidal_kinetic` is ~4× slower (and stalls
      outright with a CO₂-only ocean), and `tidal_equil` fails ("Cl has not converged", the Cl being
      the porewater's lawrencite). Give `trace_C_molal = 1e-14` (both warn). One combination is a
      known failure: `tidal_equil` + `pore_gas_fate = "to_ocean"` + cometary ocean + `c4_only` crashes
      PHREEQC 3.8.6 (segmentation fault in the ocean step, at all temperatures), whereas
      `tidal_kinetic` runs it. A 20-line reproducer is in `tests/known_failures/`, and
      `docs/MATRIX_2026-10-28.md` has the details.
    - *SIs in the `.out` file:* the printed SI tables are offset by n(H₂O) × Δlog a_w (0.143 per H₂O
      at 25 °C in the brine). Use the CSV's `si_` columns.
    - *B-dot activities above ~1 molal:* the 25 °C NaOH porewater brine reaches I ≈ 2.4 for the
      first ~3×10⁴ orbits, so treat its pH and activities as qualitative.
    - *Disk space:* `KNOBS -logfile` is off by default. Every generator takes `logfile = TRUE` to
      switch it on for convergence debugging. `phreeqc.log` then records every solver iteration:
      258–346 MB per 12-orbit production run at 150–160 °C (~300 MB), and 11.9 GB from one stalled
      run.

## 4. Worked examples

`ENTICES_runs.Rmd` has each example as code, with narrative and expected results. **The 160 °C
example case** (example 2: k = 1e-7 m², full-core layer, 12 orbits, 15 wt% IOM) is the case the rest
of this README quotes. All expected numbers below include the phosphate phases.

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

## 5. Parameters

Status: **sourced** (a reference), **placeholder** (a stand-in with a stated basis), **UNSOURCED**, or
**decided** (a modelling choice). An empty cell means not applicable.

| Name | Default | Units | What it does | Justification / source | Status | Sensitivity planned |
|---|---|---|---|---|---|---|
| `k` | none (production 1e-7) | m² | permeability, sets Darcy flow and exchange | literature ranges (§6, Flow driving), production value at the top of Liao et al.'s range | decided | yes (sweep) |
| `layer_thickness` | none (production 194200) | m | permeable layer, 194.2 km is the whole core | production choice | | yes |
| `core_rigidity`, `core_viscosity` (`EnceladusConfig`) | 17e9, 1e12 | Pa, Pa·s | η selects the layered strain within `tidal_scenario`, µ enters the skin depth (and the `"liao"` strain) | Bagheri et al. 2025 Table 5, at the geometric means of 4–70 GPa and 4.4e11–2.3e12 | placeholder | yes (η over the scenario's range) |
| `strain_model`, `tidal_scenario` (`EnceladusConfig`) | `"layered"`, `"dissipative_core"` | | the core's strain and the ocean load from a layered forward model (or `"liao"`, an isolated core) | Bagheri et al. 2025 structure, TidalPy (§6, The strain) | decided | yes (both scenarios) |
| `core_bulk_modulus` (`EnceladusConfig`), `iom_wt_percent`, `iom_vol_frac` (`Simulator`) | NULL (= K_d), 15, NULL | Pa, wt%, none | the core's bulk modulus for the strain: the drained frame K_d = K_grain (1 − φ/φc), K_grain the Hill average of the minerals and the IOM | Abers & Hacker 2016, Bass 1995, Holland & Powell (Perple_X hp634ver.dat), Delarue et al. 2026 (§6, The strain) | decided (troilite μ an estimate, IOM μ a terrestrial-kerogen bracket) | optional (μ brackets) |
| `reactive_depth` (`Simulator`) | `"skin"` | | the reacting box: the top min(L, δ), or `"layer"` | | decided | |
| `strain_scale` (`EnceladusConfig`) | 1 | | factor on the core's tidal strain (not the ocean load) | an override for sweeps | | optional |
| `phi0` (`Simulator`, the older name `porosity` still works and wins if given) | 0.32 | | initial porosity: k, α, storage, δ, strain volume, pore volume and porewater mass | PlanetProfile `phiRockMax` (one fit, Styczinski et al. 2023) | **one fit**, with an end-state check against Bagheri's core density | yes (0.20, 0.38 in the sweep) |
| `phi_c` (`Simulator`) | 0.4 | | critical porosity: α = min(1, φ/φc), φ ≥ φc is a suspension (α = 1) | Nur et al. 1998; Dvorkin & Nur 2000 (36–40% sorted grains, sandstones) | placeholder | yes (0.4–0.6) |
| `alpha` (`Simulator`) | NULL (from φ/φc) | | Biot coefficient override (must be ≥ φ) | | | |
| `eccentricity` (`EnceladusConfig`) | 0.0047 | | tidal strain | Enceladus's orbit | | |
| `max_tidal_stress` (`EnceladusConfig`) | NULL: ΔP = α·ε_v·M | Pa | tidal pore-pressure amplitude driving the Darcy flow | derived from the strain, a number overrides it | decided (derived, override documented) | through ε_v, η, φ |
| `profile_file` | the Enceladus profile | | porewater and ocean P, porewater T | PlanetProfile (Styczinski et al. 2023) | sourced | other profiles |
| ocean temperature | 274.6 K (1.45 °C) | K | ocean SOLUTION and seafloor steps | profile, top of the rock layer (274.566 K) | sourced | |
| `porewater_temp` | `"profile"` when `reactive_depth = "skin"`, else `"ocean"` | °C, `"profile"` or `"ocean"` | porewater reaction temperature, the skin's own profile mean by default | the profile, a fixed number is an explicit override or a bound (160 °C, the hot bound, in the sweep) | decided | the 160 °C bound |
| `organic_wt_percent` | 0 (production 15) | wt% | IOM in the rock | deliberate, to be swept | | yes |
| `carbon_mode` | `"c4_only"` | | carbon redox treatment and database | §6 | decided | `"coupled"` as the bound |
| `pore_gas_fate` | `"to_ocean"` | | the pore bubble's fate | §6, PROGRESS Part 7 | decided, `"retain"` is the other bound | yes (the bracket) |
| `initial_ocean` | `"pure_water"` | | ocean start, `cometary_ice_ocean()` option | attribution, Miller et al. 2026 | sourced | the cometary option |
| `cometary_ice_ocean(v)` | 0.02 | | volatile depletion of 67P ice | Miller et al. 2026 (from Titan's N budget) | sourced, not for Enceladus | v |
| `include_NH3` | TRUE | | NH₃ in the cometary ocean (0.67 per 100 H₂O × v) | Rubin et al. 2019 Table 2, and Miller et al.'s text applies v to it | decided (Miller et al.'s Table 5 omits it) | |
| `cometary_ice_ocean(f)` | 1 | | fraction of the accreted cometary volatiles retained in the ocean, applied after v | | | f |
| `include_less_volatile` | FALSE | | ices less volatile than NH₃ at 273 K, undepleted (× f) | Rubin et al. 2019 Table 2; Fray & Schmitt 2009 or boiling point; `data/comet67p_volatility.csv` | option | yes |
| `trace_C_molal` | NULL (production 1e-14) | mol/kgw | numerical C seed | smallest value tested to converge | numerical | |
| `mode`, `grouping_factor`, `substeps` | `"ocean_circulations"`, none, 1 | | orbits per step | grouped mode compounds exchange: x_step = 1 − (1 − x_tide)^G | | |
| `rec_eqb_minerals` | FALSE | | promote fast primaries to equilibrium (bulk_kinetic) | a modelling choice | decided | |
| mineral SSA (`Simulator(ssa =)`) | `"legacy"`: 0.1 (fo, en, di, teph), 5 (ab, an, Kfs, tro), or `"geometric"`: 6/(ρ_i d_i) × `roughness` from `grain_diameter` | m²/g | reactive surface area | legacy: none, geometric: grain size (§6, Permeability and surface area from grain size) | legacy **UNSOURCED**, geometric an option | yes |
| Fe, Ni metal | `dissolve_only` equilibrium | | corrodes completely in step 1 | §6 (Fe/Ni metal) | decided (an accepted idealisation) | |
| schreibersite rate | 1.370e-7 | s⁻¹ | bulk first-order | one 1-day lab measurement (RECONSTRUCTION_STATUS §4, §8) | placeholder | |
| IOM kinetics | A = 2e15 s⁻¹, distributed Ea | | IOM channels | Burnham 2019; Miller et al. 2025 | sourced (`Ea_N` a one-point fit, `IOM_S` uncalibrated) | Ea width, CHn m0 |
| IOM density | Delarue Eq. 10 at the residual IOM's H/C (1.320 at t = 0) | g/cm³ | `bulk_density()` only | Delarue et al. 2026 | sourced, a fixed number is an option | |
| `ch4_cutoff_C` | 150 | °C | only for `carbon_mode = "temperature"` | inherited default | | |

## 6. Model choices and justifications

### Carbon mode (`c4_only`)

Coupled carbon is C(+4) only. Reduced carbon enters only as decoupled methane (`Mtg`) from the IOM,
with IOM_CHn as `Mtg 1 H -2.2` (a 2.2 e⁻ reduction imposed by hand). In `"coupled"` runs equilibrium
reduces carbonate with H₂, so they are an upper bound on abiotic methanation, not the expected case
(PROGRESS Part 3).

### IOM composition

The composition is Miller et al. (2026) SI Table 9, from Alexander et al. (2007) CR-chondrite IOM, so
it combines a CR composition with CM (Murchison) kinetics. The sulfur provenance is open.

### IOM kinetics

The IOM reacts as distributed-Ea pools (the Vitrimat approach of Burnham 2019), with A = 2e15 s⁻¹,
calibrated to the heating experiments of Miller et al. (2025). `Ea_N` = 242.78 kJ/mol is a one-point
fit and its m0 0.875 a floor, and `IOM_S` has no kinetic calibration. The 350 °C residue check passes
for H, but O is under-released by 13–30% (PROGRESS Part 4).

### Rock

The rock is the Scenario R reconstruction (`reconstruction/RECONSTRUCTION_STATUS.md`) on the CI
composition of Lodders et al. (2025). Lawrencite and the 15 wt% organics are deliberate, and the
organics will be swept.

### Mineral rates

The silicates, feldspars and troilite (as hexagonal pyrrhotite) use Palandri & Kharaka (2004).
Tephroite uses log k(H⁺) = −5.4 at pH 2 and n = 0.47 (Pincus et al. 2026, Table S1, citing Casey et
al. 1993), with an Ea that falls with pH. That Ea is a fit to Casey et al.'s Fig. 4 (54, 26, 22 and 4
kJ/mol at pH 2.5–5.15, 0 taken at pH 6.2, and zero above pH 5.85). Tephroite has no neutral branch, so
the rate at pH > 7 stays at its pH-7 value.

### Surface areas (UNSOURCED)

Forsterite, enstatite, diopside and tephroite are at 0.1 m²/g, and albite, anorthite, K-feldspar and
troilite at 5 m²/g, so the feldspars react against 50× the area of the mafic silicates.

### Flow driving

The flow driving is UNSOURCED as a single value and is to be swept. The Simulator's Darcy flux is
q = (k/μ)·ΔP/L, with ΔP·sin(2πt/P) the tidal pore pressure, L = `layer_thickness` and μ = 1.8e-3 Pa·s,
so only k·ΔP/L matters. ΔP is derived from the strain (see Exchange below), whereas the table in this
section uses a fixed ΔP. The exchanged fluid per orbit is |q| integrated over the orbit times the
core's surface area (q is per unit total area). Every orbit's volume is treated as fluid that leaves
the rock and is replaced by ocean water, so **the oscillating flow is assumed to give net exchange**
(an open question, below). Production uses k = 1e-7 m² and L = 194.2 km, and the table uses ΔP = 1e5
Pa. The literature constrains ranges, each with its own limitation:
- Liao et al. (2020), poroviscoelastic tidal model of the core. They consider k from 1e-12 to
  1e-7 m². Their Fig. B3 has a maximum pore pressure of ~3×10⁴ Pa at the poles and ~0 at
  mid-latitudes, with mean pore velocities of 3.8e-6 m/s (radial) and 2.2e-6 m/s (tangential). That
  is one case (η = 1e16 Pa·s, μ = 1 GPa, α = 0.1, κ = 1e-8 m², φ = 0.2). The 1e5 Pa default is ~3×
  their plotted peak.
- Choblet et al. (2017), 3-D thermal convection in a tidally heated porous core. They explore k
  from 1e-15 to 1e-12 m² and prefer 1e-14 to 1e-13 m², with φ 0.1–0.3 and 10–30 GW of tidal power.
  Upwellings exceed 363 K, and the ocean is processed through the core above 363 K in 25–250 Myr.
  The limitation is that the convection is steady and buoyancy-driven, with no tidal pumping and no
  chemistry.
- Zandanel et al. (2021), 1-D reactive-transport columns of olivine only, with water velocities
  of 1, 10 and 100 cm/yr (which they take from Choblet et al. for k from 1e-15 to 1e-12 m²), at 1, 25,
  90 and 150 °C. Olivine alteration lasts up to ~75 Myr, and H₂ production (> 1e12 mol/yr) is far
  above the ≥ 1e9 mol/yr inferred from the plume. Their 0.1 mm olivine grains have a geometric
  surface area of ~0.01 m²/g, 10× below ours, and their short lifetimes support a spin-up.
- Fisher et al. (2024), 3-D outcrop-to-outcrop hydrothermal siphons scaled to ocean-world
  gravity. Aquifer k is 1e-13 to 1e-10 m², and a siphon at Enceladus's gravity needs ~1e-11 to 1e-10
  m². W/R is 4.9–150 per Gyr at Enceladus's gravity. The limitation is that it is an Earth
  ridge-flank analogue covering only the upper 4–5 km.

Production k = 1e-7 m² is the top of Liao et al.'s range, 3–6 orders of magnitude above Fisher et
al.'s aquifers and ~6 above Choblet et al.'s preferred range. The table gives its flow as the code
computes it (`Simulator$run_simulation()` with each k and `max_tidal_stress`):

| k (m²) | ΔP (Pa) | mean \|q\| (cm/yr) | fluid per orbit / porewater | pore displacement amplitude | ocean circulation time |
|---|---|---|---|---|---|
| 1e-15 | 1e5 | 5.8e-4 | 1.0e-12 | 1.7e-8 m | 7.8e9 yr |
| 1e-14 | 1e5 | 5.8e-3 | 1.0e-11 | 1.7e-7 m | 7.8e8 yr |
| 1e-13 | 1e5 | 5.8e-2 | 1.0e-10 | 1.7e-6 m | 7.8e7 yr |
| 1e-12 | 1e5 | 0.58 | 1.0e-9 | 1.7e-5 m | 7.8e6 yr |
| 1e-10 | 1e5 | 58 | 1.0e-7 | 1.7e-3 m | 7.8e4 yr |
| 1e-8 | 3e4 | 1.7e3 | 3.1e-6 | 5.1 cm | 2.6e3 yr |
| 1e-7 | 1e5 (production) | 5.8e4 | 1.0e-4 | 1.7 m | 78 yr |

The circulation time is the ocean's mass over the exchange per orbit, as
`sim$get_summary_metrics()$circulation_time_years` reports it. Everything scales with k·ΔP, so ΔP =
1e4 Pa is 10× slower. At k = 1e-8 m² and 3e4 Pa the displacement, 5 cm, is close to the ~7 cm implied
by Liao et al.'s 3.8e-6 m/s. Choblet et al.'s preferred k, with ΔP 1e4–1e5 Pa, gives circulation times
of 7.8e7–7.8e9 yr against their 25–250 Myr. Zandanel et al.'s 1–100 cm/yr corresponds here to k ≈
2e-12 to 2e-10 m² at 1e5 Pa (2e-11 to 2e-9 at 1e4 Pa). The 160 °C example case (k = 1e-7 m²)
circulates the ocean ~10⁵–10⁶× faster than Choblet et al. This table is the Darcy term alone, with the
gradient over the whole layer. The exchange itself is the smallest of three volumes, with the Darcy
gradient over the skin depth (see Exchange below), and at production k the poroelastic term binds.

### The sweep (designed, not run)

1. Treat the driving term as one quantity, reported three ways: k·ΔP/L, the mean Darcy flux, and
   the fluid exchanged per orbit as a fraction of the porewater. k, `max_tidal_stress` and
   `layer_thickness` already set it.
2. Ranges: ΔP 1e4–1e5 Pa (Liao et al.), and k 1e-15–1e-7 m², marking the Choblet (1e-14–1e-13),
   Fisher (1e-13–1e-10) and Liao (1e-12–1e-7) sub-ranges. This brackets Zandanel et al.'s 1–100 cm/yr.
3. Report each case's ocean circulation time, for comparison with Choblet et al.'s 25–250 Myr.
4. Open question: how the oscillating flow becomes net exchange. Liao et al.'s velocities imply
   ~7 cm of pore-fluid displacement per tide, back and forth. The code counts |q| over the whole orbit
   as one-way exchange, which is an upper bound if back-and-forth flow brings the same fluid back. L,
   the layer thickness (194.2 km), sets the pressure gradient, whereas Liao et al.'s pressure varies
   with latitude over the core's radius.

### Exchange: the smallest of three volumes per half-orbit

The core behaves as a sponge. Tides squeeze and expand the matrix, expelling and drawing in pore water
(the poroelastic coupling of Liao et al. 2020). Per half-orbit the exchanged volume is ΔV =
min(ΔV_poro, ΔV_Darcy, V_pore), over the reacting box (the skin, below):
1. **ΔV_poro = α·ε_v·V_box**, the fluid the tidal volumetric strain expels (drained limit), where α is
   the Biot coefficient (below).
2. **ΔV_Darcy**, what drains through the rock in half an orbit. This is half the |q|·A integral over
   the orbit, i.e. (k/μ)(ΔP/L_eff)·A·T/π, with the pressure gradient over L_eff = min(L, δ), the skin
   depth below. (Written with the amplitude, (k/μ)(ΔP/L)·A·T/2, it would be π/2 larger. The code uses
   the integral of the sinusoid.)
3. **V_pore = φ·V_box**, saturation.

The exchange fraction is x = ΔV/V_pore per half-orbit. Per orbit the two halves compound, x_orbit =
1 − (1 − x)², and the grouped steps compound that again, x_step = 1 − (1 − x_orbit)^G. Since x ≤ 1 by
construction, no separate cap at the pore volume is needed. `sim$exchange` has all three volumes and
the binding limit (`poro`, `darcy` or `pore_volume`), and `get_summary_metrics()` reports
`exchange_fraction_per_half_orbit`, `exchange_limit` and the circulation time (ocean mass ÷ exchange
per orbit). Every generated input stamps header lines with the reacting box and the exchange. They are
constant within one generated input and change between the chunks of a run whose φ evolves (Long runs,
below).

### Porosity and the Biot coefficient

The initial porosity is `phi0` (Simulator, with `porosity` the older name), default 0.32, the
PlanetProfile profile's `phiRockMax`. This is **one fit** (Styczinski et al. 2023) with its own
assumptions, and the sweep covers 0.1–0.4 (Choblet et al.'s 0.1–0.3 and the range the core density
implies). φ controls k (as φ³ under Kozeny–Carman), α, the storage, δ and the strain-limited volume.

The Biot coefficient comes from the critical-porosity model (Nur, Mavko, Dvorkin & Galmudi 1998, *The
Leading Edge* 17, 357–362, as cited by Dvorkin & Nur 2000). The frame's bulk modulus falls linearly to
zero at the critical porosity φc, K_frame = K_grain(1 − φ/φc), so **α = 1 − K_frame/K_grain = min(1,
φ/φc)**, which is ≥ φ as a porous solid requires. `phi_c` defaults to 0.4, the critical porosity of
sorted grain packs and sandstones in Dvorkin & Nur (2000)'s table after Nur et al. (1998) (36–40%, with
40% for sandstones and carbonates and 65% for chalk). It is still a placeholder. **φ ≥ φc means α =
1.** The solid no longer forms a frame there and the "rock" is a suspension (granular media are loose
at φc, Dvorkin & Nur 2000), so α is capped. `alpha = <number>` overrides it, and an override below φ is
an error. Bagheri et al. (2025) give a single, fixed core bulk modulus (100 GPa, their Table 2), not
frame and grain moduli, so α cannot be computed from them. At φ = 0.32, α = 0.80 (0.64 with φc = 0.5).
The box's solid mass is V·(ρ_core − φ·ρ_water), so the W/R follows φ.

### The core-density check

ρ_bulk = (1 − φ)ρ_solid + φρ_fluid, and Bagheri et al. (2025) give a present-day core density of
2160–2400 kg/m³. `phi_from_core_density(rho_solid, rho_core)` inverts it. At t = 0 with 15 wt% IOM
(ρ_solid 2.920 g/cm³), Bagheri's range implies φ = 0.27–0.40. Bagheri's value is present-day, so
**it is a check on the end state, not an input**. `bulk_density()` adds a `core_density_ok` column and
warns when the final ρ_bulk is outside 2160–2400. φ and the organic content are linked, because less
IOM means a denser solid and, for the same core density, a higher φ (a forsterite-only solid, 3.21
g/cm³, needs φ = 0.37–0.48).

### The strain: a layered forward model (`strain_model = "layered"`, default)

The core's tidal response is computed in the coupled body, not for an isolated core. The body is
Bagheri et al. (2025)'s structure, an ice shell, an ocean and a viscoelastic core, homogeneous and
compressible (bulk moduli 10, 2.1 and 100 GPa, their Table 2), with Andrade rheology in the shell and
the core (α = 0.2, β = µ^(α−1)η^(−α), their Appendix A), an inviscid static ocean, e = 0.0047, n =
5.307e-5 rad/s, R = 251.99 km and M = 1.08101e20 kg (their Table 1). The tidal radial functions come
from TidalPy (J. P. Renaud, version 0.8.1, its shooting-method radial solver).
`scripts/layered_core_strain.py` builds the models and writes the tables (`python -I
scripts/layered_core_strain.py`, ~40 s, needs numpy and TidalPy). Parameters are sampled
log-uniformly (as in Bagheri's MCMC) within the Table 5 ranges of each scenario, with the core density
set by the total mass (samples outside 2160–2400 kg/m³ rejected).

*Validation against Bagheri et al. (2025) Table 3.* The Table 5 ranges were sampled with the core
viscosity over the scenario's range, keeping, as Table 3 does, the models with 18–28 GW of tidal heat
and ≥ 90% of it in the scenario's layer. The core's share comes from the energy flux through the core
surface, r²·Im[y₁*y₂ + 6y₃*y₄ + y₅*y₆/4πG], which equals the Eq. 5 total at the surface to 1e-15 and is
continuous across the ocean. The table gives mean ± 3 sd:

| | dissipative core | Table 3 | dissipative shell | Table 3 |
|---|---|---|---|---|
| \|k₂\| | 0.0243 ± 0.0122 | 0.0240 ± 0.0102 | 0.0320 ± 0.0179 | 0.0317 ± 0.0130 |
| \|h₂\| | 0.0464 ± 0.0325 | 0.0459 ± 0.0294 | 0.0867 ± 0.0418 | 0.0848 ± 0.0359 |
| \|l₂\| | 0.0105 ± 0.0080 | 0.0104 ± 0.0072 | 0.0204 ± 0.0105 | 0.0200 ± 0.0090 |
| k₂ lag | 40.3° ± 21.5 | 39.35° ± 19.72 | 29.4° ± 18.6 | 28.91° ± 17.25 |
| h₂ lag | 3.10° ± 3.47 | 2.95° ± 2.69 | 29.3° ± 18.5 | 28.73° ± 17.18 |
| l₂ lag | 3.25° ± 3.68 | 3.10° ± 2.85 | 30.3° ± 19.1 | 29.71° ± 17.89 |

Every mean is within 5% of Table 3. The spreads are somewhat wider because the libration and J₃
constraints of their MCMC are not applied. (A static or a dynamic ocean changes k₂ by < 1%, and an
incompressible model lowers |h₂| by 5%.)

*Bagheri et al.'s Eq. 8 is not used.* As printed, Im(k̂₂^core) = 4πGρ_core R_core⁴ / (5R_E³) · Im(ΔR_core)⁴
(and their B4) is not dimensionally consistent, because the right side has a displacement to the fourth
power and is not a pure number. The core's share of the heat comes instead from the time-averaged energy
flux through a sphere of radius r, F(r) ∝ r²·Im[y₁*y₂ + 6y₃*y₄ + y₅*y₆/(4πG)], which is the power
dissipated inside r. It equals their Eq. 5 total at the surface (−5R·Im k₂/(4πG), to 1e-15), is
continuous across the inviscid ocean (to 1e-13) and vanishes for an elastic core (1e-4 of the total at
η = 1e25, Andrade's residual). The core's share is F(R_core)/F(R_E) (`scripts/layered_core_strain.py`
documents the same).

*Outputs* (`data/layered_core_strain.csv`, `data/layered_core_strain_summary.csv`, 500 samples per
viscosity). For each scenario and core viscosity (dissipative core: 8 values over 4.4e11–2.3e12 Pa·s,
dissipative shell: 3e13, 1e14 … 1e22), the tables give these, each as mean ± sd:
- the pole amplitude of the **volumetric strain** ∇·u = dy₁/dr + (2y₁ − 6y₃)/r (times U_pole =
  (3/2)(nR)²e) against depth below the core surface, 0–180 km;
- the core-surface radial displacement and the ocean-top displacement;
- the seafloor normal stress and the Love numbers.

The tables are computed at K_core = K_d = 9.798 GPa (below). Over the top 100 km:

| scenario, η (Pa·s) | ∇·u at 0 km | 50 km | 100 km | u_r core surface (pole) |
|---|---|---|---|---|
| core, 4.4e11 | 7.81e-8 ± 0.86e-8 | 4.52e-8 | 2.11e-8 | 2.15 ± 0.24 m |
| core, 1e12 (placeholder) | 8.00e-8 ± 0.87e-8 | 4.52e-8 | 2.03e-8 | 0.97 ± 0.11 m |
| core, 2.3e12 | 7.93e-8 ± 0.87e-8 | 4.49e-8 | 2.02e-8 | 0.43 ± 0.05 m |
| shell, 3e13 | 8.12e-8 ± 0.92e-8 | 4.59e-8 | 2.06e-8 | 0.035 ± 0.004 m |
| shell, 1e16–1e22 | 6.9–7.3e-8 | 3.9–4.1e-8 | 1.8e-8 | 0.002–0.003 m |

The volumetric strain does not follow η or the scenario, because the core's bulk modulus sets it. The
core-surface displacement follows η (2.2 to 0.4 m over the dissipative-core range), but the deformation
is almost all shear. The dilatation is the mean-stress change, ~900 Pa, over K_core, and scales as
1/K_core (checked over 2–100 GPa to 0.2%). The displacements and Love numbers do not change with
K_core, and the Table 3 validation is identical at 9.8 GPa and at Bagheri's 100 GPa. Liao's ε₀ =
2u_r/R is a linear strain, ~1e-5 for these displacements too, whereas the fluid is expelled by the
dilatation.

### The core's bulk modulus: the drained frame, K_d = K_grain·(1 − φ/φc)

The core's bulk modulus is the drained frame's, not Bagheri's 100 GPa (an intact-rock value the Love
numbers do not constrain). The same factor gives μ_d = μ_grain·(1 − φ/φc), checked against Bagheri's
Table 5 core rigidity (4–70 GPa). K_grain and μ_grain are Hill averages (mean of the Voigt and Reuss
bounds), by volume, over the reconstruction's minerals and the IOM (`grain_elastic_moduli()`,
`data/mineral_elastic_moduli.csv`):

| phase (vol% of the minerals) | K (GPa) | μ (GPa) | source |
|---|---|---|---|
| forsterite (26.2), enstatite (26.5), diopside (6.8), anorthite (2.3), low albite (11.9), orthoclase for K-feldspar (0.8) | 127.3, 107.9, 116.5, 81.9, 60.1, 58.3 | 81.6, 77.9, 72.8, 36.0, 36.0, 28.1 | Abers & Hacker (2016), G-cubed 17, 616–624 (their database, isothermal K) |
| α-iron (6.2) | 166.7 | 81.5 | Bass (1995), AGU Reference Shelf 2, Table 1 |
| tephroite (0.5) | 128 | 54 | Bass (1995) Table 7 |
| magnesiochromite (0.6) | 203.3 | 104.9 | **proxy:** chromite FeCr₂O₄, Bass (1995) Table 1 (MgCr₂O₄ not listed) |
| troilite (16.9) | 65.8 (K′ 4.17) | **31.5, an estimate** (27–35.5) | K₀: Holland & Powell data via Perple_X 7.2.1 `hp634ver.dat`, cited as found. The file's title is "Holland & Powell, JMG 2011 (TC-DS633, 23 Jun 2017) + additional data", entry `tro`, b6 = 658000 bar, b8 = 4.17 (trot, lot, trov identical). μ: no measurement found, so K₀ × the μ/K of Bass (1995)'s monosulfides, sphalerite 0.41 and galena 0.54 |
| IOM (33.3% of the solid at 15 wt%) | 10.2 (K′ 7.7) | **3 (1–5)**, terrestrial-kerogen analogue | K₀: Delarue et al. (2026), EPSL 691, 120227, Eq. 7 / Fig. 1 (a Murnaghan fit to an anthracite coal in a diamond-anvil cell to 7 GPa, after which their text assumes 10 GPa for low-temperature carbonaceous matter). μ: terrestrial kerogen. Yan & Han (2013), SEG Expanded Abstracts (doi:10.1190/segam2013-1319.1), invert K = 3.5–5.0 GPa for Green River kerogen and speculate μ ≈ K/2 (1.75–2.5 GPa, μ was not measured). Othman et al. (2022), ACS Omega 7, 32829–32839, tabulate nanoindentation E = 1.9–2.2 GPa (porous) to 6–15 GPa (non-porous kerogen), i.e. μ ≈ 0.7–6.3 GPa at ν 0.2–0.3, and atomistic type II kerogen μ = 4–7.5 GPa. Their table quotes Yan & Han's μ as half of Young's modulus, where Yan & Han wrote half of the bulk modulus. The 1–5 GPa bracket spans the experimental values, and the simulations reach higher |

Ni (0.7%) and schreibersite (0.5%) have no sourced modulus and are left out, with the minerals
renormalised (98.8% of their volume is covered). The IOM's volume is the remaining IOM at the Delarue
et al. density (1.320 g/cm³ for the bulk IOM), beside the reconstruction's mineral volume, as
`bulk_density()` computes it. Results at φ = 0.32, φc = 0.4 (frame factor 0.2):

| | K_grain | μ_grain | K_d | μ_d (Table 5: 4–70 GPa) |
|---|---|---|---|---|
| minerals only, without troilite (for comparison) | 107.5 | 67.8 | 21.5 | 13.6 |
| minerals with troilite | 99.0 | 59.5 (57.9–60.8) | 19.8 | 11.9 (11.6–12.2) |
| with the IOM (default) | 49.0 | 26.1 (22.9–28.9) | 9.80 | 5.21 (4.57–5.79), inside |

The ranges are the troilite and IOM shear moduli at their bracket ends. The IOM halves the grain
moduli, and with it the Voigt and Reuss bounds are far apart (K 72.8 and 25.2 GPa), so the Hill average
is a coarse estimate.

K_d follows φ(t) and the IOM. Each Simulator takes the IOM's share of the solid volume, either from
`iom_wt_percent` (default 15, the production value) or from `iom_vol_frac`. `iom_wt_percent` is
independent of PHREEQC's `organic_wt_percent`, so that a run and its organic-free blank share one
exchange. `run_tidal_evolution()` passes each chunk's remaining IOM and φ, and the chunk table reports
`K_d_Pa`, `mu_d_Pa` and `iom_vol_frac`. **At φ ≥ φc there is no frame (α = 1).** The Simulator then
stops with an error, and `run_tidal_evolution()` flags the run (`stopped = "phi reached phi_c"`) and
returns the chunks so far, so K_d never goes to zero. `config$core_bulk_modulus = <number>` overrides
K_d.

In the code, ε₀ is the pole value at the core surface (interpolated in log η between the tabulated
viscosities of `tidal_scenario`, and an error outside that range), times K_table/K_core. **ε_v =
0.19245 × the r²-weighted mean over the reacting box of the pole amplitude against depth**, kept at
its 180 km value below that. `strain_scale` multiplies it (an override, default 1).
`tests/exchange_check.R` sections 7 and 8 check the interpolation, the mean, the load, the grain moduli,
the 1/K scaling, the override and the stop at φc, and `tests/porosity_evolution_check.R` checks the
per-chunk K_d.

### The Liao strain (`strain_model = "liao"`, for reproducing earlier runs)

ε₀ is Liao et al. (2020)'s polar strain amplitude (§2.2, Fig. 2a):
ε₀ = (9/4π)·e·ω²/(ρG)·(5/3)·|k₂*|, k₂* = (3/2)/(1 + 19μ*/(2ρgR)), with the Maxwell rigidity μ* =
μ·iωτ/(1 + iωτ), τ = η/μ, g = (4/3)πGρR. Their §2.2 text and Appendix A write k₂ with an extra factor
3/5, but the Fig. 2 caption's form (without it) is the one that reproduces Fig. 2a. It gives the elastic
contours at μ = 2.1e7, 8.0e7, 8.6e8 and 8.6e9 Pa against 2.1e7, 7.5e7, 8e8 and 8.6e9 read off the
figure, and the 10⁻⁵ contour meets μ = 1e10 at η = 1.6e12 as in the figure (with their R₀ = 186 km,
33 h, e = 0.0047, ρ = 2400 kg/m³). The 3/5 version is ~1.6× off. `tests/exchange_check.R` checks this.

ε₀ is a polar *radial* strain, whereas the fluid is expelled by the *volumetric* strain. Liao et al.'s
simplified model does not give the dilatation (with a Laplace displacement potential it is zero in
their viscoelastic end-member). **Under `"liao"` it is assumed** that the dilatation has the degree-2
pattern (r/R)²·P₂(cos θ) with amplitude ε₀ at the pole. A degree-2 field has zero mean, so each
half-orbit an equal volume compresses (expelling fluid) and expands (drawing it in). The compressed
part's volume mean is ε_v = 0.19245 × (3/5)(1 − x⁵)/(1 − x³) × ε₀, x = (R − L)/R, where 0.19245 is half
the surface mean of |P₂| (0.385, the compressed half of the pattern) and the rest is the radial mean of
(r/R)² over the box. The factor follows the reacting box. With `reactive_depth = "skin"` it is
evaluated over the skin, 0.192 for a thin skin (r ≈ R, 0.192 at k ≤ 1e-11, where δ < 1 km) falling to
0.138 for the 74 km skin at k = 1e-7. With `"layer"` it is evaluated over the whole layer, 0.115 for
the full core. `sim$exchange$strain_factor` reports it. If the matrix is nearly incompressible the
dilatation is smaller, so ΔV_poro is an upper estimate.

The rigidity and viscosity come from Bagheri et al. (2025) Table 5, "dissipation in the core" (µ 4–70
GPa, η 4.4e11–2.3e12 Pa·s, the scenario consistent with tidally driven flow in the core). The defaults
are **placeholders at the geometric means**, `core_rigidity` = 17 GPa and `core_viscosity` = 1e12
Pa·s, until values are chosen. Under `"liao"`, at these viscosities the Maxwell time η/µ is seconds to
minutes against a 33 h tide, so the matrix responds viscously, |µ*| ≈ ωη, and ε₀ depends on η, not on
µ:

| η (Pa·s) | ε₀ | ε_v (full core) | ΔV_poro / V_pore per half-orbit (= ε_v/φc) |
|---|---|---|---|
| 4.4e11 | 3.8e-5 | 4.4e-6 | 1.1e-5 |
| 1e12 (placeholder) | 1.7e-5 | 2.0e-6 | 5.0e-6 |
| 2.3e12 | 7.6e-6 | 8.7e-7 | 2.2e-6 |

These apply to any µ in 4–70 GPa, with a full-layer box. With α = φ/φc the ratio α·ε_v/φ = ε_v/φc does
not depend on φ. Bagheri's rigidities lie above Liao's explored 0.01–10 GPa, and Liao's viscosities
(1e12–1e18) only touch Bagheri's at 1–2.3e12. The overlap (µ 4–10 GPa, η 1–2.3e12) is where the viscous
knee above was checked. Bagheri's core radius and density range (183–213 km, 2160–2400 kg/m³) move ε₀
between 1.4e-5 and 2.1e-5 at the placeholders, and the profile's 194.2 km and 2353 kg/m³, which the
code uses, are inside both ranges. The pore-volume limit cannot bind at real porosities, because even a
fluid matrix (k₂ → 3/2) gives ε₀ ≈ 1.5e-4, so α·ε_v ≤ 1.7e-5 of the box per half-orbit, against φ =
0.1–0.4.

### The tidal pore pressure

The pressure that drives the Darcy term is not an independent parameter. It is the undrained pore
pressure of the same strain, **ΔP = α·ε_v·M**, with the Biot modulus M = 1/S = (K_u − K_m)/α² of Liao
et al. (2020) Table A1 (their Eq. A5 with no change of fluid content). `max_tidal_stress = <number>`
overrides it, for reproducing earlier runs or for a sweep. Under the layered strain (default, K_d = 9.8
GPa), the strain's ΔP is 64 Pa in a thin skin at the placeholders (47 Pa in the 69 km skin), against
the ocean load's 26 Pa (below). At 100 GPa it would be 6 Pa. Under `"liao"` at the placeholders (φ =
0.32, φc = 0.4, α = 0.80, M = 5.2 GPa, ε₀ = 1.73e-5), the pole amplitude α·M·ε₀ is 7.2e4 Pa, and the
box mean is 1.4e4 Pa in a thin skin (ε_v/ε₀ = 0.192), 1.0e4 Pa in the 69 km skin at k = 1e-7 and 8.3e3
Pa over the whole layer (0.115).

This falls well short of Liao et al.'s ~3e4 Pa (Fig. B3a). The same formula at their Fig. B3 case (η =
1e16 Pa·s, µ = 1 GPa, α = 0.1, φ = 0.2, ε₀ = 8.6e-7) gives only 1.1e3 Pa at the pole, ~28× below their
plot. The ocean's tidal loading of the seafloor (below) is far too small to close that gap, which
remains unexplained (their full poroviscoelastic solution is not reproduced here). The strain-derived
ΔP does not change the exchange where the strain binds.

### Ocean tidal loading

The ocean column over the seafloor changes thickness with the tide, and its weight change, ΔP_load =
ρ_w·g_c·|u_ocean − u_core|, loads the pores. The calculation uses Bagheri et al. (2025)'s eccentricity
tide (their Eq. 1, U₂ = r²ω²e[−(3/2)P₂⁰(cos θ)cos ωt + P₂²(cos θ)((3/4)cos ωt cos 2φ + sin ωt sin
2φ)], ω = n the mean motion, here 5.289e-5 rad/s from the 118,800 s period (their 5.307e-5), e =
0.0047) and the radial displacement u_r = h₂U₂/g (their Eq. 2), at the pole (U₂ = (3/2)r²ω²e):
- ocean surface: taken to move with the surface, whose whole-body Love number is Bagheri et al.'s
  Table 3, "dissipative core": |h₂| = 0.0459 ± 0.0294, lag 2.95° ± 2.69°, at R = 252.1 km, g = 0.114.
  This gives u_ocean = 0.50 m at the pole (their Eq. 4's maximum over the surface,
  12√(2/7)|h₂|R²ω²e/g, is 2.16 m, and 0.78–3.54 m over the ±);
- core surface: under `"liao"`, the core's own h₂ = (5/3)k₂* of the uniform Maxwell core (Liao's
  formula, the one behind ε₀). |h₂| = 0.29, lagging 83°, at the placeholders (0.64 at η = 4.4e11, 0.13
  at 2.3e12), so u_core = 1.7 m (3.7, 0.7 m);
- the difference is taken with the phases. ΔP_load = 213 Pa at the pole (459 and 108 Pa at the η range
  ends), and its box mean is the pole value times the strain factor: 41 Pa in a thin skin, 30 Pa in the
  69 km skin, 25 Pa over the layer.

Under the layered strain (default), both displacements come from the layered model (the ocean top and
the core surface, with their phases, and |u_ocean − u_core| averaged over the samples). At the
placeholders u_ocean is 0.57 m, u_core 0.97 m and ΔP_load 133 Pa at the pole (262 and 87 Pa at η =
4.4e11 and 2.3e12, and ~130 Pa in the dissipative-shell scenario, where the core hardly moves and the
ocean top moves 1.04 m). That is 26 Pa in a thin skin, against the strain's 64 Pa at K_d (6 Pa at 100
GPa), and without it the exchange would be 29% smaller. The layered model's whole normal stress on the
core surface is larger, ~880 Pa at the pole, because it includes the tidal potential's own pressure in
the ocean. That part acts on the pore fluid as a body force too and drives no flow, so the
column-thickness form is kept. This is an approximation, since the poroelastic split of a boundary load
between the matrix and the pore fluid is not modelled.

The load is **added to the strain-derived ΔP as an amplitude**, which is an upper bound because the
phase between the loading and the matrix strain inside the core is not known. The volume it stores,
S·ΔP_load·L_eff·A, is added to the strain volume, so the "poro" limit is ΔV ≈ S·(ΔP_strain +
ΔP_load)·δ·A in the skin. `ocean_load = FALSE` drops it. Under `"liao"`, at every k and η tested the
load is 0.3% of the strain term (strain 1.38e4 + load 41 = 1.38e4 Pa in a thin skin at the placeholders,
1.00e4 + 30 Pa at k = 1e-7, and 8,257 + 25 Pa over the layer).

### The exchange under the layered strain (K_d)

Per half-orbit, in the skin, at the placeholders:

| K_core | ε₀ | k = 1e-14: exchange, circulation | k = 1e-11 | k = 1e-7 |
|---|---|---|---|---|
| 100 GPa (Bagheri) | 7.8e-9 | 1.93e-8, 6.2e8 yr | 1.93e-8, 2.0e7 yr | 1.93e-8, 2.9e5 yr |
| K_d, minerals only (19.8 GPa) | 4.0e-8 | 3.46e-8, 3.5e8 yr | 3.45e-8, 1.1e7 yr | 3.05e-8, 1.8e5 yr |
| K_d with the IOM (9.80 GPa, default) | 8.0e-8 | 5.40e-8, 2.2e8 yr | 5.39e-8, 7.1e6 yr | 4.49e-8, 1.2e5 yr |
| `"liao"` | 1.7e-5 | 8.3e-6, 1.4e6 yr | 8.3e-6, 4.6e4 yr | 6.1e-6, 908 yr |

At the default, the exchange is 6.8e-8 to 4.8e-8 over η = 4.4e11–2.3e12 (k = 1e-14) and 4.8–5.4e-8 in
the dissipative-shell scenario, and the strain limit ("poro") still binds in the skin. `strain_scale`
= 0.16 lowers it to 2.2e-8.

With ΔP derived, the three volumes are tied together. The storage volume S·ΔP·δ·A = α·ε_v·δ·A equals
the strain volume α·ε_v·V_skin to O(δ/R), and the Darcy volume is 2/f ≈ 2.3 times the storage volume
(exactly, since δ² = c/ω). In the skin the strain therefore binds at every k, by a factor 2.3–3.3, and
the exchange per half-orbit is ε_v·α/φ, independent of k. k sets only how deep the skin reaches.

### The skin depth

The tidal pressure change diffuses into the rock only a skin depth δ = √(c/ω) (Liao et al. 2020,
Table A1 and Fig. B4), with their poroelastic diffusivity c = (k/η_f)·(K_m + 4µ/3)(K_u − K_m) /
(α²(K_u + 4µ/3)), K_m = (1 − α)K_s, K_u = K_m + α²K_sK_f/(φK_s + (α − φ)K_f), K_f = 2.2 GPa and K_s =
10 GPa (their Table A1, `water_bulk_modulus`, `rock_bulk_modulus`) and µ the matrix rigidity. (K_u −
K_m)/α² is the Biot modulus M = 1/S, with S = φ/K_f + (α − φ)/K_s the storage (1.93e-10 Pa⁻¹ at φ =
0.32, α = 0.80), so c is k/(η_f S) times f = (K_m + 4µ/3)/(K_u + 4µ/3) = 0.88. Their convention has no
factor 2 (δ = √(c/ω), not √(2D/ω)), and this code follows it. δ is 22 m at k = 1e-14 m², 0.69 km at
1e-11 and 69 km at 1e-7. δ is the depth the exchange actually reaches, which is what makes the
well-mixed porewater box an approximation (below). As a cross-check, the volume the pressure change
stores in the skin, S·ΔP·δ·A, is the δ-based Darcy volume divided by 2/f = 2.27 at every k (exactly:
(ωT/π)/f).

The table shows which limit binds over the whole layer (`reactive_depth = "layer"`, φ = 0.32, α =
0.80, ΔP from the strain = 8.3e3 Pa, placeholders):

| k (m²) | δ | ΔV_poro (m³) | ΔV_Darcy (m³) | S·ΔP·δ·A (m³) | V_pore (m³) | binds | x per half-orbit | circulation time |
|---|---|---|---|---|---|---|---|---|
| 1e-14 | 22 m | 4.9e10 | 3.8e7 | 1.7e7 | 9.8e15 | darcy | 3.8e-9 | 1.1e6 yr |
| 1e-12 | 219 m | 4.9e10 | 3.8e8 | 1.7e8 | 9.8e15 | darcy | 3.8e-8 | 1.1e5 yr |
| 1e-11 | 692 m | 4.9e10 | 1.2e9 | 5.2e8 | 9.8e15 | darcy | 1.2e-7 | 3.3e4 yr |
| 1e-7 (production) | 69 km | 4.9e10 | 1.2e11 | 5.2e10 | 9.8e15 | poro | 5.0e-6 | 811 yr |

With the skin depth the Darcy volume grows as √k, not k. Over the whole layer the storage volume is only
δ·A/V_layer = 3δ/R of the strain volume, so Darcy binds until δ ≈ R/7, 28 km (k ≈ 2e-8 m²).

### The reacting box is the skin (`reactive_depth = "skin"`, default)

When δ < L, only the top δ of the layer exchanges in a tide, so the box that reacts is the top L_eff =
min(L, δ) (`sim$reactive_thickness()`). Its temperature and pressure are the profile means over the
top L_eff, its rock mass and pore volume are those of the L_eff shell, and the exchange fraction is
relative to the skin's pore volume (the poroelastic volume uses the skin's degree-2 factor).
`reactive_depth = "layer"` is the whole-layer box. In the skin the poroelastic limit binds at every k,
because the skin's pore volume grows with δ as the Darcy volume does, and the strain expels ε_v/φc ≈
8e-6 of it per half-orbit, less than what drains (2.3–3.3× less, above). Full core, α = φ/φc (φc =
0.4), placeholders:

| φ0 | α | k (m²) | L_eff | skin T | skin P | x per half-orbit | binds | circulation time |
|---|---|---|---|---|---|---|---|---|
| 0.32 | 0.8 | 1e-14 | 22 m | 1.4 °C | 66.7 atm | 8.32e-6 | poro | 1.4e6 yr |
| 0.32 | 0.8 | 1e-12 | 219 m | 1.4 °C | 66.7 atm | 8.31e-6 | poro | 1.4e5 yr |
| 0.32 | 0.8 | 1e-11 | 692 m | 1.4 °C | 66.7 atm | 8.29e-6 | poro | 4.6e4 yr |
| 0.32 | 0.8 | 1e-7 (production) | 69 km | 31.0 °C | 155.8 atm | 6.06e-6 | poro | 912 yr |
| 0.2 | 0.5 | 1e-14 | 28 m | 1.4 °C | 66.7 atm | 8.32e-6 | poro | 1.8e6 yr |
| 0.2 | 0.5 | 1e-12 | 284 m | 1.4 °C | 66.7 atm | 8.31e-6 | poro | 1.8e5 yr |
| 0.2 | 0.5 | 1e-11 | 899 m | 1.4 °C | 66.7 atm | 8.28e-6 | poro | 5.7e4 yr |
| 0.2 | 0.5 | 1e-7 | 90 km | 37.2 °C | 179.0 atm | 5.64e-6 | poro | 1,359 yr |
| 0.4 | 1 | 1e-14 | 19 m | 1.4 °C | 66.7 atm | 8.32e-6 | poro | 1.3e6 yr |
| 0.4 | 1 | 1e-12 | 192 m | 1.4 °C | 66.7 atm | 8.31e-6 | poro | 1.3e5 yr |
| 0.4 | 1 | 1e-11 | 606 m | 1.4 °C | 66.7 atm | 8.29e-6 | poro | 4.2e4 yr |
| 0.4 | 1 | 1e-7 | 61 km | 28.4 °C | 146.3 atm | 6.26e-6 | poro | 767 yr |

(The profile's grid is coarser than ~1 km near the top, so skins up to that depth take its top row.)
At a fixed k a larger φ gives a larger α and a shorter δ (the storage grows), a thinner, colder skin,
and a faster circulation (more pore water per volume exchanged at the same fraction). One-orbit
`tidal_kinetic` runs at each k (φ0 = 0.32, then α = 0.64, `porewater_temp = "profile"`, 15 wt% IOM)
completed with 0 errors. At low k, tidal pumping alone reaches only cold rock just below the seafloor
(1.4 °C, ocean pressure). Deep, hot rock is reached only as δ approaches L, i.e. at the top of the
permeability range. This is the comparison point with buoyancy-driven models (Choblet et al. 2017),
whose hot upwellings tidal pumping does not produce, and those models are out of scope here.
`porewater_temp` still overrides the profile temperature when given as a number.

### Porosity that evolves with the chemistry

`run_tidal_evolution()` (`porosity_evolution.R`) evolves the porosity, with `porosity_evolution = TRUE`
for production runs and FALSE in the quick self-tests. The frame is rigid, with no compaction: φ(t) =
1 − V_solid(t)/V_total, with V_total = V_solid(0)/(1 − φ0) fixed and V_solid from `bulk_density()`
(minerals, the remaining IOM at the Delarue density, clathrates, with `origin_pqi` passing the IOM
inventory across restarts). The run is split into chunks (`tidal_kinetic`, then
`tidal_kinetic_restart` with `stop_on_depletion = FALSE` and `allow_pressure_change = TRUE`, since the
skin's pressure moves with δ). Between chunks a new Simulator takes φ(t) and recomputes k (with
`k_of_phi`, e.g. Kozeny–Carman), α, the storage, ΔP, δ, the reacting box, the strain-limited volume and
the exchange fraction. A chunk whose φ changes by more than `phi_tol` (1% relative) is rerun with half
the orbits, down to one orbit, which is accepted and flagged (`over_tol`). The next chunk is scaled from
the last one's rate. φ ≤ `phi_min` (0.01) is capped and flagged `clogged`.
`tests/porosity_evolution_check.R` checks the bookkeeping.

### Long runs: grouped chunks and a geometric time step

With `t_end_years` (or `grouped = TRUE`), `run_tidal_evolution()` runs chunks of `steps_per_chunk` (4)
PHREEQC steps of G orbits each, with the exchange compounded as x_step = 1 − (1 − x_orbit)^G
(`tidal_kinetic`'s grouped mode), so the output is log-spaced in time. G starts at 1 and doubles per
accepted chunk, limited by accuracy:
- x_step ≤ 0.1 (G ≤ ln 0.9 / ln(1 − x_orbit), from each chunk's own exchange), which bounds the error
  of reacting the porewater for G orbits between exchanges;
- the porosity change per chunk within `phi_tol` (from the last chunk's rate), with a chunk that
  exceeds it redone at G/2.

There is no fixed cap (`g_max`, default none). The 1e5-orbit cap of the automatic
`"ocean_circulations"` grouping is not used, because the two criteria already bound the step (x_step ≤
0.1 allows ~1.7e6 orbits, 6,000 yr, at the reference). The run ends at `t_end_years`, at a **steady
state** (`ss_tol`: the ocean's C(+4) by actual oxidation state, pH and Na each change by less than
`ss_tol` relative over one ocean circulation time), at clogging or at φ = φc, and `stopped` says which.
Chunk directories have short names (`c0001`, ...) because Windows refuses paths over 260 characters, and
an output path over 120 characters is an error before any run.

`tests/grouped_validation_check.R` (~10 min) validates the grouping. The reference case to orbit 1,000
grouped (14 chunks, 56 steps, G up to 84) and ungrouped (1,000 steps) agree within 0.41% in the ocean's
C(+4), pH, H₂ and Na. At the reference the porosity sets the step (φ falls from 0.32 to ~0.17 by 5e4 yr,
by hydration), not x_step, and G reaches ~2.5e5 orbits (950 yr) by 5e4 yr, at 3–16 s per chunk.

### Step limits, set by validation

The 1,000-orbit comparison above cannot test the step limits, because there G stays below 100 orbits.
An ungrouped run long enough for x_step to reach 0.1 (G ~1.7e6 orbits) would take weeks, so each limit
was tested separately against much finer grouped steps (which the comparison above ties to ungrouped),
on the reference and on k ≈ 1e-11 (`scripts/step_limit_validation.R`, ocean C(+4), pH, H₂ and Na at
the end):
- x_step_max, with porosity kept fixed so that only x_step limits G (the regime of a clogged rock
  under `k_floor`), to 1e6 yr, against x_step_max = 0.01. Reference: 0.03 → Na −0.47%, pH −0.04%, H₂
  +0.04%, C(+4) +0.08%, then 0.05 → Na −0.96%, 0.1 → Na −2.1%, 0.3 → −6.7%, 0.5 → −11%, 1 → −23% (C(+4)
  +5.7%, pH −3.0%). The Na error is linear in the step (0.235% per 0.01 of x_step), so extrapolated to a
  zero step 0.03 is 0.7% off and 0.05 1.2%. k ≈ 1e-11: 0.03 → Na −0.65% (pH −0.03%, C(+4) +0.03%, H₂
  0.00%), 0.1 → Na −2.9%, 0.3 → −9.5%. The slope is steeper there (0.33% per 0.01), so extrapolated to
  zero 0.03 is ~1.0% off, at the limit and not inside it.
- phi_tol, with porosity evolving (x_step_max 0.1, not binding there), to 5e4 yr, against phi_tol =
  0.001. Reference: 0.003 → C(+4) +0.66% (pH, H₂, Na ≤ 0.10%), 0.005 → C(+4) +0.96%, 0.01 → +2.05%,
  0.03 → +7.6%, 0.05 → +10%. k ≈ 1e-11, against 0.003: 0.01 → C(+4) +2.2%, 0.03 → +3.7%, 0.05 → +4.3%
  (Na ≤ 0.9%). C(+4) is the sensitive output (at 3e-10 molal at the reference, and 3e-14 at k ≈ 1e-11,
  near the seed), whereas pH, H₂ and Na move by ≤ 0.9% at every phi_tol tried.

The cases use x_step_max = 0.03 and phi_tol = 0.003, the largest tested values within 1% both
against the finest run and extrapolated to a zero step (0.05 and 0.005 pass the first by 0.04% and fail
the second). The code defaults of 0.1 and 0.01 are 2% off. `run_tidal_evolution()` keeps those defaults
for the short self-tests, and production goes through the case files.

### Cost at the validated limits

Costs were measured on this machine, a Core Ultra 7 155H, while ~15 other runs shared it. The
reference to 5e4 yr took 807 steps (phi-limited, 1.5 h). Once porosity stops limiting, x_step 0.03
allows ~3.4e5 orbits (1,280 yr) per step at ~1.2 s per step, so 1e6 yr is ~800 steps (~20 min) and the
1e8 yr backstop ~78,000 steps (~1 day). At k ≈ 1e-11 a step covers ~1,060 yr at 9–20 s per step, so its
steady state (3τ ≈ 5e7 yr) takes ~5–11 days. At k ≈ 1e-14 a step covers ~1,450 yr at ~12 s, so the 1e8
yr backstop takes ~10 days. At the reference φ falls from 0.32 to 0.17 by 5e4 yr, so clogging may end
the run before a steady state. 4.5 Gyr is out of reach for every case at these limits (~40 days at the
reference alone).

With the phosphate phases (`"apatite"`), the reference to 1e4 yr took 134 chunks and 2.8 h. The first
chunk took 845 s (its .out has 109 "Numerical method failed with this set of convergence parameters"
and 112 "Maximum iterations exceeded" warnings), and later chunks took ~26 s for 4 steps (~6.5 s per
step, ~5× the cost without phosphates). Porosity still limits the step at 1e4 yr (φ = 0.20, with Δt
growing only as t^0.18). The projection to the 1e8 yr backstop is ~6 days if the step limit takes
over (x_step 0.03, ~1,280 yr per step) and ~24 days if porosity keeps limiting (the fitted Δt(t),
~81,000 chunks). 1e6 yr takes ~0.5 day either way. Clogging or a steady state may end a run earlier.
`-diagonal_scale true` in every step cut chunk 1 to 57 s but changed the converged results: 1,122 of
7,308 values differ by more than 1e-8 relative, hydroxyapatite is missing from one step (9.4e-6 against
0 mol), and dk_Diopside differs by 14%. It is therefore used only around the ocean step.

### Clogging and k_floor

The rigid-frame default ends a run when the solids fill the pores (φ ≤ `phi_min`, 0.01), and the run
reports the clog time and the state then (`clog_time_yr`, `end_reason = "clogging"`). `k_floor` (m²,
`k_floor_m2` in a case) is the alternative, a fracture-maintained permeability from reaction-driven
cracking. Kelemen & Hirth (2012, EPSL 345–348, 81–89) estimate the crystallization pressure of olivine
hydration and carbonation at > 300 MPa up to ~200 °C, above what is needed to fracture rock. They
conclude that "the volume change due to hydration and carbonation can cause fracture formation and
dilation, maintaining or increasing permeability and reactive surface", with permeability kept for tens
of thousands of years in Oman and at Lost City. They give no permeability value, so **the two floors
swept, 1e-14 and 1e-12 m², are a modelling choice, not values from the paper.** With `k_floor`, k =
max(k(φ), k_floor), and where the rock would clog the run goes on with φ kept at `phi_min` and k =
`k_floor`. The further solid volume is taken up by the cracking, not by the pores, and the `k_floor_on`
column marks those chunks. The run reopens (k from φ) if the solids shrink again. At the reference grain
size Kozeny–Carman gives 1.4e-12 m² at φ = 0.01, so both floors act only at clogging.

### Phosphate phases

Schreibersite oxidises to orthophosphate (0.040 mol P per kg of rock, all released within ~1 yr at the
reference). Without a phosphate sink, phosphate builds up to ~1 eq/kg in the porewater (hydroxyapatite
at SI +14) and, as Na phosphate, provides most of the ocean's alkalinity. Hydroxyapatite precipitates
slowly at low temperature, because P often precipitates first as amorphous phases that later
crystallise, and small organic acids inhibit Ca-phosphate precipitation (Hausrath et al. 2024,
*Minerals* 14, 591, §4.3 and refs. 205–207, 285–287). On the model's timescales (to 1e8 yr), however,
apatite belongs in the reference. Postberg et al. (2023, *Nature*) interpret Enceladus's plume phosphate
through the hydroxyapatite–calcite buffer (hydroxyapatite + 5 CO₃²⁻ + H⁺ ↔ calcite + 3 PO₄³⁻ + H₂O, with
one point fitting the whitlockite–calcite buffer better), and their experiments formed hydroxyapatite
and Ca-/Mg-whitlockite.

`phosphate_phases` is set in every generator and restart, in the case files and in RUNINFO, and a
restart refuses a parent with a different setting. `"apatite"` (default, the reference) puts
Hydroxyapatite and the precursors in the porewater and ocean lists. `"precursor"` (a sweep case) leaves
out the apatite, so the metastable phases control P (the Ostwald step rule). There is no "no phosphate"
option.

*The data.* The phosphate set comes from one source, Thermoddem V1.10 (BRGM,
`PHREEQC_ThermoddemV1.10_15Dec2020.dat`, database date 15 Dec 2020, shipped with PHREEQC 3.8.6). It is
cited as Blanc et al. 2012, *Appl. Geochem.* 27, 2107–2116, the paper BRGM's PHREEQC download page names
as describing the database. No other citation request could be found on the site, and the file's header
gives none. The set is converted to coreclath's basis (HPO₄²⁻) with Thermoddem's own H₂PO₄⁻ = HPO₄²⁻ +
H⁺ (log K −7.212, ΔH 3.6 kJ/mol, its analytic expression). `scripts/add_phosphate_phases.py` writes the
block into the coupled and uncoupled databases and comments out coreclath's own Hydroxyapatite and
Fluorapatite (`#tdP#`, kept for reference), and `make_c4only.py` passes it on to the C4-only database.
The two databases' aqueous P agrees (H₂PO₄⁻/HPO₄²⁻ within 0.007–0.06 log units over 0–150 °C), but
their **apatite differs by 3–5 log units**. Coreclath's (LLNL, from Robie, Hemingway & Fisher 1979,
USGS Bull. 1452, a compilation whose source for apatite was not traced) has log K −3.07 at 25 °C,
Thermoddem's (Blanc & Piantone 2006) −7.30, and McDowell et al. (1977) measured −7.37 to −7.56 on this
basis. The spread across the databases is −3.1 to −7.3, and coreclath's entry may represent poorly
crystalline apatite (unconfirmed). The crystalline apatite used here is a much stronger P sink
(~10⁴ times less soluble).

*Validation* (25 °C, coreclath basis, primary values converted with coreclath's pK₃ 12.3218 and pK_w
13.9951). A phase whose database value is more than 0.5 log units from the primary value is not used
until resolved. The primary value replaces it only for a clear sign or basis slip, and a phase without
a primary source is left out.

| Phase | Database (converted) | Primary | Source | Difference | Used |
|---|---|---|---|---|---|
| Hydroxyapatite | −7.30 (Thermoddem Natur) | −7.56 (measured) / −7.37 (best value) | McDowell et al. 1977, *J. Res. NBS* 81A | +0.26 / +0.07 | yes |
| Whitlockite (β-Ca₃(PO₄)₂) | −4.30 (Thermoddem "high") | −4.28 | Gregory et al. 1974, *J. Res. NBS* 78A | −0.03 | yes |
| Whitlockite (low) | −6.03 (Thermoddem "low") | −4.28 | Gregory et al. 1974 | −1.75 | no |
| Brushite | −6.61 | −6.59 | Gregory et al. 1970, *J. Res. NBS* 74A | −0.02 | yes |
| OCP, Ca₄H(PO₄)₃ | −11.52 (Thermoddem ":3H2O", NIST46.4) | −11.44 (pK_sp 48.4, 23.5 °C) | Tung et al. 1988, *J. Res. NBS* 93 | −0.08 | yes |
| OCP (":2.5H2O", 84nan) | −9.82 | −11.44 | Tung et al. 1988 | +1.61 | no |
| Bobbierite | −0.50 | −0.56 | Taylor et al. 1963, *Trans. Faraday Soc.* 59 (also Thermoddem's source) | +0.06 | yes |
| Newberyite | −5.80 | −5.82 | Taylor et al. 1963 | +0.03 | yes |
| Vivianite | −17.70 | −11.12 | Al-Borno & Tomson 1994, *GCA* 58 | −6.57, Thermoddem's sign is flipped (+3.272 gives −11.15) | yes, the primary value |
| Struvite | (none) | −10.28 (pK_sp 13.36), the fit used gives −10.47 | Bhuiyan et al. 2007, *Environ. Technol.* 28; Ohlinger et al. 1998: −10.18 | fit vs Ohlinger −0.29 | yes |
| Monetite | −6.91 | not read (McDowell et al. 1971, *Inorg. Chem.* 10, paywalled) | | | no (unverified) |
| Fluorapatite | −22.55 | −22.64 from McCann 1968 pK 59.6, seen only in a search snippet | | | no (unverified) |
| Farringtonite, Na₃PO₄, Na₂HPO₄, NaH₂PO₄ | | no primary solubility product found | | | no |

*Temperatures.* Hydroxyapatite, whitlockite and brushite have Thermoddem's analytic expressions
(0–300 °C), vivianite Al-Borno & Tomson's fit (5–90 °C), and struvite Bhuiyan et al.'s fit (10–60 °C,
pK_sp 13.57 at 25 °C against their tabulated 13.36, diverging outside the range to log K −21 at 150
°C). OCP, bobbierite and newberyite have 25 °C values only (no enthalpy), reliable to about 50 °C,
which covers the skin. Each phase is listed only where its data apply
(`PhreeqcIntegrator$phosphate_minerals()`): vivianite to 90 °C, struvite to 60 °C, and OCP, bobbierite
and newberyite to 50 °C. The ocean (1.45 °C) is a few degrees below vivianite's and struvite's measured
range. **A hot porewater sees fewer phosphates.** The 160 °C bound case (`cases/sweep/T160.R`) lists
only hydroxyapatite, whitlockite and brushite, so the 25 °C-only phases (OCP, bobbierite, newberyite)
are excluded there, and so are vivianite and struvite. The molar volumes
(`data/mineral_molar_volumes.csv`) are Core11's for hydroxyapatite and from unit cells (webmineral.com)
for the others. **OCP and whitlockite are PROVISIONAL.** OCP's comes from the Mathew, Brown, Schroeder
& Dickens (1988) cell, read through search results only, and whitlockite's is approximate, from natural
whitlockite scaled by P, not from β-Ca₃(PO₄)₂.

### Refill

Hydration consumes porewater (after the first orbit the porewater is 20–37% short of the pores).
Between chunks the porewater is brought back to the pore volume φ(t) implies, m_pw(0)·φ(t)/φ0, by
moving ocean water into it. m_pw(0) is the first chunk's porewater, so the target is on the run's own
normalised basis and is 0 at t = 0. `tidal_kinetic_restart(refill_kg =)` mixes that much of the ocean
solution into the porewater (or porewater into the ocean, if negative) before its first orbit. MIX
moves whole solutions, so every element is conserved between the two. The gas phases are untouched, and
the MIX steps run before the SELECTED_OUTPUT is defined and add no rows (such a restart has no `i_soln`
rows). `tests/restart_closure_check.R` checks that a restart with a refill of +0.02 or −0.01 kg ends
its orbit with the same whole-system element totals (porewater, ocean, rock, pore and ocean gas) as
without it, and that the water moved. Each chunk reports the refill (`refill_in_kg`, and per orbit of
the previous chunk) next to the tidal exchange (`exchange_kg_per_orbit`).

The refill dwarfs the tidal exchange. After orbit 1 it moves 0.033–0.044 kg of ocean water into the
porewater (normalised basis), against a tidal exchange of ~1e-8 kg per orbit (layered strain at K_d,
1.9e-6 kg under the Liao strain), ~10⁶ times more. After that it falls to 1e-4–1e-5 kg per orbit,
still ~10⁴× the tidal exchange, and is sometimes negative (a little porewater returned) without IOM. In
this model, then, **hydration-driven inflow, not tidal pumping, sets how much ocean water enters the
rock** while the primaries hydrate.

### The ocean per kg of rock

The PHREEQC basis is per kg of the reacting box's rock, so the ocean box per kg of rock must be the
ocean's mass ÷ the box's rock mass. The ocean's mass is the profile's `Mocean_kg` (2.131e19 kg). The
config's geometric ocean at 1000 kg/m³ (2.113e19 kg) would give 0.9914 of the correct value at every k
(k = 1e-14: 1002.2 against 1010.9, 1e-11: 31.80 against 32.08, 1e-7: 0.4619 against 0.4659), and the
Simulator's circulation time still uses it (0.9% lighter). Each restart scales the parent's ocean
water, ocean minerals and ocean gas by old box rock ÷ new box rock (`tidal_kinetic_restart(ocean_scale
=)`, 1.07–1.38 after orbit 1, ~1.00–1.07 later). Without this scaling a restart would pass on the
parent's ocean per kg of the *parent's* box, so when δ shrank (69 → 58 km) the real ocean would
silently shrink with it.

### Drift

Each chunk's whole-system element totals (porewater, ocean, rock, IOM, both gases) are checked against
its start (t = 0, or the parent's end with the ocean scaled), cumulatively (`drift_max`,
`drift_element`), with a warning above 1e-4 relative. The production tolerance is unchanged. In the
24-orbit runs below the drift is ≤ 7e-9 at 160 °C (the solver's tolerance, Al and C), and 4–6e-6 in H
at the skin's temperature. These runs were made before restarts read the solutions back from the dump
(§3, Restarts), and there the skin's temperature changes between chunks, so each restart redefined the
porewater at the new temperature from its pH, which moved its H by that much. With the dump the
per-chunk drift in H is ≤ 2.1e-9. `tests/porosity_evolution_check.R` requires < 1e-6 at 25 °C
(measured 3.9e-10).

### How much ocean water the refill takes

`refill_summary(chunks)` scales each chunk's refill to the whole core by the reacting box's rock mass
(the skin is a shell over the whole core) and compares it with the profile's ocean mass (`Mocean_kg` =
2.13e19 kg). In the 24-orbit runs below (with the ocean scaling):

| run | after orbit 1 | cumulative, 24 orbits | fraction of the ocean | ocean concentration factor (scaled / model's own ocean) |
|---|---|---|---|---|
| skin T, 15 wt% IOM | 1.42e18 kg | 1.49e18 kg | 6.97% | 1.0749 / 1.0749 |
| skin T, no IOM | 1.41e18 kg | 1.42e18 kg | 6.65% | 1.0713 / 1.0713 |
| 160 °C, 15 wt% IOM | 1.66e18 kg | 1.82e18 kg | 8.53% | 1.0933 / 1.0933 |
| 160 °C, no IOM | 1.25e18 kg | 1.27e18 kg | 5.95% | 1.0633 / 1.0633 |

Some 6–8% of the ocean is drawn into the core in the first orbit, as the 69 km skin's metal corrodes
and its primaries hydrate (both nearly instantaneous here, since instantaneous corrosion is accepted as
an idealisation of early alteration). After that the refill adds only ~0.1–0.5% of the ocean over 23
orbits, and is sometimes slightly negative without IOM. The two concentration factors agree to all four
digits in every chunk. One comes from the cumulative refill, M_ocean/(M_ocean − refilled), and the other
from the model's own ocean water, M_ocean/(ocean per kg of rock × box rock). (Without the ocean scaling
they would be 1.06–1.09 and 1.08–1.12.) The ocean per kg of rock at a chunk's end is then the geometric
one × 1/CF, the water the refill has drawn. The concentration factor follows from the water loss alone,
and the solutes the refilled water brings into the rock are in the chemistry.

### 24-orbit runs

These runs use production k = 1e-7 at φ0 = 0.32 via Kozeny–Carman with d = 1.59 cm, the skin box, φc =
0.4, ΔP from the strain and the refill on. The last row, for comparison, uses φc 0.5, a fixed ΔP of 1e5
and no refill:

| | φ at orbit 24 | k | δ | skin T | circulation | refill, orbit 1 |
|---|---|---|---|---|---|---|
| skin temperature (`"profile"`), 15 wt% IOM | 0.279 | 6.0e-8 | 58 km | 31 → 27 °C | 1,121 yr | 0.033 kg |
| skin temperature, no IOM | 0.259 | 4.5e-8 | 52 km | 31 → 26 °C | 1,260 yr | 0.034 kg |
| 160 °C override, 15 wt% IOM | 0.222 | 2.6e-8 | 43 km | 160 °C | 1,626 yr | 0.044 kg |
| 160 °C override, no IOM | 0.169 | 9.9e-9 | 31 km | 160 °C | 2,651 yr | 0.038 kg |
| *160 °C, φc 0.5, ΔP 1e5, no refill, IOM / none* | *0.219 / 0.165* | *2.5e-8 / 9.1e-9* | *45 / 32 km* | *160 °C* | *2,014 / 3,334 yr* | |

(k, δ and circulation are those of the last chunk, which starts a little above the final φ.) The 160 °C
rows override the temperature (`porewater_temp = 160`) on the skin box, whose own temperature at k =
1e-7 is 31 °C (falling to 26–27 °C as the skin thins). At the skin's temperature the core closes much
less (φ 0.32 → 0.28 with IOM, → 0.26 without), because hydration is slower. At 160 °C it closes as with
the comparison settings (→ 0.22, 0.17), so the current defaults change little (the strain binds
throughout, so ΔP does not matter). In every case the core closes and does not clog in 24 orbits,
and the strain limit binds. The first chunk takes 77 s at the skin's temperature and 190–515 s at 160
°C.

### Permeability and surface area from grain size

`Simulator$new(k = "kozeny_carman", grain_diameter = d)` sets the permeability from the grain size and
the porosity by the Kozeny–Carman relation for granular, unconsolidated media, **k =
d²φ³/(180(1 − φ)²)** (Carman 1937). The 180 is 36 × Carman's constant 5 for packed spheres, with the
specific surface 6/d, checked against secondary sources (the 1937 paper itself not seen). A number for
`k` also works, so earlier runs can be reproduced. `kozeny_carman_k()` and `kozeny_carman_d()` (the
inverse) compute it. The table (m²):

| φ | d = 1 µm | 30 µm | 1 mm | 2 cm |
|---|---|---|---|---|
| 0.2 | 6.9e-17 | 6.2e-14 | 6.9e-11 | 2.8e-8 |
| 0.3 | 3.1e-16 | 2.8e-13 | 3.1e-10 | 1.2e-7 |

At the production φ = 0.32, Choblet et al.'s preferred 1e-14–1e-13 m² needs grains of 5–16 µm, and
**the production k = 1e-7 m² implies grains of 1.6 cm, a centimetre-scale gravel**, not a
chondritic matrix. Fractures can dominate a real rock's permeability, in which case Kozeny–Carman is a
lower bound. Alteration reduces φ and k falls as φ³, so with `run_tidal_evolution()` a Kozeny–Carman
template follows φ(t) through k automatically.

The surface area follows the same grain size. `ssa = "geometric"` sets each mineral's specific surface
to 6/(ρᵢ·dᵢ) × `roughness`, where ρᵢ = molar mass/Vm from `data/mineral_molar_volumes.csv` and
`roughness`, default 1, is the BET-over-geometric factor. `grain_diameter_by_mineral` (a named vector,
m) gives per-mineral sizes, e.g. a fine matrix for feldspars and sulfides and coarser olivine and
pyroxene, and `PhreeqcIntegrator$apply_ssa()` writes the areas into every generated input. `ssa =
"legacy"` (the default) uses the fixed 0.1 and 5 m²/g (UNSOURCED). The difference is large. At d = 1.6
cm (production k) the geometric area of forsterite is 1.2e-4 m²/g, 800× below the legacy 0.1, and at 10
µm it is 0.19 m²/g. `tests/kozeny_carman_check.R` checks the table, the Simulator and the SSA lines.

### One sweep over (φ, d, strain) (design only, not run)

This sweep replaces the separate k, surface-area and stress sweeps. It covers φ 0.1–0.4 and d 1 µm–1
cm (setting both k and the surface areas). Since ΔP is derived from the strain, it also covers the
tidal scenario and core viscosity in place of ΔP: `tidal_scenario = "dissipative_core"` with
`core_viscosity` 4.4e11, 1e12, 2.3e12 Pa·s and `"dissipative_shell"` with 3e13, 1e16, 1e22 (the exchange
at k = 1e-14 is 6.8e-8, 5.4e-8, 4.8e-8 and 5.4, 5.1, 4.8e-8 per half-orbit). Finally it covers the
troilite and IOM shear moduli at their bracket ends (`grain_elastic_moduli(mu_case =)`), which move μ_d
but not K_d. The table gives k (m²) over the grid, marked C (Choblet et al.'s preferred 1e-14–1e-13), F
(Fisher et al.'s aquifers 1e-13–1e-10) and L (Liao et al.'s 1e-12–1e-7):

| φ | 1 µm | 10 µm | 100 µm | 1 mm | 1 cm |
|---|---|---|---|---|---|
| 0.1 | 6.9e-18 | 6.9e-16 | 6.9e-14 C | 6.9e-12 F L | 6.9e-10 L |
| 0.2 | 6.9e-17 | 6.9e-15 | 6.9e-13 F | 6.9e-11 F L | 6.9e-9 L |
| 0.3 | 3.1e-16 | 3.1e-14 C | 3.1e-12 F L | 3.1e-10 L | 3.1e-8 L |
| 0.4 | 9.9e-16 | 9.9e-14 C | 9.9e-12 F L | 9.9e-10 L | 9.9e-8 L |

### Assumptions stated, not swept

- Ingested water has the bulk-ocean composition. The water drawn in is not the fluid just
  expelled (bottom currents disperse that), but it is not fresh either, because the ocean box has the
  products of every earlier cycle. MIX with the evolving ocean represents this.
- The porewater is one well-mixed box. Real exchange penetrates about the displacement distance
  per tide (centimetres in Liao et al.'s case), the pressure change reaches only the skin depth δ
  (~20 m to ~70 km over k = 1e-14 to 1e-7 m²), and both spread by dispersion. Over Myr runs this is
  an approximation.

### Fe/Ni metal: instantaneous corrosion, an accepted idealisation

The metal is a `dissolve_only` equilibrium phase, so it corrodes in the first step and ~2 mol H₂ per kg
of rock arrives at t = 0 (1.998 mol of Fe metal in the corrected rock, and the Ni corrodes at 160 °C but
not at 25 °C). This is not a rate claim. It stands for the alteration and outgassing of early
differentiation, and a pure-water ocean over fresh rock is likewise an idealisation, chosen so that the
early outgassing totals can be computed cleanly. **The science is the evolution over millions of years,
and the first orbits are a transient**, so results are to be judged on the long runs
(`docs/PROGRESS.md`).

### The H₂ budget (by electron balance)

`h2_budget(csv)` gives per orbit, on the run's basis (mol per kg of the box's rock), the H₂ generated
since t = 0 by electron balance:
- H2_from_Fe = ½ the change of Σ n·z_Fe over every Fe-bearing phase in the porewater's minerals and
  kinetic reactants and the ocean's minerals, plus the dissolved Fe as Fe(II). z_Fe comes from each
  formula by charge balance (`fe_valence()`: metal and schreibersite 0, troilite and the Fe(II)
  silicates 2, magnetite 8/3, cronstedtite 2.5, Fe(III) phases 3). The CSV has no Fe(3) column, and
  `dissolved_Fe_mol` bounds that error;
- H2_from_Ni (metal to Ni(II)) and H2_from_P (schreibersite's P to phosphate, 2.5 per P);
- `H2_generated`, their sum.

It also gives where the H₂ is at the end of the orbit (dissolved in the ocean, in the ocean gas,
dissolved in the porewater, in the pore gas, with `H2_present` their sum) and `H2_other_redox` =
generated − present, the H₂ the carbon, sulfur and nitrogen couples took up or gave. The lower bounds
(metal and schreibersite to Fe(II): `H2_from_metal`, `H2_from_schreibersite`) are also reported.
`scale_to_core(budget, sim)` multiplies by the rock mass of the reacting box, the whole permeable layer
and the whole core (V (ρ_core − φ ρ_water) from the Simulator's geometry), which gives the inventory
available for early loss to space. For the 160 °C example case (§4, k = 1e-7, skin box 69 km, 12
orbits, `to_ocean`):

| per kg of rock | orbit 1 | orbit 12 | whole core, orbit 12 |
|---|---|---|---|
| H₂ generated (electron balance) | 2.349 mol | 2.462 mol | 1.54e20 mol |
| of which by Fe (to Fe(II), magnetite 8/3, Fe(III)) | 2.347 mol | 2.444 mol | 1.52e20 mol |
| of which by schreibersite's P | 0.0016 mol | 0.018 mol | 1.1e18 mol |
| H₂ in the ocean gas | 2.31 mol | 2.43 mol | 1.51e20 mol |
| H₂ dissolved in the ocean | 0.030 mol | 0.030 mol | 1.9e18 mol |
| H₂ dissolved in the porewater | 0.0094 mol | 0.0071 mol | 4.4e17 mol |
| H₂ present (sum) | 2.349 mol | 2.462 mol | 1.54e20 mol |
| generated − present | −2.0e-5 mol | −2.4e-4 mol | |
| lower bound (metal + schreibersite to Fe(II)) | 2.002 mol | 2.037 mol | 1.27e20 mol |

(The whole core's rock is 6.24e19 kg, and here the layer is the whole core. The numbers were computed
without the phosphate phases.) The electron balance closes the budget. The H₂ generated equals the H₂
present to 1e-4 relative, whereas the lower bound falls 17% short, because of the iron oxidised beyond
Fe(II) (0.42 mol of magnetite per kg of rock by orbit 12, with the dissolved Fe at 4e-12 mol). Ni metal
does not corrode in this run (H2_from_Ni = 0), and the other redox couples gave 2.4e-4 mol. So ~1.5e20
mol of H₂ (3e17 kg) is in the ocean gas within the first orbits, the inventory available for early
loss.

### Schreibersite

Schreibersite is kinetic, with a placeholder bulk rate from one 1-day measurement
(RECONSTRUCTION_STATUS §4, with the corrosion chemistry after Pasek & Lauretta 2005). Its molar volume
comes from mindat unit-cell data.

### Gas

Each reservoir has a fixed-pressure GAS_PHASE (Parkhurst & Appelo 2013, GAS_PHASE), so a bubble forms
only when the partial pressures exceed the pressure. All gases are Peng–Robinson (Peng & Robinson
1976), with the decoupled gases using their analogues' constants (phreeqc.dat). `pore_gas_fate`
brackets the pore bubble's fate: `"to_ocean"` (default) or `"retain"` (an upper bound on retention when
`Vgas_Vpore` > 1).

`"to_ocean"` is the default because under `"retain"` the pores would contain about 2.2 pore volumes
of gas, which is physically impossible, so `"retain"` is an upper bound on retention. Carbon
attribution is the same in both fates (organic − blank ocean C agrees to 3 significant figures at 150
and 160 °C). `"to_ocean"` without plume loss saturates the ocean in H₂, so it is an upper bound on
ocean H₂. A third fate, `"retain_to_pore_volume"`, keeps gas up to the pore volume and sends the
excess to the ocean (§7, Gas).

### Pressure

Pressure comes from the profile and acts only through gases. Under the LLNL aqueous model, aqueous
molar volumes are inactive, and the minerals' were stripped (`databases/README.md`), so mineral
equilibria are at the database's reference pressure (~0.1–0.2 log units at 260 bar).

### Initial ocean

The ocean is pure water by default, so all ocean carbon is attributable to the organics. The cometary
option follows Miller et al. (2026) and Rubin et al. (2019). It starts the ocean with ~0.087 molal DIC,
so attribution must use the paired difference.

### The cometary-ice ocean, extended (an option, never the default)

`cometary_ice_ocean(v = 0.02, f = 1, include_NH3 = TRUE, include_less_volatile = FALSE, include_O2 =
TRUE)`.
- NH₃ is included by default. Rubin et al. (2019) Table 2 gives 0.67 ± 0.20 per 100 H₂O, and
  Miller et al. (2026) apply v_NH3 = 0.02 to "NH₃ and other compounds that are more volatile", so the
  depleted ice has 0.0134 per 100 H₂O (7.44e-3 molal). Its absence from their Table 5 is treated as an
  omission. Rubin et al. note a factor of 10 for NH₃ between the 2015 bulk value and the earlier 3.1 au
  measurements (0.06 and 0.15 over the two hemispheres, which themselves differ 2.5×), left unexplained.
- f scales the whole cometary load after the depletion. It is the fraction of the accreted cometary
  volatiles retained in the ocean.
- Less volatile ices. Miller et al.'s rule keeps ices less volatile than NH₃ at 273 K in full, and
  their Table 5 includes none. `data/comet67p_volatility.csv` lists every species of Rubin et al.'s
  Table 2 (bulk column, with its ± errors), classified against NH₃ at 273 K. The classification uses
  the sublimation curves of Fray & Schmitt (2009) where they cover the species (their Table 5
  polynomials evaluated at 273 K, read from the PDF, with NH₃ at 14 bar there), and otherwise the normal
  boiling point against NH₃'s 239.8 K, a proxy marked as such (standard values, not checked against a
  handbook here). Fray & Schmitt make CO, O₂, N₂, CH₄, C₂H₆, CO₂, H₂S and the noble gases more volatile,
  and SO₂, HCN, CH₃OH, HC₃N, C₆H₆ and HCOOH less. Radicals and atoms (SO, S, S₂, H₂CS, PO) are left
  unclassified (they have no boiling point, and Rubin et al. tie S, S₂ and SO partly to grains). OCS,
  C₃H₈, HCl and HBr come out more volatile by the proxy (OCS and C₃H₈ are then missing from Table 5
  too). **The 21 less-volatile CHONS species add C 1.08, N 0.177 and S 0.177 per 100 H₂O, against Table
  5's C 0.174, N 0.00356, S 0.022: 6.2×, 50× and 8×** (HF, also less volatile, is not CHONS and is not
  added). They enter by their elements, as CO₂, NH₃ and H₂S with H and O corrections, so that no element
  name can match a phase (the CO-as-cobalt trap). Under `c4_only` their carbon can only become C(+4),
  and the electrons go to H₂ or reduce sulfate, and HCN's N becomes N(−3). **Caveat:** in a cold ocean
  these organics are kinetically metastable, so converting them to CO₂ plus reductant is an equilibrium
  assumption. Under `"coupled"` the equilibrium decides as well. No species is kept as itself
  (methanol, formate and acetate exist in the database but equilibrate), so in an equilibrium model
  these species cannot "stay as themselves".

The starting ocean (v = 0.02, NH₃ and O₂ included, 1.45 °C, 66.7 atm, mol/kgw, with clathrate SIs from
the initial-ocean step, which has no clathrate phases to precipitate):

| less volatile | mode | f | pH | DIC (C+4) | H₂ | sulfate | sulfide | Mtg / reduced C | gas at 66.7 atm | clathrate SI (CH₄/Mtg, CO₂, H₂S) |
|---|---|---|---|---|---|---|---|---|---|---|
| no | c4_only | 1 | 2.08 | 0.0866 | 6e-14 | 9.4e-3 | 2.8e-3 | Mtg 0.0102 | none | −0.81, −1.07, −1.88 |
| no | c4_only | 0.1 | 2.98 | 8.7e-3 | 1e-13 | 9.4e-4 | 2.8e-4 | Mtg 1.0e-3 | none | −1.81, −2.07, −2.88 |
| no | coupled | 1 | 5.43 | 0.0960 | 1e-10 | 3e-9 | 0.0122 | CH₄ etc. 2.8e-6 | none | −1.92, −1.06, −1.25 |
| no | coupled | 0.1 | 5.45 | 9.6e-3 | 1e-10 | 3e-10 | 1.2e-3 | 2.9e-8 | none | −2.91, −2.06, −2.25 |
| yes | c4_only | 1 | 5.79 | 0.557 | 0.055 | 0 | 0.099 | Mtg 1.3e-3 | 0.19 L (H₂ 0.48, CO₂ 0.064, H₂S 6e-3, Mtg 4e-3 mol) | −1.73, −0.37, −0.38 |
| yes | c4_only | 0.1 | 5.76 | 0.068 | 0.064 | 0 | 0.011 | Mtg 8e-4 | 7e-3 L (mostly H₂) | −1.91, −1.25, −1.31 |
| yes | coupled | 1 | 5.96 | 0.402 | 4e-10 | 1e-8 | 0.108 | 0.033 (CH₄) | 0.016 L (CH₄ 0.051 mol) | +0.30, −0.55, −0.34 |
| yes | coupled | 0.1 | 6.01 | 0.042 | 5e-10 | 5e-10 | 0.011 | 1.0e-3 | none | −0.39, −1.51, −1.31 |

With the less volatile ices the organics' electrons consume all the O₂ (no sulfate is left), and the
ocean is saturated in H₂ under `c4_only` (~0.06 molal at 66.7 atm, the rest a bubble). Under `coupled`
the carbon reduces to methane instead, oversaturating CH₄ hydrate at f = 1. The less volatile organics
multiply the ocean's DIC by 6–7 (0.087 → 0.56 molal at f = 1).

### Initial fluids

Initial fluids are defined at 25 °C, pe 4 (no oxidant or reductant) and heated by
REACTION_TEMPERATURE. Defined at the run temperature, pe 4 would put O₂ into the porewater (4.9 molal
O(0) at 250 °C).

### Trace N and S seeds: tested, not added

The default ocean is pure water plus the C seed. Seeding N or S as well (1e-14 mol/kgw in every
initial solution, alone or together) changed nothing that matters in the example case (§4) at 25 and
160 °C (2 orbits, 15 wt% IOM). Every variant ran with 0 errors with and without them. At 160 °C the
solver's retries were 33 in every variant, and at 25 °C they moved without pattern (4 with C alone, 18
with C + N, 12 with C + S, 4 with all three). Runtimes differed only with the machine's load. With the
C seed and the ocean step's `-diagonal_scale`, the default already converges, so `trace_N_molal` and
`trace_S_molal` were not added.

### Trace C seed

The C seed of 1e-14 mol/kgw is numerical, so the carbon mass balance can be solved when the IOM has
released almost nothing.

### Databases

`coreclath_CH4coupled.dat` is Mohit's Core11 (Neveu 2017, extended by Diab 2021 and Melwani Daswani
2025) plus CH₄/CO₂/H₂S clathrates. The uncoupled and C4only files derive from it
(`databases/README.md`).

## 7. Interpreting outputs

### Rock density

PHREEQC reports no solid volumes, so rock density is computed in R from a run's selected output:

```r
source("phreeqc_tools.R")
d <- bulk_density("phreeqc_inputs/<file>.csv", porosity = sim$params$porosity)
d$summary     # per REACT step (and t = 0 from the .pqi): rho_solid, rho_fluid, gas_saturation, rho_bulk
d$breakdown   # per phase: mass and volume, and their fractions
```

- `rho_solid` = Σ nᵢMᵢ / Σ nᵢVᵢ over every mineral with nonzero moles (equilibrium phases, kinetic
  primaries, schreibersite, Fe/Ni metal, clathrates) and the remaining IOM. Molar volumes come from
  `data/mineral_molar_volumes.csv` (sources in `data/README.md`), at 25 °C and 1 bar. Leaving out
  thermal expansion and compression costs about 1% for silicates at 160 °C and 26 MPa, less than the
  porosity and IOM-density uncertainty.
- Remaining IOM mass = initial IOM − Σ (m0ᵢ − nᵢ) × M(formulaᵢ) over the IOM channels (each
  channel's `-formula` is what the solid loses). Its density (default `iom_density =
  "delarue2026"`) is Delarue et al. (2026) Eq. 10, ρ = 1472 − 345 (H/C) + 94 / ((H/C) + 0.122)
  kg/m³, evaluated on each row at the residual IOM's atomic H/C. That H/C is the bulk IOM's C and H
  (`iom_bulk_formula()`, H/C 0.753) less what the channels have released (Mtg counted as CH₄).
  At t = 0 the density is 1.320 g/cm³. The fit is to kerogens and coals and is within 100 kg/m³ for
  H/C below 1 (Delarue et al. Fig. 5). `iom_density = 1.4` (any number) fixes it instead.
  Compression is left out, because their bulk modulus K0 = 10.2 GPa gives 0.25% at the core's 26 MPa.
- `rho_bulk` = (1 − φ) ρ_solid + φ ρ_pore, with φ the model's porosity (0.32) and
  ρ_pore = (1 − s_g) ρ_fluid + s_g ρ_gas, s_g = min(Vgas/Vpore, 1) the pore gas saturation. ρ_fluid
  is PHREEQC's solution density (`Rho_soln`). Under the LLNL aqueous model PHREEQC does not compute
  water's density from T and P, so it is ~1 g/cm³ for water at any temperature (give
  `fluid_density =` to override).
- Gas beyond the pore volume is not counted (s_g is capped at 1, and the gas is not rock volume
  either). Rows where Vgas > Vpore are flagged in the `gas_exceeds_pores` column (with `Vgas_Vpore`),
  and `bulk_density()` warns when there are any, because their ρ_bulk is for pores full of gas, not a
  physical state. In the 160 °C example case that is the first REACT row under `"to_ocean"` (2.14
  pore volumes, before the transfer) and every row under `"retain"` (up to 2.26).

For the 160 °C example case (15 wt% IOM, φ = 0.32), ρ_solid is 2.920 g/cm³ and ρ_bulk 2.305 g/cm³ at
t = 0. At orbit 12, ρ_solid is 2.759 and ρ_bulk 2.255 (`"to_ocean"`) or 1.880 (`"retain"`, the bubble
filling the pores). With a fixed IOM density of 1.4 g/cm³ these would be 2.977 and 2.344 at t = 0 and
2.806 and 2.287 or 1.912 at orbit 12, so the Delarue density lowers ρ_solid by 1.9% at t = 0 and 1.7%
at orbit 12, and ρ_bulk by 1.4–1.7%. The residual IOM's H/C hardly moves in 12 orbits (0.75298 to
0.75298, with 0.01 g of 150 g released), so its density stays 1.3197 g/cm³. Lizardite (26% of the
solid volume by orbit 12) takes up water. The profile's own silicate-layer density is 2353 kg/m³
(`rhoSil_kgm3`).

A published core density is a bulk density (porosity and pore fluid included), so compare it with
`rho_bulk`. Bagheri et al. (2025, Table 5) fit Enceladus's libration, shape, heat flux, gravity and mass
and find a core density of 2160–2400 kg/m³ and a core radius of 183–213 km (their priors were 2000–2800
kg/m³ and 120–220 km), in both their core- and shell-dissipation cases. Against that range:
- the model's ρ_bulk is 2.305 g/cm³ at t = 0 and 2.255 at orbit 12 under `"to_ocean"`, inside it;
- 1.880 under `"retain"` (pores full of gas) is below it, one more sign that `"retain"` is a bound,
  not a state;
- the profile's own silicate-layer density, 2353 kg/m³ (`rhoSil_kgm3`, the config's
  `core_density`), and its core radius, 194.2 km, are inside it.

Iess et al. (2014) measured the gravity field and inferred a moment of inertia of ~0.335 MR². Their
"relatively low core density of ~2.4 g/cm³" is an example structure compatible with that moment
of inertia (with a 60 km, 1 g/cm³ H₂O layer), not a measured core density. The profile is fitted to
C/MR² = 0.335 ± 0.001. `EnceladusConfig`'s radius, mass, core radius and core density are the
profile's (R_m, M_kg, RsilMean_m, rhoSil_kgm3), not values from either paper.

### Gas

`Vgas_Vpore` > 1 on porewater rows means the bubble is larger than the pore space, and under
`"retain"` that is an upper bound, not a physical state. `pore_gas_fate = "retain_to_pore_volume"`
(for `tidal_kinetic` and its restart, not `tidal_equil`) keeps the pore gas up to the pore volume and
sends the excess to the ocean. The porewater step `PUT`s each component's excess, gas × max(0, 1 −
Vpore/Vgas). A short step after it (not punched in the run's CSV) removes that excess from the
porewater and its gas phase (KINETICS 4), and the ocean step adds it to the ocean (KINETICS 3, as
`to_ocean`). The pore gas left after each removal is written to `<base>_poregas.csv`, which the mass
balance (`tests/gas_mass_balance_check.R`) and restarts read. The default is `"to_ocean"`.

Worked example 1 under the three fates (25 °C, 2 orbits, end of orbit 2, the layered strain at K_d,
computed without the phosphate phases, so the `to_ocean` porewater pH is 11.04 here against 10.27 in
§4):

| fate | ocean pH | ocean H₂ (molal) | ocean gas | porewater pH | pore gas |
|---|---|---|---|---|---|
| `to_ocean` | 7.41 | 0.065 (saturated) | 0.90 L | 11.04 | none |
| `retain_to_pore_volume` | 7.62 | 0.065 (saturated) | 0.58 L | 14.06 | 0.66 pore volumes |
| `retain` | 7.62 | 4.2e-9 | none | 14.06 | 2.47 pore volumes |

(Under the Liao strain, with ~100× the exchange, the ocean pH is 7.57, 9.55, 9.55 and `retain`'s ocean
H₂ 5.7e-7.)

Under `"to_ocean"` the whole pore bubble reaches the ocean every orbit, whatever the transport rate, so
the ocean's redox is set by the gas, not by the circulation. Under `"retain_to_pore_volume"` the excess
alone (the metal's H₂ beyond 1 pore volume) still saturates the ocean in H₂. Only under `"retain"` does
the ocean see the H₂ through the tidal exchange alone. The ocean gas phase (MIX 4 rows, `g_` columns)
contains H₂ and organic Mtg under `"to_ocean"`, and `organic_contribution()` reports it.

### Steady state

The water the porewater takes in has the evolving bulk-ocean composition, so each orbit's products
feed the next. The expected end state is therefore a **steady state, an ocean composition that stops
changing, but not necessarily zero net exchange**. The IOM degrades kinetically, with a long tail of
high-Ea pools, so its release may never quite stop. The porewater (100–160 °C) and the ocean (1.45 °C)
also cannot share one saturated mineral assemblage, so circulation can sustain a constant flux (e.g.
silica dissolved in the pores and precipitated in the ocean). True equilibrium with no net exchange
needs the IOM exhausted, no kinetic primaries left and equal temperatures.

`steady_state_diagnostics(csv)` shows the approach per orbit. It gives the net elemental flux
into the ocean (the change of the ocean side, water + minerals + gas, over the orbit, since the ocean
exchanges only with the porewater) and the net change of the secondary minerals in the porewater and in
the ocean (mol of formula units per kg of the box's rock). In the 160 °C example case (12 orbits), the
porewater's secondary minerals grow by 1.09 mol in orbit 1, 0.26 in orbit 2, 0.014 in orbit 6 and
0.0038 in orbit 12 (the hydration transient dying away), and the ocean's by ~3e-10 per orbit from orbit
6. The Na flux into the ocean settles at 1.75e-8 per orbit (2.1e-8 in orbit 1), Si falls from 8.7e-10
to 1.6e-10, and C is ~1-2e-11 per orbit. These numbers use the layered strain at K_d and were computed
without the phosphate phases. Under the Liao strain, with ~140× the exchange, they are Na 2.4e-6, Si
1.2e-7 to 2.1e-8 and C rising from 2.2e-10 to 2.5e-9. Twelve orbits are far from the steady state,
which is what the long runs are for.

### Carbon

Use `carbon_by_oxidation_state()` (`phreeqc_tools.R`), not PHREEQC's valence totals, to see which
carbon is reduced, because PHREEQC files formate, acetate and other organics under C(4).

Glein & Waite (2020) adopt Postberg et al. (2009)'s salt envelopes for Enceladus's ocean carbonate,
0.2 molal Cl with 0.1 molal total dissolved carbonate (TDC = HCO₃⁻ + CO₃²⁻, without CO₂(aq)) or 0.05
molal Cl with 0.01 molal TDC, and the pH range 8.5–9, which gives log a(CO₂) between −4.6 and −3.2
(their preferred estimate). The model's oceans compare as follows:
- the cometary-ice option starts with 0.087 molal DIC, but at pH 1.85 (2.08 with NH₃, the default) it
  is almost all CO₂(aq) (0.0866 molal, log a(CO₂) ≈ −1.1). TDC is 2e-6 molal and there is no Cl, and
  with the less volatile ices DIC is 0.56 molal at pH 5.8;
- the pure-water default, after 12 orbits of the 160 °C example case, has 2.9e-10 molal DIC at pH
  8.51 (layered strain at K_d, and 3.4e-8 at pH 10.63 under the Liao strain). The long runs
  (`docs/PROGRESS.md`) are for finding whether the organic contribution can approach the observed
  carbonate over the core's lifetime.

### Phosphate and calcium

The apatite is Ca-limited and no calcite forms, so the apatite–calcite buffer is inactive in this rock.
In the reference under `"apatite"` (pure-water ocean, to 1e4 yr), hydroxyapatite takes all the P
(0.0133 mol per kg of box material from ~100 yr on, with ocean P ~1.6e-6 molal) and nearly all the Ca
the rock releases. The box material has 0.286 mol Ca per kg (0.85 kg of silicate, 15 wt% IOM: diopside
0.233, anorthite 0.053 mol), but only 0.0755 mol is released by 1e4 yr. Hydroxyapatite contains 0.0667
of it and mesolite 0.0088, dissolved Ca stays ≤ 1e-7 molal, and no carbonate forms at any time. The
table shows why so little is released (by chunk, porewater SIs from SELECTED_OUTPUT):

| t (yr) | porewater pH | Diopside left | SI | rate (mol/s) | Anorthite left | SI |
|---|---|---|---|---|---|---|
| 0 | | 0.2333 | | | 0.0528 | |
| 0.26 | 11.4 | 0.2330 | −3.4 | 4.3e-11 | 6e-6 | −10.2 |
| 5.4 | 12.6 | 0.2262 | −2.2 | 4.1e-11 | 0 | −10.6 |
| 16 | 13.3 | 0.2148 | −0.40 | 2.3e-11 | 0 | −9.2 |
| 40 | 13.6 | 0.2106 | +0.33 | 0 | 0 | −8.9 |
| 1e4 | 13.9 | 0.2106 | +0.77 | 0 | 0 | −8.9 |

- Anorthite is not limiting, since it all dissolves within ~1 yr (0.053 mol).
- Diopside's slow kinetics, then the SR > 1 guard, stop the release. Diopside dissolves at ~4e-11
  mol/s (PK04, SSA 0.1 m²/g) and gives up 0.023 mol in ~30 yr. Then the porewater (pH 13.5–14) becomes
  supersaturated in diopside (SI > 0 in 85 of 134 chunks, up to +1.5), and the rate law's guard sets
  its rate to 0.
- The low release is not an accounting artefact. Each chunk's −m0 equals the previous chunk's end
  moles, for both minerals, and at most 5.5e-14 mol of anorthite is dropped where it falls below
  PHREEQC's precision.

At the same time (1e4 yr), Na is 0.051 molal under `"apatite"` against 0.063 without phosphate phases.

### What tidal-only circulation can reach

The exchange in ENTICES is tidal only (buoyancy-driven flow is out of scope, and Choblet et al. 2017
and others are for comparison). The tidal pressure reaches only the skin depth δ, so the rock that
reacts and exchanges is the top L_eff = min(L, δ) of the layer (§6). At the defaults (φ = 0.32, φc =
0.4, α = 0.8, the layered strain at K_d = 9.8 GPa, `"dissipative_core"`, placeholders µ = 17 GPa, η =
1e12 Pa·s, ΔP from the strain plus ocean loading):

| k (m²) | L_eff | skin T | exchange per half-orbit | circulation time | (at K = 100 GPa) | (under `"liao"`) |
|---|---|---|---|---|---|---|
| 1e-14 | 22 m | 1.4 °C | 5.4e-8 | 2.2e8 yr | 1.9e-8, 6.2e8 yr | 8.3e-6, 1.4e6 yr |
| 1e-12 | 219 m | 1.4 °C | 5.4e-8 | 2.2e7 yr | 1.9e-8, 6.2e7 yr | 8.3e-6, 1.4e5 yr |
| 1e-11 | 692 m | 1.4 °C | 5.4e-8 | 7.1e6 yr | 1.9e-8, 2.0e7 yr | 8.3e-6, 4.6e4 yr |
| 1e-7 | 69 km | 31.0 °C | 4.5e-8 | 1.2e5 yr | 1.9e-8, 2.9e5 yr | 6.1e-6, 908 yr |

**Tidal pumping alone reaches only cold rock just below the seafloor unless k is at the top of the
range.** Deep, hot rock is reached only as δ approaches L. The strain limit binds throughout, so the
exchange per half-orbit hardly depends on k. k sets how deep (and so how warm and how much rock) the
reacting skin is, and with it the circulation time. **With the coupled body's strain, only k near the
top of the range circulates the ocean through the core within Enceladus's age** (1.2e5 yr at 1e-7, 7e6
yr at 1e-11, 2.2e8 yr at 1e-14). Over η = 4.4e11–2.3e12 the times are 0.75–1.15× these.

## 8. Known limitations and open decisions

Open decisions:
- **The core viscosity η is the main one.** It selects the layered strain within the tidal scenario
  and so the exchange. The reference uses 1e12 Pa·s (dissipative core), and the sweep brackets it at
  4.4e11 and 2.3e12 (Bagheri et al. 2025 Table 5) and tries the dissipative-shell scenario at 1e14. The
  core rigidity µ (17 GPa, the geometric mean of 4–70 GPa) is a placeholder too.
- Surface areas: the legacy per-mineral values are UNSOURCED, and the `"geometric"` option ties them to
  the grain size.
- φc (0.4) is a placeholder from sorted grain packs.
- The k_floor values (1e-14, 1e-12 m²) are a modelling bracket, not from a measurement (§6).
- The spin-up and the IOM's sulfur provenance are open (`docs/PROGRESS.md` Part 10).
- **The OCP and whitlockite molar volumes are PROVISIONAL** (§6, Phosphate phases).

Settled choices: the metal corrodes in the first step (an accepted idealisation), the pore gas goes to
the ocean (`"retain"` and `"retain_to_pore_volume"` are the bracket), the carbon mode is `c4_only`, and
the skin reacts at its own profile temperature (a fixed temperature is a bound). For the phosphates,
`"apatite"` is the reference, the Ca-limited apatite stands (diopside's kinetics and the SR > 1 guard,
§7), and the phosphate temperature limits and the other volumes are accepted.

Limitations: the 25 °C brine exceeds the B-dot range. The pure-water ocean is unbuffered, so its pH is
not an Enceladus pH. A run cannot reach 4.5 Gyr at the validated step limits, so production runs end at
a steady state or the 1e8 yr backstop (§3). The k ≈ 1e-14 case, whose circulation time (~5e8 yr)
exceeds the backstop, is not at steady state by design.

## 9. Repository map

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
| `docs/` | Design decisions, test write-ups, progress (`docs/PROGRESS.md`), `reviews/`, `archive/`. |
| `legacy/` | Retired files (`legacy/README.md`). **Do not use for runs.** |

## 10. References

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
