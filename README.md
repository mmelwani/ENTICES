# ENTICES

Enceladus Tidally-Induced Circulation and Exchange Simulator.

## 1. What ENTICES is

ENTICES models water–rock reaction in Enceladus' porous core. An R model of tidally driven
porewater circulation (a two-layer gravity field, tidal stress over each 1.375-day orbit, Darcy
flow through the porous core) sets how much porewater exchanges with the ocean per orbit. It
writes PHREEQC inputs that react a reconstructed primordial rock — kinetic silicates and sulfide,
Fe–Ni metal, schreibersite and insoluble organic matter (IOM) — with porewater at the core's
temperature and pressure, and exchange that porewater with the ocean every orbit. PHREEQC runs
the inputs; ENTICES' R tools read the results.

The science question is whether IOM degradation in the core adds carbon to the ocean — the
carbonate Cassini saw in the plume — on a timescale that constrains how long hydrothermal
circulation has lasted. ENTICES therefore starts, by default, from pure water and attributes ocean
carbon to the organics by difference against an otherwise identical run without them
(`pair_organic_run()`).

## 2. Requirements

- **R** (tested with 4.5.2) with the packages the code loads: `tidyverse`, `conflicted`, `dplyr`,
  `patchwork`, `glue`, `R6` (`ENTICES_LF_20260917.R`); `ggplot2` for the rate plots in
  `mineral_dissolution_constants.r` and the worked examples.
- **PHREEQC 3.8.6.** `phreeqc_tools.R` looks for it in the `PHREEQC_EXE` environment variable, then
  at `C:/Program Files/USGS/phreeqc-3.8.6-17100-x64/bin/Release/phreeqc.exe`, then on the `PATH`.
  Every generated `.pqi` is stamped `# Target PHREEQC version: 3.8.6`, and `run_phreeqc()` warns
  if the executable reports another version: minor phases and trace species differ between
  versions (Lucas's 3.5.0 outputs rerun on 3.8.6 moved dissolved Mn 130× and sulfide ~11×; bulk
  chemistry was stable). Treat the warning as "not reproducible as stamped".
- **Python 3** only for the database scripts (`scripts/make_c4only.py`, `strip_mineral_vm.py`,
  `gas_eos.py`, `export_mineral_volumes.py`).
- **Databases** (`databases/`, see `databases/README.md`), chosen by `carbon_mode` and stamped in
  each `.pqi` header:

  | `carbon_mode` | Database |
  |---|---|
  | `"c4_only"` (default) | `coreclath_C4only.dat` |
  | `"coupled"` | `coreclath_CH4coupled.dat` |
  | `"ch4_uncoupled"` | `coreclath_CH4uncoupled.dat` |
  | `"temperature"` | uncoupled below `ch4_cutoff_C` (150 °C), coupled above |

## 3. How to run (regenerating production)

Run everything from the repo root.

1. **Pull `main`.**
2. **Set paths.** PHREEQC as above (`Sys.setenv(PHREEQC_EXE = "...")` if needed).
3. **Run the self-tests first.** Each prints `SUMMARY: 0 failure(s)` when passing (`tests/README.md`
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

4. **Choose the profile.** `ph$profile_file` (default
   `profiles/EnceladusProfile_Seawater_10.0ppt_Tb272.4578K_PorousRock.txt`, a PlanetProfile model,
   `profiles/README.md`) sets the porewater pressure (layer mean, 261.2 atm), the ocean pressure
   (ocean floor, 66.7 atm) and, with `porewater_temp = "profile"`, the porewater temperature. There
   is no default pressure: a missing profile is an error.
5. **Set parameters** (§5): the circulation (`k`, `layer_thickness`), the organics
   (`calculate_phreeqc_parameters(organic_wt_percent = 15)`), the temperature and the run length.
6. **Generate the `.pqi`:**

   ```r
   source("ENTICES_LF_20260917.R")
   sim <- Simulator$new(k = 1e-7, layer_thickness = 194200, config = EnceladusConfig$new())
   sim$run_simulation()
   ph <- PhreeqcIntegrator$new(sim)
   ph$calculate_phreeqc_parameters(organic_wt_percent = 15)
   f <- ph$generate_phreeqc_input_tidal_kinetic(output_dir = "phreeqc_inputs", mode = "ungrouped",
                                                n_orbital_cycles = 12, porewater_temp = 160, substeps = 1,
                                                trace_C_molal = 1e-14)
   ```

7. **Run** with the database stamped in its header (`# Database:`):

   ```r
   db <- sub("^# Database: ", "", grep("^# Database: ", readLines(f, n = 20), value = TRUE))
   r  <- run_phreeqc(f, db, output = sub("\\.pqi$", ".out", f), workdir = dirname(f))
   ```

   or `phreeqc <file>.pqi <file>.out <repo>/databases/coreclath_C4only.dat` from the input's
   directory.
8. **Outputs** (beside the `.pqi`): `<base>.csv`, the SELECTED_OUTPUT (tab-separated); `<base>.out`;
   `<base>_initial_ocean.csv` for a cometary ocean. In the CSV:
   - `state`, `soln`, `step`: `i_soln` rows are the initial solutions; then four `react` rows per
     orbit of a tidal run — REACT (soln 1: porewater + rock), MIX 4 (soln 4: ocean + exchanged
     porewater + seafloor phases + ocean gas), MIX 1 (soln 1: resampled porewater), MIX 2 (soln 2:
     the ocean at the end of the orbit);
   - `pH`, `pe`, `Alk`, `mu`, `mass_H2O` (kg);
   - element totals (molal) for all 36 database elements; `H_total_mol`, `O_total_mol` (moles,
     water included; the `H` and `O` totals leave water out);
   - `m_` / `la_` species molalities and log activities;
   - equilibrium-phase moles by name, `d_` their change; `k_` kinetic reactants (moles), `dk_` their
     change; `<mineral>_remain_g`, `IOM_*_mol`; `Fe_metal` (the `Fe` column is dissolved Fe; `Ni`
     appears twice, dissolved then the metal);
   - `si_` primary-mineral SIs: use these, not the `.out` tables (below);
   - gas: `pressure`, `total mol`, `volume`, `g_<gas>` moles of the gas phase in use, and `Vgas_L`,
     `P_gas_atm`, `Vgas_Vpore` (gas volume / pore volume);
   - `Rho_soln` (PHREEQC's solution density).
9. **Paired runs:** the organics' contribution by difference against a same-rock blank:

   ```r
   source("organic_contribution.R")
   pr <- pair_organic_run(ph, "tidal_kinetic", mode = "ungrouped", n_orbital_cycles = 12, porewater_temp = 160,
                          trace_C_molal = 1e-14, output_dir = "phreeqc_inputs/pair160")
   pr$contribution  # DIC, CH4, NH4+, H2, total C (mol/kgw), ocean gas (mol): run, blank, difference
   ```

10. **Sweeps:** `parameter_sweep.R` (generates over a grid of `k` and layer thickness).
11. **Restarts** continue a run from its CSV (`prev_csv`), with the parent's `.pqi` beside it:
    `generate_phreeqc_input_tidal_kinetic_restart()`, `generate_phreeqc_input_bulk_kinetic_restart()`.
    They carry every element at 15 digits, the equilibrium phases, kinetic reactants, IOM channels
    and both gas phases, and refuse a parent run in another carbon mode or at other pressures
    unless told (`carbon_mode`, `allow_pressure_change`).
12. **Runtimes** (this machine, PHREEQC 3.8.6): Lucas's production case (12 orbits, 15 wt% IOM)
    320–515 s at 150–160 °C; 2 orbits at 160 °C 140–290 s; equilibrium generators seconds. The
    first orbit (the metal pulse into the gas phase) is the slowest.
13. **Common failures:**
    - *No convergence, or a step running for tens of minutes:* check the input has the KNOBS
      `-convergence_tolerance 1e-8` after SELECTED_OUTPUT; at 25 °C, or with a cometary ocean, give
      `trace_C_molal = 1e-14`.
    - *The version warning:* rerun on PHREEQC 3.8.6.
    - *A cometary ocean with a carbon-free porewater:* `tidal_kinetic` stalled and `tidal_equil` failed
      ("Cl has not converged"); give `trace_C_molal = 1e-14` (both warn). Known failure:
      `tidal_equil` + `pore_gas_fate = "to_ocean"` + cometary ocean + `c4_only` crashes PHREEQC 3.8.6
      (segmentation fault in the ocean step, all temperatures); `tidal_kinetic` runs that combination.
    - *SIs in the `.out` file:* the printed SI tables are offset by n(H₂O) × Δlog a_w (0.143 per H₂O
      at 25 °C in the brine); use the CSV's `si_` columns.
    - *B-dot activities above ~1 molal:* the 25 °C NaOH porewater brine reaches I ≈ 2.4 for the
      first ~3×10⁴ orbits; treat its pH and activities as qualitative.
    - *Disk space:* `KNOBS -logfile true` (every template) writes `phreeqc.log` with every solver
      iteration, ~300 MB per 12-orbit run; set it to `false` if you don't need it.

## 4. Worked examples

`ENTICES_runs.Rmd` has each as code with narrative and expected results:

1. **A 25 °C smoke test** (minutes): `tidal_kinetic`, 2 orbits, 15 wt% IOM, `trace_C_molal = 1e-14`.
2. **Lucas's production case:** k = 1e-7 m², full-core layer, 12 orbits, 160 °C, 15 wt% IOM. Expect
   (the default `pore_gas_fate = "to_ocean"`, orbit 12): ocean pH 11.3, ocean H₂ 6.5e-2 molal
   (saturated), ocean DIC 1.3e-7 molal (the seed is 1e-14), 0.82 L of ocean gas; the pore bubble
   0.34 L after the first orbit, 3e-4 L by orbit 12.
3. **A 100 °C organic vs blank pair** with `pair_organic_run()`.
4. **The cometary-ice ocean:** `initial_ocean = cometary_ice_ocean()`; the ocean starts at pH 1.85
   with 0.087 molal DIC.
5. **The pore-gas bracket:** the production case with `pore_gas_fate = "retain"` and `"to_ocean"`:
   ocean H₂ 1.5e-5 vs 6.5e-2 molal at orbit 12, the same organic DIC.

## 5. Parameters

Status: **sourced** (a reference), placeholder (a stand-in with a stated basis), or **UNSOURCED**.

| Name | Default | Units | What it does | Justification / source | Status | Sensitivity planned |
|---|---|---|---|---|---|---|
| `k` | none (production 1e-7) | m² | permeability; sets Darcy flow and exchange | production choice | — | yes (sweeps) |
| `layer_thickness` | none (production 194200) | m | permeable layer; 194.2 km is the whole core | production choice | — | yes |
| `max_tidal_stress` (`EnceladusConfig`) | 1e5 | Pa | tidal stress amplitude driving the flow | none | **UNSOURCED** | yes |
| `porosity` (`Simulator`) | 0.32 | — | pore volume, porewater mass | profile `phiRockMax`; Styczinski et al. 2023 | sourced | — |
| `profile_file` | the Enceladus profile | — | porewater and ocean P, porewater T | PlanetProfile (Styczinski et al. 2023) | sourced | other profiles |
| ocean temperature | 274.6 K (1.45 °C) | K | ocean SOLUTION and seafloor steps | profile, top of the rock layer (274.566 K) | sourced | — |
| `porewater_temp` | `"ocean"` | °C or `"profile"` | porewater reaction temperature | production uses 150 and 160 °C | — | yes |
| `organic_wt_percent` | 0 (production 15) | wt% | IOM in the rock | deliberate, to be swept | — | yes |
| `carbon_mode` | `"c4_only"` | — | carbon redox treatment and database | §6 | decided (PI, 10-27) | `"coupled"` as the bound |
| `pore_gas_fate` | `"to_ocean"` | — | the pore bubble's fate | §6; PROGRESS Part 7 | decided (PI, 10-29); `"retain"` is the other bound | yes (the bracket) |
| `initial_ocean` | `"pure_water"` | — | ocean start; `cometary_ice_ocean()` option | attribution; Miller et al. 2026 | sourced | the cometary option |
| `cometary_ice_ocean(v)` | 0.02 | — | volatile depletion of 67P ice | Miller et al. 2026 (from Titan's N budget) | sourced, not for Enceladus | v |
| `include_NH3` | FALSE | — | NH₃ in the cometary ocean | Table 5 lists none | **pending the PI** | — |
| `trace_C_molal` | NULL (production 1e-14) | mol/kgw | numerical C seed | smallest value tested to converge | numerical | — |
| `mode`, `grouping_factor`, `substeps` | `"ocean_circulations"`, —, 1 | — | orbits per step | grouped mode compounds exchange: x_step = 1 − (1 − x_tide)^G | — | — |
| `rec_eqb_minerals` | FALSE | — | promote fast primaries to equilibrium (bulk_kinetic) | a modelling choice | — | — |
| mineral SSA | 0.1 (fo, en, di, teph); 5 (ab, an, Kfs, tro) | m²/g | reactive surface area | none (D7) | **UNSOURCED** | yes |
| Fe, Ni metal | `dissolve_only` equilibrium | — | corrodes completely in step 1 | none (D6) | **UNSOURCED** | yes (D6) |
| schreibersite rate | 1.370e-7 | s⁻¹ | bulk first-order | one 1-day lab measurement (RECONSTRUCTION_STATUS §4, §8) | placeholder | — |
| IOM kinetics | A = 2e15 s⁻¹, distributed Ea | — | IOM channels | Burnham 2019; Miller et al. 2025 | sourced; `Ea_N` a one-point fit; `IOM_S` uncalibrated | Ea width, CHn m0 |
| IOM density | 1.4 | g/cm³ | `bulk_density()` only | none | **UNSOURCED** | — |
| `ch4_cutoff_C` | 150 | °C | only for `carbon_mode = "temperature"` | the old default | — | — |

## 6. Model choices and justifications

**Carbon mode (`c4_only`).** Coupled carbon is C(+4) only; reduced carbon enters only as decoupled
methane (`Mtg`) from the IOM, with IOM_CHn as `Mtg 1 H -2.2` (a 2.2 e⁻ reduction imposed by hand).
In `"coupled"` runs equilibrium reduces carbonate with H₂, so they are an **upper bound on abiotic
methanation**, not the expected case (PROGRESS Part 3).

**IOM composition.** Miller et al. (2026) SI Table 9, from Alexander et al. (2007) CR-chondrite IOM:
a CR composition with CM (Murchison) kinetics. The sulfur provenance is open (D10).

**IOM kinetics.** Distributed-Ea pools (the Vitrimat approach of Burnham 2019), A = 2e15 s⁻¹,
calibrated to the heating experiments of Miller et al. (2025). `Ea_N` = 242.78 kJ/mol is a one-point
fit and its m0 0.875 a floor; `IOM_S` has no kinetic calibration. The 350 °C residue check passes
for H; O is under-released by 13–30% (PROGRESS Part 4).

**Rock.** The Scenario R reconstruction (`reconstruction/RECONSTRUCTION_STATUS.md`) on the CI
composition of Lodders et al. (2025). Lawrencite and the 15 wt% organics are deliberate; the
organics will be swept.

**Mineral rates.** Palandri & Kharaka (2004) for the silicates, feldspars and troilite (as hexagonal
pyrrhotite). Tephroite: log k(H⁺) = −5.4 at pH 2, n = 0.47 (Pincus et al. 2026, Table S1, citing
Casey et al. 1993), with an Ea that falls with pH (a fit to Casey et al.'s Fig. 4) and no neutral
branch: the rate at pH > 7 is held at its pH-7 value.

**Surface areas: UNSOURCED.** Forsterite, enstatite, diopside and tephroite are at 0.1 m²/g; albite,
anorthite, K-feldspar and troilite at 5 m²/g, so the feldspars react against 50× the area of the
mafic silicates (D7).

**Fe/Ni metal: UNSOURCED instantaneous corrosion** (`dissolve_only` equilibrium phases), so the
metal's ~2.6 mol H₂ arrives at t = 0 (D6). Every gas result is dominated by this pulse until D6 is
decided.

**Schreibersite:** kinetic, with a placeholder bulk rate from one 1-day measurement
(RECONSTRUCTION_STATUS §4; corrosion chemistry after Pasek & Lauretta 2005); its molar volume from
mindat unit-cell data.

**Gas.** A fixed-pressure GAS_PHASE in each reservoir (Parkhurst & Appelo 2013, GAS_PHASE): a bubble
forms only when the partial pressures exceed the pressure. All gases Peng–Robinson (Peng & Robinson
1976); decoupled gases with their analogues' constants (phreeqc.dat). `pore_gas_fate` brackets the
pore bubble's fate: `"to_ocean"` (default) or `"retain"` (an upper bound on retention when
`Vgas_Vpore` > 1).

Why `"to_ocean"` is the default (PI, 2026-10-29): under `"retain"` the pores hold about 2.2 pore
volumes of gas, which is physically impossible, so `"retain"` is an upper bound on retention.
Carbon attribution is the same in both fates (organic − blank ocean C agrees to 3 significant
figures at 150 and 160 °C). `"to_ocean"` without plume loss saturates the ocean in H₂, so it is an
upper bound on ocean H₂. A third fate, `"retain_to_pore_volume"` (keep gas up to the pore volume,
send the excess to the ocean), is planned once the metal corrosion rate (D6) is decided, since the
bubble is metal H₂.

**Pressure** from the profile. Pressure acts only through gases: under the LLNL aqueous model
aqueous molar volumes are inactive, and the minerals' were stripped (`databases/README.md`), so
mineral equilibria are at the database's reference pressure (~0.1–0.2 log units at 260 bar).

**Initial ocean:** pure water by default, so all ocean carbon is attributable to the organics. The
cometary option follows Miller et al. (2026) and Rubin et al. (2019); it starts the ocean with
~0.087 molal DIC, so attribution must use the paired difference.

**Initial fluids** are defined at 25 °C, pe 4 (no oxidant or reductant) and heated by
REACTION_TEMPERATURE; defined at the run temperature, pe 4 put O₂ into the porewater (4.9 molal
O(0) at 250 °C).

**Trace C seed** 1e-14 mol/kgw: numerical, so the carbon mass balance can be solved when the IOM has
released almost nothing.

**Databases:** `coreclath_CH4coupled.dat` is Mohit's Core11 (Neveu 2017, extended by Diab 2021 and
Melwani Daswani 2025) plus CH₄/CO₂/H₂S clathrates; the uncoupled and C4only files derive from it
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
  channel's `-formula` is what the solid loses). Its density, `IOM_DENSITY_G_CM3` = 1.4 g/cm³, is
  **UNSOURCED** [citation needed: a measured density of chondritic IOM or kerogen].
- `rho_bulk` = (1 − φ) ρ_solid + φ ρ_pore, with φ the model's porosity (0.32) and
  ρ_pore = (1 − s_g) ρ_fluid + s_g ρ_gas, s_g = min(Vgas/Vpore, 1) the pore gas saturation. ρ_fluid
  is PHREEQC's solution density (`Rho_soln`): under the LLNL aqueous model PHREEQC does not compute
  water's density from T and P, so it is ~1 g/cm³ for water at any temperature (give
  `fluid_density =` to override).

For Lucas's 160 °C case (15 wt% IOM, φ = 0.32): ρ_solid 2.977 g/cm³ and ρ_bulk 2.344 g/cm³ at t = 0;
at orbit 12, ρ_solid 2.806 and ρ_bulk 1.91 (`"retain"`, the bubble filling the pores) or 2.29
(`"to_ocean"`). Lizardite (26% of the solid volume by orbit 12) takes up water. The profile's own
silicate-layer density is 2353 kg/m³ (`rhoSil_kgm3`).

**Comparing with Enceladus.** A published core density is a bulk density (porosity and pore fluid
included), so compare it with `rho_bulk`. The usual reference is Iess et al. (2014) [citation
needed: the core density and its uncertainty, to be taken from the paper; the PDF was not available].

### Gas

`Vgas_Vpore` > 1 on porewater rows means the bubble is larger than the pore space: under `"retain"`
that is an upper bound, not a physical state. The ocean gas phase (MIX 4 rows, `g_` columns) holds
H₂ and organic Mtg under `"to_ocean"`; `organic_contribution()` reports it.

### Carbon

Use `carbon_by_oxidation_state()` (`phreeqc_tools.R`), not PHREEQC's valence totals, to see which
carbon is reduced: PHREEQC files formate, acetate and other organics under C(4).

## 8. Known limitations and open decisions

`docs/PROGRESS.md` Part 10 lists the decisions (D6 metal corrosion rate, D7 surface areas, D9 spin-up,
D10 IOM sulfur, the default `pore_gas_fate`, NH₃ in the cometary ocean) and the next actions. In
short: the H₂ pulse is the metal (D6); surface areas are unsourced (D7); the 25 °C brine exceeds the
B-dot range; the pure-water ocean is unbuffered, so its pH is not an Enceladus pH.

## 9. Repository map

| Path | What it is |
|---|---|
| `ENTICES_LF_20260917.R` | **Entry point.** Circulation model and all PHREEQC input generators (`generate_phreeqc_input_*`). |
| `ENTICES_runs.Rmd` | The worked examples (§4). |
| `iom_module.R` | IOM decomposition channels: calibrated sub-pools. |
| `mineral_dissolution_constants.r` | Mineral rate laws in R (used by `rec_eqb_minerals`) and rate plots. |
| `phreeqc_tools.R` | Runs PHREEQC; carbon re-binning; database name tools; `bulk_density()`; `system_totals()`. |
| `organic_contribution.R` | `pair_organic_run()`, `organic_contribution()`. |
| `parameter_sweep.R` | Parameter sweep helper. |
| `CHANGELOG.md` | Change log. |
| `databases/` | PHREEQC databases (`databases/README.md`). |
| `profiles/` | The PlanetProfile interior profile (`profiles/README.md`). |
| `data/` | `mineral_molar_volumes.csv` (`data/README.md`). |
| `reconstruction/` | The primordial CI-chondrite rock reconstruction. |
| `tests/` | Self-tests and PHREEQC test inputs (`tests/README.md`). |
| `examples/` | Lucas's example production run and his run notebook; the redox budget diagnostic inputs. |
| `scripts/` | Database scripts and standalone analyses. |
| `docs/` | Design decisions, test write-ups, progress (`docs/PROGRESS.md`); `reviews/`, `archive/`. |
| `legacy/` | Retired files (`legacy/README.md`). **Do not use for runs.** |

## 10. References

- Alexander, C.M.O'D., Fogel, M., Yabuta, H., Cody, G.D. (2007). The origin and evolution of
  chondrites recorded in the elemental and isotopic compositions of their macromolecular organic
  matter. *Geochim. Cosmochim. Acta* 71, 4380–4403. doi:10.1016/j.gca.2007.06.052
- Alexander, C.M.O'D., et al. (2022). Sulfur abundances and isotopic compositions in bulk
  carbonaceous chondrites and insoluble organic matter. *Meteorit. Planet. Sci.* 57, 334–351.
  doi:10.1111/maps.13746
- Burnham, A.K. (2019). Kinetic models of vitrinite, kerogen, and bitumen reflectance. *Org.
  Geochem.* 131, 50–59. doi:10.1016/j.orggeochem.2019.03.007
- Casey, W.H., et al. (1993). *Geochim. Cosmochim. Acta* 57, 785–793 [citation needed: authors and
  title; cited via Pincus et al. 2026].
- Fray, N., Schmitt, B. (2009). Sublimation of ices of astrophysical interest: a bibliographic
  review. *Planet. Space Sci.* 57, 2053–2080.
- Glein, C.R., Waite, J.H. (2020). *Geophys. Res. Lett.* 47, e2019GL085885 [citation needed: the
  title and the DIC estimate, from the paper].
- Green, E.C.R., Holland, T.J.B., Powell, R., Weller, O.M., Riel, N. (2025). *J. Petrol.*,
  doi:10.1093/petrology/egae079 (its `tc-ds636` dataset cross-checks the molar volumes).
- Iess, L., et al. (2014) [citation needed: the Enceladus gravity paper, for the core density].
- King, A.J., et al. (2015). Modal mineralogy of CI and CI-like carbonaceous chondrites by X-ray
  diffraction. *Geochim. Cosmochim. Acta* 165, 148–160.
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
- Pincus et al. (2026). *ACS Earth Space Chem.* 10, 1174–1184 [citation needed: authors and title,
  as cited in the tephroite rate law].
- Rubin, M., et al. (2019). Elemental and molecular abundances in comet 67P/Churyumov-Gerasimenko.
  *Mon. Not. R. Astron. Soc.* 489, 594–607. doi:10.1093/mnras/stz2086
- Sloan, E.D., Koh, C.A. (2007). *Clathrate Hydrates of Natural Gases*, 3rd ed. CRC Press. Table
  2.5a, p. 78.
- Styczinski, M.J., Vance, S.D., Melwani Daswani, M. (2023). PlanetProfile: Self-consistent interior
  structure modeling for ocean worlds and rocky dwarf planets in Python. *Earth Space Sci.* 10,
  e2022EA002748. doi:10.1029/2022EA002748
- White, A.F., Brantley, S.L. (2003) [citation needed: the lab–field weathering-rate discrepancy
  paper].
