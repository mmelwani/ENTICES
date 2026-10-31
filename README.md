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
    - *25 °C ocean step (fixed 10-29):* at 25 °C the ocean step (MIX 4) failed to converge
      ("Numerical method failed on all combinations"; Al, C, Ca, Cl, Fe, K, Mg, Mn, Na, Ni, P, S and
      Si "not converged" at ~1e-15–1e-8 mol) for scattered porewater fractions: of 1e-8 to 1e-3,
      3e-8, 2e-6 and 5e-6 failed and their neighbours ran (the pre-`exchange-physics` code failed at
      3.3e-6 too; `retain` fails the same way, so it is not the gas transfer). Every ocean step now
      runs with KNOBS `-diagonal_scale true`, switched off after it (`OCEAN_STEP_KNOBS_ON/OFF`): all
      fractions tested converge, and a converging 160 °C orbit is unchanged (1e-8 relative). Set for
      the whole run instead, it shifted the porewater step's m_H₂S by 0.9%. Trace seeds were not
      needed.
    - *A cometary ocean with a carbon-free porewater:* `tidal_kinetic` is ~4× slower (and stalls
      outright with a CO₂-only ocean), and `tidal_equil` fails ("Cl has not converged", the Cl being
      the porewater's lawrencite); give `trace_C_molal = 1e-14` (both warn). Known failure:
      `tidal_equil` + `pore_gas_fate = "to_ocean"` + cometary ocean + `c4_only` crashes PHREEQC 3.8.6
      (segmentation fault in the ocean step, all temperatures); `tidal_kinetic` runs that combination.
      A 20-line reproducer is in `tests/known_failures/`. Details: `docs/MATRIX_2026-10-28.md`.
    - *SIs in the `.out` file:* the printed SI tables are offset by n(H₂O) × Δlog a_w (0.143 per H₂O
      at 25 °C in the brine); use the CSV's `si_` columns.
    - *B-dot activities above ~1 molal:* the 25 °C NaOH porewater brine reaches I ≈ 2.4 for the
      first ~3×10⁴ orbits; treat its pH and activities as qualitative.
    - *Disk space:* `KNOBS -logfile` is off by default (PI, 10-29). Every generator takes
      `logfile = TRUE` to switch it on for convergence debugging: `phreeqc.log` then holds every
      solver iteration, 258–346 MB per 12-orbit production run at 150–160 °C (~300 MB), and 11.9 GB
      from one stalled run.

## 4. Worked examples

`ENTICES_runs.Rmd` has each as code with narrative and expected results:

1. **A 25 °C smoke test** (minutes): `tidal_kinetic`, 2 orbits, 15 wt% IOM, `trace_C_molal = 1e-14`.
2. **Lucas's production case:** k = 1e-7 m², full-core layer, 12 orbits, 160 °C, 15 wt% IOM. Expect
   (the default `pore_gas_fate = "to_ocean"`, orbit 12): ocean pH 11.3, ocean H₂ 6.5e-2 molal
   (saturated), ocean DIC 1.3e-7 molal (the seed is 1e-14), 0.82 L of ocean gas; the pore bubble
   0.34 L after the first orbit, 3e-4 L by orbit 12.
3. **A 100 °C organic vs blank pair** with `pair_organic_run()`.
4. **The cometary-ice ocean:** `initial_ocean = cometary_ice_ocean()` (NH₃ included by default since
   10-31); the ocean starts at pH 2.08 with 0.087 molal DIC and 7.4e-3 molal NH₄⁺.
5. **The pore-gas bracket:** the production case with `pore_gas_fate = "retain"` and `"to_ocean"`:
   ocean H₂ 1.5e-5 vs 6.5e-2 molal at orbit 12, the same organic DIC.

## 5. Parameters

Status: **sourced** (a reference), placeholder (a stand-in with a stated basis), or **UNSOURCED**.

| Name | Default | Units | What it does | Justification / source | Status | Sensitivity planned |
|---|---|---|---|---|---|---|
| `k` | none (production 1e-7) | m² | permeability; sets Darcy flow and exchange | literature ranges (§6, Flow driving) | **production value at the top of Liao et al.'s range** | yes (D12 sweep) |
| `layer_thickness` | none (production 194200) | m | permeable layer; 194.2 km is the whole core | production choice | — | yes |
| `core_rigidity`, `core_viscosity` (`EnceladusConfig`) | 17e9, 1e12 | Pa, Pa·s | the tidal strain (exchange) | Bagheri et al. 2025 Table 5, core dissipation | **placeholders** at the geometric means of 4–70 GPa and 4.4e11–2.3e12 | yes (range ends) |
| `reactive_depth` (`Simulator`) | `"skin"` | — | the reacting box: the top min(L, δ), or `"layer"` | PI 2026-10-29 | decided | — |
| `strain_scale` (`EnceladusConfig`) | 1 | — | factor on the core's tidal strain ε₀ | bracket 0.046/0.29 = 0.16 to 1 (PI 10-31) | **largest uncertainty in the exchange** | yes (0.16–1) |
| `phi0` (`Simulator`; was `porosity`) | 0.32 | — | initial porosity: k, α, storage, δ, strain volume, pore water | PlanetProfile `phiRockMax` (one fit) | **one fit**; end-state check against Bagheri's core density | yes (0.1–0.4) |
| `phi_c` (`Simulator`) | 0.4 | — | critical porosity: α = min(1, φ/φc); φ ≥ φc is a suspension (α = 1) | Nur et al. 1998; Dvorkin & Nur 2000 (36–40% sorted grains, sandstones) | PI 10-30; still a placeholder | yes (0.4–0.6) |
| `alpha` (`Simulator`) | NULL (from φ/φc) | — | Biot coefficient override (must be ≥ φ) | — | — | — |
| `eccentricity` (`EnceladusConfig`) | 0.0047 | — | tidal strain | Enceladus's orbit | — | — |
| `max_tidal_stress` (`EnceladusConfig`) | NULL: ΔP = α·ε_v·M | Pa | tidal pore-pressure amplitude driving the Darcy flow | derived from the strain (PI 10-30); a number overrides it (the old default was an UNSOURCED 1e5) | derived; override documented | through ε_v, η, φ |
| `porosity` (`Simulator`) | 0.32 | — | pore volume, porewater mass | profile `phiRockMax`; Styczinski et al. 2023 | sourced | — |
| `profile_file` | the Enceladus profile | — | porewater and ocean P, porewater T | PlanetProfile (Styczinski et al. 2023) | sourced | other profiles |
| ocean temperature | 274.6 K (1.45 °C) | K | ocean SOLUTION and seafloor steps | profile, top of the rock layer (274.566 K) | sourced | — |
| `porewater_temp` | `"ocean"` | °C or `"profile"` | porewater reaction temperature | production uses 150 and 160 °C | — | yes |
| `organic_wt_percent` | 0 (production 15) | wt% | IOM in the rock | deliberate, to be swept | — | yes |
| `carbon_mode` | `"c4_only"` | — | carbon redox treatment and database | §6 | decided (PI, 10-27) | `"coupled"` as the bound |
| `pore_gas_fate` | `"to_ocean"` | — | the pore bubble's fate | §6; PROGRESS Part 7 | decided (PI, 10-29); `"retain"` is the other bound | yes (the bracket) |
| `initial_ocean` | `"pure_water"` | — | ocean start; `cometary_ice_ocean()` option | attribution; Miller et al. 2026 | sourced | the cometary option |
| `cometary_ice_ocean(v)` | 0.02 | — | volatile depletion of 67P ice | Miller et al. 2026 (from Titan's N budget) | sourced, not for Enceladus | v |
| `include_NH3` | TRUE | — | NH₃ in the cometary ocean (0.67 per 100 H₂O × v) | Rubin et al. 2019 Table 2; Miller et al.'s text applies v to it | decided (PI 10-29): Table 5's omission | — |
| `cometary_ice_ocean(f)` | 1 | — | fraction of the accreted cometary volatiles retained in the ocean, applied after v | PI 10-29 | — | f |
| `include_less_volatile` | FALSE | — | ices less volatile than NH₃ at 273 K, undepleted (× f) | Rubin et al. 2019 Table 2; Fray & Schmitt 2009 or boiling point; `data/comet67p_volatility.csv` | option | yes |
| `trace_C_molal` | NULL (production 1e-14) | mol/kgw | numerical C seed | smallest value tested to converge | numerical | — |
| `mode`, `grouping_factor`, `substeps` | `"ocean_circulations"`, —, 1 | — | orbits per step | grouped mode compounds exchange: x_step = 1 − (1 − x_tide)^G | — | — |
| `rec_eqb_minerals` | FALSE | — | promote fast primaries to equilibrium (bulk_kinetic) | a modelling choice | — | — |
| mineral SSA | 0.1 (fo, en, di, teph); 5 (ab, an, Kfs, tro) | m²/g | reactive surface area | none (D7) | **UNSOURCED** | yes |
| Fe, Ni metal | `dissolve_only` equilibrium | — | corrodes completely in step 1 | none (D6) | **UNSOURCED** | yes (D6) |
| schreibersite rate | 1.370e-7 | s⁻¹ | bulk first-order | one 1-day lab measurement (RECONSTRUCTION_STATUS §4, §8) | placeholder | — |
| IOM kinetics | A = 2e15 s⁻¹, distributed Ea | — | IOM channels | Burnham 2019; Miller et al. 2025 | sourced; `Ea_N` a one-point fit; `IOM_S` uncalibrated | Ea width, CHn m0 |
| IOM density | Delarue Eq. 10 at the residual IOM's H/C (1.320 at t = 0) | g/cm³ | `bulk_density()` only | Delarue et al. 2026 | sourced (PI, 10-29); a fixed number is an option | — |
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
Casey et al. 1993), with an Ea that falls with pH (a fit to Casey et al.'s Fig. 4: 54, 26, 22 and
4 kJ/mol at pH 2.5–5.15, 0 taken at pH 6.2; zero above pH 5.85) and no neutral
branch: the rate at pH > 7 is held at its pH-7 value.

**Surface areas: UNSOURCED.** Forsterite, enstatite, diopside and tephroite are at 0.1 m²/g; albite,
anorthite, K-feldspar and troilite at 5 m²/g, so the feldspars react against 50× the area of the
mafic silicates (D7).

**Flow driving: UNSOURCED as a single value; to be swept (D12, N6).** The Simulator's Darcy flux is
q = (k/μ)·ΔP/L, with ΔP·sin(2πt/P) the tidal pore pressure (since 2026-10-30 derived from the strain, see
Exchange below; this paragraph's table used the old fixed ΔP), L = `layer_thickness` and μ = 1.8e-3 Pa·s,
so only k·ΔP/L matters. The exchanged fluid per orbit is |q| integrated over the orbit times the
core's surface area (q is per unit total area; until 2026-10-29 it was also multiplied by the
porosity, 3.1× too small), and every orbit's volume is treated as fluid that leaves the
rock and is replaced by ocean water: **the oscillating flow is assumed to give net exchange** (open
question below). Production uses k = 1e-7 m², ΔP = 1e5 Pa, L = 194.2 km. The literature constrains
ranges, each with its own limitation:
- **Liao et al. (2020)**, poroviscoelastic tidal model of the core. They consider k from 1e-12 to
  1e-7 m². Their Fig. B3 has a maximum pore pressure of ~3×10⁴ Pa at the poles and ~0 at
  mid-latitudes, with mean pore velocities of 3.8e-6 m/s (radial) and 2.2e-6 m/s (tangential). That
  is one case (η = 1e16 Pa·s, μ = 1 GPa, α = 0.1, κ = 1e-8 m², φ = 0.2). The 1e5 Pa default is ~3×
  their plotted peak.
- **Choblet et al. (2017)**, 3-D thermal convection in a tidally heated porous core. They explore k
  from 1e-15 to 1e-12 m² and prefer 1e-14 to 1e-13 m², with φ 0.1–0.3 and 10–30 GW of tidal power.
  Upwellings exceed 363 K, and the ocean is processed through the core above 363 K in 25–250 Myr.
  Limitation: steady buoyancy-driven convection, no tidal pumping, no chemistry.
- **Zandanel et al. (2021)**, 1-D reactive-transport columns of olivine only, with water velocities
  of 1, 10 and 100 cm/yr (which they take from Choblet et al. for k from 1e-15 to 1e-12 m²), at 1, 25,
  90 and 150 °C. Olivine alteration lasts up to ~75 Myr, and H₂ production (> 1e12 mol/yr) is far
  above the ≥ 1e9 mol/yr inferred from the plume. Their 0.1 mm olivine grains have a geometric
  surface area of ~0.01 m²/g, 10× below ours (D7), and their short lifetimes support a spin-up (D9).
- **Fisher et al. (2024)**, 3-D outcrop-to-outcrop hydrothermal siphons scaled to ocean-world
  gravity. Aquifer k is 1e-13 to 1e-10 m², and a siphon at Enceladus's gravity needs ~1e-11 to 1e-10
  m². W/R is 4.9–150 per Gyr at Enceladus's gravity. Limitation: an Earth ridge-flank analogue
  covering only the upper 4–5 km.

**Production k = 1e-7 m² is the top of Liao et al.'s range**, 3–6 orders of magnitude above Fisher
et al.'s aquifers and ~6 above Choblet et al.'s preferred range. Its flow, as the code computes it
(`Simulator$run_simulation()` with each k and `max_tidal_stress`, 2026-10-29):

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
`sim$get_summary_metrics()$circulation_time_years` reports it. Everything scales with k·ΔP: ΔP = 1e4
Pa is 10× slower. At k = 1e-8 m² and 3e4 Pa the displacement, 5 cm, is close to the ~7 cm implied by
Liao et al.'s 3.8e-6 m/s. Choblet et al.'s preferred k, with ΔP 1e4–1e5 Pa, gives circulation times of
7.8e7–7.8e9 yr against their 25–250 Myr. Zandanel et al.'s 1–100 cm/yr corresponds here to k ≈
2e-12 to 2e-10 m² at 1e5 Pa (2e-11 to 2e-9 at 1e4 Pa). The production case circulates the ocean ~10⁵–10⁶×
faster than Choblet et al. (This table is the Darcy term alone with the gradient over the whole layer; the exchange is now the
smallest of three volumes with the Darcy gradient over the skin depth, see Exchange below: at
production k the poroelastic term binds.)

**The sweep (N6; designed, not run).**
1. Treat the driving term as one quantity, reported three ways: k·ΔP/L, the mean Darcy flux, and
   the fluid exchanged per orbit as a fraction of the porewater. k, `max_tidal_stress` and
   `layer_thickness` already set it.
2. Ranges: ΔP 1e4–1e5 Pa (Liao et al.); k 1e-15–1e-7 m², marking the Choblet (1e-14–1e-13), Fisher
   (1e-13–1e-10) and Liao (1e-12–1e-7) sub-ranges. This brackets Zandanel et al.'s 1–100 cm/yr.
3. Report each case's ocean circulation time, for comparison with Choblet et al.'s 25–250 Myr.
4. **Open question for the PI:** how the oscillating flow becomes net exchange. Liao et al.'s
   velocities imply ~7 cm of pore-fluid displacement per tide, back and forth. The code counts
   |q| over the whole orbit as one-way exchange, an upper bound, if back-and-forth flow brings the
   same fluid back. (The extra ×φ in the flow rate was a bug, fixed 2026-10-29.) L, the layer
   thickness (194.2 km), sets the pressure gradient, where Liao et al.'s pressure varies with
   latitude over the core's radius.

**Exchange: the smallest of three volumes per half-orbit (PI, 2026-10-29).** The core is a sponge:
tides squeeze and expand the matrix, expelling and drawing in pore water (the poroelastic coupling of
Liao et al. 2020). Per half-orbit the exchanged volume is ΔV = min(ΔV_poro, ΔV_Darcy, V_pore), over
the reacting box (the skin, below):
1. **ΔV_poro = α·ε_v·V_box**, the fluid the tidal volumetric strain expels (drained limit); α is the
   Biot coefficient (below).
2. **ΔV_Darcy**, what drains through the rock in half an orbit: half the |q|·A integral over the
   orbit, i.e. (k/μ)(ΔP/L_eff)·A·T/π, with the pressure gradient over L_eff = min(L, δ), the skin
   depth below. (Written with the amplitude, (k/μ)(ΔP/L)·A·T/2, it would be π/2
   larger; the integral of the sinusoid is what the code has always used.)
3. **V_pore = φ·V_box**, saturation.

The exchange fraction is x = ΔV/V_pore per half-orbit; per orbit the two halves compound, x_orbit =
1 − (1 − x)², and the grouped steps compound that again, x_step = 1 − (1 − x_orbit)^G. Since x ≤ 1
by construction, the old cap at the pore volume is no longer needed. `sim$exchange` holds all three
volumes and the binding limit (`poro`, `darcy` or `pore_volume`); `get_summary_metrics()` reports
`exchange_fraction_per_half_orbit`, `exchange_limit` and the circulation time (ocean mass ÷ exchange
per orbit); every generated input stamps header lines with the reacting box and the exchange. They
are constant through a run (they will change between chunks once φ evolves, PROGRESS N13).

**Porosity and the Biot coefficient (PI, 2026-10-29).** The initial porosity is `phi0` (Simulator;
`porosity` is the older name), default 0.32, the PlanetProfile profile's `phiRockMax`: **one fit**
(Styczinski et al. 2023) with its own assumptions; the N6 sweep covers 0.1–0.4 (Choblet et al.'s
0.1–0.3 and the range the core density implies). φ now controls k (as φ³ under Kozeny–Carman, item 3),
α, the storage, δ and the strain-limited volume. The Biot coefficient comes from the critical-porosity
model (Nur, Mavko, Dvorkin & Galmudi 1998, *The Leading Edge* 17, 357–362; as cited by Dvorkin & Nur
2000): the frame's bulk modulus falls linearly to zero at the critical porosity φc, K_frame =
K_grain(1 − φ/φc), so **α = 1 − K_frame/K_grain = min(1, φ/φc)**, which is ≥ φ as a porous solid
requires. `phi_c` defaults to **0.4** (PI, 10-30), the critical porosity of sorted grain packs and
sandstones in Dvorkin & Nur (2000)'s table after Nur et al. (1998) (36–40%; 40% for sandstones and
carbonates, 65% for chalk), still a placeholder. **φ ≥ φc means α = 1:** the solid no longer forms a
frame and the "rock" is a suspension (granular media are loose at φc, Dvorkin & Nur 2000), so α is
capped there. `alpha = <number>` overrides it; an override below φ is an error.
Bagheri et al. (2025) give a single, fixed core bulk modulus (100 GPa, their Table 2), not frame and
grain moduli, so α cannot be computed from them. At φ = 0.32, α = 0.80 (it was a free 0.2 before,
below φ, and 0.64 with φc = 0.5). The box's solid mass is V·(ρ_core − φ·ρ_water), so the W/R follows φ.

**The core-density check.** ρ_bulk = (1 − φ)ρ_solid + φρ_fluid, and Bagheri et al. (2025) give a
present-day core density of 2160–2400 kg/m³. `phi_from_core_density(rho_solid, rho_core)` inverts it:
at t = 0 with 15 wt% IOM (ρ_solid 2.920 g/cm³) Bagheri's range implies φ = 0.27–0.40. Bagheri's value
is present-day, so it is **a check on the end state, not an input**: `bulk_density()` adds a
`core_density_ok` column and warns when the final ρ_bulk is outside 2160–2400. **φ and the organic
content are linked:** less IOM means a denser solid and, for the same core density, a higher φ (a
forsterite-only solid, 3.21 g/cm³, needs φ = 0.37–0.48).

**The strain.** ε₀ is Liao et al. (2020)'s polar strain amplitude (§2.2, Fig. 2a):
ε₀ = (9/4π)·e·ω²/(ρG)·(5/3)·|k₂*|, k₂* = (3/2)/(1 + 19μ*/(2ρgR)), with the Maxwell rigidity μ* =
μ·iωτ/(1 + iωτ), τ = η/μ, g = (4/3)πGρR. Their §2.2 text and Appendix A write k₂ with an extra factor
3/5; the Fig. 2 caption's form (without it) is the one that reproduces Fig. 2a: the elastic contours
at μ = 2.1e7, 8.0e7, 8.6e8 and 8.6e9 Pa against 2.1e7, 7.5e7, 8e8 and 8.6e9 read off the figure, and
the 10⁻⁵ contour meeting μ = 1e10 at η = 1.6e12 as in the figure (with their R₀ = 186 km, 33 h,
e = 0.0047, ρ = 2400 kg/m³); the 3/5 version is ~1.6× off. `tests/exchange_check.R` checks this.

ε₀ is a polar *radial* strain; the fluid is expelled by the *volumetric* strain. Liao et al.'s
simplified model does not give the dilatation (with a Laplace displacement potential it is zero in
their viscoelastic end-member). **Assumption (placeholder pending the PI):** the dilatation has the
degree-2 pattern (r/R)²·P₂(cos θ) with amplitude ε₀ at the pole. A degree-2 field has zero mean, so
each half-orbit an equal volume compresses (expelling fluid) and expands (drawing it in); the
compressed part's volume mean is ε_v = 0.19245 × (3/5)(1 − x⁵)/(1 − x³) × ε₀, x = (R − L)/R, where 0.19245
is half the surface mean of |P₂| (0.385: the compressed half of the pattern) and the rest the radial mean
of (r/R)² over the box. **The factor follows the reacting box** (PI, 10-30): with `reactive_depth =
"skin"` it is evaluated over the skin, 0.192 for a thin skin (r ≈ R; 0.192 at k ≤ 1e-11, where δ < 1 km)
falling to 0.138 for the 74 km skin at k = 1e-7; with `"layer"` over the whole layer, 0.115 for the
full core. `sim$exchange$strain_factor` reports it. If the matrix is nearly incompressible the dilatation is
smaller, so ΔV_poro is an upper estimate.

**Rigidity and viscosity** come from Bagheri et al. (2025) Table 5, "dissipation in the core" (µ 4–70
GPa, η 4.4e11–2.3e12 Pa·s; the scenario consistent with tidally driven flow in the core), with
**placeholders at the geometric means**, `core_rigidity` = 17 GPa and `core_viscosity` = 1e12 Pa·s,
until the PI picks values. At these viscosities the Maxwell time η/µ is seconds to minutes against a
33 h tide, so the matrix responds viscously, |µ*| ≈ ωη, and **ε₀ depends on η, not on µ**:

| η (Pa·s) | ε₀ | ε_v (full core) | ΔV_poro / V_pore per half-orbit (= ε_v/φc) |
|---|---|---|---|
| 4.4e11 | 3.8e-5 | 4.4e-6 | 1.1e-5 |
| 1e12 (placeholder) | 1.7e-5 | 2.0e-6 | 5.0e-6 |
| 2.3e12 | 7.6e-6 | 8.7e-7 | 2.2e-6 |

(any µ in 4–70 GPa; full-layer box; with α = φ/φc the ratio α·ε_v/φ = ε_v/φc does not depend on φ).
Bagheri's rigidities lie above Liao's explored 0.01–10 GPa, and Liao's viscosities (1e12–1e18) only
touch Bagheri's at 1–2.3e12; the overlap (µ 4–10 GPa, η 1–2.3e12) is where the viscous knee above was
checked. Bagheri's core radius and density range (183–213 km, 2160–2400 kg/m³) move ε₀ between 1.4e-5
and 2.1e-5 at the placeholders; the profile's 194.2 km and 2353 kg/m³, which the code uses, are inside
both ranges. **The pore-volume limit cannot bind at real porosities:** even a fluid matrix (k₂ → 3/2)
gives ε₀ ≈ 1.5e-4, so α·ε_v ≤ 1.7e-5 of the box per half-orbit, against φ = 0.1–0.4.

**The tidal pore pressure (PI, 2026-10-30).** The pressure that drives the Darcy term is no longer an
independent parameter: it is the undrained pore pressure of the same strain, **ΔP = α·ε_v·M**, with the
Biot modulus M = 1/S = (K_u − K_m)/α² of Liao et al. (2020) Table A1 (their Eq. A5 with no change of
fluid content). `max_tidal_stress = <number>` still overrides it (for reproducing old runs or a sweep).
At the placeholders (φ = 0.32, φc = 0.4, α = 0.80, M = 5.2 GPa, ε₀ = 1.73e-5) the pole amplitude
α·M·ε₀ is 7.2e4 Pa and the box mean 1.4e4 Pa in a thin skin (ε_v/ε₀ = 0.192), 1.0e4 Pa in the 69 km
skin at k = 1e-7 and 8.3e3 Pa over the whole layer (0.115). **Against Liao et al.'s ~3e4 Pa (Fig.
B3a):** the same formula at their Fig. B3 case (η = 1e16 Pa·s, µ = 1 GPa, α = 0.1, φ = 0.2, ε₀ = 8.6e-7)
gives only 1.1e3 Pa at the pole, ~28× below their plot. The ocean's tidal loading of the seafloor
(below) is far too small to close that gap; it remains unexplained (their full poroviscoelastic solution,
not reproduced here). The strain-derived ΔP does not change the exchange where the strain binds.

**Ocean tidal loading (PI, 2026-10-31).** The ocean column over the seafloor changes thickness with the
tide, and its weight change, ΔP_load = ρ_w·g_c·|u_ocean − u_core|, loads the pores. With Bagheri et al.
(2025)'s eccentricity tide (their Eq. 1, U₂ = r²ω²e[−(3/2)P₂⁰(cos θ)cos ωt + P₂²(cos θ)((3/4)cos ωt
cos 2φ + sin ωt sin 2φ)], ω = n the mean motion, here 5.289e-5 rad/s from the 118,800 s period (their
5.307e-5), e = 0.0047) and the radial displacement u_r = h₂U₂/g (their Eq. 2), at the pole (U₂ =
(3/2)r²ω²e):
- **ocean surface:** taken to move with the surface, whose whole-body Love number is Bagheri et al.'s
  Table 3, "dissipative core": |h₂| = 0.0459 ± 0.0294, lag 2.95° ± 2.69°, at R = 252.1 km, g = 0.114:
  **u_ocean = 0.50 m** at the pole (their Eq. 4's maximum over the surface, 12√(2/7)|h₂|R²ω²e/g, is
  2.16 m; 0.78–3.54 m over the ±);
- **core surface:** Table 3 gives no core Love number. The code uses the core's own h₂ = (5/3)k₂* of
  the uniform Maxwell core (Liao's formula, the one behind ε₀): |h₂| = 0.29 lagging 83° at the
  placeholders (0.64 at η = 4.4e11, 0.13 at 2.3e12), **u_core = 1.7 m** (3.7, 0.7 m);
- the difference is taken with the phases: **ΔP_load = 213 Pa at the pole** (459 and 108 Pa at the η
  range ends), and its box mean is the pole value times the strain factor: 41 Pa in a thin skin, 30 Pa
  in the 69 km skin, 25 Pa over the layer.

It is **added to the strain-derived ΔP as an amplitude** (an upper bound: the phase between the loading
and the matrix strain inside the core is not known), and the volume it stores, S·ΔP_load·L_eff·A, is added
to the strain volume, so the "poro" limit is ΔV ≈ S·(ΔP_strain + ΔP_load)·δ·A in the skin. `ocean_load =
FALSE` drops it. **At every k and η tested the load is 0.3% of the strain term** (strain 1.38e4 + load 41
= 1.38e4 Pa in a thin skin at the placeholders; 1.00e4 + 30 Pa at k = 1e-7; 8,257 + 25 Pa over the layer),
so it changes the exchange by 0.3%.

**The largest uncertainty in the exchange: the core's strain (PI, 2026-10-31).** ε₀ comes from the
isolated uniform Maxwell core, whose h₂ (0.29 at the placeholders) is six times Bagheri et al.'s
whole-body surface h₂ (0.046) and lags by 83° rather than 3°. Inside the real coupled body the core's
response is not the isolated one, and the exchange scales with ε_v. **The bracket is ε₀ × (0.046/0.29)
to ε₀**, i.e. `strain_scale` (config, default 1) from 0.16 to 1, for the N6 sweep. At the low end the
exchange per half-orbit in a thin skin falls from 8.3e-6 to 1.3e-6 and the strain-derived ΔP from 1.4e4
to 2.2e3 Pa; circulation times are ~6× longer. Bagheri et al.'s Table 3 gives no core Love number; the
PI may supply core-surface displacements from their model later, which would replace the bracket.

With ΔP derived, the three volumes are tied together: the storage volume S·ΔP·δ·A = α·ε_v·δ·A equals
the strain volume α·ε_v·V_skin to O(δ/R), and the Darcy volume is 2/f ≈ 2.3 times the storage volume
(exactly, since δ² = c/ω). **In the skin the strain therefore binds at every k**, by a factor 2.3–3.3,
and the exchange per half-orbit is ε_v·α/φ, independent of k: k sets only how deep the skin reaches.

**The skin depth.** The tidal pressure change diffuses into the rock only a skin depth δ = √(c/ω)
(Liao et al. 2020, Table A1 and Fig. B4), with their poroelastic diffusivity
c = (k/η_f)·(K_m + 4µ/3)(K_u − K_m) / (α²(K_u + 4µ/3)), K_m = (1 − α)K_s, K_u = K_m + α²K_sK_f/(φK_s +
(α − φ)K_f), K_f = 2.2 GPa and K_s = 10 GPa (their Table A1; `water_bulk_modulus`,
`rock_bulk_modulus`) and µ the matrix rigidity. (K_u − K_m)/α² is the Biot modulus M = 1/S, S =
φ/K_f + (α − φ)/K_s the storage (1.93e-10 Pa⁻¹ at φ = 0.32, α = 0.80), so c is k/(η_f S) times
f = (K_m + 4µ/3)/(K_u + 4µ/3) = 0.88. Their convention has no factor 2 (δ = √(c/ω), not √(2D/ω)); this
code follows it. δ is 22 m at k = 1e-14 m², 0.69 km at 1e-11 and 69 km at 1e-7: **δ is the depth the
exchange actually reaches**, which is what makes the well-mixed porewater box an approximation
(below). Cross-check: the volume the pressure change stores in the skin, S·ΔP·δ·A, is the δ-based
Darcy volume divided by 2/f = 2.27 at every k (exactly: (ωT/π)/f).

**Which limit binds over the whole layer** (`reactive_depth = "layer"`, φ = 0.32, α = 0.80, ΔP from the
strain = 8.3e3 Pa, placeholders):

| k (m²) | δ | ΔV_poro (m³) | ΔV_Darcy (m³) | S·ΔP·δ·A (m³) | V_pore (m³) | binds | x per half-orbit | circulation time |
|---|---|---|---|---|---|---|---|---|
| 1e-14 | 22 m | 4.9e10 | 3.8e7 | 1.7e7 | 9.8e15 | darcy | 3.8e-9 | 1.1e6 yr |
| 1e-12 | 219 m | 4.9e10 | 3.8e8 | 1.7e8 | 9.8e15 | darcy | 3.8e-8 | 1.1e5 yr |
| 1e-11 | 692 m | 4.9e10 | 1.2e9 | 5.2e8 | 9.8e15 | darcy | 1.2e-7 | 3.3e4 yr |
| 1e-7 (production) | 69 km | 4.9e10 | 1.2e11 | 5.2e10 | 9.8e15 | **poro** | 5.0e-6 | 811 yr |

With the skin depth the Darcy volume grows as √k, not k. Over the whole layer the storage volume is only
δ·A/V_layer = 3δ/R of the strain volume, so Darcy binds until δ ≈ R/7, 28 km (k ≈ 2e-8 m²).

**The reacting box is the skin (`reactive_depth = "skin"`, default; PI 2026-10-29).** When δ < L,
only the top δ of the layer exchanges in a tide, so the box that reacts is the top L_eff = min(L, δ)
(`sim$reactive_thickness()`): its temperature and pressure are the profile means over the top L_eff,
its rock mass and pore volume are those of the L_eff shell, and the exchange fraction is relative to
the skin's pore volume (the poroelastic volume uses the skin's degree-2 factor). `reactive_depth =
"layer"` is the old whole-layer box. In the skin the poroelastic limit binds at every k: the skin's
pore volume grows with δ as the Darcy volume does, and the strain expels ε_v/φc ≈ 8e-6 of it per
half-orbit, less than what drains (2.3–3.3× less, above). Full core, α = φ/φc (φc = 0.4), placeholders:

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
completed with 0 errors. **Result: at low k, tidal pumping alone reaches only cold rock just below
the seafloor** (1.4 °C, ocean pressure); deep, hot rock is reached only as δ approaches L, i.e. at the
top of the permeability range. This is the comparison point with buoyancy-driven models (Choblet et
al. 2017), whose hot upwellings tidal pumping does not produce; those models are out of scope here.
`porewater_temp` still overrides the profile temperature when given as a number.

**Porosity that evolves with the chemistry (`run_tidal_evolution()`, `porosity_evolution.R`; PI
2026-10-29/30; `porosity_evolution = TRUE` for production runs, FALSE in the quick self-tests).** Rigid
frame, no compaction: φ(t) = 1 − V_solid(t)/V_total, with V_total = V_solid(0)/(1 − φ0) fixed and
V_solid from `bulk_density()` (minerals, the remaining IOM at the Delarue density, clathrates;
`origin_pqi` carries the IOM inventory across restarts). The run is split into chunks
(`tidal_kinetic`, then `tidal_kinetic_restart` with `stop_on_depletion = FALSE` and
`allow_pressure_change = TRUE`, since the skin's pressure moves with δ); between chunks a new Simulator
takes φ(t) and recomputes k (with `k_of_phi`, e.g. Kozeny–Carman), α, the storage, ΔP, δ, the reacting
box, the strain-limited volume and the exchange fraction. A chunk whose φ changes by more than
`phi_tol` (1% relative) is rerun with half the orbits, down to one orbit, which is accepted and flagged
(`over_tol`); the next chunk is scaled from the last one's rate. φ ≤ `phi_min` (0.01) is capped and
flagged `clogged`. `tests/porosity_evolution_check.R` checks the bookkeeping.

**Refill (PI, 2026-10-30).** Hydration consumes porewater (after the first orbit the porewater was
20–37% short of the pores). Between chunks the porewater is now brought back to the pore volume φ(t)
implies, m_pw(0)·φ(t)/φ0 (m_pw(0) the first chunk's porewater, so the target is on the run's own
normalised basis and is 0 at t = 0), by moving ocean water into it: `tidal_kinetic_restart(refill_kg
=)` mixes that much of the ocean solution into the porewater (or porewater into the ocean, if negative)
before its first orbit. MIX moves whole solutions, so every element is conserved between the two; the
gas phases are untouched; the MIX steps run before the SELECTED_OUTPUT is defined and add no rows
(such a restart has no `i_soln` rows). `tests/restart_closure_check.R` checks that a restart with a
refill of +0.02 or −0.01 kg ends its orbit with the same whole-system element totals (porewater, ocean,
rock, pore and ocean gas) as without it, and that the water moved. Each chunk reports the refill
(`refill_in_kg`, and per orbit of the previous chunk) next to the tidal exchange (`exchange_kg_per_orbit`).
**The refill dwarfs the tidal exchange:** after orbit 1 it moves 0.033–0.044 kg of ocean water into the
porewater (normalised basis), against a tidal exchange of 1.9e-6 kg per orbit, ~2×10⁴ times more; after
that it falls to 1e-4–1e-5 kg per orbit, still 10–100× the tidal exchange, and is sometimes negative (a
little porewater returned) without IOM. In this model, then, **hydration-driven inflow, not tidal
pumping, sets how much ocean water enters the rock** while the primaries hydrate.

**The ocean per kg of rock (PI, 2026-10-31).** The PHREEQC basis is per kg of the reacting box's rock, so
the ocean box per kg of rock must be the ocean's mass ÷ the box's rock mass. It was 0.9914 of that at
every k (k = 1e-14: 1002.2 against 1010.9; 1e-11: 31.80 against 32.08; 1e-7: 0.4619 against 0.4659),
because it used the config's geometric ocean at 1000 kg/m³ (2.113e19 kg) instead of the profile's
`Mocean_kg` (2.131e19 kg); it now uses the profile's (the Simulator's circulation time keeps the config's,
0.9% lighter). And it was not updated between chunks: a restart carried the parent's ocean per kg of the
*parent's* box, so when δ shrank (69 → 58 km) the real ocean silently shrank with it. Each restart now
scales the parent's ocean water, ocean minerals and ocean gas by old box rock ÷ new box rock
(`tidal_kinetic_restart(ocean_scale =)`; 1.07–1.38 after orbit 1, ~1.00–1.07 later).

**Drift.** Each chunk's whole-system element totals (porewater, ocean, rock, IOM, both gases) are checked
against its start (t = 0; or the parent's end with the ocean scaled), cumulatively (`drift_max`,
`drift_element`), with a warning above 1e-4 relative; the production tolerance is unchanged. In the
24-orbit runs below it is ≤ 7e-9 at 160 °C (the solver's tolerance, Al and C), and 4–6e-6 in H at the
skin's temperature: there the skin's temperature changes between chunks, and a restart redefines the
porewater at the new temperature from its pH, which moves its H by that much. `tests/porosity_evolution_check.R`
requires < 1e-6 at 25 °C (measured 3.9e-10).

**How much ocean water the refill takes (PI, 2026-10-31).** `refill_summary(chunks)` scales each chunk's
refill to the whole core by the reacting box's rock mass (the skin is a shell over the whole core) and
compares it with the profile's ocean mass (`Mocean_kg` = 2.13e19 kg). In the 24-orbit runs below (with the
ocean scaling):

| run | after orbit 1 | cumulative, 24 orbits | fraction of the ocean | ocean concentration factor (scaled / model's own ocean) |
|---|---|---|---|---|
| skin T, 15 wt% IOM | 1.42e18 kg | 1.49e18 kg | 6.97% | 1.0749 / 1.0749 |
| skin T, no IOM | 1.41e18 kg | 1.42e18 kg | 6.65% | 1.0713 / 1.0713 |
| 160 °C, 15 wt% IOM | 1.66e18 kg | 1.82e18 kg | 8.53% | 1.0933 / 1.0933 |
| 160 °C, no IOM | 1.25e18 kg | 1.27e18 kg | 5.95% | 1.0633 / 1.0633 |

**Some 6–8% of the ocean is drawn into the core in the first orbit**, as the 69 km skin's metal corrodes
and its primaries hydrate (both nearly instantaneous here: D6 accepts instantaneous corrosion as an
idealisation of early alteration). After that the refill adds only ~0.1–0.5% of the ocean over 23
orbits, and is sometimes slightly negative without IOM. **The two concentration factors now agree** (to
all four digits in every chunk): the one from the cumulative refill, M_ocean/(M_ocean − refilled), and the
one from the model's own ocean water, M_ocean/(ocean per kg of rock × box rock). (Before the ocean
scaling they were 1.06–1.09 and 1.08–1.12.) The ocean per kg of rock at a chunk's end is then the
geometric one × 1/CF, the water the refill has drawn. The concentration factor follows from the water
loss alone; the solutes the refilled water carries into the rock are in the chemistry.

**24-orbit runs (production k = 1e-7 at φ0 = 0.32 via Kozeny–Carman with d = 1.59 cm, skin box, φc = 0.4,
ΔP from the strain, refill on):**

| | φ at orbit 24 | k | δ | skin T | circulation | refill, orbit 1 |
|---|---|---|---|---|---|---|
| skin temperature (`"profile"`), 15 wt% IOM | 0.279 | 6.0e-8 | 58 km | 31 → 27 °C | 1,121 yr | 0.033 kg |
| skin temperature, no IOM | 0.259 | 4.5e-8 | 52 km | 31 → 26 °C | 1,260 yr | 0.034 kg |
| 160 °C override, 15 wt% IOM | 0.222 | 2.6e-8 | 43 km | 160 °C | 1,626 yr | 0.044 kg |
| 160 °C override, no IOM | 0.169 | 9.9e-9 | 31 km | 160 °C | 2,651 yr | 0.038 kg |
| *before (10-29): 160 °C, φc 0.5, ΔP 1e5, no refill, IOM / none* | *0.219 / 0.165* | *2.5e-8 / 9.1e-9* | *45 / 32 km* | *160 °C* | *2,014 / 3,334 yr* | — |

(k, δ and circulation are those of the last chunk, which starts a little above the final φ.) **The test
of 10-29 overrode the temperature** (`porewater_temp = 160`) on the skin box, whose own temperature at
k = 1e-7 is 31 °C (falling to 26–27 °C as the skin thins). At the skin's temperature the core closes much
less (φ 0.32 → 0.28 with IOM, → 0.26 without), because hydration is slower; at 160 °C it closes as
before (→ 0.22, 0.17), the new defaults changing little (the strain binds throughout, so ΔP does not
matter). In every case **the core closes and does not clog** in 24 orbits, and the strain limit binds.
The first chunk takes 77 s at the skin's temperature and 190–515 s at 160 °C.

**Permeability and surface area from grain size (item 3, PI 2026-10-29).** `Simulator$new(k =
"kozeny_carman", grain_diameter = d)` sets the permeability from the grain size and the porosity by the
Kozeny–Carman relation for granular, unconsolidated media, **k = d²φ³/(180(1 − φ)²)** (Carman 1937; the
180 is 36 × Carman's constant 5 for packed spheres, with the specific surface 6/d; checked against
secondary sources, the 1937 paper itself not seen). A number for `k` still works, so old runs reproduce.
`kozeny_carman_k()`, `kozeny_carman_d()` (the inverse). The table (m²):

| φ | d = 1 µm | 30 µm | 1 mm | 2 cm |
|---|---|---|---|---|
| 0.2 | 6.9e-17 | 6.2e-14 | 6.9e-11 | 2.8e-8 |
| 0.3 | 3.1e-16 | 2.8e-13 | 3.1e-10 | 1.2e-7 |

**At the production φ = 0.32, Choblet et al.'s preferred 1e-14–1e-13 m² needs grains of 5–16 µm, and the
production k = 1e-7 m² implies grains of 1.6 cm: centimetre-scale gravel**, not a chondritic matrix.
Fractures can dominate a real rock's permeability, so Kozeny–Carman is then a lower bound. Alteration
reduces φ and k falls as φ³: with `run_tidal_evolution()` (item E) a Kozeny–Carman template follows
φ(t) through k automatically.

The surface area follows the same grain size: `ssa = "geometric"` sets each mineral's specific surface
to **6/(ρᵢ·dᵢ) × `roughness`** (ρᵢ = molar mass/Vm from `data/mineral_molar_volumes.csv`; `roughness`,
default 1, is the BET-over-geometric factor), with `grain_diameter_by_mineral` (a named vector, m) for
per-mineral sizes, e.g. a fine matrix for feldspars and sulfides and coarser olivine and pyroxene;
`PhreeqcIntegrator$apply_ssa()` writes it into every generated input. `ssa = "legacy"` (the default)
keeps the current 0.1 and 5 m²/g (UNSOURCED, D7). The difference is large: at d = 1.6 cm (production k)
the geometric area of forsterite is 1.2e-4 m²/g, 800× below the legacy 0.1; at 10 µm it is 0.19 m²/g.
`tests/kozeny_carman_check.R` checks the table, the Simulator and the SSA lines.

**The sweep (N6; design only, not run): one sweep over (φ, d, strain).** It replaces the separate k,
surface-area and stress sweeps: φ 0.1–0.4, d 1 µm–1 cm (setting both k and the surface areas), and,
since ΔP is now derived from the strain (F), `strain_scale` 0.16–1 (L) in place of ΔP. k (m²) over the
grid, marked C (Choblet et al.'s preferred 1e-14–1e-13), F (Fisher et al.'s aquifers 1e-13–1e-10) and
L (Liao et al.'s 1e-12–1e-7):

| φ | 1 µm | 10 µm | 100 µm | 1 mm | 1 cm |
|---|---|---|---|---|---|
| 0.1 | 6.9e-18 | 6.9e-16 | 6.9e-14 C | 6.9e-12 F L | 6.9e-10 L |
| 0.2 | 6.9e-17 | 6.9e-15 | 6.9e-13 F | 6.9e-11 F L | 6.9e-9 L |
| 0.3 | 3.1e-16 | 3.1e-14 C | 3.1e-12 F L | 3.1e-10 L | 3.1e-8 L |
| 0.4 | 9.9e-16 | 9.9e-14 C | 9.9e-12 F L | 9.9e-10 L | 9.9e-8 L |

**Two assumptions, stated, not swept.**
- **Ingested water has the bulk-ocean composition.** The water drawn in is not the fluid just
  expelled (bottom currents disperse that), but it is not fresh either: the ocean box carries the
  products of every earlier cycle. That is what MIX with the evolving ocean does.
- **The porewater is one well-mixed box.** Real exchange penetrates about the displacement distance
  per tide (centimetres in Liao et al.'s case), the pressure change reaches only the skin depth δ
  (~20 m to ~70 km over k = 1e-14 to 1e-7 m²), and both spread by dispersion. Over Myr runs this is
  an approximation.

**Fe/Ni metal: instantaneous corrosion, an accepted idealisation (D6 closed, PI 2026-10-29).** The metal
is a `dissolve_only` equilibrium phase, so it corrodes in the first step: ~2 mol H₂ per kg of rock arrives
at t = 0 (1.998 mol of Fe metal in the corrected rock; the Ni corrodes at 160 °C but not at 25 °C). This
is not a rate claim: it stands for the alteration and outgassing of early differentiation, and a
pure-water ocean over fresh rock is likewise an idealisation, chosen so that the early outgassing totals
can be computed cleanly. **The science is the evolution over millions of years: the first orbits are a
transient**, and results are to be judged on the long runs (PROGRESS N8).

**The H₂ budget (item 6).** `h2_budget(csv)` gives per orbit, on the run's basis (mol per kg of the
box's rock): the H₂ generated since t = 0 by the metal and by schreibersite, separately, counted as
electrons/2 to Fe(II), Ni(II) and phosphate (1 H₂ per Fe or Ni, 5.5 per Fe₃P entered as zero-valent;
lower bounds, since iron oxidised further, e.g. to magnetite, makes more), and where the H₂ is at the
end of the orbit: dissolved in the ocean, in the ocean gas, dissolved in the porewater, in the pore gas.
`scale_to_core(budget, sim)` multiplies by the rock mass of the reacting box, the whole permeable layer
and the whole core (V (ρ_core − φ ρ_water) from the Simulator's geometry): the inventory available for
early loss to space. Lucas's 160 °C case (k = 1e-7, skin box 69 km, 12 orbits, `to_ocean`):

| per kg of rock | orbit 1 | orbit 12 | whole core, orbit 12 |
|---|---|---|---|
| H₂ generated by the metal | 1.998 mol | 1.998 mol | 1.25e20 mol |
| H₂ generated by schreibersite | 0.0036 mol | 0.039 mol | 2.4e18 mol |
| H₂ in the ocean gas | 2.31 mol | 2.43 mol | 1.51e20 mol |
| H₂ dissolved in the ocean | 0.030 mol | 0.030 mol | 1.9e18 mol |
| H₂ dissolved in the porewater | 0.0094 mol | 0.0071 mol | 4.4e17 mol |

(the whole core's rock is 6.24e19 kg; here the layer is the whole core). The H₂ present (2.47 mol/kg)
exceeds the metal-and-schreibersite count (2.04): the rest comes from iron oxidised beyond Fe(II) (the
magnetite and Fe(III) phases), which the lower-bound count leaves out. So **~1.5e20 mol of H₂ (3e17 kg)
sits in the ocean gas within the first orbits**, the inventory available for early loss.



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
send the excess to the ocean), is implemented since 10-31 (D6 closed; README §7, Gas).

**Pressure** from the profile. Pressure acts only through gases: under the LLNL aqueous model
aqueous molar volumes are inactive, and the minerals' were stripped (`databases/README.md`), so
mineral equilibria are at the database's reference pressure (~0.1–0.2 log units at 260 bar).

**Initial ocean:** pure water by default, so all ocean carbon is attributable to the organics. The
cometary option follows Miller et al. (2026) and Rubin et al. (2019); it starts the ocean with
~0.087 molal DIC, so attribution must use the paired difference.

**The cometary-ice ocean, extended (PI, 2026-10-29; an option, never the default).**
`cometary_ice_ocean(v = 0.02, f = 1, include_NH3 = TRUE, include_less_volatile = FALSE, include_O2 =
TRUE)`.
- **NH₃** is now included by default: Rubin et al. (2019) Table 2 gives 0.67 ± 0.20 per 100 H₂O, and
  Miller et al. (2026) apply v_NH3 = 0.02 to "NH₃ and other compounds that are more volatile", so the
  depleted ice has 0.0134 per 100 H₂O (7.44e-3 molal); its absence from their Table 5 is treated as an
  omission. Rubin et al. note a factor of 10 for NH₃ between the 2015 bulk value and the earlier 3.1 au
  measurements (0.06 and 0.15 over the two hemispheres, which themselves differ 2.5×), left unexplained.
- **f** scales the whole cometary load after the depletion: the fraction of the accreted cometary
  volatiles retained in the ocean.
- **Less volatile ices.** Miller et al.'s rule keeps ices less volatile than NH₃ at 273 K in full, and
  their Table 5 includes none. `data/comet67p_volatility.csv` lists every species of Rubin et al.'s
  Table 2 (bulk column, with its ± errors), classified against NH₃ at 273 K: by the sublimation curves
  of Fray & Schmitt (2009) (their Table 5 polynomials evaluated at 273 K, read from the PDF; NH₃ is 14
  bar there) where they cover the species, else by the normal boiling point against NH₃'s 239.8 K, a
  proxy marked as such (standard values; not checked against a handbook here). Fray & Schmitt make CO,
  O₂, N₂, CH₄, C₂H₆, CO₂, H₂S and the noble gases more volatile, and SO₂, HCN, CH₃OH, HC₃N, C₆H₆ and
  HCOOH less. Radicals and atoms (SO, S, S₂, H₂CS, PO) are left unclassified (no boiling point; Rubin et
  al. tie S, S₂ and SO partly to grains); OCS, C₃H₈, HCl and HBr come out more volatile by the proxy
  (OCS and C₃H₈ are then missing from Table 5 too). **The 21 less-volatile CHONS species add C 1.08,
  N 0.177 and S 0.177 per 100 H₂O, against Table 5's C 0.174, N 0.00356, S 0.022: 6.2×, 50× and 8×**
  (HF, also less volatile, is not CHONS and is not added). They enter by their elements (as CO₂, NH₃ and
  H₂S with H and O corrections, so no element name can match a phase, the CO-as-cobalt trap): under
  `c4_only` their carbon can only become C(+4), and the electrons go to H₂ or reduce sulfate; HCN's N
  becomes N(−3). **Caveat:** in a cold ocean these organics are kinetically metastable, so converting
  them to CO₂ plus reductant is an equilibrium assumption. Under `"coupled"` the equilibrium decides as
  well: no species is held as itself (methanol, formate and acetate exist in the database but
  equilibrate), which the PI's request ("stay as themselves") cannot be in an equilibrium model.

The starting ocean (v = 0.02, NH₃ and O₂ included; 1.45 °C, 66.7 atm; mol/kgw; clathrate SIs from the
initial-ocean step, which has no clathrate phases to precipitate):

| less volatile | mode | f | pH | DIC (C+4) | H₂ | sulfate | sulfide | Mtg / reduced C | gas at 66.7 atm | clathrate SI (CH₄/Mtg, CO₂, H₂S) |
|---|---|---|---|---|---|---|---|---|---|---|
| no | c4_only | 1 | 2.08 | 0.0866 | 6e-14 | 9.4e-3 | 2.8e-3 | Mtg 0.0102 | none | −0.81, −1.07, −1.88 |
| no | c4_only | 0.1 | 2.98 | 8.7e-3 | 1e-13 | 9.4e-4 | 2.8e-4 | Mtg 1.0e-3 | none | −1.81, −2.07, −2.88 |
| no | coupled | 1 | 5.43 | 0.0960 | 1e-10 | 3e-9 | 0.0122 | CH₄ etc. 2.8e-6 | none | −1.92, −1.06, −1.25 |
| no | coupled | 0.1 | 5.45 | 9.6e-3 | 1e-10 | 3e-10 | 1.2e-3 | 2.9e-8 | none | −2.91, −2.06, −2.25 |
| yes | c4_only | 1 | 5.79 | 0.557 | 0.055 | 0 | 0.099 | Mtg 1.3e-3 | **0.19 L** (H₂ 0.48, CO₂ 0.064, H₂S 6e-3, Mtg 4e-3 mol) | −1.73, −0.37, −0.38 |
| yes | c4_only | 0.1 | 5.76 | 0.068 | 0.064 | 0 | 0.011 | Mtg 8e-4 | **7e-3 L** (mostly H₂) | −1.91, −1.25, −1.31 |
| yes | coupled | 1 | 5.96 | 0.402 | 4e-10 | 1e-8 | 0.108 | 0.033 (CH₄) | **0.016 L** (CH₄ 0.051 mol) | **+0.30**, −0.55, −0.34 |
| yes | coupled | 0.1 | 6.01 | 0.042 | 5e-10 | 5e-10 | 0.011 | 1.0e-3 | none | −0.39, −1.51, −1.31 |

With the less volatile ices the organics' electrons consume all the O₂ (no sulfate is left), the ocean is
saturated in H₂ under `c4_only` (~0.06 molal at 66.7 atm, the rest a bubble), and under `coupled` the
carbon reduces to methane instead, oversaturating CH₄ hydrate at f = 1. The less volatile organics
multiply the ocean's DIC by 6–7 (0.087 → 0.56 molal at f = 1).

**Initial fluids** are defined at 25 °C, pe 4 (no oxidant or reductant) and heated by
REACTION_TEMPERATURE; defined at the run temperature, pe 4 put O₂ into the porewater (4.9 molal
O(0) at 250 °C).

**Trace N and S seeds: tested, not added (item 5, 10-31).** The default ocean is pure water plus the C seed.
Seeding N or S as well (1e-14 mol/kgw in every initial solution, alone or together) changed nothing that
matters in Lucas's case at 25 and 160 °C (2 orbits, 15 wt% IOM): every variant ran with 0 errors with and
without them; at 160 °C the solver's retries were 33 in every variant; at 25 °C they moved without
pattern (4 with C alone, 18 with C + N, 12 with C + S, 4 with all three); runtimes differed only with the
machine's load. With the C seed and the ocean step's `-diagonal_scale` (B) the default already converges,
so `trace_N_molal` and `trace_S_molal` were not added.

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
  channel's `-formula` is what the solid loses). Its density (default `iom_density =
  "delarue2026"`) is Delarue et al. (2026) Eq. 10, ρ = 1472 − 345 (H/C) + 94 / ((H/C) + 0.122)
  kg/m³, evaluated on each row at the residual IOM's atomic H/C: the bulk IOM's C and H
  (`iom_bulk_formula()`, H/C 0.753) less what the channels have released (Mtg counted as CH₄).
  At t = 0 that is 1.320 g/cm³. The fit is to kerogens and coals and is within 100 kg/m³ for H/C
  below 1 (Delarue et al. Fig. 5). `iom_density = 1.4` (any number) fixes it instead. Compression
  is left out: their bulk modulus K0 = 10.2 GPa gives 0.25% at the core's 26 MPa.
- `rho_bulk` = (1 − φ) ρ_solid + φ ρ_pore, with φ the model's porosity (0.32) and
  ρ_pore = (1 − s_g) ρ_fluid + s_g ρ_gas, s_g = min(Vgas/Vpore, 1) the pore gas saturation. ρ_fluid
  is PHREEQC's solution density (`Rho_soln`): under the LLNL aqueous model PHREEQC does not compute
  water's density from T and P, so it is ~1 g/cm³ for water at any temperature (give
  `fluid_density =` to override).
- Gas beyond the pore volume is not counted (s_g is capped at 1; it is not rock volume either).
  Rows where Vgas > Vpore are flagged in the `gas_exceeds_pores` column (with `Vgas_Vpore`), and
  `bulk_density()` warns when there are any: their ρ_bulk is for pores full of gas, not a physical
  state. In Lucas's 160 °C case that is the first REACT row under `"to_ocean"` (2.14 pore volumes,
  before the transfer) and every row under `"retain"` (up to 2.26).

For Lucas's 160 °C case (15 wt% IOM, φ = 0.32): ρ_solid 2.920 g/cm³ and ρ_bulk 2.305 g/cm³ at t = 0;
at orbit 12, ρ_solid 2.759 and ρ_bulk 2.255 (`"to_ocean"`) or 1.880 (`"retain"`, the bubble filling
the pores). With the old fixed 1.4 g/cm³ these were 2.977 and 2.344 at t = 0 and 2.806 and 2.287 or
1.912 at orbit 12: the Delarue density lowers ρ_solid by 1.9% at t = 0 and 1.7% at orbit 12, and
ρ_bulk by 1.4–1.7%. The residual IOM's H/C hardly moves in 12 orbits (0.75298 to 0.75298; 0.01 g
of 150 g is released), so its density stays 1.3197 g/cm³. Lizardite (26% of the solid volume by orbit 12) takes up water. The profile's own
silicate-layer density is 2353 kg/m³ (`rhoSil_kgm3`).

**Comparing with Enceladus.** A published core density is a bulk density (porosity and pore fluid
included), so compare it with `rho_bulk`. Bagheri et al. (2025, Table 5) fit Enceladus's libration,
shape, heat flux, gravity and mass and find a core density of 2160–2400 kg/m³ and a core radius of
183–213 km (their priors were 2000–2800 kg/m³ and 120–220 km), in both their core- and
shell-dissipation cases. Against that range:
- the model's ρ_bulk is 2.305 g/cm³ at t = 0 and 2.255 at orbit 12 under `"to_ocean"`: inside it;
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

`Vgas_Vpore` > 1 on porewater rows means the bubble is larger than the pore space: under `"retain"`
that is an upper bound, not a physical state. **`pore_gas_fate = "retain_to_pore_volume"`** (implemented
10-31, now that D6 is closed; `tidal_kinetic` and its restart; `tidal_equil` does not take it) keeps the
pore gas up to the pore volume and sends the excess to the ocean: the porewater step `PUT`s each
component's excess, gas × max(0, 1 − Vpore/Vgas); a short step after it (not punched in the run's CSV)
removes that excess from the porewater and its gas phase (KINETICS 4); the ocean step adds it to the
ocean (KINETICS 3, as `to_ocean`). The pore gas left after each removal is written to
`<base>_poregas.csv`, which the mass balance (`tests/gas_mass_balance_check.R`) and restarts read. The
default stays `"to_ocean"`.

**Worked example 1 under the three fates** (25 °C, 2 orbits, end of orbit 2):

| fate | ocean pH | ocean H₂ (molal) | ocean gas | porewater pH | pore gas |
|---|---|---|---|---|---|
| `to_ocean` | 7.57 | 0.065 (saturated) | 0.90 L | 11.04 | none |
| `retain_to_pore_volume` | 9.55 | 0.065 (saturated) | 0.58 L | 14.06 | 0.67 pore volumes |
| `retain` | 9.55 | 5.7e-7 | none | 14.06 | 2.46 pore volumes |

Under `"to_ocean"` the whole pore bubble reaches the ocean every orbit, whatever the transport rate, so
the ocean's redox is set by the gas, not by the circulation. Under `"retain_to_pore_volume"` the excess
alone (the metal's H₂ beyond 1 pore volume) still saturates the ocean in H₂; only under `"retain"` does
the ocean see the H₂ through the tidal exchange alone. The ocean gas phase (MIX 4 rows, `g_` columns) holds
H₂ and organic Mtg under `"to_ocean"`; `organic_contribution()` reports it.

### Steady state (item 7, PI 2026-10-29)

The water the porewater takes in has the **evolving bulk-ocean composition**: each orbit's products feed
the next. The expected end state is therefore a **steady state, an ocean composition that stops
changing, but not necessarily zero net exchange**, for two reasons: the IOM degrades kinetically, with a
long tail of high-Ea pools, so its release may never quite stop; and the porewater (100–160 °C) and the
ocean (1.45 °C) cannot share one saturated mineral assemblage, so circulation can sustain a constant flux
(e.g. silica dissolved in the pores and precipitated in the ocean). True equilibrium with no net exchange
needs the IOM exhausted, no kinetic primaries left and equal temperatures.

`steady_state_diagnostics(csv)` shows the approach directly, per orbit: the net elemental flux into the
ocean (the change of the ocean side, water + minerals + gas, over the orbit; the ocean exchanges only
with the porewater) and the net change of the secondary minerals in the porewater and in the ocean (mol
of formula units per kg of the box's rock). In Lucas's 160 °C case (12 orbits): the porewater's
secondary minerals grow by 1.09 mol in orbit 1, 0.26 in orbit 2, 0.014 in orbit 6 and 0.0038 in orbit
12 (the hydration transient dying away); the ocean's by 4e-8 per orbit from orbit 6; the Na flux into
the ocean settles at 2.4e-6 per orbit (2.9e-6 in orbit 1), Si falls from 1.2e-7 to 2.1e-8, while C
rises from 2.2e-10 to 2.5e-9 per orbit as the IOM's carbon release grows. Twelve orbits are far from
the steady state; the long runs (N8) are for that.

### Carbon

Use `carbon_by_oxidation_state()` (`phreeqc_tools.R`), not PHREEQC's valence totals, to see which
carbon is reduced: PHREEQC files formate, acetate and other organics under C(4).

**Comparing with Enceladus's ocean carbonate.** Glein & Waite (2020) adopt Postberg et al.
(2009)'s salt envelopes, 0.2 molal Cl with 0.1 molal total dissolved carbonate (TDC = HCO₃⁻ +
CO₃²⁻, without CO₂(aq)) or 0.05 molal Cl with 0.01 molal TDC, and the pH range 8.5–9, which gives
log a(CO₂) between −4.6 and −3.2 (their preferred estimate). The model's oceans:
- the cometary-ice option starts with 0.087 molal DIC, but at pH 1.85 (2.08 with NH₃, the default since
  10-31) it is almost all CO₂(aq) (0.0866 molal, log a(CO₂) ≈ −1.1); TDC is 2e-6 molal and there is no
  Cl; with the less volatile ices it is 0.56 molal at pH 5.8;
- the pure-water default, after 12 orbits of Lucas's 160 °C case, holds 1.3e-7 molal DIC at pH
  11.3, of which the organics contribute 1.29e-7 molal. Whether the organic contribution can
  approach the observed carbonate over the core's lifetime is what the long runs (PROGRESS N8) are
  for.

## 8. Known limitations and open decisions

`docs/PROGRESS.md` Part 10 lists the decisions (D6 metal corrosion rate, D7 surface areas, D9 spin-up,
D10 IOM sulfur, the default `pore_gas_fate`, NH₃ in the cometary ocean) and the next actions. In
short: the H₂ pulse is the metal (D6, closed: an accepted idealisation); surface areas are unsourced (D7); the 25 °C brine exceeds the
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
- Bagheri, A., Simons, M., Park, R.S., Berne, A., Hemingway, D., Melwani Daswani, M., Vance, S.D.
  (2025). Exploring the interior structure and mode of tidal heating in Enceladus. *Planet. Sci. J.*
  6, 245. doi:10.3847/PSJ/ae0cab
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
