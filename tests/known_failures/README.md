# tests/known_failures/

Inputs that fail in PHREEQC 3.8.6 for reasons outside ENTICES' control, kept so the failures can be
re-checked on a new PHREEQC version and reported upstream. Not part of the self-tests.

## `phreeqc386_mix_kinetics_segfault.pqi`: access violation in a MIX + KINETICS step

**Where it came from.** In the 2026-10-28 generator matrix (`docs/MATRIX_2026-10-28.md`), the
seeded cometary `tidal_equil` runs with `carbon_mode = "c4_only"` and `pore_gas_fate = "to_ocean"`
crash PHREEQC 3.8.6 at 25, 100 and 160 °C (exit 139, segmentation fault). The crash is in the
ocean step (MIX 4 with KINETICS 3, the `to_ocean` gas transfer), during the solver's convergence
retries ("Trying reduced tolerance 1e-18 ..."), so PHREEQC stops without an error message.
`tidal_kinetic` runs the same combination without crashing.

**The reproducer.** It was reduced automatically from the 25 °C generated input (623 lines)
to 20 lines in 309 PHREEQC runs. Each removal was kept only if the input still crashed. It was then
moved to the stock `llnl.dat` (Chamosite-7A for the project database's Chamosite):
1. Solution 1 equilibrates with Fe metal, albite and Chamosite-7A, forming H₂ in a gas phase.
   USER_PUNCH `PUT`s the gas's H₂ moles.
2. MIX 4 takes 3.3e-5 of solution 1 with a KINETICS reactant that adds those H₂ moles (`GET`) in
   one step.

The run crashes every time. Two controls show what it needs:
- without Chamosite-7A, the run completes (exit 0);
- without KINETICS 3, the run completes (exit 0).

The reduced input has lost the ocean solution and most phases, so it shows the same symptom (an
access violation during the convergence retries of a MIX + KINETICS step). It has not been shown to
follow the same code path as the production crash.

To re-check:

```sh
phreeqc tests/known_failures/phreeqc386_mix_kinetics_segfault.pqi out.txt "<PHREEQC>/database/llnl.dat"
echo $?   # 139 under Git Bash (0xC0000005 on Windows) while the bug is present
```

To regenerate the full input (from the repo root, in R):

```r
source("ENTICES_LF_20260917.R")
sim <- Simulator$new(k = 1e-7, layer_thickness = 194200, config = EnceladusConfig$new()); sim$run_simulation()
ph <- PhreeqcIntegrator$new(sim); ph$calculate_phreeqc_parameters()
ph$generate_phreeqc_input_tidal_equil(output_dir = "phreeqc_inputs/segfault", n_orbital_cycles = 1,
                                      porewater_temp = 25, carbon_mode = "c4_only", pore_gas_fate = "to_ocean",
                                      initial_ocean = cometary_ice_ocean(), trace_C_molal = 1e-14)
```

Workaround: use `tidal_kinetic` for that combination, or `pore_gas_fate = "retain"`.
