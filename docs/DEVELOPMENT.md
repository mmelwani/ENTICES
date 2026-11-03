# Developing ENTICES

Test, comparison and restart internals, and the known failures. [README.md](../README.md) covers running
ENTICES, [PRODUCTION_GUIDE.md](../PRODUCTION_GUIDE.md) production runs, and [MODEL.md](MODEL.md) the model.
`tests/README.md` describes each self-test, and `tests/run_all.R` runs them all.

## Comparing two runs column by column

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

## Restart internals

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

## Known failures

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
