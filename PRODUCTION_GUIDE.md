# Production guide

For running ENTICES production cases. You edit case files, run scripts and make figures. You do not need to
edit model code, and should not. Every command runs from the repository root.

## 1. Setup

1. Get the code: `git clone https://github.com/mmelwani/ENTICES.git` (or `git pull` in your copy).
   **On Windows, clone near the drive root** (e.g. `C:\ENTICES`). PHREEQC cannot write files whose full
   path is over 260 characters, and the run files have long names. `run_case.R` refuses a run directory
   whose path is too long.
2. Install the software:
   - R 4.5 or later, with the packages `tidyverse`, `glue`, `R6` and `patchwork`;
   - PHREEQC 3.8.6 (USGS). It is found at its default Windows install path, on the `PATH`, or through
     the `PHREEQC_EXE` environment variable. Each run checks the version before its first chunk and
     records it in `run.log`.

   The kit is tested on Windows, Linux and macOS.
3. Run the test suite: `Rscript tests/run_all.R`.
   - It takes a few hours under "apatite" (20 tests). `--quick` skips the long tests and takes a few minutes.
   - Passing ends with `ALL n TESTS PASSED`.
   - If a test fails, it prints `n OF m TESTS FAILED`, and its log is in `tests/output/run_all/`.
     Don't run production then, and send the log (section 7).
4. For long runs, keep the machine awake. Set the power plan so it never sleeps on mains power. A
   sleeping machine pauses the runs. They continue afterwards, but the time is lost.

## 2. Run the reference case

```
Rscript run_case.R cases/reference.R
```

This runs the organic case and its **same-rock blank** (the same minerals, no organics) side by side.
The blank follows the organic run's porosity, permeability and time steps, so the two can be subtracted at
every output time.

### When it stops

A run ends at whichever comes first:
- a steady state, when the ocean's DIC as C(+4), pH and Na each change by less than 1% over one ocean
  circulation time;
- t_end, 1e8 yr;
- clogging, when the porosity falls to `phi_min`;
- φ reaching φc.

`RUNINFO.txt` and the `end_reason` column of `summary.csv` say which ended it.

### How long

Expect hours to weeks. With the phosphate phases the reference takes ~3 h to 1e4 yr, ~0.5 day to 1e6 yr
and an estimated 6–24 days to the 1e8 yr backstop ([docs/MODEL.md §1.5](docs/MODEL.md#15-long-run-stepping), Step limits). Each case runs two PHREEQC
processes. `runs/reference/organic/run.log` gets one line per chunk, so watch that to see progress.

### Resuming

Rerun the same command, and it continues from the last complete chunk. This works after a crash, a
reboot or a power cut, and also if you stop the run with Ctrl-C. The run redoes no completed chunk and
repeats none.

Rerunning also **extends** a run that ended at `t_end`. Raise `t_end_yr` in the case file and rerun, and
the run continues from its last chunk (the blank with it). A run that ended at a steady state, by
clogging or at φc stays as it is, and the rerun says so and skips it. Lowering `t_end_yr` below what a run
has reached neither truncates nor reruns it. The rerun says the run is already past that time.

## 3. Make a new case

1. Copy a sweep case, e.g. `cases/sweep/iom_3.R`, to `cases/sweep/<your_name>.R`.
2. Set `name` (it becomes `runs/<name>/`).
3. Override only the fields you change. Everything else is inherited from `cases/reference.R`.

`cases/reference.R` lists every field, with units, the allowed range and the section of
[docs/MODEL.md](docs/MODEL.md) that justifies it. `run_case.R` checks the case before running and lists every problem in plain words: out of
range, φ0 ≥ φc, an unknown field, and so on.

Some details:
- To give a fixed permeability instead of a grain size, set `k_m2 = 1e-12, grain_diameter_m = NULL`.
- `label` is free text, copied into `RUNINFO.txt` and `summary.csv`.
- Don't reuse the name of a finished run for a different case. `run_case.R` refuses to mix them, so rename
  the case or move the old run away.

## 4. Run the sweep

```
Rscript run_sweep.R                 # every case in cases/sweep/, 2 cases at a time
Rscript run_sweep.R --parallel 4    # 4 cases at a time
```

Each case uses two cores (organic and blank), so `--parallel 4` uses 8 cores. Leave a couple free. Rerun
the same command to resume. Finished cases are skipped, and unfinished ones resume. The progress log is
`runs/sweep.log`.

These cases are set apart on purpose:
- `k_1e-14` has `ss_tol = NULL` and is labelled "not at steady state by design". Its ocean circulation
  time (~5e8 yr) is longer than t_end, so it runs to the 1e8 yr backstop.
- `k_floor_1e-14` and `k_floor_1e-12` keep a fracture-maintained permeability once the rock clogs,
  instead of ending the run. See [docs/MODEL.md §1.4](docs/MODEL.md#14-porosity-evolution), *Clogging and k_floor*.

### Staged runs (recommended)

A full sweep to the 1e8 yr backstop is weeks of machine time per case. Run it in two stages instead:

1. Run every case to 1e6 yr first. In `cases/reference.R` set `t_end_yr = 1e6` (the sweep cases inherit
   it), then run the sweep. Under the default `phosphate_phases = "apatite"` this takes about half a day
   per case (measured on the reference: 2.8 h to 1e4 yr, projected ~0.5 day to 1e6 yr).
2. Then extend the cases that matter toward 1e8 yr. Add `t_end_yr = 1e8` to those case files (or put it
   back in `cases/reference.R` to extend them all) and rerun `run_case.R` on each, or rerun the sweep.
   Each run resumes from its last chunk, and nothing is redone. Expect 6–24 days per case to 1e8 yr, ~6
   days once the step limit takes over and ~24 if porosity keeps limiting the steps ([docs/MODEL.md §1.5](docs/MODEL.md#15-long-run-stepping), *Step
   limits*).

Runs that ended at a steady state, by clogging or at φc in stage 1 are finished. Stage 2 skips them, and
the sweep log says why. Look at the stage-1 figures before choosing which cases to extend.

## 5. Make figures

```
Rscript make_figures.R                               # every case in runs/, into figures/
Rscript make_figures.R reference k_1e-11 iom_3       # chosen cases
Rscript make_figures.R --species pH,DIC_C4,H2,NH4_NH3 --output-var ocean_H2
```

This writes six figures as PDF and PNG:

| Figure | Shows |
|---|---|
| 1 | skin depth, skin temperature and circulation time against permeability |
| 2 | ocean time series (organic solid, blank dashed) |
| 3 | organic − blank ocean carbon |
| 4 | φ, k and skin depth |
| 5 | the sensitivity of one output, each case against the reference |
| 6 | the H₂ budget for the whole core |

### Cases per figure

One figure shows at most 8 cases besides the reference. The tornado figure (5) shows them all, and needs
the `reference` run plus at least one other case. Without them it is skipped, with a message.

### Species

`--species` takes the names after `ocean_` in `summary.csv`: `pH`, `pe`, `DIC_C4`, `Mtg`, `H2`,
`NH4_NH3`, `SO4`, `S_minus2`, `Na`, `Cl`, `K`, `Ca`, `Mg`, `Si`, `P`, `ionic_strength`.

### In R

`source("make_figures.R")` gives one function per figure: `fig_skin_vs_k()`, `fig_ocean()`,
`fig_org_minus_blank()`, `fig_rock()`, `fig_tornado()` and `fig_h2_budget()`. Each takes `cases`, and the
relevant ones take `species` or `output`. Each also takes `overlay`, a data frame of your own bands (`ymin`,
`ymax`) or points (`x`, `y`), with optional `panel` and `label` columns. The header of `make_figures.R` gives
the details. The figures contain no observational data, so you add those through `overlay`.

## 6. Where the outputs are

```
runs/<name>/
  summary.csv          one row per output time (the end of each chunk, log-spaced in time)
  summary_units.csv    unit and meaning of every column of summary.csv
  RUNINFO.txt          the commit, PHREEQC version, database, machine, wall time, what ended the run, the case
  case_file.R          a copy of the case file as run
  organic/  blank/     per run:
      run.log          one line per chunk, warnings, errors
      chunks.csv       the chunk table (porosity, k, exchange, drift, ...)
      state.rds        the checkpoint used to resume
      c0001/ c0002/ …  each chunk's PHREEQC input (.pqi), output (.out), selected output (.csv), DUMP
```

`summary.csv` contains these columns:
- the ocean: pH, pe, DIC as C(+4), CH₄ (Mtg), H₂, NH₄⁺ + NH₃, SO₄²⁻, ΣS(−II), Na, Cl, K, Ca, Mg, Si, P
  and the ionic strength, each with `blank_` and `org_minus_blank_` versions;
- the gases: moles in the ocean and pore gas phases;
- the H₂ budget, per kg of the reacting box's rock, and `_core_mol`, scaled to the whole core;
- rock and transport: φ, k, skin depth, skin T and P, α, K_d, the exchange and refill per orbit,
  ρ_solid and ρ_bulk with the Bagheri range check, and clogging;
- diagnostics: the drift, the steady-state flag, `end_reason` and `clog_time_yr`.

Units are in `summary_units.csv`.

## 7. If something fails

### Where to look

- `run_case.R` stops with a message. The run's own log is `runs/<name>/organic/run.log` (or `blank/`).
- If a process crashed outright, see `console.log` in the same folder.
- A failed run leaves a file called `FAILED` containing the error.

### Common causes

- A path too long (Windows): clone nearer the drive root.
- PHREEQC not found: set `PHREEQC_EXE`.
- "already being run by process N": that case is still running. Don't start it twice.
- "holds a run of a different case": you changed a case that has already run. Rename it. (Changing
  `t_end_yr`, `ss_tol` or the label is allowed, see *Resuming*.)
- "not running and has not finished": the run wrote nothing for 10 minutes and its process is gone. A
  busy machine doesn't cause this. A run is given up on only when it has written nothing for 10 minutes
  and the process table, asked again and again over at least a minute, says its process is gone (no
  answer counts as running). The process crashed or was killed. See its `console.log`, then rerun to
  resume.
- Status 137 (or "killed"): the operating system stopped PHREEQC because it ran out of memory.
  Send the files below.
- A PHREEQC crash at one chunk. Rerunning resumes and retries that chunk. If it fails again at the
  same chunk, send the files below.
- Known PHREEQC 3.8.6 failures are described in `tests/known_failures/README.md`. The kit avoids them.

### What to send back

1. the case file;
2. `runs/<name>/RUNINFO.txt` (if it exists);
3. the last ~50 lines of the `run.log` that failed, plus its `FAILED` and `console.log`;
4. for a test failure, the test's log from `tests/output/run_all/`.

**Don't change model code, tolerances or the accuracy settings (`phi_tol`, `x_step_max`) to get past a
failure.** Those are decisions (section 8).

## 8. What is fixed, and why

The physics choices are documented decisions, not tuning knobs, and [docs/MODEL.md §1](docs/MODEL.md#1-the-model-in-the-order-a-run-uses-it) gives the reason for each.
They are:
- the rigid-frame porosity evolution, the critical-porosity Biot coefficient and K_d;
- the layered tidal strain and the ocean load;
- the min-of-three exchange and the reacting skin;
- the refill;
- the carbon mode;
- the gas fates;
- the metal corrosion idealisation;
- the IOM kinetics and density.

The time-step limits (`phi_tol`, `x_step_max`) were set by validation against much finer steps ([docs/MODEL.md §1.5](docs/MODEL.md#15-long-run-stepping),
*Long runs*). Loosening them changes the answer by more than 1%. The test tolerances guard mass balance
and reproducibility. A case file changes only what this guide and `cases/reference.R` list. If a result
seems to need a change to any of the above, that is a question for the PI, not an edit.
