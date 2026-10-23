# Review of Lucas's 2026-10-19 commit, reply to Claude Code, cleanup task, and the reconstruction module

Date: 2026-10-21
From: the claude.ai review session.
For: Claude Code (VS Code) and the PI. Lucas should see §1.
Read:
- the local checkout at `rock-composition` @ `2d55fbb` (fetched 2026-10-21);
- `~/Documents/Enceladus_hydrotidal/`;
- the PHREEQC v3 manual.

I can't run PHREEQC or R here. Numbers marked *(hand estimate)* are mine;
everything else is read from the files.

---

## 1. Lucas's commit `2d55fbb` ("No-precipitation guard on tephroite added")

The message understates the commit: it is +998/−269 lines across
`ENTICES_LF_20260917.R` and `ENTICES_runs.Rmd`. What it actually does:

| change | assessment |
|---|---|
| Tephroite `IF SR > 1 THEN rate = 0` guard | ✓ Fixes review item §5. |
| **Magnetite added to the porewater list (EP1) in all generators that lacked it** | ✓ **Claude Code's "4 generators missing" is now resolved.** All six have it. (The restart generators get it through their `secondary_minerals` vectors.) |
| `rec_eqb_minerals` — auto-promote fast primaries to `dissolve_only` equilibrium, using the half-life at pH 11 × margin 2 vs one orbit | Sound idea, close to review §4's per-mineral switch. Problems listed below. |
| `bulk_equil` gains a `porewater` mode (water/rock ratio, profile T) | One bug, listed below. |
| RATES blocks built from a list, only for kinetic minerals | Good refactor. The restart reads promoted minerals back from `*_remain_g`. |

### Problems to fix

**(a) `rec_eqb_minerals` only exists in `bulk_kinetic`.** The generator that
hung was `tidal_kinetic` (the 10-06 diagnostic). As committed, it doesn't touch
the convergence problem.

**(b) The pH-11 reference is not conservative for every mineral.**
- For acid-branch minerals (forsterite, enstatite, diopside, anorthite,
  tephroite, troilite) pH 11 gives the *slowest* rate, so it under-promotes.
  That is safe.
- For albite and K-feldspar (base branch) pH 11 gives the *fastest* rate.
- If the real pH stays at 9, albite is ~14× slower and K-feldspar ~44× slower
  than assumed, so the ×2 margin doesn't cover it.

The code comment says the margin guards against "the pH 11 reference being an
underestimate of the real in-run pH excursion". For the base-branch minerals,
the risk runs the other way. In the 160 °C diagnostic pH did reach 11.2–11.5,
so the result happens to be right there. Fix the comment, and either:
- evaluate base-branch minerals at a lower reference pH; or
- (better) drop the fixed pH and use the spin-up (§2.4).

**(c) Which minerals get promoted depends on unsourced surface areas.** SSA is
0.1 m²/g for forsterite, enstatite, diopside and tephroite, and 5 m²/g for
albite, K-feldspar, anorthite and troilite. That 50× split decides the outcome.

*(Script-checked estimate from the RATES constants, 160 °C, pH 11, orbit
118,800 s):*

| mineral | t½ | promoted? |
|---|---|---|
| albite | ≈ 140 s | yes |
| K-feldspar | ≈ 5,100 s | yes |
| anorthite | ≈ 7.0×10⁴ s | no |
| forsterite | ≈ 1.0×10⁵ s | no |
| enstatite, diopside | ~10⁷ s | no |

Anorthite and forsterite are both within a factor of 2 of the threshold.
Forsterite is the borderline case and the one Claude Code found still slow
after the feldspars were gone. RECONSTRUCTION_STATUS §8 item 5 already lists
surface areas as unsourced for every kinetic phase. **Record in the CHANGELOG
that the promotion outcome is driven by SSA choices**, so nobody reads it as a
kinetic result.

**(d) The tephroite rate law is still summed.** It is still
`k = (kacid + kneut)` (still doubled at pH 7), while
`mineral_dissolution_constants.r` says "pH = 7 value for all pH > 7". Use the
acid branch evaluated at `max(ACT("H+"), 1e-7)`, as agreed on 10-06. At pH 11
the two forms agree, so the half-life check isn't affected.

**(e) The tephroite pH-dependent Ea is a reasonable fit, but note what it
implies.** The formula is `Ea = −60.11·ln(pH) + 106.17` kJ/mol, floored at 0.
*(Hand check)* it reproduces Casey's own quoted temperature sensitivity:

| pH | 10 K rise, Casey | 10 K rise, this fit |
|---|---|---|
| 2 | ~2.2× | 2.3× |
| 3–4 | ~1.5× | 1.5× |
| 5 | ~1.1× | 1.1× |

(PHREEQC Basic `LOG` is natural log, consistent with R's `log`.)

But Ea reaches 0 above pH ≈ 5.85. So at Enceladus pH (9–11.5) tephroite has
**no** temperature dependence, extrapolated from 25–45 °C data to 160 °C. It is
low stakes (0.5 wt%, Mn only), but say so in the CHANGELOG. Update
RECONSTRUCTION_STATUS §4, which still says "Ea UNRESOLVED".

**(f) Bug in `bulk_equil` porewater mode.** It computes `pore_pressure_atm`
from the profile, but SOLUTION 1 still uses `pressure
{round(ocean_pressure_atm, 1)}` (line ~2604). In porewater mode it should use
`pore_pressure_atm`.

**(g) Cosmetic, in `tidal_kinetic_restart`.** `Magnetite` is now in
`secondary_minerals` *and* still appended after Lizardite in the
SELECTED_OUTPUT list (`.all_eq_so`), so it appears twice. That likely gives a
duplicate CSV column (R renames it `Magnetite.1`). Remove the `append`.

---

## 2. Reply to Claude Code's response to REVIEW_convergence

### 2.1 Peng–Robinson: the manual answers the open question

> "If the critical temperature and pressure are not defined, the ideal gas law
> will be used. **Ideal gases and Peng-Robinson gases cannot be mixed in a
> GAS_PHASE.**" (GAS_PHASE, Description of Data Input)

So a fixed-pressure GAS_PHASE containing `Mtg(g)` (PR constants active) and
`CH4(g)`/`H2(g)` (ideal) is not allowed. Below 150 °C, where decoupled methane
is the default, that is exactly the mix we'd get.

**Recommendation:** comment out `-T_c`/`-P_c`/`-Omega` on the five decoupled
gases (`Oxg(g)`, `Hdg(g)`, `Ntg(g)`, `Mtg(g)`, `H2Sg(g)`) so the whole
database is ideal. That matches its name (`idealgas_mod`) and how every
ordinary gas is already treated. Fix the `Ntg(g)` `T_c` typo in the same edit.
Do it in **all three** `.dat` files, or better, in the single canonical
database (§3).

### 2.2 Magnetite

Resolved by Lucas in `2d55fbb` (§1). Thanks for catching the scope; my 10-06
review checked only one generator.

### 2.3 "Ideal-gas database" correction

Accepted. My statement was half right, as Claude Code says.

### 2.4 Spin-up won't hand over a quiet state

Accepted. Forsterite (1.37 mol, t½ ~ 1 day at 160 °C by my estimate) is the
long tail. The spin-up needs a **rate-based stopping rule**: stop when every
remaining kinetic primary's fractional rate is below, say, 1% per orbit. The
gas fix and pruning the `0 0` list still come first, as both sessions agree.

### 2.5 Go-ahead for what Claude Code can do now

Both proposals are fine:
- **Magnetite:** already done by Lucas, so no patch is needed.
- **Fixed-pressure GAS_PHASE test:** run it at **70 atm, bubbles left in
  place**. Do the database PR edit in §2.1 first, or test only above 150 °C
  (coupled CH₄) so no PR gas is involved.

### 2.6 Decisions still owed by the PI

Unchanged from 10-06 §6: porewater pressure, bubble fate, ocean-side gas as
plume loss, spin-up, tephroite floor.

One input for the pressure decision: Lucas's new `bulk_equil` code already
computes `mean_layer_pressure_MPa(profile)`. Using the same function for the
tidal porewater would replace the hand-entered `66.7#261.2` with a
profile-derived value.

---

## 3. Cleanup task for Claude Code

**Do this after PR #3 merges.** Lucas's open PR touches the big R file and the
Rmd, so moving files now would make his merge conflict.

Work on a new branch `cleanup` and open it as a PR. Use `git mv` / `git rm`
(history keeps everything). **Delete nothing outside git without asking the
PI.** Present the plan as a table (file → keep / move / archive / remove, and
why) and get approval before executing.

### What I found (inventory as of `2d55fbb`)

- **50 tracked files, all flat at the top level** apart from three folders.
- **Seven review-exchange docs are untracked:** HANDOVER_for_review, REVIEW_*,
  REPLY_*, RESPONSE_* (09-30) and CONVERGENCE_DIAGNOSTIC (10-06).

### Specific items

1. **The tracked `Primordial_CI_chondrite_mineralogy.Rmd` is the discredited
   v0 reconstruction.** It contains the hydration-expansion factors (1.35),
   `mag_factor 0.8` and the density-ratio approach. RECONSTRUCTION_STATUS §6
   says "Do not mine v0 for numbers".
   - The current version lives in `~/Documents/Enceladus_hydrotidal/`.
   - The working copy also has an **uncommitted, half-finished edit** of the
     v0 file. It contains two incomplete R assignments (`orthoclase_density`
     and `mol_mass_or <- ` with no value), so the file won't run.
   - **Ask the PI before discarding it.** The new density citations in it may
     be worth keeping as notes.
   - Replace the file as part of §4; don't just tidy it.
2. **Exact duplicate:** `HANDOVER_IOM_validation_2026-09-17.txt` is
   byte-identical to the `.md`. Remove the `.txt`.
3. **`From Lucas 20260924/` (untracked):** all four files are byte-identical to
   blobs already in git history (`b67926f`). It is safe to delete locally;
   confirm with the PI.
4. **`.gitignore` lists 10 specific `.pqi` filenames.** Replace them with
   patterns:
   - `phreeqc_inputs/*.pqi`, `*.out`, `*.csv` under run directories;
   - `.Rhistory`, `.RData`, `*.nb.html`, `.Rproj.user/`.

   Keep an explicit `!` exception for any example run that is deliberately
   committed.
5. **Large generated or binary files in git:**
   - `PHREEQC_Input_and_Examples_v3.pdf` (4.3 MB): replace with a link to the
     USGS page in the README.
   - The example run's `.pqi.out` (2.0 MB).
   - `phreeqc_test_output/*.out` (~400 kB).

   Keep the test outputs. They are the evidence behind
   PHREEQC_TEST_RESULTS / B3, so move them into `tests/` (below) rather than
   deleting them.
6. **Merged branches:** `iom-module` and `lucas-integration` are fully merged
   into `main`. Delete them locally and on GitHub, once the PI confirms.
7. **Three diverging copies of the database:**
   - `Core11_idealgas_mod_v4.dat` (repo);
   - `Core11_idealgas_mod_v4.dat` (Enceladus_hydrotidal; differs only by a
     tephroite RATES block);
   - `coreclath_CH4coupled.dat` and `coreclath_CH4uncoupled.dat`.

   Propose **one canonical source plus the two variants it generates**, and a
   note on which generator uses which. Don't merge them yet: Lucas is mid-merge
   on the clathrate/fluorine split (09-30 handover §5.1), so ask him first.
8. **`ENTICES_LF_20260917.R`:** a date-stamped name for the living main file
   confuses everyone (three copies with similar names exist across the
   machine). Propose renaming it to `ENTICES.R`, with `ENTICES_runs.Rmd`'s
   `source()` updated. **This is Lucas's file, so it needs his agreement.**
9. **`ENTICES_v3.R` is still on `main`** and carries the REACTION dosing bug
   (35.6× under-dose; see REPLY_to_response addendum). Move it to `legacy/`
   with a README stating the bug, or remove it (history keeps it).
10. **README.md is 73 bytes.** Write a real one covering:
    - what ENTICES is;
    - the directory layout;
    - which file is the entry point;
    - which database to use when;
    - how to run the self-tests;
    - where the design docs live.

### Proposed layout

    ENTICES.R  ENTICES_runs.Rmd  iom_module.R  mineral_dissolution_constants.r   (code, top level)
    reconstruction/        <- §4
    databases/             <- .dat files (+ README: which is canonical, which generator uses which)
    tests/                 <- iom_selftest.R, test_*.pqi, tests/output/*.out
    examples/              <- the from_lucas_20260924 example run, redox_budget inputs
    docs/                  <- PROGRESS.md, OPTION3_decision.md, PATCH_F_*.md, CHANGELOG.md, PHREEQC_TEST_*.md
    docs/archive/          <- superseded handovers, NEXT_STEPS, INTEGRATION_NOTE, the 09-30/10-06 review exchange
    legacy/                <- ENTICES_v3.R (if kept)

**Paths that move need their references updated.** For example:
- `iom_module.R` is `source()`d by relative path;
- tests reference `.dat` paths;
- docs cross-reference each other by filename.

Grep for every moved filename before committing, and run `iom_selftest.R`
afterwards.

### Not part of the cleanup, but owed

`PROGRESS.md` stops at 09-20. The 09-30 and 10-06 exchanges changed several
conclusions:
- `Hdg` is a no-op in production;
- the IOM composition is CR;
- the gas-phase artefact;
- the dosing bug.

Archiving those exchange docs without folding them into PROGRESS would lose
them. The review session will write that update; Claude Code shouldn't
archive the exchange docs until it lands.

---

## 4. Should the CI reconstruction live in ENTICES? Yes — as a separate module in the same repo

### Why the same repo

**ENTICES consumes the reconstruction's output (`ROCK_MOL_PER_KG`), and the
hand transfer between the two has already failed three ways:**
1. The repo carries the discredited v0 Rmd under the canonical name.
2. Lucas's 09-17 coefficients came from a 2-decimal xlsx, renormalised
   without magnesiochromite, with lawrencite's Fe not subtracted (fixed in
   `7e14d6a`).
3. Claude Code couldn't check RECONSTRUCTION_STATUS because it wasn't in the
   repo.

### Why a separate module

The same reasoning that split out `iom_module.R`:
- different update cadence;
- independently testable;
- ENTICES should depend only on its output.

### Proposed structure, `reconstruction/`

- `Primordial_CI_chondrite_mineralogy.Rmd` (the 2026-09-16 version)
- `Primordial_CI_chondrite_mineralogy_montecarlo.Rmd`
- `RECONSTRUCTION_STATUS.md`
- `tephroite_kinetics.txt` and `tephroite_handoff.txt` (or fold them into
  STATUS)
- **a new machine-readable output:** `reconstruction/rock_mol_per_kg.csv` with
  columns phase, formula, wt%, mol/kg and role
  (kinetic/equilibrium/inert), written by the Rmd.

### Wiring into ENTICES

ENTICES should either:
- read `ROCK_MOL_PER_KG` from that CSV; or, at minimum,
- **assert in a self-test that its hard-coded table matches the CSV** (after
  the lawrencite Fe subtraction, which is ENTICES's own modelling step).

That turns the hand transfer into a checked one.

### Leave behind

- `.RData` and `.Rhistory` (workspace state);
- `.nb.html` (885 kB, regenerated from the Rmd);
- `Superseded/` (v0 material; keep it out of the repo entirely, or put it in
  `reconstruction/archive/` with a DO-NOT-USE README if the PI wants the
  record);
- the duplicate `Core11` (§3 item 7).

### Separate repo instead?

Only worth it if the reconstruction will be published or reused on its own,
for example for other bodies. It could be split out later with
`git subtree split` without losing history. Same repo for now.

### Before moving anything

The reconstruction's own open items should travel with it, still marked
open:
- §8 item 8, the King et al. forward check (which the spin-up would provide);
- surface areas (§8 item 5);
- the schreibersite rate (§8 item 1b);
- the tephroite Ea (now fitted by Lucas; update §4).
