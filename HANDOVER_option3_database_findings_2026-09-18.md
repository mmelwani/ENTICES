# Handover: Option 3 decoupled-species names are wrong in the current draft

Date: 2026-09-18
For: the Claude Chat session that authored `PROGRESS.md` Part 3 and the
`DECOUPLED_FORMULAS` sketch it contains (via `HANDOVER_distributed_Ea_and_option3.md`
Task 4), relayed through the user.
Trigger: the user added the actual project database,
`Core11_idealgas_mod_v4.dat`, to the repo and asked me to check its format
against the decoupled-species assumptions in `iom_module.R`. This is a
correction to that sketch, found by reading the real database text -- not
by running PHREEQC (not this session's job; per the user, that's Lucas's
downstream task, ours is just producing correct input text).

## Summary

`PROGRESS.md` Part 3 assumed decoupled ammonium ("Amm") and decoupled
sulfide ("S(-2)") exist in the project database the same way stock
`phreeqc.dat` provides `Mtg`/`Amm`. **`Amm` does not exist in
`Core11_idealgas_mod_v4.dat` at all, and `S(-2)` is not a decoupled
species -- it's just a valence tag inside the normal coupled sulfur
ladder.** I built `DECOUPLED_FORMULAS` in `iom_module.R` directly from that
sketch without checking it against a real database (none was available to
this session at the time), and it needs correcting.

## What the database actually has

`Core11_idealgas_mod_v4.dat`, lines 336-342, under a section explicitly
labeled `# Redox-uncoupled gases from phreeqc.dat`:

```
Hdg	Hdg	0	Hdg		2.016   # H2 gas
Oxg	Oxg	0	Oxg		32      # O2 gas
Mtg	Mtg	0	Mtg		16.032  # CH4 gas
Sg	H2Sg	1	H2Sg		34.08   # H2S gas
Ntg	Ntg	0	Ntg		28.0134 # N2 gas
```

Five decoupled pseudo-elements exist: H2 (`Hdg`), O2 (`Oxg`), CH4 (`Mtg`),
**H2S (`Sg`, species name `H2Sg`)**, and **N2 (`Ntg`)**. `Sg` has its own
small internal acid-base ladder (`H2Sg = HSg- + H+`, added by the user --
comment says "MMD 2024/07/07" -- so this decoupled-sulfide mechanism is the
user's own prior, already-tested work, not something added for this
project). There is no decoupled ammonium species anywhere in the file, and
no decoupled "CO2-only" species either.

Compare to the normal COUPLED redox ladders in the same file (lines
258-263, 303-307, 315-323):

```
C(-4)  CH4     ...   <- same "C" mass balance as:
C      HCO3-   ...
C(+4)  HCO3-   ...

N(-3)  NH3     ...   <- same "N" mass balance as:
N(0)   N2      ...
N(+3)  NO2-    ...
N(+5)  NO3-    ...

S(-2)  HS-     ...   <- same "S" mass balance as:
S(+6)  SO4-2   ...
```

`C(-4)`, `N(-3)`/`N(0)`, and `S(-2)` are valence STATES of one shared
element (`C`, `N`, `S` respectively) -- writing `-formula C(4) 1 ...` or
`-formula S(-2) 1 ...` only sets the STARTING valence of newly-added
material; PHREEQC's normal equilibrium will still move it wherever
thermodynamics favors, through the SAME redox ladder as everything else
using that element. That is precisely the coupled behavior Option 3 is
trying to avoid. Only `Mtg`, `Sg`, `Ntg`, `Hdg`, `Oxg` are genuinely
independent mass balances that do not exchange electrons with the rest of
the system.

## Consequence for each channel's draft decoupled formula

| Channel | Draft (wrong) | Problem |
|---|---|---|
| `IOM_CH4` | `Mtg 1 H 0.435` | **Correct as drafted.** `Mtg` is real, decoupled, and molar mass matches CH4 (16.032). Still need to verify the H-balance (does PHREEQC accept a bare element line `H 0.435` alongside a pseudo-element in the same `-formula`? probably yes, but untested). |
| `IOM_S` | `S(-2) 1 H 1` | **Wrong species name.** Should be `Sg`, not `S(-2)` -- `S(-2)` doesn't decouple anything. `Sg` = H2S (2 H per S), but the coupled formula only releases 1 H per S (`S 1 H 1`), so the corrected draft is probably `Sg 1 H -1` (borrowing 1 H from solution) -- same kind of H-balance question as `Mtg`, unverified. |
| `IOM_N` | `Amm 1` | **No such species exists.** The only decoupled nitrogen option is `Ntg` (N2 gas) -- the opposite end of the redox ladder from what IOM_N releases (reduced, NH2-like nitrogen). Using `Ntg` would pin released N as N2, not as NH4+/NH3, which is not obviously what "decouple N so it doesn't equilibrate away from ammonium" was meant to achieve. This needs a real decision: either (a) add a genuine decoupled-ammonium pseudo-element to the database, analogous to how `Sg` was added for sulfide, or (b) reconsider what "decoupling nitrogen" should mean given only `Ntg` is available. I have not guessed at an answer -- this is a database-design choice, not something to infer from the existing file. |
| `IOM_CO2` | `C(4) 1 H 1 O 1.3115` | **Doesn't decouple anything.** No decoupled oxidized-carbon species exists at all (`Mtg` only covers the reduced/methane end). As drafted, "decoupled" `IOM_CO2` would still fully equilibrate through the coupled `C` ladder and could still reduce to `C(-4)`/CH4 -- which is exactly the artifact Option 3 exists to prevent. Also needs a real decision: either add a decoupled oxidized-carbon species, or reconsider whether `IOM_CO2`'s released carbon should route into `Mtg` too in decoupled mode (treating "decoupled" as "carbon stays wherever it's released, full stop," rather than "oxidized and reduced carbon are each pinned separately"). |
| `IOM_CHn` | (none defined) | Already flagged as an open question in `iom_module.R`'s docstring for a different reason (H/C=1.8 doesn't map onto `Mtg`'s stoichiometry). Given `IOM_CO2` also has no clean home, CHn is not the only carbon channel with this problem. |

## What I did and didn't change

I have **not** changed `iom_module.R`'s `DECOUPLED_FORMULAS` beyond adding a
pointer to this document, because the correct fix for `IOM_N` and `IOM_CO2`
requires a design decision (possibly a real database edit, adding new
decoupled species) that isn't mine to make. Silently substituting `Ntg` for
`Amm` or leaving `C(4)` in place would both be guesses dressed up as fixes.
`IOM_CH4` (Mtg) appears genuinely fine as drafted. `IOM_S` has an easy,
low-risk correction (`S(-2)` -> `Sg`) that I'm flagging here rather than
silently applying, since the H-balance (`H -1` vs `H 1`) still needs
checking and this and `IOM_N`/`IOM_CO2` are naturally reviewed together.

## Suggested next step

Whoever owns the Option 3 database design (this session, the next VS Code
session, or directly with the user) should decide, before any more code
changes:
1. Does the project want to ADD decoupled ammonium and decoupled
   oxidized-carbon species to `Core11_idealgas_mod_v4.dat` (following the
   `Sg` precedent -- the user already did this once for sulfide in 2024)?
2. Or does "decoupled mode" for N and CO2 mean something narrower than a
   1:1 species swap -- e.g., only CH4 and H2S get pinned (since those are
   the reduced END-products a coupled model would erroneously produce),
   while N and CO2-channel carbon stay on the coupled ladder because their
   OWN equilibrium behavior (NH3 dominant when reducing; CO2/HCO3- dominant
   unless pe is extreme) is close enough to correct without decoupling?
   That would be a real, defensible simplification of PROGRESS.md Part 3's
   original four-species design, not a bug -- but it's a decision, not
   something I should assume.

Either way, `iom_module.R`'s `DECOUPLED_FORMULAS`, its docstring, and
`PROGRESS.md` Part 3's sketch table all need to be revisited together once
that's settled, so they don't drift out of sync with each other again.
