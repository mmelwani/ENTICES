# legacy/

`ENTICES_v3.R` is an earlier version of the ENTICES generator, kept for
reference only. **Do not use it for runs.** Use `ENTICES_LF_20260917.R` at the
repo root.

It has a known organic dosing bug. Its `REACTION` block multiplies each IOM
coefficient by `iom_C_moles`, so the delivered dose is ×`iom_C_moles` of the
intended one: about 36× too small (1/0.02809 = 35.6) in the redox budget
diagnostic case. This is why IOM moved to `KINETICS` in the current code. See
`docs/reviews/REPLY_to_response_2026-09-30.md` §D.3.

## Core11_idealgas_mod_v4.dat

Mohit's base PHREEQC database, retired 2026-10-23. **Superseded by
`databases/coreclath_CH4coupled.dat`**, which contains every Core11 entry (same
species, phases and RATES) plus the CH4/CO2/H2S hydrate phases. Kept for the
record; the archived docs that cite it describe runs made with this file.
