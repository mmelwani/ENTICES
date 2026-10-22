# legacy/

`ENTICES_v3.R` is an earlier version of the ENTICES generator, kept for
reference only. **Do not use it for runs.** Use `ENTICES_LF_20260917.R` at the
repo root.

It has a known organic dosing bug. Its `REACTION` block multiplies each IOM
coefficient by `iom_C_moles`, so the delivered dose is ×`iom_C_moles` of the
intended one: about 36× too small (1/0.02809 = 35.6) in the redox budget
diagnostic case. This is why IOM moved to `KINETICS` in the current code. See
`docs/reviews/REPLY_to_response_2026-09-30.md` §D.3.
