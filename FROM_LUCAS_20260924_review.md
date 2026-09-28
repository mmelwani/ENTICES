# Review of "From Lucas 20260924" — current ENTICES code + a real production run

Date: 2026-09-24
Files reviewed: `ENTICES_LF_20260917.R`, `enceladus_k1.00e-07_d194200_T150_org15_orb12.pqi`
(+ `.out`), `coreclath_CH4coupled.dat`.

**Headline: the calibration in Lucas's current code and this run is three
patch-generations behind what's on the `iom-module` branch, and the coupled-
carbon-to-methane artefact that motivated Option 3 is visible in the actual
output, not just predicted.**

---

## 1. The IOM channels are not sourced from `iom_module.R` — they're hand-copied and stale

`ENTICES_LF_20260917.R` doesn't `source()` the module at all. It has its own
inline copy of the channel definitions (his own "PATCH C/D1/D2" labels, lines
1124-1189, duplicated again at 3952-4017 for the restart path):

```r
self$iom_pools <- data.frame(
  name      = c("IOM_CO2", "IOM_CH4", "IOM_N", "IOM_S"),
  m0_per_kg = c(4.169, 1.908, 0.875, 0.462),
  formula   = c("C 1 H 1 O 1.3115", "C 1 H 4.435", "N 1 H 2", "S 1 H 1"),
  Ea_J      = c(214e3, 228e3, 215e3, 215e3),
  logA      = c(13, 13, 13, 13),
  ...
)
```

Every one of these values is the **original PATCH_C** number, before any of
the corrections made on this branch:

| | this file | should be (branch `iom-module`, current) | correction date |
|---|---|---|---|
| A (all channels) | 1e13 | **2e15** | 2026-09-18 (Burnham 2019; 1e13 is the superseded 1989 value) |
| IOM_CO2 | single Ea, 214 kJ/mol | **7 sub-pools**, Vitrimat-2018 shape, shift +3.587 kcal | 2026-09-18 |
| IOM_CH4 | single Ea, 228 kJ/mol | **13 sub-pools**, shift +3.736 kcal | 2026-09-18 |
| IOM_N | Ea 215 kJ/mol | **242.8 kJ/mol** | 2026-09-17 (216→217.3), then 2026-09-20 (basis-bug fix →242.8) |
| IOM_S | Ea 215 kJ/mol (pinned to N) | pinned to corrected N, 242.8 kJ/mol | same |
| IOM_CHn | absent | present, 7 sub-pools | 2026-09-18 |

The file's own comment is honest about this: *"CAVEAT: single Ea per
channel... There is a low-Ea tail this misses, and for COLD Enceladus runs
that tail is the only thing that can release anything at all. Sub-pools to
follow."* — so Lucas already knew this was provisional. It just means he
wrote this before seeing the work done here since 09-17, not that anything
is broken on his end.

One thing worth telling him directly: **the 09-17 single-Ea fix already
found and fixed an 18% error in this exact CH4 value** (228 → 229.4 kJ/mol,
before being superseded again by the distributed treatment on 09-18) — so
this file's CH4 channel carries a bug that was caught and fixed twice over,
not once.

Also worth noting: he independently flagged the same nitrogen-composition
question this project spent several sessions resolving (*"NOTE: v3 comments
quoted N_3.284 (46 g N) — gram values used here; confirm against source
composition"*), and landed on the same correct answer (30 g N, i.e. N=2.142)
without yet having seen `PATCH_F_addendum_N_and_CHn.md`, which closes that
question with the same value and traces the wrong one to a misread of the
hydrogen gram figure.

## 2. The predicted artefact is visible in the actual output

This matters more than the calibration staleness on its own, because it's
not hypothetical: the coupled-carbon-collapses-to-methane problem
(`PROGRESS.md` Part 2, `OPTION3_decision.md` §2) shows up directly in this
run's porewater chemistry, right after the IOM channels react (150 °C
porewater, step 1):

    pH  = 11.746
    pe  = -11.692   (strongly reducing)
    CH4 = 4.327e-07 mol/kg
    CO2 = 2.551e-16 mol/kg

**CH4/CO2 ≈ 1.7×10⁹ — essentially all released organic carbon becomes
methane, none stays as CO2.** That is the exact signature `OPTION3_decision.md`
§2 used to justify the `Hdg`-decoupling design (its own reference point,
from an earlier run, was C(−4)/C(+4) ~10¹⁷). Cassini observes both CO2 and
CH4 in the plume; a model that only ever produces one of them at the
porewater temperature cannot reproduce that, regardless of how well the
release *rate* is calibrated. This run uses `coreclath_CH4coupled.dat` in
fully coupled mode — no `Hdg`, `Sg`, or `Mtg` in any `-formula` — so Option 3
isn't implemented in this run either, consistent with its filename and
Lucas's note that it's "the coupled CH4 database."

(The database itself *does* carry the decoupled species — `Hdg`, `Sg`,
`Mtg`, `Ntg` are all present and match `Core11_idealgas_mod_v4.dat`'s
definitions closely; it's a near-identical file, missing only the
fluorine/fluorapatite section MMD added to the `_v4` copy this project has
been testing against. So there's no database blocker to using the tested
`Sg 1 H 1` / `Mtg 1 H 0.435` formulas here once the decision to do so is
made.)

## 3. What Lucas got right, independently

Worth saying plainly, not just listing problems: his integration design
converged on the same pattern `iom_module.R`'s own INTEGRATION POINTS
section recommends, without having seen it —

- `iom_mass_kg` constructed the same way (`organic_wt_percent/100 *
  total_rock`);
- restart handling ("PATCH D1" punches `IOM_*_mol` columns so remaining
  moles can be recovered; "PATCH D2" rebuilds a config with `m0_per_kg`
  replaced by recovered absolute moles and reads it back in) is exactly the
  restart mechanism the module's docstring describes as the right approach;
- the RATES/KINETICS text his code generates is structurally identical to
  `iom_rates_block()`/`iom_kinetics_block()`'s output, just built inline
  instead of by calling them.

That's a good sign for actually wiring the module in: the surrounding
scaffolding (mass calc, restart, punch headings) doesn't need to change,
only the channel table needs to become `iom_default_config()`'s output
instead of a hand-copied one four generations old.

## 4. Recommendation

1. **Give Lucas the current `iom_module.R`** (or point him at the
   `iom-module` branch once pushed) so his next run uses the distributed-Ea,
   corrected-N calibration. His own code structure needs almost no change —
   swap the hardcoded `self$iom_pools` assignment for
   `self$iom_pools <- iom_default_config()`, called once and cached exactly
   as his `is.null(self$iom_pools)` guard already does.
2. **Decide on Option 3 before or alongside that.** The artefact is real and
   visible in his own output now, not just predicted. Since `Hdg`/`Sg` are
   already in his database, this doesn't need new database work — it needs
   the `redox_mode` decision (`PROGRESS.md` Part 5) and passing
   `iom_default_config(c(C="decoupled", S="decoupled"))`'s formulas through,
   plus entering H2 as `Hdg` in the SOLUTION blocks where organics are
   active.
3. Worth flagging to Lucas separately, since it's outside IOM entirely: the
   final mixed-ocean block in this output (after 46 mixing steps) shows
   C/N/S all below 1e-11 mol/kg — essentially none of the porewater's
   organic-derived chemistry survives to the ocean box in this run. Whether
   that's expected dilution/mixing behaviour or worth a second look is his
   call, not something this review investigated further.
