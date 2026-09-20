# Option 3 — decoupling decision, per channel

Date: 2026-09-18 (§0a and §5 updated 2026-09-20)
Responds to: `HANDOVER_option3_database_findings_2026-09-18.md`
Supersedes: the `DECOUPLED_FORMULAS` sketch in `PROGRESS.md` Part 3
Status of my own claims: the database *reasoning* is verified (see §1); the
database *file* was not available to this session, so line numbers and species
spellings are taken on trust from the VS Code session's reading. Two items at
the end need a PHREEQC run and cannot be settled by argument.

---

## 0. The handover was right, and right to stop where it did

The VS Code session found that two of four draft decoupled formulas were
wrong, not merely unverified, and declined to patch `IOM_N`/`IOM_CO2` because
the fix needs a design decision. That was the correct call — substituting
`Ntg` for `Amm` would have been actively wrong (see §2), and it would have
looked like a fix.

It also correctly identified the general principle, which I verified
independently: in PHREEQC's `SOLUTION_MASTER_SPECIES`, the **element name**
(first column) defines the mass/electron-balance group. `C(-4)` and `C(+4)`
are redox *states* of the single element `C` and share one electron balance;
`Mtg`, `Sg`, `Ntg`, `Hdg`, `Oxg` are separate elements and do not. So writing
`-formula C(4) 1 ...` sets only the *starting* valence of added material and
PHREEQC's solver then moves it down the shared `C` ladder — exactly the
behaviour Option 3 exists to prevent. This is standard PHREEQC semantics, not
an inference about this particular database.

## 0a. Why `Mtg 1` works when `CH4` would not — read this with PROGRESS Part 1

**Added 2026-09-20**, because `PROGRESS.md` Part 1 was corrected and the two
documents now need to be read together without tripping the reader.

Part 1 now states (correctly, per the PHREEQC v3 manual) that `-formula CO2` is
**legal syntax** and that naming a molecule buys you its **stoichiometry, not
its identity** — `-formula CO2` adds 1 C + 2 O as elements, which then speciate
by equilibrium exactly as `C 1 O 2` would.

Someone who has just read that could reasonably ask why the decoupled formulas
in §3 — `Mtg 1 H 0.435`, `Sg 1 H -1` — are expected to behave any differently.
The answer is the distinction Part 1 draws, applied in the other direction:

- `CH4` is a **molecule name used as stoichiometric shorthand**. PHREEQC expands
  it to 1 C + 4 H, and that carbon joins the shared `C` mass balance, free to
  move anywhere on the coupled ladder.
- `Mtg` is **not a molecule name at all — it is a distinct element.** Writing
  `Mtg 1` adds one mole of element `Mtg`, which has its own mass balance and no
  electron exchange with `C`. It cannot leave that balance, so it persists as
  methane.

So the two documents are consistent, and the mechanism is the same fact seen
from both sides: identity survives only when it is carried by an *element*, not
by a formula. This is exactly why Option 3 requires database-level decoupled
species and cannot be achieved by writing product molecules into `-formula`.

---

## 1. Decision, per channel

The key realisation is that **the four channels are not in the same
situation**, so a uniform "decouple everything" design was wrong from the
start. What matters per element is: *does coupled equilibrium put the released
material where kinetics would have left it?*

| channel | decision | decoupled formula | why |
|---|---|---|---|
| `IOM_CH4` | **decouple** | `Mtg 1 H 0.435` | Correct as drafted. `Mtg` exists and is genuinely decoupled. |
| `IOM_S` | **decouple** | `Sg 1 H -1` (pending §4) | Species name corrected from `S(-2)`. Sulfide/sulfate equilibration is genuinely inhibited at low T. |
| `IOM_N` | **leave coupled** | *(none)* | See below — coupled equilibrium already gives the right answer. |
| `IOM_CO2` | **decouple, via `Hdg`** | *(no new species; see §3)* | Coupled equilibrium gives exactly the wrong answer, but the fix is not a carbon species. |
| `IOM_CHn` | follows `IOM_CO2` | *(same mechanism)* | Its carbon has the same problem; no separate solution needed. |

### Why `IOM_N` should stay coupled — this simplifies the design

The handover offered this as option (2) and was right to. Checking the
chemistry:

- `IOM_N` releases `N 1 H 2`, i.e. nitrogen at oxidation state **−2**
  (amide/amine-like).
- NH₃/NH₄⁺ is N(−3); N₂ is N(0); NO₃⁻ is N(+5).
- At Enceladus pore conditions (pe ≈ −11.5, pH ≈ 11.5, H₂-rich),
  **NH₃/NH₄⁺ is the thermodynamically stable nitrogen species.** N₂ and NOₓ
  are favoured only under oxidising conditions.

So coupled equilibrium moves released N from −2 to −3: a one-electron step, in
the same direction kinetics would take it, to the same species. **Decoupling N
would change almost nothing at low pe.** There is no artefact to fix.

This also disposes of the `Ntg` question: pinning released N as N₂ would
*oxidise* it from −2 to 0 — moving it away from the correct answer, not
toward it. `Ntg` is worse than leaving N coupled. **Do not add a decoupled
ammonium species.** It would be database work for no change in result.

(Caveat worth keeping: this argument depends on the pore fluid staying
strongly reducing. If an oxidising case is ever modelled, revisit — under
oxidising conditions coupled N *would* run away to N₂/NO₃⁻ and decoupling
would start to matter.)

---

## 2. Why `IOM_CO2` is the opposite case

Unlike nitrogen, coupled equilibrium does **not** leave released carbon where
kinetics would. From this project's own PHREEQC output (the
`k1.00e-07_d194200_orb12` run): at pe = −11.5,
**C(−4)/C(+4) = 3.1×10⁻¹⁶ / 7.9×10⁻³³** — i.e. essentially all carbon becomes
methane, seventeen orders of magnitude of it.

`IOM_CO2` releases carbon at about **+1.6** (only ~1.3 O per C). Coupled
equilibrium drives that to −4, a 5.6-electron reduction per mole, drawing
those electrons from H₂ and from oxidising Fe(II) minerals. That is the
artefact in `PROGRESS.md` Part 2b, and for carbon it is real and large.

So option (2) works for N but must **not** be extended to CO₂ carbon.

### And the option the handover didn't consider: `Hdg`

Adding a decoupled oxidised-carbon element is unattractive:

- Carbonate has three protonation states (CO₂/HCO₃⁻/CO₃²⁻), so it is a
  substantial database edit, not a one-line addition like `Sg`.
- **More seriously: decoupled carbon would no longer see calcite, dolomite or
  any other carbonate mineral.** Carbonate precipitation would silently stop
  working for exactly the carbon we care about — a worse artefact than the one
  being fixed, and a quiet one.

Routing `IOM_CO2` carbon into `Mtg` instead is also wrong: that imposes the
same 5.6-electron reduction by hand rather than letting the solver do it.

**The better option uses a species already in the database.** The mechanism
driving carbon to CH₄ is the large H₂ reservoir setting pe ≈ −11.5.
`Hdg` (decoupled H₂) already exists in `Core11_idealgas_mod_v4.dat`. If H₂ is
decoupled, it no longer participates in the electron balance, so it cannot act
as the electron donor that reduces released carbon — and pe is then set by the
mineral assemblage rather than by dissolved H₂.

This is attractive for three reasons:

1. **No new database species.** Uses `Hdg`, already present and presumably
   already exercised.
2. It targets the actual mechanism (H₂ as electron donor) rather than
   symptomatically pinning each product.
3. It directly addresses `PROGRESS.md` Part 2b's spurious-electron-demand
   finding, since that demand was being met by H₂ plus mineral oxidation.

**This is a hypothesis, not a verified result.** It predicts that decoupling
H₂ alone should largely stop the carbon runaway. It must be tested (§4) before
being relied on. It may also prove too blunt — decoupling H₂ affects *every*
redox couple in the system, including mineral reactions we may want left
coupled. If that turns out to be a problem, the fallback is the carbonate
database edit with its calcite caveat accepted and documented.

---

## 3. Revised decoupled-mode design

    coupled mode   : all channels as in Part 1 (elements; solver decides)
    decoupled mode : IOM_CH4 -> Mtg 1 H 0.435
                     IOM_S   -> Sg 1 H -1
                     IOM_N   -> unchanged (N 1 H 2), deliberately coupled
                     IOM_CO2 -> unchanged formula; decoupling achieved by
                                Hdg in the database, not by a carbon species
                     IOM_CHn -> unchanged formula; same mechanism as IOM_CO2

Three of five channels keep one formula for both modes. Only `IOM_CH4` and
`IOM_S` need a second formula column — so the config change is smaller than
Part 3 assumed, but `redox_mode` must still be **per-element**, because H₂
decoupling is a database-level switch while Mtg/Sg are formula-level.

### Element conservation — verified for both changed channels

Part 3's invariant (identical elements released in both modes) holds:

    IOM_CH4  coupled  C 1 H 4.435       H total = 4.435
             decoupled Mtg 1 H 0.435    Mtg carries 4 H + 0.435 free = 4.435  OK
    IOM_S    coupled  S 1 H 1           H total = 1
             decoupled Sg 1 H -1        H2Sg carries 2 H − 1 free  = 1        OK

Negative coefficients are legal: the PHREEQC manual (KINETICS, `-formula`)
states that "a negative stoichiometric coefficient and a positive value for
reaction progress gives a negative mole transfer, which removes reactants from
the aqueous solution," and its own `Fe_di_ox` example uses −1.0. Verified
against the manual text in the project, not from memory.

---

## 4. The two things that need a PHREEQC run, not more argument

Both are for Lucas / the VS Code session; neither can be settled by reading.

**(a) Does `Sg 1` implicitly demand H₂Sg's two hydrogens?** The database line
is `Sg  H2Sg  1  H2Sg  34.08`, so the element `Sg` has master species H₂Sg.
Writing `Sg 1` adds one mole of element Sg. Whether PHREEQC also pulls the two
H of H₂Sg from solution automatically, or expects the formula to supply them,
determines whether the correct form is `Sg 1 H -1`, `Sg 1 H 1`, or `Sg 1`
alone. Test: add a known amount via each variant to a simple solution and
check the resulting total H and total S against hand arithmetic.

**(b) Does decoupling H₂ (`Hdg`) actually stop the carbon runaway?** Test in
isolation, no ENTICES: a reducing solution with a known H₂ inventory, add
oxidised carbon, run once with H₂ as `H2` and once as `Hdg`, and compare the
resulting C(+4)/C(−4) split and pe. This is the cheap decisive test of §2's
hypothesis, and it should be done before any module or config changes.

---

## 5. Relationship to `PROGRESS.md` — status as of 2026-09-20

**This section previously listed four edits `PROGRESS.md` needed. All four have
been applied**, in the consolidating edit of 2026-09-18/20. Recorded here so
nobody repeats them:

- ~~Part 3's sketch table replaced with §3 above~~ — **done.**
- ~~Part 3 referred to `Core11_idealgas_mod_v2.dat`; actual file is v4~~ —
  **done.**
- ~~Part 4 showed A = 1e13 and the single-Ea values~~ — **done.** Part 4 now
  carries A = 2e15, the Vitrimat 2018 distributions, and the 250 C validation
  result. See `claude/PATCH_F_vitrimat2018.md`.
- ~~Part 5 listed the nitrogen discrepancy as unresolved~~ — **done.** N = 2.142
  confirmed correct; 3.284 traced to misreading the hydrogen gram figure (46 g
  H) as nitrogen. See `claude/PATCH_F_addendum_N_and_CHn.md`.

One further `PROGRESS.md` change, made 2026-09-20, **strengthens** this
document rather than conflicting with it: Part 1's first justification ("PHREEQC
mechanically cannot add a persisting molecule") was wrong and has been
corrected. Phase and molecule names *are* legal in `-formula`, and decoupled
pseudo-elements *do* persist. The old wording contradicted this document's
entire premise; the corrected wording is consistent with it. §0a above spells
out the distinction so the two documents read together cleanly.

**Nothing in §1–§4 changes as a result.** The per-channel decisions, the
element-conservation table, the negative-coefficient justification, and the two
outstanding PHREEQC tests are all unaffected — none of them ever depended on
the incorrect claim.
