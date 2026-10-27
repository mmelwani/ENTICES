# profiles/

Interior-structure profiles that set the generators' temperatures and pressures.

## `EnceladusProfile_Seawater_10.0ppt_Tb272.4578K_PorousRock.txt`

A PlanetProfile output (Styczinski et al. 2023) for Enceladus. The file's
header records the run settings:
- ocean and pore fluid: Seawater, 10 g/kg;
- ice-shell base temperature Tb = 272.4578 K;
- porous silicate layer, maximum porosity 0.32 (`phiRockMax`), silicate EOS
  `Comet_67P_CG_v7_excluding_fluid_properties.tab`, no iron core;
- moment of inertia C/MR² = 0.335 ± 0.001.

Columns after the 81 header lines: P (MPa), T (K), r (m), phase ID, density,
and the other PlanetProfile outputs. The generators read phase 50 (porous
rock) only:
- `mean_layer_pressure_MPa()`: mean pressure over the permeable layer,
  26.47 MPa = **261.2 atm** for the full-core layer (`layer_thickness` 194.2 km).
  This is the porewater pressure.
- `ocean_floor_pressure_MPa()`: the pressure at the top of phase 50,
  6.76 MPa = **66.7 atm**. This is the ocean pressure.
- `mean_layer_temp_C()`: mean layer temperature, used when
  `porewater_temp = "profile"`.

The header also gives a mean silicate density `rhoSil_kgm3 = 2353.19`, which
is the core bulk density (porosity and pore fluid included) of this
structure model.

Provenance: this is the profile Lucas Fifer's production runs used (their
inputs carry the same 66.7 / 261.2 atm). Which PlanetProfile version and
input file produced it is not recorded here [citation needed: PlanetProfile
version and run configuration].

Reference: Styczinski, M.J., Vance, S.D., Melwani Daswani, M. (2023).
PlanetProfile: Self-consistent interior structure modeling for ocean worlds
and rocky dwarf planets in Python. *Earth and Space Science*.
[citation needed: volume, article number and DOI, not verified here]
