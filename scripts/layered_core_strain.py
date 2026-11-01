"""Tidal strain of Enceladus's core from a layered forward model (PI item N, 2026-10-31).

Structure: Bagheri et al. (2025, PSJ 6:245): ice shell, ocean and viscoelastic core, homogeneous
layers, compressible (bulk moduli shell 10 GPa and ocean 2.1 GPa: their Table 2), Andrade
rheology for shell and core (alpha = 0.2, beta = mu^(alpha-1) eta^(-alpha): their Appendix A), an
inviscid ocean (static, compressible), eccentricity tide e = 0.0047, n = 5.307e-5 rad/s (their Table 1),
R = 251.99 km, M = 1.08101e20 kg. Parameters are sampled log-uniformly (their Table 2 note) within the
Table 5 ranges of each scenario; the core density is set by the total mass (samples outside the Table 5
core-density range are rejected).

Core bulk modulus (PI 2026-11-01): not Bagheri's fixed 100 GPa (an intact-rock value the Love numbers do
not constrain) but the drained frame modulus of ENTICES's critical-porosity model, K_d = K_grain (1 - phi/phi_c),
with K_grain the Hill average of the reconstruction's minerals and the IOM (data/mineral_elastic_moduli.csv;
grain_elastic_moduli() in ENTICES_LF_20260917.R; 15 wt% IOM, 33% of the solid volume): 48.99 GPa, so K_d =
9.798 GPa at phi = 0.32, phi_c = 0.4.
--k-core sets it (Pa); it is written to both tables (K_core_Pa). The volumetric strain scales as 1/K_core
(checked here over 2-100 GPa, to 0.2%), which is how the R code follows K_d(phi(t)).

The core's share of the tidal heat. Bagheri et al. (2025) Eq. 8 (and B4), as printed,
Im(k2_core) = 4 pi G rho_core R_core^4 / (5 R_E^3) Im(Delta R_core)^4, is not dimensionally consistent (the
right side has units of m^4 kg s^-2 / m^... times a displacement to the fourth power, not a pure number), so
it is not used. Instead the time-averaged energy flux through a sphere of radius r,
F(r) = r^2 Im[conj(y1) y2 + l(l+1) conj(y3) y4 + conj(y5) y6 / (4 pi G)] (times a constant common to all r),
gives the power dissipated inside r: it equals -5 R Im(k2) / (4 pi G), i.e. their Eq. 5 total, at the surface
(to 1e-15), is continuous across the inviscid ocean (to 1e-13), and vanishes for an elastic core (1e-4 of the
total at eta_core = 1e25, Andrade's residual). The core's share is F(R_core) / F(R_E).

Forward model: TidalPy's radial solver (J. P. Renaud, TidalPy 0.8.1; shooting method, Takeuchi & Saito
1972 starting conditions). Radial functions y1..y6 are per unit tidal potential at the surface
(h2 = y1(R) g(R), 1 + k2 = y5(R)). The physical fields at the pole are those times U_pole = (3/2)(n R)^2 e
(Bagheri Eq. 1 with P20 = 1, P22 = 0).

Outputs:
  data/layered_core_strain.csv          per scenario, core viscosity and depth below the core surface:
                                        mean and sd over the samples of the pole amplitude of the
                                        volumetric strain (dilatation) div u = dy1/dr + (2 y1 - 6 y3)/r.
  data/layered_core_strain_summary.csv  per scenario and core viscosity: mean and sd of the core-surface
                                        radial displacement, the ocean-top radial displacement, the
                                        seafloor normal-stress amplitude |y2(Rc)| U_pole, the Love
                                        numbers and the tidal heat.
  stdout                                the Table 3 validation and the sensitivities.

Run: python -I scripts/layered_core_strain.py [out_dir] [--k-core Pa]  (needs numpy and TidalPy; ~40 s).
"""
import math, sys, os
import numpy as np
from TidalPy.RadialSolver import radial_solver

G = 6.6743e-11
RE = 251.99e3; M = 1.08101e20; W = 5.307e-5; ECC = 0.0047
U_POLE = 1.5 * (W * RE) ** 2 * ECC
K_SHELL, K_OCEAN = 10e9, 2.1e9
K_CORE = 9.798075e9   # K_d at phi 0.32 with 15 wt% IOM (see above); --k-core overrides
ANDRADE_ALPHA = 0.2
N_CORE, N_OTHER = 400, 60
DEPTHS_KM = np.arange(0, 181)

SCENARIOS = {
    # Bagheri et al. (2025) Table 5
    "dissipative_core": dict(D=(22e3, 34e3), eta_s=(6.3e14, 1e20), mu_s=(2e9, 4e9), rho_s=(700, 950), rho_o=(920, 1250),
                             Rc=(183e3, 213e3), mu_c=(4e9, 70e9), rho_c=(2160, 2400), eta_c=(4.4e11, 2.3e12)),
    "dissipative_shell": dict(D=(22e3, 34e3), eta_s=(2.5e13, 1.3e14), mu_s=(2e9, 4e9), rho_s=(700, 950), rho_o=(920, 1250),
                              Rc=(183e3, 213e3), mu_c=(5e9, 80e9), rho_c=(2160, 2400), eta_c=(3.16e13, 1e22)),
}
# Core viscosities of the strain table: the PI's sweep ranges (4.4e11-2.3e12; 3e13-1e22)
ETA_GRID = {"dissipative_core": np.geomspace(4.4e11, 2.3e12, 8),
            "dissipative_shell": np.array([3e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19, 1e20, 1e21, 1e22])}


def andrade(mu, eta, a=ANDRADE_ALPHA):
    beta = mu ** (a - 1) * eta ** (-a)
    return 1 / (1 / mu + beta * math.gamma(1 + a) * (1j * W) ** (-a) - 1j / (eta * W))


def solve(p, k_core):
    Ro = RE - p["D"]
    layers = [(0.0, p["Rc"], p["rho_c"], k_core, andrade(p["mu_c"], p["eta_c"]), N_CORE),
              (p["Rc"], Ro, p["rho_o"], K_OCEAN, 0j, N_OTHER),
              (Ro, RE, p["rho_s"], K_SHELL, andrade(p["mu_s"], p["eta_s"]), N_OTHER)]
    r = np.concatenate([np.linspace(lo, hi, n) for lo, hi, _, _, _, n in layers])
    rho = np.concatenate([np.full(n, float(d)) for _, _, d, _, _, n in layers])
    K = np.concatenate([np.full(n, k + 0j) for _, _, _, k, _, n in layers])
    mu = np.concatenate([np.full(n, m) for _, _, _, _, m, n in layers])
    sol = radial_solver(r, rho, K, mu, W, 1611.0, ("solid", "liquid", "solid"), (False, True, False), (False, False, False),
                        np.array([p["Rc"], Ro, RE]), degree_l=2, integration_rtol=1e-9, integration_atol=1e-12)
    if not sol.success:
        return None
    y = sol.result
    i = slice(0, N_CORE); rc = r[i]; rc[0] = max(rc[0], 1.0)
    muc = mu[0]; lam = k_core - 2 * muc / 3
    A = (2 * y[0, i] - 6 * y[2, i]) / rc
    div = (y[1, i] - lam * A) / (lam + 2 * muc) + A            # dy1/dr from the y2 constitutive relation
    flux = lambda j: r[j] ** 2 * np.imag(np.conj(y[0, j]) * y[1, j] + 6 * np.conj(y[2, j]) * y[3, j] + np.conj(y[4, j]) * y[5, j] / (4 * np.pi * G))
    depth = p["Rc"] - rc
    eps = U_POLE * np.abs(np.interp(DEPTHS_KM * 1e3, depth[::-1], div.real[::-1]) + 1j * np.interp(DEPTHS_KM * 1e3, depth[::-1], div.imag[::-1]))
    return dict(k=sol.k, h=sol.h, l=sol.l, heat=-21 / 2 * (W * RE) ** 5 * ECC ** 2 * np.imag(sol.k) / G,
                core_frac=flux(N_CORE - 1) / flux(len(r) - 1),
                # the ocean top's displacement is read on the shell side (the static ocean does not carry y1)
                u_core=abs(y[0, N_CORE - 1]) * U_POLE, u_ocean_top=abs(y[0, N_CORE + N_OTHER]) * U_POLE,
                u_load=abs(y[0, N_CORE + N_OTHER] - y[0, N_CORE - 1]) * U_POLE,
                sigma_rr_core=abs(y[1, N_CORE - 1]) * U_POLE,
                eps_top=eps[0], eps=eps, lin_core=2 * abs(y[0, N_CORE - 1]) * U_POLE / p["Rc"])


def sample(rng, sc, eta_c=None):
    """One parameter set, log-uniform in the Table 5 ranges; rho_c from the total mass (None if outside its range)."""
    lu = lambda lo, hi: math.exp(rng.uniform(math.log(lo), math.log(hi)))
    p = {k: lu(*v) for k, v in sc.items() if k != "rho_c"}
    if eta_c is not None:
        p["eta_c"] = eta_c
    Ro = RE - p["D"]
    rho_c = (3 * M / (4 * np.pi) - p["rho_o"] * (Ro ** 3 - p["Rc"] ** 3) - p["rho_s"] * (RE ** 3 - Ro ** 3)) / p["Rc"] ** 3
    if not sc["rho_c"][0] <= rho_c <= sc["rho_c"][1]:
        return None
    p["rho_c"] = rho_c
    return p


def draw(rng, sc, n, k_core, eta_c=None):
    out = []
    while len(out) < n:
        p = sample(rng, sc, eta_c)
        if p is None:
            continue
        s = solve(p, k_core)
        if s is not None:
            out.append(s)
    return out


def stats(v):
    v = np.asarray(v); return v.mean(), v.std(ddof=1)


def main(out_dir, k_core):
    rng = np.random.default_rng(20261031)
    print(f"K_core = {k_core / 1e9:.4g} GPa")
    # 1. Validation against Bagheri et al. (2025) Table 3 (mean +- 3 sigma).
    table3 = {"dissipative_core": dict(k2=(0.0240, 0.0102), h2=(0.0459, 0.0294), l2=(0.0104, 0.0072), lag_k=(39.35, 19.72), lag_h=(2.95, 2.69), lag_l=(3.10, 2.85)),
              "dissipative_shell": dict(k2=(0.0317, 0.0130), h2=(0.0848, 0.0359), l2=(0.0200, 0.0090), lag_k=(28.91, 17.25), lag_h=(28.73, 17.18), lag_l=(29.71, 17.89))}
    print("Validation: Table 5 sampling, eta_core log-uniform over the scenario's range; filtered as in Table 3: heat 18-28 GW and >= 90% of it in the scenario's layer")
    for name, sc in SCENARIOS.items():
        ss = draw(rng, sc, 3000, k_core)
        keep = [s for s in ss if 18e9 <= s["heat"] <= 28e9 and ((s["core_frac"] >= 0.9) if name == "dissipative_core" else (s["core_frac"] <= 0.1))]
        print(f"  {name}: {len(keep)} of {len(ss)} samples pass the filter")
        for lab, f in [("|k2|", lambda s: abs(s["k"])), ("|h2|", lambda s: abs(s["h"])), ("|l2|", lambda s: abs(s["l"])),
                       ("k2 lag deg", lambda s: -np.angle(s["k"], deg=True)), ("h2 lag deg", lambda s: -np.angle(s["h"], deg=True)),
                       ("l2 lag deg", lambda s: -np.angle(s["l"], deg=True))]:
            key = {"|k2|": "k2", "|h2|": "h2", "|l2|": "l2", "k2 lag deg": "lag_k", "h2 lag deg": "lag_h", "l2 lag deg": "lag_l"}[lab]
            m_all, s_all = stats([f(s) for s in ss]); m, s = stats([f(s) for s in keep])
            t = table3[name][key]
            print(f"    {lab:10s} filtered {m:.4g} +- {3 * s:.3g} (3 sd) | unfiltered {m_all:.4g} +- {3 * s_all:.3g} | Table 3 {t[0]} +- {t[1]}")
    # 2. Strain table.
    rows, summ = [], []
    for name, sc in SCENARIOS.items():
        for eta in ETA_GRID[name]:
            ss = draw(rng, sc, 500, k_core, eta_c=eta)
            E = np.array([s["eps"] for s in ss])
            for j, d in enumerate(DEPTHS_KM):
                rows.append((name, eta, d, E[:, j].mean(), E[:, j].std(ddof=1)))
            q = {key: stats([s[key] for s in ss]) for key in ("u_core", "u_ocean_top", "u_load", "sigma_rr_core", "eps_top", "lin_core", "heat", "core_frac")}
            lv = {lab: stats([abs(s[lab]) for s in ss]) for lab in ("k", "h", "l")}
            summ.append((name, eta, len(ss)) + sum((q[k] for k in q), ()) + sum((lv[k] for k in lv), ()))
            print(f"  {name} eta {eta:.3g}: dilatation at the core surface {q['eps_top'][0]:.3e} +- {q['eps_top'][1]:.2e}, "
                  f"at 50 km {E[:, 50].mean():.3e}, at 100 km {E[:, 100].mean():.3e}; u_core {q['u_core'][0]:.3f} +- {q['u_core'][1]:.3f} m; "
                  f"2 u_core/Rc {q['lin_core'][0]:.3e}; seafloor sigma_rr {q['sigma_rr_core'][0]:.1f} Pa; heat {q['heat'][0] / 1e9:.1f} GW, core {q['core_frac'][0]:.2f}")
    with open(os.path.join(out_dir, "layered_core_strain.csv"), "w", newline="\n") as f:
        f.write("scenario,eta_core_Pa_s,depth_km,eps_v_pole_mean,eps_v_pole_sd,K_core_Pa\n")
        for r in rows:
            f.write(f"{r[0]},{r[1]:.6g},{r[2]},{r[3]:.6e},{r[4]:.6e},{k_core:.7g}\n")
    cols = ["u_core_m", "u_ocean_top_m", "u_ocean_minus_core_m", "seafloor_sigma_rr_Pa", "eps_v_surface", "two_u_core_over_Rc", "heat_W", "core_heat_fraction"]
    with open(os.path.join(out_dir, "layered_core_strain_summary.csv"), "w", newline="\n") as f:
        f.write("scenario,eta_core_Pa_s,n_samples," + ",".join(f"{c}_mean,{c}_sd" for c in cols + ["k2_abs", "h2_abs", "l2_abs"]) + ",K_core_Pa\n")
        for s in summ:
            f.write(f"{s[0]},{s[1]:.6g},{s[2]}," + ",".join(f"{v:.6e}" for v in s[3:]) + f",{k_core:.7g}\n")
    # 3. Sensitivity of the dilatation to the core's bulk modulus (one mid-range model).
    p = dict(D=28e3, eta_s=1e17, mu_s=3e9, rho_s=930, rho_o=1050, Rc=195e3, mu_c=17e9, eta_c=1e12)
    Ro = RE - p["D"]; p["rho_c"] = (3 * M / (4 * np.pi) - p["rho_o"] * (Ro ** 3 - p["Rc"] ** 3) - p["rho_s"] * (RE ** 3 - Ro ** 3)) / p["Rc"] ** 3
    for kc in (100e9, k_core, 20e9, 2e9):
        s = solve(p, kc)
        print(f"  sensitivity: K_core {kc / 1e9:g} GPa: dilatation at the core surface {s['eps_top']:.3e}, u_core {s['u_core']:.3f} m, |k2| {abs(s['k']):.4f}")


if __name__ == "__main__":
    args = sys.argv[1:]
    kc = K_CORE
    if "--k-core" in args:
        i = args.index("--k-core"); kc = float(args[i + 1]); del args[i:i + 2]
    main(args[0] if args else "data", kc)
