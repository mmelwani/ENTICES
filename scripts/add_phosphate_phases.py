"""Adds the validated phosphate phases to databases/coreclath_CH4coupled.dat and coreclath_CH4uncoupled.dat
(production kit, PI 2026-11-02; docs/MODEL.md, Phosphate phases), and comments out coreclath's own
Hydroxyapatite and Fluorapatite. Then rerun scripts/make_c4only.py to carry the change into
coreclath_C4only.dat.

Sources, all converted to the coreclath basis (HPO4-2 for P(5), NH3 for N(-3)):
  - Thermoddem V1.10 (BRGM; PHREEQC_ThermoddemV1.10_15Dec2020.dat, "BDD Date 12/15/2020", shipped with
    PHREEQC 3.8.6 and in github.com/usgs-coupled/phreeqc3, database/), read from that file:
      Hydroxyapatite(Natur), Whitlockite(high), Brushite (log_k, delta_h, -analytic) and
      Ca4H(PO4)3:3H2O (OCP), Bobbierite, Newberyite (log_k at 25 C only).
    Each H2PO4- of Thermoddem's reaction becomes HPO4-2 + H+, with Thermoddem's own
    H2PO4- = HPO4-2 + H+ (log_k -7.212, delta_h 3.600 kJ/mol, -analytic): log_k, delta_h and the analytic
    coefficients add n times that reaction's.
  - Vivianite: Al-Borno & Tomson (1994, GCA 58, 5373), pKsp(T) = -234.205 + 12242.6/T + 92.510 log10 T
    (5-90 C; Ksp = (Fe+2)^3 (PO4-3)^2), with coreclath's HPO4-2 = PO4-3 + H+. Thermoddem's -3.272 is the
    same value with the sign flipped (README).
  - Struvite: Bhuiyan, Mavinic & Beckie (2007, Environ. Technol. 28, 1015), log10 Ksp = -1157.45 - 0.784 T -
    63.86/T + 556.83 log10 T + 19.54/T^2 (10-60 C; Ksp = (Mg+2)(NH4+)(PO4-3)), with coreclath's
    HPO4-2 = PO4-3 + H+ and NH3 + H+ = NH4+.
Each phase was checked against a primary solubility product before it was added (docs/MODEL.md, the table).

Idempotent: the inserted block sits between marker lines and is replaced on a rerun; the coreclath apatites
are commented with the prefix "#tdP# " once.

Usage (from the repo root):
    python scripts/add_phosphate_phases.py [--thermoddem <path to PHREEQC_ThermoddemV1.10_15Dec2020.dat>]
"""
import argparse, math, re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DBS = ["databases/coreclath_CH4coupled.dat", "databases/coreclath_CH4uncoupled.dat"]
DEFAULT_TD = "C:/Program Files/USGS/phreeqc-3.8.6-17100-x64/database/PHREEQC_ThermoddemV1.10_15Dec2020.dat"
BEGIN = "# ==== Phosphate phases (scripts/add_phosphate_phases.py; docs/MODEL.md, Phosphate phases) ===="
END = "# ==== end of phosphate phases ===="
PREFIX = "#tdP# "
# coreclath's (LLNL) H2PO4-/PO4-3/NH4+ reactions, as written in the coreclath databases
CC_HPO4_PO4 = [-7.6170e1, -3.3574e-2, 1.3405e2, 2.9658e1, 2.1140]     # HPO4-2 = PO4-3 + H+      (log_k -12.3218)
CC_NH3_NH4 = [-1.4527e1, -5.0518e-3, 3.0447e3, 6.0865, 4.7515e1]      # NH3 + H+ = NH4+           (log_k 9.2410)

# Thermoddem phases taken: (Thermoddem name, coreclath name, temperature note)
TD_PHASES = [
    ("Hydroxyapatite(Natur)", "Hydroxyapatite", "0-300 C (Thermoddem analytic)"),
    ("Whitlockite(high)", "Whitlockite", "0-300 C (Thermoddem analytic); beta-Ca3(PO4)2 as measured by Gregory et al. 1974"),
    ("Brushite", "Brushite", "0-300 C (Thermoddem analytic); measured 5-37.5 C (Gregory et al. 1970)"),
    ("Ca4H(PO4)3:3H2O", "OCP", "25 C only (no enthalpy): treat as reliable to ~50 C"),
    ("Bobbierite", "Bobbierite", "25 C only (no enthalpy): treat as reliable to ~50 C"),
    ("Newberyite", "Newberyite", "25 C only (no enthalpy): treat as reliable to ~50 C"),
]


def an(c, T):
    return c[0] + c[1] * T + c[2] / T + c[3] * math.log10(T) + c[4] / T ** 2


def read_td(path):
    t = Path(path).read_text(encoding="utf-8-sig", errors="replace").replace("\ufffd", "").split("\n")
    # H2PO4- = HPO4-2 + H+
    i = next(k for k, l in enumerate(t) if re.match(r"\s*1\.000H2PO4-\s*=\s*HPO4-2 \+ 1\.000H\+", l))
    sp = parse_entry(t, i + 1)
    ph = next(k for k, l in enumerate(t) if l.strip() == "PHASES")
    out = {}
    for name, _, _ in TD_PHASES:
        k = next(j for j in range(ph, len(t)) if t[j].strip() == name)
        e = parse_entry(t, k + 2); e["reaction"] = t[k + 1].strip(); e["refs"] = next((l.strip() for l in t[k + 2:k + 8] if "#References" in l), "")
        out[name] = e
    return sp, out


def parse_entry(t, i):
    e = {"log_k": None, "delta_h": None, "analytic": None}
    for l in t[i:i + 8]:
        s = l.strip()
        if s.startswith("log_k"): e["log_k"] = float(s.split()[1])
        elif s.startswith("delta_h"): e["delta_h"] = float(s.split()[1])
        elif s.startswith("-analytic"): e["analytic"] = [float(x) for x in s.split()[1:6]]
        elif s and not s.startswith(("#", "-", "log", "delta")): break
    return e


def convert(reaction, n_h2po4):
    """Rewrite a Thermoddem reaction 'X + aH+ = ... + nH2PO4- + ...' on the HPO4-2 basis."""
    lhs, rhs = [s.strip() for s in reaction.split("=")]
    terms = [x.strip() for x in lhs.split(" + ")]
    solid = terms[0]
    a = 0.0
    for x in terms[1:]:
        m = re.match(r"([0-9.]*)H\+$", x)
        if m: a = float(m.group(1) or 1)
    prods = []
    for x in [x.strip() for x in rhs.split(" + ")]:
        m = re.match(r"([0-9.]*)(.+)$", x)
        c, s = float(m.group(1) or 1), m.group(2)
        prods.append((c, "HPO4-2" if s == "H2PO4-" else s))
    h = a - n_h2po4
    fmt = lambda c, s: (s if abs(c - 1) < 1e-9 else f"{c:g} {s}")
    left = solid + (f" + {fmt(h, 'H+')}" if h > 1e-9 else "")
    right = " + ".join(fmt(c, s) for c, s in prods) + (f" + {fmt(-h, 'H+')}" if h < -1e-9 else "")
    return f"{left} = {right}"


def block(sp, td):
    L = [BEGIN,
         "# Converted to the HPO4-2 basis from Thermoddem V1.10 (BRGM, 15 Dec 2020) with Thermoddem's own",
         f"# H2PO4- = HPO4-2 + H+ (log_k {sp['log_k']}, delta_h {sp['delta_h']} kJ/mol); Vivianite and Struvite from the primary",
         "# literature. Each was checked against a primary solubility product (docs/MODEL.md). Molar volumes",
         "# are in data/mineral_molar_volumes.csv, not here (mineral -Vm lines are stripped in these databases).", ""]
    for name, cname, trange in TD_PHASES:
        e = td[name]
        n = float(re.search(r"([0-9.]+)H2PO4-", e["reaction"].split("=")[1]).group(1))
        L += [cname,
              f"\t{convert(e['reaction'], n)}",
              f"\tlog_k\t\t{e['log_k'] + n * sp['log_k']:.4f}"]
        if e["delta_h"] is not None: L.append(f"\t-delta_H\t{e['delta_h'] + n * sp['delta_h']:.3f}\tkJ/mol")
        if e["analytic"] is not None:
            L.append("\t-analytic\t" + " ".join(f"{a + n * b:.7e}" for a, b in zip(e["analytic"], sp["analytic"])))
        L += [f"#\tThermoddem V1.10 {name}: {e['reaction']}; log_k {e['log_k']}" +
              (f", delta_h {e['delta_h']} kJ/mol" if e["delta_h"] is not None else "") + f" ({e['refs'].replace('#References = ', '')})",
              f"#\tRange: {trange}", ""]
    # Vivianite (Al-Borno & Tomson 1994): log K = -pKsp(T) - 2 log K(HPO4-2 = PO4-3 + H+)
    viv = [234.205, 0.0, -12242.6, -92.510, 0.0]
    viv = [a - 2 * b for a, b in zip(viv, CC_HPO4_PO4)]
    L += ["Vivianite",
          "\tFe3(PO4)2:8H2O + 2 H+ = 3 Fe+2 + 2 HPO4-2 + 8 H2O",
          f"\tlog_k\t\t{an(viv, 298.15):.4f}",
          "\t-analytic\t" + " ".join(f"{x:.7e}" for x in viv),
          "#\tAl-Borno & Tomson (1994, GCA 58, 5373): pKsp = -234.205 + 12242.6/T + 92.510 log10 T, Ksp = (Fe+2)^3 (PO4-3)^2,",
          "#\tlog Ksp(25 C) = -35.767; with coreclath's HPO4-2 = PO4-3 + H+. Thermoddem's Vivianite (log_k -3.272 on H2PO4-, 94alb/tom)",
          "#\tis this value with the sign flipped (+3.272 reproduces it to 0.03).",
          "#\tRange: 5-90 C (the measurements)", ""]
    # Struvite (Bhuiyan et al. 2007): log K = log Ksp(T) - log K(HPO4-2 = PO4-3 + H+) - log K(NH3 + H+ = NH4+)
    stv = [-1157.45, -0.784, -63.86, 556.83, 19.54]
    stv = [a - b - c for a, b, c in zip(stv, CC_HPO4_PO4, CC_NH3_NH4)]
    L += ["Struvite",
          "\tMgNH4PO4:6H2O = Mg+2 + NH3 + HPO4-2 + 6 H2O",
          f"\tlog_k\t\t{an(stv, 298.15):.4f}",
          "\t-analytic\t" + " ".join(f"{x:.7e}" for x in stv),
          "#\tBhuiyan, Mavinic & Beckie (2007, Environ. Technol. 28, 1015): log10 Ksp = -1157.45 - 0.784 T - 63.86/T",
          "#\t+ 556.83 log10 T + 19.54/T^2, Ksp = (Mg+2)(NH4+)(PO4-3) (their fit; it gives pKsp 13.57 at 25 C against the",
          "#\ttabulated 13.36 +- 0.07); with coreclath's HPO4-2 = PO4-3 + H+ and NH3 + H+ = NH4+. Ohlinger et al. 1998: pKso 13.26.",
          "#\tRange: 10-60 C (the measurements)", "", END]
    return L


def comment_out(lines, name, note):
    """Comment out the PHASES entry `name` (its name line through the line before the next blank line); only
    within PHASES (RATES has a Hydroxyapatite too), and only once."""
    p0 = next(k for k, l in enumerate(lines) if l.strip() == "PHASES")
    p1 = next(k for k, l in enumerate(lines) if l.strip() == "RATES")
    if any(l.startswith(PREFIX + name) for l in lines[p0:p1]): return lines   # already done
    for i in range(p0, p1):
        l = lines[i]
        if l.rstrip() == name or l.startswith(name + " #") or l.startswith(name + "\t#"):
            j = i
            while j < len(lines) and lines[j].strip(): j += 1
            return lines[:i] + [PREFIX + note] + [PREFIX + x for x in lines[i:j]] + lines[j:]
    sys.exit(f"{name}: entry not found")


NOTE_HA = ("Replaced by Thermoddem Hydroxyapatite(Natur) (PI 2026-11-02). This LLNL entry (Ref RHF79: Robie, Hemingway & "
           "Fisher 1979, USGS Bulletin 1452, a compilation; what its apatite value rests on was not traced) is ~4 log units "
           "more soluble (log_k -3.07 vs -7.30 at 25 C); it may represent poorly crystalline apatite (unconfirmed). "
           "McDowell et al. 1977 measured the equivalent of -7.37 to -7.56 on this basis.")
NOTE_FA = ("Replaced (PI 2026-11-02): the phosphate set is Thermoddem's. Thermoddem Fluorapatite(Natur) is not added yet: "
           "its primary check (McCann 1968, pKs0 59.6 at 25 C) was not read; this LLNL entry is 2.6 log units less soluble.")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--thermoddem", default=DEFAULT_TD); a = ap.parse_args()
    sp, td = read_td(a.thermoddem)
    new = block(sp, td)
    for db in DBS:
        p = ROOT / db
        raw = p.read_bytes().decode("utf-8", errors="surrogateescape")
        lines = raw.replace("\r\n", "\n").split("\n")
        if BEGIN in lines:                                            # replace a previous block
            i, j = lines.index(BEGIN), lines.index(END)
            lines = lines[:i] + lines[j + 2:]                         # the block and the blank line after it
        lines = comment_out(lines, "Hydroxyapatite", NOTE_HA)
        lines = comment_out(lines, "Fluorapatite", NOTE_FA)
        r = next(k for k, l in enumerate(lines) if l.strip() == "RATES")
        lines = lines[:r] + new + [""] + lines[r:]
        p.write_bytes("\r\n".join(lines).encode("utf-8", errors="surrogateescape"))
        print(f"{db}: {len([l for l in new if l and not l[0] in '#\t'])} phases inserted before RATES")
    print("now run: python scripts/make_c4only.py")


if __name__ == "__main__":
    main()
