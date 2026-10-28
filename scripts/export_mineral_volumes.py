"""Write data/mineral_molar_volumes.csv: a molar volume (cm3/mol) for every non-gas PHASES entry of the
three coreclath databases, by the source priority the PI set on 2026-10-28:

  1. legacy/Core11_idealgas_mod_v4.dat -Vm (identical to coreclath for shared phases; Mohit's base
     database). The coreclath-only phases (the CH4/Mtg/CO2/H2S hydrates) take their coreclath value
     at this level, sourced in their own comments.
  2. Green et al. (2025) dataset tc-ds636.txt, from the bundle tc-thermoinput-igG25-H18corr-2024-12-21
     (THERMOCALC format; V in kJ/kbar = J/bar, x10 -> cm3/mol).
  3. Perple_X hp636ver.dat (TC-DS636, May 3 2023), raw file from
     https://github.com/jadconnolly/Perple_X/blob/main/datafiles/hp636ver.dat (V0 in J/bar, x10).
  4. Phase-specific values the PI supplied (Schreibersite, from mindat unit-cell data).
  Anything else is left blank and listed for the PI to source by hand.

Levels 2-3 are matched by an explicit name map, and a match is used only if its stoichiometry is
proportional to the PHREEQC formula (the formula-unit ratio scales the volume). Where level 2 or 3
also has a phase Core11 covers, a difference > 2% is reported (Core11's value is kept).

Molar volumes are read from both "-Vm" and unhyphenated "Vm" lines (PHREEQC accepts either), and also
from lines that make_databases.py's strip_mineral_vm pass has commented out ("#vm# -Vm ..."), so the
export gives the same table before and after the strip.

Usage (from the repo root):
  python scripts/export_mineral_volumes.py --ds636 <path to tc-ds636.txt> --hp636 <path to hp636ver.dat>
Without --ds636 / --hp636 those levels are skipped (and the script says so).
"""
import argparse, csv, re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DBS = {"coupled": "databases/coreclath_CH4coupled.dat", "uncoupled": "databases/coreclath_CH4uncoupled.dat",
       "C4only": "databases/coreclath_C4only.dat"}
CORE11 = "legacy/Core11_idealgas_mod_v4.dat"
VM_RX = re.compile(r"^\s*(?:#vm#\s*)?-?[Vv][Mm]\s+([-0-9.eE+]+)")
KEYWORDS = {"SOLUTION_MASTER_SPECIES", "SOLUTION_SPECIES", "PHASES", "RATES", "EXCHANGE_MASTER_SPECIES",
            "EXCHANGE_SPECIES", "SURFACE_MASTER_SPECIES", "SURFACE_SPECIES", "LLNL_AQUEOUS_MODEL_PARAMETERS",
            "NAMED_EXPRESSIONS", "CALCULATE_VALUES", "END", "GAS_BINARY_PARAMETERS"}

# PI-supplied values (level 4).
PI_VALUES = {
    "Schreibersite": (0.60221 * 366.10 / 8,
                      "mindat.org (min-3582), unit cell; Fe3P end-member assumed (mindat formula (Fe,Ni)3P; "
                      "ENTICES uses Fe3P). Tetragonal I-4, a = 9.05 A, c = 4.47 A, V = 366.10 A3, Z = 8 "
                      "[citation needed: primary crystal-structure reference listed by mindat]"),
}
# Clathrate cross-check against Sloan & Koh (2007), Clathrate Hydrates of Natural Gases, 3rd ed.,
# Table 2.5a (p. 78): structure I lattice parameter a (0 C unless marked), 46 H2O per cell, so for
# the database formula X:6H2O, Vm = 0.60221 a^3 / (46/6). The database values are kept (all within 2%).
SK_A = {"CH4_hydrate": (11.981, "Ogienko et al. 2006, via S&K"), "Mtg_hydrate": (11.981, "CH4 value; Ogienko et al. 2006, via S&K"),
        "CO2_hydrate": (12.07, "S&K"), "H2S_hydrate": (12.02, "S&K")}
def sk_note(p, vm):
    a, who = SK_A[p]; v = 0.60221 * a ** 3 / (46 / 6)
    return (f"; cross-check Sloan & Koh (2007) Table 2.5a p. 78, sI a = {a} A ({who}): {v:.2f} cm3/mol "
            f"({100 * (vm - v) / v:+.1f}% database vs S&K). Formula fixes the hydration number at 6; real sI "
            f"hydrates run ~5.75-6.4")
# Minerals the generators use that no database defines as a PHASE (Schreibersite enters as a KINETICS
# -formula; see ENTICES_LF_20260917.R).
EXTRA_FORMULA = {"Schreibersite": "Fe3P"}
# PHREEQC phase -> THERMOCALC / Perple_X endmember name (levels 2-3).
HP_MAP = {
    "Forsterite": "fo", "Fayalite": "fa", "Tephroite": "teph", "Larnite": "lrn", "Monticellite": "mont",
    "Enstatite": "en", "Ferrosilite": "fs", "Diopside": "di", "Hedenbergite": "hed", "Wollastonite": "wo",
    "Pseudowollastonite": "pswo", "Rhodonite": "rhod", "Jadeite": "jd", "Albite": "ab", "Anorthite": "an",
    "K-Feldspar": "san", "Sanidine_high": "san", "Microcline": "mic", "Quartz": "q", "Tridymite": "trd",
    "Cristobalite": "crst", "Coesite": "coe", "Stishovite": "stv", "Magnetite": "mt", "Hematite": "hem",
    "Troilite": "tro", "Pyrite": "pyr", "Fe": "iron", "Ni": "Ni", "Calcite": "cc", "Aragonite": "arag",
    "Dolomite": "dol", "Magnesite": "mag", "Siderite": "sid", "Rhodochrosite": "rhc", "Brucite": "br",
    "Periclase": "per", "Lime": "lime", "Corundum": "cor", "Diaspore": "dsp", "Goethite": "gth",
    "Chrysotile": "chr", "Lizardite": "liz", "Antigorite": "atg", "Talc": "ta", "Minnesotaite": "minn",
    "Greenalite": "glt", "Kaolinite": "kao", "Pyrophyllite": "prl", "Clinochlore-14A": "clin",
    "Daphnite-14A": "daph", "Clinozoisite": "cz", "Zoisite": "zo", "Epidote": "ep", "Epidote-ord": "ep",
    "Prehnite": "pre", "Laumontite": "lmt", "Wairakite": "wrk", "Stilbite": "stlb", "Heulandite": "heu",
    "Analcime": "anl", "Tremolite": "tr", "Anthophyllite": "anth", "Muscovite": "mu", "Paragonite": "pa",
    "Phlogopite": "phl", "Annite": "ann", "Anhydrite": "anh", "Halite": "hlt", "Sylvite": "syv",
    "Grossular": "gr", "Andradite": "andr", "Pyrope": "py", "Almandine": "alm", "Spinel": "sp",
    "Hercynite": "herc", "Rutile": "ru", "Ilmenite": "ilm", "Titanite": "sph", "Nepheline": "ne",
    "Kalsilite": "kls", "Leucite": "lc", "Akermanite": "ak", "Gehlenite": "geh", "Merwinite": "merw",
    "Rankinite": "rnk", "Spurrite": "spu", "Tilleyite": "ty", "Cordierite": "crd", "Eskolaite": "esk",
    "Bunsenite": "NiO", "Graphite": "gph", "Lawsonite": "law", "Margarite": "ma", "Riebeckite": "rieb",
    "Glaucophane": "gl", "Pargasite": "parg", "Tschermakite": "ts", "Andalusite": "and", "Kyanite": "ky",
    "Sillimanite": "sill", "Manganosite": "mang", "Bixbyite": "bix", "Pyrophanite": "pnt",
    "Geikielite": "geik", "Baddeleyite": "bdy", "Zircon": "zrc", "Ulvospinel": "usp",
    "Magnesioferrite": "mft", "Magnesiochromite": "picr", "Tenorite": "ten", "Cuprite": "cup", "Cu": "Cu",
    "Sulfur": "S", "Covellite": "Cv", "Chalcocite": "Cc", "Chalcopyrite": "Ccp", "Bornite": "Bn",
    "Azurite": "Azr", "Malachite": "Mal", "Celadonite": "cel", "Sudoite": "sud", "Deerite": "deer",
}
TC_ELEMENTS = {1: "Si", 2: "Ti", 3: "Al", 4: "Fe", 5: "Mg", 6: "Mn", 7: "Ca", 8: "Na", 9: "K", 10: "O",
               11: "H", 12: "C", 13: "Cl", 15: "Ni", 17: "S", 19: "Cr"}
OXIDES = {"Na2O": {"Na": 2, "O": 1}, "MgO": {"Mg": 1, "O": 1}, "Al2O3": {"Al": 2, "O": 3}, "SiO2": {"Si": 1, "O": 2},
          "K2O": {"K": 2, "O": 1}, "CaO": {"Ca": 1, "O": 1}, "TiO2": {"Ti": 1, "O": 2}, "MnO": {"Mn": 1, "O": 1},
          "FeO": {"Fe": 1, "O": 1}, "NiO": {"Ni": 1, "O": 1}, "ZrO2": {"Zr": 1, "O": 2}, "Cl2": {"Cl": 2},
          "O2": {"O": 2}, "H2O": {"H": 2, "O": 1}, "CO2": {"C": 1, "O": 2}, "CuO": {"Cu": 1, "O": 1},
          "Cr2O3": {"Cr": 2, "O": 3}, "S2": {"S": 2}, "F2": {"F": 2}, "N2": {"N": 2}, "ZnO": {"Zn": 1, "O": 1}}

def parse_formula(f):
    """Element counts of a PHREEQC formula: parentheses, ':' hydrates, decimal coefficients."""
    def group(s, i):
        out = {}
        while i < len(s) and s[i] != ")":
            if s[i] == "(":
                sub, i = group(s, i + 1); i += 1
                m = re.match(r"[0-9.]+", s[i:]); k = float(m.group()) if m else 1.0; i += m.end() if m else 0
                for e, n in sub.items(): out[e] = out.get(e, 0) + n * k
                continue
            m = re.match(r"([A-Z][a-z]*)([0-9.]*)", s[i:])
            if not m: raise ValueError(f"cannot parse {f!r} at {s[i:]!r}")
            out[m.group(1)] = out.get(m.group(1), 0) + (float(m.group(2)) if m.group(2) else 1.0); i += m.end()
        return out, i
    total = {}
    for part in f.split(":"):
        m = re.match(r"^([0-9.]+)(.*)", part); k, body = (float(m.group(1)), m.group(2)) if m else (1.0, part)
        for e, n in group(body, 0)[0].items(): total[e] = total.get(e, 0) + k * n
    return total

def read_db(path):
    """{phase: dict(formula, vm, vm_hyphen, comment)} and {element: gfw} from a PHREEQC database."""
    lines = (ROOT / path).read_text(encoding="utf-8", errors="replace").splitlines()
    sec = None; cur = None; phases = {}; gfw = {}
    for l in lines:
        code = l.split("#", 1)[0]
        first = code.split()[0] if code.split() else ""
        if first in KEYWORDS and len(code.split()) == 1: sec = first; cur = None; continue
        if sec == "SOLUTION_MASTER_SPECIES" and code.strip():
            t = code.split()
            if len(t) >= 5 and "(" not in t[0]:
                try: gfw[t[0]] = float(t[4])
                except ValueError: pass
        if sec != "PHASES": continue
        m = VM_RX.match(l)
        if m and cur:
            if phases[cur]["vm"] is None:
                phases[cur].update(vm=float(m.group(1)), vm_hyphen=bool(re.match(r"^\s*(?:#vm#\s*)?-", l)),
                                   comment=l.split("#", 1)[1].strip() if "#" in l.replace("#vm#", "", 1) else "")
            continue
        if re.match(r"^[A-Za-z(\[]", code) and "=" not in code:
            cur = first; phases[cur] = dict(formula=None, vm=None, vm_hyphen=None, comment="")
        elif cur and phases[cur]["formula"] is None and "=" in code:
            phases[cur]["formula"] = code.split("=")[0].split("+")[0].split()[0]
    return phases, gfw

def read_ds636(path):
    """{name: (V cm3/mol, element counts)} from a THERMOCALC tc-ds file."""
    L = Path(path).read_text(encoding="utf-8", errors="replace").splitlines(); out = {}
    for i, l in enumerate(L):
        t = l.split()
        if len(t) >= 4 and re.match(r"^[A-Za-z]", l) and t[2] == "0" and i + 1 < len(L):
            nums = t[3:]; comp = {}
            for j in range(0, len(nums) - 1, 2):
                if nums[j] == "0": break
                k = int(nums[j]); comp[TC_ELEMENTS.get(k, f"?{k}")] = float(nums[j + 1])
            v = L[i + 1].split()
            if len(v) >= 4 and "?14" not in comp:
                out[t[0]] = (float(v[3]) * 10.0, comp)
    return out

def read_hp636(path):
    """{name: (V cm3/mol, element counts)} from a Perple_X data file (V0 in J/bar)."""
    L = Path(path).read_text(encoding="utf-8", errors="replace").splitlines(); out = {}
    for i, l in enumerate(L):
        if "EoS" in l and not l.lstrip().startswith("|") and i + 2 < len(L):
            name = l.split()[0]; comp = {}
            for ox, n in re.findall(r"([A-Za-z0-9]+)\(([-0-9.]+)\)", L[i + 1]):
                for e, k in OXIDES.get(ox, {f"?{ox}": 1}).items(): comp[e] = comp.get(e, 0) + k * float(n)
            m = re.search(r"V0\s*=\s*([-0-9.eE+]+)", L[i + 2])
            if m: out[name] = (float(m.group(1)) * 10.0, {e: n for e, n in comp.items() if abs(n) > 1e-12})
    return out

def scale(phreeqc, hp):
    """Formula-unit ratio PHREEQC/HP if the stoichiometries are proportional, else None."""
    if set(phreeqc) != set(hp): return None
    r = [phreeqc[e] / hp[e] for e in phreeqc]
    return r[0] if max(r) - min(r) <= 1e-6 * abs(r[0]) else None

def production_phases():
    """Minerals in the generators' EQUILIBRIUM_PHASES and KINETICS lists (gases excluded)."""
    src = (ROOT / "ENTICES_LF_20260917.R").read_text(encoding="utf-8")
    names = set(re.findall(r"^   ([A-Za-z][\w\-\(\):.]*)\s+0\s+(?:0|\{)", src, re.M))
    for vec in re.findall(r"secondary_minerals <- c\((.*?)\)\s*\n", src, re.S):
        names |= set(re.findall(r'"([^"]+)"', vec))
    names |= {"Enstatite", "Forsterite", "Troilite", "Albite", "Diopside", "Anorthite", "K-Feldspar", "Tephroite",
              "Schreibersite", "Fe", "Ni", "Lawrencite", "Quartz", "Antigorite", "Foshagite", "Hatrurite", "Epidote",
              "Epidote-ord", "CH4_hydrate", "Mtg_hydrate", "CO2_hydrate", "H2S_hydrate"}
    return {n for n in names if not n.endswith("(g)")}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ds636"); ap.add_argument("--hp636"); ap.add_argument("--out", default="data/mineral_molar_volumes.csv")
    a = ap.parse_args()
    core, gfw = read_db(CORE11)
    dbs = {k: read_db(p)[0] for k, p in DBS.items()}
    ds = read_ds636(a.ds636) if a.ds636 else {}
    hp = read_hp636(a.hp636) if a.hp636 else {}
    if not a.ds636: print("note: --ds636 not given; level 2 skipped")
    if not a.hp636: print("note: --hp636 not given; level 3 skipped")
    prod = production_phases()
    phases = sorted({p for d in dbs.values() for p, v in d.items() if not p.endswith("(g)")} | set(EXTRA_FORMULA), key=str.lower)
    rows, conflicts, diffs, unfilled, unhyph = [], [], [], [], {}
    for k, d in dbs.items():
        unhyph[k] = sorted(p for p, v in d.items() if v["vm_hyphen"] is False)
    unhyph["Core11"] = sorted(p for p, v in core.items() if v["vm_hyphen"] is False)
    for p in phases:
        ent = next((d[p] for d in dbs.values() if p in d), {"formula": None})
        formula = ent["formula"] or (core.get(p) or {}).get("formula") or EXTRA_FORMULA.get(p)
        try: comp = parse_formula(formula)
        except Exception: comp = None
        mm = sum(gfw.get(e, float("nan")) * n for e, n in comp.items()) if comp else float("nan")
        in_dbs = ";".join(k for k, d in dbs.items() if p in d) or "none (KINETICS -formula only)"
        vals = {k: d[p]["vm"] for k, d in dbs.items() if p in d and d[p]["vm"] is not None}
        if len(set(vals.values())) > 1: conflicts.append((p, vals))
        c11 = core.get(p, {}).get("vm")
        alt = {}
        hn = HP_MAP.get(p)
        for lvl, src, tag in ((2, ds, "Green et al. (2025) tc-ds636.txt (tc-thermoinput-igG25-H18corr-2024-12-21)"),
                              (3, hp, "Perple_X hp636ver.dat (TC-DS636, 2023-05-03), github.com/jadconnolly/Perple_X")):
            if hn and hn in src and comp:
                s = scale(comp, src[hn][1])
                if s is not None: alt[lvl] = (src[hn][0] * s, f"{tag}, endmember {hn}" + (f" x{s:g}" if s != 1 else ""))
                else: print(f"warning: {p} -> {hn}: stoichiometry differs in level {lvl}; not used")
        if c11 is not None:
            vm, source, lvl = c11, "legacy/Core11_idealgas_mod_v4.dat -Vm" + (f" ({core[p]['comment']})" if core[p]["comment"] else ""), 1
        elif vals and p not in core:
            v = next(iter(vals.values())); c = next(d[p]["comment"] for d in dbs.values() if p in d and d[p]["vm"] is not None)
            vm, source, lvl = v, "coreclath database (not in Core11)" + (f": {c}" if c else ""), 1
        elif 2 in alt: vm, source, lvl = alt[2][0], alt[2][1], 2
        elif 3 in alt: vm, source, lvl = alt[3][0], alt[3][1], 3
        elif p in PI_VALUES: vm, source, lvl = PI_VALUES[p][0], PI_VALUES[p][1], 4
        else: vm, source, lvl = None, "UNFILLED: for the PI to source", None
        if lvl == 1:
            for l2, (v2, s2) in alt.items():
                if abs(v2 - vm) > 0.02 * vm: diffs.append((p, vm, l2, v2, s2))
        if p in SK_A and vm is not None: source += sk_note(p, vm)
        if vm is None: unfilled.append(p)
        rows.append({"phase": p, "formula": formula, "molar_mass_g_mol": f"{mm:.4f}" if mm == mm else "",
                     "Vm_cm3_mol": f"{vm:.4f}" if vm is not None else "", "source": source,
                     "priority_level": lvl if lvl is not None else "", "databases": in_dbs,
                     "production": "yes" if p in prod else "no"})
    out = ROOT / a.out; out.parent.mkdir(exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n"); w.writeheader(); w.writerows(rows)
    n_by = {l: sum(1 for r in rows if r["priority_level"] == l) for l in (1, 2, 3, 4)}
    print(f"wrote {out.relative_to(ROOT)}: {len(rows)} phases; by level {n_by}; unfilled {len(unfilled)}")
    print("production phases without a volume:", sorted(set(unfilled) & prod) or "none")
    print("all phases without a volume:", unfilled)
    print("database conflicts:", conflicts or "none")
    print("level 2/3 differences > 2% from Core11 (Core11 kept):")
    for p, v1, l2, v2, s2 in diffs: print(f"  {p}: Core11 {v1:.3f}, level {l2} {v2:.3f} ({100*(v2-v1)/v1:+.1f}%), {s2}")
    print("unhyphenated Vm entries:", {k: v for k, v in unhyph.items()})
    missing = sorted(prod - set(phases))
    if missing: print("production names not found as PHASES:", missing)

if __name__ == "__main__":
    main()
