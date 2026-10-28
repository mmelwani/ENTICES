"""Regenerates databases/coreclath_C4only.dat from databases/coreclath_CH4uncoupled.dat.

coreclath_C4only.dat (carbon_mode "c4_only", the ENTICES default) is the uncoupled
database with every coupled carbon species, phase and master species whose carbon is
not +4 by actual oxidation state commented out (each line prefixed "#c4only# "), and a
header listing each removal. Nothing else is changed. Edit the uncoupled database, then
rerun this script; tests/c4only_db_check.R checks that the committed file is what this
script produces.

Usage (from the repo root):
    python scripts/make_c4only.py                       # rewrites databases/coreclath_C4only.dat
    python scripts/make_c4only.py <src.dat> <out.dat>
    python scripts/make_c4only.py <src.dat> <out.dat> --keep HCOOH HCOO-
        (--keep re-admits the named entries, for sensitivity runs; not for production)
    python scripts/make_c4only.py --no-strip-mineral-vm
        (keeps mineral -Vm lines; by default strip_mineral_vm is on, and the output gets the same
        pass that scripts/strip_mineral_vm.py gives the other coreclath databases)

Line endings: the source is read with any line ending and the output is written with
CRLF, like the other databases. The source hash in the header is of the source with LF
line endings, so it does not depend on how git checked the file out.
"""
import hashlib
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from strip_mineral_vm import strip_lines
args = [a for a in sys.argv[1:]]
strip_mineral_vm = "--no-strip-mineral-vm" not in args
args = [a for a in args if a != "--no-strip-mineral-vm"]
keep_extra = set()
if "--keep" in args:
    i = args.index("--keep")
    keep_extra = set(args[i + 1:])
    args = args[:i]
src = args[0] if len(args) > 0 else os.path.join(ROOT, "databases", "coreclath_CH4uncoupled.dat")
out = args[1] if len(args) > 1 else os.path.join(ROOT, "databases", "coreclath_C4only.dat")

text = open(src, "rb").read().decode("utf-8").replace("\r\n", "\n")
L = text.split("\n")
PREFIX = "#c4only# "

# name -> (what, mean C oxidation state). Oxidation states as in phreeqc_tools.R
# CARBON_MOIETIES (H +1, O -2, N -3, S -2).
ACETATE_COMPLEX = ("acetate complex", "0")
REMOVE_SPECIES = {
    "C2H6": ("ethane", "-3"), "C2H4": ("ethene", "-2"), "C3H8": ("propane", "-8/3"), "CO": ("carbon monoxide", "+2"),
    "HCOOH": ("formic acid", "+2"), "HCOO-": ("formate", "+2"),
    "CH3COOH": ("acetic acid", "0"), "CH3COO-": ("acetate", "0"),
    "CH3OH": ("methanol", "-2"), "CH3CH2OH": ("ethanol", "-2"), "CH2O": ("formaldehyde", "0"),
    "C3H6O2": ("propanoic acid", "-2/3"), "C3H5O2-": ("propanoate", "-2/3"),
    "CH3COCOOH": ("pyruvic acid", "+2/3"), "CH3COCOO-": ("pyruvate", "+2/3"),
    "C2H5NO2": ("glycine", "+1"), "C2H4NO2-": ("glycinate", "+1"), "C2H6NO2+": ("glycinium", "+1"),
    "Na(C2H4NO2)": ("Na-glycinate", "+1"), "C4H8N2O3": ("glycylglycine", "+1"),
    "C3H7NO2": ("alanine", "0"), "C2H5NO": ("acetamide", "0"), "C6H8O7": ("citric acid", "+1"),
    "CH3NH2": ("methylamine", "-2"), "CN-": ("cyanide", "+2"), "HCN": ("hydrogen cyanide", "+2"),
}
REMOVE_PHASES = {
    "Alanine": ("alanine (s)", "0"), "Glycine": ("glycine (s)", "+1"), "Citric_Acid": ("citric acid (s)", "+1"),
    "C": ("graphite", "0"), "Pyridine": ("pyridine (s)", "-0.4"),
    "KerogenC128": ("kerogen C128H68O7", "-0.42"), "KerogenC292": ("kerogen C292H288O12", "-0.90"),
    "KerogenC406": ("kerogen C406H528O19", "-1.21"), "KerogenC415": ("kerogen C415H698O22", "-1.58"),
    "KerogenC515": ("kerogen C515H596O72", "-0.88"),
    "TiC": ("titanium carbide", "-4"), "UC": ("uranium carbide", "-4"),
    "C2H4(g)": ("ethene gas", "-2"), "C2H6(g)": ("ethane gas", "-3"), "C3H8(g)": ("propane gas", "-8/3"),
    "CO(g)": ("carbon monoxide gas", "+2"),
}
REMOVE_MASTERS = {"C(-3)": "C2H6", "C(-2)": "C2H4", "C(+2)": "CO"}
C_RX = re.compile(r"C(?![a-z])")          # C, not Ca/Cl/Co/...


def has_c(s):
    return bool(C_RX.search(s)) or "Mtg" in s


# Parse the entries (line ranges) of the three sections that hold carbon.
KEYS = {"SOLUTION_MASTER_SPECIES", "SOLUTION_SPECIES", "PHASES", "RATES", "EXCHANGE_MASTER_SPECIES",
        "EXCHANGE_SPECIES", "SURFACE_MASTER_SPECIES", "SURFACE_SPECIES", "LLNL_AQUEOUS_MODEL_PARAMETERS", "END"}
sec, cur, entries = None, None, []
for i, line in enumerate(L):
    s = line.split("#")[0].rstrip()
    if s.strip() in KEYS:
        sec, cur = s.strip(), None
        continue
    if sec == "SOLUTION_MASTER_SPECIES" and s and not s[0].isspace():
        entries.append({"sec": "MASTER", "name": s.split()[0], "rxn": s.strip(), "lines": [i]})
    elif sec == "SOLUTION_SPECIES":
        if s and not s[0].isspace() and "=" in s:
            cur = {"sec": "SPECIES", "name": s.split("=")[1].split()[0], "rxn": s.strip(), "lines": [i]}
            entries.append(cur)
        elif cur and s.strip():
            cur["lines"].append(i)
    elif sec == "PHASES":
        if s and not s[0].isspace() and "=" not in s:
            cur = {"sec": "PHASE", "name": s.split()[0], "rxn": "", "lines": [i]}
            entries.append(cur)
        elif cur and s.strip():
            if not cur["rxn"] and "=" in s:
                cur["rxn"] = s.strip()
            cur["lines"].append(i)

# Select the removals: the named entries, then every complex of acetic acid.
removed, kept_c, seen = [], [], set()
for e in entries:
    n = e["name"]
    if n in keep_extra:
        continue
    if e["sec"] == "MASTER":
        if n in REMOVE_MASTERS:
            removed.append((e, "master species of " + REMOVE_MASTERS[n], ""))
        continue
    tab = REMOVE_SPECIES if e["sec"] == "SPECIES" else REMOVE_PHASES
    if n in tab:
        removed.append((e,) + tab[n])
        seen.add((e["sec"], n))
    elif e["sec"] == "SPECIES" and has_c(n):
        kept_c.append(n)
    elif e["sec"] == "PHASE" and e["rxn"] and has_c(e["rxn"].split("=")[0].split()[0]):
        kept_c.append(n)
for e in entries:
    n = e["name"]
    if e["sec"] == "SPECIES" and n not in keep_extra and ("SPECIES", n) not in seen and \
            ("(CH3COO)" in n or re.search(r"CH3COO\+|CH3COO$", n)):
        removed.append((e,) + ACETATE_COMPLEX)
        seen.add(("SPECIES", n))
        if n in kept_c:
            kept_c.remove(n)
missing = [n for t, d in (("SPECIES", REMOVE_SPECIES), ("PHASE", REMOVE_PHASES))
           for n in d if (t, n) not in seen and n not in keep_extra]
if missing:
    raise SystemExit("not found in source: " + ", ".join(missing))
gone = {e["name"] for e, *_ in removed}
for e in entries:  # a kept entry must not be written from a removed one
    if e["name"] in gone or e["sec"] == "MASTER":
        continue
    toks = set(re.split(r"\s+", e["rxn"].replace("=", " ")))
    if toks & gone:
        raise SystemExit(f"kept {e['name']} uses removed {toks & gone}: {e['rxn']}")

new = list(L)
for e, *_ in removed:
    for i in e["lines"]:
        new[i] = PREFIX + new[i]

sha = hashlib.sha1(text.encode("utf-8")).hexdigest()[:12]
H = ["# ============================================================================",
     "# coreclath_C4only.dat -- carbon_mode = \"c4_only\" (ENTICES)",
     f"# Derived 2026-10-26 from {os.path.basename(src)} (sha1 {sha}, LF line endings)",
     "# by scripts/make_c4only.py; regenerate it with that script, do not edit by hand.",
     "#",
     "# Coupled carbon is C(+4) only. Every coupled carbon species, phase and master",
     "# species whose carbon is not +4 by actual oxidation state (H +1, O -2, N -3,",
     "# S -2; phreeqc_tools.R CARBON_MOIETIES) is commented out below, each line",
     f"# prefixed '{PREFIX.strip()}' (grep for it). Nothing else is changed: every C(+4)",
     "# species and phase and every Mtg species and phase (Mtg, Mtg(g), Mtg_hydrate)",
     "# is kept, and reduced carbon enters only as decoupled Mtg.",
     "# Kept although they hold non-carbonate carbon: SCN- and CO(NH2)2 (urea) --",
     "# their carbon is +4 (S -2, N -3), so they are C(+4) by this rule.",
     "#",
     f"# Removed ({len(removed)} entries; source line numbers):"]
for e, what, ox in removed:
    lab = {"MASTER": "master ", "SPECIES": "species", "PHASE": "phase  "}[e["sec"]]
    oxs = f"C({ox})" if ox else ""
    H.append(f"#   {lab} {e['name']:<16} {what:<28} {oxs:<9} line {e['lines'][0] + 1}: {e['rxn'][:60]}")
if keep_extra:
    H.append("# SENSITIVITY VARIANT -- re-admitted: " + ", ".join(sorted(keep_extra)))
H.append(f"# Kept carbon species/phases ({len(kept_c)}):")
w = "#   "
for n in kept_c:
    if len(w) + len(n) > 80:
        H.append(w.rstrip(", "))
        w = "#   "
    w += n + ", "
H.append(w.rstrip(", "))
H += ["# ============================================================================", "#"]
new, n_vm, _ = strip_lines(new, strip_mineral_vm)  # a no-op if the source was already stripped
open(out, "wb").write("\r\n".join(H + new).encode("utf-8"))
print(f"{out}: removed {len(removed)} entries, kept {len(kept_c)} carbon species/phases")
