"""Gives every gas in the coreclath databases its Peng-Robinson critical constants (-T_c, -P_c, -Omega).

PHREEQC treats a GAS_PHASE with Peng-Robinson (Peng & Robinson 1976) only if every component has
all three constants, and cannot mix ideal and PR gases in one phase. On 2026-10-23 the decoupled
gases' constants were commented out so they could share a phase with the (then ideal) ordinary
gases. The PI's decision of 2026-10-27 is the reverse: every gas PR.

- Ordinary gases: their own commented constants are restored (sources as given in each entry:
  eGAS tables, NIST webbook, phreeqc.dat).
- Decoupled gases take their analogue's constants, as in phreeqc.dat distributed with PHREEQC 3.8.6:
  Mtg(g) <- CH4, Hdg(g) <- H2, H2Sg(g) <- H2S, Ntg(g) <- N2, Oxg(g) <- O2.
- H2(g) Omega: the entry's commented value was +0.225, labelled "phreeqc.dat"; phreeqc.dat (3.8.6)
  has -0.225 (H2's acentric factor is negative), which is used here, as Hdg(g) already had.

Each restored line ends with TAG, so "--ideal" comments those lines out again exactly.
coreclath_C4only.dat is regenerated from the uncoupled database by scripts/make_c4only.py.

Usage (from the repo root):
  python scripts/gas_eos.py           # restore the constants (Peng-Robinson gases)
  python scripts/gas_eos.py --ideal   # comment them out again (ideal gases)
"""
import os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FILES = ["databases/coreclath_CH4coupled.dat", "databases/coreclath_CH4uncoupled.dat"]
TAG = "  # PR constants restored 2026-10-28 (scripts/gas_eos.py)"
CONST_RX = re.compile(r"^#(\s+-(T_c|P_c|Omega)\b.*)$")
KEYS = {"SOLUTION_MASTER_SPECIES", "SOLUTION_SPECIES", "PHASES", "RATES", "EXCHANGE_MASTER_SPECIES",
        "EXCHANGE_SPECIES", "SURFACE_MASTER_SPECIES", "SURFACE_SPECIES", "LLNL_AQUEOUS_MODEL_PARAMETERS", "END"}


def entries(L):
    """(name, [line indices]) of the active PHASES entries."""
    out, sec, cur = [], None, None
    for i, l in enumerate(L):
        code = l.split("#")[0].rstrip()
        if code.strip() in KEYS: sec, cur = code.strip(), None; continue
        if sec != "PHASES": continue
        if code and not code[0].isspace() and "=" not in code:
            cur = (code.split()[0], [i]); out.append(cur)
        elif cur and not (l[:1] and not l[0].isspace() and not l.startswith("#")):
            if l.startswith("#") and re.match(r"^#[A-Za-z]", l): cur = None   # a commented-out entry starts
            elif cur: cur[1].append(i)
    return out


def restore(L):
    n = 0
    for name, idx in entries(L):
        if not name.endswith("(g)"): continue
        cand = [i for i in idx if CONST_RX.match(L[i])]
        if not cand: continue
        new23 = [i for i in cand if "commented out 2026-10-23" in L[i]]
        use = new23 if new23 else cand
        for i in use:
            line = L[i][1:]
            line = re.sub(r"\s*# commented out 2026-10-23:.*$", "", line)
            if name == "H2(g)" and "-Omega" in line:
                line = re.sub(r"-Omega(\s+)0\.225", r"-Omega\g<1>-0.225", line) + "  # -0.225 as in phreeqc.dat (was +0.225)"
            L[i] = line + TAG; n += 1
    return n


def make_ideal(L):
    n = 0
    for i, l in enumerate(L):
        if l.endswith(TAG):
            L[i] = "#" + l[: -len(TAG)]; n += 1
    return n


def main():
    ideal = "--ideal" in sys.argv[1:]
    for f in FILES:
        p = os.path.join(ROOT, f)
        L = open(p, "rb").read().decode("utf-8").replace("\r\n", "\n").split("\n")
        n = make_ideal(L) if ideal else restore(L)
        open(p, "wb").write("\r\n".join(L).encode("utf-8"))
        print(f"{f}: {'commented out' if ideal else 'restored'} {n} critical-constant lines")


if __name__ == "__main__":
    main()
