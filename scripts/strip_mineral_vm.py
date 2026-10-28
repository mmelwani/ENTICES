"""Comments out the molar volume (-Vm / Vm) of every non-gas PHASES entry in the coreclath databases.

Why (PI decision, 2026-10-28). These databases use the LLNL aqueous model, under which PHREEQC sets
every aqueous species' molar volume to zero. A mineral's reaction volume is then minus its whole
molar volume, so its log K moves with pressure by Vm dP / (ln10 RT): antigorite by -4.16 log units
between 66.7 and 261.2 atm at 160 C. That is not physical. Without the term, mineral equilibria are
evaluated at the database's reference pressure (1 bar to 100 C, the saturation curve above), which
with real dissolution volumes of tens of cm3/mol costs about 0.1-0.2 log units at 260 bar. Gases keep
their entries: pressure acts on them through the GAS_PHASE.

Each stripped line is prefixed "#vm# " (unhyphenated "Vm" lines included: PHREEQC accepts either), so
--restore reverses it exactly. The volumes are kept in data/mineral_molar_volumes.csv
(scripts/export_mineral_volumes.py reads the prefixed lines too) and in legacy/Core11_idealgas_mod_v4.dat.
coreclath_C4only.dat is not processed here: scripts/make_c4only.py derives it from the uncoupled
database and applies the same pass.

Usage (from the repo root):
  python scripts/strip_mineral_vm.py            # strip (STRIP_MINERAL_VM = True) the coupled and uncoupled databases
  python scripts/strip_mineral_vm.py --restore  # put the volumes back (strip_mineral_vm = False)
"""
import os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FILES = ["databases/coreclath_CH4coupled.dat", "databases/coreclath_CH4uncoupled.dat"]
STRIP_MINERAL_VM = True
PREFIX = "#vm# "
VM_RX = re.compile(r"^\s*-?[Vv][Mm]\s")
MARK = ("# Mineral molar volumes (-Vm) commented out with '#vm# ' by scripts/strip_mineral_vm.py "
        "(LLNL aqueous model: see that script); --restore reverses it.")
KEYS = {"SOLUTION_MASTER_SPECIES", "SOLUTION_SPECIES", "PHASES", "RATES", "EXCHANGE_MASTER_SPECIES",
        "EXCHANGE_SPECIES", "SURFACE_MASTER_SPECIES", "SURFACE_SPECIES", "LLNL_AQUEOUS_MODEL_PARAMETERS",
        "NAMED_EXPRESSIONS", "CALCULATE_VALUES", "END", "GAS_BINARY_PARAMETERS"}


def strip_lines(L, strip_mineral_vm=STRIP_MINERAL_VM):
    """The lines L with non-gas PHASES Vm lines commented out (or, with strip_mineral_vm False,
    restored). Returns (lines, n_changed, n_unhyphenated)."""
    out, sec, phase, n, unh = [], None, None, 0, 0
    for line in L:
        code = line.split("#")[0].rstrip()
        if code.strip() in KEYS:
            sec, phase = code.strip(), None
        elif sec == "PHASES" and code and not code[0].isspace() and "=" not in code:
            phase = code.split()[0]
        if strip_mineral_vm:
            if sec == "PHASES" and phase and not phase.endswith("(g)") and VM_RX.match(line):
                n += 1; unh += not re.match(r"^\s*-", line)
                line = PREFIX + line
        elif line.startswith(PREFIX):
            n += 1; line = line[len(PREFIX):]
        out.append(line)
    if strip_mineral_vm and n and MARK not in out:
        out.insert(0, MARK)
    if not strip_mineral_vm:
        out = [l for l in out if l != MARK]
    return out, n, unh


def main():
    restore = "--restore" in sys.argv[1:]
    for f in FILES:
        p = os.path.join(ROOT, f)
        text = open(p, "rb").read().decode("utf-8").replace("\r\n", "\n")
        lines, n, unh = strip_lines(text.split("\n"), strip_mineral_vm=not restore)
        open(p, "wb").write("\r\n".join(lines).encode("utf-8"))
        print(f"{f}: {'restored' if restore else 'stripped'} {n} mineral Vm lines" + ("" if restore else f" ({unh} unhyphenated)"))


if __name__ == "__main__":
    main()
