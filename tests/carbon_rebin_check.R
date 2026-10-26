# carbon_by_oxidation_state() (phreeqc_tools.R) on the N1 leak test at 100 C
# (tests/N1_leak_T100.pqi, run 2026-10-26 on main 255981f with PHREEQC 3.8.6 and
# coreclath_CH4uncoupled.dat; output in tests/output/N1_leak_T100.out and .csv).
#
# Run from the repo root:  Rscript tests/carbon_rebin_check.R

source("phreeqc_tools.R")

failures <- 0
check <- function(ok, what) {
  cat(if (isTRUE(ok)) "[PASS]" else "[FAIL]", what, "\n")
  if (!isTRUE(ok)) failures <<- failures + 1
}
close_to <- function(a, b, rel = 1e-4) isTRUE(abs(a - b) <= rel * abs(b))

out <- carbon_by_oxidation_state("tests/output/N1_leak_T100.out")
st <- out$by_state

# 1. All 10 distribution blocks parsed, every carbon species classified.
check(max(st$block) == 10, "10 'Distribution of species' blocks")
check(!anyNA(out$by_species$ox), "every carbon species has an oxidation state")

# 2. Re-binned totals reproduce PHREEQC's own total C in every block. The output
#    prints group totals to 4 significant figures, so the tolerance is 5e-4.
lines <- readLines("tests/output/N1_leak_T100.out")
starts <- grep("Distribution of species", lines, fixed = TRUE)
ctx_of <- function(b) { x <- out$by_species$context[out$by_species$block == b]; if (length(x)) x[1] else "no carbon" }
for (b in seq_along(starts)) {
  end <- starts[b] + grep("Saturation indices", lines[starts[b]:length(lines)], fixed = TRUE)[1] - 1
  g <- grep("^C\\([^)]*\\)\\s+[-0-9.e+]+\\s*$", lines[starts[b]:end], value = TRUE)
  phreeqc_C <- sum(as.numeric(sub("^\\S+\\s+", "", g)))
  ours <- sum(st$C_molal[st$block == b & !st$decoupled])
  check(if (phreeqc_C == 0) ours == 0 else close_to(ours, phreeqc_C, 5e-4),
        sprintf("block %2d (%s): re-binned C %.4e vs PHREEQC total C %.4e", b, ctx_of(b), ours, phreeqc_C))
}

# 3. Orbit-2 porewater (block 7) against an independent implementation
#    (the Python parser used for the 2026-10-26 N1 report).
ref <- c(`-3` = 8.5120000000e-15, `-2` = 3.9800006944e-17, `0` = 2.0884428201e-15,
         `2` = 3.3390000198e-11, `4` = 4.5559442186e-08)
b7 <- st[st$block == 7, ]
for (ox in names(ref))
  check(close_to(sum(b7$C_molal[abs(b7$ox - as.numeric(ox)) < 1e-9]), ref[[ox]]),
        sprintf("block 7 C(%s) = %.6e", ox, ref[[ox]]))
# PHREEQC files formate, acetate, methanol under C(4); the re-binning must not.
sp7 <- out$by_species[out$by_species$block == 7, ]
check(all(sp7$ox[sp7$species %in% c("HCOO-", "CH3COO-", "CH3OH")] < 4), "formate, acetate, methanol are not C(+4)")

# 4. Selected-output mode: only the punched m_ species, so a subset. Ethane is
#    punched; it must match the output file (5 vs 4 printed figures: 5e-4).
so <- carbon_by_oxidation_state("tests/output/N1_leak_T100.csv")
e_so  <- sum(so$by_species$C_molal[so$by_species$block == 7 & so$by_species$species == "C2H6"])
e_out <- sum(sp7$C_molal[sp7$species == "C2H6"])
check(close_to(e_so, e_out, 5e-4), sprintf("selected output row 7 ethane C %.4e vs output %.4e", e_so, e_out))
check(sum(so$by_state$C_molal[so$by_state$block == 7]) <= sum(b7$C_molal) * (1 + 5e-4),
      "selected-output total is a subset of the output-file total")

# 5. Ligand multiplier, Mtg, and unknown species, on a synthetic block.
tmp <- tempfile(fileext = ".out")
writeLines(c("Using solution 1.", "----------------------------Distribution of species----------------------------",
             "C(4)          1.000e-06",
             "   Mn(CH3COO)3-    1.000e-08   1.000e-08    -8.000    -8.000     0.000     (0)  ",
             "   C9H9X           1.000e-09   1.000e-09    -9.000    -9.000     0.000     (0)  ",
             "Mtg           2.000e-07",
             "   Mtg             2.000e-07   2.000e-07    -6.699    -6.699     0.000     (0)  ",
             "------------------------------Saturation indices-------------------------------"), tmp)
w <- NULL
syn <- withCallingHandlers(carbon_by_oxidation_state(tmp), warning = function(c) { w <<- conditionMessage(c); invokeRestart("muffleWarning") })
sp <- syn$by_species
check(sp$nC[sp$species == "Mn(CH3COO)3-"] == 6 && sp$ox[sp$species == "Mn(CH3COO)3-"] == 0, "Mn(CH3COO)3- counts 6 C at C(0)")
check(isTRUE(sp$decoupled[sp$species == "Mtg"]) && sp$ox[sp$species == "Mtg"] == -4, "Mtg counted as decoupled C(-4)")
check(!is.null(w) && grepl("C9H9X", w), "unclassified species is warned about")

cat(sprintf("\nSUMMARY: %d failure(s).\n", failures))
if (failures > 0) quit(status = 1)
