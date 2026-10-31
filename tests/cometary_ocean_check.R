# The cometary-ice initial ocean (cometary_ice_ocean(), an option; the default ocean is pure water).
#
# Generates a tidal_kinetic input with the cometary ocean (c4_only), runs only its first two
# simulations (the pure-water solutions, then the ices reacting into the ocean at 1.45 C and
# 66.7 atm) and reads <base>_initial_ocean.csv. Checks:
#  1. Leak check: every coupled carbon species is C(+4) (carbon_by_oxidation_state()).
#  2. Element closure: C, Mtg, S and Ntg in the ocean equal the ices' inputs (to 1e-9), and no
#     element enters that should not (REACTION matches phase names first, case-insensitively, so
#     "CO" would be the cobalt phase Co: CO is entered as CO2 minus O).
#  3. Expected numbers (v = 0.02, O2 included, NH3 not): coupled C(+4) 0.0866 molal; Mtg 0.0102;
#     pH ~2; H2 ~0; no bubble. Sulfate: the hand estimate 0.0086 molal S(+6) / 0.0036 S(-2) leaves
#     out the ethane term: "Mtg 2 H -2" reduces ethane's C from -3 to -4, drawing 2 e- per C2H6
#     from the coupled system, 0.0064 eq, i.e. 0.0008 molal more sulfate: 0.0094 S(+6), 0.0028 S(-2).
#     Checked within 10%.
#  4. include_O2 = FALSE: no sulfate; H2 from the CO shift (0.0344) less the ethane term (0.0032),
#     ~0.031 molal; pH ~3.8. include_NH3 = TRUE: N(-3) 0.0074 molal enters.
#  5. (2026-10-31) The defaults: include_NH3 = TRUE. f scales every input linearly (and the ocean's C, N, S).
#     include_less_volatile = TRUE adds the 21 less-volatile CHONS ices of data/comet67p_volatility.csv by
#     their elements (C 1.08, N 0.177, S 0.177 per 100 H2O x f): element closure, and under c4_only all
#     coupled C is C(+4).
# Checks 1-4 use the Table 5 set as printed (include_NH3 = FALSE).
#
# Run from the repo root:  Rscript tests/cometary_ocean_check.R

suppressPackageStartupMessages({ library(tidyverse); library(glue); library(R6) })
invisible(capture.output(source("ENTICES_LF_20260917.R")))

failures <- 0
check <- function(ok, what) {
  cat(if (isTRUE(ok)) "[PASS]" else "[FAIL]", what, "\n")
  if (!isTRUE(ok)) failures <<- failures + 1
}
quiet <- function(expr) { invisible(capture.output(r <- suppressWarnings(suppressMessages(expr)))); r }
sim <- Simulator$new(k = 1e-07, layer_thickness = 194200, config = EnceladusConfig$new())
quiet(sim$run_simulation())
ph <- quiet(PhreeqcIntegrator$new(sim))
od <- tempfile("cometary_"); dir.create(od)

ocean0 <- function(ice, tag) {
  f <- quiet(ph$generate_phreeqc_input_tidal_kinetic(output_dir = file.path(od, tag), mode = "ungrouped", n_orbital_cycles = 1,
                                                     porewater_temp = 100, initial_ocean = ice))
  l <- readLines(f); cut <- file.path(dirname(f), "ocean0.pqi")
  writeLines(l[seq_len(which(l == "END")[2])], cut)                 # the solutions, then the cometary step
  db <- sub("^# Database: ", "", grep("^# Database: ", l, value = TRUE)[1])
  r <- run_phreeqc(cut, db, workdir = dirname(f), warn_version = FALSE)
  csv <- list.files(dirname(f), "_initial_ocean\\.csv$", full.names = TRUE)[1]
  x <- utils::read.delim(csv, check.names = FALSE, strip.white = TRUE); names(x) <- trimws(names(x))
  list(x = x[nrow(x), ], csv = csv, errors = r$errors, water = ph$phreeqc_params$water_mass)
}
s6 <- function(x) x$`m_SO4-2` + x$`m_HSO4-`
s2 <- function(x) x$m_H2S + x$`m_HS-`

# 1-3. Default: v = 0.02, O2 included, NH3 not.
ice <- cometary_ice_ocean(include_NH3 = FALSE)
o <- ocean0(ice, "default"); x <- o$x; m <- ice$mol_per_kgw
check(o$errors == 0, "cometary step runs (0 errors)")
st <- carbon_by_oxidation_state(o$csv)$by_state
cpl <- st[!st$decoupled, ]
check(all(cpl$ox[cpl$C_molal > 1e-15] == 4), sprintf("leak check: coupled C %.4e molal, all C(+4)", sum(cpl$C_molal)))
w0 <- o$water; w <- x$mass_H2O
inp <- c(C = m[["CO2"]] + m[["CO"]], Mtg = m[["CH4"]] + 2 * m[["C2H6"]], S = m[["H2S"]], Ntg = m[["N2"]]) * w0
got <- c(C = x$C, Mtg = x$Mtg, S = x$S, Ntg = x$Ntg) * w
other <- setdiff(intersect(database_elements("databases/coreclath_C4only.dat"), names(x)), c("H", "O", names(got)))
stray <- other[vapply(other, function(e) abs(as.numeric(x[[e]])) > 1e-20, logical(1))]
check(all(abs(got / inp - 1) < 1e-9) && !length(stray),
      sprintf("closure: C %.6e / %.6e, Mtg %.6e / %.6e, S %.6e / %.6e, Ntg %.6e / %.6e mol (ocean / input); other elements: %s",
              got[["C"]], inp[["C"]], got[["Mtg"]], inp[["Mtg"]], got[["S"]], inp[["S"]], got[["Ntg"]], inp[["Ntg"]],
              if (length(stray)) paste(stray, collapse = " ") else "none"))
check(abs(x$C / 0.0866 - 1) < 0.01 && abs(x$Mtg / 0.0102 - 1) < 0.01 && abs(x$pH - 2) < 0.5 && x$m_H2 < 1e-6 &&
        abs(s6(x) / 0.0094 - 1) < 0.1 && abs(s2(x) / 0.0028 - 1) < 0.1 && x$Vgas_L == 0,
      sprintf("v = 0.02: C(+4) %.4f, Mtg %.4f molal; pH %.2f; H2 %.1e; S(+6) %.4f, S(-2) %.4f molal; bubble %g L",
              x$C, x$Mtg, x$pH, x$m_H2, s6(x), s2(x), x$Vgas_L))

# 4. Toggles.
x <- ocean0(cometary_ice_ocean(include_O2 = FALSE, include_NH3 = FALSE), "noO2")$x
check(s6(x) < 1e-10 && abs(x$m_H2 / 0.031 - 1) < 0.15 && abs(x$pH - 3.8) < 0.5,
      sprintf("include_O2 = FALSE: S(+6) %.1e, H2 %.4f molal, pH %.2f", s6(x), x$m_H2, x$pH))
ice3 <- cometary_ice_ocean(include_NH3 = TRUE)
x <- ocean0(ice3, "NH3")$x
check(abs(x$N / ice3$mol_per_kgw[["NH3"]] - 1) < 1e-3 && abs(ice3$mol_per_kgw[["NH3"]] - 7.44e-3) < 1e-5,
      sprintf("include_NH3 = TRUE: N(-3) %.4e molal in (input %.4e); pH %.2f", x$N, ice3$mol_per_kgw[["NH3"]], x$pH))

# 5. Defaults, f, the less volatile ices.
d <- cometary_ice_ocean()
check(isTRUE(d$include_NH3) && !isTRUE(d$include_less_volatile) && d$f == 1, "defaults: include_NH3 = TRUE, include_less_volatile = FALSE, f = 1")
a1 <- cometary_ice_ocean(include_less_volatile = TRUE); a01 <- cometary_ice_ocean(f = 0.1, include_less_volatile = TRUE)
check(all(abs(a01$mol_per_kgw / a1$mol_per_kgw - 0.1) < 1e-12) && all(abs(a01$less_volatile_mol_per_kgw / a1$less_volatile_mol_per_kgw - 0.1) < 1e-12) &&
        length(a1$less_volatile_species) == 21 && abs(a1$less_volatile_mol_per_kgw[["C"]] * 100 * 0.0180153 - 1.0806) < 1e-3,
      sprintf("f = 0.1 scales every input by 0.1; less volatile: %d species, C %.4f per 100 H2O", length(a1$less_volatile_species),
              a1$less_volatile_mol_per_kgw[["C"]] * 100 * 0.0180153))
for (ff in c(1, 0.1)) {
  ice <- cometary_ice_ocean(f = ff, include_less_volatile = TRUE); o <- ocean0(ice, sprintf("lv_f%g", ff)); x <- o$x
  m <- ice$mol_per_kgw; lv <- ice$less_volatile_mol_per_kgw; w0 <- o$water; w <- x$mass_H2O
  inp <- c(C = m[["CO2"]] + m[["CO"]] + lv[["C"]], Mtg = m[["CH4"]] + 2 * m[["C2H6"]], S = m[["H2S"]] + lv[["S"]], N = m[["NH3"]] + lv[["N"]],
           Ntg = m[["N2"]]) * w0
  got <- c(C = x$C, Mtg = x$Mtg, S = x$S, N = x$N, Ntg = x$Ntg) * w
  for (k in grep("^g_", names(x), value = TRUE)) { fe <- formula_elements(sub("\\(g\\)$", "", sub("^g_", "", k)))
    for (e in intersect(names(fe), names(got))) got[e] <- got[e] + as.numeric(x[[k]]) * fe[[e]] }   # the ocean's gas bubble
  st <- carbon_by_oxidation_state(o$csv)$by_state; cpl <- st[!st$decoupled, ]
  check(o$errors == 0 && all(abs(got / inp - 1) < 1e-9) && all(cpl$ox[cpl$C_molal > 1e-15] == 4),
        sprintf("less volatile, f = %g: closure (solution + gas) C %.4e, S %.4e, N %.4e mol (worst %.1e); coupled C all C(+4); pH %.2f; bubble %.3g L", ff,
                got[["C"]], got[["S"]], got[["N"]], max(abs(got / inp - 1)), x$pH, x$Vgas_L))
}

unlink(od, recursive = TRUE)
cat(sprintf("\nSUMMARY: %d failure(s).\n", failures))
if (failures > 0) quit(status = 1)
