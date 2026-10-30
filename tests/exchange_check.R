# Exchange physics (exchange-physics branch, 2026-10-29): Simulator$calculate_exchange().
#
#  1. The polar strain amplitude FlowUtils$tidal_strain_amplitude() (Liao et al. 2020, Fig. 2 caption form)
#     reproduces Liao's Fig. 2a: with their R0 = 186 km, a 33 h period, e = 0.0047 and rho = 2400 kg/m3, the
#     elastic-limit (eta = 1e18) contours 1e-4.5, 1e-5, 1e-6, 1e-7 sit at mu = 2.1e7, 7.5e7, 8e8, 8.6e9 Pa
#     as read off the figure (within 15%); and at mu = 1e10 the 1e-5 contour meets eta ~ 1.6e12 (within 25%).
#     The 3/5 factor written in their section 2.2 text and Appendix A misses these by ~1.6x.
#  2. The degree-2 dilatation factor: 0.19245 x 3/5 for a full core, -> 0.19245 for a thin layer.
#  3. Exchange per half-orbit = min(poro, darcy, pore_volume) at k = 1e-14, 1e-13, 1e-12, 1e-11, 1e-7 m2
#     (production phi, full core, Darcy over the skin depth): darcy binds below k ~ 1.2e-11, poro at 1e-7; x_half <= 1; x_orbit = 1 - (1 - x_half)^2;
#     exchanging_fluid_mass = x_orbit * porewater mass; circulation time = ocean mass / exchange per orbit.
#  4. The pore-volume limit: the poroelastic term can never reach it at real porosities (a fluid matrix gives
#     k2 -> 3/2, so eps0 <= ~1.5e-4 and alpha eps_v <= ~1.7e-5 of V_layer), so it is exercised with phi = 1e-7
#     and a huge k: x_half = 1, x_orbit = 1.
#  5. The generated tidal_kinetic input stamps the exchange line.
#
# Generates inputs only; nothing is run in PHREEQC.
# Run from the repo root:  Rscript tests/exchange_check.R

suppressPackageStartupMessages({ library(tidyverse); library(glue); library(R6) })
invisible(capture.output(source("ENTICES_LF_20260917.R")))

failures <- 0
check <- function(ok, what) {
  cat(if (isTRUE(ok)) "[PASS]" else "[FAIL]", what, "\n")
  if (!isTRUE(ok)) failures <<- failures + 1
}
quiet <- function(expr) { invisible(capture.output(r <- suppressWarnings(suppressMessages(expr)))); r }

# 1. Liao et al. (2020) Fig. 2a.
eps_liao <- function(mu, eta) FlowUtils$tidal_strain_amplitude(mu, eta, 2400, 186e3, 0.0047, 33 * 3600)
fig <- data.frame(leps = c(-4.5, -5, -6, -7), mu = c(2.1e7, 7.5e7, 8e8, 8.6e9))
mu_fit <- vapply(fig$leps, function(l) 10^uniroot(function(lm) log10(eps_liao(10^lm, 1e18)) - l, c(4, 13))$root, 0)
check(all(abs(mu_fit / fig$mu - 1) < 0.15),
      sprintf("Liao Fig. 2a, elastic limit: contours at mu = %s Pa (figure: %s)",
              paste(sprintf("%.2g", mu_fit), collapse = ", "), paste(sprintf("%.2g", fig$mu), collapse = ", ")))
knee <- 10^uniroot(function(le) log10(eps_liao(1e10, 10^le)) + 5, c(9, 18))$root
check(abs(knee / 1.6e12 - 1) < 0.25, sprintf("Liao Fig. 2a, viscous knee: the 1e-5 contour meets mu = 1e10 at eta = %.3g (figure ~1.6e12)", knee))

# 2. Geometric factor.
f_full <- FlowUtils$degree2_dilatation_factor(194.2e3, 194.2e3); f_thin <- FlowUtils$degree2_dilatation_factor(194.2e3, 1)
check(abs(f_full - 0.19245 * 0.6) < 1e-5 && abs(f_thin - 0.19245) < 1e-4,
      sprintf("degree-2 dilatation factor: %.5f (full core), %.5f (thin layer)", f_full, f_thin))

# 3. The three limits.
cfg <- EnceladusConfig$new()
expect <- c("1e-14" = "darcy", "1e-13" = "darcy", "1e-12" = "darcy", "1e-11" = "darcy", "1e-07" = "poro")
for (k in c(1e-14, 1e-13, 1e-12, 1e-11, 1e-7)) {
  sim <- quiet({ s <- Simulator$new(k = k, layer_thickness = 194200, config = cfg, reactive_depth = "layer"); s$run_simulation(); s })
  x <- sim$exchange; m <- quiet(sim$get_summary_metrics())
  pw <- 1000 * sim$params$porosity * x$V_layer
  circ <- cfg$get_derived_properties()$ocean_mass / m$exchanging_fluid_mass * cfg$orbital_period / (365.25 * 24 * 3600)
  ok <- x$binding == expect[[format(k)]] && x$dV == min(x$dV_poro, x$dV_darcy, x$dV_pore_volume) && x$x_half <= 1 &&
    abs(x$x_orbit - (1 - (1 - x$x_half)^2)) < 1e-15 && abs(m$exchanging_fluid_mass / (x$x_orbit * pw) - 1) < 1e-12 &&
    abs(m$circulation_time_years / circ - 1) < 1e-9
  check(ok, sprintf("k = %.0e: poro %.3e, darcy %.3e, pore_volume %.3e m3 -> %s; x_half %.3e; circulation %.3g yr",
                    k, x$dV_poro, x$dV_darcy, x$dV_pore_volume, x$binding, x$x_half, m$circulation_time_years))
}

# 3b. The skin depth: delta = sqrt(c / omega) with Liao et al. (2020) Table A1's c, recomputed here from the
#     moduli; the Darcy volume over L_eff = min(L, delta), (k/mu)(dP/L_eff) A T/pi, equals the code's; and the
#     storage cross-check S dP delta A is the Darcy volume / (2/f), f = (K_m + 4/3 mu)/(K_u + 4/3 mu).
for (k in c(1e-14, 1e-7)) {
  sim <- quiet({ s <- Simulator$new(k = k, layer_thickness = 194200, config = cfg, reactive_depth = "layer"); s$run_simulation(); s })
  a <- sim$biot_alpha(); phi <- 0.32; Ks <- 10e9; Kf <- 2.2e9; mu <- cfg$core_rigidity
  Km <- (1 - a) * Ks; Ku <- Km + a^2 * Ks * Kf / (phi * Ks + (a - phi) * Kf)
  c_liao <- k / 1.8e-3 * (Km + 4/3 * mu) * (Ku - Km) / (a^2 * (Ku + 4/3 * mu))
  delta <- sqrt(c_liao / (2 * pi / cfg$orbital_period)); Leff <- min(194200, delta)
  darcy <- k / 1.8e-3 * sim$tidal_pressure() / Leff * 4 * pi * cfg$core_radius^2 * cfg$orbital_period / pi
  f <- (Km + 4/3 * mu) / (Ku + 4/3 * mu); x <- sim$exchange
  check(abs(x$delta / delta - 1) < 1e-9 && abs(x$dV_darcy / darcy - 1) < 2e-3 && abs(x$dV_darcy / x$dV_storage - 2 / f) < 0.01,
        sprintf("k = %.0e: skin depth %.4g m (Table A1: %.4g); darcy %.4g m3 (by hand %.4g); darcy / storage %.3f (2/f = %.3f)",
                k, x$delta, delta, x$dV_darcy, darcy, x$dV_darcy / x$dV_storage, 2 / f))
}

# 3c. The default skin box: the reacting box is the top L_eff = min(L, delta); V_pore = phi * V_skin, and the
#     poroelastic volume uses the skin's degree-2 factor. Poro binds at every k (the skin's pore volume scales
#     with delta as the Darcy volume does, and alpha eps_v / phi ~ 2e-6 is below what drains).
for (k in c(1e-14, 1e-12, 1e-11, 1e-7)) {
  sim <- quiet({ s <- Simulator$new(k = k, layer_thickness = 194200, config = cfg); s$run_simulation(); s })
  x <- sim$exchange; R <- cfg$core_radius; Le <- min(194200, x$delta)
  Vs <- 4 / 3 * pi * (R^3 - (R - Le)^3)
  ok <- abs(sim$reactive_thickness() / Le - 1) < 1e-12 && abs(x$V_pore / (0.32 * Vs) - 1) < 1e-12 &&
    abs(x$dV_poro / (sim$biot_alpha() * x$eps0 * FlowUtils$degree2_dilatation_factor(R, Le) * Vs) - 1) < 1e-12 && x$binding == "poro"
  check(ok, sprintf("skin box, k = %.0e: L_eff %.4g m; x_half %.3e (%s)", k, Le, x$x_half, x$binding))
}

# 3d. The pore pressure from the strain (F, 2026-10-30): dP = alpha eps_v M, M = 1/S; the strain factor is the
#     box's (0.192 for a thin skin, 0.115 for the full layer, G); in the skin, the storage volume S dP delta A
#     equals the strain volume alpha eps_v V_skin to O(delta/R), so Darcy / strain = 2/f there; an override of
#     max_tidal_stress is used as given.
s1 <- quiet({ s <- Simulator$new(k = 1e-14, layer_thickness = 194200, config = cfg); s$run_simulation(); s })
sl <- quiet({ s <- Simulator$new(k = 1e-14, layer_thickness = 194200, config = cfg, reactive_depth = "layer"); s$run_simulation(); s })
x1 <- s1$exchange
ok <- abs(x1$dP / (s1$biot_alpha() * x1$eps_v / s1$skin_depth()$S) - 1) < 1e-12 && x1$dP_source == "strain" &&
  abs(x1$strain_factor - FlowUtils$degree2_dilatation_factor(cfg$core_radius, x1$delta)) < 1e-15 && abs(x1$strain_factor - 0.19245) < 1e-4 &&
  abs(sl$exchange$strain_factor - 0.19245 * 0.6) < 1e-5 && abs(x1$dV_storage / x1$dV_poro - 1) < 1e-3
ov <- quiet({ s <- Simulator$new(k = 1e-14, layer_thickness = 194200, config = EnceladusConfig$new(list(max_tidal_stress = 1e5))); s$run_simulation(); s })
check(ok && ov$exchange$dP == 1e5 && ov$exchange$dP_source == "override",
      sprintf("dP = alpha eps_v M: %.4g Pa (skin, factor %.4f), %.4g Pa (layer, factor %.4f); storage / strain %.4f in the skin; override %g Pa",
              x1$dP, x1$strain_factor, sl$exchange$dP, sl$exchange$strain_factor, x1$dV_storage / x1$dV_poro, ov$exchange$dP))

# 4. Pore-volume limit.
eps_fluid <- FlowUtils$tidal_strain_amplitude(1e10, 1, cfg$core_density, cfg$core_radius, cfg$eccentricity, cfg$orbital_period)
cat(sprintf("       (fluid-matrix limit: eps0 = %.3g; alpha = 1 gives %.3g of V_layer per half-orbit)\n",
            eps_fluid, eps_fluid * FlowUtils$degree2_dilatation_factor(194.2e3, 194.2e3)))
big <- EnceladusConfig$new(list(core_viscosity = 1))
sim <- quiet({ s <- Simulator$new(k = 1, layer_thickness = 194200, porosity = 1e-7, config = big, alpha = 1); s$run_simulation(); s })
check(sim$exchange$binding == "pore_volume" && sim$exchange$x_half == 1 && sim$exchange$x_orbit == 1,
      sprintf("pore-volume limit: %s, x_half %g, x_orbit %g", sim$exchange$binding, sim$exchange$x_half, sim$exchange$x_orbit))

# 6. Porosity and the Biot coefficient (D, 2026-10-29): alpha = min(1, phi/phi_c) (critical porosity, Nur et al.
#    1998); an override below phi is an error; the box's solid mass is V (rho_core - phi rho_water).
s <- Simulator$new(k = 1e-7, layer_thickness = 194200, config = cfg)
s2 <- Simulator$new(k = 1e-7, layer_thickness = 194200, config = cfg, phi0 = 0.2)
s4 <- Simulator$new(k = 1e-7, layer_thickness = 194200, config = cfg, phi0 = 0.45)
so <- Simulator$new(k = 1e-7, layer_thickness = 194200, config = cfg, alpha = 0.5)
e <- tryCatch({ Simulator$new(k = 1e-7, layer_thickness = 194200, config = cfg, alpha = 0.2); "no error" }, error = function(e) conditionMessage(e))
check(s$params$phi_c == 0.4 && abs(s$biot_alpha() - 0.8) < 1e-12 && abs(s2$biot_alpha() - 0.5) < 1e-12 && s4$biot_alpha() == 1 && so$biot_alpha() == 0.5 &&
        grepl("must be >= the porosity", e),
      sprintf("alpha (phi_c 0.4): %.3g (phi 0.32), %.3g (phi 0.2), %.3g (phi 0.45 >= phi_c: a suspension), override 0.5; override 0.2 < phi errors", s$biot_alpha(), s2$biot_alpha(), s4$biot_alpha()))
p2 <- quiet({ s2$run_simulation(); ph2 <- PhreeqcIntegrator$new(s2); ph2$calculate_phreeqc_parameters(organic_wt_percent = 15) })
Vb <- 4 / 3 * pi * (cfg$core_radius^3 - (cfg$core_radius - s2$reactive_thickness())^3)
check(abs(p2$layer_mass / (Vb * (cfg$core_density - 0.2 * 1000)) - 1) < 1e-12 && abs(p2$porewater_mass / (1000 * 0.2 * Vb) - 1) < 1e-12,
      sprintf("phi0 = 0.2: box solid mass %.4g kg = V (rho_core - 0.2 rho_w); porewater %.4g kg = 0.2 V rho_w", p2$layer_mass, p2$porewater_mass))

# 5. The header.
od <- file.path(tempdir(), "exchange_check"); unlink(od, recursive = TRUE)
f <- quiet({ s <- Simulator$new(k = 1e-7, layer_thickness = 194200, config = cfg); s$run_simulation()
  ph <- PhreeqcIntegrator$new(s); ph$calculate_phreeqc_parameters(organic_wt_percent = 15)
  ph$generate_phreeqc_input_tidal_kinetic(output_dir = od, mode = "ungrouped", n_orbital_cycles = 1, porewater_temp = 25, trace_C_molal = 1e-14) })
check(any(grepl("^# Exchange per half-orbit .*limited by poro", readLines(f))) && any(grepl("^# Reacting box: top ", readLines(f))),
      "tidal_kinetic header stamps the reacting box and the exchange (limited by poro)")
unlink(od, recursive = TRUE)

cat(sprintf("\nSUMMARY: %d failure(s).\n", failures))
if (failures > 0) quit(status = 1)
