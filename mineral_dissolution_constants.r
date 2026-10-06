# From Core11_idealgas_mod_v2.dat, and/or ENTICES input files
# Constants originally from Palandri and Kharaka 2004
# H   = H+ activity (= 10^-pH for ideal solution)
# TK  = temperature in Kelvin
# fe3 = Fe3+ activity (only used for Pyrrhotite)

k_primary_mineral <- function(mineral, mol=1.0, H=1e-7, TK = 274.6, fe3 = 1e-30) {
    if (mineral == "Magnetite") {
        kacid <- 10^(-8.59)  * exp(-18.6e3/8.314 * (1/TK-1/298.15)) * H^0.279
        kneut <- 10^(-10.78) * exp(-18.6e3/8.314 * (1/TK-1/298.15))
        k <- kacid + kneut #mol/m^2/s
        mw <- 231.517 #g/mol
        SSA <- 5 #m^2/g #additional factor of 10 in Haller et al, not included here yet
        k <- k*SSA*mw*mol
    } else if (mineral == "Fayalite") {
        kacid <- 10^(-4.80)  * exp(-94.4e3/8.314 * (1/TK-1/298.15))* H^1.0
        kneut <- 10^(-12.80) * exp(-94.4e3/8.314 * (1/TK-1/298.15))
        k <- kacid + kneut
        mw <- 203.7555
        SSA <- 0.1
        k <- k*SSA*mw*mol
    } else if (mineral == "Anorthite") {
        kacid <- 10^(-3.50) * exp(-16.6e3/8.314 * (1/TK-1/298.15)) * H^1.411
        kneut <- 10^(-9.12) * exp(-17.8e3/8.314 * (1/TK-1/298.15))
        k <- kacid + kneut
        mw <- 278.164
        SSA <- 5
        k <- k*SSA*mw*mol
    } else if (mineral == "Albite") {
        kacid <- 10^(-10.16) * exp(-65.0e3/8.314 * (1/TK-1/298.15)) * H^0.457
        kneut <- 10^(-12.56) * exp(-69.8e3/8.314 * (1/TK-1/298.15))
        kbase <- 10^(-15.60) * exp(-71.0e3/8.314 * (1/TK-1/298.15)) * H^(-0.572)
        k <- kacid + kneut + kbase
        mw <- 262.1798
        SSA <- 5
        k <- k*SSA*mw*mol
    } else if (mineral == "Enstatite") {
        kacid <- 10^(-9.02)  * exp(-80.0e3/8.314 * (1/TK-1/298.15)) * H^0.6
        kneut <- 10^(-12.72) * exp(-80.0e3/8.314 * (1/TK-1/298.15))
        k <- kacid + kneut
        mw <- 100.3725
        SSA <- 0.1 #additional factor of 10 in Haller et al, not included here yet
        k <- k*SSA*mw*mol
    } else if (mineral == "Ferrosilite") {
        kacid <- 10^(-8.30)  * exp(-47.2e3/8.314 * (1/TK-1/298.15)) * H^0.650
        kneut <- 10^(-11.70) * exp(-66.1e3/8.314 * (1/TK-1/298.15))
        k <- kacid + kneut
        mw <- 131.9145
        SSA <- 0.1 #additional factor of 10 in Haller et al, not included here yet
        k <- k*SSA*mw*mol
    } else if (mineral == "Forsterite") {
        kacid <- 10^(-6.85)  * exp(-67.2e3/8.314 * (1/TK-1/298.15)) * H^0.47
        kneut <- 10^(-10.64) * exp(-79.0e3/8.314 * (1/TK-1/298.15))
        k <- kacid + kneut
        mw<- 140.6715
        SSA <- 0.1
        k <- k*SSA*mw*mol
    } else if (mineral == "Pyrrhotite") {
        # Hexagonal pyrrhotite (PK04); rate depends on Fe3+ as oxidant
        kacid <- 10^(-6.79) * exp(-63.0e3/8.314 * (1/TK-1/298.15)) * H^(-0.090) * fe3^0.356
        k <- kacid
        mw <- 87.913
        SSA <- 5
        k <- k*SSA*mw*mol
    } else if (mineral == "Troilite") {
        # No separate PK04 troilite entry; hexagonal pyrrhotite parameters used as proxy (same FeS chemistry)
        kacid <- 10^(-6.79) * exp(-63.0e3/8.314 * (1/TK-1/298.15)) * H^(-0.090) * fe3^0.356
        k <- kacid
        mw <- 87.913
        SSA <- 5
        k <- k*SSA*mw*mol
    } else if (mineral == "Diopside") {
        kacid <- 10^(-6.36) * exp(-96.1e3/8.314 * (1/TK-1/298.15)) * H^0.71
        kneut <- 10^(-11.11) * exp(-40.6e3/8.314 * (1/TK-1/298.15))
        k <- kacid + kneut
        mw <- 216.55
        SSA <- 0.1
        k <- k*SSA*mw*mol
    } else if (mineral == "K-Feldspar") {
        kacid <- 10^(-10.06) * exp(-51.7e3/8.314 * (1/TK-1/298.15)) * H^0.5
        kneut <- 10^(-12.41) * exp(-38.0e3/8.314 * (1/TK-1/298.15))
        kbase <- 10^(-21.20) * exp(-94.1e3/8.314 * (1/TK-1/298.15)) * H^(-0.823)
        k <- kacid + kneut + kbase
        mw <- 278.33
        SSA <- 5
        k <- k*SSA*mw*mol
    } else if (mineral == "Tephroite") {
        # Pincus et al. (2026) ACS Earth Space Chem 10, 1174-1184, Table S1,
        # citing Casey et al. (1993) GCA 57, 785-793 (synthetic tephroite,
        # anoxic). Ea decreases with increasing pH (Casey Fig 4), assumed
        # Ea=0 at pH~6.2, clamped at 0 above that. No base branch: pH>7
        # source data not fittable, so neutral rate is held flat at its
        # pH=7 value for all pH>7. See ENTICES_LF_*.R's Tephroite RATES
        # block (ported 1:1 here) for the full derivation.
        Ea <- -60.11 * log(-log10(H)) + 106.17
        if (Ea < 0) Ea <- 0
        kacid <- 10^(-4.46) * exp(-Ea*1e3/8.314 * (1/TK-1/298.15)) * H^0.47
        kneut <- 10^(-4.46) * exp(-Ea*1e3/8.314 * (1/TK-1/298.15)) * 1e-7^0.47
        k <- kacid + kneut
        mw <- 201.96
        SSA <- 0.1
        k <- k*SSA*mw*mol
    } else if (mineral == "IOM_labile") {
        k <- 10^(13) * exp(-180000 / (8.314 * TK)) * mol
    } else if (mineral == "IOM_mid") {
        k <- 10^(13) * exp(-215000 / (8.314 * TK)) * mol
    } else if (mineral == "IOM_refract") {
        k <- 10^(13) * exp(-285000 / (8.314 * TK)) * mol
    } else if (mineral == "IOM_CO2") {
        k <- 10^(13) * exp(-214000 / (8.314 * TK)) * mol
    } else if (mineral == "IOM_CH4") {
        k <- 10^(13) * exp(-228000 / (8.314 * TK)) * mol
    } else if (mineral == "IOM_N") {
        k <- 10^(13) * exp(-215000 / (8.314 * TK)) * mol
    } else if (mineral == "IOM_S") {
        k <- 10^(13) * exp(-215000 / (8.314 * TK)) * mol
    } else {
        stop("Unknown mineral: ", mineral)
    }
    k
}

dissolution_per_step <- function(mineral,mol=1.0, step_orbits=1.0, H=1e-7, TK = 274.6, fe3=1e-30) {
    k <- k_primary_mineral(mineral, mol, H, TK, fe3)
    step_seconds = step_orbits*118800 
    moles <- k*step_seconds

    moles
}

check_step_size_minerals <- function(mineral_list, step_orbits, H=1e-7, TK = 274.6) {
    for (mineral in mineral_list) {}
    mineral_list
}

# Read the last dk_* value for each primary mineral from a tidal kinetic CSV and
# return dissolution rates per second.  rate_mol_s is mol/s; rate_frac_s is 1/s
# (mol/s divided by remaining moles).  Missing dk/k columns and depleted minerals
# (k = 0) return NULL.
last_dissolution_rate <- function(csv_path) {
    primary_minerals <- c("Forsterite", "Fayalite", "Enstatite", "Ferrosilite",
                          "Pyrrhotite", "Anorthite", "Albite", "Magnetite",
                          "Troilite", "Diopside", "K-Feldspar", "Tephroite",
                          "IOM_CO2", "IOM_CH4","IOM_N", "IOM_S")
                          #"IOM_labile", "IOM_mid", "IOM_refract")

    dat <- read.delim(csv_path, sep = "\t", header = TRUE,
                      strip.white = TRUE, check.names = FALSE)
    dat <- dat[trimws(dat[["step"]]) != "-99", ]

    # Last REACT row: soln=1 with at least one k_* > 0
    k_cols <- paste0("k_", primary_minerals)
    present_k <- k_cols[k_cols %in% names(dat)]
    if (length(present_k) == 0)
        stop("No k_* columns found in: ", csv_path)

    is_react <- apply(dat[, present_k, drop = FALSE], 1,
                      function(r) any(suppressWarnings(as.numeric(r)) > 0, na.rm = TRUE))
    react_rows <- which(is_react)
    if (length(react_rows) == 0)
        stop("No REACT rows (k_* > 0) found in: ", csv_path)

    last_row <- dat[tail(react_rows, 1), ]

    # Time step from difference between last two REACT rows' cumulative times;
    # fall back to the first REACT row's Time_Years if only one REACT row exists.
    time_step_yr <- if (length(react_rows) >= 2) {
        t_last <- suppressWarnings(as.numeric(dat[tail(react_rows, 1), "Time_Years"]))
        t_prev <- suppressWarnings(as.numeric(dat[tail(react_rows, 2)[1], "Time_Years"]))
        t_last - t_prev
    } else {
        suppressWarnings(as.numeric(dat[react_rows[1], "Time_Years"]))
    }
    if (is.na(time_step_yr) || time_step_yr <= 0)
        stop("Could not determine a positive time step from Time_Years column")

    time_step_s <- time_step_yr * 365.25 * 24 * 3600

    rate_mol_s <- setNames(lapply(primary_minerals, function(m) {
        col <- paste0("dk_", m)
        if (!col %in% names(last_row)) return(NA_real_)
        -suppressWarnings(as.numeric(last_row[[col]])) / time_step_s
    }), primary_minerals)

    rate_frac_s <- setNames(lapply(primary_minerals, function(m) {
        r <- rate_mol_s[[m]]
        if (is.na(r)) return(NA_real_)
        k_col  <- paste0("k_", m)
        dk_col <- paste0("dk_", m)
        if (!k_col %in% names(last_row) || !dk_col %in% names(last_row)) return(NA_real_)
        k_curr <- suppressWarnings(as.numeric(last_row[[k_col]]))
        dk     <- suppressWarnings(as.numeric(last_row[[dk_col]]))
        # k_curr is the remaining moles AFTER this step; dk = k_curr - k_prior (<=0 as
        # the mineral dissolves), so the moles available at the START of this step
        # (the correct denominator for a fractional rate) is k_curr - dk, not k_curr.
        k_prior <- k_curr - dk
        if (is.na(k_prior) || k_prior <= 0) return(NA_real_)
        r / k_prior
    }), primary_minerals)

    list(rate_mol_s = rate_mol_s, rate_frac_s = rate_frac_s)
}

# Latest remaining moles (k_<mineral>) for each primary mineral, from the same
# last REACT row that last_dissolution_rate() grabs dk_<mineral> from — i.e. the
# moles remaining AFTER that row's dk change (k_curr, not k_prior). Missing k_*
# columns return NA.
latest_mineral_moles <- function(csv_path) {
    primary_minerals <- c("Forsterite", "Fayalite", "Enstatite", "Ferrosilite",
                          "Pyrrhotite", "Anorthite", "Albite", "Magnetite",
                          "Troilite", "Diopside", "K-Feldspar", "Tephroite",
                          "IOM_CO2", "IOM_CH4","IOM_N", "IOM_S")

    dat <- read.delim(csv_path, sep = "\t", header = TRUE,
                      strip.white = TRUE, check.names = FALSE)
    dat <- dat[trimws(dat[["step"]]) != "-99", ]

    # Last REACT row: soln=1 with at least one k_* > 0
    k_cols <- paste0("k_", primary_minerals)
    present_k <- k_cols[k_cols %in% names(dat)]
    if (length(present_k) == 0)
        stop("No k_* columns found in: ", csv_path)

    is_react <- apply(dat[, present_k, drop = FALSE], 1,
                      function(r) any(suppressWarnings(as.numeric(r)) > 0, na.rm = TRUE))
    react_rows <- which(is_react)
    if (length(react_rows) == 0)
        stop("No REACT rows (k_* > 0) found in: ", csv_path)

    last_row <- dat[tail(react_rows, 1), ]

    setNames(lapply(primary_minerals, function(m) {
        k_col <- paste0("k_", m)
        if (!k_col %in% names(last_row)) return(NA_real_)
        suppressWarnings(as.numeric(last_row[[k_col]]))
    }), primary_minerals)
}




# ---- Plot k vs pH and vs temperature, saving both PNGs ----
# Called explicitly (not gated by sys.nframe()) so it works the same whether
# this file is run via Rscript, sourced from a console, or sourced from an
# Rmd/knitr chunk -- sys.nframe() is never 0 inside knitr's evaluation, so a
# frame-depth guard here could never pass from a notebook.
plot_mineral_rates <- function(output_dir = ".") {
library(dplyr)
library(ggplot2)

minerals  <- c("Forsterite", "Fayalite", "Enstatite", "Ferrosilite",
               "Pyrrhotite", "Anorthite", "Albite", "Magnetite",
               "Troilite", "Diopside", "K-Feldspar", "Tephroite",
               "IOM_labile", "IOM_mid", "IOM_refract",
               "IOM_CO2", "IOM_CH4", "IOM_N", "IOM_S")

pH_range  <- seq(0, 14, length.out = 300)
H_vals    <- 10^(-pH_range)
TK_val    <- 274.6    # Enceladus ocean floor (~1.45 C)
fe3_val   <- 1e-30    # fixed Fe3+ activity for Pyrrhotite

rate_data <- do.call(rbind, lapply(minerals, function(m) {
    k_vals <- sapply(H_vals, function(H)
        k_primary_mineral(m, H = H, TK = TK_val, fe3 = fe3_val))
    data.frame(mineral = m, pH = pH_range, k = k_vals)
}))

mineral_color_map <- c(
  "Forsterite"  = "#E69F00", "Fayalite"    = "#56B4E9",
  "Enstatite"   = "#009E73", "Ferrosilite" = "#F0E442",
  "Pyrrhotite"  = "#0072B2", "Anorthite"   = "#D55E00",
  "Albite"      = "#CC79A7", "Magnetite"   = "#882255",
  # New primary assemblage (2026-09-18): grouped by mineral family with the
  # existing entries above -- Troilite (Pyrrhotite's FeS proxy/replacement),
  # Diopside (pyroxene, alongside Enstatite), K-Feldspar (feldspar, alongside
  # Albite/Anorthite), Tephroite (Mn-olivine, alongside Forsterite/Fayalite).
  "Troilite"    = "#A6761D", "Diopside"    = "#4DAC26",
  "K-Feldspar"  = "#762A83", "Tephroite"   = "#053061",
  # IOM pools: k here is a first-order decay constant (s^-1), not a surface-
  # area-normalised dissolution rate (mol m^-2 s^-1) like the rock minerals
  # above — shown on the same axis for convenience, not a like-for-like value.
  "IOM_labile"  = "#000000", "IOM_mid"     = "#666666", "IOM_refract" = "#B2B2B2",
  "IOM_CO2"     = "#7570B3", "IOM_CH4"     = "#E41A1C",
  "IOM_N"       = "#556B2F", "IOM_S"       = "#DAA520"
)

# Kerogen pools (labile/mid/refract): dotted. Product/decay-channel pools
# (CO2/CH4/N/S): dashed. Rock minerals: solid (default).
mineral_linetype_map <- c(
  "Forsterite"  = "solid", "Fayalite"    = "solid",
  "Enstatite"   = "solid", "Ferrosilite" = "solid",
  "Pyrrhotite"  = "solid", "Anorthite"   = "solid",
  "Albite"      = "solid", "Magnetite"   = "solid",
  "Troilite"    = "solid", "Diopside"    = "solid",
  "K-Feldspar"  = "solid", "Tephroite"   = "solid",
  "IOM_labile"  = "dotted", "IOM_mid"    = "dotted", "IOM_refract" = "dotted",
  "IOM_CO2"     = "dashed", "IOM_CH4"    = "dashed",
  "IOM_N"       = "dashed", "IOM_S"      = "dashed"
)

p_k_pH <- ggplot(rate_data, aes(x = pH, y = k, colour = mineral, linetype = mineral)) +
    geom_line(linewidth = 1.0) +
    scale_linetype_manual(values = mineral_linetype_map) +
    scale_y_log10(
        labels = scales::label_scientific(),
        sec.axis = sec_axis(
            transform = ~ 1 / (. * 365.25 * 24 * 3600),
            name   = expression(tau==1/k[total]~~"(yr)"),
            labels = scales::label_scientific()
        )
    ) +
    scale_color_manual(values=mineral_color_map) +
    labs(
        x       = "pH",
        y       = expression(k[total]~~"(mol" ~~ s^{-1} * ")"),
        title   = sprintf("Mineral & IOM dissolution/decomposition rate constants vs pH  [T = %.1f K]", TK_val),
        caption = sprintf("k assumes 1 mol present (mol/s, directly comparable).  Pyrrhotite Fe³⁺ fixed at %.0e.  IOM k is pH-independent.", fe3_val)
    ) +
    theme_minimal(base_size = 11) +
    theme(
        aspect.ratio      = 0.65,
        axis.line         = element_line(colour = "grey40"),
        axis.ticks        = element_line(colour = "grey40"),
        axis.ticks.length = unit(0.15, "cm"),
        plot.caption      = element_text(size = 7, colour = "grey45")
    )

ggsave(file.path(output_dir, "mineral_rate_constants_pH.png"), plot = p_k_pH, dpi = 300, width = 10, height = 5.5)

# ---- Plot k vs temperature (constant pH) ----
pH_val    <- 7
H_val     <- 10^(-pH_val)
TC_range  <- seq(0, 300, length.out = 300)   # degrees C
TK_range  <- TC_range + 273.15

rate_data_T <- do.call(rbind, lapply(minerals, function(m) {
    k_vals <- sapply(TK_range, function(TK)
        k_primary_mineral(m, H = H_val, TK = TK, fe3 = fe3_val))
    data.frame(mineral = m, TC = TC_range, k = k_vals)
}))

# k floor so that tau = 1/k never exceeds 1e10 years on the secondary axis.
# Gridlines are always drawn from the PRIMARY scale's breaks, so to get a
# gridline at every power-of-ten tau value, compute primary breaks by
# inverting the transform at each power-of-ten year value, rather than using
# round powers of ten in k directly (the two can't both be round, since the
# mol/s <-> year conversion factor isn't itself a power of ten).
sec_per_year  <- 365.25 * 24 * 3600
tau_max_years <- 1e10
k_floor_T     <- 1 / (tau_max_years * sec_per_year)
k_max_T       <- max(rate_data_T$k, na.rm = TRUE)
tau_min_years <- 1 / (k_max_T * sec_per_year)
tau_breaks_T  <- 10^seq(floor(log10(tau_min_years)), ceiling(log10(tau_max_years)))
y_breaks_T    <- sort(1 / (tau_breaks_T * sec_per_year))

p_k_T <- ggplot(rate_data_T, aes(x = TC, y = k, colour = mineral, linetype = mineral)) +
    geom_line(linewidth = 1.0) +
    scale_linetype_manual(values = mineral_linetype_map) +
    scale_y_log10(
        limits = c(k_floor_T, NA),
        breaks = y_breaks_T,
        labels = scales::label_scientific(),
        sec.axis = sec_axis(
            transform = ~ 1 / (. * sec_per_year),
            name   = expression(tau==1/k[total]~~"(yr)"),
            breaks = tau_breaks_T,
            labels = scales::label_scientific()
        )
    ) +
    scale_color_manual(values=mineral_color_map) +
    labs(
        x       = "Temperature (°C)",
        y       = expression(k[total]~~"(mol" ~~ s^{-1} * ")"),
        title   = sprintf("Mineral & IOM dissolution/decomposition rate constants vs temperature  [pH = %.1f]", pH_val),
        caption = sprintf("k assumes 1 mol present (mol/s, directly comparable).  Pyrrhotite Fe³⁺ fixed at %.0e.", fe3_val)
    ) +
    theme_minimal(base_size = 11) +
    theme(
        aspect.ratio      = 0.65,
        axis.line         = element_line(colour = "grey40"),
        axis.ticks        = element_line(colour = "grey40"),
        axis.ticks.length = unit(0.15, "cm"),
        plot.caption      = element_text(size = 7, colour = "grey45")
    )

ggsave(file.path(output_dir, "mineral_rate_constants_T.png"), plot = p_k_T, dpi = 300, width = 10, height = 5.5)

invisible(list(k_vs_pH = p_k_pH, k_vs_T = p_k_T))
}

# Run automatically only when this file is executed as a top-level script
# (Rscript, or "Source" in RStudio/R console) -- never true inside knitr, so
# notebooks must call plot_mineral_rates() explicitly after sourcing.
if (sys.nframe() == 0) {
  plot_mineral_rates()
}