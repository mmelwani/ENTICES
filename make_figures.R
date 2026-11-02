# Figures from production runs (production kit, item 4; PRODUCTION_GUIDE.md). From the repo root:
#
#   Rscript make_figures.R                          (every runs/*/summary.csv; all six figures to figures/)
#   Rscript make_figures.R reference k_1e-11 iom_3  (these cases only)
#   Rscript make_figures.R --species pH,DIC_C4,H2,Na --output-var ocean_DIC_C4 --runs-dir runs --out figures
#
# Or source("make_figures.R") in R and call the figure functions, one per figure, each returning a ggplot:
#   fig_skin_vs_k(), fig_ocean(), fig_org_minus_blank(), fig_rock(), fig_tornado(), fig_h2_budget().
# Every function takes `summaries` (from read_summaries()), `cases` (names; NULL = all) and `overlay`: an
# optional data frame to draw on top, e.g. observational bands or points that you supply. Its columns:
#   y, or ymin and ymax (in the panel's units; ymin/ymax without x draws a horizontal band across the panel);
#   x, or xmin and xmax (optional; time in yr, or k in m2 for fig_skin_vs_k; a band without xmin spans the panel);
#   `panel` (optional: draw only on that panel; the species code for fig_ocean and fig_org_minus_blank, e.g. "pH"
#   or "DIC_C4"; phi, k_m2 or delta_m for fig_rock; delta_m, skin_T_C or circulation_yr for fig_skin_vs_k;
#   h2_generated_cum, h2_present, h2_ocean_dissolved, ... for fig_h2_budget); `label` (optional text, placed at
#   x (or xmin) and y (or ymax)).
# No observational data are included here.
suppressPackageStartupMessages({ library(ggplot2) })

# The reference is near-black; other cases take the validated categorical slots in a fixed order (the dataviz
# reference palette: CVD-checked for adjacent lines), the same colour for a case in every figure. Organic
# runs are solid lines, blanks dashed.
FIG_REFERENCE_COLOUR <- "#0b0b0b"
FIG_SLOTS <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")

# summary.csv of every run under runs_dir (or the named cases), bound into one data frame.
read_summaries <- function(runs_dir = "runs", cases = NULL) {
  f <- Sys.glob(file.path(runs_dir, "*", "summary.csv"))
  nm <- basename(dirname(f))
  if (!is.null(cases)) { miss <- setdiff(cases, nm); if (length(miss)) stop("no summary.csv for: ", paste(miss, collapse = ", ")); f <- f[nm %in% cases] }
  if (!length(f)) stop("no summary.csv found under ", runs_dir)
  s <- lapply(f, utils::read.csv, stringsAsFactors = FALSE)
  cols <- Reduce(union, lapply(s, names))
  do.call(rbind, lapply(s, function(d) { for (cc in setdiff(cols, names(d))) d[[cc]] <- NA; d[cols] }))
}

# One colour per case: the reference near-black, the others in slot order (sorted names, so the mapping does not
# depend on which figure or selection). More than 8 non-reference cases in one figure: an error (use fewer
# cases per figure, or the tornado figure, which needs no colours).
case_colours <- function(cases) {
  cases <- unique(cases); others <- sort(setdiff(cases, "reference"))
  if (length(others) > length(FIG_SLOTS))
    stop(sprintf("%d cases in one figure: at most %d besides the reference can be told apart; pass fewer cases (fig_tornado() shows all of them)",
                 length(others), length(FIG_SLOTS)))
  c(if ("reference" %in% cases) c(reference = FIG_REFERENCE_COLOUR), stats::setNames(FIG_SLOTS[seq_along(others)], others))
}

fig_theme <- function() theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
        axis.line = element_line(colour = "grey40", linewidth = 0.3), strip.text = element_text(face = "bold", hjust = 0),
        legend.position = "bottom", plot.title = element_text(face = "bold"))

log_lab <- scales::label_scientific(digits = 2)
log_x_time <- function() scale_x_log10("Time (yr)", breaks = scales::breaks_log(6), labels = log_lab)

pick <- function(summaries, cases) {
  if (is.null(cases)) return(summaries)
  miss <- setdiff(cases, unique(summaries$case)); if (length(miss)) stop("not in the summaries: ", paste(miss, collapse = ", "))
  summaries[summaries$case %in% cases, , drop = FALSE]
}

# Draw an overlay data frame (see the header) on a plot. panels: the facet levels, named by the codes an
# overlay's `panel` may use (e.g. c(pH = "pH", Na = "Na (mol/kgw)")); rows without a panel go on every panel.
add_overlay <- function(p, overlay, panels = NULL, log_x = TRUE) {
  if (is.null(overlay) || !nrow(overlay)) return(p)
  ov <- overlay
  if (!is.null(panels)) {
    if (!"panel" %in% names(ov)) ov$panel <- NA
    lab <- ifelse(ov$panel %in% names(panels), panels[as.character(ov$panel)], as.character(ov$panel))
    bad <- setdiff(stats::na.omit(lab), panels)
    if (length(bad)) stop("overlay panel(s) not in this figure: ", paste(bad, collapse = ", "), " (use ", paste(names(panels), collapse = ", "), ")")
    every <- is.na(lab)
    ov <- rbind(transform(ov[!every, , drop = FALSE], panel = lab[!every]),
                do.call(rbind, lapply(panels, function(pp) transform(ov[every, , drop = FALSE], panel = rep(pp, sum(every))))))
    ov$panel <- factor(ov$panel, levels = panels)
  }
  band <- all(c("ymin", "ymax") %in% names(ov))
  if (band) {
    b <- ov[!is.na(ov$ymin) & !is.na(ov$ymax), , drop = FALSE]
    if (!"xmin" %in% names(b)) b$xmin <- NA
    if (!"xmax" %in% names(b)) b$xmax <- NA
    b$xmin[is.na(b$xmin)] <- if (log_x) 0 else -Inf; b$xmax[is.na(b$xmax)] <- Inf   # 0 is -Inf after a log10 axis
    p <- p + geom_rect(data = b, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE,
                       fill = "grey50", alpha = 0.15, colour = NA)
  }
  if (all(c("x", "y") %in% names(ov))) {
    pt <- ov[!is.na(ov$x) & !is.na(ov$y), , drop = FALSE]
    p <- p + geom_point(data = pt, aes(x = x, y = y), inherit.aes = FALSE, shape = 23, size = 2.5, fill = "white", colour = "grey20")
  }
  if ("label" %in% names(ov)) {
    lb <- ov[!is.na(ov$label), , drop = FALSE]
    if (nrow(lb)) {
      lb$lx <- if ("x" %in% names(lb)) lb$x else if ("xmin" %in% names(lb)) lb$xmin else NA
      if ("x" %in% names(lb) && "xmin" %in% names(lb)) lb$lx[is.na(lb$lx)] <- lb$xmin[is.na(lb$lx)]
      lb$ly <- if ("y" %in% names(lb)) lb$y else lb$ymax
      if ("y" %in% names(lb) && "ymax" %in% names(lb)) lb$ly[is.na(lb$ly)] <- lb$ymax[is.na(lb$ly)]
      lb <- lb[is.finite(lb$lx) & is.finite(lb$ly), , drop = FALSE]   # a label needs a position (x or xmin, y or ymax)
      if (nrow(lb)) p <- p + geom_text(data = lb, aes(x = lx, y = ly, label = label), inherit.aes = FALSE, size = 2.6, hjust = -0.05, vjust = -0.4, colour = "grey20")
    }
  }
  p
}

# The species fig_ocean() can show: the ocean_* columns of summary.csv without the prefix.
FIG_SPECIES_LABELS <- c(pH = "pH", pe = "pe", DIC_C4 = "DIC, C(+4) (mol/kgw)", Mtg = "CH4 (Mtg) (mol/kgw)", H2 = "H2 (mol/kgw)",
                        NH4_NH3 = "NH4+ + NH3 (mol/kgw)", SO4 = "SO4-2 (mol/kgw)", S_minus2 = "HS- + H2S (mol/kgw)",
                        Na = "Na (mol/kgw)", Cl = "Cl (mol/kgw)", K = "K (mol/kgw)", Ca = "Ca (mol/kgw)", Mg = "Mg (mol/kgw)",
                        Si = "Si (mol/kgw)", P = "P (mol/kgw)", ionic_strength = "Ionic strength (mol/kgw)")
FIG_LINEAR_SPECIES <- c("pH", "pe")

# 1. Skin depth, skin temperature and circulation time against permeability, each case's trajectory (k falls
# as the pores fill; the start marked with a circle).
fig_skin_vs_k <- function(summaries, cases = NULL, overlay = NULL) {
  s <- pick(summaries, cases); cols <- case_colours(s$case)
  s <- s[order(s$case, s$time_yr), ]
  start <- s[!duplicated(s$case), ]
  one <- function(col, title, log_y) {
    p <- ggplot(s, aes(k_m2, .data[[col]], colour = case)) + geom_path(linewidth = 0.6) +
      geom_point(data = start, size = 2.2, shape = 21, fill = "white", stroke = 0.8) +
      scale_x_log10(expression(Permeability ~ k ~ (m^2)), breaks = scales::breaks_log(5), labels = log_lab) +
      (if (log_y) scale_y_log10(NULL, breaks = scales::breaks_log(6), labels = log_lab) else scale_y_continuous(NULL)) +
      scale_colour_manual("Case", values = cols) + fig_theme() + labs(subtitle = title)
    add_overlay(p, if (is.null(overlay) || !"panel" %in% names(overlay)) overlay else overlay[overlay$panel %in% col, , drop = FALSE])
  }
  patchwork::wrap_plots(one("delta_m", "Skin depth (m)", TRUE), one("skin_T_C", "Skin temperature (degC)", FALSE),
                        one("circulation_yr", "Ocean circulation time (yr)", TRUE), ncol = 3, guides = "collect") +
    patchwork::plot_annotation(title = "Skin depth, skin temperature and circulation time against permeability",
                               subtitle = "Each case's trajectory as k falls with porosity; circle: t = 0") &
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
}

# 2. Ocean time series of the chosen species, organic (solid) and blank (dashed), one panel per species.
fig_ocean <- function(summaries, cases = NULL, species = c("pH", "DIC_C4", "H2", "Na"), blank = TRUE, overlay = NULL) {
  s <- pick(summaries, cases); cols <- case_colours(s$case)
  bad <- setdiff(species, names(FIG_SPECIES_LABELS)); if (length(bad)) stop("unknown species: ", paste(bad, collapse = ", "), " (choose from ", paste(names(FIG_SPECIES_LABELS), collapse = ", "), ")")
  long <- do.call(rbind, lapply(species, function(sp) {
    o <- data.frame(case = s$case, time_yr = s$time_yr, panel = sp, run = "organic", value = s[[paste0("ocean_", sp)]])
    b <- if (blank && paste0("blank_ocean_", sp) %in% names(s)) data.frame(case = s$case, time_yr = s$time_yr, panel = sp, run = "blank", value = s[[paste0("blank_ocean_", sp)]]) else NULL
    rbind(o, b)
  }))
  long <- long[is.finite(long$value) & long$time_yr > 0, ]
  long$panel <- factor(long$panel, levels = species, labels = FIG_SPECIES_LABELS[species])
  p <- ggplot(long, aes(time_yr, value, colour = case, linetype = run, group = interaction(case, run))) +
    geom_line(linewidth = 0.6) + facet_wrap(~panel, scales = "free_y") + log_x_time() +
    scale_y_continuous(NULL) + scale_colour_manual("Case", values = cols) +
    scale_linetype_manual("Run", values = c(organic = "solid", blank = "22")) + fig_theme() +
    labs(title = "Ocean composition", subtitle = "Concentrations on a linear axis per panel; pH and pe unitless")
  add_overlay(p, overlay, stats::setNames(unname(FIG_SPECIES_LABELS[species]), species))
}

# 3. The organics' contribution: organic minus blank ocean carbon (DIC as C(+4), and methane) against time.
fig_org_minus_blank <- function(summaries, cases = NULL, species = c("DIC_C4", "Mtg"), overlay = NULL) {
  s <- pick(summaries, cases); cols <- case_colours(s$case)
  need <- paste0("org_minus_blank_ocean_", species); miss <- setdiff(need, names(s))
  if (length(miss)) stop("no organic minus blank columns (", paste(miss, collapse = ", "), "): were the blanks run?")
  long <- do.call(rbind, lapply(species, function(sp) data.frame(case = s$case, time_yr = s$time_yr, panel = sp, value = s[[paste0("org_minus_blank_ocean_", sp)]])))
  long <- long[is.finite(long$value) & long$time_yr > 0, ]
  long$panel <- factor(long$panel, levels = species, labels = paste("Organic - blank:", FIG_SPECIES_LABELS[species]))
  p <- ggplot(long, aes(time_yr, value, colour = case)) + geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_line(linewidth = 0.6) + facet_wrap(~panel, scales = "free_y") + log_x_time() +
    scale_y_continuous(NULL) + scale_colour_manual("Case", values = cols) + fig_theme() +
    labs(title = "The organics' contribution to ocean carbon", subtitle = "Organic run minus its same-rock blank (mol/kgw)")
  add_overlay(p, overlay, stats::setNames(paste("Organic - blank:", FIG_SPECIES_LABELS[species]), species))
}

# 4. Porosity, permeability and skin depth against time (clogging and k_floor marked).
fig_rock <- function(summaries, cases = NULL, overlay = NULL) {
  s <- pick(summaries, cases); cols <- case_colours(s$case)
  long <- rbind(data.frame(case = s$case, time_yr = s$time_yr, panel = "Porosity", value = s$phi),
                data.frame(case = s$case, time_yr = s$time_yr, panel = "Permeability (m2)", value = s$k_m2),
                data.frame(case = s$case, time_yr = s$time_yr, panel = "Skin depth (m)", value = s$delta_m))
  long$panel <- factor(long$panel, levels = c("Porosity", "Permeability (m2)", "Skin depth (m)"))
  clog <- unique(s[is.finite(s$clog_time_yr), c("case", "clog_time_yr")])
  p <- ggplot(long[long$time_yr > 0, ], aes(time_yr, value, colour = case)) + geom_line(linewidth = 0.6) +
    facet_wrap(~panel, scales = "free_y", ncol = 3) + log_x_time() + scale_y_log10(NULL, breaks = scales::breaks_log(6), labels = log_lab) +
    scale_colour_manual("Case", values = cols) + fig_theme() +
    labs(title = "Porosity, permeability and skin depth", subtitle = if (nrow(clog)) "Dotted lines: clog time (phi reaches phi_min)" else NULL)
  if (nrow(clog)) p <- p + geom_vline(data = clog, aes(xintercept = clog_time_yr, colour = case), linetype = "dotted", linewidth = 0.4, show.legend = FALSE)
  add_overlay(p, overlay, c(phi = "Porosity", k_m2 = "Permeability (m2)", delta_m = "Skin depth (m)"))
}

# 5. Sweep sensitivity: each case's end value of one output against the reference's (tornado style): bars of
# log10(case / reference) for concentrations, or the difference for pH/pe. Cases that did not end at steady
# state are marked with their end reason.
fig_tornado <- function(summaries, output = "ocean_DIC_C4", cases = NULL, reference = "reference", overlay = NULL) {
  s <- pick(summaries, cases)
  if (!output %in% names(s)) stop("no column ", output, " in the summaries")
  if (!reference %in% s$case) stop("the reference case \"", reference, "\" is not among the summaries")
  last <- do.call(rbind, lapply(split(s, s$case), function(d) d[which.max(d$time_yr), ]))
  ref <- last[last$case == reference, output]
  log_scale <- !(sub("^(org_minus_blank_|blank_)?ocean_", "", output) %in% FIG_LINEAR_SPECIES) && all(last[[output]] > 0, na.rm = TRUE)
  d <- last[last$case != reference, ]
  d$effect <- if (log_scale) log10(d[[output]] / ref) else d[[output]] - ref
  d$note <- ifelse(d$end_reason == "steady state", "", paste0(" (", d$end_reason, ")"))
  d$case_lab <- paste0(d$case, d$note)
  d <- d[order(abs(d$effect)), ]; d$case_lab <- factor(d$case_lab, levels = d$case_lab)
  d$sign <- ifelse(d$effect >= 0, "above the reference", "below the reference")
  p <- ggplot(d, aes(effect, case_lab, fill = sign)) + geom_col(width = 0.7) + geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.4) +
    scale_fill_manual(NULL, values = c(`above the reference` = "#2a78d6", `below the reference` = "#eb6834")) +
    scale_x_continuous(if (log_scale) sprintf("log10(case / reference), %s", output) else sprintf("case - reference, %s", output)) +
    scale_y_discrete(NULL) + fig_theme() +
    labs(title = sprintf("Sensitivity of %s at the end of each run", output),
         subtitle = sprintf("Reference: %.4g at t = %.3g yr (%s). In brackets: runs not ended by steady state",
                            ref, last$time_yr[last$case == reference], last$end_reason[last$case == reference]))
  add_overlay(p, overlay, log_x = FALSE)
}

# 6. The H2 budget against time, scaled to the whole core (mol): generated (cumulative), present, and where it is.
fig_h2_budget <- function(summaries, cases = NULL, overlay = NULL) {
  s <- pick(summaries, cases); cols <- case_colours(s$case)
  comp <- c(h2_generated_cum_core_mol = "Generated (cumulative)", h2_present_core_mol = "Present (all reservoirs)",
            h2_ocean_dissolved_core_mol = "Ocean, dissolved", h2_ocean_gas_core_mol = "Ocean gas",
            h2_porewater_dissolved_core_mol = "Porewater, dissolved", h2_pore_gas_core_mol = "Pore gas")
  comp <- comp[names(comp) %in% names(s)]
  long <- do.call(rbind, lapply(names(comp), function(cc) data.frame(case = s$case, time_yr = s$time_yr, panel = comp[[cc]], value = s[[cc]])))
  long <- long[long$time_yr > 0 & is.finite(long$value) & long$value > 0, ]
  long$panel <- factor(long$panel, levels = unname(comp))
  p <- ggplot(long, aes(time_yr, value, colour = case)) + geom_line(linewidth = 0.6) + facet_wrap(~panel, ncol = 3) +
    log_x_time() + scale_y_log10("H2 (mol, whole core)", breaks = scales::breaks_log(6), labels = log_lab) +
    scale_colour_manual("Case", values = cols) + fig_theme() +
    labs(title = "Hydrogen budget, scaled to the whole core", subtitle = "Positive values only on the log axis")
  add_overlay(p, overlay, stats::setNames(unname(comp), sub("_core_mol$", "", names(comp))))
}

# Write a figure as PDF and PNG.
save_figure <- function(p, name, out = "figures", width = 7.2, height = 4.2) {
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(out, paste0(name, ".pdf")), p, width = width, height = height, device = grDevices::pdf)
  ggsave(file.path(out, paste0(name, ".png")), p, width = width, height = height, dpi = 300, bg = "white")
  invisible(file.path(out, paste0(name, c(".pdf", ".png"))))
}

# All six figures for a set of cases.
make_all_figures <- function(summaries, cases = NULL, species = c("pH", "DIC_C4", "H2", "Na"), output = "ocean_DIC_C4", out = "figures") {
  line_cases <- if (is.null(cases)) unique(summaries$case) else cases
  if (length(setdiff(line_cases, "reference")) > length(FIG_SLOTS)) {
    line_cases <- c(intersect("reference", line_cases), utils::head(sort(setdiff(line_cases, "reference")), length(FIG_SLOTS)))
    message("more than 8 cases: the line figures show ", paste(line_cases, collapse = ", "), "; the tornado shows all")
  }
  f <- c(save_figure(fig_skin_vs_k(summaries, line_cases), "fig1_skin_vs_k", out, height = 3.4),
         save_figure(fig_ocean(summaries, line_cases, species), "fig2_ocean", out, height = 5),
         if (any(grepl("^org_minus_blank_", names(summaries)))) save_figure(fig_org_minus_blank(summaries, line_cases), "fig3_org_minus_blank", out, height = 3.4),
         save_figure(fig_rock(summaries, line_cases), "fig4_rock", out, height = 3.4),
         if ("reference" %in% summaries$case && length(unique(summaries$case)) > 1)
           save_figure(fig_tornado(summaries, output, cases), "fig5_tornado", out, height = 1.2 + 0.28 * length(unique(summaries$case))),
         save_figure(fig_h2_budget(summaries, line_cases), "fig6_h2_budget", out, height = 4.6))
  invisible(f)
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  opt <- function(flag, default) { i <- match(flag, args); if (is.na(i)) default else args[i + 1] }
  flags <- c("--species", "--output-var", "--runs-dir", "--out")
  cases <- args[!(args %in% flags) & !(seq_along(args) - 1) %in% match(flags, args)]
  s <- read_summaries(opt("--runs-dir", "runs"), if (length(cases)) cases else NULL)
  f <- make_all_figures(s, NULL, species = strsplit(opt("--species", "pH,DIC_C4,H2,Na"), ",")[[1]],
                        output = opt("--output-var", "ocean_DIC_C4"), out = opt("--out", "figures"))
  cat("wrote", length(f), "files to", opt("--out", "figures"), "\n")
}
