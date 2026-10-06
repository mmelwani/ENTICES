## Compares organics-on vs organics-off final state for the matched ENTICES
## case (enceladus_k1.00e-12_d600_t7.50e+04_circ1). Run from this directory,
## or `Rscript phreeqc_inputs/redox_budget/compare.R` from the repo root.
this_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
if (!file.exists(file.path(this_dir, "org_on.csv"))) this_dir <- "phreeqc_inputs/redox_budget"
setwd(this_dir)

rd <- function(f) {
  d <- read.delim(f, check.names = FALSE, strip.white = TRUE)
  names(d) <- trimws(names(d))
  d[nrow(d), , drop = FALSE]      # final state
}
on  <- rd("org_on.csv")
off <- rd("org_off.csv")

cat("=== REDOX INDICATORS (final state) ===\n")
key <- c("pH","pe","m_H2","m_CH4","m_CO2","m_HCO3-","m_CO3-2","m_HS-","m_SO4-2",
         "m_NH4+","m_NH3","m_N2","m_Fe+2","m_Fe+3","C","S","N","Fe","H","O")
cat(sprintf("%-12s %14s %14s %12s\n","quantity","organics OFF","organics ON","ratio on/off"))
for (k in key) {
  if (!(k %in% names(on))) next
  a <- suppressWarnings(as.numeric(off[[k]])); b <- suppressWarnings(as.numeric(on[[k]]))
  if (is.na(a) || is.na(b)) next
  r <- if (a != 0) b/a else NA
  cat(sprintf("%-12s %14.6g %14.6g %12s\n", k, a, b,
              if (is.na(r)) "--" else sprintf("%.4g", r)))
}

cat("\n=== KINETIC PRIMARY MINERALS (moles remaining, k_ prefix) ===\n")
for (base in c("Forsterite","Fayalite","Enstatite","Ferrosilite","Pyrrhotite","Anorthite","Albite","Magnetite")) {
  k <- paste0("k_", base)
  if (!(k %in% names(on))) { cat(sprintf("%-14s (column not found)\n", base)); next }
  a <- as.numeric(off[[k]]); b <- as.numeric(on[[k]])
  pct <- if (a != 0) (b/a-1)*100 else NA
  cat(sprintf("%-14s off=%14.8g  on=%14.8g  delta=%+14.8g (%+.4f%%)\n", base, a, b, b-a, pct))
}

cat("\n=== SECONDARY PHASES: all with a non-trivial difference ===\n")
dcols <- grep("^d_", names(on), value = TRUE)
rows <- list()
for (k in dcols) {
  a <- suppressWarnings(as.numeric(off[[k]])); b <- suppressWarnings(as.numeric(on[[k]]))
  if (is.na(a) || is.na(b)) next
  if (max(abs(a), abs(b)) < 1e-12) next
  rows[[length(rows)+1]] <- data.frame(phase = sub("^d_","",k), off = a, on = b,
                                       diff = b - a, stringsAsFactors = FALSE)
}
if (length(rows)) {
  df <- do.call(rbind, rows)
  df <- df[order(-abs(df$diff)), ]
  for (i in seq_len(nrow(df)))
    cat(sprintf("%-20s off=%12.6g  on=%12.6g  delta=%+12.6g\n",
                df$phase[i], df$off[i], df$on[i], df$diff[i]))
} else cat("(none above 1e-12)\n")
