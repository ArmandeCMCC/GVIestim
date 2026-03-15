# export_paris_appendix_effectmod_tables.R
#
# Build Overleaf-ready appendix tables:
# - one table per heat metric
# - three greenness columns (GVI, NDVI, IMU total)
# - rows for primary (p10->p90) and companion (min->max) contrasts at p99
#
# Inputs (run-tagged):
#   paris_effectmod_heat_green_summary_<tag>.csv
#   paris_effectmod_paper_table_p99_minmax_<tag>.csv
#
# Outputs:
#   tables/appendix/effectmod_<heat>_<tag>.tex
#   tables/appendix/effectmod_tables_all_<tag>.tex
#   tables/appendix/effectmod_appendix_table_source_<tag>.csv

library(data.table)

base_dir <- "/Users/armandeaboudrar-meda/Desktop/GVIestim/paper"

norm_tag <- function(x) {
  x <- trimws(x)
  if (!nzchar(x)) return("")
  if (!startsWith(x, "_")) x <- paste0("_", x)
  x
}

run_tag <- norm_tag(Sys.getenv("RUN_TAG", "rerun_rule1_sum_ghsheat_ghsgreen"))

in_summary <- file.path(base_dir, paste0("paris_effectmod_heat_green_summary", run_tag, ".csv"))
in_minmax  <- file.path(base_dir, paste0("paris_effectmod_paper_table_p99_minmax", run_tag, ".csv"))

out_dir <- file.path(base_dir, "tables", "appendix")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_all_tex <- file.path(out_dir, paste0("effectmod_tables_all", run_tag, ".tex"))
out_src_csv <- file.path(out_dir, paste0("effectmod_appendix_table_source", run_tag, ".csv"))

if (!file.exists(in_summary)) stop("Missing summary file: ", in_summary)
if (!file.exists(in_minmax))  stop("Missing minmax file: ", in_minmax)

summary_dt <- fread(in_summary)
minmax_dt  <- fread(in_minmax)

main_green <- c("gvi_popw_points", "ndvi_popw_ghs", "imu_veg_total_popw_ghs")
main_heat <- c(
  "t2m_min", "t2m_mean", "t2m_max",
  "lst_min", "lst_mean", "lst_max",
  "wbgt_min", "wbgt_mean", "wbgt_max"
)

green_labels <- c(
  gvi_popw_points = "GVI (point pop-w., GHS)",
  ndvi_popw_ghs = "NDVI (pop-w., GHS)",
  imu_veg_total_popw_ghs = "IMU total (pop-w., GHS)"
)

heat_labels <- c(
  t2m_min = "T2M min",
  t2m_mean = "T2M mean",
  t2m_max = "T2M max",
  lst_min = "LST min",
  lst_mean = "LST mean",
  lst_max = "LST max",
  wbgt_min = "WBGT min",
  wbgt_mean = "WBGT mean",
  wbgt_max = "WBGT max"
)

esc_latex <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([_%&#$])", "\\\\\\1", x, perl = TRUE)
  x
}

fmt_est <- function(est, lo, hi, digits = 2) {
  ok <- is.finite(est) & is.finite(lo) & is.finite(hi)
  out <- rep("NA", length(est))
  if (any(ok)) {
    out[ok] <- sprintf(
      paste0("%.", digits, "f (", "%.", digits, "f, ", "%.", digits, "f)"),
      est[ok], lo[ok], hi[ok]
    )
  }
  out
}

primary_dt <- summary_dt[
  heat_metric %in% main_heat & green_metric_raw %in% main_green,
  .(
    heat_metric,
    green_metric_raw,
    row_label = "Primary (p10$\\rightarrow$p90)",
    value = fmt_est(
      as.numeric(cr_mod_p99_p10p90),
      as.numeric(cr_mod_p99_p10p90_low),
      as.numeric(cr_mod_p99_p10p90_high)
    ),
    heat_mmt = as.numeric(heat_mmt),
    heat_p99 = as.numeric(heat_p99)
  )
]

companion_dt <- minmax_dt[
  heat_metric %in% main_heat & green_metric_raw %in% main_green,
  .(
    heat_metric,
    green_metric_raw,
    row_label = "Secondary (min$\\rightarrow$max)",
    value = fmt_est(
      as.numeric(cr_mod_paper),
      as.numeric(cr_mod_paper_low),
      as.numeric(cr_mod_paper_high)
    ),
    heat_mmt = as.numeric(heat_mmt),
    heat_p99 = as.numeric(heat_p99)
  )
]

plot_dt <- rbindlist(list(primary_dt, companion_dt), fill = TRUE)

# Save source used for appendix tables for auditability
fwrite(plot_dt, out_src_csv)

build_one_table <- function(h) {
  sub <- plot_dt[heat_metric == h]
  if (nrow(sub) == 0) return(NULL)

  sub[, row_order := factor(
    row_label,
    levels = c("Primary (p10$\\rightarrow$p90)", "Secondary (min$\\rightarrow$max)")
  )]
  setorder(sub, row_order)
  sub[, row_order := NULL]
  wide <- dcast(sub, row_label ~ green_metric_raw, value.var = "value")
  wide[, row_order := factor(
    row_label,
    levels = c("Primary (p10$\\rightarrow$p90)", "Secondary (min$\\rightarrow$max)")
  )]
  setorder(wide, row_order)
  wide[, row_order := NULL]

  # Ensure all three expected columns exist
  for (g in main_green) {
    if (!(g %in% names(wide))) wide[, (g) := "NA"]
  }
  setcolorder(wide, c("row_label", main_green))

  hmmt <- unique(na.omit(sub$heat_mmt))
  hp99 <- unique(na.omit(sub$heat_p99))
  hmmt_txt <- if (length(hmmt) > 0) sprintf("%.2f", hmmt[1]) else "NA"
  hp99_txt <- if (length(hp99) > 0) sprintf("%.2f", hp99[1]) else "NA"

  h_lab <- heat_labels[[h]]
  cap <- paste0(
    "Effect-modification estimates (\\% change in RR) for ", h_lab,
    " at p99 (MMT = ", hmmt_txt, ", p99 = ", hp99_txt, ")."
  )

  lines <- c(
    "\\begin{table}[H]",
    "\\centering",
    paste0("\\caption{", cap, "}"),
    paste0("\\label{tab:effectmod_", h, esc_latex(run_tag), "}"),
    "\\small",
    "\\begin{tabular}{lccc}",
    "\\toprule",
    paste(
      "Contrast",
      green_labels[main_green[1]],
      green_labels[main_green[2]],
      green_labels[main_green[3]],
      sep = " & "
    ),
    "\\\\",
    "\\midrule"
  )

  for (i in seq_len(nrow(wide))) {
    lines <- c(
      lines,
      paste(
        wide$row_label[i],
        wide[[main_green[1]]][i],
        wide[[main_green[2]]][i],
        wide[[main_green[3]]][i],
        sep = " & "
      ),
      "\\\\"
    )
  }

  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    ""
  )

  lines
}

all_lines <- character()
for (h in main_heat) {
  tb <- build_one_table(h)
  if (is.null(tb)) next
  out_tex <- file.path(out_dir, paste0("effectmod_", h, run_tag, ".tex"))
  writeLines(tb, out_tex)
  all_lines <- c(all_lines, tb)
}

writeLines(all_lines, out_all_tex)

cat("\nSaved appendix table source:\n - ", out_src_csv, "\n", sep = "")
cat("Saved combined LaTeX tables:\n - ", out_all_tex, "\n", sep = "")
cat("Saved one-file-per-heat tables in:\n - ", out_dir, "\n", sep = "")
