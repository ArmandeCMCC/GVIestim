# plot_paris_main_insight_figures.R
#
# Main paper figures: compact exposure-response curves, direct
# contrasts (GVI vs NDVI, GVI vs IMU total), and burden heatmap
# Default RUN_TAG: _rerun_rule1_sum_ghsheat_ghsgreen.
#
# Inputs (RUN_TAG-specific):
#   - paris_effectmod_results_harmonized{RUN_TAG}.csv
#   - paris_effectmod_heat_green_curves{RUN_TAG}.csv
#   - paris_attr_scenario_keytable{RUN_TAG}.csv
#

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

`%||%` <- function(x, y) if (is.null(x)) y else x
args_all <- commandArgs(trailingOnly = FALSE)
file_flag <- "--file="
script_arg <- args_all[grepl(file_flag, args_all)]
script_path <- if (length(script_arg) > 0) sub(file_flag, "", script_arg[1]) else "."
this_dir <- normalizePath(dirname(script_path), mustWork = FALSE)
setwd(this_dir)

tag_name <- function(base, run_tag) {
  if (nzchar(run_tag)) {
    stem <- sub("\\.[^.]+$", "", base)
    ext <- sub("^.*(\\.[^.]+)$", "\\1", base)
    if (identical(stem, ext)) ext <- ""
    paste0(stem, run_tag, ext)
  } else {
    base
  }
}

save_plot_dual <- function(plot_obj, png_path, pdf_path, width, height, dpi = 320) {
  ggsave(png_path, plot_obj, width = width, height = height, dpi = dpi)
  ggsave(pdf_path, plot_obj, width = width, height = height)
}

# config
run_tag_raw <- Sys.getenv("RUN_TAG", "rerun_rule1_sum_ghsheat_ghsgreen")
run_tag <- if (nzchar(run_tag_raw)) {
  if (startsWith(run_tag_raw, "_")) run_tag_raw else paste0("_", run_tag_raw)
} else {
  ""
}

fig_dir <- file.path(getwd(), "figs")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

harm_file <- tag_name("paris_effectmod_results_harmonized.csv", run_tag)
curves_file <- tag_name("paris_effectmod_heat_green_curves.csv", run_tag)
burden_file <- tag_name("paris_attr_scenario_keytable.csv", run_tag)

if (!file.exists(harm_file)) stop("Missing harmonized file: ", harm_file)
if (!file.exists(curves_file)) stop("Missing curves file: ", curves_file)
if (!file.exists(burden_file)) stop("Missing burden key table: ", burden_file)

harm <- fread(harm_file)
curves <- fread(curves_file)
burden <- fread(burden_file)

heat_labels <- c(
  t2m_min = "T2M min", t2m_mean = "T2M mean", t2m_max = "T2M max",
  lst_min = "LST min", lst_mean = "LST mean", lst_max = "LST max",
  wbgt_min = "WBGT min", wbgt_mean = "WBGT mean", wbgt_max = "WBGT max"
)
green_labels <- c(
  gvi_popw_points = "GVI (point pop-w., GHS)",
  ndvi_popw_ghs = "NDVI (pop-w., GHS)",
  imu_veg_total_popw_ghs = "IMU total (pop-w., GHS)",
  imu_high_popw_ghs = "IMU high (pop-w., GHS)",
  imu_low_popw_ghs = "IMU low (pop-w., GHS)"
)

# 1) compact curves
build_compact_curves <- function(heat_set, out_stem, subtitle_txt) {
  green_set <- c("gvi_popw_points", "ndvi_popw_ghs", "imu_veg_total_popw_ghs")
  level_set <- c("low_p10", "high_p90")

  dt <- copy(curves)[
    heat_metric %in% heat_set &
      green_metric_raw %in% green_set &
      green_level %in% level_set
  ]
  if (nrow(dt) == 0) stop("No curve rows for requested set: ", out_stem)

  dt[, heat_label := factor(heat_labels[heat_metric], levels = heat_labels[heat_set])]
  dt[, green_label := factor(green_labels[green_metric_raw], levels = green_labels[green_set])]
  dt[, green_level := factor(green_level, levels = level_set, labels = c("Low greenness (p10)", "High greenness (p90)"))]

  guides_dt <- unique(dt[, .(heat_metric, green_metric_raw, heat_mmt, heat_p99, heat_label, green_label)])

  p <- ggplot(dt, aes(x = temp, y = rr, linetype = green_level, group = green_level)) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.3) +
    geom_vline(data = guides_dt, aes(xintercept = heat_mmt), linetype = "dotted", linewidth = 0.3, show.legend = FALSE) +
    geom_vline(data = guides_dt, aes(xintercept = heat_p99), linetype = "dotdash", linewidth = 0.3, show.legend = FALSE) +
    geom_ribbon(aes(ymin = rr_low, ymax = rr_high, fill = green_level), alpha = 0.14, colour = NA) +
    geom_line(linewidth = 0.7) +
    facet_grid(green_label ~ heat_label, scales = "free_x") +
    labs(
      title = "Compact exposure-response comparison",
      subtitle = subtitle_txt,
      x = "Heat metric value",
      y = "Relative risk (RR)",
      linetype = NULL,
      fill = NULL
    ) +
    theme_bw(base_size = 10.5) +
    theme(
      strip.background = element_rect(fill = "grey95", colour = "grey75"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )

  out_png <- file.path(fig_dir, tag_name(paste0(out_stem, ".png"), run_tag))
  out_pdf <- file.path(fig_dir, tag_name(paste0(out_stem, ".pdf"), run_tag))
  save_plot_dual(p, out_png, out_pdf, width = 10.8, height = 8.6)
  c(out_png, out_pdf)
}

curve_main_paths <- build_compact_curves(
  heat_set = c("wbgt_mean", "t2m_mean"),
  out_stem = "paris_curves_effectmod_compact_main",
  subtitle_txt = "Heat metrics: WBGT mean and T2M mean; greenness: GVI, NDVI, IMU total (p10 vs p90)"
)

curve_strong_paths <- build_compact_curves(
  heat_set = c("wbgt_min", "t2m_min"),
  out_stem = "paris_curves_effectmod_compact_strong",
  subtitle_txt = "Stronger-signal view: WBGT min and T2M min; greenness: GVI, NDVI, IMU total (p10 vs p90)"
)

# one-per-family best-signal panel for map-choice (for paper):
# WBGT min, T2M min, LST min x (GVI, NDVI, IMU total)
curve_bestsignal_paths <- build_compact_curves(
  heat_set = c("wbgt_min", "t2m_min", "lst_min"),
  out_stem = "paris_curves_effectmod_bestsignal_trio",
  subtitle_txt = "Best-signal family trio: WBGT min, T2M min, LST min; greenness: GVI, NDVI, IMU total (p10 vs p90)"
)

# 2) direct contrast figure
heat_order <- c("t2m_min", "t2m_mean", "t2m_max", "lst_min", "lst_mean", "lst_max", "wbgt_min", "wbgt_mean", "wbgt_max")
base_green <- "gvi_popw_points"
comparators <- c("ndvi_popw_ghs", "imu_veg_total_popw_ghs")

eff <- copy(harm)[
  heat_metric %in% heat_order &
    green_metric_raw %in% c(base_green, comparators),
  .(
    heat_metric, green_metric_raw,
    est = cr_mod_p99_p10p90,
    lo = cr_mod_p99_p10p90_low,
    hi = cr_mod_p99_p10p90_high
  )
]

eff[, `:=`(
  theta = log1p(est / 100),
  theta_lo = log1p(lo / 100),
  theta_hi = log1p(hi / 100)
)]
eff[, se_theta := (theta_hi - theta_lo) / (2 * 1.96)]

base_dt <- eff[green_metric_raw == base_green, .(heat_metric, theta_base = theta, se_base = se_theta)]

contrast_dt <- rbindlist(lapply(comparators, function(cmp) {
  cmp_dt <- eff[green_metric_raw == cmp, .(heat_metric, theta_cmp = theta, se_cmp = se_theta)]
  z <- merge(base_dt, cmp_dt, by = "heat_metric", all = FALSE)
  z[, `:=`(
    delta_theta = theta_base - theta_cmp,
    se_delta = sqrt(se_base^2 + se_cmp^2) # independence approximation
  )]
  z[, `:=`(
    delta_lo = delta_theta - 1.96 * se_delta,
    delta_hi = delta_theta + 1.96 * se_delta
  )]
  z[, `:=`(
    delta_pct = (exp(delta_theta) - 1) * 100,
    delta_pct_lo = (exp(delta_lo) - 1) * 100,
    delta_pct_hi = (exp(delta_hi) - 1) * 100,
    comparator = cmp
  )]
  z[]
}), fill = TRUE)

contrast_dt[, contrast_label := fifelse(
  comparator == "ndvi_popw_ghs",
  "GVI minus NDVI",
  "GVI minus IMU total"
)]
contrast_dt[, heat_label := factor(heat_labels[heat_metric], levels = rev(heat_labels[heat_order]))]

contrast_plot <- ggplot(
  contrast_dt,
  aes(x = delta_pct, y = heat_label)
) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35) +
  geom_errorbarh(aes(xmin = delta_pct_lo, xmax = delta_pct_hi), height = 0.22, linewidth = 0.45, colour = "#6b6b6b") +
  geom_point(size = 2.0, colour = "#1f4e79") +
  facet_wrap(~ contrast_label, ncol = 2) +
  labs(
    title = "Direct contrast of map-choice effects at p99",
    subtitle = "Primary estimand (p10→p90). Positive values mean stronger attenuation for GVI.",
    x = "Difference in attenuation (%): GVI minus comparator",
    y = NULL,
    caption = "CIs use an independent-SE approximation from reported model contrasts."
  ) +
  theme_bw(base_size = 10.8) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = "grey75")
  )

out_contrast_png <- file.path(fig_dir, tag_name("paris_effectmod_direct_contrasts.png", run_tag))
out_contrast_pdf <- file.path(fig_dir, tag_name("paris_effectmod_direct_contrasts.pdf", run_tag))
save_plot_dual(contrast_plot, out_contrast_png, out_contrast_pdf, width = 11.2, height = 5.8)

out_contrast_csv <- tag_name("paris_effectmod_direct_contrasts.csv", run_tag)
fwrite(contrast_dt, out_contrast_csv)

# 3) two-panel burden figure
greens_keep <- c(
  "GVI (point pop-w., GHS)",
  "NDVI (pop-w., GHS)",
  "IMU total (pop-w., GHS)"
)
burden <- copy(burden)[green %in% greens_keep]
burden[, panel := paste0(heat, " (", threshold, ")")]
panel_order <- c("WBGT mean (MMT)", "T2M mean (MMT)", "LST max (MMT)")
burden <- burden[panel %in% panel_order]

# add ref observed burden (average across metrics within panel) as per-100k/JJAS
obs_ref <- burden[, .(ad_obs_ref = mean(ad_obs, na.rm = TRUE)), by = panel]
obs_ref[, ad_obs_ref_per100k := ad_obs_ref / 2200000 * 100000 / 10]
burden <- merge(burden, obs_ref[, .(panel, ad_obs_ref_per100k)], by = "panel", all.x = TRUE)
burden[, panel_lab := sprintf("%s | obs %.2f/100k", panel, ad_obs_ref_per100k)]
burden[, panel_lab := factor(panel_lab, levels = sprintf("%s | obs %.2f/100k", panel_order, obs_ref[match(panel_order, panel), ad_obs_ref_per100k]))]
burden[, green := factor(green, levels = greens_keep)]

make_tile <- function(dt, val_col, title_txt, fill_name) {
  ggplot(dt, aes(x = green, y = panel_lab, fill = get(val_col))) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", get(val_col))), size = 3.0) +
    scale_fill_gradient2(
      low = "#b2182b", mid = "#f7f7f7", high = "#2166ac", midpoint = 0,
      name = fill_name
    ) +
    labs(x = NULL, y = NULL, title = title_txt) +
    theme_bw(base_size = 10.2) +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1),
      panel.grid = element_blank(),
      plot.title = element_text(size = 10.5)
    )
}

burden_per100k <- make_tile(
  burden,
  "ad_avoided_per100k_per_summer",
  "Avoided deaths per 100k (max-all)",
  "Avoided deaths\nper 100k"
)
burden_pct <- make_tile(
  burden,
  "pct_reduction_vs_observed_attr",
  "Percent reduction vs observed attributable burden (max-all)",
  "% reduction"
)

burden_two_panel <- burden_per100k + burden_pct +
  plot_layout(ncol = 1, guides = "keep") +
  plot_annotation(
    title = "Burden comparison across greenness metrics",
    subtitle = "Direct GHS-weighted heat and harmonized GHS-weighted greenness; interpret rows (endpoints) separately.",
    caption = "Scenario outputs (not causal claims). 'obs' in row labels is the reference attributable burden per 100k/JJAS."
  )

out_burden2_png <- file.path(fig_dir, tag_name("paris_burden_key_twopanel.png", run_tag))
out_burden2_pdf <- file.path(fig_dir, tag_name("paris_burden_key_twopanel.pdf", run_tag))
save_plot_dual(burden_two_panel, out_burden2_png, out_burden2_pdf, width = 8.5, height = 10)

cat("Saved focused figures:\n")
cat(" -", curve_main_paths[1], "\n")
cat(" -", curve_main_paths[2], "\n")
cat(" -", curve_strong_paths[1], "\n")
cat(" -", curve_strong_paths[2], "\n")
cat(" -", curve_bestsignal_paths[1], "\n")
cat(" -", curve_bestsignal_paths[2], "\n")
cat(" -", out_contrast_png, "\n")
cat(" -", out_contrast_pdf, "\n")
cat(" -", out_contrast_csv, "\n")
cat(" -", out_burden2_png, "\n")
cat(" -", out_burden2_pdf, "\n")
