# plot_paris_effectmod_results.R
#
# Produce main figures (forest plots + curve plots).
# Main paper run tag: _rerun_rule1_sum_ghsheat_ghsgreen.
#
# plotting / results-packaging script for Paris heat × greenness CTS/DLNM
# results, updated for recent ggplot2 versions.

# Inputs:
#   - paris_effectmod_results_harmonized.csv
#   - paris_effectmod_heat_green_curves.csv

# Main outputs:
#   1) Forest plot of PRIMARY effect modification at p99 (p10 -> p90 contrast)
#   2) Companion forest plots for harmonized (-1 IQR -> +1 IQR) and
#      Paris-paper (min -> max) contrasts
#   3) Exposure-response curve plot for selected heat × greenness combinations

# Output files:
#   - figs/paris_forest_effectmod_p99.png
#   - figs/paris_forest_effectmod_p99.pdf
#   - figs/paris_forest_effectmod_p99_native.png
#   - figs/paris_forest_effectmod_p99_native.pdf
#   - figs/paris_forest_effectmod_p99_by_family.png
#   - figs/paris_forest_effectmod_p99_by_family.pdf
#   - figs/paris_forest_effectmod_p99_iqr.png
#   - figs/paris_forest_effectmod_p99_iqr.pdf
#   - figs/paris_forest_effectmod_p99_minmax.png
#   - figs/paris_forest_effectmod_p99_minmax.pdf
#   - figs/paris_curves_effectmod_selected.png           (PRIMARY p10->p90)
#   - figs/paris_curves_effectmod_selected.pdf
#   - figs/paris_curves_effectmod_selected_iqr.png       (sensitivity)
#   - figs/paris_curves_effectmod_selected_iqr.pdf
#   - figs/paris_curves_effectmod_selected_native.png
#   - figs/paris_curves_effectmod_selected_native.pdf
#   - figs/paris_curves_effectmod_selected_iqr_native.png
#   - figs/paris_curves_effectmod_selected_iqr_native.pdf

library(data.table)
library(ggplot2)

base_dir <- "/Users/armandeaboudrar-meda/Desktop/GVIestim/paper" # relative : base_dir <- here::here() 

norm_tag <- function(x) {
  x <- trimws(x)
  if (!nzchar(x)) return("")
  if (!startsWith(x, "_")) x <- paste0("_", x)
  x
}

tag_name <- function(filename, tag) {
  if (!nzchar(tag)) return(filename)
  sub("\\.([A-Za-z0-9]+)$", paste0(tag, ".\\1"), filename)
}

tag_file <- function(stem, tag, ext = ".csv") {
  file.path(base_dir, paste0(stem, tag, ext))
}

run_tag <- norm_tag(Sys.getenv("RUN_TAG", "rerun_rule1_sum_ghsheat_ghsgreen"))
harm_tag <- norm_tag(Sys.getenv("HARM_TAG", run_tag))
curves_tag <- norm_tag(Sys.getenv("CURVES_TAG", run_tag))

in_harm   <- tag_file("paris_effectmod_results_harmonized", harm_tag)
in_curves <- tag_file("paris_effectmod_heat_green_curves", curves_tag)

fig_dir <- file.path(base_dir, "figs")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

out_forest_png    <- file.path(fig_dir, tag_name("paris_forest_effectmod_p99.png", run_tag))
out_forest_pdf    <- file.path(fig_dir, tag_name("paris_forest_effectmod_p99.pdf", run_tag))
out_forest_native_png <- file.path(fig_dir, tag_name("paris_forest_effectmod_p99_native.png", run_tag))
out_forest_native_pdf <- file.path(fig_dir, tag_name("paris_forest_effectmod_p99_native.pdf", run_tag))
out_forestfam_png <- file.path(fig_dir, tag_name("paris_forest_effectmod_p99_by_family.png", run_tag))
out_forestfam_pdf <- file.path(fig_dir, tag_name("paris_forest_effectmod_p99_by_family.pdf", run_tag))
out_forestiqr_png <- file.path(fig_dir, tag_name("paris_forest_effectmod_p99_iqr.png", run_tag))
out_forestiqr_pdf <- file.path(fig_dir, tag_name("paris_forest_effectmod_p99_iqr.pdf", run_tag))
out_forestpaper_png <- file.path(fig_dir, tag_name("paris_forest_effectmod_p99_minmax.png", run_tag))
out_forestpaper_pdf <- file.path(fig_dir, tag_name("paris_forest_effectmod_p99_minmax.pdf", run_tag))
out_forestpaper_native_png <- file.path(fig_dir, tag_name("paris_forest_effectmod_p99_minmax_native.png", run_tag))
out_forestpaper_native_pdf <- file.path(fig_dir, tag_name("paris_forest_effectmod_p99_minmax_native.pdf", run_tag))
out_curves_png    <- file.path(fig_dir, tag_name("paris_curves_effectmod_selected.png", run_tag))
out_curves_pdf    <- file.path(fig_dir, tag_name("paris_curves_effectmod_selected.pdf", run_tag))
out_curvesiqr_png <- file.path(fig_dir, tag_name("paris_curves_effectmod_selected_iqr.png", run_tag))
out_curvesiqr_pdf <- file.path(fig_dir, tag_name("paris_curves_effectmod_selected_iqr.pdf", run_tag))
out_curves_native_png <- file.path(fig_dir, tag_name("paris_curves_effectmod_selected_native.png", run_tag))
out_curves_native_pdf <- file.path(fig_dir, tag_name("paris_curves_effectmod_selected_native.pdf", run_tag))
out_curvesiqr_native_png <- file.path(fig_dir, tag_name("paris_curves_effectmod_selected_iqr_native.png", run_tag))
out_curvesiqr_native_pdf <- file.path(fig_dir, tag_name("paris_curves_effectmod_selected_iqr_native.pdf", run_tag))

cat("\nPlotting config:\n")
cat(" - HARM_TAG:", ifelse(nzchar(harm_tag), harm_tag, "<none>"), "\n")
cat(" - CURVES_TAG:", ifelse(nzchar(curves_tag), curves_tag, "<none>"), "\n")
cat(" - RUN_TAG:", ifelse(nzchar(run_tag), run_tag, "<none>"), "\n")
cat(" - input harmonized:", in_harm, "\n")
cat(" - input curves:", in_curves, "\n")

# 1) load data
harm <- fread(in_harm)
curves <- fread(in_curves)

stopifnot(all(c(
  "heat_metric", "green_metric_raw",
  "cr_mod_p99", "cr_mod_p99_low", "cr_mod_p99_high",
  "cr_mod_p99_p10p90", "cr_mod_p99_p10p90_low", "cr_mod_p99_p10p90_high",
  "cr_mod_p99_minmax", "cr_mod_p99_minmax_low", "cr_mod_p99_minmax_high"
) %in% names(harm)))

stopifnot(all(c(
  "heat_metric", "green_metric_raw", "green_level",
  "temp", "rr", "rr_low", "rr_high"
) %in% names(curves)))

# 2) save helper 
save_plot_dual <- function(plot_obj, png_path, pdf_path, width, height, dpi = 320) {
  # remove stale files if they exist
  if (file.exists(png_path)) file.remove(png_path)
  if (file.exists(pdf_path)) file.remove(pdf_path)
  
  # PNG: prefer ragg
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(
      filename = png_path,
      width = width,
      height = height,
      units = "in",
      res = dpi
    )
    print(plot_obj)
    dev.off()
  } else if (capabilities("aqua")) {
    png(
      filename = png_path,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      type = "quartz"
    )
    print(plot_obj)
    dev.off()
  } else {
    stop("Need either the 'ragg' package or macOS Quartz support to save PNG files.")
  }
  
  # PDF: base pdf() avoids cairo
  pdf(
    file = pdf_path,
    width = width,
    height = height,
    onefile = FALSE,
    useDingbats = FALSE
  )
  print(plot_obj)
  dev.off()
}

# 3) labels / ordering helpers
heat_order <- c(
  "t2m_min", "t2m_mean", "t2m_max",
  "lst_min", "lst_mean", "lst_max",
  "wbgt_min", "wbgt_mean", "wbgt_max"
)

green_order <- c(
  "gvi_popw_points",
  "ndvi_popw_ghs",
  "imu_veg_total_popw_ghs",
  "imu_high_popw_ghs",
  "imu_low_popw_ghs"
)

green_order_native <- c(
  "gvi_mean",
  "ndvi_areaw_arr",
  "imu_veg_total"
)

heat_labels <- c(
  t2m_min   = "T2M min",
  t2m_mean  = "T2M mean",
  t2m_max   = "T2M max",
  lst_min   = "LST min",
  lst_mean  = "LST mean",
  lst_max   = "LST max",
  wbgt_min  = "WBGT min",
  wbgt_mean = "WBGT mean",
  wbgt_max  = "WBGT max"
)

green_labels <- c(
  gvi_mean           = "GVI (native mean)",
  gvi_popw_points    = "GVI (point pop-w.)",
  ndvi_native_jjas   = "NDVI (native JJAS)",
  ndvi_popw_ghs      = "NDVI (GHS pop-w.)",
  ndvi_areaw_arr     = "NDVI (arr area-w.)",
  imu_veg_total      = "IMU total (native)",
  imu_veg_total_popw = "IMU total (pop-w.)",
  imu_high_popw      = "IMU high veg (pop-w.)",
  imu_low_popw       = "IMU low veg (pop-w.)",
  imu_veg_total_popw_ghs = "IMU total (GHS pop-w.)",
  imu_high_popw_ghs      = "IMU high veg (GHS pop-w.)",
  imu_low_popw_ghs       = "IMU low veg (GHS pop-w.)"
)

heat_family_labels <- c(
  T2M = "T2M",
  LST = "LST",
  WBGT = "WBGT"
)

level_labels <- c(
  low_minus_1iqr = "Low greenness (-1 IQR)",
  high_plus_1iqr = "High greenness (+1 IQR)",
  low_p10 = "Low greenness (p10)",
  high_p90 = "High greenness (p90)"
)

get_heat_family <- function(x) {
  fifelse(grepl("^t2m_", x), "T2M",
          fifelse(grepl("^lst_", x), "LST",
                  fifelse(grepl("^wbgt_", x), "WBGT", "Other")))
}

get_heat_stat <- function(x) sub("^[^_]+_", "", x)

# 4) clean / relabel harmonised results
plot_dt <- copy(harm)[
  heat_metric %in% heat_order & green_metric_raw %in% green_order
]

plot_dt_native <- copy(harm)[
  heat_metric %in% heat_order & green_metric_raw %in% green_order_native
]

plot_dt[, heat_metric := factor(heat_metric, levels = rev(heat_order))]
plot_dt[, green_metric_raw := factor(green_metric_raw, levels = green_order)]

plot_dt_native[, heat_metric := factor(heat_metric, levels = rev(heat_order))]
plot_dt_native[, green_metric_raw := factor(green_metric_raw, levels = green_order_native)]

plot_dt[, heat_label := factor(
  heat_labels[as.character(heat_metric)],
  levels = rev(heat_labels[heat_order])
)]

plot_dt_native[, heat_label := factor(
  heat_labels[as.character(heat_metric)],
  levels = rev(heat_labels[heat_order])
)]

plot_dt[, green_label := factor(
  green_labels[as.character(green_metric_raw)],
  levels = green_labels[green_order]
)]

plot_dt_native[, green_label := factor(
  green_labels[as.character(green_metric_raw)],
  levels = green_labels[green_order_native]
)]

plot_dt[, heat_family := factor(
  get_heat_family(as.character(heat_metric)),
  levels = c("T2M", "LST", "WBGT"),
  labels = heat_family_labels[c("T2M", "LST", "WBGT")]
)]

plot_dt_native[, heat_family := factor(
  get_heat_family(as.character(heat_metric)),
  levels = c("T2M", "LST", "WBGT"),
  labels = heat_family_labels[c("T2M", "LST", "WBGT")]
)]

plot_dt[, heat_stat := factor(
  get_heat_stat(as.character(heat_metric)),
  levels = c("min", "mean", "max"),
  labels = c("min", "mean", "max")
)]

plot_dt_native[, heat_stat := factor(
  get_heat_stat(as.character(heat_metric)),
  levels = c("min", "mean", "max"),
  labels = c("min", "mean", "max")
)]

plot_dt[, heat_y := factor(
  paste(heat_family, heat_stat),
  levels = rev(c(
    "T2M min", "T2M mean", "T2M max",
    "LST min", "LST mean", "LST max",
    "WBGT min", "WBGT mean", "WBGT max"
  ))
)]

plot_dt_native[, heat_y := factor(
  paste(heat_family, heat_stat),
  levels = rev(c(
    "T2M min", "T2M mean", "T2M max",
    "LST min", "LST mean", "LST max",
    "WBGT min", "WBGT mean", "WBGT max"
  ))
)]

# Estimands at p99
plot_dt[, `:=`(
  cr_mod_primary = cr_mod_p99_p10p90,
  cr_mod_primary_low = cr_mod_p99_p10p90_low,
  cr_mod_primary_high = cr_mod_p99_p10p90_high,
  cr_mod_iqr = cr_mod_p99,
  cr_mod_iqr_low = cr_mod_p99_low,
  cr_mod_iqr_high = cr_mod_p99_high,
  cr_mod_paper = cr_mod_p99_minmax,
  cr_mod_paper_low = cr_mod_p99_minmax_low,
  cr_mod_paper_high = cr_mod_p99_minmax_high
)]

plot_dt_native[, `:=`(
  cr_mod_primary = cr_mod_p99_p10p90,
  cr_mod_primary_low = cr_mod_p99_p10p90_low,
  cr_mod_primary_high = cr_mod_p99_p10p90_high,
  cr_mod_iqr = cr_mod_p99,
  cr_mod_iqr_low = cr_mod_p99_low,
  cr_mod_iqr_high = cr_mod_p99_high,
  cr_mod_paper = cr_mod_p99_minmax,
  cr_mod_paper_low = cr_mod_p99_minmax_low,
  cr_mod_paper_high = cr_mod_p99_minmax_high
)]

# 5) forest plot: PRIMARY p10 -> p90 contrast at p99
forest_plot <- ggplot(
  plot_dt,
  aes(x = cr_mod_primary, y = heat_y, xmin = cr_mod_primary_low, xmax = cr_mod_primary_high)
) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(orientation = "y", width = 0.18, linewidth = 0.45) +
  geom_point(size = 1.9) +
  facet_wrap(~ green_label, ncol = 2, scales = "fixed") +
  labs(
    title = "Effect modification of heat-related mortality by greenness",
    subtitle = "PRIMARY contrast: greenness from p10 to p90, evaluated at each heat metric's p99",
    x = "% change in RR at p99 for greener vs less green districts (p10 -> p90)",
    y = NULL,
    caption = "Negative values indicate lower heat-mortality risk in greener districts."
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title.position = "plot"
  )

save_plot_dual(forest_plot, out_forest_png, out_forest_pdf, width = 12, height = 8.5)

# native-main forest plot
forest_plot_native <- ggplot(
  plot_dt_native,
  aes(x = cr_mod_primary, y = heat_y, xmin = cr_mod_primary_low, xmax = cr_mod_primary_high)
) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(orientation = "y", width = 0.18, linewidth = 0.45) +
  geom_point(size = 1.9) +
  facet_wrap(~ green_label, ncol = 2, scales = "fixed") +
  labs(
    title = "Effect modification of heat-related mortality by greenness",
    subtitle = "PRIMARY native-metric set: greenness from p10 to p90 at each heat metric's p99",
    x = "% change in RR at p99 for greener vs less green districts (p10 -> p90)",
    y = NULL,
    caption = "Negative values indicate lower heat-mortality risk in greener districts."
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title.position = "plot"
  )

save_plot_dual(forest_plot_native, out_forest_native_png, out_forest_native_pdf, width = 11, height = 8.2)

# 5b) forest plot: harmonised -1 IQR -> +1 IQR at p99
forest_iqr_plot <- ggplot(
  plot_dt,
  aes(x = cr_mod_iqr, y = heat_y, xmin = cr_mod_iqr_low, xmax = cr_mod_iqr_high)
) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(orientation = "y", width = 0.18, linewidth = 0.45) +
  geom_point(size = 1.9) +
  facet_wrap(~ green_label, ncol = 2, scales = "fixed") +
  labs(
    title = "Effect modification of heat-related mortality by greenness",
    subtitle = "Harmonized contrast: greenness from -1 IQR to +1 IQR, evaluated at each heat metric's p99",
    x = "% change in RR at p99 (harmonized -1 IQR -> +1 IQR)",
    y = NULL,
    caption = "Negative values indicate lower heat-mortality risk in greener districts."
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title.position = "plot"
  )

save_plot_dual(forest_iqr_plot, out_forestiqr_png, out_forestiqr_pdf, width = 12, height = 8.5)

# 5c) forest plot: Paris-paper min -> max contrast at p99
forest_paper_plot <- ggplot(
  plot_dt,
  aes(x = cr_mod_paper, y = heat_y, xmin = cr_mod_paper_low, xmax = cr_mod_paper_high)
) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(orientation = "y", width = 0.18, linewidth = 0.45) +
  geom_point(size = 1.9) +
  facet_wrap(~ green_label, ncol = 2, scales = "fixed") +
  labs(
    title = "Effect modification of heat-related mortality by greenness",
    subtitle = "Paris-paper contrast: greenness from observed minimum to observed maximum at p99",
    x = "% change in RR at p99 (min -> max greenness)",
    y = NULL,
    caption = "Negative values indicate lower heat-mortality risk in greener districts."
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title.position = "plot"
  )

save_plot_dual(forest_paper_plot, out_forestpaper_png, out_forestpaper_pdf, width = 12, height = 8.5)

# native-main Paris-paper min -> max forest
forest_paper_plot_native <- ggplot(
  plot_dt_native,
  aes(x = cr_mod_paper, y = heat_y, xmin = cr_mod_paper_low, xmax = cr_mod_paper_high)
) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(orientation = "y", width = 0.18, linewidth = 0.45) +
  geom_point(size = 1.9) +
  facet_wrap(~ green_label, ncol = 2, scales = "fixed") +
  labs(
    title = "Effect modification of heat-related mortality by greenness",
    subtitle = "Paris-paper contrast (native-metric set): observed minimum to observed maximum at p99",
    x = "% change in RR at p99 (min -> max greenness)",
    y = NULL,
    caption = "Negative values indicate lower heat-mortality risk in greener districts."
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title.position = "plot"
  )

save_plot_dual(
  forest_paper_plot_native,
  out_forestpaper_native_png,
  out_forestpaper_native_pdf,
  width = 11,
  height = 8.2
)

# 6) forest plot by heat family
forest_family_plot <- ggplot(
  plot_dt,
  aes(x = cr_mod_primary, y = heat_stat, xmin = cr_mod_primary_low, xmax = cr_mod_primary_high)
) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(orientation = "y", width = 0.18, linewidth = 0.45) +
  geom_point(size = 1.9) +
  facet_grid(green_label ~ heat_family, scales = "fixed") +
  labs(
    title = "Effect modification at p99 by heat family and greenness metric",
    subtitle = "PRIMARY contrast: greenness from p10 to p90",
    x = "% change in RR (p10 -> p90)",
    y = "Heat summary statistic",
    caption = "Negative values indicate lower heat-mortality risk in greener districts."
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    panel.grid.minor = element_blank(),
    plot.title.position = "plot"
  )

save_plot_dual(forest_family_plot, out_forestfam_png, out_forestfam_pdf, width = 13, height = 9)

# 7) curve plot subset
curve_heat_metrics <- c("t2m_mean", "lst_max", "wbgt_mean")
curve_green_metrics <- c(
  "gvi_popw_points",
  "ndvi_popw_ghs",
  "imu_veg_total_popw_ghs",
  "imu_high_popw_ghs",
  "imu_low_popw_ghs"
)
curve_green_metrics_native <- c(
  "gvi_mean",
  "ndvi_areaw_arr",
  "imu_veg_total"
)

curve_base <- copy(curves)[
  heat_metric %in% curve_heat_metrics &
    green_metric_raw %in% curve_green_metrics
]

curve_base_native <- copy(curves)[
  heat_metric %in% curve_heat_metrics &
    green_metric_raw %in% curve_green_metrics_native
]

build_curve_plot <- function(dt, level_order, title_txt, subtitle_txt, green_metrics_levels) {
  dt <- copy(dt)
  dt[, heat_metric := factor(heat_metric, levels = curve_heat_metrics)]
  dt[, green_metric_raw := factor(green_metric_raw, levels = green_metrics_levels)]
  dt[, green_level := factor(
    green_level,
    levels = level_order,
    labels = level_labels[level_order]
  )]
  
  dt[, heat_label := factor(
    heat_labels[as.character(heat_metric)],
    levels = heat_labels[curve_heat_metrics]
  )]
  
  dt[, green_label := factor(
    green_labels[as.character(green_metric_raw)],
    levels = green_labels[green_metrics_levels]
  )]
  
  guides <- unique(dt[, .(
    heat_metric, green_metric_raw, heat_label, green_label, heat_mmt, heat_p99
  )])
  
  ggplot(
    dt,
    aes(x = temp, y = rr, group = green_level, linetype = green_level)
  ) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.35) +
    geom_vline(
      data = guides,
      aes(xintercept = heat_mmt),
      linetype = "dotted",
      linewidth = 0.35,
      show.legend = FALSE
    ) +
    geom_vline(
      data = guides,
      aes(xintercept = heat_p99),
      linetype = "dotdash",
      linewidth = 0.35,
      show.legend = FALSE
    ) +
    geom_ribbon(
      aes(ymin = rr_low, ymax = rr_high, fill = green_level),
      alpha = 0.14,
      colour = NA
    ) +
    geom_line(linewidth = 0.7) +
    facet_grid(green_label ~ heat_label, scales = "free_x") +
    labs(
      title = title_txt,
      subtitle = subtitle_txt,
      x = "Heat metric value",
      y = "Relative risk (RR)",
      linetype = NULL,
      fill = NULL
    ) +
    theme_bw(base_size = 10.5) +
    theme(
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      panel.grid.minor = element_blank(),
      plot.title.position = "plot",
      legend.position = "bottom"
    )
}

# PRIMARY curve plot: p10/p90 greenness
curve_primary_levels <- c("low_p10", "high_p90")
if (!all(curve_primary_levels %in% unique(curve_base$green_level))) {
  stop("Missing p10/p90 curve levels in curves input. Re-run estimate_paris_effect_modification_cts_heat_green.R.")
}
curve_primary_dt <- curve_base[green_level %in% curve_primary_levels]
  curve_primary_plot <- build_curve_plot(
    curve_primary_dt,
    level_order = curve_primary_levels,
    title_txt = "Exposure-response curves by heat metric and greenness level",
    subtitle_txt = "PRIMARY contrast: low/high greenness as p10/p90; curves centered at metric-specific MMT",
    green_metrics_levels = curve_green_metrics
  )
save_plot_dual(curve_primary_plot, out_curves_png, out_curves_pdf, width = 13, height = 13)

# Sensitivity curve plot: -1/+1 IQR greenness
curve_iqr_levels <- c("low_minus_1iqr", "high_plus_1iqr")
if (all(curve_iqr_levels %in% unique(curve_base$green_level))) {
  curve_iqr_dt <- curve_base[green_level %in% curve_iqr_levels]
  curve_iqr_plot <- build_curve_plot(
    curve_iqr_dt,
    level_order = curve_iqr_levels,
    title_txt = "Exposure-response curves by heat metric and greenness level",
    subtitle_txt = "Sensitivity contrast: low/high greenness as -1/+1 IQR; curves centered at metric-specific MMT",
    green_metrics_levels = curve_green_metrics
  )
  save_plot_dual(curve_iqr_plot, out_curvesiqr_png, out_curvesiqr_pdf, width = 13, height = 13)
}

# main curve plots
curve_primary_dt_native <- curve_base_native[green_level %in% curve_primary_levels]
if (nrow(curve_primary_dt_native) > 0) {
  curve_primary_plot_native <- build_curve_plot(
    curve_primary_dt_native,
    level_order = curve_primary_levels,
    title_txt = "Exposure-response curves by heat metric and greenness level",
    subtitle_txt = "PRIMARY native-metric set: low/high greenness as p10/p90; curves centered at metric-specific MMT",
    green_metrics_levels = curve_green_metrics_native
  )
  save_plot_dual(curve_primary_plot_native, out_curves_native_png, out_curves_native_pdf, width = 12.5, height = 10.8)
}

if (all(curve_iqr_levels %in% unique(curve_base_native$green_level))) {
  curve_iqr_dt_native <- curve_base_native[green_level %in% curve_iqr_levels]
  curve_iqr_plot_native <- build_curve_plot(
    curve_iqr_dt_native,
    level_order = curve_iqr_levels,
    title_txt = "Exposure-response curves by heat metric and greenness level",
    subtitle_txt = "Sensitivity native-metric set: low/high greenness as -1/+1 IQR; curves centered at metric-specific MMT",
    green_metrics_levels = curve_green_metrics_native
  )
  save_plot_dual(curve_iqr_plot_native, out_curvesiqr_native_png, out_curvesiqr_native_pdf, width = 12.5, height = 10.8)
}

# 8) console preview
cat("\nSaved figures:\n")
cat(" -", out_forest_png, "\n")
cat(" -", out_forest_pdf, "\n")
cat(" -", out_forest_native_png, "\n")
cat(" -", out_forest_native_pdf, "\n")
cat(" -", out_forestfam_png, "\n")
cat(" -", out_forestfam_pdf, "\n")
cat(" -", out_forestiqr_png, "\n")
cat(" -", out_forestiqr_pdf, "\n")
cat(" -", out_forestpaper_png, "\n")
cat(" -", out_forestpaper_pdf, "\n")
cat(" -", out_forestpaper_native_png, "\n")
cat(" -", out_forestpaper_native_pdf, "\n")
cat(" -", out_curves_png, "\n")
cat(" -", out_curves_pdf, "\n")
if (file.exists(out_curvesiqr_png)) cat(" -", out_curvesiqr_png, "\n")
if (file.exists(out_curvesiqr_pdf)) cat(" -", out_curvesiqr_pdf, "\n")
if (file.exists(out_curves_native_png)) cat(" -", out_curves_native_png, "\n")
if (file.exists(out_curves_native_pdf)) cat(" -", out_curves_native_pdf, "\n")
if (file.exists(out_curvesiqr_native_png)) cat(" -", out_curvesiqr_native_png, "\n")
if (file.exists(out_curvesiqr_native_pdf)) cat(" -", out_curvesiqr_native_pdf, "\n")

cat("\nMost protective p99 results (PRIMARY p10 -> p90):\n")
print(
  plot_dt[order(cr_mod_p99_p10p90)][1:12, .(
    heat_metric = as.character(heat_metric),
    green_metric = as.character(green_metric_raw),
    cr_mod_primary,
    cr_mod_primary_low,
    cr_mod_primary_high
  )]
)

cat("\nMost protective p99 results (Paris-paper min -> max):\n")
print(
  plot_dt[order(cr_mod_paper)][1:12, .(
    heat_metric = as.character(heat_metric),
    green_metric = as.character(green_metric_raw),
    cr_mod_paper,
    cr_mod_paper_low,
    cr_mod_paper_high
  )]
)

cat("\nDone.\n")
