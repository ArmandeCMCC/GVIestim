# plot_paris_weighted_vs_native_comparison.R
#
# Paired forest plot comparing the full weighted pipeline
# (GHS-weighted heat + GHS-weighted greenness) with the full native pipeline
# (area-average heat + unweighted greenness).
#
# For each heat metric, two dots are shown (weighted vs native),
# faceted by greenness family (GVI, NDVI, IMU total).

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

base_dir <- "/Users/armandeaboudrar-meda/Desktop/GVIestim/paper"
setwd(base_dir)

# load weighted (main) results
w <- fread("paris_effectmod_primary_table_p99_rerun_rule1_sum_ghsheat_ghsgreen.csv")

# load native/unweighted results
n <- fread("paris_effectmod_primary_table_p99_native_nativeheat.csv")

# green metric mapping: weighted -> native
green_map <- data.table(
  green_w = c("gvi_popw_points", "ndvi_popw_ghs", "imu_veg_total_popw_ghs"),
  green_n = c("gvi_mean", "ndvi_areaw_arr", "imu_veg_total"),
  green_family = c("GVI", "NDVI", "IMU total")
)

# subset weighted
w3 <- w[green_metric_raw %in% green_map$green_w,
         .(heat_metric, green_metric_raw,
           cr = cr_mod_primary, cr_lo = cr_mod_primary_low, cr_hi = cr_mod_primary_high)]
w3 <- merge(w3, green_map[, .(green_w, green_family)],
            by.x = "green_metric_raw", by.y = "green_w")
w3[, pipeline := "GHS-weighted heat + green"]

# subset native/unweighted
n3 <- n[green_metric_raw %in% green_map$green_n,
         .(heat_metric, green_metric_raw,
           cr = cr_mod_primary, cr_lo = cr_mod_primary_low, cr_hi = cr_mod_primary_high)]
n3 <- merge(n3, green_map[, .(green_n, green_family)],
            by.x = "green_metric_raw", by.y = "green_n")
n3[, pipeline := "Native heat + green"]

# combine
dt <- rbind(
  w3[, .(heat_metric, green_family, cr, cr_lo, cr_hi, pipeline)],
  n3[, .(heat_metric, green_family, cr, cr_lo, cr_hi, pipeline)]
)

# heat metric labels and ordering
heat_labels <- c(
  "t2m_min" = "T2M min", "t2m_mean" = "T2M mean", "t2m_max" = "T2M max",
  "lst_min" = "LST min", "lst_mean" = "LST mean", "lst_max" = "LST max",
  "wbgt_min" = "WBGT min", "wbgt_mean" = "WBGT mean", "wbgt_max" = "WBGT max"
)
heat_order <- c("T2M min", "T2M mean", "T2M max",
                "LST min", "LST mean", "LST max",
                "WBGT min", "WBGT mean", "WBGT max")

dt[, heat_label := heat_labels[heat_metric]]
dt[, heat_label := factor(heat_label, levels = rev(heat_order))]
dt[, green_family := factor(green_family, levels = c("GVI", "NDVI", "IMU total"))]
dt[, pipeline := factor(pipeline,
                        levels = c("GHS-weighted heat + green", "Native heat + green"))]

p <- ggplot(dt, aes(x = cr, y = heat_label, colour = pipeline, shape = pipeline)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.35) +
  geom_pointrange(
    aes(xmin = cr_lo, xmax = cr_hi),
    position = position_dodge(width = 0.55),
    size = 0.35, linewidth = 0.45
  ) +
  facet_wrap(~ green_family, ncol = 3) +
  scale_colour_manual(
    values = c("GHS-weighted heat + green" = "#1f4e79",
               "Native heat + green" = "#c44e52"),
    name = NULL
  ) +
  scale_shape_manual(
    values = c("GHS-weighted heat + green" = 16,
               "Native heat + green" = 17),
    name = NULL
  ) +
  labs(
    x = "Effect-modification estimate (% change in RR at p99, p10 \u2192 p90)",
    y = NULL,
    title = "Population-weighted vs native pipeline: effect-modification comparison",
    subtitle = "Primary estimand (p10\u2192p90 at p99). Each pair shares the same heat metric but differs in heat support and greenness weighting."
  ) +
  theme_bw(base_size = 10.5) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = "grey75")
  )

fig_dir <- file.path(base_dir, "figs")
out_png <- file.path(fig_dir, "paris_weighted_vs_native_forest.png")
out_pdf <- file.path(fig_dir, "paris_weighted_vs_native_forest.pdf")

ggsave(out_png, p, width = 12, height = 6.5, dpi = 300)
ggsave(out_pdf, p, width = 12, height = 6.5)

cat("Saved:\n -", out_png, "\n -", out_pdf, "\n")
