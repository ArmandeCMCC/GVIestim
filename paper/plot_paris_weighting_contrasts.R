# plot_paris_weighting_contrasts.R
# Generates the "weighted vs unweighted" contrast figure:
#   For each heat x green combination, shows the difference in
#   effect-modification estimate (p10->p90 at p99) between
#   the population-weighted and the unweighted (native) metric.
#   Positive = stronger attenuation for the population-weighted metric.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

base_dir <- "/Users/armandeaboudrar-meda/Desktop/GVIestim/paper" # relative 
setwd(base_dir)

run_tag <- "_rerun_rule1_sum_ghsheat_ghsgreen"

weighted_file <- paste0("paris_effectmod_primary_table_p99", run_tag, ".csv")
native_file   <- paste0("paris_effectmod_primary_table_p99_native", run_tag, ".csv")

w <- fread(weighted_file)
n <- fread(native_file)

# Map weighted green metrics to their native counterparts
green_map <- data.table(
  green_w = c("gvi_popw_points", "ndvi_popw_ghs", "imu_veg_total_popw_ghs"),
  green_n = c("gvi_mean", "ndvi_native_jjas", "imu_veg_total"),
  green_label = c("GVI", "NDVI", "IMU total")
)

# Subset to primary 3 green metrics
w3 <- w[green_metric_raw %in% green_map$green_w,
         .(heat_metric, green_metric_raw,
           cr_w = cr_mod_primary, cr_w_low = cr_mod_primary_low, cr_w_high = cr_mod_primary_high)]
n3 <- n[green_metric_raw %in% green_map$green_n,
         .(heat_metric, green_metric_raw,
           cr_n = cr_mod_primary, cr_n_low = cr_mod_primary_low, cr_n_high = cr_mod_primary_high)]

# Add mapping keys
w3 <- merge(w3, green_map, by.x = "green_metric_raw", by.y = "green_w")
n3 <- merge(n3, green_map[, .(green_n, green_label)], by.x = "green_metric_raw", by.y = "green_n")

# Merge weighted and native
dt <- merge(
  w3[, .(heat_metric, green_label, cr_w, cr_w_low, cr_w_high)],
  n3[, .(heat_metric, green_label, cr_n, cr_n_low, cr_n_high)],
  by = c("heat_metric", "green_label")
)

# Contrast: native - weighted (positive = weighted shows stronger attenuation)
dt[, contrast := cr_n - cr_w]
# Approximate SE: sqrt(se_w^2 + se_n^2), from CI widths
dt[, se_w := (cr_w_high - cr_w_low) / (2 * 1.96)]
dt[, se_n := (cr_n_high - cr_n_low) / (2 * 1.96)]
dt[, se_contrast := sqrt(se_w^2 + se_n^2)]
dt[, contrast_low := contrast - 1.96 * se_contrast]
dt[, contrast_high := contrast + 1.96 * se_contrast]

# Heat metric labels
heat_labels <- c(
  "t2m_min" = "T2M min", "t2m_mean" = "T2M mean", "t2m_max" = "T2M max",
  "lst_min" = "LST min", "lst_mean" = "LST mean", "lst_max" = "LST max",
  "wbgt_min" = "WBGT min", "wbgt_mean" = "WBGT mean", "wbgt_max" = "WBGT max"
)
dt[, heat_label := heat_labels[heat_metric]]
dt[, heat_family := sub("_.*", "", heat_metric)]
dt[, heat_family := toupper(heat_family)]
dt[, heat_family := factor(heat_family, levels = c("T2M", "LST", "WBGT"))]

# Order heat metrics: min, mean, max within each family
heat_order <- c("T2M min", "T2M mean", "T2M max",
                "LST min", "LST mean", "LST max",
                "WBGT min", "WBGT mean", "WBGT max")
dt[, heat_label := factor(heat_label, levels = rev(heat_order))]

dt[, green_label := factor(green_label, levels = c("GVI", "NDVI", "IMU total"))]

p <- ggplot(dt, aes(x = contrast, y = heat_label, colour = green_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(
    aes(xmin = contrast_low, xmax = contrast_high),
    position = position_dodge(width = 0.6),
    size = 0.4
  ) +
  facet_wrap(~green_label, ncol = 3) +
  scale_colour_manual(values = c("GVI" = "#2ca02c", "NDVI" = "#1f77b4", "IMU total" = "#d62728")) +
  labs(
    x = "Difference in attenuation (pp): population-weighted minus unweighted",
    y = NULL,
    title = "Effect of population-weighting on effect-modification estimates",
    subtitle = "Primary estimand (p10\u2192p90 at p99). Positive = stronger attenuation for pop-weighted metric.",
    caption = "CIs use independent-SE approximation. Corrected run (rule=1, sum projection)."
  ) +
  theme_bw(base_size = 10.8) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = "grey75")
  )

fig_dir <- file.path(base_dir, "figs")

out_png <- file.path(fig_dir, paste0("paris_weighting_contrasts", run_tag, ".png"))
out_pdf <- file.path(fig_dir, paste0("paris_weighting_contrasts", run_tag, ".pdf"))

ggsave(out_png, p, width = 12, height = 5.5, dpi = 300)
ggsave(out_pdf, p, width = 12, height = 5.5)

cat("Saved:\n")
cat(" -", out_png, "\n")
cat(" -", out_pdf, "\n")

# Also print summary
cat("\nContrast summary (positive = pop-weighted shows stronger attenuation):\n")
print(dt[, .(heat_label, green_label, contrast = round(contrast, 2),
             CI = paste0("[", round(contrast_low, 2), ", ", round(contrast_high, 2), "]"))])
