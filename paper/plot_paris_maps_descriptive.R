# plot_paris_maps_descriptive.R
#
# Four-panel choropleth map of Paris arrondissements:
#   (a) GVI (GHS-POP weighted, point-level)
#   (b) NDVI (GHS-POP weighted)
#   (c) IMU total vegetation (GHS-POP weighted)
#   (d) Mean daily maximum WBGT (GHS-POP sum-weighted, JJAS 2008-2017)

library(sf)
library(data.table)
library(ggplot2)
library(patchwork)

base_dir <- "/Users/armandeaboudrar-meda/Desktop/GVIestim/paper" # to modify if ran somewhere else 

# 1) Load arrondissement boundaries
arr_shp_path <- file.path(base_dir,
  "NDVI_Paris_2014_2017/_tmp_paris_shp/arrondissements.shp")
arr_sf <- st_read(arr_shp_path, quiet = TRUE)
arr_sf$arr <- as.character(arr_sf$c_arinsee)
paris_arr <- sprintf("751%02d", 1:20)
arr_sf <- arr_sf[arr_sf$arr %in% paris_arr, c("arr", "geometry")]
stopifnot(nrow(arr_sf) == 20)

# 2) Load greenness metrics
green <- fread(file.path(base_dir, "paris_arr_greenness_metrics.csv"))
green[, arr := as.character(arr)]

# 3) Compute mean WBGT max by arrondissement
heat <- fread(file.path(base_dir,
  "paris_panel_residence_urbclim_daily_2008_2017_ghspopweighted_sumproj_JJAS.csv"))
heat[, arr := as.character(arr)]
wbgt_mean <- heat[, .(wbgt_max_mean = mean(wbgt_max, na.rm = TRUE)), by = arr]

# 4) Merge into spatial object
arr_sf <- merge(arr_sf, green[, .(arr, gvi_popw_points, ndvi_popw_ghs,
                                   imu_veg_total_popw_ghs)],
                by = "arr", all.x = TRUE)
arr_sf <- merge(arr_sf, wbgt_mean, by = "arr", all.x = TRUE)

# Arrondissement labels (1-20)
arr_sf$arr_num <- as.integer(substr(arr_sf$arr, 4, 5))
centroids <- st_centroid(st_geometry(arr_sf))
arr_sf$cx <- st_coordinates(centroids)[, 1]
arr_sf$cy <- st_coordinates(centroids)[, 2]

# 5) Common theme
map_theme <- theme_minimal(base_size = 10) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.key.width = unit(1.2, "cm"),
    legend.key.height = unit(0.3, "cm"),
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    plot.margin = margin(2, 2, 2, 2)
  )

make_map <- function(sf_data, fill_var, title, fill_label, palette = "YlGn",
                     direction = 1) {
  ggplot(sf_data) +
    geom_sf(aes(fill = .data[[fill_var]]), colour = "grey30", linewidth = 0.3) +
    geom_text(aes(x = cx, y = cy, label = arr_num),
              size = 2.5, colour = "grey25", fontface = "bold") +
    scale_fill_distiller(palette = palette, direction = direction,
                         name = fill_label) +
    ggtitle(title) +
    map_theme
}

# 6) Build panels
p_gvi <- make_map(arr_sf, "gvi_popw_points",
                  "(a) Green View Index",
                  "GVI (%)")

p_ndvi <- make_map(arr_sf, "ndvi_popw_ghs",
                   "(b) NDVI",
                   "NDVI")

p_imu <- make_map(arr_sf, "imu_veg_total_popw_ghs",
                  "(c) Planimetric vegetation",
                  "IMU total (%)")

p_wbgt <- make_map(arr_sf, "wbgt_max_mean",
                   "(d) Mean daily max WBGT",
                   "WBGT max (\u00B0C)",
                   palette = "YlOrRd", direction = 1)

# 7) Combine and save
p_all <- (p_gvi | p_ndvi) / (p_imu | p_wbgt) +
  plot_annotation(
    theme = theme(plot.margin = margin(5, 5, 5, 5))
  )

fig_dir <- file.path(base_dir, "figs")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

out_pdf <- file.path(fig_dir, "paris_maps_descriptive.pdf")
out_png <- file.path(fig_dir, "paris_maps_descriptive.png")

ggsave(out_pdf, p_all, width = 9, height = 9, device = cairo_pdf)
ggsave(out_png, p_all, width = 9, height = 9, dpi = 300)

cat("Saved:\n -", out_pdf, "\n -", out_png, "\n")
