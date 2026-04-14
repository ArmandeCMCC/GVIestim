# build_paris_arr_ndvi_metrics_esa10m.R
#
# Replacement for build_paris_arr_ndvi_metrics.R
# Uses ESA WorldCover Sentinel-2 10m NDVI composites (2020-2021)
# instead of Copernicus CLMS 300m 10-daily NDVI (2014-2017).
#
# Input tiles (already downloaded):
#   ESA_WC_NDVI_2020_N48E002.tif  (3 bands: p90, p50, p10)
#   ESA_WC_NDVI_2021_N48E002.tif  (3 bands: p90, p50, p10)
#
# Band 1 (p90) represents annual peak greenness — analogous to
# the summer JJAS mean used in the Copernicus 300m pipeline.
# We average p90 across 2020 and 2021 for temporal smoothing.
#
# Produces both GHS-POP weighted (main paper) and unweighted/native:
#   - ndvi_popw_ghs  : GHS-POP weighted mean per arrondissement
#   - ndvi_areaw_arr : area-average per arrondissement
#
# Also produces IMU-support variants for robustness:
#   - ndvi_popw_imu  : NDVI mean per IMU polygon, then pop-weighted to arr
#   - ndvi_areaw_imu : NDVI mean per IMU polygon, then area-weighted to arr
#
# Output: paris_arr_ndvi_metrics.csv (same format as original script)
#
# Run BEFORE build_paris_arr_greenness_metrics.R.

if (Sys.info()[["sysname"]] == "Darwin") {
  if (!nzchar(Sys.getenv("PROJ_LIB"))) Sys.setenv(PROJ_LIB = "/opt/homebrew/opt/proj/share/proj")
  if (!nzchar(Sys.getenv("GDAL_DATA"))) Sys.setenv(GDAL_DATA = "/opt/homebrew/opt/gdal/share/gdal")
}

library(data.table)
library(sf)
library(terra)
library(exactextractr)

sf::sf_use_s2(TRUE)
terraOptions(progress = 1)

base_dir <- "/Users/armandeaboudrar-meda/Desktop/GVIestim/paper"
analysis_crs <- 2154  # Lambert-93

# Paths
ndvi_2020_path <- file.path(base_dir, "ESA_WC_NDVI_2020_N48E002.tif")
ndvi_2021_path <- file.path(base_dir, "ESA_WC_NDVI_2021_N48E002.tif")

arr_shp_candidates <- c(
  file.path(base_dir, "arrondissements.shp"),
  file.path(base_dir, "NDVI_Paris_2014_2017", "_tmp_paris_shp", "arrondissements.shp")
)

imu_candidates <- c(
  file.path(base_dir, "Îlots_morphologiques_urbains_d_Île-de-France.shp"),
  file.path(base_dir, "Îlots_morphologiques_urbains_d_Île-de-France.shp"),
  file.path(base_dir, "IMU2022_O.geojson")
)

pop_rast_candidates <- c(
  file.path(base_dir, "GHS_POP_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif"),
  "GHS_POP_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif"
)

out_arr_path <- file.path(base_dir, "paris_arr_ndvi_metrics.csv")

# Helpers
read_first_valid_sf <- function(paths) {
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0) return(NULL)
  for (p in existing) {
    obj <- try(suppressWarnings(st_read(p, quiet = TRUE)), silent = TRUE)
    if (!inherits(obj, "try-error")) return(list(path = p, data = obj))
  }
  NULL
}

safe_wmean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  weighted.mean(x[ok], w[ok])
}

safe_share_nonmissing <- function(x, w) {
  ok_w <- is.finite(w) & w > 0
  ok_x <- is.finite(x)
  denom <- sum(w[ok_w], na.rm = TRUE)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  sum(w[ok_w & ok_x], na.rm = TRUE) / denom
}


# 1) Load ESA WorldCover NDVI tiles and compute pixel-wise average
stopifnot(file.exists(ndvi_2020_path), file.exists(ndvi_2021_path))

cat("Loading ESA WorldCover NDVI tiles...\n")
r2020 <- rast(ndvi_2020_path)
r2021 <- rast(ndvi_2021_path)

cat("  2020:", nlyr(r2020), "bands,", ncol(r2020), "x", nrow(r2020), "pixels\n")
cat("  2021:", nlyr(r2021), "bands,", ncol(r2021), "x", nrow(r2021), "pixels\n")

# Extract Band 1 (p90 = annual peak greenness)
p90_2020 <- r2020[[1]]
p90_2021 <- r2021[[1]]

# Free multi-band rasters
rm(r2020, r2021); gc()

# Pixel-wise average of the two years
cat("Computing pixel-wise average of 2020 and 2021 p90 NDVI...\n")
ndvi_avg <- terra::mean(p90_2020, p90_2021, na.rm = TRUE)
rm(p90_2020, p90_2021); gc()

# 2) Read arrondissement polygons
arr_path <- arr_shp_candidates[file.exists(arr_shp_candidates)][1]
if (is.na(arr_path)) {
  stop("No arrondissement shapefile found. Checked: ",
       paste(arr_shp_candidates, collapse = " ; "))
}

arr_sf <- st_read(arr_path, quiet = TRUE)
stopifnot("c_arinsee" %in% names(arr_sf))
arr_sf$arr <- as.character(arr_sf$c_arinsee)
paris_arr <- sprintf("751%02d", 1:20)
arr_sf <- arr_sf[arr_sf$arr %in% paris_arr, c("arr", "geometry")]
stopifnot(nrow(arr_sf) == 20)

arr_sf_prj <- st_transform(arr_sf, analysis_crs)
arr_sf_prj <- st_make_valid(arr_sf_prj)

# Crop NDVI to Paris bounding box + buffer (saves memory)
paris_bbox_wgs <- st_bbox(st_transform(arr_sf_prj, 4326))
crop_ext <- ext(
  paris_bbox_wgs["xmin"] - 0.02,
  paris_bbox_wgs["xmax"] + 0.02,
  paris_bbox_wgs["ymin"] - 0.02,
  paris_bbox_wgs["ymax"] + 0.02
)
cat("Cropping NDVI raster to Paris extent...\n")
ndvi_paris <- crop(ndvi_avg, crop_ext)
rm(ndvi_avg); gc()

cat("  Cropped raster:", ncol(ndvi_paris), "x", nrow(ndvi_paris), "pixels\n")
cat("  Resolution:", round(res(ndvi_paris)[1] * 111000, 1), "m (approx)\n")

# Mask out invalid values (NDVI should be in [-1, 1])
ndvi_paris[ndvi_paris < -1 | ndvi_paris > 1] <- NA

# Quick sanity check
gstats <- global(ndvi_paris, c("mean", "min", "max"), na.rm = TRUE)
cat("  NDVI stats after crop: mean=", round(gstats$mean, 4),
    " min=", round(gstats$min, 4),
    " max=", round(gstats$max, 4), "\n")

# 3) Read GHS-POP raster and reproject to NDVI grid
pop_rast_path <- pop_rast_candidates[file.exists(pop_rast_candidates)][1]
if (is.na(pop_rast_path)) {
  stop("No GHS population raster found. Checked: ",
       paste(pop_rast_candidates, collapse = " ; "))
}

cat("Reprojecting GHS-POP to NDVI grid (method=sum)...\n")
pop_rast <- rast(pop_rast_path)
pop_on_ndvi <- terra::project(pop_rast, ndvi_paris, method = "sum")
pop_on_ndvi[pop_on_ndvi < 0] <- 0
rm(pop_rast); gc()

# 4) Arrondissement-level NDVI extraction
# Transform arrondissements to raster CRS for exact_extract
arr_extract <- st_transform(arr_sf_prj, crs(ndvi_paris))
arr_extract <- st_make_valid(arr_extract)

cat("Extracting NDVI by arrondissement (area-weighted + GHS-POP weighted)...\n")

ndvi_areaw <- exact_extract(ndvi_paris, arr_extract, "mean", progress = FALSE)
ndvi_popw  <- exact_extract(
  ndvi_paris, arr_extract, "weighted_mean",
  weights = pop_on_ndvi, progress = FALSE
)

arr_ndvi <- data.table(
  arr = arr_extract$arr,
  ndvi_popw_ghs = as.numeric(ndvi_popw),
  ndvi_areaw_arr = as.numeric(ndvi_areaw),
  ndvi_n_years_jjas_ghs = 2L,    # 2020 + 2021
  ndvi_n_dates_jjas_ghs = 2L     # 2 annual composites
)

cat("\nArrondissement NDVI (10m ESA WorldCover, 2020-2021 p90 avg):\n")
print(arr_ndvi[, .(
  ndvi_popw_ghs_min = min(ndvi_popw_ghs, na.rm = TRUE),
  ndvi_popw_ghs_max = max(ndvi_popw_ghs, na.rm = TRUE),
  ndvi_popw_ghs_mean = mean(ndvi_popw_ghs, na.rm = TRUE),
  ndvi_areaw_arr_min = min(ndvi_areaw_arr, na.rm = TRUE),
  ndvi_areaw_arr_max = max(ndvi_areaw_arr, na.rm = TRUE),
  ndvi_areaw_arr_mean = mean(ndvi_areaw_arr, na.rm = TRUE)
)])

# 5) IMU-level NDVI extraction (for robustness metrics)
cat("\nReading IMU polygons...\n")
imu_pick <- read_first_valid_sf(imu_candidates)
if (is.null(imu_pick)) {
  stop("No readable IMU file found. Checked: ",
       paste(imu_candidates, collapse = " ; "))
}
imu <- imu_pick$data
names(imu) <- tolower(names(imu))

if (!("popmen_imu" %in% names(imu))) imu$popmen_imu <- NA_real_
imu <- imu[, c("popmen_imu", "geometry")]

imu_zm <- try(suppressWarnings(st_zm(imu, drop = TRUE, what = "ZM")), silent = TRUE)
if (!inherits(imu_zm, "try-error")) imu <- imu_zm

if (is.na(st_crs(imu))) {
  if (grepl("\\.geojson$", tolower(imu_pick$path))) {
    st_crs(imu) <- 4326
  } else {
    st_crs(imu) <- 2154
  }
}

imu$popmen_imu <- suppressWarnings(as.numeric(imu$popmen_imu))
imu_prj <- st_transform(imu, analysis_crs)
imu_prj <- imu_prj[!st_is_empty(imu_prj), ]
imu_prj <- st_filter(imu_prj, arr_sf_prj, .predicate = st_intersects)
imu_prj <- st_make_valid(imu_prj)
imu_prj$imu_area <- as.numeric(st_area(imu_prj))
imu_prj$imu_id <- seq_len(nrow(imu_prj))

cat("IMU polygons in Paris support:", nrow(imu_prj), "\n")

# Build IMU x arrondissement intersection weights
cat("Building IMU x arrondissement intersection weights...\n")
imu_arr_int <- suppressWarnings(st_intersection(
  imu_prj[, c("imu_id", "popmen_imu", "imu_area", "geometry")],
  arr_sf_prj["arr"]
))
imu_arr_int <- imu_arr_int[!st_is_empty(imu_arr_int), ]
imu_arr_int$part_area <- as.numeric(st_area(imu_arr_int))

imu_arr_key <- as.data.table(st_drop_geometry(imu_arr_int))
imu_arr_key <- imu_arr_key[is.finite(part_area) & part_area > 0]
imu_arr_key[, area_share := fifelse(
  is.finite(imu_area) & imu_area > 0,
  pmin(part_area / imu_area, 1),
  NA_real_
)]
imu_arr_key[, pop_w_part := fifelse(
  is.finite(popmen_imu) & popmen_imu > 0 & is.finite(area_share),
  popmen_imu * area_share,
  NA_real_
)]

imu_arr_key <- imu_arr_key[, .(
  area_w = sum(part_area, na.rm = TRUE),
  pop_w = if (any(is.finite(pop_w_part) & pop_w_part > 0)) {
    sum(pop_w_part[is.finite(pop_w_part) & pop_w_part > 0], na.rm = TRUE)
  } else {
    NA_real_
  }
), by = .(imu_id, arr)]

# Extract mean NDVI per IMU polygon
imu_extract <- st_transform(imu_prj, crs(ndvi_paris))
imu_extract <- st_make_valid(imu_extract)

cat("Extracting NDVI by IMU polygon...\n")
imu_ndvi_vals <- exact_extract(ndvi_paris, imu_extract, "mean", progress = FALSE)

imu_mean <- data.table(
  imu_id = imu_prj$imu_id,
  ndvi_jjas_2020_2021 = as.numeric(imu_ndvi_vals)
)

# Aggregate IMU to arrondissement
imu_for_arr <- merge(imu_arr_key, imu_mean, by = "imu_id", all.x = TRUE)

arr_ndvi_imu <- imu_for_arr[, .(
  ndvi_popw_imu = safe_wmean(ndvi_jjas_2020_2021, pop_w),
  ndvi_areaw_imu = safe_wmean(ndvi_jjas_2020_2021, area_w),
  ndvi_popw_coverage = safe_share_nonmissing(ndvi_jjas_2020_2021, pop_w),
  ndvi_areaw_coverage = safe_share_nonmissing(ndvi_jjas_2020_2021, area_w),
  ndvi_n_years_jjas_imu = 2L,
  ndvi_n_dates_jjas_imu = 2L
), by = arr]

# 6) Merge and save
arr_ndvi <- merge(
  data.table(arr = paris_arr),
  arr_ndvi,
  by = "arr", all.x = TRUE
)
arr_ndvi <- merge(arr_ndvi, arr_ndvi_imu, by = "arr", all.x = TRUE)
setorder(arr_ndvi, arr)

fwrite(arr_ndvi, out_arr_path)

cat("\nSaved:", out_arr_path, "\n")

cat("\nArrondissement NDVI final summary:\n")
print(arr_ndvi[, .(
  ndvi_popw_ghs_min = min(ndvi_popw_ghs, na.rm = TRUE),
  ndvi_popw_ghs_max = max(ndvi_popw_ghs, na.rm = TRUE),
  ndvi_areaw_arr_min = min(ndvi_areaw_arr, na.rm = TRUE),
  ndvi_areaw_arr_max = max(ndvi_areaw_arr, na.rm = TRUE),
  ndvi_n_years_jjas_ghs_min = min(ndvi_n_years_jjas_ghs, na.rm = TRUE),
  ndvi_n_years_jjas_ghs_max = max(ndvi_n_years_jjas_ghs, na.rm = TRUE),
  ndvi_popw_imu_min = min(ndvi_popw_imu, na.rm = TRUE),
  ndvi_popw_imu_max = max(ndvi_popw_imu, na.rm = TRUE),
  ndvi_areaw_imu_min = min(ndvi_areaw_imu, na.rm = TRUE),
  ndvi_areaw_imu_max = max(ndvi_areaw_imu, na.rm = TRUE)
)])

cat("\nAll 20 arrondissements:\n")
print(arr_ndvi)

cat("\nDone. This replaces the output of build_paris_arr_ndvi_metrics.R.\n")
cat("Next step: run build_paris_arr_greenness_metrics.R to merge into unified file.\n")
