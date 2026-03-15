# build_paris_arr_greenness_metrics.R
#
# Build arrondissement-level greenness metrics from GVI points, IMU
# polygons, and NDVI summaries. Produces BOTH the GHS-POP pop-weighted
# metrics (main paper) and the native/unweighted metrics (appendix),
# all in a single output file.
#
# Main paper metrics (GHS-POP weighted):
#   - gvi_popw_points  : GVI point-level weighted by GHS-POP pixel value
#   - ndvi_popw_ghs    : NDVI (from build_paris_arr_ndvi_metrics.R)
#   - imu_*_popw_ghs   : IMU vegetation weighted by GHS-POP within
#                         arrondissement-allocated IMU polygon fragments
#
# Native/appendix metrics (area-weighted or unweighted):
#   - gvi_mean         : simple point mean per arrondissement
#   - ndvi_native_jjas : area-average NDVI
#   - imu_veg_total    : area-weighted IMU vegetation
#
# Output:
#   - paris_arr_greenness_metrics.csv

if (Sys.info()[["sysname"]] == "Darwin") {
  if (!nzchar(Sys.getenv("PROJ_LIB"))) Sys.setenv(PROJ_LIB = "/opt/homebrew/opt/proj/share/proj")
  if (!nzchar(Sys.getenv("GDAL_DATA"))) Sys.setenv(GDAL_DATA = "/opt/homebrew/opt/gdal/share/gdal")
}

library(data.table)
library(sf)
library(terra)
library(exactextractr)

sf::sf_use_s2(TRUE)

base_dir <- "/Users/armandeaboudrar-meda/Desktop/GVIestim/paper" # relative:   base_dir <- here::here()
analysis_crs <- 2154   # Lambert-93

# 0) Paths
gvi_csv_path <- file.path(base_dir, "gvi_356_cities.csv")
ndvi_native_path <- file.path(base_dir, "paris_ndvi_2014_2017_by_arr_mean.csv")
ndvi_imu_path <- file.path(base_dir, "paris_arr_ndvi_metrics.csv")  

imu_candidates <- c(
  file.path(base_dir, "Îlots_morphologiques_urbains_d_Île-de-France.shp"),
  file.path(base_dir, "Îlots_morphologiques_urbains_d_Île-de-France.shp"),
  file.path(base_dir, "IMU2022_O.geojson")
)

read_first_valid_sf <- function(paths, label = "sf object") {
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0) return(NULL)
  for (p in existing) {
    obj <- try(suppressWarnings(st_read(p, quiet = TRUE)), silent = TRUE)
    if (!inherits(obj, "try-error")) {
      return(list(path = p, data = obj))
    }
  }
  NULL
}

imu_pick <- read_first_valid_sf(imu_candidates, "IMU")
if (is.null(imu_pick)) {
  stop("No readable IMU file found. Checked: ", paste(imu_candidates, collapse = " ; "))
}
imu_path <- imu_pick$path

out_path <- file.path(base_dir, "paris_arr_greenness_metrics.csv")

pop_rast_candidates <- c(
  file.path(base_dir, "GHS_POP_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif"),
  "GHS_POP_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif"
)

pop_rast_path <- pop_rast_candidates[file.exists(pop_rast_candidates)][1]
if (is.na(pop_rast_path)) {
  stop("No population raster found for point weighting. Put GHS_POP_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif in base_dir or working directory.")
}
pop_rast <- rast(pop_rast_path)

# 0b) Optional NDVI inputs
ndvi_native <- NULL
if (file.exists(ndvi_native_path)) {
  ndvi_native <- fread(ndvi_native_path)
  ndvi_native[, arr := as.character(arr)]
  ndvi_native <- ndvi_native[, .(
    arr,
    ndvi_native_all = ndvi_mean_2014_2017_all,
    ndvi_native_jjas = ndvi_mean_2014_2017_jjas,
    ndvi_native_sd_all = ndvi_sd_2014_2017_all,
    ndvi_native_sd_jjas = ndvi_sd_2014_2017_jjas,
    ndvi_n_dates_all = n_dates_all,
    ndvi_n_dates_jjas = n_dates_jjas
  )]
  setorder(ndvi_native, arr)
}

ndvi_imu <- NULL
if (file.exists(ndvi_imu_path)) {
  ndvi_imu <- fread(ndvi_imu_path)
  ndvi_imu[, arr := as.character(arr)]
  keep_ndvi_imu <- intersect(
    c("arr",
      "ndvi_popw_ghs", "ndvi_areaw_arr",
      "ndvi_n_dates_jjas_ghs", "ndvi_n_years_jjas_ghs",
      "ndvi_popw_imu", "ndvi_areaw_imu",
      "ndvi_n_dates_jjas_imu", "ndvi_n_years_jjas_imu"),
    names(ndvi_imu)
  )
  ndvi_imu <- ndvi_imu[, ..keep_ndvi_imu]
  setorder(ndvi_imu, arr)
}

# helper
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

# 1) Read Paris arrondissement polygons 
arr_shp_candidates <- c(
  file.path(base_dir, "arrondissements.shp"),
  file.path(base_dir, "NDVI_Paris_2014_2017", "_tmp_paris_shp", "arrondissements.shp")
)

arr_shp_path <- arr_shp_candidates[file.exists(arr_shp_candidates)][1]
if (is.na(arr_shp_path)) {
  stop("No local arrondissement shapefile found. Expected one of: ",
       paste(arr_shp_candidates, collapse = " ; "))
}

arr_sf <- st_read(arr_shp_path, quiet = TRUE)
if (!("c_arinsee" %in% names(arr_sf))) {
  stop("Local arrondissement shapefile must contain 'c_arinsee'. File: ", arr_shp_path)
}

arr_sf$arr <- as.character(arr_sf$c_arinsee)
paris_arr <- sprintf("751%02d", 1:20)
arr_sf <- arr_sf[arr_sf$arr %in% paris_arr, c("arr", "geometry")]
stopifnot(nrow(arr_sf) == 20)

arr_sf_prj <- st_transform(arr_sf, analysis_crs)
arr_sf_prj <- st_make_valid(arr_sf_prj)

# 2) Read IMU and prepare vegetation variables
imu <- imu_pick$data
names(imu) <- tolower(names(imu))

needed_cols <- c("iv_haute", "iv_basse", "geometry")
missing_cols <- setdiff(needed_cols, names(imu))
if (length(missing_cols) > 0) {
  stop("IMU file is missing required columns: ", paste(missing_cols, collapse = ", "))
}

if (!("popmen_imu" %in% names(imu))) {
  imu$popmen_imu <- NA_real_
}

imu <- imu[, c("iv_haute", "iv_basse", "popmen_imu", "geometry")]

# Drop Z/M if present
imu_zm <- try(suppressWarnings(st_zm(imu, drop = TRUE, what = "ZM")), silent = TRUE)
if (!inherits(imu_zm, "try-error")) imu <- imu_zm

# If CRS is missing, infer it:
# - shapefile IMU is usually projected (Lambert-93 / EPSG:2154)
# - geojson is usually lon/lat (EPSG:4326)
if (is.na(st_crs(imu))) {
  if (grepl("\\.geojson$", tolower(imu_path))) {
    st_crs(imu) <- 4326
  } else {
    st_crs(imu) <- 2154
  }
}

imu_high_raw <- suppressWarnings(as.numeric(imu$iv_haute))
imu_low_raw  <- suppressWarnings(as.numeric(imu$iv_basse))
imu_pop_raw  <- suppressWarnings(as.numeric(imu$popmen_imu))

# IMU values are usually encoded as 0-1 fractions even though conceptually they are %
scale_factor <- if (max(c(imu_high_raw, imu_low_raw), na.rm = TRUE) <= 1.5) 100 else 1

imu$imu_high <- imu_high_raw * scale_factor
imu$imu_low  <- imu_low_raw  * scale_factor
imu$imu_veg_total <- imu$imu_high + imu$imu_low
imu$popmen_imu <- imu_pop_raw

imu_prj <- st_transform(imu, analysis_crs)
imu_prj <- imu_prj[!st_is_empty(imu_prj), ]

# Keep only IMU polygons near/intersecting Paris
imu_prj <- st_filter(imu_prj, arr_sf_prj, .predicate = st_intersects)

# IMU area for area-weighted averages
imu_prj$imu_area <- as.numeric(st_area(imu_prj))

# Stable polygon id for downstream joins
imu_prj$imu_id <- seq_len(nrow(imu_prj))

# 3) Allocate IMU polygons to arrondissements by area intersection
cat("Building IMU x arrondissement area intersections...\n")

imu_prj <- st_make_valid(imu_prj)
imu_arr_int <- suppressWarnings(st_intersection(
  imu_prj[, c("imu_id", "popmen_imu", "imu_area", "geometry")],
  arr_sf_prj["arr"]
))
imu_arr_int <- imu_arr_int[!st_is_empty(imu_arr_int), ]
imu_arr_int$part_area <- as.numeric(st_area(imu_arr_int))

# Alternative population weights based on the same GHS raster used for GVI points.
# This keeps IMU vegetation on IMU support, but harmonizes the population source.
imu_arr_int_wgs <- st_transform(imu_arr_int, 4326)
imu_arr_int$pop_ghs_part <- suppressWarnings(as.numeric(
  exact_extract(pop_rast, imu_arr_int_wgs, "sum")
))
imu_arr_int$pop_ghs_part[!is.finite(imu_arr_int$pop_ghs_part) | imu_arr_int$pop_ghs_part <= 0] <- NA_real_

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

# Collapse multipart intersections to one row per IMU x arrondissement
imu_arr_key <- imu_arr_key[, .(
  area_w = sum(part_area, na.rm = TRUE),
  pop_w = if (any(is.finite(pop_w_part) & pop_w_part > 0)) {
    sum(pop_w_part[is.finite(pop_w_part) & pop_w_part > 0], na.rm = TRUE)
  } else {
    NA_real_
  },
  pop_w_ghs = if (any(is.finite(pop_ghs_part) & pop_ghs_part > 0)) {
    sum(pop_ghs_part[is.finite(pop_ghs_part) & pop_ghs_part > 0], na.rm = TRUE)
  } else {
    NA_real_
  }
), by = .(imu_id, arr)]

imu_attr_dt <- as.data.table(st_drop_geometry(imu_prj))[, .(
  imu_id, imu_high, imu_low, imu_veg_total, popmen_imu, imu_area
)]

imu_poly_dt <- merge(imu_arr_key, imu_attr_dt, by = "imu_id", all.x = TRUE)
setorder(imu_poly_dt, arr, imu_id)

# Arrondissement IMU metrics (area- and population-weighted after area allocation)
imu_arr <- imu_poly_dt[, .(
  imu_high = safe_wmean(imu_high, area_w),
  imu_low = safe_wmean(imu_low, area_w),
  imu_veg_total = safe_wmean(imu_veg_total, area_w),
  
  imu_high_popw = safe_wmean(imu_high, pop_w),
  imu_low_popw = safe_wmean(imu_low, pop_w),
  imu_veg_total_popw = safe_wmean(imu_veg_total, pop_w),

  imu_high_popw_ghs = safe_wmean(imu_high, pop_w_ghs),
  imu_low_popw_ghs = safe_wmean(imu_low, pop_w_ghs),
  imu_veg_total_popw_ghs = safe_wmean(imu_veg_total, pop_w_ghs),

  imu_popw_ghs_coverage = safe_share_nonmissing(imu_veg_total, pop_w_ghs),
  
  imu_n_parts = .N
), by = arr]

setorder(imu_arr, arr)

# 4) GVI aggregation
#    4A) keep legacy arrondissement point-mean GVI
#    4B) + harmonised GVI on IMU support, then area/pop aggregation to arr
cat("Reading GVI CSV...\n")
gvi <- fread(gvi_csv_path, showProgress = TRUE)
names(gvi) <- tolower(names(gvi))

required_gvi <- c("city", "x", "y", "gvi")
missing_gvi <- setdiff(required_gvi, names(gvi))
if (length(missing_gvi) > 0) {
  stop("GVI CSV is missing required columns: ", paste(missing_gvi, collapse = ", "))
}

gvi <- gvi[tolower(city) == "paris"]
gvi <- gvi[!is.na(x) & !is.na(y) & !is.na(gvi)]

if ("lcz_filter_v3" %in% names(gvi)) {
  gvi <- gvi[lcz_filter_v3 %in% 1:10]
}

# Collapse duplicate coordinates
gvi <- gvi[, .(gvi = mean(gvi, na.rm = TRUE)), by = .(x, y)]

cat("Paris GVI points after filtering:", nrow(gvi), "\n")

gvi_sf <- st_as_sf(gvi, coords = c("x", "y"), crs = 4326, remove = FALSE)
gvi_sf_prj <- st_transform(gvi_sf, analysis_crs)

# GHS raster is in lon/lat, so extract using WGS84 points
gvi_sf_wgs <- st_transform(gvi_sf_prj, 4326)

gvi_sf_prj$pt_pop_ghs <- terra::extract(pop_rast, vect(gvi_sf_wgs))[, 2]
gvi_sf_prj$pt_pop_ghs <- suppressWarnings(as.numeric(gvi_sf_prj$pt_pop_ghs))
gvi_sf_prj$pt_pop_ghs[!is.finite(gvi_sf_prj$pt_pop_ghs) | gvi_sf_prj$pt_pop_ghs < 0] <- NA_real_

# 4A) legacy arrondissement-level point mean (keep for robustness)
gvi_join_arr <- st_join(gvi_sf_prj, arr_sf_prj["arr"], left = FALSE)

gvi_arr_legacy <- as.data.table(st_drop_geometry(gvi_join_arr))[, .(
  gvi_mean = mean(gvi, na.rm = TRUE),
  gvi_n_points_legacy = .N
), by = arr]

setorder(gvi_arr_legacy, arr)

# 4B) NEW PRIMARY: true point-level population-weighted GVI
# Join points to arrondissements and weight by point-level GHS population
gvi_join_arr_pop <- st_join(
  gvi_sf_prj[, c("gvi", "pt_pop_ghs")],
  arr_sf_prj["arr"],
  left = FALSE
)

gvi_arr_pointpopw <- as.data.table(st_drop_geometry(gvi_join_arr_pop))[, .(
  gvi_popw_points = safe_wmean(gvi, pt_pop_ghs),
  gvi_n_points_popw = sum(is.finite(pt_pop_ghs) & pt_pop_ghs > 0),
  gvi_n_points_total = .N,
  gvi_points_popweight_coverage = mean(is.finite(pt_pop_ghs) & pt_pop_ghs > 0),
  gvi_pointpop_sum = sum(pt_pop_ghs[is.finite(pt_pop_ghs) & pt_pop_ghs > 0], na.rm = TRUE)
), by = arr]

setorder(gvi_arr_pointpopw, arr)

# 4C) HARMONISATION ROBUSTNESS: GVI -> IMU polygon mean -> arrondissement
gvi_join_imu <- st_join(
  gvi_sf_prj,
  imu_prj[, "imu_id"],
  left = FALSE,
  largest = TRUE
)

gvi_point_dt <- as.data.table(st_drop_geometry(gvi_join_imu))

gvi_poly_dt <- gvi_point_dt[, .(
  gvi_poly = mean(gvi, na.rm = TRUE),
  gvi_n_points = .N
), by = imu_id]

imu_gvi_dt <- merge(
  imu_poly_dt[, .(imu_id, arr, area_w, pop_w)],
  gvi_poly_dt,
  by = "imu_id",
  all.x = TRUE
)

gvi_arr_harmonized <- imu_gvi_dt[, .(
  gvi_areaw_imu = safe_wmean(gvi_poly, area_w),
  gvi_popw_imu  = safe_wmean(gvi_poly, pop_w),
  
  gvi_n_imu_polys_with_points = uniqueN(imu_id[is.finite(gvi_poly)]),
  gvi_n_points = sum(gvi_n_points, na.rm = TRUE),
  
  gvi_areaw_coverage = safe_share_nonmissing(gvi_poly, area_w),
  gvi_popw_coverage  = safe_share_nonmissing(gvi_poly, pop_w)
), by = arr]

setorder(gvi_arr_harmonized, arr)

# 5) Merge and save
out <- data.table(arr = paris_arr)
out <- merge(out, imu_arr, by = "arr", all.x = TRUE)
out <- merge(out, gvi_arr_legacy, by = "arr", all.x = TRUE)
out <- merge(out, gvi_arr_pointpopw, by = "arr", all.x = TRUE)
out <- merge(out, gvi_arr_harmonized, by = "arr", all.x = TRUE)
if (!is.null(ndvi_native)) {
  out <- merge(out, ndvi_native, by = "arr", all.x = TRUE)
}
if (!is.null(ndvi_imu)) {
  out <- merge(out, ndvi_imu, by = "arr", all.x = TRUE)
}
setorder(out, arr)

fwrite(out, out_path)

cat("\nSaved:\n -", out_path, "\n")

# 6) Sanity checks
cat("\nSummary ranges:\n")
print(out[, .(
  imu_high_min = min(imu_high, na.rm = TRUE),
  imu_high_max = max(imu_high, na.rm = TRUE),
  imu_low_min = min(imu_low, na.rm = TRUE),
  imu_low_max = max(imu_low, na.rm = TRUE),
  imu_veg_total_min = min(imu_veg_total, na.rm = TRUE),
  imu_veg_total_max = max(imu_veg_total, na.rm = TRUE),
  
  imu_high_popw_min = min(imu_high_popw, na.rm = TRUE),
  imu_high_popw_max = max(imu_high_popw, na.rm = TRUE),
  imu_low_popw_min = min(imu_low_popw, na.rm = TRUE),
  imu_low_popw_max = max(imu_low_popw, na.rm = TRUE),
  imu_veg_total_popw_min = min(imu_veg_total_popw, na.rm = TRUE),
  imu_veg_total_popw_max = max(imu_veg_total_popw, na.rm = TRUE),

  imu_veg_total_popw_ghs_min = min(imu_veg_total_popw_ghs, na.rm = TRUE),
  imu_veg_total_popw_ghs_max = max(imu_veg_total_popw_ghs, na.rm = TRUE),
  
  gvi_mean_min = min(gvi_mean, na.rm = TRUE),
  gvi_mean_max = max(gvi_mean, na.rm = TRUE),
  
  gvi_popw_points_min = min(gvi_popw_points, na.rm = TRUE),
  gvi_popw_points_max = max(gvi_popw_points, na.rm = TRUE),
  
  gvi_areaw_imu_min = min(gvi_areaw_imu, na.rm = TRUE),
  gvi_areaw_imu_max = max(gvi_areaw_imu, na.rm = TRUE),
  gvi_popw_imu_min = min(gvi_popw_imu, na.rm = TRUE),
  gvi_popw_imu_max = max(gvi_popw_imu, na.rm = TRUE),
  
  ndvi_native_jjas_min = min(ndvi_native_jjas, na.rm = TRUE),
  ndvi_native_jjas_max = max(ndvi_native_jjas, na.rm = TRUE),
  ndvi_native_all_min = min(ndvi_native_all, na.rm = TRUE),
  ndvi_native_all_max = max(ndvi_native_all, na.rm = TRUE),
  
  gvi_points_popweight_coverage_min = min(gvi_points_popweight_coverage, na.rm = TRUE),
  gvi_points_popweight_coverage_max = max(gvi_points_popweight_coverage, na.rm = TRUE),
  gvi_pointpop_sum_min = min(gvi_pointpop_sum, na.rm = TRUE),
  gvi_pointpop_sum_max = max(gvi_pointpop_sum, na.rm = TRUE)
)])

if ("ndvi_popw_imu" %in% names(out)) {
  cat("\nNDVI weighted ranges:\n")
  print(out[, .(
    ndvi_popw_ghs_min = min(ndvi_popw_ghs, na.rm = TRUE),
    ndvi_popw_ghs_max = max(ndvi_popw_ghs, na.rm = TRUE),
    ndvi_areaw_arr_min = min(ndvi_areaw_arr, na.rm = TRUE),
    ndvi_areaw_arr_max = max(ndvi_areaw_arr, na.rm = TRUE),
    ndvi_popw_imu_min = min(ndvi_popw_imu, na.rm = TRUE),
    ndvi_popw_imu_max = max(ndvi_popw_imu, na.rm = TRUE),
    ndvi_areaw_imu_min = min(ndvi_areaw_imu, na.rm = TRUE),
    ndvi_areaw_imu_max = max(ndvi_areaw_imu, na.rm = TRUE)
  )])
}

cat("\nFirst rows:\n")
print(out)

cat("\nCorrelation checks:\n")
if (all(is.finite(out$gvi_mean)) && all(is.finite(out$imu_veg_total_popw))) {
  cat("  cor(gvi_mean, imu_veg_total_popw) =",
      round(cor(out$gvi_mean, out$imu_veg_total_popw), 3), "\n")
}
if (all(is.finite(out$gvi_popw_imu)) && all(is.finite(out$imu_veg_total_popw))) {
  cat("  cor(gvi_popw_imu, imu_veg_total_popw) =",
      round(cor(out$gvi_popw_imu, out$imu_veg_total_popw), 3), "\n")
}
if (all(is.finite(out$gvi_areaw_imu)) && all(is.finite(out$imu_veg_total))) {
  cat("  cor(gvi_areaw_imu, imu_veg_total) =",
      round(cor(out$gvi_areaw_imu, out$imu_veg_total), 3), "\n")
}

if (all(is.finite(out$gvi_popw_points)) && all(is.finite(out$imu_veg_total_popw))) {
  cat("  cor(gvi_popw_points, imu_veg_total_popw) =",
      round(cor(out$gvi_popw_points, out$imu_veg_total_popw), 3), "\n")
}
if (all(is.finite(out$gvi_popw_points)) && all(is.finite(out$imu_veg_total_popw_ghs))) {
  cat("  cor(gvi_popw_points, imu_veg_total_popw_ghs) =",
      round(cor(out$gvi_popw_points, out$imu_veg_total_popw_ghs), 3), "\n")
}
