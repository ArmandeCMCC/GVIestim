# build_paris_arr_ndvi_metrics.R
# compute NDVI-derived arrondissement metrics before
#   build_paris_arr_greenness_metrics.R merges all greenness families

#   Outputs from this script feed ndvi_* fields used in the
#   _final_nativeheat_ghsgreen pipeline.
#
# Build NDVI metrics for Paris from JJAS 2014-2017 TIFFs:
#   - ndvi_popw_ghs  : NDVI extracted to arrondissement and weighted with GHS POP
#                      (same population source as gvi_popw_points and imu_*_popw_ghs)
#   - ndvi_areaw_arr : NDVI extracted to arrondissement (area-weighted / unweighted mean)
#   - ndvi_popw_imu  : NDVI extracted to IMU polygons, then IMU-pop weighted
#                      to arrondissement (IMU support robustness)
#   - ndvi_areaw_imu : NDVI extracted to IMU polygons, then area-weighted
#                      to arrondissement (IMU support robustness)
#
# Design choice:
#   - Keep "native NDVI" as main NDVI metric in the paper
#   - Add population-harmonized NDVI with GHS POP at arrondissement level
#   - Keep IMU-supported NDVI as support-sensitivity
#
# Inputs:
#   - cdse_ndvi_paris/paris_ndvi_2014_2017_cog_filelist.csv
#   - cdse_ndvi_paris/ndvi300_paris_jjas_2014_2017_tifs/ndvi300_paris_YYYYMMDD.tif
#     (JJAS local TIFFs; preferred)
#   - GHS_POP_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif
#   - IMU polygons
#   - local arrondissement shapefile
#
# Optional:
#   - If local clips are missing and download_missing = TRUE, uses AWS CLI to
#     stream missing dates from CDSE object storage.
#
# Outputs:
#   - paris_arr_ndvi_metrics.csv
#   - paris_ndvi_arr_jjas_2014_2017_by_arr_date.csv
#   - paris_ndvi_imu_jjas_2014_2017_by_imu.csv
#   - paris_ndvi_imu_jjas_2014_2017_partial.csv (checkpoint for resume)
#   - paris_ndvi_arr_jjas_2014_2017_partial.csv (checkpoint for resume)
#   - paris_ndvi_imu_jjas_2014_2017_failed_dates.csv

if (Sys.info()[["sysname"]] == "Darwin") {
  if (!nzchar(Sys.getenv("PROJ_LIB"))) Sys.setenv(PROJ_LIB = "/opt/homebrew/opt/proj/share/proj")
  if (!nzchar(Sys.getenv("GDAL_DATA"))) Sys.setenv(GDAL_DATA = "/opt/homebrew/opt/gdal/share/gdal")
}

library(data.table)
library(sf)
library(terra)
library(exactextractr)
library(lubridate)

sf::sf_use_s2(TRUE)
terraOptions(progress = 1)

base_dir <- "/Users/armandeaboudrar-meda/Desktop/GVIestim/paper"
analysis_crs <- 2154

# runtime controls
download_missing <- FALSE
endpoint_url <- "https://eodata.dataspace.copernicus.eu"
min_dates_required <- 30

# paths
filelist_path <- file.path(base_dir, "cdse_ndvi_paris", "paris_ndvi_2014_2017_cog_filelist.csv")
local_clip_candidates <- c(
  file.path(base_dir, "cdse_ndvi_paris", "ndvi300_paris_jjas_2014_2017_tifs"),
  file.path(base_dir, "cdse_ndvi_paris", "ndvi_paris_clips_2014_2017")
)
local_clip_dir <- local_clip_candidates[dir.exists(local_clip_candidates)][1]
if (is.na(local_clip_dir)) {
  stop("No local NDVI TIFF directory found. Checked: ",
       paste(local_clip_candidates, collapse = " ; "))
}

checkpoint_path <- file.path(base_dir, "paris_ndvi_imu_jjas_2014_2017_partial.csv")
checkpoint_arr_path <- file.path(base_dir, "paris_ndvi_arr_jjas_2014_2017_partial.csv")
failed_path <- file.path(base_dir, "paris_ndvi_imu_jjas_2014_2017_failed_dates.csv")

out_arr_path <- file.path(base_dir, "paris_arr_ndvi_metrics.csv")
out_imu_path <- file.path(base_dir, "paris_ndvi_imu_jjas_2014_2017_by_imu.csv")
out_arr_date_path <- file.path(base_dir, "paris_ndvi_arr_jjas_2014_2017_by_arr_date.csv")

arr_shp_candidates <- c(
  file.path(base_dir, "arrondissements.shp"),
  file.path(base_dir, "NDVI_Paris_2014_2017", "_tmp_paris_shp", "arrondissements.shp")
)

imu_candidates <- c(
  file.path(base_dir, "Îlots_morphologiques_urbains_d_Île-de-France.shp"),
  file.path(base_dir, "Îlots_morphologiques_urbains_d_Île-de-France.shp"),
  file.path(base_dir, "IMU2022_O.geojson")
)

pop_rast_candidates <- c(
  file.path(base_dir, "GHS_POP_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif"),
  "GHS_POP_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif"
)

# helpers
read_first_valid_sf <- function(paths) {
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

ndvi_object_from_prefix <- function(x) {
  x <- sub("^s3://EODATA/", "s3://eodata/", x)
  x <- sub("/$", "", x)
  prod_dir <- dirname(x)
  prod_name <- basename(x)
  tif_name <- sub("^c_gls_NDVI300_", "c_gls_NDVI300-NDVI_", prod_name)
  tif_name <- sub("_cog$", ".tiff", tif_name)
  paste0(prod_dir, "/", prod_name, "/", tif_name)
}

load_ndvi_raster <- function(date_chr, s3_prefix, download_missing, endpoint_url) {
  ymd <- gsub("-", "", date_chr)
  local_path <- file.path(local_clip_dir, paste0("ndvi300_paris_", ymd, ".tif"))
  if (file.exists(local_path)) {
    return(list(r = rast(local_path), tmp = NULL, source = "local_clip"))
  }

  if (!download_missing) {
    return(NULL)
  }

  aws_bin <- Sys.which("aws")
  if (!nzchar(aws_bin)) {
    stop("AWS CLI not found. Install aws cli or set download_missing <- FALSE.")
  }
  if (!nzchar(Sys.getenv("AWS_ACCESS_KEY_ID")) || !nzchar(Sys.getenv("AWS_SECRET_ACCESS_KEY"))) {
    stop("Missing AWS credentials in environment for download_missing = TRUE.")
  }

  object_path <- ndvi_object_from_prefix(s3_prefix)
  tmp_src <- tempfile(pattern = paste0("ndvi_", ymd, "_"), fileext = ".tiff")

  status <- system2(
    aws_bin,
    args = c(
      "s3", "cp",
      object_path,
      tmp_src,
      "--endpoint-url", endpoint_url,
      "--only-show-errors"
    )
  )

  if (!identical(status, 0L) || !file.exists(tmp_src)) {
    if (file.exists(tmp_src)) unlink(tmp_src)
    return(NULL)
  }

  list(r = rast(tmp_src), tmp = tmp_src, source = "remote_tmp")
}

normalize_ndvi <- function(r) {
  gmax <- suppressWarnings(global(r, "max", na.rm = TRUE)[1, 1])
  gmin <- suppressWarnings(global(r, "min", na.rm = TRUE)[1, 1])

  # If raster is DN-coded (0..255), convert to physical NDVI
  if (is.finite(gmax) && gmax > 2) {
    r <- classify(r, matrix(c(251, 255, NA), ncol = 3, byrow = TRUE), include.lowest = TRUE)
    r <- r * 0.004 - 0.08
  }

  r[r < -1 | r > 1] <- NA
  r
}

# 1) read target JJAS date list
if (!file.exists(filelist_path)) {
  stop("Missing NDVI file list: ", filelist_path)
}

files <- fread(filelist_path)
if (!all(c("date", "s3_path_cli") %in% names(files))) {
  stop("NDVI file list must contain columns: date, s3_path_cli")
}

files[, date := as.Date(date)]
files <- files[year(date) %in% 2014:2017 & month(date) %in% 6:9]
setorder(files, date)
if (nrow(files) == 0) stop("No JJAS dates found in file list.")

# 2) read arrondissement polygons 
arr_path <- arr_shp_candidates[file.exists(arr_shp_candidates)][1]
if (is.na(arr_path)) {
  stop("No arrondissement shapefile found. Checked: ", paste(arr_shp_candidates, collapse = " ; "))
}

arr_sf <- st_read(arr_path, quiet = TRUE)
if (!("c_arinsee" %in% names(arr_sf))) {
  stop("Arrondissement shapefile must contain c_arinsee: ", arr_path)
}
arr_sf$arr <- as.character(arr_sf$c_arinsee)
paris_arr <- sprintf("751%02d", 1:20)
arr_sf <- arr_sf[arr_sf$arr %in% paris_arr, c("arr", "geometry")]
stopifnot(nrow(arr_sf) == 20)

arr_sf_prj <- st_transform(arr_sf, analysis_crs)
arr_sf_prj <- st_make_valid(arr_sf_prj)

# oopulation raster used for harmonised NDVI weighting (same source as GVI pop-w.)
pop_rast_path <- pop_rast_candidates[file.exists(pop_rast_candidates)][1]
if (is.na(pop_rast_path)) {
  stop("No GHS population raster found. Checked: ",
       paste(pop_rast_candidates, collapse = " ; "))
}
pop_rast <- rast(pop_rast_path)

# 3) read IMU and build arrondissement allocation weights (same support as IMU)
imu_pick <- read_first_valid_sf(imu_candidates)
if (is.null(imu_pick)) {
  stop("No readable IMU file found. Checked: ", paste(imu_candidates, collapse = " ; "))
}
imu <- imu_pick$data
names(imu) <- tolower(names(imu))

needed_cols <- c("iv_haute", "iv_basse", "geometry")
missing_cols <- setdiff(needed_cols, names(imu))
if (length(missing_cols) > 0) {
  stop("IMU file is missing required columns: ", paste(missing_cols, collapse = ", "))
}
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

# 4) resume checkpoint (IMU x date NDVI)
if (file.exists(checkpoint_path)) {
  ndvi_imu_date <- fread(checkpoint_path)
  ndvi_imu_date[, date := as.Date(date)]
} else {
  ndvi_imu_date <- data.table(
    imu_id = integer(),
    date = as.Date(character()),
    mean_ndvi_imu = numeric()
  )
}

if (file.exists(checkpoint_arr_path)) {
  ndvi_arr_date <- fread(checkpoint_arr_path)
  ndvi_arr_date[, date := as.Date(date)]
} else {
  ndvi_arr_date <- data.table(
    arr = character(),
    date = as.Date(character()),
    ndvi_areaw_arr_date = numeric(),
    ndvi_popw_ghs_date = numeric()
  )
}

done_dates_imu <- unique(as.character(ndvi_imu_date$date))
done_dates_arr <- unique(as.character(ndvi_arr_date$date))
done_dates <- intersect(done_dates_imu, done_dates_arr)
files_todo <- files[!(as.character(date) %in% done_dates)]

cat("Already done dates (both IMU + arrondissement):", length(done_dates), "\n")
cat("Remaining dates:", nrow(files_todo), "\n")

failed_dates <- character()
imu_extract <- NULL
imu_extract_crs <- NA_character_
arr_extract <- NULL
arr_extract_crs <- NA_character_
pop_on_ndvi <- NULL
pop_on_ndvi_key <- NA_character_

# 5) loop over missing dates
for (i in seq_len(nrow(files_todo))) {
  date_i <- as.character(files_todo$date[i])
  s3_i <- files_todo$s3_path_cli[i]

  cat(sprintf("[%03d/%03d] %s\n", i, nrow(files_todo), date_i))

  rr <- load_ndvi_raster(date_i, s3_i, download_missing, endpoint_url)
  if (is.null(rr)) {
    cat("  failed: missing local clip and no successful remote fetch\n")
    failed_dates <- c(failed_dates, date_i)
    next
  }

  r <- normalize_ndvi(rr$r)

  # cache IMU/arrondissement transformations in the raster CRS
  r_crs <- terra::crs(r, proj = TRUE)
  if (!is.character(r_crs) || !nzchar(r_crs)) {
    if (!is.null(rr$tmp) && file.exists(rr$tmp)) unlink(rr$tmp)
    failed_dates <- c(failed_dates, date_i)
    next
  }

  if (is.na(imu_extract_crs) || !identical(r_crs, imu_extract_crs)) {
    imu_extract <- st_transform(imu_prj, r_crs)
    imu_extract <- st_make_valid(imu_extract)
    imu_extract_crs <- r_crs
  }

  if (is.na(arr_extract_crs) || !identical(r_crs, arr_extract_crs)) {
    arr_extract <- st_transform(arr_sf_prj, r_crs)
    arr_extract <- st_make_valid(arr_extract)
    arr_extract_crs <- r_crs
  }

  # Reproject population raster to the NDVI grid once per grid definition
  grid_key <- paste(
    terra::crs(r, proj = TRUE),
    paste(round(terra::res(r), 12), collapse = ","),
    paste(round(c(terra::xmin(r), terra::xmax(r), terra::ymin(r), terra::ymax(r)), 8), collapse = ","),
    terra::ncol(r),
    terra::nrow(r),
    sep = "|"
  )
  if (is.na(pop_on_ndvi_key) || !identical(pop_on_ndvi_key, grid_key)) {
    pop_on_ndvi <- try(terra::project(pop_rast, r, method = "sum"), silent = TRUE)
    if (inherits(pop_on_ndvi, "try-error")) {
      if (!is.null(rr$tmp) && file.exists(rr$tmp)) unlink(rr$tmp)
      failed_dates <- c(failed_dates, date_i)
      next
    }
    # Zero-weight (not NA) for unpopulated cells: exact_extract's built-in
    # "weighted_mean" does NOT skip NA weights, so NA here would propagate NaN
    # for any polygon overlapping unpopulated cells (e.g. Bois de Vincennes/Boulogne).
    pop_on_ndvi[pop_on_ndvi < 0] <- 0
    pop_on_ndvi_key <- grid_key
  }

  vals_imu <- try(exactextractr::exact_extract(r, imu_extract, "mean", progress = FALSE), silent = TRUE)
  vals_arr_areaw <- try(exactextractr::exact_extract(r, arr_extract, "mean", progress = FALSE), silent = TRUE)
  vals_arr_popw <- try(
    exactextractr::exact_extract(
      r, arr_extract, "weighted_mean",
      weights = pop_on_ndvi,
      progress = FALSE
    ),
    silent = TRUE
  )

  if (inherits(vals_imu, "try-error") ||
      inherits(vals_arr_areaw, "try-error") ||
      inherits(vals_arr_popw, "try-error")) {
    if (!is.null(rr$tmp) && file.exists(rr$tmp)) unlink(rr$tmp)
    failed_dates <- c(failed_dates, date_i)
    next
  }

  vals_imu_dt <- data.table(
    imu_id = imu_prj$imu_id,
    date = as.Date(date_i),
    mean_ndvi_imu = as.numeric(vals_imu)
  )

  vals_arr_dt <- data.table(
    arr = as.character(arr_sf_prj$arr),
    date = as.Date(date_i),
    ndvi_areaw_arr_date = as.numeric(vals_arr_areaw),
    ndvi_popw_ghs_date = as.numeric(vals_arr_popw)
  )

  ndvi_imu_date <- rbind(ndvi_imu_date, vals_imu_dt, fill = TRUE)
  ndvi_imu_date <- unique(ndvi_imu_date, by = c("imu_id", "date"))
  fwrite(ndvi_imu_date, checkpoint_path)

  ndvi_arr_date <- rbind(ndvi_arr_date, vals_arr_dt, fill = TRUE)
  ndvi_arr_date <- unique(ndvi_arr_date, by = c("arr", "date"))
  fwrite(ndvi_arr_date, checkpoint_arr_path)

  if (!is.null(rr$tmp) && file.exists(rr$tmp)) unlink(rr$tmp)
}

# Keep failed dates log
if (length(failed_dates) > 0) {
  fwrite(data.table(date = unique(failed_dates)), failed_path)
}

if (nrow(ndvi_imu_date) == 0) {
  stop("No NDVI IMU-date rows available. Provide local clips or set download_missing <- TRUE.")
}
if (nrow(ndvi_arr_date) == 0) {
  stop("No NDVI arrondissement-date rows available.")
}

n_dates_used_imu <- uniqueN(ndvi_imu_date$date)
n_dates_used_arr <- uniqueN(ndvi_arr_date$date)
n_dates_used <- min(n_dates_used_imu, n_dates_used_arr)
if (n_dates_used < min_dates_required) {
  stop(
    "Only ", n_dates_used, " JJAS dates available across checkpoints; ",
    "minimum required is ", min_dates_required, ". ",
    "Enable download_missing <- TRUE or provide more local NDVI clips."
  )
}

setorder(ndvi_imu_date, imu_id, date)
setorder(ndvi_arr_date, arr, date)

# 6) year-balanced JJAS aggregation at IMU level
ndvi_imu_date[, year := as.integer(format(date, "%Y"))]

imu_year <- ndvi_imu_date[, .(
  ndvi_jjas_year = mean(mean_ndvi_imu, na.rm = TRUE),
  n_dates_year = sum(!is.na(mean_ndvi_imu))
), by = .(imu_id, year)]

imu_mean <- imu_year[, .(
  ndvi_jjas_2014_2017 = mean(ndvi_jjas_year, na.rm = TRUE),
  ndvi_n_years_jjas = sum(!is.na(ndvi_jjas_year)),
  ndvi_n_dates_jjas = sum(n_dates_year, na.rm = TRUE)
), by = imu_id]

# 7) year-balanced JJAS aggregation at arrondissement level 
ndvi_arr_date[, year := as.integer(format(date, "%Y"))]

arr_year <- ndvi_arr_date[, .(
  ndvi_popw_ghs_year = mean(ndvi_popw_ghs_date, na.rm = TRUE),
  ndvi_areaw_arr_year = mean(ndvi_areaw_arr_date, na.rm = TRUE),
  ndvi_n_dates_year_ghs = sum(!is.na(ndvi_popw_ghs_date))
), by = .(arr, year)]

arr_mean_ghs <- arr_year[, .(
  ndvi_popw_ghs = mean(ndvi_popw_ghs_year, na.rm = TRUE),
  ndvi_areaw_arr = mean(ndvi_areaw_arr_year, na.rm = TRUE),
  ndvi_n_years_jjas_ghs = sum(is.finite(ndvi_popw_ghs_year)),
  ndvi_n_dates_jjas_ghs = sum(ndvi_n_dates_year_ghs, na.rm = TRUE)
), by = arr]

# 8) arrondissement aggregation using IMU support weights (robustness)
imu_for_arr <- merge(
  imu_arr_key,
  imu_mean,
  by = "imu_id",
  all.x = TRUE
)

arr_ndvi <- imu_for_arr[, .(
  ndvi_popw_imu = safe_wmean(ndvi_jjas_2014_2017, pop_w),
  ndvi_areaw_imu = safe_wmean(ndvi_jjas_2014_2017, area_w),
  ndvi_popw_coverage = safe_share_nonmissing(ndvi_jjas_2014_2017, pop_w),
  ndvi_areaw_coverage = safe_share_nonmissing(ndvi_jjas_2014_2017, area_w),
  ndvi_n_years_jjas_imu = if (any(is.finite(ndvi_n_years_jjas))) {
    min(ndvi_n_years_jjas[is.finite(ndvi_n_years_jjas)], na.rm = TRUE)
  } else {
    NA_real_
  },
  ndvi_n_dates_jjas_imu = if (any(is.finite(ndvi_n_dates_jjas))) {
    min(ndvi_n_dates_jjas[is.finite(ndvi_n_dates_jjas)], na.rm = TRUE)
  } else {
    NA_real_
  }
), by = arr]

arr_ndvi <- merge(data.table(arr = paris_arr), arr_ndvi, by = "arr", all.x = TRUE)
arr_ndvi <- merge(arr_ndvi, arr_mean_ghs, by = "arr", all.x = TRUE)
setorder(arr_ndvi, arr)

# 9) save
fwrite(imu_mean, out_imu_path)
fwrite(ndvi_arr_date, out_arr_date_path)
fwrite(arr_ndvi, out_arr_path)

cat("\nSaved:\n")
cat(" -", checkpoint_path, "\n")
cat(" -", checkpoint_arr_path, "\n")
cat(" -", out_arr_date_path, "\n")
cat(" -", out_imu_path, "\n")
cat(" -", out_arr_path, "\n")
if (file.exists(failed_path)) cat(" -", failed_path, "\n")

cat("\nArrondissement NDVI range (GHS-pop + IMU-supported):\n")
print(arr_ndvi[, .(
  ndvi_popw_ghs_min = min(ndvi_popw_ghs, na.rm = TRUE),
  ndvi_popw_ghs_max = max(ndvi_popw_ghs, na.rm = TRUE),
  ndvi_areaw_arr_min = min(ndvi_areaw_arr, na.rm = TRUE),
  ndvi_areaw_arr_max = max(ndvi_areaw_arr, na.rm = TRUE),
  ndvi_n_years_jjas_ghs_min = min(ndvi_n_years_jjas_ghs, na.rm = TRUE),
  ndvi_n_years_jjas_ghs_max = max(ndvi_n_years_jjas_ghs, na.rm = TRUE),
  ndvi_n_dates_jjas_ghs_min = min(ndvi_n_dates_jjas_ghs, na.rm = TRUE),
  ndvi_n_dates_jjas_ghs_max = max(ndvi_n_dates_jjas_ghs, na.rm = TRUE),
  ndvi_popw_imu_min = min(ndvi_popw_imu, na.rm = TRUE),
  ndvi_popw_imu_max = max(ndvi_popw_imu, na.rm = TRUE),
  ndvi_areaw_imu_min = min(ndvi_areaw_imu, na.rm = TRUE),
  ndvi_areaw_imu_max = max(ndvi_areaw_imu, na.rm = TRUE),
  ndvi_n_years_jjas_imu_min = min(ndvi_n_years_jjas_imu, na.rm = TRUE),
  ndvi_n_years_jjas_imu_max = max(ndvi_n_years_jjas_imu, na.rm = TRUE),
  ndvi_n_dates_jjas_imu_min = min(ndvi_n_dates_jjas_imu, na.rm = TRUE),
  ndvi_n_dates_jjas_imu_max = max(ndvi_n_dates_jjas_imu, na.rm = TRUE)
)])
