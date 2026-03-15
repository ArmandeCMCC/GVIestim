# build_paris_panel_deaths_heat_urbclim.R
#
# APPENDIX comparison: native (area-average) heat panel.
# Extracts UrbClim daily heat metrics as simple arrondissement-level
# area means (no population weighting). This is the counterpart of
# build_paris_panel_deaths_heat_ghspopweighted_sumproj.R (main paper).
#
# Both panels cover the same 20 arrondissements, same 2008-2017 period,
# same 9 heat metrics; the only difference is weighting.
#
# Outputs:
#   - paris_panel_residence_urbclim_daily_2008_2017.csv
#   - paris_panel_residence_urbclim_daily_2008_2017_JJAS.csv

if (Sys.info()[["sysname"]] == "Darwin") {
  if (!nzchar(Sys.getenv("PROJ_LIB"))) Sys.setenv(PROJ_LIB = "/opt/homebrew/opt/proj/share/proj")
  if (!nzchar(Sys.getenv("GDAL_DATA"))) Sys.setenv(GDAL_DATA = "/opt/homebrew/opt/gdal/share/gdal")
}

library(data.table)
library(lubridate)
library(sf)
library(terra)
library(exactextractr)
library(ncdf4)

terra::gdalCache(4000)

base_dir <- "/Users/armandeaboudrar-meda/Desktop/GVIestim/paper"

deaths_full_path <- file.path(base_dir, "paris_daily_deaths_arr_residence_fullgrid_2008_2017.csv")
arr_shp_candidates <- c(
  file.path(base_dir, "arrondissements.shp"),
  file.path(base_dir, "NDVI_Paris_2014_2017", "_tmp_paris_shp", "arrondissements.shp")
)
arr_shp_path <- arr_shp_candidates[file.exists(arr_shp_candidates)][1]
if (is.na(arr_shp_path)) {
  stop("No local arrondissement shapefile found. Expected one of: ",
       paste(arr_shp_candidates, collapse = " ; "))
}

urbclim_base <- file.path(base_dir, "UrbClim")
metric_dirs <- list(
  LST  = file.path(urbclim_base, "LST"),
  T2M  = file.path(urbclim_base, "T2M"),
  WBGT = file.path(urbclim_base, "WBGT")
)

years <- 2008:2017
stats <- c("min","mean","max")  

# 1) read deaths panel (full 2008-2017 grid by arrondissement)
deaths <- fread(deaths_full_path)
stopifnot(all(c("arr","date","deaths") %in% names(deaths)))
deaths[, date := as.Date(date)]
deaths[, arr := as.character(arr)]
stopifnot(uniqueN(deaths$arr) == 20)

# 2) read arrondissement polygons and project to UrbClim CRS (EPSG:3035)
arr_sf <- st_read(arr_shp_path, quiet = TRUE)

if (!("c_arinsee" %in% names(arr_sf))) {
  stop("Local arrondissement shapefile must contain 'c_arinsee'. File: ", arr_shp_path)
}
arr_sf$arr <- as.character(arr_sf$c_arinsee)
paris_arr <- sprintf("751%02d", 1:20)
arr_sf <- arr_sf[arr_sf$arr %in% paris_arr, c("arr","geometry")]
stopifnot(nrow(arr_sf) == 20)

arr_sf_3035 <- st_transform(arr_sf, 3035)

# 3) Helper functions:
#    - parse date names from NetCDF layers
#    - validate NetCDF readability
#    - extract arrondissement mean per day for one file
parse_date_from_colname <- function(x) {
  s <- as.character(x)
  s <- sub("^[A-Za-z]+\\.", "", s)   # drop "mean."
  s <- sub("^X", "", s)
  s <- sub("\\.$", "", s)
  s <- sub("\\.[0-9]+$", "", s)
  s <- gsub("\\.", "-", s)
  s <- substr(s, 1, 10)
  as.Date(s, format = "%Y-%m-%d")
}

# check that it reads UrbClim files, and that they are complete, not corrupted
assert_nc_ok <- function(nc_path) {
  ok <- TRUE
  tryCatch({
    nc <- nc_open(nc_path)
    nc_close(nc)
  }, error = function(e) ok <<- FALSE)
  if (!ok) stop("NetCDF file is not readable (likely corrupted/partial download): ", nc_path)
}

extract_arr_daily_mean <- function(nc_path, value_name, arr_polygons_3035) {
  assert_nc_ok(nc_path) 
  
  r <- rast(nc_path)    
  # some UrbClim NetCDF files don't have CRS
  # UrbClim Paris grids are in LAEA Europe (EPSG:3035).
  if (!nzchar(crs(r))) {
    crs(r) <- "EPSG:3035"
  }
  # Kelvin -> Celsius if needed
  mx <- global(r, "max", na.rm = TRUE)[1,1]
  if (is.finite(mx) && mx > 200) r <- r - 273.15
  
  dts <- as.Date(terra::time(r))
  if (is.null(dts)) stop("No time dimension detected in: ", nc_path)
  names(r) <- format(dts, "%Y-%m-%d")
  
  vals <- exact_extract(r, arr_polygons_3035, "mean")
  dt <- as.data.table(vals)
  dt[, arr := arr_polygons_3035$arr]
  
  long <- melt(dt, id.vars = "arr", variable.name = "date_raw", value.name = value_name)
  long[, date := parse_date_from_colname(date_raw)]
  if (any(is.na(long$date))) {
    bad <- unique(long[is.na(date), date_raw])[1:20]
    stop("Failed to parse some date column names. Examples: ", paste(bad, collapse=", "))
  }
  long[, date_raw := NULL]
  setorder(long, arr, date)
  long
}

# 4) extract daily heat exposures for each metric/stat/year and merge to one
#    arrondissement-day exposure table.
expo_list <- list()

for (metric in names(metric_dirs)) {
  for (stt in stats) {
    colname <- paste0(tolower(metric), "_", stt)
    cat("\n=== Processing", metric, stt, "===\n")
    
    one_list <- list()
    for (yy in years) {
      nc_path <- file.path(metric_dirs[[metric]],
                           sprintf("%s_year_daily_%s_%d.nc", metric, stt, yy))
      if (!file.exists(nc_path)) stop("Missing file: ", nc_path)
      cat("  -", basename(nc_path), "\n")
      
      one_list[[as.character(yy)]] <- extract_arr_daily_mean(nc_path, colname, arr_sf_3035)
    }
    expo_list[[paste(metric, stt, sep="_")]] <- rbindlist(one_list)
  }
}

expo <- Reduce(function(a, b) merge(a, b, by=c("arr","date"), all=TRUE), expo_list)

# 5) merge deaths + exposures, add calendar fields, save full panel and JJAS.
panel <- merge(deaths, expo, by=c("arr","date"), all.x=TRUE)

expo_cols <- setdiff(names(expo), c("arr","date"))
panel <- panel[rowSums(is.na(panel[, ..expo_cols])) < length(expo_cols)]

panel[, dow := wday(date, week_start = 1)]
panel[, year := year(date)]
panel[, month := month(date)]
panel[, ym := format(date, "%Y-%m")]
setorder(panel, arr, date)

out_panel <- file.path(base_dir, "paris_panel_residence_urbclim_daily_2008_2017.csv")
out_panel_jjas <- file.path(base_dir, "paris_panel_residence_urbclim_daily_2008_2017_JJAS.csv")

fwrite(panel, out_panel)
fwrite(panel[month %in% 6:9], out_panel_jjas)

cat("\nSaved:\n -", out_panel, "\n -", out_panel_jjas, "\n")
