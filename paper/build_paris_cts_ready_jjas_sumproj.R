# build_paris_cts_ready_jjas_sumproj.R: Script 5
# Merge deaths + heat panel + greenness into one CTS-ready JJAS dataset
#
# Final run configuration:
#   HEAT_SUPPORT=ghspopw_sum  RUN_TAG=_rerun_rule1_sum_ghsheat_ghsgreen
#
# Inputs:
#   - Heat panel built with GHS-POP sum-preserving reprojection
#     (paris_panel_residence_urbclim_daily_2008_2017_ghspopweighted_sumproj_JJAS.csv)
#   - Arrondissement greenness metrics (paris_arr_greenness_metrics.csv)
#
# This script merges the heat and greenness panels, creates CTS-ready
# time variables (district-year-month strata, day of week, day of season),
# and IQR-standardises all greenness metrics across the 20 arrondissements.
#
# MMT-centered heat variables are NOT created here; they depend on the
# baseline DLNM fit and are handled downstream.
#
# Outputs (tagged with RUN_TAG):
#   - paris_cts_ready_jjas_2008_2017{RUN_TAG}.csv
#   - paris_greenness_metric_summary{RUN_TAG}.csv
#   - paris_heat_metric_summary{RUN_TAG}.csv

library(data.table)
library(lubridate)

base_dir <- "/Users/armandeaboudrar-meda/Desktop/GVIestim/paper" # relative:   base_dir <- here::here()

norm_tag <- function(x) {
  x <- trimws(x)
  if (!nzchar(x)) return("")
  if (!startsWith(x, "_")) x <- paste0("_", x)
  x
}

tag_file <- function(stem, tag, ext = ".csv") {
  file.path(base_dir, paste0(stem, tag, ext))
}

# HEAT_SUPPORT:
#   - ghspopw_sum : GHS-POP-weighted panel built with sum-preserving reprojection (default)
#   - ghspopw  : GHS-POP-weighted panel with bilinear reprojection
#   - native   : native arrondissement-area heat panel
#   - popw     : IMU-population-weighted heat panel
heat_support <- tolower(trimws(Sys.getenv("HEAT_SUPPORT", "ghspopw_sum")))
if (heat_support %in% c("popw", "weighted", "imuweighted")) {
  panel_file <- "paris_panel_residence_urbclim_daily_2008_2017_imuweighted_JJAS.csv"
  default_run_tag <- ""
} else if (heat_support %in% c("native", "areaw", "area")) {
  panel_file <- "paris_panel_residence_urbclim_daily_2008_2017_JJAS.csv"
  default_run_tag <- "_nativeheat"
} else if (heat_support %in% c("ghspopw", "ghs", "ghsweighted", "ghs_popw")) {
  panel_file <- "paris_panel_residence_urbclim_daily_2008_2017_ghspopweighted_JJAS.csv"
  default_run_tag <- "_ghsheat"
} else if (heat_support %in% c("ghspopw_sum", "ghs_sum", "ghsweighted_sum")) {
  panel_file <- "paris_panel_residence_urbclim_daily_2008_2017_ghspopweighted_sumproj_JJAS.csv"
  default_run_tag <- "_ghsheat_sumproj"
} else {
  stop("Invalid HEAT_SUPPORT='", heat_support, "'. Use 'popw', 'native', 'ghspopw', or 'ghspopw_sum'.")
}

run_tag <- norm_tag(Sys.getenv("RUN_TAG", default_run_tag))

panel_path <- file.path(base_dir, panel_file)
green_path <- file.path(base_dir, "paris_arr_greenness_metrics.csv")

out_cts_path  <- tag_file("paris_cts_ready_jjas_2008_2017", run_tag)
out_green_sum <- tag_file("paris_greenness_metric_summary", run_tag)
out_heat_sum  <- tag_file("paris_heat_metric_summary", run_tag)

cat("\nCTS-ready build config:\n")
cat(" - HEAT_SUPPORT:", heat_support, "\n")
cat(" - RUN_TAG:", ifelse(nzchar(run_tag), run_tag, "<none>"), "\n")
cat(" - panel input:", panel_path, "\n")

# 1) load inputs
panel <- fread(panel_path)
green <- fread(green_path)

stopifnot(all(c("arr", "date", "deaths") %in% names(panel)))
stopifnot("arr" %in% names(green))

panel[, date := as.Date(date)]
panel[, arr := as.character(arr)]
green[, arr := as.character(arr)]

# restrict to JJAS just in case
panel <- panel[month(date) %in% 6:9]

# one row per arrondissement-day
dup_check <- panel[, .N, by = .(arr, date)][N > 1]
if (nrow(dup_check) > 0) {
  stop("Duplicate arrondissement-date rows found in panel.")
}

# one row per arrondissement in greenness table
dup_green <- green[, .N, by = arr][N > 1]
if (nrow(dup_green) > 0) {
  stop("Duplicate arrondissement rows found in greenness table.")
}

# 2) merge panel + greenness
dt <- merge(panel, green, by = "arr", all.x = TRUE)

# keep arrondissement ordering stable
paris_arr <- sprintf("751%02d", 1:20)
dt <- dt[arr %in% paris_arr]
setorder(dt, arr, date)

# 3) basic CTS-ready time variables
if (!("year" %in% names(dt)))  dt[, year := year(date)]
if (!("month" %in% names(dt))) dt[, month := month(date)]
if (!("dow" %in% names(dt)))   dt[, dow := wday(date, week_start = 1)]   # 1=Mon,...,7=Sun
if (!("ym" %in% names(dt)))    dt[, ym := format(date, "%Y-%m")]

dt[, arr_num := as.integer(sub("^751", "", arr))]
dt[, year_f := factor(year)]
dt[, month_f := factor(month)]
dt[, dow_f := factor(dow, levels = 1:7, labels = c("Mon","Tue","Wed","Thu","Fri","Sat","Sun"))]
dt[, arr_f := factor(arr, levels = paris_arr)]

# district-year-month stratum as in CTS design
dt[, stratum := paste(arr, year, sprintf("%02d", month), sep = "_")]
dt[, stratum_id := .GRP, by = stratum]

# day of summer season: June 1 = 1
dt[, season_start := as.Date(sprintf("%d-06-01", year))]
dt[, day_of_season := as.integer(date - season_start) + 1L]
dt[, season_start := NULL]

# useful helpers
dt[, weekend := as.integer(dow %in% c(6, 7))]
dt[, deaths_nonzero := as.integer(deaths > 0)]

# 4) define metric groups
heat_cols <- intersect(
  c("t2m_min", "t2m_mean", "t2m_max",
    "lst_min", "lst_mean", "lst_max",
    "wbgt_min", "wbgt_mean", "wbgt_max"),
  names(dt)
)

green_cols <- intersect(
  c("gvi_mean",
    "gvi_popw_points",
    "gvi_areaw_imu",
    "gvi_popw_imu",
    "ndvi_native_all",
    "ndvi_native_jjas",
    "ndvi_popw_ghs",
    "ndvi_areaw_arr",
    "ndvi_popw_imu",
    "ndvi_areaw_imu",
    "imu_high", "imu_low", "imu_veg_total",
    "imu_high_popw", "imu_low_popw", "imu_veg_total_popw",
    "imu_high_popw_ghs", "imu_low_popw_ghs", "imu_veg_total_popw_ghs"),
  names(dt)
)

if (length(heat_cols) == 0) {
  stop("No expected heat columns found.")
}
if (length(green_cols) == 0) {
  stop("No expected greenness columns found.")
}

# 5) IQR-standardise greenness metrics
# Standardisation:  G_std_iqr = (G - median(G)) / IQR(G)
# Computed over the N=20 arrondissement-level values (each arrondissement
# contributes one value because greenness is time-invariant).
# IQR standardisation (rather than SD) makes effect-modification
# coefficients comparable across greenness metrics that have very different
# natural scales (GVI 0-1, NDVI 0-1, IMU 0-100%) and is robust to
# the small-N, potentially non-normal distribution of district-level
# greenness. type=7 is R's default quantile algorithm (Hyndman & Fan
# type 7), used throughout for consistency.
std_iqr <- function(x) {
  med <- median(x, na.rm = TRUE)
  iq  <- IQR(x, na.rm = TRUE, type = 7)
  if (!is.finite(iq) || iq <= 0) return(rep(NA_real_, length(x)))
  (x - med) / iq
}

for (v in green_cols) {
  dt[, paste0(v, "_std_iqr") := std_iqr(get(v))]
}

# 6) metric summary tables
summarise_metric <- function(x) {
  q <- quantile(x, probs = c(0, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99, 1.00), na.rm = TRUE, type = 7)
  list(
    n = sum(!is.na(x)),
    min = unname(q[1]),
    p10 = unname(q[2]),
    q1 = unname(q[3]),
    median = unname(q[4]),
    q3 = unname(q[5]),
    p90 = unname(q[6]),
    p99 = unname(q[7]),
    max = unname(q[8]),
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    iqr = IQR(x, na.rm = TRUE, type = 7)
  )
}

green_summary <- rbindlist(lapply(green_cols, function(v) {
  out <- as.data.table(summarise_metric(dt[[v]]))
  out[, metric := v]
  out[, type := "greenness"]
  out
}), fill = TRUE)

setcolorder(green_summary, c("type", "metric", setdiff(names(green_summary), c("type", "metric"))))

heat_summary <- rbindlist(lapply(heat_cols, function(v) {
  out <- as.data.table(summarise_metric(dt[[v]]))
  out[, metric := v]
  out[, type := "heat"]
  out
}), fill = TRUE)

setcolorder(heat_summary, c("type", "metric", setdiff(names(heat_summary), c("type", "metric"))))

# 7) final checks
# all 20 arrondissements present?
n_arr <- uniqueN(dt[["arr"]])
if (n_arr != 20) {
  warning("Expected 20 arrondissements, found ", n_arr, ". Continuing with available arrondissements.")
}

# JJAS only?
stopifnot(all(dt[["month"]] %in% 6:9))

# day-of-season range should be 1..122 depending on leap year / date coverage
cat("day_of_season range:", min(dt$day_of_season), "to", max(dt$day_of_season), "\n")

# missing greenness check
green_missing <- dt[, lapply(.SD, function(x) mean(is.na(x))), .SDcols = green_cols]
cat("\nShare missing in greenness columns:\n")
print(green_missing)

# 8) save outputs
fwrite(dt, out_cts_path)
fwrite(green_summary, out_green_sum)
fwrite(heat_summary, out_heat_sum)

cat("\nSaved:\n")
cat(" -", out_cts_path, "\n")
cat(" -", out_green_sum, "\n")
cat(" -", out_heat_sum, "\n")

# 9) preview
cat("\nCTS-ready columns preview:\n")
print(names(dt))

cat("\nFirst rows:\n")
print(dt[1:10, .(
  arr, date, deaths,
  t2m_mean, t2m_max, lst_mean, lst_max, wbgt_mean, wbgt_max,
  gvi_mean, gvi_mean_std_iqr,
  ndvi_native_jjas, ndvi_native_jjas_std_iqr,
  imu_veg_total_popw, imu_veg_total_popw_std_iqr,
  stratum, stratum_id, dow, day_of_season
)])
