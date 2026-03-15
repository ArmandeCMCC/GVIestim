# estimate_paris_baseline_cts_heat_metrics_rule1.R: Script 6
#
# Fit baseline heat-only CTS/DLNM models to recover MMT, p99, and IQR
# for each of the 9 heat metrics (T2M/LST/WBGT x min/mean/max).
# Used for BOTH the main paper run (_rerun_rule1_sum_ghsheat_ghsgreen)
# and the native/unweighted appendix comparison — only the CTS input differs.
#
# Key methodological choice: approx(rule=1) returns NA outside the
# observed range instead of extrapolating (rule=2). This is the
# corrected version used for all final results.
#
# Model:
#   Conditional quasi-Poisson (gnm with eliminate) on arrondissement-
#   year-month strata, with day-of-week FE, ns(day_of_season, 4 df)
#   interacted with year, and a DLNM crossbasis (ns with 1 knot at
#   p90, integer lag 0-1).
#
# Outputs (tagged with RUN_TAG):
#   - paris_baseline_heat_metric_summary{RUN_TAG}.csv
#   - paris_baseline_heat_metric_curves{RUN_TAG}.csv

library(data.table)
library(lubridate)
library(splines)
library(gnm)
library(dlnm)

base_dir <- "/Users/armandeaboudrar-meda/Desktop/GVIestim/paper" # relative : base_dir <- here::here()
 
norm_tag <- function(x) {
  x <- trimws(x)
  if (!nzchar(x)) return("")
  if (!startsWith(x, "_")) x <- paste0("_", x)
  x
}

tag_file <- function(stem, tag, ext = ".csv") {
  file.path(base_dir, paste0(stem, tag, ext))
}

run_tag <- norm_tag(Sys.getenv("RUN_TAG", ""))
cts_tag <- norm_tag(Sys.getenv("CTS_TAG", run_tag))

in_path <- tag_file("paris_cts_ready_jjas_2008_2017", cts_tag)

out_summary <- tag_file("paris_baseline_heat_metric_summary", run_tag)
out_curves  <- tag_file("paris_baseline_heat_metric_curves", run_tag)

cat("\nBaseline run config:\n")
cat(" - CTS_TAG:", ifelse(nzchar(cts_tag), cts_tag, "<none>"), "\n")
cat(" - RUN_TAG:", ifelse(nzchar(run_tag), run_tag, "<none>"), "\n")
cat(" - input CTS:", in_path, "\n")

# 1) load data
dt <- fread(in_path)
dt[, date := as.Date(date)]

required_cols <- c("arr", "date", "deaths", "dow_f", "year_f", "day_of_season", "stratum_id")
missing_req <- setdiff(required_cols, names(dt))
if (length(missing_req) > 0) {
  stop("Missing required columns in CTS-ready file: ", paste(missing_req, collapse = ", "))
}

heat_metrics <- intersect(
  c("t2m_min", "t2m_mean", "t2m_max",
    "lst_min", "lst_mean", "lst_max",
    "wbgt_min", "wbgt_mean", "wbgt_max"),
  names(dt)
)

if (length(heat_metrics) == 0) {
  stop("No heat metric columns found.")
}

# restrict to JJAS just in case
dt <- dt[month(date) %in% 6:9]

# factors needed by gnm
dt[, stratum_f := factor(stratum_id)]
dt[, dow_f := factor(dow_f)]
dt[, year_f := factor(year_f)]

# 2) helpers
# rule=1: return NA for temperatures outside the observed range rather than
# extrapolating (rule=2). This is the conservative choice — the dose-response
# curve should not be evaluated beyond the support of the data.
pred_extract <- function(cp, x) {
  rr  <- approx(cp$predvar, cp$allRRfit,  xout = x, rule = 1)$y
  low <- approx(cp$predvar, cp$allRRlow,  xout = x, rule = 1)$y
  hi  <- approx(cp$predvar, cp$allRRhigh, xout = x, rule = 1)$y
  list(rr = rr, low = low, high = hi)
}

# 3) loop over heat metrics
summary_list <- list()
curve_list <- list()

for (metric in heat_metrics) {
  cat("\n============================================================\n")
  cat("Fitting baseline CTS/DLNM for:", metric, "\n")
  cat("============================================================\n")
  
  sub <- copy(dt[!is.na(get(metric))])
  
  if (nrow(sub) == 0) {
    warning("Skipping ", metric, ": all values are NA.")
    next
  }
  
  x <- sub[[metric]]
  
  metric_min <- min(x, na.rm = TRUE)
  metric_max <- max(x, na.rm = TRUE)
  # type=7 is R's default quantile algorithm (Hyndman & Fan type 7),
  # used throughout for consistency across all scripts.
  metric_p90 <- as.numeric(quantile(x, 0.90, na.rm = TRUE, type = 7))
  metric_p99 <- as.numeric(quantile(x, 0.99, na.rm = TRUE, type = 7))
  metric_iqr <- IQR(x, na.rm = TRUE, type = 7)
  metric_med <- median(x, na.rm = TRUE)
  
  cat("Range:", round(metric_min, 3), "to", round(metric_max, 3), "\n")
  cat("p90:", round(metric_p90, 3), " p99:", round(metric_p99, 3), " IQR:", round(metric_iqr, 3), "\n")
  
  # prediction grid
  at_grid <- seq(
    floor(metric_min * 10) / 10,
    ceiling(metric_max * 10) / 10,
    by = 0.1
  )
  
  # DLNM crossbasis
  # - Exposure-response: natural cubic B-spline with 1 internal knot at
  #   p90 of the summer temperature distribution, following Paris Paper.
  #   The knot at p90 concentrates
  #   flexibility in the upper tail where the heat-mortality relationship
  #   steepens, while keeping the curve smooth elsewhere.
  # - Lag-response: integer function over lag 0-1 (same-day + next-day),
  #   capturing the effect of heat on mortality. Short lags are
  #   standard in summer heat-mortality studies because heat effects on
  #   all-cause mortality are immediate 
  cb <- crossbasis(
    x,
    lag = 1,
    argvar = list(fun = "ns", knots = metric_p90),
    arglag = list(fun = "integer")
  )
  
  # Baseline CTS model following the case time series (CTS) design
  # - quasipoisson(): accounts for overdispersion in daily death counts.
  # - eliminate = stratum_f: conditions out district-year-month intercepts
  #   via exact Poisson conditioning, controlling all time-invariant
  #   confounders within each stratum (equivalent to conditional logistic
  #   regression but for count data). This is central to the CTS design.
  # - ns(day_of_season, df=4):year_f: smooth intra-seasonal trend with
  #   4 df, interacted with year to allow the seasonal shape to vary
  #   across years (captures e.g. varying summer onset).
  # - dow_f: day-of-week fixed effects for within-week mortality patterns.
  fit <- gnm(
    deaths ~ cb + dow_f + ns(day_of_season, df = 4):year_f,
    eliminate = stratum_f,
    family = quasipoisson(),
    data = sub,
    trace = FALSE
  )
  
  # first prediction just to locate MMT
  cp0 <- crosspred(
    cb,
    fit,
    at = at_grid,
    cen = metric_med
  )
  
  mmt_idx <- which.min(cp0$allRRfit)
  mmt <- cp0$predvar[mmt_idx]
  
  cat("MMT:", round(mmt, 3), "\n")
  
  # re-center predictions at MMT
  cp <- crosspred(
    cb,
    fit,
    at = at_grid,
    cen = mmt
  )
  
  # comparable evaluation points
  target_1iqr <- mmt + metric_iqr
  target_2iqr <- mmt + 2 * metric_iqr
  
  p99_res  <- pred_extract(cp, metric_p99)
  iqr1_res <- pred_extract(cp, target_1iqr)
  iqr2_res <- pred_extract(cp, target_2iqr)
  
  # save curve
  curve_dt <- data.table(
    metric = metric,
    temp = cp$predvar,
    rr = cp$allRRfit,
    rr_low = cp$allRRlow,
    rr_high = cp$allRRhigh,
    mmt = mmt,
    p90 = metric_p90,
    p99 = metric_p99,
    iqr = metric_iqr
  )
  
  curve_list[[metric]] <- curve_dt
  
  # save summary
  summary_dt <- data.table(
    metric = metric,
    n_obs = nrow(sub),
    total_deaths = sum(sub$deaths, na.rm = TRUE),
    
    min = metric_min,
    p90 = metric_p90,
    p99 = metric_p99,
    max = metric_max,
    median = metric_med,
    iqr = metric_iqr,
    mmt = mmt,
    
    target_1iqr = target_1iqr,
    target_2iqr = target_2iqr,
    
    rr_p99 = p99_res$rr,
    rr_p99_low = p99_res$low,
    rr_p99_high = p99_res$high,
    cr_p99 = (p99_res$rr - 1) * 100,
    cr_p99_low = (p99_res$low - 1) * 100,
    cr_p99_high = (p99_res$high - 1) * 100,
    
    rr_1iqr = iqr1_res$rr,
    rr_1iqr_low = iqr1_res$low,
    rr_1iqr_high = iqr1_res$high,
    cr_1iqr = (iqr1_res$rr - 1) * 100,
    cr_1iqr_low = (iqr1_res$low - 1) * 100,
    cr_1iqr_high = (iqr1_res$high - 1) * 100,
    
    rr_2iqr = iqr2_res$rr,
    rr_2iqr_low = iqr2_res$low,
    rr_2iqr_high = iqr2_res$high,
    cr_2iqr = (iqr2_res$rr - 1) * 100,
    cr_2iqr_low = (iqr2_res$low - 1) * 100,
    cr_2iqr_high = (iqr2_res$high - 1) * 100
  )
  
  summary_list[[metric]] <- summary_dt
}

# 4) bind and save
summary_out <- rbindlist(summary_list, fill = TRUE)
curves_out  <- rbindlist(curve_list, fill = TRUE)

setorder(summary_out, metric)
setorder(curves_out, metric, temp)

fwrite(summary_out, out_summary)
fwrite(curves_out, out_curves)

cat("\nSaved:\n")
cat(" -", out_summary, "\n")
cat(" -", out_curves, "\n")

# 5) preview
cat("\nDone.\n")
