# estimate_paris_attributable_deaths_scenarios.R
#
# Run ONE attributable-burden scenario (one heat x one green metric).
# Used for both the main paper run (_rerun_rule1_sum_ghsheat_ghsgreen)
# and the native appendix comparison, controlled by RUN_TAG / CTS_TAG.
#
# For each selected heat x green pair, refits the effect-modification
# model and estimates observed, counterfactual, and avoided heat-
# attributable deaths under three greening scenarios:
#   1) max_all      — all arrondissements set to observed max greenness
#   2) p90_floor    — all below p90 raised to p90
#   3) plus1iqr_cap — each gets +1 IQR, capped at observed max
#
# Uses rule=1 in approx() for consistency with baseline and effectmod.
# Attribution threshold defaults to MMT.
#
# Inputs (tagged with CTS_TAG / BASELINE_TAG):
#   - paris_cts_ready_jjas_2008_2017{CTS_TAG}.csv
#   - paris_baseline_heat_metric_summary{BASELINE_TAG}.csv
#
# Outputs (tagged with RUN_TAG):
#   - paris_attr_scenario_rowlevel_<stub>{RUN_TAG}.csv
#   - paris_attr_scenario_by_arr_<stub>{RUN_TAG}.csv
#   - paris_attr_scenario_summary_<stub>{RUN_TAG}.csv

library(data.table)
library(lubridate)
library(splines)
library(gnm)
library(dlnm)

base_dir <- "/Users/armandeaboudrar-meda/Desktop/GVIestim/paper" # relative 

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
baseline_tag <- norm_tag(Sys.getenv("BASELINE_TAG", run_tag))

in_cts  <- tag_file("paris_cts_ready_jjas_2008_2017", cts_tag)
in_base <- tag_file("paris_baseline_heat_metric_summary", baseline_tag)

cat("\nAttributable single-run config:\n")
cat(" - CTS_TAG:", ifelse(nzchar(cts_tag), cts_tag, "<none>"), "\n")
cat(" - BASELINE_TAG:", ifelse(nzchar(baseline_tag), baseline_tag, "<none>"), "\n")
cat(" - RUN_TAG:", ifelse(nzchar(run_tag), run_tag, "<none>"), "\n")
cat(" - input CTS:", in_cts, "\n")
cat(" - input baseline:", in_base, "\n")

# 0) USER CHOICES
selected_heat <- Sys.getenv("SELECTED_HEAT", unset = "t2m_mean")
selected_green_std <- Sys.getenv("SELECTED_GREEN_STD", unset = "imu_veg_total_popw_ghs_std_iqr")
selected_green_raw <- sub("_std_iqr$", "", selected_green_std)

# Attribution threshold:
#   "mmt" = count heat-attributable deaths only above the fitted MMT
#   or set a numeric threshold manually, e.g. 25
heat_threshold_mode <- Sys.getenv("HEAT_THRESHOLD_MODE", unset = "mmt")

sanitize_token <- function(x) gsub("[^a-zA-Z0-9]+", "", tolower(x))
format_threshold_token <- function(threshold_mode) {
  if (is.character(threshold_mode) && length(threshold_mode) == 1 && threshold_mode == "mmt") {
    return("mmt")
  }
  thr_num <- suppressWarnings(as.numeric(threshold_mode))
  if (!is.finite(thr_num)) stop("heat_threshold_mode must be 'mmt' or numeric-like.")
  thr_chr <- format(round(thr_num, 3), trim = TRUE, scientific = FALSE)
  thr_chr <- gsub("\\.", "p", thr_chr)
  paste0("p", thr_chr)
}

run_stub <- paste(
  sanitize_token(selected_heat),
  sanitize_token(selected_green_raw),
  format_threshold_token(heat_threshold_mode),
  sep = "_"
)
out_stub <- if (nzchar(run_tag)) paste0(run_stub, run_tag) else run_stub
out_row  <- file.path(base_dir, paste0("paris_attr_scenario_rowlevel_", out_stub, ".csv"))
out_arr  <- file.path(base_dir, paste0("paris_attr_scenario_by_arr_", out_stub, ".csv"))
out_sum  <- file.path(base_dir, paste0("paris_attr_scenario_summary_", out_stub, ".csv"))

# 1) load inputs
dt <- fread(in_cts)
base <- fread(in_base)

dt[, date := as.Date(date)]

required_cols <- c("arr", "date", "deaths", "day_of_season", "stratum_id")
missing_req <- setdiff(required_cols, names(dt))
if (length(missing_req) > 0) {
  stop("Missing required columns in CTS-ready file: ", paste(missing_req, collapse = ", "))
}

if (!(selected_heat %in% names(dt))) {
  stop("Selected heat metric not found in CTS-ready file: ", selected_heat)
}
if (!(selected_green_std %in% names(dt))) {
  stop("Selected greenness metric not found in CTS-ready file: ", selected_green_std)
}
if (!(selected_green_raw %in% names(dt))) {
  stop("Selected raw greenness metric not found in CTS-ready file: ", selected_green_raw)
}

# Restrict to JJAS just in case
dt <- dt[month(date) %in% 6:9]

# User-set denominator for standardized impact metrics
city_population_ref <- 2200000
n_summers <- uniqueN(year(dt$date))

# Factors for gnm
dt[, stratum_f := factor(stratum_id)]

if ("dow_f" %in% names(dt)) {
  dt[, dow_f := factor(as.character(dow_f),
                       levels = c("Mon","Tue","Wed","Thu","Fri","Sat","Sun"))]
} else if ("dow" %in% names(dt)) {
  dt[, dow_f := factor(dow,
                       levels = 1:7,
                       labels = c("Mon","Tue","Wed","Thu","Fri","Sat","Sun"))]
} else {
  stop("Need either dow_f or dow in CTS-ready file.")
}

if ("year" %in% names(dt)) {
  dt[, year_f := factor(year)]
} else if ("year_f" %in% names(dt)) {
  dt[, year_f := factor(as.character(year_f))]
} else {
  stop("Need either year or year_f in CTS-ready file.")
}

brow <- base[metric == selected_heat]
if (nrow(brow) != 1) {
  stop("Could not find exactly one baseline row for heat metric: ", selected_heat)
}

heat_min <- brow$min[1]
heat_max <- brow$max[1]
heat_p90 <- brow$p90[1]
heat_p99 <- brow$p99[1]
heat_iqr <- brow$iqr[1]
heat_mmt <- brow$mmt[1]

cat("Selected heat metric: ", selected_heat, "\n", sep = "")
cat("Selected greenness metric: ", selected_green_std, "\n", sep = "")
cat("Raw greenness metric: ", selected_green_raw, "\n", sep = "")
cat("Heat MMT: ", round(heat_mmt, 3), "\n", sep = "")
cat("Heat p99: ", round(heat_p99, 3), "\n", sep = "")

# 2) helpers
resolve_block_names <- function(fit, object_name, mat_obj) {
  beta_names <- names(coef(fit))
  cn <- colnames(mat_obj)
  
  cand1 <- paste0(object_name, cn)
  if (all(cand1 %in% beta_names)) return(cand1)
  
  cand2 <- cn
  if (all(cand2 %in% beta_names)) return(cand2)
  
  cand3 <- beta_names[grepl(paste0("^", object_name), beta_names)]
  if (length(cand3) == ncol(mat_obj)) return(cand3)
  
  stop(
    "Could not match coefficient names for block ", object_name,
    ". Example fitted names: ", paste(head(beta_names, 20), collapse = ", ")
  )
}

make_effective_block <- function(beta_h, beta_g, V_hh, V_hg, V_gh, V_gg, g_value) {
  list(
    coef = beta_h + g_value * beta_g,
    vcov = V_hh + g_value * V_hg + g_value * V_gh + (g_value^2) * V_gg
  )
}

# 3) fit selected heat × green model
sub <- copy(dt[!is.na(get(selected_heat)) & !is.na(get(selected_green_std))])

x <- sub[[selected_heat]]
g <- sub[[selected_green_std]]

at_grid <- seq(
  floor(heat_min * 10) / 10,
  ceiling(heat_max * 10) / 10,
  by = 0.1
)

cbh <- crossbasis(
  x,
  lag = 1,
  argvar = list(fun = "ns", knots = heat_p90),
  arglag = list(fun = "integer")
)
colnames(cbh) <- paste0("h", seq_len(ncol(cbh)))

cbg <- cbh
cbg[] <- cbh[] * g
colnames(cbg) <- paste0("g", seq_len(ncol(cbg)))

fit <- gnm(
  deaths ~ cbh + cbg + dow_f + ns(day_of_season, df = 4):year_f,
  eliminate = stratum_f,
  family = quasipoisson(),
  data = sub,
  trace = FALSE
)

beta <- coef(fit)
V <- vcov(fit)

h_names <- resolve_block_names(fit, "cbh", cbh)
g_names <- resolve_block_names(fit, "cbg", cbg)

beta_h <- beta[h_names]
beta_g <- beta[g_names]

V_hh <- V[h_names, h_names, drop = FALSE]
V_gg <- V[g_names, g_names, drop = FALSE]
V_hg <- V[h_names, g_names, drop = FALSE]
V_gh <- V[g_names, h_names, drop = FALSE]

# 4) build arrondissement-level greenness scenarios
arr_green <- unique(sub[, .(
  arr,
  green_raw_obs = get(selected_green_raw),
  green_std_obs = get(selected_green_std)
)], by = "arr")

if (nrow(arr_green) != 20) {
  warning("Expected 20 arrondissement greenness rows; found ", nrow(arr_green))
}

green_median <- median(arr_green$green_raw_obs, na.rm = TRUE)
green_iqr <- IQR(arr_green$green_raw_obs, na.rm = TRUE, type = 7)
green_p90 <- as.numeric(quantile(arr_green$green_raw_obs, 0.90, na.rm = TRUE, type = 7))
green_max <- max(arr_green$green_raw_obs, na.rm = TRUE)

if (!is.finite(green_iqr) || green_iqr <= 0) {
  stop("Greenness IQR is not positive for raw metric: ", selected_green_raw)
}

to_std <- function(raw_value) {
  (raw_value - green_median) / green_iqr
}

scenario_long <- rbindlist(list(
  arr_green[, .(
    arr,
    scenario = "observed",
    green_raw_target = green_raw_obs,
    green_std_target = green_std_obs
  )],
  arr_green[, .(
    arr,
    scenario = "max_all",
    green_raw_target = green_max,
    green_std_target = to_std(green_max)
  )],
  arr_green[, .(
    arr,
    scenario = "p90_floor",
    green_raw_target = pmax(green_raw_obs, green_p90),
    green_std_target = to_std(pmax(green_raw_obs, green_p90))
  )],
  arr_green[, .(
    arr,
    scenario = "plus1iqr_cap",
    green_raw_target = pmin(green_raw_obs + green_iqr, green_max),
    green_std_target = to_std(pmin(green_raw_obs + green_iqr, green_max))
  )]
), use.names = TRUE)

# 5) RR prediction function for arbitrary greenness level
cp_cache <- new.env(parent = emptyenv())

get_cp_for_g <- function(g_value_std) {
  key <- format(round(g_value_std, 8), scientific = FALSE, trim = TRUE)
  
  if (exists(key, envir = cp_cache, inherits = FALSE)) {
    return(get(key, envir = cp_cache, inherits = FALSE))
  }
  
  eff <- make_effective_block(
    beta_h = beta_h,
    beta_g = beta_g,
    V_hh = V_hh,
    V_hg = V_hg,
    V_gh = V_gh,
    V_gg = V_gg,
    g_value = g_value_std
  )
  
  cp <- crosspred(
    cbh,
    coef = eff$coef,
    vcov = eff$vcov,
    model.link = "log",
    at = at_grid,
    cen = heat_mmt
  )
  
  assign(key, cp, envir = cp_cache)
  cp
}

predict_rr_vec <- function(temp_vec, g_std_vec) {
  out <- rep(NA_real_, length(temp_vec))
  u <- sort(unique(round(g_std_vec, 8)))
  
  for (gu in u) {
    idx <- which(round(g_std_vec, 8) == gu)
    cp <- get_cp_for_g(gu)
    out[idx] <- approx(cp$predvar, cp$allRRfit, xout = temp_vec[idx], rule = 1)$y
  }
  
  out
}

# 6) expand data to scenario panel and compute attributable deaths
row_dt <- merge(
  sub[, .(
    arr, date, deaths,
    heat_value = get(selected_heat),
    green_raw_obs = get(selected_green_raw),
    green_std_obs = get(selected_green_std)
  )],
  scenario_long,
  by = "arr",
  all.x = TRUE,
  allow.cartesian = TRUE
)

# predicted RR under observed and scenario greenness
row_dt[, rr_obs := predict_rr_vec(heat_value, green_std_obs)]
row_dt[, rr_scn := predict_rr_vec(heat_value, green_std_target)]

# Heat attribution threshold
if (is.character(heat_threshold_mode) && length(heat_threshold_mode) == 1 &&
    tolower(heat_threshold_mode) == "mmt") {
  threshold_value <- heat_mmt
} else {
  thr_num <- suppressWarnings(as.numeric(heat_threshold_mode))
  if (length(thr_num) == 1 && is.finite(thr_num)) {
    threshold_value <- thr_num
  } else {
    stop("heat_threshold_mode must be 'mmt' or a single numeric value.")
  }
}

row_dt[, above_threshold := heat_value > threshold_value]

# Heat-attributable fraction only above threshold
row_dt[, af_obs := fifelse(
  above_threshold,
  pmax((rr_obs - 1) / rr_obs, 0),
  0
)]

row_dt[, af_scn := fifelse(
  above_threshold,
  pmax((rr_scn - 1) / rr_scn, 0),
  0
)]

row_dt[, ad_obs := af_obs * deaths]
row_dt[, ad_scn := af_scn * deaths]
row_dt[, ad_avoided := ad_obs - ad_scn]

# 7) summaries
by_arr <- row_dt[, .(
  deaths_total = sum(deaths, na.rm = TRUE),
  deaths_hot_days = sum(deaths[above_threshold], na.rm = TRUE),
  ad_obs = sum(ad_obs, na.rm = TRUE),
  ad_scn = sum(ad_scn, na.rm = TRUE),
  ad_avoided = sum(ad_avoided, na.rm = TRUE),
  
  green_raw_obs = unique(green_raw_obs)[1],
  green_raw_target = unique(green_raw_target)[1]
), by = .(scenario, arr)]

setorder(by_arr, scenario, arr)

summary_dt <- row_dt[, .(
  n_rows = .N,
  n_hot_days_rows = sum(above_threshold, na.rm = TRUE),
  
  total_deaths = sum(deaths, na.rm = TRUE),
  total_deaths_hot_days = sum(deaths[above_threshold], na.rm = TRUE),
  
  ad_obs = sum(ad_obs, na.rm = TRUE),
  ad_scn = sum(ad_scn, na.rm = TRUE),
  ad_avoided = sum(ad_avoided, na.rm = TRUE),
  
  pct_reduction_vs_observed_attr = 100 * (sum(ad_avoided, na.rm = TRUE) / sum(ad_obs, na.rm = TRUE)),
  ad_obs_per100k_per_summer =
    1e5 * sum(ad_obs, na.rm = TRUE) / city_population_ref / n_summers,
  ad_scn_per100k_per_summer =
    1e5 * sum(ad_scn, na.rm = TRUE) / city_population_ref / n_summers,
  ad_avoided_per100k_per_summer =
    1e5 * sum(ad_avoided, na.rm = TRUE) / city_population_ref / n_summers
), by = scenario]

summary_dt[, `:=`(
  selected_heat = selected_heat,
  selected_green_std = selected_green_std,
  selected_green_raw = selected_green_raw,
  heat_mmt = heat_mmt,
  heat_p99 = heat_p99,
  heat_iqr = heat_iqr,
  attribution_threshold = threshold_value,
  green_median = green_median,
  green_iqr = green_iqr,
  green_p90 = green_p90,
  green_max = green_max,
  population_denominator = city_population_ref,
  n_summers = n_summers
)]

setcolorder(summary_dt, c(
  "selected_heat", "selected_green_std", "selected_green_raw",
  "scenario",
  "heat_mmt", "heat_p99", "heat_iqr", "attribution_threshold",
  "green_median", "green_iqr", "green_p90", "green_max",
  "population_denominator", "n_summers",
  setdiff(names(summary_dt), c(
    "selected_heat", "selected_green_std", "selected_green_raw",
    "scenario",
    "heat_mmt", "heat_p99", "heat_iqr", "attribution_threshold",
    "green_median", "green_iqr", "green_p90", "green_max",
    "population_denominator", "n_summers"
  ))
))

# 8) save
fwrite(row_dt, out_row)
fwrite(by_arr, out_arr)
fwrite(summary_dt, out_sum)

cat("\nSaved:\n")
cat(" -", out_row, "\n")
cat(" -", out_arr, "\n")
cat(" -", out_sum, "\n")

cat("\nScenario summary:\n")
print(summary_dt[, .(
  scenario,
  total_deaths,
  ad_obs,
  ad_scn,
  ad_avoided,
  pct_reduction_vs_observed_attr,
  ad_avoided_per100k_per_summer
)])

cat("\nTop arrondissement gains under max_all:\n")
print(
  by_arr[scenario == "max_all"][order(-ad_avoided)][1:10, .(
    arr,
    green_raw_obs,
    green_raw_target,
    ad_obs,
    ad_scn,
    ad_avoided
  )]
)

cat("\nDone.\n")
