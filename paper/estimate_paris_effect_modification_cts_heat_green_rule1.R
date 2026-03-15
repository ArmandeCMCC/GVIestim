# estimate_paris_effect_modification_cts_heat_green_rule1.R
#
# Fit heat x greenness interaction CTS/DLNM models for all 9 heat metrics
# crossed with every available IQR-standardised greenness metric.
# Used for BOTH the main paper run (_rerun_rule1_sum_ghsheat_ghsgreen)
# and the native appendix comparison — only the CTS and baseline inputs differ.
#
# Key methodological choices:
#   - approx(rule=1): NA outside observed range (corrected version)
#   - Linear interaction: cbg = cbh * greenness_std_iqr
#   - Primary contrast: p10 -> p90 in standardised greenness at heat p99
#
# The script fits ALL available green metrics in the CTS file. Downstream
# summarize_paris_effectmod_results.R separates them into:
#   - GHS-POP weighted set (main paper): gvi_popw_points, ndvi_popw_ghs,
#     imu_veg_total_popw_ghs, imu_high_popw_ghs, imu_low_popw_ghs
#   - Native/unweighted set (appendix): gvi_mean, ndvi_native_jjas, imu_veg_total
#
# Model:
#   Conditional quasi-Poisson (gnm with eliminate) on arrondissement-
#   year-month strata, with day-of-week FE, ns(day_of_season, 4 df)
#   interacted with year, DLNM crossbasis for heat (ns, 1 knot at p90,
#   integer lag 0-1), and linear interaction with greenness.
#
# Outputs (tagged with RUN_TAG):
#   - paris_effectmod_heat_green_summary{RUN_TAG}.csv
#   - paris_effectmod_heat_green_curves{RUN_TAG}.csv

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

in_path <- tag_file("paris_cts_ready_jjas_2008_2017", cts_tag)
baseline_path <- tag_file("paris_baseline_heat_metric_summary", baseline_tag)

out_summary <- tag_file("paris_effectmod_heat_green_summary", run_tag)
out_curves  <- tag_file("paris_effectmod_heat_green_curves", run_tag)

cat("\nEffect-modification run config:\n")
cat(" - CTS_TAG:", ifelse(nzchar(cts_tag), cts_tag, "<none>"), "\n")
cat(" - BASELINE_TAG:", ifelse(nzchar(baseline_tag), baseline_tag, "<none>"), "\n")
cat(" - RUN_TAG:", ifelse(nzchar(run_tag), run_tag, "<none>"), "\n")
cat(" - input CTS:", in_path, "\n")
cat(" - input baseline:", baseline_path, "\n")

# 1) load data
dt <- fread(in_path)
base <- fread(baseline_path)

dt[, date := as.Date(date)]

required_cols <- c("arr", "date", "deaths", "day_of_season", "stratum_id")
missing_req <- setdiff(required_cols, names(dt))
if (length(missing_req) > 0) {
  stop("Missing required columns in CTS-ready file: ", paste(missing_req, collapse = ", "))
}

# Restrict to JJAS just in case
dt <- dt[month(date) %in% 6:9]

# Factors for gnm
dt[, stratum_f := factor(stratum_id)]

# Rebuild clean factors after fread()
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

heat_metrics <- intersect(
  c("t2m_min", "t2m_mean", "t2m_max",
    "lst_min", "lst_mean", "lst_max",
    "wbgt_min", "wbgt_mean", "wbgt_max"),
  names(dt)
)

green_metrics <- intersect(
  c("gvi_mean_std_iqr",
    "gvi_popw_points_std_iqr",
    "ndvi_native_jjas_std_iqr",
    "ndvi_native_all_std_iqr",
    "ndvi_popw_ghs_std_iqr",
    "ndvi_areaw_arr_std_iqr",
    "imu_veg_total_std_iqr",
    "imu_high_std_iqr",
    "imu_low_std_iqr",
    "ndvi_popw_imu_std_iqr",
    "ndvi_areaw_imu_std_iqr",
    "imu_veg_total_popw_std_iqr",
    "imu_high_popw_std_iqr",
    "imu_low_popw_std_iqr",
    "imu_veg_total_popw_ghs_std_iqr",
    "imu_high_popw_ghs_std_iqr",
    "imu_low_popw_ghs_std_iqr",
    "gvi_areaw_imu_std_iqr",
    "gvi_popw_imu_std_iqr"),
  names(dt)
)

if (length(heat_metrics) == 0) stop("No heat metrics found.")
if (length(green_metrics) == 0) stop("No greenness metrics found.")

base <- base[metric %in% heat_metrics]
setkey(base, metric)

# 2) helpers
pred_extract <- function(cp, x) {
  rr  <- approx(cp$predvar, cp$allRRfit,  xout = x, rule = 1)$y
  low <- approx(cp$predvar, cp$allRRlow,  xout = x, rule = 1)$y
  hi  <- approx(cp$predvar, cp$allRRhigh, xout = x, rule = 1)$y
  list(rr = rr, low = low, high = hi)
}

# interaction-only contrast:
# if cp_int is the RR curve for a +1 unit increase in standardised green,
# then moving from -1 to +1 corresponds to delta_g = 2
contrast_from_cp_int <- function(cp_int, x, delta_g = 2) {
  one <- pred_extract(cp_int, x)
  
  log_rr1 <- log(one$rr)
  se_log1 <- (log(one$high) - log(one$low)) / (2 * 1.96)
  
  log_rr_d <- delta_g * log_rr1
  se_log_d <- abs(delta_g) * se_log1
  
  rr_d   <- exp(log_rr_d)
  low_d  <- exp(log_rr_d - 1.96 * se_log_d)
  high_d <- exp(log_rr_d + 1.96 * se_log_d)
  
  list(
    rr = rr_d,
    low = low_d,
    high = high_d,
    cr = (rr_d - 1) * 100,
    cr_low = (low_d - 1) * 100,
    cr_high = (high_d - 1) * 100
  )
}

# robustly recover coefficient names for matrix terms inside gnm
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

# effective coefficients / vcov for a given greenness level g
make_effective_block <- function(beta_h, beta_g, V_hh, V_hg, V_gh, V_gg, g_value) {
  list(
    coef = beta_h + g_value * beta_g,
    vcov = V_hh + g_value * V_hg + g_value * V_gh + (g_value^2) * V_gg
  )
}

# 3) loop over heat × greenness models
summary_list <- list()
curve_list <- list()
model_counter <- 0L

# start clean output files
if (file.exists(out_summary)) file.remove(out_summary)
if (file.exists(out_curves)) file.remove(out_curves)

for (h in heat_metrics) {
  
  brow <- base[.(h)]
  if (nrow(brow) != 1) {
    warning("Skipping ", h, ": baseline summary row not found or duplicated.")
    next
  }
  
  h_min <- brow$min[1]
  h_max <- brow$max[1]
  h_p90 <- brow$p90[1]
  h_p99 <- brow$p99[1]
  h_iqr <- brow$iqr[1]
  h_mmt <- brow$mmt[1]
  
  at_grid <- seq(
    floor(h_min * 10) / 10,
    ceiling(h_max * 10) / 10,
    by = 0.1
  )
  
  target_1iqr <- h_mmt + h_iqr
  target_2iqr <- h_mmt + 2 * h_iqr
  
  for (gstd in green_metrics) {
    
    graw <- sub("_std_iqr$", "", gstd)
    
    cat("\n============================================================\n")
    cat("Heat:", h, " | Green:", gstd, "\n")
    cat("============================================================\n")
    
    sub <- copy(dt[!is.na(get(h)) & !is.na(get(gstd))])
    
    if (nrow(sub) == 0) {
      warning("Skipping ", h, " x ", gstd, ": no usable rows.")
      next
    }
    
    # Raw greenness summary across arrondissements
    if (graw %in% names(sub)) {
      g_arr <- unique(sub[, .(arr, g_raw = get(graw), g_std = get(gstd))], by = "arr")
      g_raw_median <- median(g_arr$g_raw, na.rm = TRUE)
      g_raw_iqr    <- IQR(g_arr$g_raw, na.rm = TRUE, type = 7)
      g_raw_low    <- g_raw_median - g_raw_iqr
      g_raw_high   <- g_raw_median + g_raw_iqr
      g_raw_p10    <- as.numeric(quantile(g_arr$g_raw, 0.10, na.rm = TRUE, type = 7))
      g_raw_p90    <- as.numeric(quantile(g_arr$g_raw, 0.90, na.rm = TRUE, type = 7))

      if (is.finite(g_raw_iqr) && g_raw_iqr > 0) {
        g_std_p10 <- (g_raw_p10 - g_raw_median) / g_raw_iqr
        g_std_p90 <- (g_raw_p90 - g_raw_median) / g_raw_iqr
      } else {
        g_std_p10 <- as.numeric(quantile(g_arr$g_std, 0.10, na.rm = TRUE, type = 7))
        g_std_p90 <- as.numeric(quantile(g_arr$g_std, 0.90, na.rm = TRUE, type = 7))
      }
    } else {
      g_raw_median <- NA_real_
      g_raw_iqr    <- NA_real_
      g_raw_low    <- NA_real_
      g_raw_high   <- NA_real_
      g_raw_p10    <- NA_real_
      g_raw_p90    <- NA_real_
      g_std_p10    <- NA_real_
      g_std_p90    <- NA_real_
    }

    # Fallback if standardized p10/p90 is not finite
    if (!is.finite(g_std_p10)) g_std_p10 <- -1
    if (!is.finite(g_std_p90)) g_std_p90 <- 1
    delta_std_p10p90 <- g_std_p90 - g_std_p10
    
    x <- sub[[h]]
    g <- sub[[gstd]]
    
    # Heat crossbasis
    cbh <- crossbasis(
      x,
      lag = 1,
      argvar = list(fun = "ns", knots = h_p90),
      arglag = list(fun = "integer")
    )
    colnames(cbh) <- paste0("h", seq_len(ncol(cbh)))
    
    # Interaction block: same crossbasis columns multiplied by green
    cbg <- cbh
    cbg[] <- cbh[] * g
    colnames(cbg) <- paste0("g", seq_len(ncol(cbg)))
    
    fit <- tryCatch(
      gnm(
        deaths ~ cbh + cbg + dow_f + ns(day_of_season, df = 4):year_f,
        eliminate = stratum_f,
        family = quasipoisson(),
        data = sub,
        trace = FALSE
      ),
      error = function(e) {
        warning("Model failed for ", h, " x ", gstd, ": ", e$message)
        NULL
      }
    )
    
    if (is.null(fit)) next
    
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
    
    # Curves for low/high greenness under two contrasts:
    #  - sensitivity: -1/+1 IQR
    #  - primary visual contrast: p10/p90 in standardized greenness
    eff_low_iqr  <- make_effective_block(beta_h, beta_g, V_hh, V_hg, V_gh, V_gg, g_value = -1)
    eff_high_iqr <- make_effective_block(beta_h, beta_g, V_hh, V_hg, V_gh, V_gg, g_value =  1)
    eff_p10      <- make_effective_block(beta_h, beta_g, V_hh, V_hg, V_gh, V_gg, g_value = g_std_p10)
    eff_p90      <- make_effective_block(beta_h, beta_g, V_hh, V_hg, V_gh, V_gg, g_value = g_std_p90)

    cp_low_iqr <- crosspred(
      cbh,
      coef = eff_low_iqr$coef,
      vcov = eff_low_iqr$vcov,
      model.link = "log",
      at = at_grid,
      cen = h_mmt
    )

    cp_high_iqr <- crosspred(
      cbh,
      coef = eff_high_iqr$coef,
      vcov = eff_high_iqr$vcov,
      model.link = "log",
      at = at_grid,
      cen = h_mmt
    )

    cp_p10 <- crosspred(
      cbh,
      coef = eff_p10$coef,
      vcov = eff_p10$vcov,
      model.link = "log",
      at = at_grid,
      cen = h_mmt
    )

    cp_p90 <- crosspred(
      cbh,
      coef = eff_p90$coef,
      vcov = eff_p90$vcov,
      model.link = "log",
      at = at_grid,
      cen = h_mmt
    )
    
    # Interaction-only curve: RR per +1 standardized-green unit
    cp_int <- crosspred(
      cbh,
      coef = beta_g,
      vcov = V_gg,
      model.link = "log",
      at = at_grid,
      cen = h_mmt
    )
    
    # Predictions at target temperatures
    low_p99   <- pred_extract(cp_low_iqr,  h_p99)
    high_p99  <- pred_extract(cp_high_iqr, h_p99)
    mod_p99   <- contrast_from_cp_int(cp_int, h_p99, delta_g = 2)

    low_1iqr  <- pred_extract(cp_low_iqr,  target_1iqr)
    high_1iqr <- pred_extract(cp_high_iqr, target_1iqr)
    mod_1iqr  <- contrast_from_cp_int(cp_int, target_1iqr, delta_g = 2)

    low_2iqr  <- pred_extract(cp_low_iqr,  target_2iqr)
    high_2iqr <- pred_extract(cp_high_iqr, target_2iqr)
    mod_2iqr  <- contrast_from_cp_int(cp_int, target_2iqr, delta_g = 2)

    # Additional support-safe p10->p90 contrasts (for traceability)
    mod_p99_p10p90  <- contrast_from_cp_int(cp_int, h_p99,      delta_g = delta_std_p10p90)
    mod_1iqr_p10p90 <- contrast_from_cp_int(cp_int, target_1iqr, delta_g = delta_std_p10p90)
    mod_2iqr_p10p90 <- contrast_from_cp_int(cp_int, target_2iqr, delta_g = delta_std_p10p90)
    
    # Summary row
    srow <- data.table(
      heat_metric = h,
      green_metric_std = gstd,
      green_metric_raw = graw,
      
      n_obs = nrow(sub),
      total_deaths = sum(sub$deaths, na.rm = TRUE),
      
      heat_mmt = h_mmt,
      heat_p90 = h_p90,
      heat_p99 = h_p99,
      heat_iqr = h_iqr,
      heat_target_1iqr = target_1iqr,
      heat_target_2iqr = target_2iqr,
      
      green_low_std = -1,
      green_high_std = 1,
      green_raw_median = g_raw_median,
      green_raw_iqr = g_raw_iqr,
      green_raw_low = g_raw_low,
      green_raw_high = g_raw_high,
      green_raw_p10 = g_raw_p10,
      green_raw_p90 = g_raw_p90,
      green_std_p10 = g_std_p10,
      green_std_p90 = g_std_p90,
      delta_std_p10p90 = delta_std_p10p90,
      
      rr_low_p99 = low_p99$rr,
      rr_low_p99_low = low_p99$low,
      rr_low_p99_high = low_p99$high,
      rr_high_p99 = high_p99$rr,
      rr_high_p99_low = high_p99$low,
      rr_high_p99_high = high_p99$high,
      rr_ratio_p99 = mod_p99$rr,
      rr_ratio_p99_low = mod_p99$low,
      rr_ratio_p99_high = mod_p99$high,
      cr_mod_p99 = mod_p99$cr,
      cr_mod_p99_low = mod_p99$cr_low,
      cr_mod_p99_high = mod_p99$cr_high,
      
      rr_low_1iqr = low_1iqr$rr,
      rr_low_1iqr_low = low_1iqr$low,
      rr_low_1iqr_high = low_1iqr$high,
      rr_high_1iqr = high_1iqr$rr,
      rr_high_1iqr_low = high_1iqr$low,
      rr_high_1iqr_high = high_1iqr$high,
      rr_ratio_1iqr = mod_1iqr$rr,
      rr_ratio_1iqr_low = mod_1iqr$low,
      rr_ratio_1iqr_high = mod_1iqr$high,
      cr_mod_1iqr = mod_1iqr$cr,
      cr_mod_1iqr_low = mod_1iqr$cr_low,
      cr_mod_1iqr_high = mod_1iqr$cr_high,
      
      rr_low_2iqr = low_2iqr$rr,
      rr_low_2iqr_low = low_2iqr$low,
      rr_low_2iqr_high = low_2iqr$high,
      rr_high_2iqr = high_2iqr$rr,
      rr_high_2iqr_low = high_2iqr$low,
      rr_high_2iqr_high = high_2iqr$high,
      rr_ratio_2iqr = mod_2iqr$rr,
      rr_ratio_2iqr_low = mod_2iqr$low,
      rr_ratio_2iqr_high = mod_2iqr$high,
      cr_mod_2iqr = mod_2iqr$cr,
      cr_mod_2iqr_low = mod_2iqr$cr_low,
      cr_mod_2iqr_high = mod_2iqr$cr_high,

      # support-safe p10->p90 contrasts
      rr_ratio_p99_p10p90 = mod_p99_p10p90$rr,
      rr_ratio_p99_p10p90_low = mod_p99_p10p90$low,
      rr_ratio_p99_p10p90_high = mod_p99_p10p90$high,
      cr_mod_p99_p10p90 = mod_p99_p10p90$cr,
      cr_mod_p99_p10p90_low = mod_p99_p10p90$cr_low,
      cr_mod_p99_p10p90_high = mod_p99_p10p90$cr_high,

      rr_ratio_1iqr_p10p90 = mod_1iqr_p10p90$rr,
      rr_ratio_1iqr_p10p90_low = mod_1iqr_p10p90$low,
      rr_ratio_1iqr_p10p90_high = mod_1iqr_p10p90$high,
      cr_mod_1iqr_p10p90 = mod_1iqr_p10p90$cr,
      cr_mod_1iqr_p10p90_low = mod_1iqr_p10p90$cr_low,
      cr_mod_1iqr_p10p90_high = mod_1iqr_p10p90$cr_high,

      rr_ratio_2iqr_p10p90 = mod_2iqr_p10p90$rr,
      rr_ratio_2iqr_p10p90_low = mod_2iqr_p10p90$low,
      rr_ratio_2iqr_p10p90_high = mod_2iqr_p10p90$high,
      cr_mod_2iqr_p10p90 = mod_2iqr_p10p90$cr,
      cr_mod_2iqr_p10p90_low = mod_2iqr_p10p90$cr_low,
      cr_mod_2iqr_p10p90_high = mod_2iqr_p10p90$cr_high
    )
    
    model_counter <- model_counter + 1L
    summary_list[[model_counter]] <- srow
    
    # Curves for plotting
    curve_low_iqr <- data.table(
      heat_metric = h,
      green_metric_std = gstd,
      green_metric_raw = graw,
      green_level = "low_minus_1iqr",
      green_value_std = -1,
      temp = cp_low_iqr$predvar,
      rr = cp_low_iqr$allRRfit,
      rr_low = cp_low_iqr$allRRlow,
      rr_high = cp_low_iqr$allRRhigh,
      heat_mmt = h_mmt,
      heat_p99 = h_p99,
      heat_iqr = h_iqr
    )

    curve_high_iqr <- data.table(
      heat_metric = h,
      green_metric_std = gstd,
      green_metric_raw = graw,
      green_level = "high_plus_1iqr",
      green_value_std = 1,
      temp = cp_high_iqr$predvar,
      rr = cp_high_iqr$allRRfit,
      rr_low = cp_high_iqr$allRRlow,
      rr_high = cp_high_iqr$allRRhigh,
      heat_mmt = h_mmt,
      heat_p99 = h_p99,
      heat_iqr = h_iqr
    )

    curve_p10 <- data.table(
      heat_metric = h,
      green_metric_std = gstd,
      green_metric_raw = graw,
      green_level = "low_p10",
      green_value_std = g_std_p10,
      temp = cp_p10$predvar,
      rr = cp_p10$allRRfit,
      rr_low = cp_p10$allRRlow,
      rr_high = cp_p10$allRRhigh,
      heat_mmt = h_mmt,
      heat_p99 = h_p99,
      heat_iqr = h_iqr
    )

    curve_p90 <- data.table(
      heat_metric = h,
      green_metric_std = gstd,
      green_metric_raw = graw,
      green_level = "high_p90",
      green_value_std = g_std_p90,
      temp = cp_p90$predvar,
      rr = cp_p90$allRRfit,
      rr_low = cp_p90$allRRlow,
      rr_high = cp_p90$allRRhigh,
      heat_mmt = h_mmt,
      heat_p99 = h_p99,
      heat_iqr = h_iqr
    )

    curve_list[[model_counter]] <- rbind(
      curve_low_iqr,
      curve_high_iqr,
      curve_p10,
      curve_p90
    )
    
    # Incremental save so progress is not lost if R crashes
    summary_sofar <- rbindlist(summary_list, fill = TRUE)
    curves_sofar  <- rbindlist(curve_list, fill = TRUE)
    
    setorder(summary_sofar, heat_metric, green_metric_std)
    setorder(curves_sofar, heat_metric, green_metric_std, green_level, temp)
    
    fwrite(summary_sofar, out_summary)
    fwrite(curves_sofar, out_curves)
    
    cat("Saved progress for:", h, "x", gstd, "\n")
    cat("  cr_mod_p99 =", round(mod_p99$cr, 2), "\n")
    
    rm(sub, x, g, cbh, cbg, fit, beta, V, cp_low_iqr, cp_high_iqr, cp_p10, cp_p90, cp_int,
       curve_low_iqr, curve_high_iqr, curve_p10, curve_p90, summary_sofar, curves_sofar)
    gc()
  }
}

# 4) preview
summary_out <- fread(out_summary)

cat("\nSaved:\n")
cat(" -", out_summary, "\n")
cat(" -", out_curves, "\n")

cat("\nPreview:\n")
print(summary_out[, .(
  heat_metric, green_metric_raw,
  heat_mmt, heat_p99,
  cr_mod_p99, cr_mod_1iqr, cr_mod_2iqr
)])
