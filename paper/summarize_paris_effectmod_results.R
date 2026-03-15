# summarize_paris_effectmod_results.R
#
# Convert raw effectmod outputs into paper-ready estimands/tables.
# Main paper run tag: _rerun_rule1_sum_ghsheat_ghsgreen.
#
# Post-process Paris heat × greenness effect-modification results

# Inputs:
#   - paris_effectmod_heat_green_summary.csv
#   - paris_greenness_metric_summary.csv

# Main goals:
#   1) Keep the PRIMARY estimand:
#        greenness from p10 to p90
#   2) Keep the harmonized sensitivity estimand:
#        greenness from -1 IQR to +1 IQR
#   3) Add a paper-style companion estimand:
#        greenness from observed min to observed max
#   4) Create clean tables for:
#        - p99 endpoint
#        - MMT + 1 heat-IQR endpoint
#        - MMT + 2 heat-IQR endpoint

# Outputs:
#   - paris_effectmod_results_harmonized.csv
#   - weighted/extended main set:
#       paris_effectmod_primary_table_p99.csv
#       paris_effectmod_primary_table_1iqr.csv
#       paris_effectmod_primary_table_2iqr.csv
#       paris_effectmod_ranked_p99.csv
#       paris_effectmod_paper_table_p99_minmax.csv
#       paris_effectmod_ranked_p99_minmax.csv
#   - native main set (Hong Kong-style comparison):
#       paris_effectmod_primary_table_p99_native.csv
#       paris_effectmod_primary_table_1iqr_native.csv
#       paris_effectmod_primary_table_2iqr_native.csv
#       paris_effectmod_ranked_p99_native.csv
#       paris_effectmod_paper_table_p99_minmax_native.csv
#       paris_effectmod_ranked_p99_minmax_native.csv

library(data.table)

base_dir <- "/Users/armandeaboudrar-meda/Desktop/GVIestim/paper"

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
effmod_tag <- norm_tag(Sys.getenv("EFFMOD_TAG", run_tag))
green_tag <- norm_tag(Sys.getenv("GREEN_TAG", effmod_tag))

in_effmod <- tag_file("paris_effectmod_heat_green_summary", effmod_tag)
in_green_try <- tag_file("paris_greenness_metric_summary", green_tag)
in_green_default <- file.path(base_dir, "paris_greenness_metric_summary.csv")
in_green <- if (file.exists(in_green_try)) in_green_try else in_green_default

out_all   <- tag_file("paris_effectmod_results_harmonized", run_tag)
out_p99   <- tag_file("paris_effectmod_primary_table_p99", run_tag)
out_1iqr  <- tag_file("paris_effectmod_primary_table_1iqr", run_tag)
out_2iqr  <- tag_file("paris_effectmod_primary_table_2iqr", run_tag)
out_rank  <- tag_file("paris_effectmod_ranked_p99", run_tag)
out_p99_paper  <- tag_file("paris_effectmod_paper_table_p99_minmax", run_tag)
out_rank_paper <- tag_file("paris_effectmod_ranked_p99_minmax", run_tag)

out_p99_native <- tag_file("paris_effectmod_primary_table_p99_native", run_tag)
out_1iqr_native <- tag_file("paris_effectmod_primary_table_1iqr_native", run_tag)
out_2iqr_native <- tag_file("paris_effectmod_primary_table_2iqr_native", run_tag)
out_rank_native <- tag_file("paris_effectmod_ranked_p99_native", run_tag)
out_p99_paper_native <- tag_file("paris_effectmod_paper_table_p99_minmax_native", run_tag)
out_rank_paper_native <- tag_file("paris_effectmod_ranked_p99_minmax_native", run_tag)

cat("\nSummarize effect-modification config:\n")
cat(" - EFFMOD_TAG:", ifelse(nzchar(effmod_tag), effmod_tag, "<none>"), "\n")
cat(" - GREEN_TAG:", ifelse(nzchar(green_tag), green_tag, "<none>"), "\n")
cat(" - RUN_TAG:", ifelse(nzchar(run_tag), run_tag, "<none>"), "\n")
cat(" - input effectmod:", in_effmod, "\n")
cat(" - input greenness summary:", in_green, "\n")

# 1) load
eff <- fread(in_effmod)
green <- fread(in_green)

stopifnot(all(c(
  "heat_metric", "green_metric_raw",
  "rr_ratio_p99", "rr_ratio_p99_low", "rr_ratio_p99_high",
  "rr_ratio_1iqr", "rr_ratio_1iqr_low", "rr_ratio_1iqr_high",
  "rr_ratio_2iqr", "rr_ratio_2iqr_low", "rr_ratio_2iqr_high"
) %in% names(eff)))

stopifnot(all(c("metric", "min", "p10", "median", "p90", "max", "iqr") %in% names(green)))

green <- green[type == "greenness"]
green <- green[, .(
  green_metric_raw = metric,
  green_min = min,
  green_p10 = p10,
  green_median = median,
  green_p90 = p90,
  green_max = max,
  green_iqr_check = iqr
)]

green_order_main <- c(
  "gvi_popw_points",
  "ndvi_popw_ghs",
  "imu_veg_total_popw_ghs",
  "imu_high_popw_ghs",
  "imu_low_popw_ghs"
)

green_order_native_main <- c(
  "gvi_mean",
  "ndvi_native_jjas",
  "imu_veg_total"
)

# 2) join greenness distribution summaries
dt <- merge(eff, green, by = "green_metric_raw", all.x = TRUE)

if (any(is.na(dt$green_iqr_check))) {
  stop("Some greenness metrics in effectmod results were not found in greenness summary.")
}

# These models were fit on standardized greenness:
#   G_std = (G - median(G)) / IQR(G)
#
# The fitted rr_ratio_* in eff correspond to a move from -1 std-IQR to +1 std-IQR,
# i.e. delta_std = 2.
#
# We now rescale that to:
#   - min -> max
#   - p10 -> p90
#
# delta_std(min->max) = (max - min) / IQR
# delta_std(p10->p90) = (p90 - p10) / IQR

dt[, delta_std_minmax := (green_max - green_min) / green_iqr_check]
dt[, delta_std_p10p90 := (green_p90 - green_p10) / green_iqr_check]
dt[, iqr_low_below_observed_min := green_raw_low < green_min]
dt[, iqr_high_above_observed_max := green_raw_high > green_max]

# 3) helper to rescale RR contrasts from delta_std = 2 to arbitrary delta_std
rescale_rr_contrast <- function(rr2, low2, high2, delta_to, delta_from = 2) {
  # rr2 is the RR ratio for delta_from standardized units (here 2)
  # We assume log-linearity in the modifier contrast:
  #   log(RR_delta) = (delta_to / delta_from) * log(RR_delta_from)
  #
  # Recover SE on log scale from the reported CI
  log_rr2 <- log(rr2)
  se2 <- (log(high2) - log(low2)) / (2 * 1.96)
  
  scale_factor <- delta_to / delta_from
  
  log_rr_to <- scale_factor * log_rr2
  se_to <- abs(scale_factor) * se2
  
  rr_to <- exp(log_rr_to)
  low_to <- exp(log_rr_to - 1.96 * se_to)
  high_to <- exp(log_rr_to + 1.96 * se_to)
  
  list(
    rr = rr_to,
    low = low_to,
    high = high_to,
    cr = (rr_to - 1) * 100,
    cr_low = (low_to - 1) * 100,
    cr_high = (high_to - 1) * 100
  )
}

# 4) compute paper-style and p10-p90 sensitivities for each endpoint
add_rescaled_endpoint <- function(DT, prefix_in, prefix_out, delta_var) {
  rr_col   <- paste0("rr_ratio_", prefix_in)
  low_col  <- paste0("rr_ratio_", prefix_in, "_low")
  high_col <- paste0("rr_ratio_", prefix_in, "_high")
  
  out_rr   <- paste0("rr_ratio_", prefix_out)
  out_low  <- paste0("rr_ratio_", prefix_out, "_low")
  out_high <- paste0("rr_ratio_", prefix_out, "_high")
  out_cr   <- paste0("cr_mod_", prefix_out)
  out_cr_l <- paste0("cr_mod_", prefix_out, "_low")
  out_cr_h <- paste0("cr_mod_", prefix_out, "_high")
  
  DT[, c(out_rr, out_low, out_high, out_cr, out_cr_l, out_cr_h) := {
    res <- rescale_rr_contrast(
      rr2 = get(rr_col),
      low2 = get(low_col),
      high2 = get(high_col),
      delta_to = get(delta_var),
      delta_from = 2
    )
    list(res$rr, res$low, res$high, res$cr, res$cr_low, res$cr_high)
  }]
}

# min -> max
add_rescaled_endpoint(dt, "p99",  "p99_minmax",  "delta_std_minmax")
add_rescaled_endpoint(dt, "1iqr", "1iqr_minmax", "delta_std_minmax")
add_rescaled_endpoint(dt, "2iqr", "2iqr_minmax", "delta_std_minmax")

# p10 -> p90
add_rescaled_endpoint(dt, "p99",  "p99_p10p90",  "delta_std_p10p90")
add_rescaled_endpoint(dt, "1iqr", "1iqr_p10p90", "delta_std_p10p90")
add_rescaled_endpoint(dt, "2iqr", "2iqr_p10p90", "delta_std_p10p90")

# 5) significance / direction flags
flag_effect <- function(low, high) {
  out <- rep("uncertain", length(low))
  out[is.finite(low) & is.finite(high) & high < 0] <- "protective"
  out[is.finite(low) & is.finite(high) & low > 0]  <- "harmful"
  out
}

dt[, mod_flag_p99 := flag_effect(cr_mod_p99_low, cr_mod_p99_high)]
dt[, mod_flag_p99_p10p90 := flag_effect(cr_mod_p99_p10p90_low, cr_mod_p99_p10p90_high)]
dt[, mod_flag_1iqr := flag_effect(cr_mod_1iqr_low, cr_mod_1iqr_high)]
dt[, mod_flag_2iqr := flag_effect(cr_mod_2iqr_low, cr_mod_2iqr_high)]

# ordering helpers
dt[, abs_cr_mod_p99_p10p90 := abs(cr_mod_p99_p10p90)]
dt[, abs_cr_mod_p99 := abs(cr_mod_p99)]
dt[, abs_cr_mod_2iqr := abs(cr_mod_2iqr)]

# clean tables
make_endpoint_table <- function(DT, endpoint = c("p99", "1iqr", "2iqr")) {
  endpoint <- match.arg(endpoint)
  
  cr_col      <- paste0("cr_mod_", endpoint)
  cr_low_col  <- paste0("cr_mod_", endpoint, "_low")
  cr_high_col <- paste0("cr_mod_", endpoint, "_high")
  flag_col    <- paste0("mod_flag_", endpoint)
  
  out <- copy(DT)[, .(
    heat_metric,
    green_metric_raw,
    
    heat_mmt,
    heat_p99,
    heat_iqr,
    
    green_raw_median,
    green_raw_iqr,
    green_raw_low,
    green_raw_high,
    green_min,
    green_p10,
    green_p90,
    green_max,
    iqr_low_below_observed_min,
    iqr_high_above_observed_max,
    
    # harmonised contrast used in model output: -1 IQR -> +1 IQR
    cr_mod_iqr = get(cr_col),
    cr_mod_iqr_low = get(cr_low_col),
    cr_mod_iqr_high = get(cr_high_col),
    mod_flag_iqr = get(flag_col),
    
    # sensitivity: min -> max
    cr_mod_minmax = get(paste0("cr_mod_", endpoint, "_minmax")),
    cr_mod_minmax_low = get(paste0("cr_mod_", endpoint, "_minmax_low")),
    cr_mod_minmax_high = get(paste0("cr_mod_", endpoint, "_minmax_high")),
    
    # robust raw-contrast: p10 -> p90 (PRIMARY estimand for cross-metric comparison)
    cr_mod_p10p90 = get(paste0("cr_mod_", endpoint, "_p10p90")),
    cr_mod_p10p90_low = get(paste0("cr_mod_", endpoint, "_p10p90_low")),
    cr_mod_p10p90_high = get(paste0("cr_mod_", endpoint, "_p10p90_high"))
  )]
  
  out[, `:=`(
    primary_estimand = "p10_to_p90",
    cr_mod_primary = cr_mod_p10p90,
    cr_mod_primary_low = cr_mod_p10p90_low,
    cr_mod_primary_high = cr_mod_p10p90_high
  )]
  
  out[, mod_flag_primary := flag_effect(cr_mod_primary_low, cr_mod_primary_high)]
  
  setorder(out, heat_metric, green_metric_raw)
  out
}

add_fdr_primary <- function(tab) {
  tab <- copy(tab)
  tab[, `:=`(
    log_rr_primary = log1p(cr_mod_primary / 100),
    se_log_rr_primary = (log1p(cr_mod_primary_high / 100) - log1p(cr_mod_primary_low / 100)) / (2 * 1.96)
  )]
  tab[, z_primary := fifelse(is.finite(se_log_rr_primary) & se_log_rr_primary > 0,
                             log_rr_primary / se_log_rr_primary, NA_real_)]
  tab[, p_primary := fifelse(is.finite(z_primary), 2 * pnorm(-abs(z_primary)), NA_real_)]
  tab[, q_primary_bh := p.adjust(p_primary, method = "BH")]
  tab
}

tab_p99_all  <- make_endpoint_table(dt, "p99")
tab_1iqr_all <- make_endpoint_table(dt, "1iqr")
tab_2iqr_all <- make_endpoint_table(dt, "2iqr")

tab_p99  <- tab_p99_all[green_metric_raw %in% green_order_main]
tab_1iqr <- tab_1iqr_all[green_metric_raw %in% green_order_main]
tab_2iqr <- tab_2iqr_all[green_metric_raw %in% green_order_main]

tab_p99[, green_order_id := match(green_metric_raw, green_order_main)]
tab_1iqr[, green_order_id := match(green_metric_raw, green_order_main)]
tab_2iqr[, green_order_id := match(green_metric_raw, green_order_main)]

setorder(tab_p99, heat_metric, green_order_id)
setorder(tab_1iqr, heat_metric, green_order_id)
setorder(tab_2iqr, heat_metric, green_order_id)

tab_p99[, green_order_id := NULL]
tab_1iqr[, green_order_id := NULL]
tab_2iqr[, green_order_id := NULL]

tab_p99 <- add_fdr_primary(tab_p99)
tab_1iqr <- add_fdr_primary(tab_1iqr)
tab_2iqr <- add_fdr_primary(tab_2iqr)

# ranked p99 table: strongest protective first under PRIMARY estimand
rank_p99 <- copy(tab_p99)
setorder(rank_p99, cr_mod_p10p90)

# dedicated Paris-paper estimand table at p99: min -> max
tab_p99_paper <- copy(tab_p99)[, .(
  heat_metric,
  green_metric_raw,
  heat_mmt,
  heat_p99,
  heat_iqr,
  green_min,
  green_p10,
  green_p90,
  green_max,
  cr_mod_paper = cr_mod_minmax,
  cr_mod_paper_low = cr_mod_minmax_low,
  cr_mod_paper_high = cr_mod_minmax_high
)]
tab_p99_paper[, mod_flag_paper := flag_effect(cr_mod_paper_low, cr_mod_paper_high)]
tab_p99_paper[, green_order_id := match(green_metric_raw, green_order_main)]
setorder(tab_p99_paper, heat_metric, green_order_id)
tab_p99_paper[, green_order_id := NULL]

rank_p99_paper <- copy(tab_p99_paper)
setorder(rank_p99_paper, cr_mod_paper)

# native-main counterparts (clean map-choice comparison)
tab_p99_native  <- tab_p99_all[green_metric_raw %in% green_order_native_main]
tab_1iqr_native <- tab_1iqr_all[green_metric_raw %in% green_order_native_main]
tab_2iqr_native <- tab_2iqr_all[green_metric_raw %in% green_order_native_main]

tab_p99_native[, green_order_id := match(green_metric_raw, green_order_native_main)]
tab_1iqr_native[, green_order_id := match(green_metric_raw, green_order_native_main)]
tab_2iqr_native[, green_order_id := match(green_metric_raw, green_order_native_main)]

setorder(tab_p99_native, heat_metric, green_order_id)
setorder(tab_1iqr_native, heat_metric, green_order_id)
setorder(tab_2iqr_native, heat_metric, green_order_id)

tab_p99_native[, green_order_id := NULL]
tab_1iqr_native[, green_order_id := NULL]
tab_2iqr_native[, green_order_id := NULL]

tab_p99_native <- add_fdr_primary(tab_p99_native)
tab_1iqr_native <- add_fdr_primary(tab_1iqr_native)
tab_2iqr_native <- add_fdr_primary(tab_2iqr_native)

rank_p99_native <- copy(tab_p99_native)
setorder(rank_p99_native, cr_mod_p10p90)

tab_p99_paper_native <- copy(tab_p99_native)[, .(
  heat_metric,
  green_metric_raw,
  heat_mmt,
  heat_p99,
  heat_iqr,
  green_min,
  green_p10,
  green_p90,
  green_max,
  cr_mod_paper = cr_mod_minmax,
  cr_mod_paper_low = cr_mod_minmax_low,
  cr_mod_paper_high = cr_mod_minmax_high
)]
tab_p99_paper_native[, mod_flag_paper := flag_effect(cr_mod_paper_low, cr_mod_paper_high)]
tab_p99_paper_native[, green_order_id := match(green_metric_raw, green_order_native_main)]
setorder(tab_p99_paper_native, heat_metric, green_order_id)
tab_p99_paper_native[, green_order_id := NULL]

rank_p99_paper_native <- copy(tab_p99_paper_native)
setorder(rank_p99_paper_native, cr_mod_paper)

# 7) save
setorder(dt, heat_metric, green_metric_raw)

fwrite(dt, out_all)
fwrite(tab_p99, out_p99)
fwrite(tab_1iqr, out_1iqr)
fwrite(tab_2iqr, out_2iqr)
fwrite(rank_p99, out_rank)
fwrite(tab_p99_paper, out_p99_paper)
fwrite(rank_p99_paper, out_rank_paper)
fwrite(tab_p99_native, out_p99_native)
fwrite(tab_1iqr_native, out_1iqr_native)
fwrite(tab_2iqr_native, out_2iqr_native)
fwrite(rank_p99_native, out_rank_native)
fwrite(tab_p99_paper_native, out_p99_paper_native)
fwrite(rank_p99_paper_native, out_rank_paper_native)

cat("\nSaved:\n")
cat(" -", out_all, "\n")
cat(" -", out_p99, "\n")
cat(" -", out_1iqr, "\n")
cat(" -", out_2iqr, "\n")
cat(" -", out_rank, "\n")
cat(" -", out_p99_paper, "\n")
cat(" -", out_rank_paper, "\n")
cat(" -", out_p99_native, "\n")
cat(" -", out_1iqr_native, "\n")
cat(" -", out_2iqr_native, "\n")
cat(" -", out_rank_native, "\n")
cat(" -", out_p99_paper_native, "\n")
cat(" -", out_rank_paper_native, "\n")

# 8) preview
cat("\nTop protective results at p99 (PRIMARY p10 -> p90 contrast):\n")
print(rank_p99[1:12, .(
  heat_metric, green_metric_raw,
  cr_mod_primary, cr_mod_primary_low, cr_mod_primary_high,
  q_primary_bh,
  mod_flag_primary,
  cr_mod_iqr, cr_mod_iqr_low, cr_mod_iqr_high,
  cr_mod_minmax, cr_mod_minmax_low, cr_mod_minmax_high,
  iqr_low_below_observed_min
)])

cat("\nTop protective results at p99 (PRIMARY p10 -> p90, native-main set):\n")
print(rank_p99_native[1:12, .(
  heat_metric, green_metric_raw,
  cr_mod_primary, cr_mod_primary_low, cr_mod_primary_high,
  q_primary_bh,
  mod_flag_primary
)])

cat("\nTop protective results at p99 (Paris-paper min -> max contrast):\n")
print(rank_p99_paper[1:12, .(
  heat_metric, green_metric_raw,
  cr_mod_paper, cr_mod_paper_low, cr_mod_paper_high,
  mod_flag_paper
)])

cat("\nDone.\n")
