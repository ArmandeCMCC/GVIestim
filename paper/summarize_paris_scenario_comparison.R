# summarize_paris_scenario_comparison.R
#
# Combine all scenario summary files into one comparison table for paper
# reporting and ranking. Main paper run tag: _rerun_rule1_sum_ghsheat_ghsgreen.
#
# Stack multiple Paris attributable-deaths scenario runs into one
# comparison table.
#
# What it does:
#   1) Finds all files like:
#        paris_attr_scenario_summary*.csv
#   2) Reads and stacks them
#   3) Deduplicates repeated runs of the same specification by keeping
#      the most recently modified file
#   4) Creates:
#        - a long master table
#        - a wide publication-style comparison table
#        - a ranked table (best avoided-deaths scenarios first)
#
# Expected inputs:
#   files saved by estimate_paris_attributable_deaths_scenarios.R, e.g.
#     paris_attr_scenario_summary.csv
#     paris_attr_scenario_summary_t2mmean_imuvegtotal_p22.csv
#     paris_attr_scenario_summary_wbgtmean_imuhigh_mmt.csv
#
# Outputs:
#   - paris_attr_scenario_comparison_long.csv
#   - paris_attr_scenario_comparison_wide.csv
#   - paris_attr_scenario_comparison_ranked.csv

library(data.table)

base_dir <- "/Users/armandeaboudrar-meda/Desktop/GVIestim/paper" # relative : base_dir <- here::here() 

norm_tag <- function(x) {
  x <- trimws(x)
  if (!nzchar(x)) return("")
  if (!startsWith(x, "_")) x <- paste0("_", x)
  x
}

run_tag <- norm_tag(Sys.getenv("RUN_TAG", ""))

tag_file <- function(stem, tag, ext = ".csv") {
  file.path(base_dir, paste0(stem, tag, ext))
}

out_long <- tag_file("paris_attr_scenario_comparison_long", run_tag)
out_wide <- tag_file("paris_attr_scenario_comparison_wide", run_tag)
out_rank <- tag_file("paris_attr_scenario_comparison_ranked", run_tag)

# Backward-compatibility defaults for legacy scenario files that predate
# standardised per-100k outputs.
default_population_denominator <- 2200000
default_n_summers <- 10

# 1) find scenario summary files (filtered by RUN_TAG if set)
if (nzchar(run_tag)) {
  # Escape special regex chars in tag, then match files ending with tag.csv
  tag_escaped <- gsub("([.\\\\|()[{^$*+?])", "\\\\\\1", run_tag)
  file_pattern <- paste0("^paris_attr_scenario_summary.*", tag_escaped, "\\.csv$")
} else {
  file_pattern <- "^paris_attr_scenario_summary.*\\.csv$"
}

files <- list.files(
  path = base_dir,
  pattern = file_pattern,
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No files matching ", file_pattern, " were found in base_dir.")
}

cat("RUN_TAG filter:", ifelse(nzchar(run_tag), run_tag, "<none>"), "\n")
cat("Found", length(files), "scenario summary files:\n")
print(basename(files))

# 2) read and stack
need_cols <- c(
  "selected_heat", "selected_green_std", "selected_green_raw",
  "scenario",
  "heat_mmt", "heat_p99", "heat_iqr", "attribution_threshold",
  "green_median", "green_iqr", "green_p90", "green_max",
  "population_denominator", "n_summers",
  "n_rows", "n_hot_days_rows",
  "total_deaths", "total_deaths_hot_days",
  "ad_obs", "ad_scn", "ad_avoided",
  "pct_reduction_vs_observed_attr",
  "ad_obs_per100k_per_summer",
  "ad_scn_per100k_per_summer",
  "ad_avoided_per100k_per_summer"
)

new_cols <- c(
  "population_denominator", "n_summers",
  "ad_obs_per100k_per_summer",
  "ad_scn_per100k_per_summer",
  "ad_avoided_per100k_per_summer"
)

read_one <- function(path) {
  dt <- fread(path)
  
  miss <- setdiff(setdiff(need_cols, new_cols), names(dt))
  if (length(miss) > 0) {
    stop(
      "File is missing required columns:\n  ",
      basename(path), "\nMissing: ", paste(miss, collapse = ", ")
    )
  }

  if (!("population_denominator" %in% names(dt))) {
    dt[, population_denominator := default_population_denominator]
  }
  if (!("n_summers" %in% names(dt))) {
    dt[, n_summers := default_n_summers]
  }

  dt[, population_denominator := suppressWarnings(as.numeric(population_denominator))]
  dt[, n_summers := suppressWarnings(as.numeric(n_summers))]
  dt[!is.finite(population_denominator) | population_denominator <= 0, population_denominator := default_population_denominator]
  dt[!is.finite(n_summers) | n_summers <= 0, n_summers := default_n_summers]

  if (!("ad_obs_per100k_per_summer" %in% names(dt))) {
    dt[, ad_obs_per100k_per_summer := 1e5 * ad_obs / population_denominator / n_summers]
  }
  if (!("ad_scn_per100k_per_summer" %in% names(dt))) {
    dt[, ad_scn_per100k_per_summer := 1e5 * ad_scn / population_denominator / n_summers]
  }
  if (!("ad_avoided_per100k_per_summer" %in% names(dt))) {
    dt[, ad_avoided_per100k_per_summer := 1e5 * ad_avoided / population_denominator / n_summers]
  }
  
  info <- file.info(path)
  
  dt[, source_file := basename(path)]
  dt[, source_path := path]
  dt[, source_mtime := as.POSIXct(info$mtime, origin = "1970-01-01")]
  
  dt
}

all_dt <- rbindlist(lapply(files, read_one), fill = TRUE, use.names = TRUE)

# 3) deduplicate repeated runs
#    keep the most recent file for each specification × scenario
key_cols <- c(
  "selected_heat",
  "selected_green_std",
  "selected_green_raw",
  "attribution_threshold",
  "scenario"
)

setorderv(all_dt, c(key_cols, "source_mtime"), order = c(rep(1L, length(key_cols)), -1L))
all_dt <- all_dt[, .SD[1], by = key_cols]

# 4) clean labels / helpers
heat_labels <- c(
  t2m_min   = "T2M min",
  t2m_mean  = "T2M mean",
  t2m_max   = "T2M max",
  lst_min   = "LST min",
  lst_mean  = "LST mean",
  lst_max   = "LST max",
  wbgt_min  = "WBGT min",
  wbgt_mean = "WBGT mean",
  wbgt_max  = "WBGT max"
)

green_labels <- c(
  gvi_mean           = "GVI",
  gvi_popw_points    = "GVI (point pop-w.)",
  imu_veg_total      = "IMU total",
  imu_veg_total_popw = "IMU total (pop-w.)",
  imu_high           = "IMU high vegetation",
  imu_high_popw      = "IMU high vegetation (pop-w.)",
  imu_low            = "IMU low vegetation",
  imu_low_popw       = "IMU low vegetation (pop-w.)"
)

scenario_labels <- c(
  observed      = "Observed",
  max_all       = "Set all arr. to max greenness",
  p90_floor     = "Raise all arr. below p90 to p90",
  plus1iqr_cap  = "Add +1 greenness IQR (capped at max)"
)

all_dt[, heat_label := fifelse(
  selected_heat %in% names(heat_labels),
  heat_labels[selected_heat],
  selected_heat
)]

all_dt[, green_label := fifelse(
  selected_green_raw %in% names(green_labels),
  green_labels[selected_green_raw],
  selected_green_raw
)]

all_dt[, scenario_label := fifelse(
  scenario %in% names(scenario_labels),
  scenario_labels[scenario],
  scenario
)]

all_dt[, threshold_label := fifelse(
  abs(attribution_threshold - heat_mmt) < 1e-8,
  "MMT",
  paste0(">", format(round(attribution_threshold, 3), nsmall = 0, trim = TRUE), "C")
)]

all_dt[, spec_id := paste(
  selected_heat,
  selected_green_raw,
  paste0("thr", format(round(attribution_threshold, 3), trim = TRUE, scientific = FALSE)),
  sep = "__"
)]

main_green_metrics <- c(
  "imu_high_popw",
  "imu_veg_total_popw",
  "gvi_popw_points"
)

all_dt_main <- all_dt[selected_green_raw %in% main_green_metrics]

# 5) long master table
long_dt <- copy(all_dt_main)

setcolorder(long_dt, c(
  "spec_id",
  "heat_label", "green_label", "threshold_label",
  "selected_heat", "selected_green_std", "selected_green_raw",
  "scenario", "scenario_label",
  "heat_mmt", "heat_p99", "heat_iqr", "attribution_threshold",
  "green_median", "green_iqr", "green_p90", "green_max",
  "population_denominator", "n_summers",
  "n_rows", "n_hot_days_rows",
  "total_deaths", "total_deaths_hot_days",
  "ad_obs", "ad_scn", "ad_avoided", "pct_reduction_vs_observed_attr",
  "ad_obs_per100k_per_summer",
  "ad_scn_per100k_per_summer",
  "ad_avoided_per100k_per_summer",
  "source_file", "source_mtime"
))

setorder(long_dt, heat_label, green_label, attribution_threshold, scenario)

# 6) wide publication-style table
#    One row per specification, columns for each scenario
id_cols <- c(
  "spec_id",
  "heat_label", "green_label", "threshold_label",
  "selected_heat", "selected_green_std", "selected_green_raw",
  "heat_mmt", "heat_p99", "heat_iqr", "attribution_threshold",
  "green_median", "green_iqr", "green_p90", "green_max",
  "population_denominator", "n_summers"
)

wide_base <- unique(long_dt[, ..id_cols])

get_scenario_cols <- function(scn) {
  tmp <- long_dt[scenario == scn, .(
    spec_id,
    total_deaths,
    total_deaths_hot_days,
    ad_obs,
    ad_scn,
    ad_avoided,
    pct_reduction_vs_observed_attr,
    ad_obs_per100k_per_summer,
    ad_scn_per100k_per_summer,
    ad_avoided_per100k_per_summer
  )]
  
  setnames(
    tmp,
    old = c(
      "total_deaths", "total_deaths_hot_days", "ad_obs", "ad_scn", "ad_avoided",
      "pct_reduction_vs_observed_attr",
      "ad_obs_per100k_per_summer",
      "ad_scn_per100k_per_summer",
      "ad_avoided_per100k_per_summer"
    ),
    new = paste0(
      c(
        "total_deaths", "total_deaths_hot_days", "ad_obs", "ad_scn", "ad_avoided", "pct_reduction",
        "ad_obs_per100k_per_summer",
        "ad_scn_per100k_per_summer",
        "ad_avoided_per100k_per_summer"
      ),
      "_", scn
    )
  )
  
  tmp
}

wide_dt <- copy(wide_base)
for (scn in c("observed", "max_all", "p90_floor", "plus1iqr_cap")) {
  wide_dt <- merge(wide_dt, get_scenario_cols(scn), by = "spec_id", all.x = TRUE)
}

# Add a few especially useful compact columns
wide_dt[, avoided_best := pmax(
  ad_avoided_max_all,
  ad_avoided_p90_floor,
  ad_avoided_plus1iqr_cap,
  na.rm = TRUE
)]

wide_dt[, pct_reduction_best := pmax(
  pct_reduction_max_all,
  pct_reduction_p90_floor,
  pct_reduction_plus1iqr_cap,
  na.rm = TRUE
)]

setorder(wide_dt, heat_label, green_label, attribution_threshold)

# 7) ranked table
#    rank non-observed scenarios by avoided deaths
rank_dt <- long_dt[scenario != "observed", .(
  spec_id,
  heat_label,
  green_label,
  threshold_label,
  scenario,
  scenario_label,
  heat_mmt,
  heat_p99,
  heat_iqr,
  attribution_threshold,
  green_median,
  green_iqr,
  green_p90,
  green_max,
  population_denominator,
  n_summers,
  total_deaths,
  total_deaths_hot_days,
  ad_obs,
  ad_scn,
  ad_avoided,
  pct_reduction_vs_observed_attr,
  ad_obs_per100k_per_summer,
  ad_scn_per100k_per_summer,
  ad_avoided_per100k_per_summer,
  source_file,
  source_mtime
)]

setorder(rank_dt, -ad_avoided, -pct_reduction_vs_observed_attr)

# 8) save
fwrite(long_dt, out_long)
fwrite(wide_dt, out_wide)
fwrite(rank_dt, out_rank)

cat("\nSaved:\n")
cat(" -", out_long, "\n")
cat(" -", out_wide, "\n")
cat(" -", out_rank, "\n")

# 9) console preview
cat("\nSpecifications found:\n")
print(unique(long_dt[, .(
  heat_label,
  green_label,
  threshold_label,
  selected_heat,
  selected_green_raw,
  attribution_threshold
)]))

cat("\nWide comparison preview:\n")
print(wide_dt[, .(
  heat_label,
  green_label,
  threshold_label,
  ad_obs_observed,
  ad_avoided_max_all,
  pct_reduction_max_all,
  ad_avoided_per100k_per_summer_max_all,
  ad_avoided_p90_floor,
  pct_reduction_p90_floor,
  ad_avoided_per100k_per_summer_p90_floor,
  ad_avoided_plus1iqr_cap,
  pct_reduction_plus1iqr_cap,
  ad_avoided_per100k_per_summer_plus1iqr_cap
)])

cat("\nTop scenario results:\n")
print(rank_dt[1:min(15, .N), .(
  heat_label,
  green_label,
  threshold_label,
  scenario_label,
  ad_obs,
  ad_scn,
  ad_avoided,
  pct_reduction_vs_observed_attr,
  ad_avoided_per100k_per_summer
)])

cat("\nDone.\n")
