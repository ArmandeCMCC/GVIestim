# diagnose_residual_autocorrelation.R
#
# checking whether quasi-Poisson model-based SEs are adequate by testing
# for residual within-stratum serial correlation in deviance residuals.
#
# If the median lag-1 ACF is small (< ~0.10), the quasi-Poisson SEs are
# adequate and no clustering adjustment is needed.
#
# Runs on the MAIN CTS file by default. Set RUN_TAG to check the native
# appendix run.
#
# Outputs:
#   - Printed summary of lag-1 to lag-3 ACF across strata
#   - figs/paris_residual_acf_diagnostic{RUN_TAG}.pdf

library(data.table)
library(lubridate)
library(splines)
library(gnm)
library(dlnm)
library(ggplot2)

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

run_tag <- norm_tag(Sys.getenv("RUN_TAG", "rerun_rule1_sum_ghsheat_ghsgreen"))
cts_tag <- norm_tag(Sys.getenv("CTS_TAG", run_tag))

in_path <- tag_file("paris_cts_ready_jjas_2008_2017", cts_tag)
cat("Input CTS:", in_path, "\n")

fig_dir <- file.path(base_dir, "figs")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# 1) Load and prepare
dt <- fread(in_path)
dt[, date := as.Date(date)]
dt <- dt[month(date) %in% 6:9]
dt[, stratum_f := factor(stratum_id)]
dt[, dow_f := factor(dow_f)]
dt[, year_f := factor(year_f)]

heat_metrics <- intersect(
  c("t2m_min", "t2m_mean", "t2m_max",
    "lst_min", "lst_mean", "lst_max",
    "wbgt_min", "wbgt_mean", "wbgt_max"),
  names(dt)
)

# 2) Fit baseline models and extract residual ACF
acf_results <- list()

for (metric in heat_metrics) {
  cat("\n--- Diagnostic for:", metric, "---\n")

  sub <- copy(dt[!is.na(get(metric))])
  if (nrow(sub) == 0) next

  x <- sub[[metric]]
  metric_p90 <- as.numeric(quantile(x, 0.90, na.rm = TRUE, type = 7))

  cb <- crossbasis(
    x,
    lag = 1,
    argvar = list(fun = "ns", knots = metric_p90),
    arglag = list(fun = "integer")
  )

  fit <- gnm(
    deaths ~ cb + dow_f + ns(day_of_season, df = 4):year_f,
    eliminate = stratum_f,
    family = quasipoisson(),
    data = sub,
    trace = FALSE
  )

  # Deviance residuals — gnm may drop observations (e.g. singleton strata),
  # so match back via the named index vector rather than positional assignment.
  res <- residuals(fit, type = "deviance")
  idx <- as.integer(names(res))
  sub[, dev_resid := NA_real_]
  sub[idx, dev_resid := res]

  # Dispersion parameter
  phi <- summary(fit)$dispersion
  cat("  Dispersion (phi):", round(phi, 3), "\n")

  # Compute lag-1 to lag-3 ACF within each stratum
  # (strata are arrondissement x year x month, so consecutive rows
  #  within a stratum are consecutive days)
  setorder(sub, stratum_id, date)

  acf_by_stratum <- sub[!is.na(dev_resid), {
    n <- .N
    if (n > 6) {
      a <- acf(dev_resid, lag.max = 3, plot = FALSE)$acf
      list(
        lag1 = a[2],
        lag2 = a[3],
        lag3 = a[4],
        n_days = n
      )
    } else {
      list(lag1 = NA_real_, lag2 = NA_real_, lag3 = NA_real_, n_days = n)
    }
  }, by = stratum_id]

  acf_by_stratum <- acf_by_stratum[!is.na(lag1)]
  acf_by_stratum[, heat_metric := metric]

  cat("  Strata checked:", nrow(acf_by_stratum), "\n")
  cat("  Lag-1 ACF — median:", round(median(acf_by_stratum$lag1), 4),
      " mean:", round(mean(acf_by_stratum$lag1), 4),
      " p95:", round(quantile(acf_by_stratum$lag1, 0.95), 4), "\n")
  cat("  Lag-2 ACF — median:", round(median(acf_by_stratum$lag2), 4),
      " mean:", round(mean(acf_by_stratum$lag2), 4), "\n")
  cat("  Lag-3 ACF — median:", round(median(acf_by_stratum$lag3), 4),
      " mean:", round(mean(acf_by_stratum$lag3), 4), "\n")

  acf_results[[metric]] <- acf_by_stratum
}

# 3) Combine and summarize
acf_all <- rbindlist(acf_results, fill = TRUE)

cat("\n============================================================\n")
cat("OVERALL SUMMARY across all heat metrics and strata\n")
cat("============================================================\n")

overall <- acf_all[, .(
  n_strata = .N,
  lag1_median = median(lag1),
  lag1_mean = mean(lag1),
  lag1_p05 = quantile(lag1, 0.05),
  lag1_p95 = quantile(lag1, 0.95),
  lag2_median = median(lag2),
  lag3_median = median(lag3)
), by = heat_metric]

print(overall)

cat("\nPooled across all metrics:\n")
cat("  Lag-1 median:", round(median(acf_all$lag1), 4), "\n")
cat("  Lag-1 mean:  ", round(mean(acf_all$lag1), 4), "\n")
cat("  Lag-1 [5%, 95%]: [",
    round(quantile(acf_all$lag1, 0.05), 4), ",",
    round(quantile(acf_all$lag1, 0.95), 4), "]\n")

# Verdict
med_lag1 <- median(acf_all$lag1)
if (abs(med_lag1) < 0.05) {
  cat("\n>> VERDICT: Negligible residual autocorrelation (median |lag-1| < 0.05).\n")
  cat("   Quasi-Poisson SEs are adequate.\n")
} else if (abs(med_lag1) < 0.10) {
  cat("\n>> VERDICT: Small residual autocorrelation (median |lag-1| < 0.10).\n")
  cat("   Quasi-Poisson SEs are likely adequate; mention in limitations.\n")
} else {
  cat("\n>> VERDICT: Non-trivial residual autocorrelation (median |lag-1| >= 0.10).\n")
  cat("   Consider cluster-robust SEs or sensitivity analysis.\n")
}

# 4) Diagnostic plot
acf_long <- melt(
  acf_all,
  id.vars = c("stratum_id", "heat_metric", "n_days"),
  measure.vars = c("lag1", "lag2", "lag3"),
  variable.name = "lag",
  value.name = "acf"
)
acf_long[, lag := factor(lag, levels = c("lag1", "lag2", "lag3"),
                         labels = c("Lag 1", "Lag 2", "Lag 3"))]

heat_labels <- c(
  t2m_min = "T2M min", t2m_mean = "T2M mean", t2m_max = "T2M max",
  lst_min = "LST min", lst_mean = "LST mean", lst_max = "LST max",
  wbgt_min = "WBGT min", wbgt_mean = "WBGT mean", wbgt_max = "WBGT max"
)
acf_long[, heat_label := heat_labels[heat_metric]]

p <- ggplot(acf_long, aes(x = acf)) +
  geom_histogram(bins = 40, fill = "steelblue", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = c(-0.10, 0.10), linetype = "dotted",
             linewidth = 0.3, colour = "red") +
  facet_grid(lag ~ heat_label) +
  labs(
    title = "Residual autocorrelation diagnostic (deviance residuals within CTS strata)",
    subtitle = paste0("Run: ", run_tag, ". Red dotted lines at +/- 0.10."),
    x = "Autocorrelation",
    y = "Count (strata)"
  ) +
  theme_bw(base_size = 9.5) +
  theme(
    strip.background = element_rect(fill = "grey95", colour = "grey75"),
    panel.grid.minor = element_blank()
  )

out_pdf <- file.path(fig_dir, paste0("paris_residual_acf_diagnostic", run_tag, ".pdf"))
out_png <- file.path(fig_dir, paste0("paris_residual_acf_diagnostic", run_tag, ".png"))
ggsave(out_pdf, p, width = 14, height = 6)
ggsave(out_png, p, width = 14, height = 6, dpi = 300)

cat("\nSaved diagnostic figure:\n -", out_pdf, "\n -", out_png, "\n")
cat("\nDone.\n")
