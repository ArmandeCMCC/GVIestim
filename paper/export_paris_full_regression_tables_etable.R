# export_paris_full_regression_tables_etable.R
#
# Full regression tables using fixest::etable
# - One table per heat metric
# - Three columns per table: GVI, NDVI, IMU total
# - Run-tag aware (default: final_ghsheat_ghsgreen)

library(data.table)
library(lubridate)
library(splines)
library(dlnm)
library(fixest)

base_dir <- "/Users/armandeaboudrar-meda/Desktop/GVIestim/paper"

norm_tag <- function(x) {
  x <- trimws(x)
  if (!nzchar(x)) return("")
  if (!startsWith(x, "_")) x <- paste0("_", x)
  x
}

run_tag <- norm_tag(Sys.getenv("RUN_TAG", "final_ghsheat_ghsgreen"))
cts_tag <- norm_tag(Sys.getenv("CTS_TAG", run_tag))

in_cts <- file.path(base_dir, paste0("paris_cts_ready_jjas_2008_2017", cts_tag, ".csv"))
if (!file.exists(in_cts)) stop("Missing CTS file: ", in_cts)

out_dir <- file.path(base_dir, "tables", "appendix", "full_regression_etable")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

to_longtable <- function(tex_path) {
  x <- readLines(tex_path, warn = FALSE)

  i_tab_begin <- grep("^\\s*\\\\begin\\{tabular\\}\\{", x)[1]
  i_tab_end <- tail(grep("^\\s*\\\\end\\{tabular\\}", x), 1)
  i_cap <- grep("^\\s*\\\\caption\\{", x)[1]

  if (!length(i_tab_begin) || is.na(i_tab_begin) || !length(i_tab_end) || is.na(i_tab_end)) {
    warning("Could not convert to longtable (tabular block not found): ", tex_path)
    return(invisible(FALSE))
  }

  tab_open <- x[i_tab_begin]
  tab_spec <- sub("^\\s*\\\\begin\\{tabular\\}(\\{.*\\})\\s*$", "\\1", tab_open)
  if (!nzchar(tab_spec) || identical(tab_spec, tab_open)) {
    warning("Could not parse tabular column spec: ", tex_path)
    return(invisible(FALSE))
  }

  cap_line <- if (length(i_cap) && !is.na(i_cap)) x[i_cap] else NULL
  body <- x[(i_tab_begin + 1):(i_tab_end - 1)]

  out <- c(
    paste0("\\begin{longtable}", tab_spec),
    if (!is.null(cap_line)) cap_line else NULL,
    if (!is.null(cap_line)) "\\\\" else NULL,
    body,
    "\\end{longtable}",
    ""
  )

  writeLines(out, tex_path)
  invisible(TRUE)
}

heat_metrics_all <- c(
  "t2m_min", "t2m_mean", "t2m_max",
  "lst_min", "lst_mean", "lst_max",
  "wbgt_min", "wbgt_mean", "wbgt_max"
)

heat_filter_raw <- trimws(Sys.getenv("HEAT_FILTER", ""))
if (nzchar(heat_filter_raw)) {
  heat_metrics <- trimws(strsplit(heat_filter_raw, ",", fixed = TRUE)[[1]])
  heat_metrics <- intersect(heat_metrics, heat_metrics_all)
  if (length(heat_metrics) == 0) stop("HEAT_FILTER provided but no valid heat metrics matched.")
} else {
  heat_metrics <- heat_metrics_all
}

green_metrics <- c(
  "gvi_popw_points_std_iqr",
  "ndvi_popw_ghs_std_iqr",
  "imu_veg_total_popw_ghs_std_iqr"
)

green_labels <- c(
  gvi_popw_points_std_iqr = "GVI",
  ndvi_popw_ghs_std_iqr = "NDVI",
  imu_veg_total_popw_ghs_std_iqr = "IMU total"
)

heat_labels <- c(
  t2m_min = "T2M min",
  t2m_mean = "T2M mean",
  t2m_max = "T2M max",
  lst_min = "LST min",
  lst_mean = "LST mean",
  lst_max = "LST max",
  wbgt_min = "WBGT min",
  wbgt_mean = "WBGT mean",
  wbgt_max = "WBGT max"
)

dt <- fread(in_cts)
dt[, date := as.Date(date)]
dt <- dt[month(date) %in% 6:9]

if (!("stratum_id" %in% names(dt))) stop("Missing stratum_id in CTS file.")
dt[, stratum_f := factor(stratum_id)]

if ("dow_f" %in% names(dt)) {
  dt[, dow_f := factor(as.character(dow_f), levels = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"))]
} else if ("dow" %in% names(dt)) {
  dt[, dow_f := factor(dow, levels = 1:7, labels = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"))]
} else {
  stop("Need dow_f or dow in CTS file.")
}

if ("year" %in% names(dt)) {
  dt[, year_f := factor(year)]
} else if ("year_f" %in% names(dt)) {
  dt[, year_f := factor(as.character(year_f))]
} else {
  stop("Need year or year_f in CTS file.")
}

need_cols <- unique(c("deaths", "day_of_season", heat_metrics_all, green_metrics))
miss_cols <- setdiff(need_cols, names(dt))
if (length(miss_cols) > 0) stop("CTS file missing columns: ", paste(miss_cols, collapse = ", "))

tex_files <- character()

for (h in heat_metrics) {
  cat("Running etable export for heat metric: ", h, "\n", sep = "")

  fits <- list()

  for (gstd in green_metrics) {
    sub <- copy(dt[!is.na(get(h)) & !is.na(get(gstd))])
    if (nrow(sub) == 0) {
      warning("No rows for ", h, " x ", gstd)
      next
    }

    h_p90 <- as.numeric(quantile(sub[[h]], probs = 0.90, na.rm = TRUE, type = 7))
    if (!is.finite(h_p90)) {
      warning("Skipping model (h_p90 not finite): ", h, " x ", gstd)
      next
    }

    cbh <- crossbasis(
      sub[[h]],
      lag = 1,
      argvar = list(fun = "ns", knots = h_p90),
      arglag = list(fun = "integer")
    )
    colnames(cbh) <- paste0("h", seq_len(ncol(cbh)))

    cbg <- cbh
    cbg[] <- cbh[] * sub[[gstd]]
    colnames(cbg) <- paste0("g", seq_len(ncol(cbg)))

    model_dt <- cbind(
      sub[, .(deaths, day_of_season, year_f, dow_f, stratum_f)],
      as.data.frame(cbh),
      as.data.frame(cbg)
    )

    rhs <- paste(c(colnames(cbh), colnames(cbg), "dow_f", "ns(day_of_season, 4):year_f"), collapse = " + ")
    fml <- as.formula(paste0("deaths ~ ", rhs, " | stratum_f"))

    fit <- feglm(
      fml,
      data = model_dt,
      family = quasipoisson()
    )

    fits[[green_labels[[gstd]]]] <- fit
  }

  if (length(fits) == 0) {
    warning("No fitted models for heat metric: ", h)
    next
  }

  tex_file <- file.path(out_dir, paste0("etable_full_regression_", h, run_tag, ".tex"))
  tex_files <- c(tex_files, tex_file)

  # Keep all coefficients; FE are absorbed and reported in fit stats
  etable(
    fits,
    tex = TRUE,
    file = tex_file,
    replace = TRUE,
    title = paste0("Full regression coefficients (", heat_labels[[h]], ", run ", run_tag, ")"),
    label = paste0("tab:etable_full_reg_", h, run_tag),
    depvar = FALSE,
    headers = list("^Greenness metric" = c("GVI", "NDVI", "IMU total")),
    fixef_sizes = TRUE,
    fitstat = c("n"),
    digits = 4,
    signif.code = c("***" = 0.001, "**" = 0.01, "*" = 0.05)
  )

  # etable writes a table+tabular float; convert to longtable so content can
  # split across pages in the appendix instead of being clipped.
  to_longtable(tex_file)
}

master_file <- file.path(out_dir, paste0("etable_full_regression_tables_all", run_tag, ".tex"))
master_lines <- paste0("\\input{tables/appendix/full_regression_etable/", basename(tex_files), "}")
writeLines(c(master_lines, ""), master_file)

cat("\nSaved etable full regression tables in:\n - ", out_dir, "\n", sep = "")
cat("Master include file:\n - ", master_file, "\n", sep = "")
