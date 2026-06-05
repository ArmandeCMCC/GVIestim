# Does greenness metric choice matter? Replication code

Replication code for the analysis of how urban greenness metric choice
(GVI, NDVI, IMU) affects the assessment of heat-mortality attenuation
in Paris.

This work extends the data and baseline model of
[Achebak et al. (2026)](https://doi.org/10.1038/s42949-025-00334-5)
(*npj Urban Sustainability* 6:29), which characterised contextual
factors modifying the heat-mortality association in Paris using a case
time series (CTS) design with distributed lag non-linear models (DLNM).

## Required R packages

```r
# Core data handling
install.packages(c("data.table", "lubridate"))

# Spatial
install.packages(c("sf", "terra", "exactextractr"))

# NetCDF (for UrbClim processing)
install.packages("ncdf4")

# Statistical modelling
install.packages(c("splines", "gnm", "dlnm"))

# fixest (for appendix full regression tables only)
install.packages("fixest")

# Plotting
install.packages(c("ggplot2", "patchwork"))
```

Tested with R >= 4.3. The `gnm` package (generalised nonlinear models)
provides `eliminate` for efficient conditional Poisson estimation in the
CTS design. The `dlnm` package provides `crossbasis` and `crosspred`
for distributed lag non-linear modelling.

## Required input data

| File | Description | Source |
|------|-------------|--------|
| `data.xlsx` | Daily mortality by arrondissement (2008-2017), residence-based | CepiDc / Inserm (request required) |
| `UrbClim/` | Hourly 100m x 100m air temperature NetCDF files for Paris (T2M, LST, WBGT) | UrbClim model (De Ridder et al. 2015) |
| `GHS_POP_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif` | GHS-POP 2020 population raster (3 arc-second, ~100m) | [GHSL](https://ghsl.jrc.ec.europa.eu/) |
| `arrondissements.shp` | Paris 20-arrondissement boundary polygons | [data.gouv.fr](https://www.data.gouv.fr/) |
| `gvi_356_cities.csv` | Street-level Green View Index points (LCZ-filtered) | Provided dataset |
| IMU shapefile or `IMU2022_O.geojson` | Ilots Morphologiques Urbains (urban islands) with vegetation fractions (`iv_haute`, `iv_basse`) | [IAU IdF](https://www.institutparisregion.fr/) |
| ESA WorldCover Sentinel-2 NDVI 10m | Annual peak-greenness (p90) composites for 2020 and 2021 (main paper) | [ESA WorldCover](https://registry.opendata.aws/esa-worldcover-vito-composites/) |
| NDVI Copernicus CLMS 300m COGs (optional) | NDVI 10-daily composites (2014-2017, JJAS) — used only by the legacy NDVI script for appendix robustness checks | [Copernicus Land](https://land.copernicus.eu/) |

## Execution order

All scripts are driven by environment variables for run-tag configuration.
The main paper results use tag `_rerun_rule1_sum_ghsheat_ghsgreen`.
The native (appendix) comparison uses tag `_nativeheat`.

### Step 1: Data preparation

These scripts can run independently of each other.

```bash
# 1a. Aggregate individual mortality records into daily arrondissement counts
Rscript build_paris_daily_deaths_fullgrid.R

# 1b. MAIN: Build NDVI arrondissement metrics from ESA WorldCover Sentinel-2
#     at 10m resolution (annual peak-greenness p90 composites, 2020-2021)
Rscript build_paris_arr_ndvi_metrics_esa10m.R

# 1b'. OPTIONAL (legacy / appendix sensitivity): rebuild NDVI metrics from the
#      older Copernicus CLMS 300m (2014-2017 JJAS) product. Not used by the
#      main paper, but populates the `ndvi_native_*` columns referenced in some
#      sensitivity outputs.
# Rscript build_paris_arr_ndvi_metrics.R

# 1c. Build all greenness metrics (GVI, NDVI, IMU — both weighted and native)
Rscript build_paris_arr_greenness_metrics.R

# 1d. MAIN: Build GHS-POP sum-weighted heat panel (method="sum" reprojection)
Rscript build_paris_panel_deaths_heat_ghspopweighted_sumproj.R

# 1e. APPENDIX: Build native (area-average) heat panel
Rscript build_paris_panel_deaths_heat_urbclim.R
```

### Step 2: CTS assembly

Merge deaths + heat panel + greenness into analysis-ready datasets.

```bash
# MAIN paper CTS (GHS-POP sum-weighted heat)
HEAT_SUPPORT=ghspopw_sum \
  RUN_TAG=_rerun_rule1_sum_ghsheat_ghsgreen \
  Rscript build_paris_cts_ready_jjas_sumproj.R

# APPENDIX native CTS (area-average heat)
HEAT_SUPPORT=native \
  RUN_TAG=_nativeheat \
  Rscript build_paris_cts_ready_jjas_sumproj.R
```

### Step 3: Baseline heat-mortality estimation

Fit heat-only CTS/DLNM models to recover MMT, p99, and IQR for each
heat metric.

```bash
# MAIN
RUN_TAG=_rerun_rule1_sum_ghsheat_ghsgreen \
  Rscript estimate_paris_baseline_cts_heat_metrics_rule1.R

# APPENDIX
RUN_TAG=_nativeheat \
  Rscript estimate_paris_baseline_cts_heat_metrics_rule1.R
```

### Step 4: Effect modification (heat x greenness interaction)

Fit interaction models for all 9 heat x all available greenness metrics.

```bash
# MAIN
RUN_TAG=_rerun_rule1_sum_ghsheat_ghsgreen \
  Rscript estimate_paris_effect_modification_cts_heat_green_rule1.R

# APPENDIX
RUN_TAG=_nativeheat \
  Rscript estimate_paris_effect_modification_cts_heat_green_rule1.R
```

### Step 5: Summarize effect modification results

Convert raw model outputs to paper-ready estimands with BH FDR
correction, separating GHS-weighted and native metric sets.

```bash
# MAIN
RUN_TAG=_rerun_rule1_sum_ghsheat_ghsgreen \
  Rscript summarize_paris_effectmod_results.R

# APPENDIX
RUN_TAG=_nativeheat \
  Rscript summarize_paris_effectmod_results.R
```

### Step 6: Attributable burden scenarios

Estimate deaths avoided under greening scenarios (max_all, p90_floor,
plus1iqr_cap).

```bash
# MAIN: 3 heat x 4 GHS-weighted green, all at MMT
RUN_TAG=_rerun_rule1_sum_ghsheat_ghsgreen \
  Rscript run_paris_attributable_scenarios_batch.R

# APPENDIX: 3 heat x 3 native green, at MMT (+ T2M mean at 22 deg C)
RUN_TAG=_nativeheat \
  Rscript run_paris_attributable_scenarios_batch_native.R
```

### Step 7: Summarize scenario comparison

Stack scenario summary CSVs into comparison tables.

```bash
RUN_TAG=_rerun_rule1_sum_ghsheat_ghsgreen \
  Rscript summarize_paris_scenario_comparison.R
```

### Step 8: Figures and tables

```bash
TAG=_rerun_rule1_sum_ghsheat_ghsgreen

# Forest plots + exposure-response curves (all 9x5 + 9x3 grids)
RUN_TAG=$TAG Rscript plot_paris_effectmod_results.R

# Compact curves, direct contrasts, burden heatmap (main figures)
RUN_TAG=$TAG Rscript plot_paris_main_insight_figures.R

# Weighted vs native comparison figure (appendix)
Rscript plot_paris_weighting_contrasts.R

# Appendix effect-modification LaTeX tables
RUN_TAG=$TAG Rscript export_paris_appendix_effectmod_tables.R

# Full regression coefficient tables (fixest/etable)
RUN_TAG=$TAG Rscript export_paris_full_regression_tables_etable.R
```

## Expected outputs

### Intermediate data (from Steps 1-2)

| File | Description |
|------|-------------|
| `paris_daily_deaths_arr_residence_fullgrid_2008_2017.csv` | Daily deaths by arrondissement (20 x 3652 days) |
| `paris_arr_greenness_metrics.csv` | All greenness metrics (20 arrondissements) |
| `paris_panel_residence_urbclim_daily_2008_2017_ghspopweighted_sumproj_JJAS.csv` | GHS-POP weighted heat panel (JJAS) |
| `paris_panel_residence_urbclim_daily_2008_2017_JJAS.csv` | Native heat panel (JJAS) |
| `paris_cts_ready_jjas_2008_2017{RUN_TAG}.csv` | Analysis-ready CTS dataset |

### Estimation outputs (from Steps 3-6)

| File pattern | Description |
|-------------|-------------|
| `paris_baseline_heat_metric_summary{RUN_TAG}.csv` | MMT, p99, IQR for each heat metric |
| `paris_baseline_heat_metric_curves{RUN_TAG}.csv` | Baseline dose-response curves |
| `paris_effectmod_heat_green_summary{RUN_TAG}.csv` | All effect-modification contrasts |
| `paris_effectmod_heat_green_curves{RUN_TAG}.csv` | Stratified exposure-response curves |
| `paris_effectmod_results_harmonized{RUN_TAG}.csv` | Paper-ready estimands with FDR |
| `paris_effectmod_primary_table_p99{RUN_TAG}.csv` | Primary estimand table (p10->p90 at p99) |
| `paris_attr_scenario_summary_*{RUN_TAG}.csv` | Burden scenario summaries |

### Figures (from Step 8)

| File | Description |
|------|-------------|
| `figs/paris_curves_effectmod_compact_main{RUN_TAG}.pdf` | Compact exposure-response curves (WBGT mean, T2M mean) |
| `figs/paris_effectmod_direct_contrasts{RUN_TAG}.pdf` | GVI vs NDVI / GVI vs IMU direct contrasts |
| `figs/paris_burden_key_twopanel{RUN_TAG}.pdf` | Avoided deaths heatmap |
| `figs/paris_weighting_contrasts{RUN_TAG}.pdf` | Pop-weighted vs native comparison |
| `figs/paris_forest_effectmod_*{RUN_TAG}.pdf` | Forest plots (one per heat metric) |

### Tables (from Step 8)

| File | Description |
|------|-------------|
| `tables/appendix/effectmod_tables_all{RUN_TAG}.tex` | Combined appendix effect-modification tables |
| `tables/appendix/full_regression_etable/etable_full_regression_*{RUN_TAG}.tex` | Full regression coefficient tables |

## Key modeling choices

All modeling choices follow [Achebak et al. (2026)](https://doi.org/10.1038/s42949-025-00334-5).

| Choice | Value | Rationale |
|--------|-------|-----------|
| **DLNM lag structure** | Integer function, lag 0-1 | Captures same-day + next-day acute heat effects on mortality. Short lags are standard in summer heat studies because heat-mortality effects are immediate. |
| **DLNM knot placement** | 1 internal knot at p90 | Concentrates flexibility in the upper tail where the dose-response relationship steepens, while keeping the curve smooth in the comfort range. |
| **Seasonal trend df** | `ns(day_of_season, df=4):year_f` | 4-df natural spline for smooth intra-seasonal control, interacted with year to allow varying summer shapes across years. |
| **Stratum definition** | Arrondissement x year x month | CTS design: conditions out all time-invariant confounders within each stratum via `gnm::eliminate`. |
| **Error distribution** | Quasi-Poisson | Accounts for overdispersion in daily death counts while maintaining log-link structure. |
| **Greenness standardisation** | `(G - median) / IQR` across N=20 arrondissements | Makes coefficients comparable across metrics with different natural scales (GVI 0-1, NDVI 0-1, IMU 0-100%). IQR is robust to non-normality with small N. |
| **Quantile type** | `type = 7` (R default) | Hyndman & Fan type 7; used throughout for consistency. |
| **Interpolation rule** | `approx(rule = 1)` | Returns NA outside the observed temperature range rather than extrapolating. Conservative choice to avoid out-of-sample inference. |
| **Population denominator** | 2,200,000 | Average Paris intra-muros population 2008-2017 (Achebak et al. 2026). |
| **FDR correction** | Benjamini-Hochberg | Applied across 45 primary comparisons (9 heat x 5 green metrics) to control expected false discovery rate. |
| **Heat population weighting** | `terra::project(method = "sum")` | Sum-preserving reprojection of GHS-POP to UrbClim grid conserves total population mass, unlike bilinear interpolation which can drop arrondissements with narrow coverage. |

## Environment variables reference

| Variable | Description | Default |
|----------|-------------|---------|
| `RUN_TAG` | Suffix appended to all output filenames | `""` |
| `CTS_TAG` | Tag for the input CTS file (defaults to `RUN_TAG`) | `= RUN_TAG` |
| `BASELINE_TAG` | Tag for the baseline summary file (defaults to `RUN_TAG`) | `= RUN_TAG` |
| `HEAT_SUPPORT` | Heat panel type: `ghspopw_sum` (default), `native`, `ghspopw`, `popw` | `ghspopw_sum` |
| `HEAT_FILTER` | Comma-separated heat metrics to run (estimation scripts) | all 9 |
| `SELECTED_HEAT` | Single heat metric (scenario runner) | `t2m_mean` |
| `SELECTED_GREEN_STD` | Single greenness metric (scenario runner) | `imu_veg_total_popw_ghs_std_iqr` |
| `HEAT_THRESHOLD_MODE` | Attribution threshold: `mmt` or numeric | `mmt` |

## Script inventory

| # | Script | Role |
|---|--------|------|
| 1 | `build_paris_daily_deaths_fullgrid.R` | Aggregate mortality records into daily arrondissement counts |
| 2 | `build_paris_panel_deaths_heat_ghspopweighted_sumproj.R` | **MAIN**: GHS-POP sum-weighted UrbClim heat panel |
| 3 | `build_paris_panel_deaths_heat_urbclim.R` | **APPENDIX**: Native (area-average) UrbClim heat panel |
| 4 | `build_paris_arr_greenness_metrics.R` | All greenness metrics (GVI, NDVI, IMU; weighted + native) |
| 5 | `build_paris_arr_ndvi_metrics_esa10m.R` | **MAIN**: NDVI from ESA WorldCover Sentinel-2 10m (annual p90 composites, 2020-2021). Replaces the older 300m source. |
| 5b | `build_paris_arr_ndvi_metrics.R` | **LEGACY** (optional, appendix sensitivity): NDVI from Copernicus CLMS 300m (JJAS 2014-2017). Populates the `ndvi_native_*` columns; not used by the published main analysis. |
| 6 | `build_paris_cts_ready_jjas_sumproj.R` | Merge deaths + heat + greenness into CTS-ready dataset |
| 7 | `estimate_paris_baseline_cts_heat_metrics_rule1.R` | Baseline CTS/DLNM (MMT, p99, dose-response curves) |
| 8 | `estimate_paris_effect_modification_cts_heat_green_rule1.R` | Heat x greenness interaction models |
| 9 | `estimate_paris_attributable_deaths_scenarios.R` | Single-spec attributable burden scenario runner |
| 10 | `run_paris_attributable_scenarios_batch.R` | **MAIN**: Batch burden scenarios (3 heat x 4 green) |
| 11 | `run_paris_attributable_scenarios_batch_native.R` | **APPENDIX**: Batch burden scenarios (3 heat x 3 green) |
| 12 | `summarize_paris_effectmod_results.R` | Paper-ready estimands + BH FDR correction |
| 13 | `summarize_paris_scenario_comparison.R` | Stack scenario summaries into comparison tables |
| 14 | `plot_paris_effectmod_results.R` | Forest plots + full exposure-response curve grids |
| 15 | `plot_paris_main_insight_figures.R` | Compact curves, direct contrasts, burden heatmap |
| 16 | `plot_paris_weighting_contrasts.R` | Pop-weighted vs native comparison figure |
| 17 | `export_paris_appendix_effectmod_tables.R` | Appendix LaTeX effect-modification tables |
| 18 | `export_paris_full_regression_tables_etable.R` | Full regression coefficient tables (fixest/etable) |
| 19 | `plot_paris_weighted_vs_native_comparison.R` | Generates figs/paris_weighted_vs_native_forest.png |
| 20 | `plot_paris_maps_descriptive.R` | Generates descriptive maps of Paris |
