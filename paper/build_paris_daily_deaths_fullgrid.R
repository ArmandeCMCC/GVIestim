# build_paris_daily_deaths_fullgrid.R: Script 1 
#
# Create the complete arrondissement-day mortality panel (place of
# residence) from data.xlsx. Produces a full 2008-01-01 to 2017-12-31
# grid with zeros for missing arrondissement-days.
#
# Shared input for both the main paper heat panel
# (build_paris_panel_deaths_heat_ghspopweighted_sumproj.R) and the
# native/unweighted appendix panel (build_paris_panel_deaths_heat_urbclim.R).
# Run this script first.
#
# Outputs:
#   - paris_daily_deaths_arr_residence_fullgrid_2008_2017.csv
#   - paris_daily_deaths_arr_residence_fullgrid_2008_2017_JJAS.csv
#   - paris_daily_deaths_arr_residence_fullgrid_2008_2017_byage.csv
#   - paris_daily_deaths_arr_residence_fullgrid_2008_2017_byage_JJAS.csv

library(data.table)
library(lubridate)
library(readxl)

base_dir  <- "/Users/armandeaboudrar-meda/Desktop/GVIestim/paper" # need to make it relative :   base_dir <- here::here()
xlsx_path <- file.path(base_dir, "data.xlsx")

paris_arr <- sprintf("751%02d", 1:20)

# 1) load raw mortality extract and keep required columns
dt <- as.data.table(read_excel(xlsx_path, sheet = 1))
stopifnot(all(c("date_dc","classe_age","arrondissement_de_residence","Effectifs") %in% names(dt)))

# 2) clean types, keep Paris arrondissements, and build standardised arr codes.
dt[, date := as.Date(date_dc)]
dt[, age_group := as.character(classe_age)]
dt[, arr_num := as.integer(arrondissement_de_residence)]
dt[, deaths := as.integer(Effectifs)]

dt <- dt[!is.na(date) & !is.na(arr_num)]
dt <- dt[arr_num >= 1 & arr_num <= 20]
dt[, arr := sprintf("751%02d", arr_num)]
dt <- dt[arr %chin% paris_arr]

# 3) build total-deaths daily series 
is_total <- grepl("^Tous", dt$age_group)

if (any(is_total)) {
  daily_total <- dt[is_total, .(deaths = sum(deaths, na.rm = TRUE)), by = .(arr, date)]
} else {
  daily_total <- dt[, .(deaths = sum(deaths, na.rm = TRUE)), by = .(arr, date)]
}
setorder(daily_total, arr, date)

# 4) build age-specific daily series
daily_age <- dt[!is_total, .(deaths = sum(deaths, na.rm = TRUE)), by = .(arr, date, age_group)]
setorder(daily_age, arr, date, age_group)

# 5) expand to full arrondissement-day grid (2008-01-01 to 2017-12-31)
start_date <- as.Date("2008-01-01")
end_date   <- as.Date("2017-12-31")
all_dates  <- seq(start_date, end_date, by = "day")

grid <- CJ(arr = paris_arr, date = all_dates)
daily_full <- merge(grid, daily_total, by = c("arr","date"), all.x = TRUE)
daily_full[is.na(deaths), deaths := 0L]
setorder(daily_full, arr, date)

stopifnot(nrow(daily_full) == length(paris_arr) * length(all_dates))  # 73060

daily_full_JJAS <- daily_full[month(date) %in% 6:9]

# 6) expand age-specific table to full arrondissement-day-age grid
age_levels <- sort(unique(daily_age$age_group))
grid_age <- CJ(arr = paris_arr, date = all_dates, age_group = age_levels)

daily_age_full <- merge(grid_age, daily_age, by = c("arr","date","age_group"), all.x = TRUE)
daily_age_full[is.na(deaths), deaths := 0L]
setorder(daily_age_full, arr, date, age_group)

daily_age_full_JJAS <- daily_age_full[month(date) %in% 6:9]

# 7) save outputs (full + JJAS and total + by-age).
out_full <- file.path(base_dir, "paris_daily_deaths_arr_residence_fullgrid_2008_2017.csv")
out_jjas <- file.path(base_dir, "paris_daily_deaths_arr_residence_fullgrid_2008_2017_JJAS.csv")
out_age  <- file.path(base_dir, "paris_daily_deaths_arr_residence_fullgrid_2008_2017_byage.csv")
out_age_jjas <- file.path(base_dir, "paris_daily_deaths_arr_residence_fullgrid_2008_2017_byage_JJAS.csv")

fwrite(daily_full, out_full)
fwrite(daily_full_JJAS, out_jjas)
fwrite(daily_age_full, out_age)
fwrite(daily_age_full_JJAS, out_age_jjas)

cat("Saved:\n -", out_full, "\n -", out_jjas, "\n -", out_age, "\n -", out_age_jjas, "\n")

# 8) sanity checks 
cat("\nSanity:\n")
cat("  Date range:", as.character(min(daily_full$date)), "to", as.character(max(daily_full$date)), "\n")
cat("  Unique arr:", uniqueN(daily_full$arr), "\n")
cat("  Total deaths:", sum(daily_full$deaths), "\n")
