library(readr)
library(dplyr)
library(forecast)
library(ggplot2)
library(tseries)
library(zoo)

source("forecast.R")
source("acf_pacf.R")
source("residuals_diagnostic.R")

# ---------------------------------------------------------------------------
# Load and filter to LRT Ampang Line, Jan 2019 - Jun 2026
#
# NOTE ON MCO HANDLING: this version starts from Jan 2019 (rather than the
# original Jan 2022 window that excluded MCO entirely). The Mar 2020 - Dec
# 2021 MCO lockdown period is NOT modeled on its raw crashed values -
# instead it is replaced with a known-intervention imputation (linear
# interpolation, then STL trend+seasonal reconstruction) so the model sees
# a smoothed, seasonally-consistent path through that window rather than
# the real shock. This trades off some realism for a longer training
# series without letting the MCO crash distort ACF/PACF identification or
# residual diagnostics. See mco_imputation.R for the standalone exhibit
# version of this same imputation logic.
# ---------------------------------------------------------------------------
df <- read_csv(
  "ridership_headline.csv",
  show_col_types = FALSE
)

stopifnot("rail_lrt_ampang" %in% names(df))

daily_ampang <- df %>%
  mutate(date = as.Date(date)) %>%
  filter(date >= as.Date("2019-01-01")) %>%
  arrange(date)

# Detect granularity instead of assuming daily: if every date already
# falls on the 1st of its month AND there's exactly one row per month,
# treat the source as already-monthly and skip aggregation. Otherwise
# aggregate daily rows to monthly totals as before.
is_already_monthly <- all(format(daily_ampang$date, "%d") == "01") &&
  !any(duplicated(format(daily_ampang$date, "%Y-%m")))

if (is_already_monthly) {
  cat("Source data is already monthly (one row per month) - skipping daily aggregation.\n")
  df_ampang <- daily_ampang %>%
    mutate(date = as.Date(format(date, "%Y-%m-01"))) %>%
    arrange(date)
} else {
  days_in_month_of <- function(month_start) {
    next_month <- seq(month_start, by = "month", length.out = 2)[2]
    as.integer(next_month - month_start)
  }
  
  df_ampang <- daily_ampang %>%
    mutate(date = as.Date(format(date, "%Y-%m-01"))) %>%
    group_by(date) %>%
    summarise(
      rail_lrt_ampang = sum(rail_lrt_ampang, na.rm = TRUE),
      n_days_observed = n(),
      .groups = "drop"
    ) %>%
    arrange(date) %>%
    mutate(n_days_expected = sapply(date, days_in_month_of))
  
  incomplete <- df_ampang %>% filter(n_days_observed != n_days_expected)
  
  if (nrow(incomplete) > 0) {
    cat("Warning: incomplete month(s) detected (partial daily data):\n")
    print(incomplete %>% select(date, n_days_observed, n_days_expected))
    
    last_row <- df_ampang[nrow(df_ampang), ]
    if (last_row$n_days_observed != last_row$n_days_expected && nrow(incomplete) == 1) {
      cat("Dropping the trailing partial month (", format(last_row$date, "%b %Y"),
          ") so it isn't treated as a full month's total.\n", sep = "")
      df_ampang <- df_ampang[-nrow(df_ampang), ]
    } else {
      stop("Incomplete month(s) found outside the trailing month - inspect ",
           "daily_ampang for gaps before continuing.")
    }
  }
}

n_months_expected <- length(
  seq(
    as.Date(format(min(df_ampang$date), "%Y-%m-01")),
    as.Date(format(max(df_ampang$date), "%Y-%m-01")),
    by = "month"
  )
)

if (nrow(df_ampang) != n_months_expected) {
  stop(sprintf(
    paste0(
      "Expected %d monthly rows between %s and %s, but got %d rows ",
      "after aggregating to monthly. Inspect df_ampang$date for gaps."
    ),
    n_months_expected,
    format(min(df_ampang$date), "%b %Y"),
    format(max(df_ampang$date), "%b %Y"),
    nrow(df_ampang)
  ))
}

cat("Aggregated to", nrow(df_ampang), "monthly rows:",
    format(min(df_ampang$date), "%b %Y"), "-",
    format(max(df_ampang$date), "%b %Y"), "\n")

# ---------------------------------------------------------------------------
# Known-intervention imputation for the MCO window (2020-03-01 to 2021-12-01)
# Step 1: linear interpolation across the masked MCO gap.
# Step 2: STL decomposition on the interpolated series; MCO-window values
#         are replaced with STL trend + seasonal (drops the remainder),
#         giving a smoother, seasonally-consistent "expected" path.
# ---------------------------------------------------------------------------
mco_start <- as.Date("2020-03-01")
mco_end   <- as.Date("2021-12-01")

df_ampang <- df_ampang %>%
  mutate(
    is_mco = date >= mco_start & date <= mco_end,
    rail_lrt_ampang_raw = rail_lrt_ampang,
    value_masked = ifelse(is_mco, NA, rail_lrt_ampang)
  ) %>%
  mutate(value_linear = na.approx(value_masked, x = date, na.rm = FALSE))

value_ts_for_stl <- ts(
  df_ampang$value_linear,
  start = c(as.integer(format(min(df_ampang$date), "%Y")),
            as.integer(format(min(df_ampang$date), "%m"))),
  frequency = 12
)

stl_fit <- stl(value_ts_for_stl, s.window = "periodic", robust = TRUE)
stl_components <- as.data.frame(stl_fit$time.series)

df_ampang <- df_ampang %>%
  mutate(
    stl_trend    = stl_components$trend,
    stl_seasonal = stl_components$seasonal,
    stl_reconstructed = stl_trend + stl_seasonal,
    rail_lrt_ampang = ifelse(is_mco, stl_reconstructed, rail_lrt_ampang_raw)
  )

cat("\nMCO window (", format(mco_start, "%b %Y"), " to ", format(mco_end, "%b %Y"),
    ") replaced with STL-imputed values:\n", sep = "")
print(
  df_ampang %>%
    filter(is_mco) %>%
    select(date, rail_lrt_ampang_raw, value_linear, rail_lrt_ampang)
)

# ---------------------------------------------------------------------------
ampang_vals  <- df_ampang$rail_lrt_ampang
ampang_dates <- df_ampang$date
n <- length(ampang_vals)

ampang_ts <- ts(
  ampang_vals,
  start = c(2019, 1),
  frequency = 12
)

# ---------------------------------------------------------------------------
# Train/test split: 11 months held out for testing (Aug 2025 - Jun 2026),
# remaining months (now including the imputed 2019-2021 window) for training
# ---------------------------------------------------------------------------
n_test  <- 11
n_train <- n - n_test

train_vals  <- ampang_vals[1:n_train]
train_dates <- ampang_dates[1:n_train]

test_vals  <- ampang_vals[(n_train + 1):n]
test_dates <- ampang_dates[(n_train + 1):n]

train_start_year  <- as.integer(format(train_dates[1], "%Y"))
train_start_month <- as.integer(format(train_dates[1], "%m"))

train_ts <- ts(
  train_vals,
  start = c(train_start_year, train_start_month),
  frequency = 12
)

test_start_year  <- as.integer(format(test_dates[1], "%Y"))
test_start_month <- as.integer(format(test_dates[1], "%m"))

test_ts <- ts(
  test_vals,
  start = c(test_start_year, test_start_month),
  frequency = 12
)

train <- list(date = train_dates, value = train_vals)
test  <- list(date = test_dates,  value = test_vals)

# ---------------------------------------------------------------------------
# Stationarity tests: ADF and KPSS on the level series vs. the differenced
# series. This is order-selection evidence for d = 1, run BEFORE any model
# is fit - same purpose as the ACF/PACF block below, different tool.
# Both tests are run together because they test OPPOSITE null hypotheses:
#   - ADF null:  series HAS a unit root (non-stationary)
#                -> large p-value  = fail to reject = non-stationary
#   - KPSS null: series IS stationary
#                -> small p-value  = reject         = non-stationary
# Agreement between the two is stronger evidence than either test alone.
# ---------------------------------------------------------------------------
diff_train_ts <- diff(train_ts)

cat("\n--- Stationarity tests (training series) ---\n")

adf_level <- adf.test(train_ts)
cat(sprintf(
  "ADF on level series:        statistic = %.3f, p-value = %.4f  (%s)\n",
  adf_level$statistic, adf_level$p.value,
  ifelse(adf_level$p.value > 0.05, "non-stationary - fail to reject unit root",
         "stationary - reject unit root")
))

adf_diff <- adf.test(diff_train_ts)
cat(sprintf(
  "ADF on differenced series:  statistic = %.3f, p-value = %.4f  (%s)\n",
  adf_diff$statistic, adf_diff$p.value,
  ifelse(adf_diff$p.value > 0.05, "non-stationary - fail to reject unit root",
         "stationary - reject unit root")
))

kpss_level <- kpss.test(train_ts, null = "Level")
cat(sprintf(
  "KPSS on level series:       statistic = %.3f, p-value = %.4f  (%s)\n",
  kpss_level$statistic, kpss_level$p.value,
  ifelse(kpss_level$p.value < 0.05, "non-stationary - reject null of stationarity",
         "stationary - fail to reject null")
))

kpss_diff <- kpss.test(diff_train_ts, null = "Level")
cat(sprintf(
  "KPSS on differenced series: statistic = %.3f, p-value = %.4f  (%s)\n",
  kpss_diff$statistic, kpss_diff$p.value,
  ifelse(kpss_diff$p.value < 0.05, "non-stationary - reject null of stationarity",
         "stationary - fail to reject null")
))

# ---------------------------------------------------------------------------
# ACF/PACF of the DIFFERENCED training series. This is order-selection
# evidence, computed BEFORE any model is fit, to justify d = 1 and the
# AR(1) choice - not the same thing as a post-fit residual check.
#
# IMPORTANT: with 2019 data now in the training set, re-examine this chart
# carefully - the correlogram used to justify ARIMA(1,1,0) was originally
# identified on the 2022-2026 window only. Confirm the same AR(1) cutoff
# pattern still holds before assuming order = c(1,1,0) is still correct.
# ---------------------------------------------------------------------------
generate_acf_pacf_chart(
  diff_train_ts,
  model_name = "Differenced Series (d = 1)",
  acf_path  = "acf.png",
  pacf_path = "pacf.png"
)

# ---------------------------------------------------------------------------
# Model: ARIMA(1,1,0), fit on TRAIN, evaluated on the 11-month hold-out.
# Compared against a naive baseline as the sanity check.
# ---------------------------------------------------------------------------
model_arima <- Arima(train_ts, order = c(1, 1, 0))
fc_arima    <- forecast(model_arima, h = length(test_vals))
mape_arima  <- mean(abs((test_vals - as.numeric(fc_arima$mean)) / test_vals)) * 100

cat("\n--- Hold-out evaluation (11-month test window) ---\n")
cat(sprintf(
  "%-24s AIC = %.1f   Hold-out MAPE = %.2f%%\n",
  "ARIMA(1,1,0)", model_arima$aic, mape_arima
))

naive_fc   <- naive(train_ts, h = length(test_vals))
naive_mape <- mean(abs((test_vals - as.numeric(naive_fc$mean)) / test_vals)) * 100
cat(sprintf(
  "%-24s              Hold-out MAPE = %.2f%%\n",
  "Naive baseline", naive_mape
))

print(summary(model_arima))

lb6 <- Box.test(
  residuals(model_arima),
  lag = 6,
  type = "Ljung-Box",
  fitdf = length(model_arima$coef)
)

lb12 <- Box.test(
  residuals(model_arima),
  lag = 12,
  type = "Ljung-Box",
  fitdf = length(model_arima$coef)
)

print(lb6)
print(lb12)

# ---------------------------------------------------------------------------
# Final forecast on FULL data (train + test), Jul-Dec 2026
# ---------------------------------------------------------------------------
final_model <- Arima(
  ampang_ts,
  order = c(1, 1, 0)
)

forecast_2026h2 <- forecast(
  final_model,
  h = 6
)

print(forecast_2026h2$mean)

last_month_start <- as.Date(
  format(max(ampang_dates), "%Y-%m-01")
)

forecast_dates_h2 <- seq(
  last_month_start,
  by = "month",
  length.out = 7
)[-1]

# Forecast chart: full history (train + test) plus the Jul-Dec 2026
# forecast, annotated with the hold-out MAPE for credibility.
generate_forecast_chart(
  train, test,
  forecast_vals  = as.numeric(forecast_2026h2$mean),
  forecast_dates = forecast_dates_h2,
  mape = mape_arima,
  output_path = "forecast.png",
  title_suffix = "",
  model_name = "ARIMA"
)

# Appendix - residual diagnostics on the final model
generate_residuals_chart(
  final_model,
  output_path = "residuals_diagnostic.png",
  model_name = "ARIMA"
)