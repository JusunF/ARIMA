library(readr)
library(dplyr)
library(forecast)
library(ggplot2)
library(tseries)

source("forecast.R")
source("acf_pacf.R")
source("residuals_diagnostic.R")

# ---------------------------------------------------------------------------
# Load and filter to LRT Ampang Line, Jan 2022 - Jun 2026
# (Jan 2022 start excludes the 2020-2021 MCO disruption period.)
# ---------------------------------------------------------------------------
df <- read_csv(
  "ridership_headline.csv",
  show_col_types = FALSE
)

stopifnot("rail_lrt_ampang" %in% names(df))

# The raw CSV is DAILY (confirmed: ~1,642 rows for a ~54-month window is
# ~30.4 rows/month). Aggregate to monthly totals - sum, not mean, since
# ridership is a count of riders per day and the model needs monthly volume.
daily_ampang <- df %>%
  mutate(date = as.Date(date)) %>%
  filter(date >= as.Date("2022-01-01")) %>%
  arrange(date)

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

# --- Safety check ------------------------------------------------------
# Everything below assumes ONE row per calendar month, each summing a
# COMPLETE set of days. Flag any month where observed days != calendar
# days in that month - most likely the most recent month if the data
# extract was pulled mid-month (a partial month would understate that
# month's total and distort both training and the hold-out MAPE).
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

ampang_vals  <- df_ampang$rail_lrt_ampang
ampang_dates <- df_ampang$date
n <- length(ampang_vals)

ampang_ts <- ts(
  ampang_vals,
  start = c(2022, 1),
  frequency = 12
)

# ---------------------------------------------------------------------------
# 43/11 split: 11 months held out for testing
# ---------------------------------------------------------------------------
n_test  <- 11
n_train <- n - n_test

train_vals  <- ampang_vals[1:n_train]
train_dates <- ampang_dates[1:n_train]

test_vals  <- ampang_vals[(n_train + 1):n]
test_dates <- ampang_dates[(n_train + 1):n]

train_ts <- ts(
  train_vals,
  start = c(2022, 1),
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