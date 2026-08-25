library(readr)
library(dplyr)
library(forecast)
library(ggplot2)

source("forecast.R")
source("acf_pacf.R")
source("residuals_diagnostic.R")

df <- read_csv(
  "ridership_headline.csv",
  show_col_types = FALSE
)

df_ampang <- df %>%
  mutate(date = as.Date(date)) %>%
  filter(date >= as.Date("2022-01-01")) %>%
  arrange(date)

ampang_vals <- df_ampang$rail_lrt_ampang
ampang_dates <- df_ampang$date

n <- length(ampang_vals)

ampang_ts <- ts(
  ampang_vals,
  start = c(2022, 1),
  frequency = 12
)

n_test <- 11
n_train <- n - n_test

train_vals <- ampang_vals[1:n_train]
train_dates <- ampang_dates[1:n_train]

test_vals <- ampang_vals[(n_train + 1):n]
test_dates <- ampang_dates[(n_train + 1):n]

train_ts <- ts(
  train_vals,
  start = c(2022, 1),
  frequency = 12
)

test_start_year <- as.integer(
  format(test_dates[1], "%Y")
)

test_start_month <- as.integer(
  format(test_dates[1], "%m")
)

test_ts <- ts(
  test_vals,
  start = c(test_start_year, test_start_month),
  frequency = 12
)

train <- list(
  date = train_dates,
  value = train_vals
)

test <- list(
  date = test_dates,
  value = test_vals
)

auto_model <- auto.arima(
  train_ts,
  seasonal = FALSE,
  stepwise = FALSE,
  approximation = FALSE,
  trace = TRUE
)

print(summary(auto_model))

compare_models <- function(
    order,
    include_drift,
    label) {

  fit <- Arima(
    train_ts,
    order = order,
    include.drift = include_drift
  )

  fc <- forecast(
    fit,
    h = length(test_vals)
  )

  mape <- mean(
    abs(
      (test_vals - as.numeric(fc$mean)) /
        test_vals
    )
  ) * 100

  lb6 <- Box.test(
    residuals(fit),
    lag = 6,
    type = "Ljung-Box",
    fitdf = length(fit$coef)
  )

  lb12 <- Box.test(
    residuals(fit),
    lag = 12,
    type = "Ljung-Box",
    fitdf = length(fit$coef)
  )

  data.frame(
    model = label,
    order = paste0(
      "(",
      paste(order, collapse = ","),
      ")"
    ),
    drift = include_drift,
    AIC = round(AIC(fit), 2),
    MAPE = round(mape, 2),
    LB_p_lag6 = round(lb6$p.value, 4),
    LB_p_lag12 = round(lb12$p.value, 4),
    passes_LB =
      lb6$p.value > 0.05 &
      lb12$p.value > 0.05
  )
}

comparison_table <- rbind(
  compare_models(
    c(1, 1, 0),
    FALSE,
    "ARIMA(1,1,0) no-drift"
  ),

  compare_models(
    c(1, 1, 0),
    TRUE,
    "ARIMA(1,1,0) with drift"
  ),

  compare_models(
    c(1, 1, 1),
    FALSE,
    "ARIMA(1,1,1) no-drift"
  ),

  compare_models(
    c(1, 1, 1),
    TRUE,
    "ARIMA(1,1,1) with drift"
  )
)

print(comparison_table)

model_nodrift <- Arima(
  train_ts,
  order = c(1, 1, 0),
  include.drift = FALSE
)

fc_nodrift <- forecast(
  model_nodrift,
  h = length(test_vals)
)

mape_nodrift <- mean(
  abs(
    (test_vals - as.numeric(fc_nodrift$mean)) /
      test_vals
  )
) * 100

cat(
  "No-drift hold-out MAPE:",
  round(mape_nodrift, 2),
  "%\n"
)

print(
  Box.test(
    residuals(model_nodrift),
    lag = 6,
    type = "Ljung-Box"
  )
)

print(
  Box.test(
    residuals(model_nodrift),
    lag = 12,
    type = "Ljung-Box"
  )
)

model_drift <- Arima(
  train_ts,
  order = c(1, 1, 0),
  include.drift = TRUE
)

fc_drift <- forecast(
  model_drift,
  h = length(test_vals)
)

mape_drift <- mean(
  abs(
    (test_vals - as.numeric(fc_drift$mean)) /
      test_vals
  )
) * 100

cat(
  "With-drift hold-out MAPE:",
  round(mape_drift, 2),
  "%\n"
)

print(
  Box.test(
    residuals(model_drift),
    lag = 6,
    type = "Ljung-Box"
  )
)

print(
  Box.test(
    residuals(model_drift),
    lag = 12,
    type = "Ljung-Box"
  )
)

final_nodrift <- Arima(
  ampang_ts,
  order = c(1, 1, 0),
  include.drift = FALSE
)

forecast_2026h2_nodrift <- forecast(
  final_nodrift,
  h = 6
)

print(forecast_2026h2_nodrift$mean)

final_drift <- Arima(
  ampang_ts,
  order = c(1, 1, 0),
  include.drift = TRUE
)

forecast_2026h2_drift <- forecast(
  final_drift,
  h = 6
)

print(forecast_2026h2_drift$mean)

last_month_start <- as.Date(
  format(
    max(ampang_dates),
    "%Y-%m-01"
  )
)

forecast_dates_h2 <- seq(
  last_month_start,
  by = "month",
  length.out = 7
)[-1]

generate_acf_pacf_chart(
  residuals(final_nodrift),
  model_name = "LRT Ampang ARIMA(1,1,0) No-Drift",
  acf_path = "acf.png",
  pacf_path = "pacf.png"
)

generate_forecast_chart(
  train,
  test,
  forecast_vals = as.numeric(
    forecast_2026h2_nodrift$mean
  ),
  forecast_dates = forecast_dates_h2,
  mape = mape_nodrift,
  output_path = "forecast.png",
  title_suffix = " (No Drift)"
)

generate_residuals_chart(
  final_nodrift,
  output_path = "residuals.png",
  model_name = "LRT Ampang ARIMA(1,1,0) No-Drift"
)