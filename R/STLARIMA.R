library(readr)
library(dplyr)
library(forecast)
library(ggplot2)
library(tseries)
library(zoo)

source("stl_decompose.R")
source("stl_arima_grid_search.R")
source("stl_arima_select.R")
source("forecast.R")
source("acf_pacf.R")
source("residuals_diagnostic.R")

df <- read_csv(
  "ridership_headline.csv",
  show_col_types = FALSE
)

stopifnot("rail_lrt_ampang" %in% names(df))

daily_ampang <- df %>%
  mutate(date = as.Date(date)) %>%
  filter(date >= as.Date("2019-01-01")) %>%
  arrange(date)

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
    mutate(
      n_days_expected = sapply(date, days_in_month_of)
    )
  
  incomplete <- df_ampang %>%
    filter(n_days_observed != n_days_expected)
  
  if (nrow(incomplete) > 0) {
    
    cat("Warning: incomplete month(s) detected:\n")
    
    print(
      incomplete %>%
        select(date, n_days_observed, n_days_expected)
    )
    
    last_row <- df_ampang[nrow(df_ampang), ]
    
    if (
      last_row$n_days_observed != last_row$n_days_expected &&
      nrow(incomplete) == 1
    ) {
      
      cat(
        "Dropping the trailing partial month (",
        format(last_row$date, "%b %Y"),
        ") so it isn't treated as a full month's total.\n",
        sep = ""
      )
      
      df_ampang <- df_ampang[-nrow(df_ampang), ]
      
    } else {
      
      stop(
        "Incomplete month(s) found outside the trailing month - inspect ",
        "daily_ampang for gaps before continuing."
      )
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
  
  stop(
    sprintf(
      paste0(
        "Expected %d monthly rows between %s and %s, but got %d rows ",
        "after aggregating to monthly. Inspect df_ampang$date for gaps."
      ),
      n_months_expected,
      format(min(df_ampang$date), "%b %Y"),
      format(max(df_ampang$date), "%b %Y"),
      nrow(df_ampang)
    )
  )
}

cat(
  "Aggregated to",
  nrow(df_ampang),
  "monthly rows:",
  format(min(df_ampang$date), "%b %Y"),
  "-",
  format(max(df_ampang$date), "%b %Y"),
  "\n"
)

mco_start <- as.Date("2020-03-01")
mco_end   <- as.Date("2021-12-01")

df_ampang <- df_ampang %>%
  mutate(
    is_mco = date >= mco_start & date <= mco_end,
    rail_lrt_ampang_raw = rail_lrt_ampang,
    value_masked = ifelse(is_mco, NA, rail_lrt_ampang)
  ) %>%
  mutate(
    value_linear = na.approx(
      value_masked,
      x = date,
      na.rm = FALSE
    )
  )

value_ts_for_stl <- ts(
  df_ampang$value_linear,
  start = c(
    as.integer(format(min(df_ampang$date), "%Y")),
    as.integer(format(min(df_ampang$date), "%m"))
  ),
  frequency = 12
)

mco_impute_decomp <- perform_stl_decomposition(value_ts_for_stl)
stl_components_impute <- as.data.frame(mco_impute_decomp$stl_fit$time.series)

df_ampang <- df_ampang %>%
  mutate(
    stl_trend = stl_components_impute$trend,
    stl_seasonal = stl_components_impute$seasonal,
    stl_reconstructed = stl_trend + stl_seasonal,
    rail_lrt_ampang = ifelse(
      is_mco,
      stl_reconstructed,
      rail_lrt_ampang_raw
    )
  )

cat(
  "\nMCO window (",
  format(mco_start, "%b %Y"),
  " to ",
  format(mco_end, "%b %Y"),
  ") replaced with STL-imputed values:\n",
  sep = ""
)

print(
  df_ampang %>%
    filter(is_mco) %>%
    select(
      date,
      rail_lrt_ampang_raw,
      value_linear,
      rail_lrt_ampang
    )
)

ampang_vals <- df_ampang$rail_lrt_ampang
ampang_dates <- df_ampang$date
n <- length(ampang_vals)

ampang_ts <- ts(
  ampang_vals,
  start = c(2019, 1),
  frequency = 12
)

n_test <- 12
n_train <- n - n_test

train_vals <- ampang_vals[1:n_train]
train_dates <- ampang_dates[1:n_train]

test_vals <- ampang_vals[(n_train + 1):n]
test_dates <- ampang_dates[(n_train + 1):n]

train_start_year <- as.integer(format(train_dates[1], "%Y"))
train_start_month <- as.integer(format(train_dates[1], "%m"))

train_ts <- ts(
  train_vals,
  start = c(train_start_year, train_start_month),
  frequency = 12
)

test_start_year <- as.integer(format(test_dates[1], "%Y"))
test_start_month <- as.integer(format(test_dates[1], "%m"))

test_ts <- ts(
  test_vals,
  start = c(test_start_year, test_start_month),
  frequency = 12
)

train <- list(date = train_dates, value = train_vals)
test  <- list(date = test_dates,  value = test_vals)

intervention_dates <- as.Date(c("2022-06-01", "2023-04-01"))

make_intervention_xreg <- function(dates_vec) {
  m <- sapply(intervention_dates, function(d) as.numeric(dates_vec == d))
  colnames(m) <- paste0("pulse_", format(intervention_dates, "%Y%m"))
  m
}

xreg_train  <- make_intervention_xreg(train_dates)
xreg_test   <- make_intervention_xreg(test_dates)
xreg_full   <- make_intervention_xreg(ampang_dates)
xreg_future <- matrix(0, nrow = 6, ncol = 2,
                      dimnames = list(NULL, colnames(xreg_train)))

train_decomp <- perform_stl_decomposition(train_ts)

seasadj_train_ts <- train_decomp$seasadj
seasonal_fc_test <- extrapolate_seasonal(train_decomp$seasonal, h = length(test_vals))

cat("\n--- STL decomposition (training series) ---\n")
cat(sprintf(
  "Seasonal component range: %.0f to %.0f\n",
  min(train_decomp$seasonal), max(train_decomp$seasonal)
))

diff_seasadj_train_ts <- diff(seasadj_train_ts)

cat("\n--- Stationarity tests (seasonally-adjusted training series) ---\n")

adf_level <- adf.test(seasadj_train_ts)
cat(sprintf(
  "ADF on level seas-adj series:       statistic = %.3f, p-value = %.4f  (%s)\n",
  adf_level$statistic, adf_level$p.value,
  ifelse(adf_level$p.value > 0.05, "non-stationary", "stationary")
))

adf_diff <- adf.test(diff_seasadj_train_ts)
cat(sprintf(
  "ADF on differenced seas-adj series: statistic = %.3f, p-value = %.4f  (%s)\n",
  adf_diff$statistic, adf_diff$p.value,
  ifelse(adf_diff$p.value > 0.05, "non-stationary", "stationary")
))

kpss_level <- kpss.test(seasadj_train_ts, null = "Level")
cat(sprintf(
  "KPSS on level seas-adj series:       statistic = %.3f, p-value = %.4f  (%s)\n",
  kpss_level$statistic, kpss_level$p.value,
  ifelse(kpss_level$p.value < 0.05, "non-stationary", "stationary")
))

kpss_diff <- kpss.test(diff_seasadj_train_ts, null = "Level")
cat(sprintf(
  "KPSS on differenced seas-adj series: statistic = %.3f, p-value = %.4f  (%s)\n",
  kpss_diff$statistic, kpss_diff$p.value,
  ifelse(kpss_diff$p.value < 0.05, "non-stationary", "stationary")
))

generate_acf_pacf_chart(
  diff_seasadj_train_ts,
  model_name = "Differenced Seasonally-Adjusted Series (d = 1)",
  output_path = "stl_acf_pacf_comparison.png"
)

cat("\n============================================================\n")
cat("STARTING STLARIMA GRID SEARCH\n")
cat("Testing p = 0:9, d = 0:9, q = 0:9 on the seasonally-adjusted series\n")
cat("============================================================\n")

stl_arima_results <- run_stl_arima_grid_search(
  seasadj_train_ts,
  test_vals,
  seasonal_fc_test,
  p_range = 0:9,
  d_range = 0:9,
  q_range = 0:9
)

cat("\n============================================================\n")
cat("TOP 10 STLARIMA MODELS BY HOLD-OUT MAPE\n")
cat("============================================================\n")
print(head(stl_arima_results, 10), row.names = FALSE)

best_result <- select_stl_arima_model(stl_arima_results, mape_gap_threshold = 1.0)

selected_order <- c(best_result$p, best_result$d, best_result$q)

model_stl_arima <- Arima(seasadj_train_ts, order = selected_order, xreg = xreg_train)

cat("\n--- Coefficient significance (selected model, training fit) ---\n")
print(summary(model_stl_arima))
if (length(model_stl_arima$coef) > 0) {
  coef_se <- sqrt(diag(model_stl_arima$var.coef))
  z_vals <- model_stl_arima$coef / coef_se
  p_vals <- 2 * (1 - pnorm(abs(z_vals)))
  print(data.frame(
    coefficient = names(model_stl_arima$coef),
    estimate = as.numeric(model_stl_arima$coef),
    std_error = coef_se,
    p_value = p_vals
  ))
} else {
  cat("Model has no AR/MA coefficients (pure differenced series).\n")
}

fc_stl_arima <- forecast(model_stl_arima, h = length(test_vals), xreg = xreg_test)
fc_stl_arima_final <- as.numeric(fc_stl_arima$mean) + seasonal_fc_test

mape_stl_arima <- mean(abs((test_vals - fc_stl_arima_final) / test_vals)) * 100

cat("\n============================================================\n")
cat("FINAL HOLD-OUT EVALUATION\n")
cat("============================================================\n")
cat(sprintf(
  "STLARIMA(%d,%d,%d) AIC = %.2f   Hold-out MAPE = %.2f%%\n",
  selected_order[1], selected_order[2], selected_order[3],
  model_stl_arima$aic, mape_stl_arima
))

naive_fc <- naive(train_ts, h = length(test_vals))
naive_mape <- mean(abs((test_vals - as.numeric(naive_fc$mean)) / test_vals)) * 100

cat(sprintf(
  "Naive baseline              Hold-out MAPE = %.2f%%\n",
  naive_mape
))

lb6 <- Box.test(residuals(model_stl_arima), lag = 6, type = "Ljung-Box",
                fitdf = sum(selected_order[c(1, 3)]))
lb12 <- Box.test(residuals(model_stl_arima), lag = 12, type = "Ljung-Box",
                 fitdf = sum(selected_order[c(1, 3)]))

cat("\n--- Final Model Ljung-Box Tests (training fit) ---\n")
print(lb6)
print(lb12)

full_decomp <- perform_stl_decomposition(ampang_ts)
seasadj_full_ts <- full_decomp$seasadj
seasonal_fc_h2 <- extrapolate_seasonal(full_decomp$seasonal, h = 6)

final_model <- Arima(seasadj_full_ts, order = selected_order, xreg = xreg_full)

forecast_seasadj_h2 <- forecast(final_model, h = 6, xreg = xreg_future)
forecast_2026h2_mean <- as.numeric(forecast_seasadj_h2$mean) + seasonal_fc_h2

cat("\n============================================================\n")
cat("JUL-DEC 2026 FORECAST (STLARIMA, full-data refit)\n")
cat("============================================================\n")
print(round(forecast_2026h2_mean, 0))

cat("\n--- Full-data refit Ljung-Box check ---\n")
lb_full_12 <- Box.test(residuals(final_model), lag = 12, type = "Ljung-Box",
                       fitdf = sum(selected_order[c(1, 3)]))
lb_full_18 <- Box.test(residuals(final_model), lag = 18, type = "Ljung-Box",
                       fitdf = sum(selected_order[c(1, 3)]))
lb_full_24 <- Box.test(residuals(final_model), lag = 24, type = "Ljung-Box",
                       fitdf = sum(selected_order[c(1, 3)]))
print(lb_full_12)
print(lb_full_18)
print(lb_full_24)

last_month_start <- as.Date(format(max(ampang_dates), "%Y-%m-01"))

forecast_dates_h2 <- seq(
  last_month_start,
  by = "month",
  length.out = 7
)[-1]

model_label <- paste0(
  "STLARIMA(", selected_order[1], ",", selected_order[2], ",", selected_order[3], ")",
  " + Interventions"
)

generate_forecast_chart(
  train,
  test,
  forecast_vals = forecast_2026h2_mean,
  forecast_dates = forecast_dates_h2,
  mape = mape_stl_arima,
  output_path = "forecast.png",
  title_suffix = "",
  model_name = model_label
)

generate_residuals_chart(
  final_model,
  output_path = "residuals_diagnostic.png",
  model_name = model_label
)

png("stl_decomposition_full.png", width = 3300, height = 2400, res = 300)
plot(full_decomp$stl_fit, main = "STL Decomposition - LRT Ampang Ridership (Full Series)")
dev.off()
cat("\nSTL decomposition chart saved as: stl_decomposition_full.png\n")

print(summary(final_model))
