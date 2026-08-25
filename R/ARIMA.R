# ARIMA.R
# ARIMA(1,1,0) analysis for LRT Ampang Line ridership (BMMS2094)
# R equivalent of ARIMA.py — run in Posit Cloud

# --- 0. Packages ------------------------------------------------------------
# Uncomment the line below the FIRST time you run this in a fresh Posit Cloud project:
# install.packages(c("forecast", "ggplot2", "scales", "patchwork", "readr", "dplyr"))

library(readr)
library(dplyr)
library(forecast)   # Arima(), auto.arima(), forecast()
library(ggplot2)

source("forecast.R")             # generate_forecast_chart()
source("acf_pacf.R")             # generate_acf_pacf_chart()
source("residuals_diagnostic.R") # generate_residuals_chart()

# --- 1. Load and filter to LRT Ampang Line, Jan 2022 - Jun 2026 -------------
df <- read_csv("ridership_headline.csv", show_col_types = FALSE)

df_ampang <- df %>%
  mutate(date = as.Date(date)) %>%
  filter(date >= as.Date("2022-01-01")) %>%
  arrange(date)

# rail_lrt_ampang column -> monthly ts
ampang_vals  <- df_ampang$rail_lrt_ampang
ampang_dates <- df_ampang$date

n <- length(ampang_vals)
ampang_ts <- ts(ampang_vals, start = c(2022, 1), frequency = 12)

# --- 2. Fixed 80/20 split: 43 train (79.6%) / 11 test (20.4%) of 54 months --
n_test <- 11
n_train <- n - n_test

train_vals  <- ampang_vals[1:n_train]
train_dates <- ampang_dates[1:n_train]
test_vals   <- ampang_vals[(n_train + 1):n]
test_dates  <- ampang_dates[(n_train + 1):n]

train_ts <- ts(train_vals, start = c(2022, 1), frequency = 12)
test_start_year  <- as.integer(format(test_dates[1], "%Y"))
test_start_month <- as.integer(format(test_dates[1], "%m"))
test_ts  <- ts(test_vals, start = c(test_start_year, test_start_month), frequency = 12)

train <- list(date = train_dates, value = train_vals)
test  <- list(date = test_dates,  value = test_vals)

# --- 3. AIC-based search (auto.arima ~ auto_arima) ---------------------------
auto_model <- auto.arima(train_ts, seasonal = FALSE, stepwise = FALSE,
                          approximation = FALSE, trace = TRUE)
print(summary(auto_model))

# --- 3b. Model comparison: (1,1,0) vs (1,1,1), drift vs no-drift --------------
# Fits all 4 combinations on the SAME train/test split so results are
# directly comparable. Use this table to decide both the order AND drift.

compare_models <- function(order, include_drift, label) {
  fit <- Arima(train_ts, order = order, include.drift = include_drift)
  fc  <- forecast(fit, h = length(test_vals))
  mape <- mean(abs((test_vals - as.numeric(fc$mean)) / test_vals)) * 100
  lb6  <- Box.test(residuals(fit), lag = 6,  type = "Ljung-Box", fitdf = length(fit$coef))
  lb12 <- Box.test(residuals(fit), lag = 12, type = "Ljung-Box", fitdf = length(fit$coef))
  data.frame(
    model        = label,
    order        = paste0("(", paste(order, collapse = ","), ")"),
    drift        = include_drift,
    AIC          = round(AIC(fit), 2),
    MAPE         = round(mape, 2),
    LB_p_lag6    = round(lb6$p.value, 4),
    LB_p_lag12   = round(lb12$p.value, 4),
    passes_LB    = lb6$p.value > 0.05 & lb12$p.value > 0.05
  )
}

comparison_table <- rbind(
  compare_models(c(1, 1, 0), FALSE, "ARIMA(1,1,0) no-drift"),
  compare_models(c(1, 1, 0), TRUE,  "ARIMA(1,1,0) with drift"),
  compare_models(c(1, 1, 1), FALSE, "ARIMA(1,1,1) no-drift"),
  compare_models(c(1, 1, 1), TRUE,  "ARIMA(1,1,1) with drift")
)

cat("\n=== Model comparison: order x drift ===\n")
print(comparison_table)
write_csv(comparison_table, "arima_model_comparison.csv")
cat("\nSaved comparison table to arima_model_comparison.csv\n")
cat("Decision rule: prefer the model with passes_LB = TRUE (residuals are white noise)\n")
cat("and the lowest MAPE among those that pass. Only fall back to AIC as a tiebreaker.\n\n")

# --- 3c. Seasonality diagnostic: is a hidden seasonal pattern causing the ----
# --- Ljung-Box failures seen across all 4 non-seasonal models? --------------

# 1) Residual ACF out to 24 lags for the best-AIC model (1,1,1 with drift),
#    looking specifically for spikes at lag 12 (and 6, 24) which would
#    indicate leftover monthly seasonality the non-seasonal model can't absorb
best_fit <- Arima(train_ts, order = c(1, 1, 1), include.drift = TRUE)

png("diagnostic_residual_acf_seasonal_check.png", width = 3000, height = 1500, res = 300)
Acf(residuals(best_fit), lag.max = 24,
    main = "Residual ACF (1,1,1 drift) - checking lag 6/12/24 for seasonality")
abline(v = c(6, 12, 18, 24), col = "red", lty = 2)
dev.off()
cat("\nSaved diagnostic_residual_acf_seasonal_check.png\n")
cat("Look for bars crossing the dashed blue CI lines AT the red dashed verticals\n")
cat("(lags 6, 12, 18, 24) - that pattern would indicate residual seasonality.\n\n")

# 2) Let auto.arima search WITH seasonality allowed, to see if it picks a
#    seasonal order - if it does, that's strong evidence seasonality is real
#    and not just noise (even though your assigned model stays non-seasonal)
cat("=== Seasonal auto.arima search (diagnostic only - not your final model) ===\n")
seasonal_check <- auto.arima(train_ts, seasonal = TRUE, stepwise = FALSE,
                              approximation = FALSE, trace = TRUE)
print(summary(seasonal_check))
cat("\nIf this picked a seasonal (P,D,Q)[12] order above (1,1,1)/(1,1,0),\n")
cat("that confirms seasonality is present in the LRT Ampang series.\n\n")

# --- 4. No-drift model: ARIMA(1,1,0), trend='n' -> include.drift = FALSE ----
model_nodrift <- Arima(train_ts, order = c(1, 1, 0), include.drift = FALSE)
fc_nodrift <- forecast(model_nodrift, h = length(test_vals))
mape_nodrift <- mean(abs((test_vals - as.numeric(fc_nodrift$mean)) / test_vals)) * 100
cat("No-drift hold-out MAPE:", mape_nodrift, "\n")
cat("No-drift Ljung-Box (lag 6):\n");  print(Box.test(residuals(model_nodrift), lag = 6,  type = "Ljung-Box"))
cat("No-drift Ljung-Box (lag 12):\n"); print(Box.test(residuals(model_nodrift), lag = 12, type = "Ljung-Box"))

# --- 5. With-drift model: ARIMA(1,1,0), trend='t' -> include.drift = TRUE ---
model_drift <- Arima(train_ts, order = c(1, 1, 0), include.drift = TRUE)
fc_drift <- forecast(model_drift, h = length(test_vals))
mape_drift <- mean(abs((test_vals - as.numeric(fc_drift$mean)) / test_vals)) * 100
cat("With-drift hold-out MAPE:", mape_drift, "\n")
cat("With-drift Ljung-Box (lag 6):\n");  print(Box.test(residuals(model_drift), lag = 6,  type = "Ljung-Box"))
cat("With-drift Ljung-Box (lag 12):\n"); print(Box.test(residuals(model_drift), lag = 12, type = "Ljung-Box"))

# --- 6. Final forecasts on full data, Jul-Dec 2026 ---------------------------
final_nodrift <- Arima(ampang_ts, order = c(1, 1, 0), include.drift = FALSE)
forecast_2026h2_nodrift <- forecast(final_nodrift, h = 6)
cat("No-drift final forecast:\n"); print(forecast_2026h2_nodrift$mean)

final_drift <- Arima(ampang_ts, order = c(1, 1, 0), include.drift = TRUE)
forecast_2026h2_drift <- forecast(final_drift, h = 6)
cat("With-drift final forecast:\n"); print(forecast_2026h2_drift$mean)

# Dates for the 6-month final forecast horizon (Jul-Dec 2026)
# (base R only, no lubridate needed): step from the 1st of the month after the last observed date
last_month_start <- as.Date(format(max(ampang_dates), "%Y-%m-01"))
forecast_dates_h2 <- seq(last_month_start, by = "month", length.out = 7)[-1]

# --- 7. Generate charts -------------------------------------------------------
generate_forecast_chart(train, test,
                         forecast_vals  = as.numeric(forecast_2026h2_nodrift$mean),
                         forecast_dates = forecast_dates_h2,
                         mape = mape_nodrift,
                         output_path = "arima_forecast_nodrift.png",
                         title_suffix = " (No Drift)")

generate_forecast_chart(train, test,
                         forecast_vals  = as.numeric(forecast_2026h2_drift$mean),
                         forecast_dates = forecast_dates_h2,
                         mape = mape_drift,
                         output_path = "arima_forecast_drift.png",
                         title_suffix = " (With Drift)")

# --- 8. Diagnostics for no-drift final model -----------------------------------
generate_acf_pacf_chart(residuals(final_nodrift),
                         model_name = "LRT Ampang ARIMA(1,1,0) No-Drift",
                         output_path = "arima_acf_pacf_nodrift.png")

generate_residuals_chart(final_nodrift,
                          output_path = "arima_residuals_nodrift.png",
                          model_name = "LRT Ampang ARIMA(1,1,0) No-Drift",
                          pacf_path = "arima_pacf_nodrift.png")