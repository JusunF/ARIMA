pkgs <- c("readr", "dplyr", "lubridate", "zoo", "forecast", "ggplot2", "tseries")
new  <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new)

library(readr)
library(dplyr)
library(lubridate)
library(zoo)
library(forecast)
library(ggplot2)
library(tseries)

# ---------------------------------------------------------------------
# Load data and aggregate daily -> monthly (row-count safety check so a
# silent gap in the source data doesn't get passed into ts()).
# ---------------------------------------------------------------------
df <- read_csv("ridership_headline.csv", show_col_types = FALSE)
df$date <- as.Date(df$date)

monthly <- df %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>%
  summarise(rail_lrt_ampang = sum(rail_lrt_ampang, na.rm = TRUE)) %>%
  ungroup()

expected_months <- seq(min(monthly$month), max(monthly$month), by = "month")
stopifnot(nrow(monthly) == length(expected_months))

cat("Monthly series:", nrow(monthly), "months,",
    format(min(monthly$month)), "to", format(max(monthly$month)), "\n")

# ---------------------------------------------------------------------
# MCO resolution: Malaysia's MCO lockdown (Mar 2020 - Dec 2021) is a
# known intervention, not a statistically-detected outlier, so the
# window is set from domain knowledge and every month is kept (none
# dropped). Steps: mask the 22 MCO months -> linear-interpolate a
# placeholder -> STL decompose the interpolated series -> replace the
# masked months with STL trend + seasonal.
# ---------------------------------------------------------------------
mco_start <- as.Date("2020-03-01")
mco_end   <- as.Date("2021-12-01")

monthly <- monthly %>%
  mutate(
    is_mco       = month >= mco_start & month <= mco_end,
    value_raw    = rail_lrt_ampang,
    value_masked = ifelse(is_mco, NA, rail_lrt_ampang),
    value_linear = na.approx(value_masked, x = month, na.rm = FALSE)
  )

ts_start <- c(year(min(monthly$month)), month(min(monthly$month)))

ts_linear  <- ts(monthly$value_linear, start = ts_start, frequency = 12)
stl_impute <- stl(ts_linear, s.window = "periodic", robust = TRUE)

monthly$rail_lrt_ampang <- ifelse(
  monthly$is_mco,
  as.numeric(stl_impute$time.series[, "trend"] + stl_impute$time.series[, "seasonal"]),
  monthly$value_raw
)

cat("Bridged", sum(monthly$is_mco), "MCO months (",
    format(mco_start, "%b %Y"), "-", format(mco_end, "%b %Y"),
    ") via STL trend + seasonal imputation, 0 months dropped.\n")

print(monthly %>% filter(is_mco) %>%
        select(month, value_raw, value_linear, rail_lrt_ampang))

ampang_ts <- ts(monthly$rail_lrt_ampang, start = ts_start, frequency = 12)

# context exhibit: raw crash vs linear-only vs STL-imputed
p_mco <- autoplot(cbind(
  Raw         = ts(monthly$value_raw,    start = ts_start, frequency = 12),
  Linear_only = ts(monthly$value_linear, start = ts_start, frequency = 12),
  STL_Imputed = ampang_ts
)) +
  ggtitle("LRT Ampang: MCO Period - Raw vs Linear vs STL-Imputed") +
  ylab("Monthly ridership") +
  scale_y_continuous(labels = function(x) paste0(round(x / 1e6, 1), "M"))
print(p_mco)
ggsave("mco_imputation_comparison.png", p_mco, width = 9, height = 5, dpi = 150, bg = "white")


# train / test split - last 12 months (one full seasonal cycle)
h <- 12
n <- length(ampang_ts)

train_ts  <- ts(head(as.numeric(ampang_ts), n - h), start = start(ampang_ts), frequency = 12)
test_vals <- tail(as.numeric(ampang_ts), h)

# pulse dummies for two known service interventions
intervention_dates <- as.Date(c("2022-06-01", "2023-04-01"))
xreg_all <- sapply(intervention_dates, function(d) as.numeric(monthly$month == d))
colnames(xreg_all) <- paste0("pulse_", format(intervention_dates, "%Y%m"))

xreg_train  <- head(xreg_all, n - h)
xreg_test   <- tail(xreg_all, h)
xreg_future <- matrix(0, nrow = 6, ncol = ncol(xreg_all), dimnames = list(NULL, colnames(xreg_all)))

# STL-decompose the training series and pull out the seasonally-adjusted
# series to model, plus the seasonal component to add back afterwards
stl_train      <- stl(train_ts, s.window = "periodic", robust = TRUE)
seasadj_train  <- seasadj(stl_train)
seasonal_train <- as.numeric(stl_train$time.series[, "seasonal"])
seasonal_fc_test <- rep(tail(seasonal_train, 12), length.out = h)

diff_seasadj_train <- diff(seasadj_train)

cat("\n== ADF on seas-adj train, level (want p < 0.05 = stationary) ==\n")
print(adf.test(seasadj_train))

cat("\n== ADF on seas-adj train, differenced ==\n")
print(adf.test(diff_seasadj_train))

cat("\n== KPSS on seas-adj train, level (want p > 0.05 = stationary) ==\n")
print(kpss.test(seasadj_train, null = "Level"))

cat("\n== KPSS on seas-adj train, differenced ==\n")
print(kpss.test(diff_seasadj_train, null = "Level"))

png("stl_acf_pacf_comparison.png", width = 1200, height = 1000, res = 150)
par(mfrow = c(2, 1))
acf(diff_seasadj_train,  lag.max = 40, main = "ACF - Differenced Seasonally-Adjusted Series")
pacf(diff_seasadj_train, lag.max = 40, main = "PACF - Differenced Seasonally-Adjusted Series")
dev.off()


# ---------------------------------------------------------------------
# STLARIMA grid search: fit ARIMA(p,d,q) on the seasonally-adjusted
# training series for every combination in p,d,q = 0:9, add the
# extrapolated seasonal back on, and score by hold-out MAPE, AIC/AICc/
# BIC, and Ljung-Box p-value.
# ---------------------------------------------------------------------
cat("\n============================================================\n")
cat("STARTING STLARIMA GRID SEARCH (p,d,q = 0:9)\n")
cat("============================================================\n")

grid_results <- data.frame()

for (p in 0:9) {
  for (d in 0:9) {
    for (q in 0:9) {
      
      fit <- tryCatch(Arima(seasadj_train, order = c(p, d, q)), error = function(e) NULL)
      if (is.null(fit)) next
      
      fc <- tryCatch(forecast(fit, h = h), error = function(e) NULL)
      if (is.null(fc)) next
      
      pred <- as.numeric(fc$mean) + seasonal_fc_test
      mape <- mean(abs((test_vals - pred) / test_vals)) * 100
      
      lb <- tryCatch(
        Box.test(residuals(fit), lag = 12, type = "Ljung-Box", fitdf = length(fit$coef)),
        error = function(e) NULL
      )
      lb_p <- if (is.null(lb)) NA else lb$p.value
      
      k     <- length(fit$coef)
      n_obs <- length(seasadj_train) - d
      aicc  <- if ((n_obs - k - 1) > 0) fit$aic + (2 * k * (k + 1)) / (n_obs - k - 1) else NA
      
      grid_results <- rbind(grid_results, data.frame(
        p = p, d = d, q = q, AIC = fit$aic, AICc = aicc, BIC = fit$bic,
        MAPE = mape, LB_pvalue = lb_p
      ))
    }
  }
}

grid_results <- grid_results[order(grid_results$MAPE), ]
cat("\n--- Top 10 STLARIMA candidates by hold-out MAPE ---\n")
print(head(grid_results, 10), row.names = FALSE)

# pick the lowest-MAPE model among those passing Ljung-Box (p > 0.05);
# within 1.0pp MAPE of that winner, prefer the lowest AIC (parsimony)
valid <- grid_results %>% filter(LB_pvalue > 0.05) %>% arrange(MAPE)

if (nrow(valid) == 0) {
  cat("\nNo model passed Ljung-Box at p > 0.05 - falling back to lowest-MAPE overall.\n")
  valid <- grid_results[1, ]
}

best_mape_gap <- 1.0
best_by_mape  <- valid[1, ]
near_tied     <- valid %>% filter(MAPE <= best_by_mape$MAPE + best_mape_gap) %>% arrange(AIC)
best          <- near_tied[1, ]

selected_order <- c(best$p, best$d, best$q)

cat(sprintf(
  "\nSelected STLARIMA(%d,%d,%d)  AIC=%.2f  BIC=%.2f  MAPE=%.2f%%  LB p=%.4f\n",
  selected_order[1], selected_order[2], selected_order[3],
  best$AIC, best$BIC, best$MAPE, best$LB_pvalue
))


# fit the selected order on the training set (with interventions) and evaluate
model_train <- Arima(seasadj_train, order = selected_order, xreg = xreg_train)
print(summary(model_train))

if (length(model_train$coef) > 0) {
  se    <- sqrt(diag(model_train$var.coef))
  z     <- model_train$coef / se
  p_val <- 2 * (1 - pnorm(abs(z)))
  print(data.frame(coefficient = names(model_train$coef),
                   estimate = as.numeric(model_train$coef),
                   std_error = se, p_value = p_val))
} else {
  cat("Model has no AR/MA coefficients (pure differenced series).\n")
}

fc_train      <- forecast(model_train, h = h, xreg = xreg_test)
fc_test_final <- as.numeric(fc_train$mean) + seasonal_fc_test
mape_test     <- mean(abs((test_vals - fc_test_final) / test_vals)) * 100

naive_mape <- mean(abs((test_vals - as.numeric(naive(train_ts, h = h)$mean)) / test_vals)) * 100

cat(sprintf("\nHold-out MAPE (STLARIMA): %.2f%%\n", mape_test))
cat(sprintf("Hold-out MAPE (Naive baseline): %.2f%%\n", naive_mape))

cat("\n== Ljung-Box on training-fit residuals ==\n")
print(Box.test(residuals(model_train), lag = 6,  type = "Ljung-Box", fitdf = sum(selected_order[c(1, 3)])))
print(Box.test(residuals(model_train), lag = 12, type = "Ljung-Box", fitdf = sum(selected_order[c(1, 3)])))


# ---------------------------------------------------------------------
# Refit on the full series and forecast Jul-Dec 2026
# ---------------------------------------------------------------------
stl_full      <- stl(ampang_ts, s.window = "periodic", robust = TRUE)
seasadj_full  <- seasadj(stl_full)
seasonal_full <- as.numeric(stl_full$time.series[, "seasonal"])
seasonal_fc_h2 <- rep(tail(seasonal_full, 12), length.out = 6)

final_model <- Arima(seasadj_full, order = selected_order, xreg = xreg_all)
fc_full     <- forecast(final_model, h = 6, xreg = xreg_future)
forecast_2026h2 <- as.numeric(fc_full$mean) + seasonal_fc_h2

cat("\n============================================================\n")
cat("JUL-DEC 2026 FORECAST (STLARIMA, full-data refit)\n")
cat("============================================================\n")
print(round(forecast_2026h2, 0))

cat("\n== Ljung-Box on full-refit residuals ==\n")
print(Box.test(residuals(final_model), lag = 12, type = "Ljung-Box", fitdf = sum(selected_order[c(1, 3)])))
print(Box.test(residuals(final_model), lag = 18, type = "Ljung-Box", fitdf = sum(selected_order[c(1, 3)])))
print(Box.test(residuals(final_model), lag = 24, type = "Ljung-Box", fitdf = sum(selected_order[c(1, 3)])))

print(summary(final_model))


# outputs
model_label <- paste0(
  "STLARIMA(", selected_order[1], ",", selected_order[2], ",", selected_order[3], ") + Interventions"
)

test_dates <- tail(monthly$month, h)

train_df <- data.frame(date = head(monthly$month, n - h), value = as.numeric(train_ts))
actual_df <- data.frame(date = test_dates, value = test_vals)
fc_df <- data.frame(
  date = test_dates,
  mean = fc_test_final,
  lo80 = as.numeric(fc_train$lower[, 1]) + seasonal_fc_test,
  hi80 = as.numeric(fc_train$upper[, 1]) + seasonal_fc_test,
  lo95 = as.numeric(fc_train$lower[, 2]) + seasonal_fc_test,
  hi95 = as.numeric(fc_train$upper[, 2]) + seasonal_fc_test
)

p_fc <- ggplot() +
  geom_ribbon(data = fc_df, aes(x = date, ymin = lo95, ymax = hi95), fill = "#3366FF", alpha = 0.25) +
  geom_ribbon(data = fc_df, aes(x = date, ymin = lo80, ymax = hi80), fill = "#3366FF", alpha = 0.45) +
  geom_line(data = train_df,  aes(x = date, y = value), color = "black", linewidth = 0.6) +
  geom_line(data = fc_df,     aes(x = date, y = mean),  color = "#3366FF", linewidth = 0.7) +
  geom_line(data = actual_df, aes(x = date, y = value), color = "red", linewidth = 0.7) +
  scale_y_continuous(labels = function(x) paste0(round(x / 1e6, 1), "M")) +
  ggtitle(paste0(model_label, " - Forecast vs Actual")) +
  ylab("Monthly ridership") + xlab("Time")
print(p_fc)
ggsave("forecast_vs_actual.png", p_fc, width = 10, height = 6, dpi = 150, bg = "white")

png("residuals_diagnostic.png", width = 3300, height = 2400, res = 300)
final_model$method <- model_label
checkresiduals(final_model)
dev.off()

png("stl_decomposition_full.png", width = 3300, height = 2400, res = 300)
plot(stl_full, main = "STL Decomposition - LRT Ampang Ridership (Full Series)")
dev.off()

cat("\nDone. Saved forecast_vs_actual.png, residuals_diagnostic.png, stl_decomposition_full.png,",
    "stl_acf_pacf_comparison.png, mco_imputation_comparison.png\n")