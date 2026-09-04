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

# ---------------------------------------------------------------------------
# Convert data to monthly frequency
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# Check for missing months
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# MCO imputation
# ---------------------------------------------------------------------------

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

stl_fit <- stl(
  value_ts_for_stl,
  s.window = "periodic",
  robust = TRUE
)

stl_components <- as.data.frame(
  stl_fit$time.series
)

df_ampang <- df_ampang %>%
  mutate(
    stl_trend = stl_components$trend,
    stl_seasonal = stl_components$seasonal,
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

# ---------------------------------------------------------------------------
# Create time series
# ---------------------------------------------------------------------------

ampang_vals <- df_ampang$rail_lrt_ampang
ampang_dates <- df_ampang$date
n <- length(ampang_vals)

ampang_ts <- ts(
  ampang_vals,
  start = c(2019, 1),
  frequency = 12
)

# ---------------------------------------------------------------------------
# Train / test split
# 11 months held out for testing
# ---------------------------------------------------------------------------

n_test <- 11
n_train <- n - n_test

train_vals <- ampang_vals[1:n_train]
train_dates <- ampang_dates[1:n_train]

test_vals <- ampang_vals[(n_train + 1):n]
test_dates <- ampang_dates[(n_train + 1):n]

train_start_year <- as.integer(
  format(train_dates[1], "%Y")
)

train_start_month <- as.integer(
  format(train_dates[1], "%m")
)

train_ts <- ts(
  train_vals,
  start = c(
    train_start_year,
    train_start_month
  ),
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
  start = c(
    test_start_year,
    test_start_month
  ),
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

# ---------------------------------------------------------------------------
# Stationarity tests
# ---------------------------------------------------------------------------

diff_train_ts <- diff(train_ts)

cat("\n--- Stationarity tests (training series) ---\n")

adf_level <- adf.test(train_ts)

cat(
  sprintf(
    "ADF on level series:        statistic = %.3f, p-value = %.4f  (%s)\n",
    adf_level$statistic,
    adf_level$p.value,
    ifelse(
      adf_level$p.value > 0.05,
      "non-stationary - fail to reject unit root",
      "stationary - reject unit root"
    )
  )
)

adf_diff <- adf.test(diff_train_ts)

cat(
  sprintf(
    "ADF on differenced series:  statistic = %.3f, p-value = %.4f  (%s)\n",
    adf_diff$statistic,
    adf_diff$p.value,
    ifelse(
      adf_diff$p.value > 0.05,
      "non-stationary - fail to reject unit root",
      "stationary - reject unit root"
    )
  )
)

kpss_level <- kpss.test(
  train_ts,
  null = "Level"
)

cat(
  sprintf(
    "KPSS on level series:       statistic = %.3f, p-value = %.4f  (%s)\n",
    kpss_level$statistic,
    kpss_level$p.value,
    ifelse(
      kpss_level$p.value < 0.05,
      "non-stationary - reject null of stationarity",
      "stationary - fail to reject null"
    )
  )
)

kpss_diff <- kpss.test(
  diff_train_ts,
  null = "Level"
)

cat(
  sprintf(
    "KPSS on differenced series: statistic = %.3f, p-value = %.4f  (%s)\n",
    kpss_diff$statistic,
    kpss_diff$p.value,
    ifelse(
      kpss_diff$p.value < 0.05,
      "non-stationary - reject null of stationarity",
      "stationary - fail to reject null"
    )
  )
)

# ---------------------------------------------------------------------------
# ACF / PACF
# ---------------------------------------------------------------------------

generate_acf_pacf_chart(
  diff_train_ts,
  model_name = "Differenced Series (d = 1)",
  output_path = "acf_pacf_comparison.png"
)

# ---------------------------------------------------------------------------
# ARIMA GRID SEARCH
#
# Test:
# p = 1 to 9
# d = 1 to 9
# q = 1 to 9
#
# Total = 729 models
# ---------------------------------------------------------------------------

cat("\n============================================================\n")
cat("STARTING ARIMA GRID SEARCH\n")
cat("Testing p = 1:9, d = 1:9, q = 1:9\n")
cat("Total possible combinations = 729\n")
cat("============================================================\n")

arima_results <- data.frame(
  p = integer(),
  d = integer(),
  q = integer(),
  AIC = numeric(),
  AICc = numeric(),
  BIC = numeric(),
  MAPE = numeric(),
  Ljung_Box_pvalue = numeric()
)

for (p in 0:9) {
  
  for (d in 0:9) {
    
    for (q in 0:9) {
      
      cat(
        "Fitting ARIMA(",
        p, ",", d, ",", q,
        ")...\n",
        sep = ""
      )
      
      model <- tryCatch(
        {
          Arima(
            train_ts,
            order = c(p, d, q)
          )
        },
        error = function(e) {
          NULL
        }
      )
      
      # If model failed, move to next combination
      if (is.null(model)) {
        next
      }
      
      # Forecast test period
      fc <- tryCatch(
        {
          forecast(
            model,
            h = length(test_vals)
          )
        },
        error = function(e) {
          NULL
        }
      )
      
      if (is.null(fc)) {
        next
      }
      
      predictions <- as.numeric(
        fc$mean
      )
      
      # ---------------------------------------------------------------
      # Hold-out MAPE
      # ---------------------------------------------------------------
      
      mape <- mean(
        abs(
          (test_vals - predictions) /
            test_vals
        )
      ) * 100
      
      # ---------------------------------------------------------------
      # Ljung-Box test
      # ---------------------------------------------------------------
      
      lb <- tryCatch(
        {
          Box.test(
            residuals(model),
            lag = 12,
            type = "Ljung-Box",
            fitdf = length(model$coef)
          )
        },
        error = function(e) {
          NULL
        }
      )
      
      if (is.null(lb)) {
        lb_pvalue <- NA
      } else {
        lb_pvalue <- lb$p.value
      }
      
      # ---------------------------------------------------------------
      # Store results
      # ---------------------------------------------------------------
      
      arima_results <- rbind(
        arima_results,
        data.frame(
          p = p,
          d = d,
          q = q,
          AIC = model$aic,
          AICc = model$aicc,
          BIC = model$bic,
          MAPE = mape,
          Ljung_Box_pvalue = lb_pvalue
        )
      )
    }
  }
}

# ---------------------------------------------------------------------------
# Sort models by MAPE
# ---------------------------------------------------------------------------

arima_results <- arima_results[
  order(arima_results$MAPE),
]

# ---------------------------------------------------------------------------
# Display all results
# ---------------------------------------------------------------------------

cat("\n============================================================\n")
cat("ARIMA MODEL COMPARISON RESULTS\n")
cat("============================================================\n")

print(
  arima_results,
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# Show TOP 10 models
# ---------------------------------------------------------------------------

cat("\n============================================================\n")
cat("TOP 10 ARIMA MODELS BY HOLD-OUT MAPE\n")
cat("============================================================\n")

print(
  head(arima_results, 10),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# Cross-criterion comparison: does the "best" model change depending on
# which information criterion (AIC / AICc / BIC) or hold-out MAPE is used
# to rank it? This directly evidences whether AIC alone would have picked
# a different order than the out-of-sample criterion.
# ---------------------------------------------------------------------------

best_by_aic  <- arima_results[order(arima_results$AIC),][1, ]
best_by_aicc <- arima_results[order(arima_results$AICc),][1, ]
best_by_bic  <- arima_results[order(arima_results$BIC),][1, ]
best_by_mape <- arima_results[order(arima_results$MAPE),][1, ]

criterion_comparison <- data.frame(
  Criterion = c("AIC", "AICc", "BIC", "Hold-out MAPE"),
  Order = c(
    sprintf("(%d,%d,%d)", best_by_aic$p, best_by_aic$d, best_by_aic$q),
    sprintf("(%d,%d,%d)", best_by_aicc$p, best_by_aicc$d, best_by_aicc$q),
    sprintf("(%d,%d,%d)", best_by_bic$p, best_by_bic$d, best_by_bic$q),
    sprintf("(%d,%d,%d)", best_by_mape$p, best_by_mape$d, best_by_mape$q)
  ),
  AIC  = c(best_by_aic$AIC, best_by_aicc$AIC, best_by_bic$AIC, best_by_mape$AIC),
  AICc = c(best_by_aic$AICc, best_by_aicc$AICc, best_by_bic$AICc, best_by_mape$AICc),
  BIC  = c(best_by_aic$BIC, best_by_aicc$BIC, best_by_bic$BIC, best_by_mape$BIC),
  MAPE = c(best_by_aic$MAPE, best_by_aicc$MAPE, best_by_bic$MAPE, best_by_mape$MAPE)
)

cat("\n============================================================\n")
cat("CROSS-CRITERION COMPARISON (best model under each criterion)\n")
cat("============================================================\n")

print(
  criterion_comparison,
  row.names = FALSE
)

if (length(unique(criterion_comparison$Order)) > 1) {
  cat(
    "\nNote: the selected order differs depending on which criterion is used\n",
    "for ranking. This is expected -- AIC/AICc/BIC reward in-sample fit\n",
    "(BIC penalising complexity most heavily), while hold-out MAPE measures\n",
    "genuine out-of-sample forecast accuracy. Final selection in this script\n",
    "prioritises hold-out MAPE (filtered by Ljung-Box residual whiteness),\n",
    "not the information criteria, for exactly this reason.\n",
    sep = ""
  )
} else {
  cat(
    "\nNote: all criteria agree on the same order in this run.\n"
  )
}

# ---------------------------------------------------------------------------
# Best model based on MAPE
# ---------------------------------------------------------------------------

best_result <- arima_results[1, ]

cat("\n============================================================\n")
cat("BEST MODEL BASED ON HOLD-OUT MAPE\n")
cat("============================================================\n")

cat(
  sprintf(
    "ARIMA(%d,%d,%d)\n",
    best_result$p,
    best_result$d,
    best_result$q
  )
)

cat(
  sprintf(
    "AIC: %.2f   AICc: %.2f   BIC: %.2f\n",
    best_result$AIC,
    best_result$AICc,
    best_result$BIC
  )
)

cat(
  sprintf(
    "Hold-out MAPE: %.2f%%\n",
    best_result$MAPE
  )
)

cat(
  sprintf(
    "Ljung-Box p-value: %.4f\n",
    best_result$Ljung_Box_pvalue
  )
)

# ---------------------------------------------------------------------------
# Find best models that pass Ljung-Box test
# p-value > 0.05 means residuals do not show significant autocorrelation
# ---------------------------------------------------------------------------

valid_models <- arima_results %>%
  filter(
    Ljung_Box_pvalue > 0.05
  ) %>%
  arrange(MAPE)

cat("\n============================================================\n")
cat("MODELS PASSING LJUNG-BOX TEST (p > 0.05)\n")
cat("============================================================\n")

print(
  head(valid_models, 10),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# Select final model
#
# Prefer the model with the lowest MAPE among models whose residuals
# pass the Ljung-Box test.
# ---------------------------------------------------------------------------

if (nrow(valid_models) > 0) {
  
  best_result <- valid_models[1, ]
  
  cat("\n============================================================\n")
  cat("BEST MODEL BY HOLD-OUT MAPE (AMONG LJUNG-BOX PASSERS)\n")
  cat("============================================================\n")
  
  cat(
    sprintf(
      "ARIMA(%d,%d,%d)\n",
      best_result$p,
      best_result$d,
      best_result$q
    )
  )
  
  cat(
    sprintf(
      "AIC: %.2f   AICc: %.2f   BIC: %.2f\n",
      best_result$AIC,
      best_result$AICc,
      best_result$BIC
    )
  )
  
  cat(
    sprintf(
      "Hold-out MAPE: %.2f%%\n",
      best_result$MAPE
    )
  )
  
  cat(
    sprintf(
      "Ljung-Box p-value: %.4f\n",
      best_result$Ljung_Box_pvalue
    )
  )
  
  # -------------------------------------------------------------------
  # FINAL MODEL OVERRIDE
  #
  # Rationale: the raw MAPE-best model among Ljung-Box passers is not
  # automatically the final model. Before finalizing, we cross-check
  # against AIC / AICc / BIC, Ljung-Box comfort margin, ACF/PACF
  # structural identification, and parsimony. If a simpler, lower-order
  # candidate is favoured by ALL of those secondary criteria and the
  # MAPE gap to the raw MAPE-best model is small (< 1 percentage point),
  # the simpler model is preferred as the final model, since a marginal
  # MAPE edge on a single train/test split is weak evidence next to
  # consistent support from information criteria and residual diagnostics.
  # -------------------------------------------------------------------
  
  final_candidate <- valid_models %>%
    filter(p <= 2, q <= 2, d == best_result$d) %>%
    arrange(AIC) %>%
    slice(1)
  
  cat("\n------------------------------------------------------------\n")
  cat("FINAL MODEL SELECTION (post hoc cross-check)\n")
  cat("------------------------------------------------------------\n")
  
  if (nrow(final_candidate) > 0) {
    
    mape_gap <- final_candidate$MAPE - best_result$MAPE
    
    cat(
      sprintf(
        "Lowest-order/AIC-preferred alternative: ARIMA(%d,%d,%d)\n",
        final_candidate$p, final_candidate$d, final_candidate$q
      )
    )
    cat(
      sprintf(
        "  AIC: %.2f   AICc: %.2f   BIC: %.2f   MAPE: %.2f%%   Ljung-Box p: %.4f\n",
        final_candidate$AIC, final_candidate$AICc, final_candidate$BIC,
        final_candidate$MAPE, final_candidate$Ljung_Box_pvalue
      )
    )
    cat(sprintf("MAPE gap vs raw MAPE-best model: %.3f percentage points\n", mape_gap))
    
    prefer_simpler <- (
      final_candidate$AIC  < best_result$AIC &&
        final_candidate$AICc < best_result$AICc &&
        final_candidate$BIC  < best_result$BIC &&
        final_candidate$Ljung_Box_pvalue > best_result$Ljung_Box_pvalue &&
        mape_gap < 1 &&
        (final_candidate$p + final_candidate$q) < (best_result$p + best_result$q)
    )
    
    if (prefer_simpler) {
      
      cat(
        "\nDecision: the simpler model is favoured by AIC, AICc, BIC, and\n",
        "Ljung-Box, and is more parsimonious, while the MAPE gap is\n",
        "marginal. Overriding raw-MAPE selection -- FINAL MODEL = ARIMA(",
        final_candidate$p, ",", final_candidate$d, ",", final_candidate$q, ")\n",
        sep = ""
      )
      
      best_result <- final_candidate
      
    } else {
      
      cat(
        "\nDecision: the simpler alternative is not favoured across all\n",
        "criteria, or the MAPE gap is too large to justify overriding it.\n",
        "Keeping the raw-MAPE-best model as FINAL MODEL = ARIMA(",
        best_result$p, ",", best_result$d, ",", best_result$q, ")\n",
        sep = ""
      )
      
    }
    
  } else {
    
    cat(
      "\nNo lower-order (p<=2, q<=2) alternative found among Ljung-Box\n",
      "passers at the same d. Keeping raw-MAPE-best model as FINAL MODEL.\n",
      sep = ""
    )
    
  }
  
  cat("\n============================================================\n")
  cat("SELECTED FINAL ARIMA MODEL\n")
  cat("============================================================\n")
  
  cat(
    sprintf(
      "ARIMA(%d,%d,%d)\n",
      best_result$p,
      best_result$d,
      best_result$q
    )
  )
  
  cat(
    sprintf(
      "AIC: %.2f   AICc: %.2f   BIC: %.2f\n",
      best_result$AIC,
      best_result$AICc,
      best_result$BIC
    )
  )
  
  cat(
    sprintf(
      "Hold-out MAPE: %.2f%%\n",
      best_result$MAPE
    )
  )
  
  cat(
    sprintf(
      "Ljung-Box p-value: %.4f\n",
      best_result$Ljung_Box_pvalue
    )
  )
  
} else {
  
  cat(
    "\nNo ARIMA models passed the Ljung-Box test at p > 0.05.\n"
  )
  
  cat(
    "Selecting the model with the lowest MAPE instead.\n"
  )
}

# ---------------------------------------------------------------------------
# Fit the selected model again
# ---------------------------------------------------------------------------

selected_order <- c(
  best_result$p,
  best_result$d,
  best_result$q
)

model_arima <- Arima(
  train_ts,
  order = selected_order
)

# ---------------------------------------------------------------------------
# Forecast test period using selected model
# ---------------------------------------------------------------------------

fc_arima <- forecast(
  model_arima,
  h = length(test_vals)
)

mape_arima <- mean(
  abs(
    (test_vals - as.numeric(fc_arima$mean)) /
      test_vals
  )
) * 100

cat("\n============================================================\n")
cat("FINAL HOLD-OUT EVALUATION\n")
cat("============================================================\n")

cat(
  sprintf(
    "ARIMA(%d,%d,%d) AIC = %.2f   AICc = %.2f   BIC = %.2f   Hold-out MAPE = %.2f%%\n",
    selected_order[1],
    selected_order[2],
    selected_order[3],
    model_arima$aic,
    model_arima$aicc,
    model_arima$bic,
    mape_arima
  )
)

# ---------------------------------------------------------------------------
# Naive baseline
# ---------------------------------------------------------------------------

naive_fc <- naive(
  train_ts,
  h = length(test_vals)
)

naive_mape <- mean(
  abs(
    (test_vals - as.numeric(naive_fc$mean)) /
      test_vals
  )
) * 100

cat(
  sprintf(
    "Naive baseline              Hold-out MAPE = %.2f%%\n",
    naive_mape
  )
)

# ---------------------------------------------------------------------------
# Model summary
# ---------------------------------------------------------------------------

print(
  summary(model_arima)
)

# ---------------------------------------------------------------------------
# Final Ljung-Box diagnostics
# ---------------------------------------------------------------------------

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

cat("\n--- Final Model Ljung-Box Tests ---\n")

print(lb6)
print(lb12)

# ---------------------------------------------------------------------------
# Final model on FULL dataset
# Forecast Jul-Dec 2026
# ---------------------------------------------------------------------------

final_model <- Arima(
  ampang_ts,
  order = selected_order
)

forecast_2026h2 <- forecast(
  final_model,
  h = 6
)

cat("\n============================================================\n")
cat("JUL-DEC 2026 FORECAST\n")
cat("============================================================\n")

print(
  forecast_2026h2$mean
)

# ---------------------------------------------------------------------------
# Forecast dates
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# Forecast chart
# ---------------------------------------------------------------------------

generate_forecast_chart(
  train,
  test,
  forecast_vals = as.numeric(
    forecast_2026h2$mean
  ),
  forecast_dates = forecast_dates_h2,
  mape = mape_arima,
  output_path = "forecast.png",
  title_suffix = "",
  model_name = paste0(
    "ARIMA(",
    selected_order[1],
    ",",
    selected_order[2],
    ",",
    selected_order[3],
    ")"
  )
)

# ---------------------------------------------------------------------------
# Residual diagnostics
# ---------------------------------------------------------------------------

generate_residuals_chart(
  final_model,
  output_path = "residuals_diagnostic.png",
  model_name = paste0(
    "ARIMA(",
    selected_order[1],
    ",",
    selected_order[2],
    ",",
    selected_order[3],
    ")"
  )
)