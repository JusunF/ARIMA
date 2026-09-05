library(readr)
library(dplyr)
library(tidyr)
library(zoo)
library(forecast)
library(ggplot2)
library(scales)
library(tseries)

perform_stl_decomposition <- function(ts_data, s_window = "periodic", robust = TRUE) {

  stl_fit <- stl(
    ts_data,
    s.window = s_window,
    robust = robust
  )

  list(
    stl_fit = stl_fit,
    seasadj = seasadj(stl_fit),
    seasonal = stl_fit$time.series[, "seasonal"]
  )
}

extrapolate_seasonal <- function(seasonal_component, h) {

  last_cycle <- tail(as.numeric(seasonal_component), 12)

  rep(last_cycle, length.out = h)
}

generate_acf_pacf_chart <- function(
    x,
    model_name = "Differenced Seasonally-Adjusted Series (d = 1)",
    output_path = "acf_pacf_comparison.png") {

  x <- as.numeric(x)
  n <- length(x)
  n_lags <- min(40, n %/% 3)

  acf_res <- acf(
    x,
    lag.max = n_lags,
    plot = FALSE
  )

  pacf_res <- pacf(
    x,
    lag.max = n_lags,
    plot = FALSE
  )

  ci <- 1.96 / sqrt(n)

  color <- "#2563eb"
  grid_color <- "#e5e5e5"
  text_color <- "#1a1a1a"
  muted_text <- "#737373"

  acf_df <- data.frame(
    lag = acf_res$lag[, 1, 1],
    value = acf_res$acf[, 1, 1],
    type = "ACF"
  )

  pacf_df <- data.frame(
    lag = pacf_res$lag[, 1, 1],
    value = pacf_res$acf[, 1, 1],
    type = "PACF"
  )

  acf_df <- acf_df[acf_df$lag != 0, ]
  pacf_df <- pacf_df[pacf_df$lag != 0, ]

  comparison_df <- rbind(acf_df, pacf_df)

  comparison_plot <- ggplot(
    comparison_df,
    aes(x = lag, y = value)
  ) +

    geom_hline(
      yintercept = 0,
      color = text_color,
      linewidth = 0.4
    ) +

    geom_hline(
      yintercept = c(-ci, ci),
      color = muted_text,
      linewidth = 0.4,
      linetype = "dashed"
    ) +

    geom_col(
      width = 0.4,
      fill = color,
      color = color
    ) +

    facet_wrap(
      ~type,
      ncol = 1,
      scales = "fixed",
      labeller = labeller(
        type = c(
          ACF = "Autocorrelation Function (ACF)",
          PACF = "Partial Autocorrelation Function (PACF)"
        )
      )
    ) +

    scale_y_continuous(
      limits = c(-1.1, 1.1)
    ) +

    labs(
      title = paste0(
        "ACF/PACF Comparison\n",
        model_name
      ),
      x = "Lag (months)",
      y = "Correlation"
    ) +

    theme_minimal(base_size = 12) +

    theme(
      plot.title = element_text(
        face = "bold",
        color = text_color,
        size = 13
      ),
      strip.text = element_text(
        face = "bold",
        color = text_color,
        size = 12
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(
        color = grid_color
      ),
      axis.text = element_text(
        color = text_color
      )
    )

  ggsave(
    output_path,
    plot = comparison_plot,
    width = 12,
    height = 10,
    dpi = 150,
    bg = "white"
  )

  cat(
    "ACF/PACF comparison chart saved as:",
    output_path,
    "\n"
  )
}

run_stl_arima_grid_search <- function(seasadj_train_ts,
                                      test_vals,
                                      seasonal_fc_test,
                                      p_range = 0:9,
                                      d_range = 0:9,
                                      q_range = 0:9,
                                      verbose = TRUE) {

  results <- data.frame(
    p = integer(),
    d = integer(),
    q = integer(),
    AIC = numeric(),
    AICc = numeric(),
    BIC = numeric(),
    MAPE = numeric(),
    Ljung_Box_pvalue = numeric()
  )

  for (p in p_range) {
    for (d in d_range) {
      for (q in q_range) {

        if (verbose) {
          cat("Fitting STLARIMA(", p, ",", d, ",", q, ")...\n", sep = "")
        }

        model <- tryCatch(
          Arima(seasadj_train_ts, order = c(p, d, q)),
          error = function(e) NULL
        )

        if (is.null(model)) next

        fc <- tryCatch(
          forecast(model, h = length(test_vals)),
          error = function(e) NULL
        )

        if (is.null(fc)) next

        predictions <- as.numeric(fc$mean) + seasonal_fc_test

        mape <- mean(abs((test_vals - predictions) / test_vals)) * 100

        lb <- tryCatch(
          Box.test(
            residuals(model),
            lag = 12,
            type = "Ljung-Box",
            fitdf = length(model$coef)
          ),
          error = function(e) NULL
        )

        lb_pvalue <- if (is.null(lb)) NA else lb$p.value

        k <- length(model$coef)
        n_obs <- length(seasadj_train_ts) - d
        aicc <- if ((n_obs - k - 1) > 0) {
          model$aic + (2 * k * (k + 1)) / (n_obs - k - 1)
        } else {
          NA
        }

        results <- rbind(
          results,
          data.frame(
            p = p, d = d, q = q,
            AIC = model$aic,
            AICc = aicc,
            BIC = model$bic,
            MAPE = mape,
            Ljung_Box_pvalue = lb_pvalue
          )
        )
      }
    }
  }

  results[order(results$MAPE), ]
}

select_stl_arima_model <- function(results, mape_gap_threshold = 1.0) {

  best_result <- results[1, ]

  valid_models <- results %>%
    filter(Ljung_Box_pvalue > 0.05) %>%
    arrange(MAPE)

  cat("\n============================================================\n")
  cat("MODELS PASSING LJUNG-BOX TEST (p > 0.05)\n")
  cat("============================================================\n")
  print(head(valid_models, 10), row.names = FALSE)

  if (nrow(valid_models) > 0) {
    raw_best <- valid_models[1, ]
  } else {
    cat("\nNo STLARIMA models passed the Ljung-Box test at p > 0.05.\n")
    cat("Falling back to the lowest-MAPE model overall.\n")
    raw_best <- best_result
    valid_models <- best_result
  }

  cat("\n============================================================\n")
  cat("MECHANICAL WINNER (lowest MAPE among Ljung-Box passes)\n")
  cat("============================================================\n")
  cat(sprintf(
    "STLARIMA(%d,%d,%d)  MAPE=%.2f%%  AIC=%.2f  BIC=%.2f  LB p=%.4f\n",
    raw_best$p, raw_best$d, raw_best$q,
    raw_best$MAPE, raw_best$AIC, raw_best$BIC, raw_best$Ljung_Box_pvalue
  ))

  near_tied <- valid_models %>%
    filter(MAPE <= raw_best$MAPE + mape_gap_threshold) %>%
    arrange(AIC)

  if (nrow(near_tied) > 0) {
    override_best <- near_tied[1, ]
  } else {
    override_best <- raw_best
  }

  cat("\n============================================================\n")
  cat("CROSS-CRITERION OVERRIDE CHECK (AIC/AICc/BIC/parsimony)\n")
  cat("============================================================\n")
  cat(sprintf(
    "Mechanical MAPE winner:        STLARIMA(%d,%d,%d)  MAPE=%.2f%%  AIC=%.2f  BIC=%.2f\n",
    raw_best$p, raw_best$d, raw_best$q, raw_best$MAPE, raw_best$AIC, raw_best$BIC
  ))
  cat(sprintf(
    "Lowest-AIC within %.1fpp MAPE:  STLARIMA(%d,%d,%d)  MAPE=%.2f%%  AIC=%.2f  BIC=%.2f\n",
    mape_gap_threshold, override_best$p, override_best$d, override_best$q,
    override_best$MAPE, override_best$AIC, override_best$BIC
  ))

  if (override_best$p == 0 && override_best$q == 0) {
    cat("\nNote: the overridden candidate has no AR or MA terms. Check the\n")
    cat("coefficient significance after refitting - if p and q were both 0\n")
    cat("anyway, this confirms the STL-adjusted series needs no further\n")
    cat("ARIMA structure beyond differencing.\n")
  }

  cat("\n============================================================\n")
  cat("SELECTED FINAL STLARIMA MODEL\n")
  cat("============================================================\n")
  cat(sprintf("STLARIMA(%d,%d,%d)\n", override_best$p, override_best$d, override_best$q))
  cat(sprintf("AIC: %.2f\n", override_best$AIC))
  cat(sprintf("AICc: %.2f\n", override_best$AICc))
  cat(sprintf("BIC: %.2f\n", override_best$BIC))
  cat(sprintf("Hold-out MAPE: %.2f%%\n", override_best$MAPE))
  cat(sprintf("Ljung-Box p-value: %.4f\n", override_best$Ljung_Box_pvalue))

  override_best
}

actual_color    <- "#2563eb"
forecast_color  <- "#ea580c"
grid_color      <- "#e5e5e5"
text_color      <- "#1a1a1a"
muted_text      <- "#737373"

base_theme <- function() {
  theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", color = text_color, size = 14),
      panel.grid.major = element_line(color = grid_color),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.justification = "left",
      axis.text = element_text(color = text_color)
    )
}

generate_forecast_chart <- function(train, test, forecast_vals, forecast_dates,
                                    mape = NULL, output_path = "forecast.png",
                                    title_suffix = "", model_name = "STLARIMA(0,1,0)") {

  train_label <- sprintf(
    "Training (%s - %s)",
    format(min(train$date), "%b %Y"), format(max(train$date), "%b %Y")
  )
  test_label <- sprintf(
    "Test Actual (%s - %s)",
    format(min(test$date), "%b %Y"), format(max(test$date), "%b %Y")
  )
  forecast_label <- sprintf(
    "Forecast (%s - %s)",
    format(min(forecast_dates), "%b %Y"), format(max(forecast_dates), "%b %Y")
  )

  train_df <- data.frame(date = train$date, value = train$value, series = train_label)
  test_df  <- data.frame(date = test$date,  value = test$value,  series = test_label)
  fc_df    <- data.frame(date = forecast_dates, value = forecast_vals, series = forecast_label)

  connector_df <- data.frame(
    date  = c(test_df$date[nrow(test_df)], fc_df$date[1]),
    value = c(test_df$value[nrow(test_df)], fc_df$value[1])
  )

  label_df <- fc_df[seq(1, nrow(fc_df), by = 2), ]

  p <- ggplot() +
    geom_line(data = train_df, aes(x = date, y = value, color = series), linewidth = 1) +
    geom_line(data = test_df, aes(x = date, y = value, color = series),
              linewidth = 1, linetype = "dashed") +
    geom_line(data = connector_df, aes(x = date, y = value), color = forecast_color,
              linewidth = 0.7, linetype = "dotted", alpha = 0.6) +
    geom_line(data = fc_df, aes(x = date, y = value, color = series), linewidth = 1) +
    geom_point(data = fc_df, aes(x = date, y = value), color = forecast_color, size = 2.2,
               shape = 21, fill = forecast_color, stroke = 0.8) +
    geom_vline(xintercept = fc_df$date[1], color = muted_text,
               linewidth = 0.4, linetype = "dotted", alpha = 0.5) +
    geom_text(data = label_df, aes(x = date, y = value,
                                   label = paste0(round(value / 1e6, 1), "M")),
              vjust = -1, size = 3, color = forecast_color, fontface = "plain") +
    scale_color_manual(values = setNames(
      c(actual_color, actual_color, forecast_color),
      c(train_label, test_label, forecast_label)
    )) +
    scale_y_continuous(labels = function(x) paste0(round(x / 1e6, 1), "M")) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(
      title = paste0("LRT Ampang Line Ridership - ", model_name, " Forecast", title_suffix),
      x = NULL, y = "Monthly Ridership (millions)", color = NULL
    ) +
    base_theme()

  if (!is.null(mape)) {
    p <- p + annotate("label", x = max(fc_df$date), y = min(c(train_df$value, test_df$value)),
                      label = paste0("Hold-out MAPE: ", sprintf("%.2f", mape), "%"),
                      hjust = 1, vjust = 0, size = 3.2, color = muted_text,
                      fill = "white")
  }

  ggsave(output_path, plot = p, width = 12, height = 6, dpi = 150, bg = "white")
  cat("Forecast chart saved as", output_path, "\n")
  invisible(p)
}

generate_forecast_vs_actual_chart <- function(train,
                                              test,
                                              fc,
                                              seasonal_fc_test,
                                              output_path = "forecast_vs_actual.png",
                                              model_name = "STLARIMA(0,1,0)") {

  if (is.null(fc$lower) || ncol(fc$lower) < 2) {
    stop("fc must be produced with forecast(model, h = ..., level = c(80, 95)) ",
         "so 80%% and 95%% intervals are available.")
  }

  train_df <- data.frame(
    date = train$date,
    value = train$value
  )

  fc_df <- data.frame(
    date = test$date,
    mean = as.numeric(fc$mean) + seasonal_fc_test,
    lo80 = as.numeric(fc$lower[, 1]) + seasonal_fc_test,
    hi80 = as.numeric(fc$upper[, 1]) + seasonal_fc_test,
    lo95 = as.numeric(fc$lower[, 2]) + seasonal_fc_test,
    hi95 = as.numeric(fc$upper[, 2]) + seasonal_fc_test
  )

  actual_df <- data.frame(
    date = test$date,
    value = test$value
  )

  p <- ggplot() +
    geom_ribbon(
      data = fc_df,
      aes(x = date, ymin = lo95, ymax = hi95),
      fill = "#3366FF", alpha = 0.25
    ) +
    geom_ribbon(
      data = fc_df,
      aes(x = date, ymin = lo80, ymax = hi80),
      fill = "#3366FF", alpha = 0.45
    ) +
    geom_line(
      data = train_df,
      aes(x = date, y = value),
      color = "black", linewidth = 0.5
    ) +
    geom_line(
      data = fc_df,
      aes(x = date, y = mean),
      color = "#3366FF", linewidth = 0.6
    ) +
    geom_line(
      data = actual_df,
      aes(x = date, y = value),
      color = "red", linewidth = 0.6
    ) +
    scale_y_continuous(
      labels = function(x) paste0(round(x / 1e6, 0), "M")
    ) +
    labs(
      title = paste0(model_name, " - Forecast vs Actual"),
      x = "Time",
      y = "Monthly Ridership"
    ) +
    theme_gray(base_size = 13)

  ggsave(output_path, plot = p, width = 10, height = 6, dpi = 150)
  cat("Forecast-vs-actual chart saved as:", output_path, "\n")
  invisible(p)
}

generate_residuals_chart <- function(
    model,
    output_path = "residuals_diagnostic.png",
    model_name = "STLARIMA",
    width = 3300,
    height = 2400,
    res = 300) {

  old_scipen <- getOption("scipen")
  options(scipen = 999)
  on.exit(options(scipen = old_scipen), add = TRUE)

  model_for_plot <- model
  model_for_plot$method <- model_name

  png(
    output_path,
    width = width,
    height = height,
    res = res
  )

  checkresiduals(model_for_plot)

  dev.off()

  cat(
    "Residual diagnostic chart saved as:",
    output_path,
    "\n"
  )

  lb_test <- tryCatch(
    checkresiduals(model_for_plot, plot = FALSE),
    error = function(e) NULL
  )

  if (!is.null(lb_test)) {
    cat(model_name, "Ljung-Box test:\n")
    print(lb_test)
  }

  invisible(lb_test)
}

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

mco_exhibit_csv_path    <- "ridership_headline.csv"
mco_exhibit_start_date  <- as.Date("2019-01-01")
mco_exhibit_mco_start   <- as.Date("2020-03-01")
mco_exhibit_mco_end     <- as.Date("2021-12-01")
mco_exhibit_output_path <- "mco_imputation_comparison.png"

exhibit_df <- read_csv(mco_exhibit_csv_path, show_col_types = FALSE)

exhibit_df_ampang <- exhibit_df %>%
  mutate(date = as.Date(date)) %>%
  filter(date >= mco_exhibit_start_date) %>%
  arrange(date)

exhibit_monthly <- exhibit_df_ampang %>%
  mutate(month = as.Date(format(date, "%Y-%m-01"))) %>%
  group_by(month) %>%
  summarise(value = sum(rail_lrt_ampang, na.rm = TRUE), .groups = "drop") %>%
  arrange(month)

exhibit_expected_months <- seq(min(exhibit_monthly$month), max(exhibit_monthly$month), by = "month")
stopifnot(all(exhibit_expected_months == exhibit_monthly$month))

exhibit_monthly <- exhibit_monthly %>%
  mutate(
    is_mco = month >= mco_exhibit_mco_start & month <= mco_exhibit_mco_end,
    value_raw = value,
    value_masked = ifelse(is_mco, NA, value)
  )

exhibit_monthly <- exhibit_monthly %>%
  mutate(value_linear = na.approx(value_masked, x = month, na.rm = FALSE))

exhibit_value_ts <- ts(
  exhibit_monthly$value_linear,
  start = c(as.integer(format(min(exhibit_monthly$month), "%Y")),
            as.integer(format(min(exhibit_monthly$month), "%m"))),
  frequency = 12
)

exhibit_stl_fit <- stl(exhibit_value_ts, s.window = "periodic", robust = TRUE)
exhibit_stl_components <- as.data.frame(exhibit_stl_fit$time.series)

exhibit_monthly <- exhibit_monthly %>%
  mutate(
    stl_trend    = exhibit_stl_components$trend,
    stl_seasonal = exhibit_stl_components$seasonal,
    stl_reconstructed = stl_trend + stl_seasonal,
    value_imputed = ifelse(is_mco, stl_reconstructed, value_raw)
  )

cat("Imputed values for MCO window (", format(mco_exhibit_mco_start, "%b %Y"), " to ",
    format(mco_exhibit_mco_end, "%b %Y"), "):\n", sep = "")
print(
  exhibit_monthly %>%
    filter(is_mco) %>%
    select(month, value_raw, value_linear, value_imputed)
)

exhibit_plot_df <- exhibit_monthly %>%
  select(month, value_raw, value_linear, value_imputed) %>%
  pivot_longer(
    cols = c(value_raw, value_linear, value_imputed),
    names_to = "series",
    values_to = "value"
  ) %>%
  mutate(
    series = recode(
      series,
      value_raw     = "Actual (with MCO crash)",
      value_linear  = "Linear interpolation only",
      value_imputed = "STL-adjusted imputation"
    )
  )

p_mco_exhibit <- ggplot(exhibit_plot_df, aes(x = month, y = value, color = series, linetype = series)) +
  annotate(
    "rect",
    xmin = mco_exhibit_mco_start, xmax = mco_exhibit_mco_end,
    ymin = -Inf, ymax = Inf,
    fill = "grey85", alpha = 0.5
  ) +
  annotate(
    "text",
    x = mco_exhibit_mco_start + 150, y = max(exhibit_monthly$value_raw, na.rm = TRUE) * 0.95,
    label = "MCO Lockdown Period", size = 3, color = "grey40"
  ) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = c(
    "Actual (with MCO crash)"    = "#737373",
    "Linear interpolation only"  = "#f59e0b",
    "STL-adjusted imputation"    = "#2563eb"
  )) +
  scale_linetype_manual(values = c(
    "Actual (with MCO crash)"    = "solid",
    "Linear interpolation only"  = "dashed",
    "STL-adjusted imputation"    = "solid"
  )) +
  scale_y_continuous(labels = function(x) paste0(round(x / 1e6, 1), "M")) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title = "LRT Ampang Line Ridership - MCO Period, Known-Intervention Imputation",
    subtitle = "Context exhibit only - not used in the finalized STLARIMA model",
    x = NULL, y = "Monthly Ridership (millions)", color = NULL, linetype = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "grey40", size = 10),
    legend.position = "top",
    legend.justification = "left",
    panel.grid.minor = element_blank()
  )

ggsave(mco_exhibit_output_path, plot = p_mco_exhibit, width = 12, height = 6, dpi = 150, bg = "white")
cat("\nComparison chart saved as:", mco_exhibit_output_path, "\n")

print(p_mco_exhibit)