# forecast_vs_actual.R
#
# Recreates the "SARIMA(...) - Forecast vs Actual" style chart:
#   - default ggplot2 grey theme (no custom theme/palette)
#   - training history in black
#   - forecast mean in blue, with 80%/95% shaded confidence ribbons
#   - actual hold-out values overlaid in red on top of the forecast
#
# Usage (inside STLARIMA.R, AFTER model_stl_arima has been fit):
#
#   fc_stl_arima <- forecast(model_stl_arima, h = length(test_vals), level = c(80, 95))
#
#   generate_forecast_vs_actual_chart(
#     train, test, fc_stl_arima, seasonal_fc_test,
#     output_path = "forecast_vs_actual.png",
#     model_name  = model_label   # e.g. "STLARIMA(2,1,2)"
#   )
#
# Note: `fc_stl_arima$mean/$lower/$upper` are on the seasonally-ADJUSTED
# scale. seasonal_fc_test (the deterministic seasonal-naive extrapolation)
# must be added back to the mean AND to the lower/upper bounds, since it
# doesn't contribute any extra uncertainty of its own - it's a constant
# shift per period, not a random component.

library(ggplot2)

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
      y = "train"
    ) +
    theme_gray(base_size = 13)
  
  ggsave(output_path, plot = p, width = 10, height = 6, dpi = 150)
  cat("Forecast-vs-actual chart saved as:", output_path, "\n")
  invisible(p)
}