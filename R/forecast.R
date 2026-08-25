# forecast.R
# Forecast-chart generation for LRT Ampang ARIMA analysis (R / ggplot2 version)
# Mirrors forecast.py's generate_forecast_chart()

library(ggplot2)
library(scales)

generate_forecast_chart <- function(train, test, forecast_vals, forecast_dates,
                                     mape, output_path = "arima_forecast.png",
                                     title_suffix = "") {
  # train, test: named numeric vectors/ts with Date-coercible names (or pass dates separately)
  # forecast_vals: numeric vector of forecast point estimates
  # forecast_dates: Date vector matching forecast_vals

  actual_color   <- "#2563eb"
  forecast_color <- "#ea580c"
  grid_color     <- "#e5e5e5"
  text_color     <- "#1a1a1a"
  muted_text     <- "#737373"

  train_df <- data.frame(date = train$date, value = train$value, series = "Training (Jan 2022-Jul 2025)")
  test_df  <- data.frame(date = test$date,  value = test$value,  series = "Test Actual (Aug 2025-Jun 2026)")
  fc_df    <- data.frame(date = forecast_dates, value = forecast_vals, series = "Forecast (Jul-Dec 2026)")

  # connector segment from last test point to first forecast point
  connector_df <- data.frame(
    date  = c(test_df$date[nrow(test_df)], fc_df$date[1]),
    value = c(test_df$value[nrow(test_df)], fc_df$value[1])
  )

  label_df <- fc_df[seq(1, nrow(fc_df), by = 2), ]

  p <- ggplot() +
    geom_line(data = train_df, aes(x = date, y = value, color = series), linewidth = 1) +
    geom_line(data = test_df, aes(x = date, y = value, color = series), linewidth = 1, linetype = "dashed") +
    geom_line(data = connector_df, aes(x = date, y = value), color = forecast_color,
              linewidth = 0.7, linetype = "dotted", alpha = 0.6) +
    geom_line(data = fc_df, aes(x = date, y = value, color = series), linewidth = 1) +
    geom_point(data = fc_df, aes(x = date, y = value), color = forecast_color, size = 2.2,
               shape = 21, fill = forecast_color, stroke = 0.8) +
    geom_vline(xintercept = as.numeric(fc_df$date[1]), color = muted_text,
               linewidth = 0.4, linetype = "dotted", alpha = 0.5) +
    geom_text(data = label_df, aes(x = date, y = value,
              label = paste0(round(value / 1e6, 1), "M")),
              vjust = -1, size = 3, color = forecast_color, fontface = "plain") +
    annotate("label", x = max(fc_df$date), y = min(c(train_df$value, test_df$value)),
             label = paste0("Hold-out MAPE: ", sprintf("%.2f", mape), "%"),
             hjust = 1, vjust = 0, size = 3.2, color = muted_text,
             fill = "white", label.size = 0.3) +
    scale_color_manual(values = setNames(
      c(actual_color, actual_color, forecast_color),
      c("Training (Jan 2022-Jul 2025)", "Test Actual (Aug 2025-Jun 2026)", "Forecast (Jul-Dec 2026)")
    )) +
    scale_y_continuous(labels = function(x) paste0(round(x / 1e6, 1), "M")) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(
      title = paste0("LRT Ampang Line Ridership - ARIMA(1,1,0) Forecast", title_suffix),
      x = NULL, y = "Monthly Ridership (millions)", color = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", color = text_color, size = 14),
      panel.grid.major = element_line(color = grid_color),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.justification = "left",
      axis.text = element_text(color = text_color)
    )

  ggsave(output_path, plot = p, width = 12, height = 6, dpi = 150, bg = "white")
  cat("\nForecast chart saved as", output_path, "\n")
  invisible(p)
}