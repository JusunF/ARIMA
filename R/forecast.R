# forecast.R
# Forecast chart generation for LRT Ampang Line ARIMA analysis.
#
# generate_forecast_chart() -> full history + the future forecast (model
# fit on ALL data), annotated with the hold-out MAPE for credibility.

library(ggplot2)
library(scales)

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

# ---------------------------------------------------------------------------
# Final forecast: full history (train + test) plus the future forecast,
# with a dotted connector and the hold-out MAPE annotated. All series
# labels/date ranges are derived from the data, never hardcoded.
# ---------------------------------------------------------------------------
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