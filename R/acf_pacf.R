library(ggplot2)

# ---------------------------------------------------------------------------
# Generate combined ACF/PACF comparison chart.
#
# Purpose:
#   Used BEFORE fitting the ARIMA model, on the differenced training series,
#   to provide evidence for selecting the AR and MA orders.
#
# Input:
#   x = time series to analyse
#   model_name = label shown in the chart title
#   output_path = output PNG filename
# ---------------------------------------------------------------------------

generate_acf_pacf_chart <- function(
    x,
    model_name = "Differenced Series (d = 1)",
    output_path = "acf_pacf_comparison.png") {
  
  x <- as.numeric(x)
  n <- length(x)
  n_lags <- min(40, n %/% 3)
  
  # Calculate ACF and PACF
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
  
  # 95% confidence interval
  ci <- 1.96 / sqrt(n)
  
  color <- "#2563eb"
  grid_color <- "#e5e5e5"
  text_color <- "#1a1a1a"
  muted_text <- "#737373"
  
  # Convert results to data frames
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
  
  # Remove lag 0
  acf_df <- acf_df[acf_df$lag != 0, ]
  pacf_df <- pacf_df[pacf_df$lag != 0, ]
  
  # Combine ACF and PACF
  comparison_df <- rbind(acf_df, pacf_df)
  
  # Combined chart
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
  
  # Save combined chart
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