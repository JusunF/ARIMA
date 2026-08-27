library(ggplot2)

# Generic ACF/PACF chart pair. Used for two DIFFERENT purposes in this
# project - keep that distinction straight when calling it:
#   1) Fig. 1 (order-selection evidence): call on diff(train_ts), BEFORE
#      fitting, to justify the AR/MA order.
#   2) Residual diagnostics (appendix): call on residuals(model), AFTER
#      fitting, to confirm the residuals are white noise.
generate_acf_pacf_chart <- function(
    x,
    model_name = "ARIMA(1,1,0)",
    acf_path = "acf.png",
    pacf_path = "pacf.png") {
  
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
    value = acf_res$acf[, 1, 1]
  )
  
  pacf_df <- data.frame(
    lag = pacf_res$lag[, 1, 1],
    value = pacf_res$acf[, 1, 1]
  )
  
  # Remove lag 0
  acf_df <- acf_df[acf_df$lag != 0, ]
  pacf_df <- pacf_df[pacf_df$lag != 0, ]
  
  # Function for drawing a single chart
  draw_chart <- function(
    df,
    title,
    ylab) {
    
    ggplot(df, aes(x = lag, y = value)) +
      
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
      
      scale_y_continuous(
        limits = c(-1.1, 1.1)
      ) +
      
      labs(
        title = paste0(title, "\n", model_name),
        x = "Lag (months)",
        y = ylab
      ) +
      
      theme_minimal(base_size = 12) +
      
      theme(
        plot.title = element_text(
          face = "bold",
          color = text_color,
          size = 13
        ),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(
          color = grid_color
        ),
        axis.text = element_text(
          color = text_color
        )
      )
  }
  
  # ACF
  acf_plot <- draw_chart(
    acf_df,
    "Autocorrelation Function (ACF)",
    "Autocorrelation"
  )
  
  # PACF
  pacf_plot <- draw_chart(
    pacf_df,
    "Partial Autocorrelation Function (PACF)",
    "Partial Autocorrelation"
  )
  
  # Save separately
  ggsave(
    acf_path,
    plot = acf_plot,
    width = 12,
    height = 6,
    dpi = 150,
    bg = "white"
  )
  
  ggsave(
    pacf_path,
    plot = pacf_plot,
    width = 12,
    height = 6,
    dpi = 150,
    bg = "white"
  )
  
  cat("ACF chart saved as:", acf_path, "\n")
  cat("PACF chart saved as:", pacf_path, "\n")
}