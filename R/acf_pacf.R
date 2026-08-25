# acf_pacf.R
# ACF/PACF residual-diagnostic chart generation for ARIMA analysis (R version)
# Mirrors acf_pacf.py's generate_acf_pacf_chart()

library(ggplot2)
library(patchwork)

generate_acf_pacf_chart <- function(residuals, model_name = "ARIMA(1,1,0)",
                                     output_path = "arima_acf_pacf.png") {
  residuals <- as.numeric(residuals)
  n <- length(residuals)
  n_lags <- min(40, n %/% 3)

  acf_res  <- acf(residuals, lag.max = n_lags, plot = FALSE)
  pacf_res <- pacf(residuals, lag.max = n_lags, plot = FALSE)

  ci <- 1.96 / sqrt(n)
  color      <- "#2563eb"
  grid_color <- "#e5e5e5"
  text_color <- "#1a1a1a"
  muted_text <- "#737373"

  acf_df  <- data.frame(lag = acf_res$lag[, 1, 1],  value = acf_res$acf[, 1, 1])
  pacf_df <- data.frame(lag = pacf_res$lag[, 1, 1], value = pacf_res$acf[, 1, 1])

  draw_chart <- function(df, title, ylab, drop_lag0 = FALSE) {
    if (drop_lag0) df <- df[df$lag != 0, ]
    sig_df <- df[abs(df$value) > ci & df$lag != 0, ]

    ggplot(df, aes(x = lag, y = value)) +
      geom_hline(yintercept = 0, color = text_color, linewidth = 0.4) +
      geom_ribbon(aes(x = lag, ymin = -ci, ymax = ci), inherit.aes = FALSE,
                  fill = grid_color, alpha = 0.5) +
      geom_col(width = 0.4, fill = color, color = color, alpha = ifelse(df$lag == 0, 0.3, 1)) +
      { if (nrow(sig_df) > 0)
          geom_text(data = sig_df, aes(x = lag, y = value,
                    label = paste0("lag ", lag)),
                    vjust = ifelse(sig_df$value > 0, -0.6, 1.4),
                    size = 2.6, color = muted_text)
      } +
      scale_y_continuous(limits = c(-1.1, 1.1)) +
      labs(title = title, x = NULL, y = ylab) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", color = text_color, size = 13),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = text_color)
      )
  }

  p1 <- draw_chart(acf_df,  "Autocorrelation Function (ACF) of Residuals", "Autocorrelation")
  p2 <- draw_chart(pacf_df, "Partial Autocorrelation Function (PACF) of Residuals", "Partial Autocorrelation") +
    labs(x = "Lag (months)")

  caption_text <- paste0(model_name, " residuals | n=", n,
                          " | 95% CI: \u00b1", sprintf("%.3f", ci))
  p1 <- p1 + annotate("label", x = -Inf, y = -Inf, label = caption_text,
                       hjust = -0.02, vjust = -0.3, size = 3, color = muted_text,
                       fill = "white", label.size = 0.3)

  combined <- p1 / p2

  ggsave(output_path, plot = combined, width = 12, height = 8, dpi = 150, bg = "white")
  cat("\nACF/PACF chart saved as", output_path, "\n")
  invisible(combined)
}