# residuals_diagnostic.R
# Reusable residual-diagnostic chart function (mirrors Steven's checkresiduals() usage)
# Source this file, then call generate_residuals_chart() from ARIMA.R

library(forecast)

#' Save a checkresiduals() diagnostic panel (residuals-over-time, ACF, histogram)
#' for a fitted Arima() model, and optionally a separate PACF chart.
#'
#' @param model        A fitted Arima() model object (e.g. model_nodrift)
#' @param output_path  File path for the checkresiduals PNG
#' @param model_name   Label used in console output / plot title context
#' @param width,height,res  png() device settings (defaults match Steven's script)
#' @param include_pacf Logical; if TRUE also saves a companion PACF chart
#' @param pacf_path    File path for the PACF PNG (used only if include_pacf = TRUE)
#' @param lag_max      Max lag for the companion PACF chart
generate_residuals_chart <- function(model,
                                      output_path = "arima_residuals.png",
                                      model_name = "ARIMA",
                                      width = 3300, height = 2400, res = 300,
                                      include_pacf = TRUE,
                                      pacf_path = "arima_pacf.png",
                                      lag_max = 24) {

  # Main 3-panel diagnostic: residuals over time, ACF, histogram + normal curve
  png(output_path, width = width, height = height, res = res)
  checkresiduals(model)
  dev.off()
  cat("\nResidual diagnostic chart saved as", output_path, "\n")

  # checkresiduals() prints Ljung-Box test results to console automatically;
  # capture and return them so they can be logged / written to file if needed
  lb_test <- tryCatch(
    checkresiduals(model, plot = FALSE),
    error = function(e) NULL
  )
  if (!is.null(lb_test)) {
    cat(model_name, "Ljung-Box test:\n")
    print(lb_test)
  }

  # Optional companion PACF chart (checkresiduals already bundles an ACF panel;
  # standalone combined ACF+PACF lives in acf_pacf.R's generate_acf_pacf_chart())
  if (include_pacf) {
    png(pacf_path, width = 3000, height = 1500, res = 300)
    Pacf(residuals(model), lag.max = lag_max,
         main = paste(model_name, "Residuals PACF"))
    dev.off()
    cat("PACF chart saved as", pacf_path, "\n")
  }

  invisible(lb_test)
}