library(forecast)

generate_residuals_chart <- function(
    model,
    output_path = "residuals_diagnostic.png",
    model_name = "STLARIMA",
    width = 3300,
    height = 2400,
    res = 300) {
  
  # checkresiduals() draws a base-R plot, which switches to scientific
  # notation (1e+06) by default for large values. scipen strongly biases
  # base R toward plain/fixed notation instead. Scoped with on.exit() so
  # it reverts to the previous setting once this function finishes, and
  # doesn't affect any other chart in the project (e.g. the ggplot2
  # charts in forecast.R / acf_pacf.R already handle their own formatting).
  old_scipen <- getOption("scipen")
  options(scipen = 999)
  on.exit(options(scipen = old_scipen), add = TRUE)
  
  # checkresiduals() does NOT forward a `main =` argument through to its
  # plot title - `...` is only passed to the internal Box.test() call.
  # Both the plot title AND the "Residuals from ..." line printed with the
  # Ljung-Box test are built internally from object$method, which for an
  # Arima-derived fit is auto-generated as e.g. "ARIMA(2,1,2)" and knows
  # nothing about the STL wrapping. Overwrite $method on a copy of the
  # model (stats/residuals are untouched - only the label changes) and use
  # that copy everywhere checkresiduals() is called, so both the chart
  # title and the printed test output say STLARIMA instead of ARIMA.
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