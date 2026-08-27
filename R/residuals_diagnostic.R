library(forecast)

generate_residuals_chart <- function(
    model,
    output_path = "arima_residuals.png",
    model_name = "ARIMA",
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
  
  png(
    output_path,
    width = width,
    height = height,
    res = res
  )
  
  checkresiduals(model)
  
  dev.off()
  
  cat(
    "Residual diagnostic chart saved as:",
    output_path,
    "\n"
  )
  
  lb_test <- tryCatch(
    checkresiduals(model, plot = FALSE),
    error = function(e) NULL
  )
  
  if (!is.null(lb_test)) {
    cat(model_name, "Ljung-Box test:\n")
    print(lb_test)
  }
  
  invisible(lb_test)
}