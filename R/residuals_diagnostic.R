library(forecast)

generate_residuals_chart <- function(
    model,
    output_path = "arima_residuals.png",
    model_name = "ARIMA",
    width = 3300,
    height = 2400,
    res = 300) {

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