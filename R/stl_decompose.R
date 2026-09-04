library(forecast)

# ---------------------------------------------------------------------------
# perform_stl_decomposition()
#
# Runs STL decomposition on a monthly ts() object and returns the fitted
# stl object alongside the two pieces needed downstream:
#   - seasadj: the seasonally-adjusted series (trend + remainder), which is
#     what the ARIMA(p,d,q) component is fit to in STLARIMA.
#   - seasonal: the raw seasonal component, used to extrapolate the
#     seasonal pattern forward for forecasting.
# ---------------------------------------------------------------------------
perform_stl_decomposition <- function(ts_data, s_window = "periodic", robust = TRUE) {

  stl_fit <- stl(
    ts_data,
    s.window = s_window,
    robust = robust
  )

  list(
    stl_fit = stl_fit,
    seasadj = seasadj(stl_fit),
    seasonal = stl_fit$time.series[, "seasonal"]
  )
}

# ---------------------------------------------------------------------------
# extrapolate_seasonal()
#
# Seasonal-naive extrapolation of an STL seasonal component: repeats the
# last observed 12-month seasonal cycle for `h` future periods. This is
# added back to the ARIMA forecast (fit on the seasonally-adjusted series)
# to reconstruct the final STLARIMA forecast on the original scale.
# ---------------------------------------------------------------------------
extrapolate_seasonal <- function(seasonal_component, h) {

  last_cycle <- tail(as.numeric(seasonal_component), 12)

  rep(last_cycle, length.out = h)
}