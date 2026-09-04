library(forecast)

# ---------------------------------------------------------------------------
# run_stl_arima_grid_search()
#
# Fits ARIMA(p,d,q) for every combination in p_range x d_range x q_range on
# the seasonally-adjusted TRAINING series. For each model, forecasts the
# test horizon on the seas-adj scale, adds back seasonal_fc_test (the
# seasonal-naive extrapolation of the STL seasonal component), and scores
# the result against the ACTUAL test values.
#
# Returns a data frame of results sorted by hold-out MAPE (ascending),
# with AIC / AICc / BIC / Ljung-Box p-value (lag 12) for each model.
# ---------------------------------------------------------------------------
run_stl_arima_grid_search <- function(seasadj_train_ts,
                                      test_vals,
                                      seasonal_fc_test,
                                      p_range = 0:9,
                                      d_range = 0:9,
                                      q_range = 0:9,
                                      verbose = TRUE) {
  
  results <- data.frame(
    p = integer(),
    d = integer(),
    q = integer(),
    AIC = numeric(),
    AICc = numeric(),
    BIC = numeric(),
    MAPE = numeric(),
    Ljung_Box_pvalue = numeric()
  )
  
  for (p in p_range) {
    for (d in d_range) {
      for (q in q_range) {
        
        if (verbose) {
          cat("Fitting STLARIMA(", p, ",", d, ",", q, ")...\n", sep = "")
        }
        
        model <- tryCatch(
          Arima(seasadj_train_ts, order = c(p, d, q)),
          error = function(e) NULL
        )
        
        if (is.null(model)) next
        
        fc <- tryCatch(
          forecast(model, h = length(test_vals)),
          error = function(e) NULL
        )
        
        if (is.null(fc)) next
        
        predictions <- as.numeric(fc$mean) + seasonal_fc_test
        
        mape <- mean(abs((test_vals - predictions) / test_vals)) * 100
        
        lb <- tryCatch(
          Box.test(
            residuals(model),
            lag = 12,
            type = "Ljung-Box",
            fitdf = length(model$coef)
          ),
          error = function(e) NULL
        )
        
        lb_pvalue <- if (is.null(lb)) NA else lb$p.value
        
        k <- length(model$coef)
        n_obs <- length(seasadj_train_ts) - d
        aicc <- if ((n_obs - k - 1) > 0) {
          model$aic + (2 * k * (k + 1)) / (n_obs - k - 1)
        } else {
          NA
        }
        
        results <- rbind(
          results,
          data.frame(
            p = p, d = d, q = q,
            AIC = model$aic,
            AICc = aicc,
            BIC = model$bic,
            MAPE = mape,
            Ljung_Box_pvalue = lb_pvalue
          )
        )
      }
    }
  }
  
  results[order(results$MAPE), ]
}