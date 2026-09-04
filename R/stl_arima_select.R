library(dplyr)

# ---------------------------------------------------------------------------
# select_stl_arima_model()
#
# Applies the same two-stage selection methodology used for the
# non-decomposed ARIMA model in this project:
#   1. Mechanical winner: lowest hold-out MAPE among models that pass the
#      Ljung-Box test (p > 0.05).
#   2. Cross-criterion override: if that mechanical winner beats a
#      simpler, lower-AIC candidate by less than `mape_gap_threshold`
#      percentage points, that gap is within noise for this sample size,
#      so the lower-AIC/more parsimonious candidate is selected instead.
#
# Prints the comparison so the reasoning is visible in the console output,
# and returns the single-row data frame for the final selected model.
# ---------------------------------------------------------------------------
select_stl_arima_model <- function(results, mape_gap_threshold = 1.0) {
  
  best_result <- results[1, ]
  
  valid_models <- results %>%
    filter(Ljung_Box_pvalue > 0.05) %>%
    arrange(MAPE)
  
  cat("\n============================================================\n")
  cat("MODELS PASSING LJUNG-BOX TEST (p > 0.05)\n")
  cat("============================================================\n")
  print(head(valid_models, 10), row.names = FALSE)
  
  if (nrow(valid_models) > 0) {
    raw_best <- valid_models[1, ]
  } else {
    cat("\nNo STLARIMA models passed the Ljung-Box test at p > 0.05.\n")
    cat("Falling back to the lowest-MAPE model overall.\n")
    raw_best <- best_result
    valid_models <- best_result
  }
  
  cat("\n============================================================\n")
  cat("MECHANICAL WINNER (lowest MAPE among Ljung-Box passes)\n")
  cat("============================================================\n")
  cat(sprintf(
    "STLARIMA(%d,%d,%d)  MAPE=%.2f%%  AIC=%.2f  BIC=%.2f  LB p=%.4f\n",
    raw_best$p, raw_best$d, raw_best$q,
    raw_best$MAPE, raw_best$AIC, raw_best$BIC, raw_best$Ljung_Box_pvalue
  ))
  
  near_tied <- valid_models %>%
    filter(MAPE <= raw_best$MAPE + mape_gap_threshold) %>%
    arrange(AIC)
  
  if (nrow(near_tied) > 0) {
    override_best <- near_tied[1, ]
  } else {
    # No Ljung-Box-passing candidate within the MAPE margin (should not
    # normally happen, since raw_best itself always qualifies) - fall
    # back to the mechanical winner.
    override_best <- raw_best
  }
  
  cat("\n============================================================\n")
  cat("CROSS-CRITERION OVERRIDE CHECK (AIC/AICc/BIC/parsimony)\n")
  cat("============================================================\n")
  cat(sprintf(
    "Mechanical MAPE winner:        STLARIMA(%d,%d,%d)  MAPE=%.2f%%  AIC=%.2f  BIC=%.2f\n",
    raw_best$p, raw_best$d, raw_best$q, raw_best$MAPE, raw_best$AIC, raw_best$BIC
  ))
  cat(sprintf(
    "Lowest-AIC within %.1fpp MAPE:  STLARIMA(%d,%d,%d)  MAPE=%.2f%%  AIC=%.2f  BIC=%.2f\n",
    mape_gap_threshold, override_best$p, override_best$d, override_best$q,
    override_best$MAPE, override_best$AIC, override_best$BIC
  ))
  
  if (override_best$p == 0 && override_best$q == 0) {
    cat("\nNote: the overridden candidate has no AR or MA terms. Check the\n")
    cat("coefficient significance after refitting - if p and q were both 0\n")
    cat("anyway, this confirms the STL-adjusted series needs no further\n")
    cat("ARIMA structure beyond differencing.\n")
  }
  
  cat("\n============================================================\n")
  cat("SELECTED FINAL STLARIMA MODEL\n")
  cat("============================================================\n")
  cat(sprintf("STLARIMA(%d,%d,%d)\n", override_best$p, override_best$d, override_best$q))
  cat(sprintf("AIC: %.2f\n", override_best$AIC))
  cat(sprintf("AICc: %.2f\n", override_best$AICc))
  cat(sprintf("BIC: %.2f\n", override_best$BIC))
  cat(sprintf("Hold-out MAPE: %.2f%%\n", override_best$MAPE))
  cat(sprintf("Ljung-Box p-value: %.4f\n", override_best$Ljung_Box_pvalue))
  
  override_best
}