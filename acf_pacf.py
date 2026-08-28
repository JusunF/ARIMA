import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from statsmodels.tsa.stattools import acf, pacf, adfuller, kpss
from statsmodels.graphics.tsaplots import plot_acf, plot_pacf
from statsmodels.tsa.seasonal import STL

# ============================================================
# 1. LOAD DATA (monthly aggregated)
# ============================================================

FILE_PATH = "ridership_ampang_monthly.csv"
DATE_COLUMN = "date"
VALUE_COLUMN = "rail_lrt_ampang"

df = pd.read_csv(FILE_PATH)
df[DATE_COLUMN] = pd.to_datetime(df[DATE_COLUMN])
df[VALUE_COLUMN] = pd.to_numeric(df[VALUE_COLUMN], errors="coerce")
df = df.dropna(subset=[DATE_COLUMN, VALUE_COLUMN]).sort_values(DATE_COLUMN)

monthly = df.set_index(DATE_COLUMN)[VALUE_COLUMN].resample("MS").sum()
monthly = monthly.loc["2019-01-01":"2026-06-01"].asfreq("MS")

# ============================================================
# 2. MCO IMPUTATION (matching R script logic)
#    - Linear interpolation across MCO window (2020-03 to 2021-12)
#    - STL decomposition on interpolated series
#    - Replace MCO window with STL trend + seasonal (dropping remainder)
# ============================================================

mco_start = pd.Timestamp("2020-03-01")
mco_end = pd.Timestamp("2021-12-01")

# Convert to float to allow STL imputation (produces float values)
monthly_raw = monthly.astype(float).copy()

# Step 1: Linear interpolation across MCO window
is_mco = (monthly.index >= mco_start) & (monthly.index <= mco_end)
monthly_masked = monthly_raw.copy()
monthly_masked[is_mco] = np.nan
monthly_linear = monthly_masked.interpolate(method="time")

# Step 2: STL on interpolated series, reconstruct MCO window from trend + seasonal
# STL requires a pandas Series with a DatetimeIndex with known frequency
monthly_linear.index.freq = "MS"  # ensure frequency is set
stl_fit = STL(monthly_linear, period=12, robust=True).fit()
stl_components = pd.DataFrame({
    'trend': stl_fit.trend,
    'seasonal': stl_fit.seasonal,
    'resid': stl_fit.resid
}, index=monthly_linear.index)

stl_reconstructed = stl_components['trend'] + stl_components['seasonal']

# Replace MCO window with STL reconstruction
monthly_imputed = monthly_raw.copy()
monthly_imputed[is_mco] = stl_reconstructed[is_mco]

# Print imputed values for verification
print("=" * 60)
print("MCO IMPUTATION RESULTS")
print("=" * 60)
print(f"MCO window: {mco_start.strftime('%b %Y')} to {mco_end.strftime('%b %Y')}")
imputed_df = pd.DataFrame({
    'raw': monthly_raw[is_mco],
    'linear_interpolated': monthly_linear[is_mco],
    'stl_imputed': monthly_imputed[is_mco]
})
print(imputed_df.to_string())
print()

# ============================================================
# 3. CONFIRM d (ORDER OF DIFFERENCING) VIA ADF + KPSS
#    Using the IMPUTED series
# ============================================================

TEST_SIZE = 11  # matches the finalized 43 train / 11 test split
train = monthly_imputed.iloc[:-TEST_SIZE]

print("=" * 60)
print("STEP 1: STATIONARITY CHECK (determines d) - USING IMPUTED DATA")
print("=" * 60)


def stationarity_report(series, label):
    adf_p = adfuller(series)[1]
    kpss_p = kpss(series, nlags="auto")[1]
    print(f"--- {label} ---")
    print(f"ADF p-value:  {adf_p:.4f}  ({'stationary' if adf_p < 0.05 else 'non-stationary'} per ADF)")
    print(f"KPSS p-value: {kpss_p:.4f}  ({'stationary' if kpss_p > 0.05 else 'non-stationary'} per KPSS)")
    print()


stationarity_report(train, "Level series (d=0)")

diff1 = train.diff().dropna()
stationarity_report(diff1, "First-differenced series (d=1)")

# Use the differenced series for ACF/PACF identification (d = 1)
series_to_use = diff1
d = 1


# ============================================================
# 4. COMPUTE ACF / PACF VALUES + SIGNIFICANCE BAND
# ============================================================

print("=" * 60)
print("STEP 2: ACF / PACF ON DIFFERENCED IMPUTED SERIES (identifies p, q)")
print("=" * 60)

n = len(series_to_use)
ci = 1.96 / np.sqrt(n)
n_lags = min(14, n // 3)  # keep lags well below n to avoid unstable estimates

acf_vals = acf(series_to_use, nlags=n_lags)
pacf_vals = pacf(series_to_use, nlags=n_lags)

print(f"n = {n}, 95% significance bound = +/- {ci:.4f}\n")
print("Lag | ACF     | Sig? | PACF    | Sig?")
acf_sig_lags = []
pacf_sig_lags = []
for lag in range(1, n_lags + 1):
    a = acf_vals[lag]
    p = pacf_vals[lag]
    a_sig = abs(a) > ci
    p_sig = abs(p) > ci
    if a_sig:
        acf_sig_lags.append(lag)
    if p_sig:
        pacf_sig_lags.append(lag)
    print(f"{lag:3d} | {a:7.4f} | {'YES' if a_sig else '':<4} | {p:7.4f} | {'YES' if p_sig else '':<4}")


# ============================================================
# 5. AUTO-SUGGEST (p, d, q) FROM THE CUTOFF PATTERN
# ============================================================
# Rule of thumb (Box-Jenkins):
#   PACF cuts off at lag p, ACF tails off      -> AR(p):  ARIMA(p, d, 0)
#   ACF cuts off at lag q, PACF tails off      -> MA(q):  ARIMA(0, d, q)
#   Both cut off quickly (lag 1 only)          -> compare which is "cleaner"
#   Both tail off gradually                    -> mixed ARMA, harder to read directly

def first_gap_after_start(sig_lags, max_lag):
    """Returns the lag of the last significant spike before the first
    non-significant lag (i.e. where the pattern 'cuts off')."""
    if not sig_lags or 1 not in sig_lags:
        return 0
    cutoff = 1
    for lag in range(1, max_lag + 1):
        if lag in sig_lags:
            cutoff = lag
        else:
            break
    return cutoff


p_suggested = first_gap_after_start(pacf_sig_lags, n_lags)
q_suggested = first_gap_after_start(acf_sig_lags, n_lags)

print("\n" + "-" * 60)
print("INTERPRETATION")
print("-" * 60)
print(f"PACF cuts off after lag: {p_suggested}  -> suggests p = {p_suggested}")
print(f"ACF cuts off after lag:  {q_suggested}  -> suggests q = {q_suggested}")

# Decide which pattern is cleaner: AR signature vs MA signature
pacf_tail = len(pacf_sig_lags) > 1 and pacf_sig_lags != list(range(1, len(pacf_sig_lags) + 1))
acf_tail = len(acf_sig_lags) > 1 and acf_sig_lags != list(range(1, len(acf_sig_lags) + 1))

print()
if p_suggested >= 1 and q_suggested >= 1 and p_suggested <= 2 and q_suggested >= 3:
    # PACF cuts off quickly, ACF keeps oscillating/decaying -> AR signature
    recommendation = (p_suggested, d, 0)
    print(f"PACF cuts off sharply while ACF decays gradually -> AR({p_suggested}) signature")
elif q_suggested >= 1 and p_suggested >= 3:
    recommendation = (0, d, q_suggested)
    print(f"ACF cuts off sharply while PACF decays gradually -> MA({q_suggested}) signature")
elif p_suggested == 1 and q_suggested == 1:
    recommendation = (1, d, 0)
    print("Both ACF and PACF are significant only at lag 1.")
    print("Defaulting to the simpler AR(1) reading: ARIMA(1,{},0) -- ".format(d)
          + "cross-check against ACF decay shape in the plot.")
else:
    recommendation = (p_suggested, d, q_suggested)
    print("No clean single cutoff pattern -- consider this a candidate ARMA(p,q) "
          "to validate further (AIC / hold-out MAPE / Ljung-Box).")

print(f"\n>>> SUGGESTED ORDER FROM ACF/PACF (on imputed data): ARIMA{recommendation} <<<")


# ============================================================
# 6. PLOT ACF / PACF
# ============================================================

fig, axes = plt.subplots(2, 1, figsize=(10, 8))

plot_acf(series_to_use, lags=n_lags, ax=axes[0], color="#2563eb")
axes[0].set_title(f"ACF of Differenced Training Series (d={d}, MCO-imputed)", fontsize=13, fontweight="bold")
axes[0].set_xlabel("Lag (months)")
axes[0].set_ylabel("Autocorrelation")

plot_pacf(series_to_use, lags=n_lags, ax=axes[1], method="ywm", color="#ea580c")
axes[1].set_title(f"PACF of Differenced Training Series (d={d}, MCO-imputed)", fontsize=13, fontweight="bold")
axes[1].set_xlabel("Lag (months)")
axes[1].set_ylabel("Partial Autocorrelation")

plt.tight_layout()
plt.savefig("acf_pacf_identification.png", dpi=150, bbox_inches="tight")
print("\nSaved chart: acf_pacf_identification.png")