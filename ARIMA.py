import matplotlib
matplotlib.use("TkAgg")  # Explicit backend for interactive display

import pandas as pd
from pmdarima import auto_arima
from sklearn.metrics import mean_absolute_percentage_error
from statsmodels.stats.diagnostic import acorr_ljungbox
from statsmodels.tsa.arima.model import ARIMA

from forecast import generate_forecast_chart
from acf_pacf import generate_acf_pacf_chart

# Load and filter to MRT Kajang Line, Jan 2022 - Jun 2026
df = pd.read_csv("ridership_headline.csv", parse_dates=["date"])
df_ampang = df[df["date"] >= "2022-01-01"].copy()
ampang = df_ampang.set_index("date")["rail_lrt_ampang"].asfreq("MS")

# Fixed 80/20 split: 43 train (79.6%) / 11 test (20.4%) of 54 months
n_test = int(len(ampang) * 0.2)  # 10.8 -> adjust to 11 for whole months
n_test = 11
train, test = ampang[:-n_test], ampang[-n_test:]

# AIC-based search
auto_model = auto_arima(train, seasonal=False, trace=True)
print(auto_model.summary())

# No-drift model (trend="n")
model_nodrift = ARIMA(train, order=(1, 1, 0), trend="n").fit()
fc_nodrift = model_nodrift.forecast(steps=len(test))
mape_nodrift = mean_absolute_percentage_error(test, fc_nodrift) * 100
print("No-drift hold-out MAPE:", mape_nodrift)
print("No-drift Ljung-Box:", acorr_ljungbox(model_nodrift.resid, lags=[6, 12]))

# With-drift model (trend="t" for linear trend - needed when d=1)
model_drift = ARIMA(train, order=(1, 1, 0), trend="t").fit()
fc_drift = model_drift.forecast(steps=len(test))
mape_drift = mean_absolute_percentage_error(test, fc_drift) * 100
print("With-drift hold-out MAPE:", mape_drift)
print("With-drift Ljung-Box:", acorr_ljungbox(model_drift.resid, lags=[6, 12]))

# Final forecasts on full data, Jul-Dec 2026
final_nodrift = ARIMA(ampang, order=(1, 1, 0), trend="n").fit()
forecast_2026h2_nodrift = final_nodrift.forecast(steps=6)
print("No-drift final forecast:")
print(forecast_2026h2_nodrift)

final_drift = ARIMA(ampang, order=(1, 1, 0), trend="t").fit()
forecast_2026h2_drift = final_drift.forecast(steps=6)
print("With-drift final forecast:")
print(forecast_2026h2_drift)

# Generate charts
generate_forecast_chart(train, test, forecast_2026h2_nodrift, mape_nodrift,
                        output_path="arima_forecast_nodrift.png",
                        title_suffix=" (No Drift)")

generate_forecast_chart(train, test, forecast_2026h2_drift, mape_drift,
                        output_path="arima_forecast_drift.png",
                        title_suffix=" (With Drift)")

# ACF/PACF of no-drift final model residuals
generate_acf_pacf_chart(final_nodrift.resid)