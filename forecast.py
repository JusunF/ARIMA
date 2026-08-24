"""Forecast-chart generation for LRT Ampang ARIMA analysis — with drift comparison."""

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd


def generate_forecast_chart(train, test, forecast, mape, output_path="arima_forecast.png", show=True, title_suffix=""):
    """Create and save the training, test, and forecast chart."""
    forecast_dates = pd.date_range(start=test.index[-1] + pd.DateOffset(months=1), periods=len(forecast), freq="MS")
    actual_color, forecast_color = "#2563eb", "#ea580c"
    grid_color, text_color, muted_text = "#e5e5e5", "#1a1a1a", "#737373"
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(train.index, train.values, color=actual_color, linewidth=2, label="Training (Jan 2022–Dec 2025)", zorder=3)
    ax.plot(test.index, test.values, color=actual_color, linewidth=2, linestyle="--", label="Test Actual (Jan–Jun 2026)", zorder=3)
    ax.plot(forecast_dates, forecast.values, color=forecast_color, linewidth=2, label="Forecast (Jul–Dec 2026)", zorder=3)
    ax.scatter(forecast_dates, forecast.values, color=forecast_color, s=40, zorder=4, edgecolors="white", linewidths=1.5)
    ax.plot([test.index[-1], forecast_dates[0]], [test.values[-1], forecast.values[0]], color=forecast_color, linewidth=1.5, linestyle=":", alpha=0.6, zorder=2)
    ax.axvline(x=forecast_dates[0], color=muted_text, linewidth=1, linestyle=":", alpha=0.5, zorder=1)
    ax.set_title(f"MRT Kajang Line Ridership — ARIMA(1,1,0) Forecast{title_suffix}", fontsize=14, fontweight="600", color=text_color, pad=16)
    ax.set_ylabel("Monthly Ridership (millions)", fontsize=11, color=text_color, labelpad=10)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda value, _: f"{value / 1e6:.1f}M"))
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.xaxis.set_minor_locator(mdates.MonthLocator((1, 7)))
    ax.tick_params(axis="x", which="major", labelsize=10, color=text_color)
    ax.tick_params(axis="y", labelsize=10, color=text_color)
    ax.grid(True, color=grid_color, linewidth=1, alpha=0.5)
    ax.set_axisbelow(True)
    legend = ax.legend(loc="upper left", frameon=True, facecolor="white", edgecolor=grid_color, fontsize=10, labelcolor=text_color)
    legend.get_frame().set_linewidth(0.5)
    for index, (date, value) in enumerate(zip(forecast_dates, forecast.values)):
        if index % 2 == 0:
            ax.annotate(f"{value / 1e6:.1f}M", xy=(date, value), xytext=(0, 12), textcoords="offset points", ha="center", va="bottom", fontsize=9, color=forecast_color, fontweight="500", bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor=forecast_color, linewidth=0.5))
    ax.annotate(f"Hold-out MAPE: {mape:.2f}%", xy=(0.98, 0.02), xycoords="axes fraction", ha="right", va="bottom", fontsize=10, color=muted_text, bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor=grid_color, linewidth=0.5))
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"\nForecast chart saved as '{output_path}'")
    if show:
        plt.show()
    else:
        plt.close(fig)