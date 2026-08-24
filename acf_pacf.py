"""ACF/PACF residual-diagnostic chart generation for ARIMA analysis."""

import numpy as np
import matplotlib.pyplot as plt
from statsmodels.tsa.stattools import acf, pacf


def generate_acf_pacf_chart(residuals, model_name="ARIMA(1,1,0)", output_path="arima_acf_pacf.png", show=True):
    """Create and save ACF and PACF charts for model residuals."""
    residuals = np.asarray(residuals)
    n_lags = min(40, len(residuals) // 3)
    acf_values = acf(residuals, nlags=n_lags, fft=True)
    pacf_values = pacf(residuals, nlags=n_lags, method="ywm")
    confidence_interval = 1.96 / np.sqrt(len(residuals))
    color, grid_color, text_color, muted_text = "#2563eb", "#e5e5e5", "#1a1a1a", "#737373"
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=False)

    def draw_chart(axis, values, title, ylabel):
        for lag, value in enumerate(values):
            axis.bar(lag, value, width=0.4, color=color, edgecolor=color, linewidth=1.5, alpha=0.3 if lag == 0 else 1, zorder=3)
        axis.axhspan(-confidence_interval, confidence_interval, color=grid_color, alpha=0.5, zorder=1)
        axis.axhline(y=0, color=text_color, linewidth=0.8, zorder=2)
        axis.set_title(title, fontsize=13, fontweight="600", color=text_color, pad=12)
        axis.set_ylabel(ylabel, fontsize=11, color=text_color, labelpad=10)
        axis.set_xlim(-0.5, n_lags + 0.5)
        axis.set_ylim(-1.1, 1.1)
        axis.grid(True, color=grid_color, linewidth=1, alpha=0.5)
        axis.set_axisbelow(True)
        axis.tick_params(axis="both", labelsize=10, color=text_color)
        for lag, value in enumerate(values[1:], start=1):
            if abs(value) > confidence_interval:
                axis.annotate(f"lag {lag}", xy=(lag, value), xytext=(0, 8 if value > 0 else -12), textcoords="offset points", ha="center", va="bottom" if value > 0 else "top", fontsize=8, color=muted_text, fontweight="500")

    draw_chart(axes[0], acf_values, "Autocorrelation Function (ACF) of Residuals", "Autocorrelation")
    draw_chart(axes[1], pacf_values, "Partial Autocorrelation Function (PACF) of Residuals", "Partial Autocorrelation")
    axes[1].set_xlabel("Lag (months)", fontsize=11, color=text_color, labelpad=10)
    axes[0].annotate(f"{model_name} residuals | n={len(residuals)} | 95% CI: ±{confidence_interval:.3f}", xy=(0.02, 0.02), xycoords="axes fraction", ha="left", va="bottom", fontsize=9, color=muted_text, bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=grid_color, linewidth=0.5))
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"\nACF/PACF chart saved as '{output_path}'")
    if show:
        plt.show()
    else:
        plt.close(fig)
