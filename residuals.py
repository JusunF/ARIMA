"""Residuals diagnostic chart: time plot, histogram, and Q-Q plot."""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats


def generate_residuals_chart(residuals, model_name="ARIMA(1,1,0)", output_path="arima_residuals_diagnostic.png", show=True):
    """Create and save a 3-panel residuals diagnostic chart: time plot, histogram, Q-Q plot."""
    residuals = np.asarray(residuals)
    color, grid_color, text_color, muted_text = "#2563eb", "#e5e5e5", "#1a1a1a", "#737373"

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Panel 1: Residuals over time
    ax1 = axes[0]
    ax1.plot(residuals, color=color, linewidth=1.2, zorder=3)
    ax1.scatter(range(len(residuals)), residuals, color=color, s=15, zorder=4, alpha=0.6)
    ax1.axhline(y=0, color=text_color, linewidth=0.8, linestyle="--", zorder=2)
    std = residuals.std()
    ax1.axhspan(-2 * std, 2 * std, color=grid_color, alpha=0.4, zorder=1)
    ax1.set_title("Residuals Over Time", fontsize=12, fontweight="600", color=text_color, pad=10)
    ax1.set_xlabel("Time Index", fontsize=10, color=text_color)
    ax1.set_ylabel("Residual", fontsize=10, color=text_color)
    ax1.grid(True, color=grid_color, linewidth=1, alpha=0.5)
    ax1.set_axisbelow(True)
    ax1.tick_params(labelsize=9, color=text_color)

    # Panel 2: Histogram with normal curve overlay
    ax2 = axes[1]
    n, bins, patches = ax2.hist(residuals, bins=12, color=color, alpha=0.7, edgecolor="white", zorder=3, density=True)
    x = np.linspace(residuals.min(), residuals.max(), 100)
    mu, sigma = residuals.mean(), residuals.std()
    ax2.plot(x, stats.norm.pdf(x, mu, sigma), color="#ea580c", linewidth=2, label="Normal fit", zorder=4)
    ax2.set_title("Residual Distribution", fontsize=12, fontweight="600", color=text_color, pad=10)
    ax2.set_xlabel("Residual", fontsize=10, color=text_color)
    ax2.set_ylabel("Density", fontsize=10, color=text_color)
    ax2.grid(True, color=grid_color, linewidth=1, alpha=0.5)
    ax2.set_axisbelow(True)
    ax2.tick_params(labelsize=9, color=text_color)
    ax2.legend(loc="upper right", fontsize=9, frameon=True, facecolor="white", edgecolor=grid_color)

    # Panel 3: Q-Q plot
    ax3 = axes[2]
    (osm, osr), (slope, intercept, r) = stats.probplot(residuals, dist="norm")
    ax3.scatter(osm, osr, color=color, s=20, zorder=3, alpha=0.7)
    ax3.plot(osm, slope * osm + intercept, color="#ea580c", linewidth=1.5, zorder=4)
    ax3.set_title(f"Normal Q-Q Plot (R²={r**2:.3f})", fontsize=12, fontweight="600", color=text_color, pad=10)
    ax3.set_xlabel("Theoretical Quantiles", fontsize=10, color=text_color)
    ax3.set_ylabel("Sample Quantiles", fontsize=10, color=text_color)
    ax3.grid(True, color=grid_color, linewidth=1, alpha=0.5)
    ax3.set_axisbelow(True)
    ax3.tick_params(labelsize=9, color=text_color)

    fig.suptitle(f"{model_name} Residual Diagnostics", fontsize=14, fontweight="600", color=text_color, y=1.02)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"\nResiduals diagnostic chart saved as '{output_path}'")
    if show:
        plt.show()
    else:
        plt.close(fig)