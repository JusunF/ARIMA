library(readr)
library(dplyr)
library(zoo)
library(ggplot2)
library(tidyr)

# ---------------------------------------------------------------------------
# MCO known-intervention imputation - CONTEXT / EXHIBIT ONLY.
# Not used by the finalized ARIMA(1,1,0) model. Shows what LRT Ampang
# ridership would likely have looked like without the MCO lockdown shock:
#   1. Linear interpolation across the MCO window (initial fill).
#   2. STL decomposition on that interpolated series; the MCO window is
#      replaced with STL trend + seasonal (dropping the remainder).
# ---------------------------------------------------------------------------

csv_path    <- "ridership_headline.csv"
start_date  <- as.Date("2019-01-01")
mco_start   <- as.Date("2020-03-01")
mco_end     <- as.Date("2021-12-01")
output_path <- "mco_imputation_comparison.png"

# --- 1. Load full history and aggregate to monthly totals -----------------
df <- read_csv(csv_path, show_col_types = FALSE)

df_ampang <- df %>%
  mutate(date = as.Date(date)) %>%
  filter(date >= start_date) %>%
  arrange(date)

monthly <- df_ampang %>%
  mutate(month = as.Date(format(date, "%Y-%m-01"))) %>%
  group_by(month) %>%
  summarise(value = sum(rail_lrt_ampang, na.rm = TRUE), .groups = "drop") %>%
  arrange(month)

expected_months <- seq(min(monthly$month), max(monthly$month), by = "month")
stopifnot(all(expected_months == monthly$month))

# --- 2. Mask the MCO window as NA ------------------------------------------
monthly <- monthly %>%
  mutate(
    is_mco = month >= mco_start & month <= mco_end,
    value_raw = value,
    value_masked = ifelse(is_mco, NA, value)
  )

# --- 3. Step 1: linear interpolation across the masked gap -----------------
monthly <- monthly %>%
  mutate(value_linear = na.approx(value_masked, x = month, na.rm = FALSE))

# --- 4. Step 2: STL on the interpolated series, reconstruct MCO window
#        from trend + seasonal only -----------------------------------------
value_ts <- ts(
  monthly$value_linear,
  start = c(as.integer(format(min(monthly$month), "%Y")),
            as.integer(format(min(monthly$month), "%m"))),
  frequency = 12
)

stl_fit <- stl(value_ts, s.window = "periodic", robust = TRUE)
stl_components <- as.data.frame(stl_fit$time.series)

monthly <- monthly %>%
  mutate(
    stl_trend    = stl_components$trend,
    stl_seasonal = stl_components$seasonal,
    stl_reconstructed = stl_trend + stl_seasonal,
    value_imputed = ifelse(is_mco, stl_reconstructed, value_raw)
  )

cat("Imputed values for MCO window (", format(mco_start, "%b %Y"), " to ",
    format(mco_end, "%b %Y"), "):\n", sep = "")
print(
  monthly %>%
    filter(is_mco) %>%
    select(month, value_raw, value_linear, value_imputed)
)

# --- 5. Comparison chart ----------------------------------------------------
plot_df <- monthly %>%
  select(month, value_raw, value_linear, value_imputed) %>%
  pivot_longer(
    cols = c(value_raw, value_linear, value_imputed),
    names_to = "series",
    values_to = "value"
  ) %>%
  mutate(
    series = recode(
      series,
      value_raw     = "Actual (with MCO crash)",
      value_linear  = "Linear interpolation only",
      value_imputed = "STL-adjusted imputation"
    )
  )

p <- ggplot(plot_df, aes(x = month, y = value, color = series, linetype = series)) +
  annotate(
    "rect",
    xmin = mco_start, xmax = mco_end,
    ymin = -Inf, ymax = Inf,
    fill = "grey85", alpha = 0.5
  ) +
  annotate(
    "text",
    x = mco_start + 150, y = max(monthly$value_raw, na.rm = TRUE) * 0.95,
    label = "MCO Lockdown Period", size = 3, color = "grey40"
  ) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = c(
    "Actual (with MCO crash)"    = "#737373",
    "Linear interpolation only"  = "#f59e0b",
    "STL-adjusted imputation"    = "#2563eb"
  )) +
  scale_linetype_manual(values = c(
    "Actual (with MCO crash)"    = "solid",
    "Linear interpolation only"  = "dashed",
    "STL-adjusted imputation"    = "solid"
  )) +
  scale_y_continuous(labels = function(x) paste0(round(x / 1e6, 1), "M")) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title = "LRT Ampang Line Ridership - MCO Period, Known-Intervention Imputation",
    subtitle = "Context exhibit only - not used in the finalized ARIMA(1,1,0) model",
    x = NULL, y = "Monthly Ridership (millions)", color = NULL, linetype = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "grey40", size = 10),
    legend.position = "top",
    legend.justification = "left",
    panel.grid.minor = element_blank()
  )

ggsave(output_path, plot = p, width = 12, height = 6, dpi = 150, bg = "white")
cat("\nComparison chart saved as:", output_path, "\n")

print(p)