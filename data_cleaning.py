import pandas as pd

# ---------------------------------------------------------------------------
# Clean raw ridership_headline.csv (DAILY granularity) down to a clean
# MONTHLY LRT Ampang Line series, Jan 2022 - Jun 2026.
#
# Why this is needed:
#   ridership_headline.csv has one row PER DAY, not per month. Passing the
#   filtered daily rows straight into a time-series model would incorrectly
#   treat each DAY as a monthly observation. Therefore, daily values are
#   aggregated into true monthly totals first.
# ---------------------------------------------------------------------------

# Read raw daily dataset
raw = pd.read_csv("ridership_headline.csv")

# Convert date column to datetime
raw["date"] = pd.to_datetime(raw["date"], errors="coerce")

# Filter from 2022-01-01 and remove missing LRT Ampang values
monthly_ampang = (
    raw[
        (raw["date"] >= "2019-01-01") &
        (raw["rail_lrt_ampang"].notna())
    ]
    .assign(
        month=lambda x: x["date"].dt.to_period("M")
    )
    .groupby("month", as_index=False)
    .agg(
        rail_lrt_ampang=("rail_lrt_ampang", "sum"),
        n_days=("rail_lrt_ampang", "count")
    )
    .sort_values("month")
)

# Convert month to YYYY-MM-01
monthly_ampang["month"] = (
    monthly_ampang["month"]
    .dt.to_timestamp()
)

# Sanity checks -------------------------------------------------------------

print("Number of months:", len(monthly_ampang))
print("First month:", monthly_ampang["month"].iloc[0].strftime("%Y-%m-%d"))
print("Last month:", monthly_ampang["month"].iloc[-1].strftime("%Y-%m-%d"))

# Flag months with fewer than 28 days of data
incomplete = monthly_ampang[
    monthly_ampang["n_days"] < 28
]

if len(incomplete) > 0:
    print("\nWARNING: months with fewer than 28 days of data:")
    print(incomplete)
else:
    print("All months have full daily coverage. No missing values.")

# Write cleaned monthly CSV -----------------------------------------------

monthly_out = (
    monthly_ampang[
        ["month", "rail_lrt_ampang"]
    ]
    .rename(columns={"month": "date"})
)

monthly_out["date"] = monthly_out["date"].dt.strftime("%Y-%m-%d")

monthly_out.to_csv(
    "ridership_ampang_monthly.csv",
    index=False
)

print(
    "\nCleaned monthly file written to "
    "ridership_ampang_monthly.csv"
)