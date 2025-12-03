"""
Compute a compact invariant table from GKA summary CSVs.

For each summary CSV in Results/summaries, aggregate a few statistics:
- median/90th pct parity_odd
- median eta_slope, eta_knee, knee/radius
- median borgt_slope, purity, coverage
- counts (total, valid rows with finite eta_slope/eta_knee)

Outputs: Results/summaries/invariant_meta.csv
"""
from pathlib import Path
import numpy as np
import pandas as pd


def summarize_one(csv_path: Path) -> dict:
    df = pd.read_csv(csv_path)
    out = {"dataset": csv_path.name}

    # Total counts
    out["n_total"] = len(df)

    # Valid mask for eta-based stats
    mask = (
        df["eta_slope"].replace([np.inf, -np.inf], np.nan).notna()
        & df["eta_knee"].replace([np.inf, -np.inf], np.nan).notna()
        & df["radius"].replace([np.inf, -np.inf], np.nan).notna()
    )
    out["n_valid"] = int(mask.sum())

    def med(series, m=None):
        if m is not None:
            series = series[m]
        series = series.replace([np.inf, -np.inf], np.nan).dropna()
        return float(series.median()) if len(series) else np.nan

    def pct(series, q, m=None):
        if m is not None:
            series = series[m]
        series = series.replace([np.inf, -np.inf], np.nan).dropna()
        return float(np.percentile(series, q)) if len(series) else np.nan

    out["parity_odd_med"] = med(df["parity_odd"])
    out["parity_odd_p90"] = pct(df["parity_odd"], 90)
    out["eta_slope_med"] = med(df["eta_slope"], mask)
    out["eta_knee_med"] = med(df["eta_knee"], mask)
    out["knee_over_radius_med"] = med(df["eta_knee"] / df["radius"], mask)
    out["borgt_slope_med"] = med(df["borgt_slope"])
    out["purity_med"] = med(df["purity"])
    out["coverage_med"] = med(df["coverage"])
    return out


def main():
    root = Path("Results/summaries")
    if not root.exists():
        raise SystemExit(f"Missing summaries folder: {root}")
    rows = []
    for csv in sorted(root.glob("gka_oam_image_summary_*.csv")):
        rows.append(summarize_one(csv))
    if not rows:
        raise SystemExit("No summary CSVs found in Results/summaries")
    out = pd.DataFrame(rows)
    out_path = root / "invariant_meta.csv"
    out.to_csv(out_path, index=False)
    print(f"Wrote {out_path}")
    print(out)


if __name__ == "__main__":
    main()
