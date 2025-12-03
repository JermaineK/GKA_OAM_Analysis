"""
Agent: quick roll-up + plots for GKA/OAM summaries.

Reads the CSV produced by gka_oam_extended.py, aggregates by folder,
and emits a few sanity plots (histograms / scatter) for fast inspection.
"""
from pathlib import Path
import matplotlib

# headless backend
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402


def main():
    root = Path(__file__).resolve().parent
    # prefer latest v4, then v3, then v2
    for name in ["gka_oam_image_summary_v4.csv", "gka_oam_image_summary_v3.csv", "gka_oam_image_summary_v2.csv"]:
        csv_path = root / name
        if csv_path.exists():
            break
    else:
        raise SystemExit("Missing summary CSV (v4/v3/v2). Run gka_oam_extended.py first.")

    df = pd.read_csv(csv_path)
    if "source" not in df.columns:
        raise SystemExit("CSV missing 'source' column; run gka_oam_extended.py first.")

    df["folder"] = df["source"].apply(lambda p: str(Path(p).parent))

    plot_dir = root / "gka_oam_plots_rollup"
    plot_dir.mkdir(parents=True, exist_ok=True)

    # m_star histogram
    plt.figure(figsize=(8, 4))
    counts = df["m_star"].value_counts().sort_index()
    counts.plot(kind="bar")
    plt.xlabel("m_star")
    plt.ylabel("count")
    plt.title("Dominant mode counts")
    plt.tight_layout()
    plt.savefig(plot_dir / "m_star_counts.png", dpi=150)
    plt.close()

    # purity vs coverage scatter, color by top folders
    top_folders = df["folder"].value_counts().index[:8]
    df["folder_top"] = df["folder"].where(df["folder"].isin(top_folders), "other")

    plt.figure(figsize=(7, 5))
    for folder, sub in df.groupby("folder_top"):
        plt.scatter(sub["coverage"], sub["purity"], s=18, alpha=0.7, label=folder)
    plt.xlabel("coverage")
    plt.ylabel("purity")
    plt.title("Purity vs coverage by folder")
    plt.legend(fontsize=8, loc="best", framealpha=0.6)
    plt.tight_layout()
    plt.savefig(plot_dir / "purity_vs_coverage.png", dpi=150)
    plt.close()

    # mean purity by folder (top 10)
    purity_by_folder = (
        df.groupby("folder")["purity"].mean().sort_values(ascending=False).head(10)
    )
    plt.figure(figsize=(9, 4))
    purity_by_folder.plot(kind="bar")
    plt.ylabel("mean purity")
    plt.title("Top folders by mean purity")
    plt.tight_layout()
    plt.savefig(plot_dir / "mean_purity_by_folder.png", dpi=150)
    plt.close()

    # coherence vs m_star scatter
    plt.figure(figsize=(7, 5))
    plt.scatter(df["m_star"], df["coherence_q"], s=20, alpha=0.7)
    plt.xlabel("m_star")
    plt.ylabel("coherence_q (p75/p25)")
    plt.title("Coherence vs dominant mode")
    plt.tight_layout()
    plt.savefig(plot_dir / "coherence_vs_mstar.png", dpi=150)
    plt.close()

    print(f"Plots written to {plot_dir}")


if __name__ == "__main__":
    main()
