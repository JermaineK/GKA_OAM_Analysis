"""
Quick visuals for the Excel-derived OAM datasets (Fig.2/3/4).
Outputs PNGs to Additional data/plots_excel.
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402


def plot_fig2(base: Path, out_dir: Path):
    xls = pd.ExcelFile(base / "Figure2_data.xlsx")
    sheets = [s for s in xls.sheet_names if s.startswith("Fig.2(") and s != "Fig.2(g)"]
    rates = [
        "Transmission under 0.2 us^-1",
        "Transmission under 0.8 us^-1",
        "Transmission under 3.4 us^-1",
        "Transmission under 7.1us^-1",
    ]
    colors = ["C0", "C1", "C2", "C3"]

    for sh in sheets:
        df = xls.parse(sh)
        if "Detuning\\2pi * MHz" in df.columns:
            x = df["Detuning\\2pi * MHz"]
        else:
            x = df.iloc[:, 0]
        plt.figure(figsize=(7, 4))
        for r, c in zip(rates, colors):
            if r not in df:
                continue
            plt.plot(x, df[r], label=r, color=c, alpha=0.9)
        plt.xlabel("Detuning / 2π (MHz)")
        plt.ylabel("Transmission")
        plt.title(sh)
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(out_dir / f"fig2_{sh}.png", dpi=150)
        plt.close()

    # Fig 2(g)
    if "Fig.2(g)" in xls.sheet_names:
        df = xls.parse("Fig.2(g)")
        plt.figure(figsize=(6, 4))
        plt.errorbar(df["L"], df["Transmission"], yerr=df["Transmission error bar"], fmt="o-", label="Transmission")
        plt.xlabel("L")
        plt.ylabel("Transmission")
        plt.twinx()
        plt.plot(df["L"], df["VRR / 2pi * MHz"], "s--", color="C3", label="VRR / 2π")
        plt.ylabel("VRR / 2π (MHz)")
        plt.title("Fig.2(g): Transmission and VRR vs L")
        plt.tight_layout()
        plt.savefig(out_dir / "fig2_g.png", dpi=150)
        plt.close()


def plot_fig3(base: Path, out_dir: Path):
    xls = pd.ExcelFile(base / "Figure3_data.xlsx")
    plt.figure(figsize=(6, 4))
    for L, color in zip(["L=0", "L=1", "L=2"], ["C0", "C1", "C2"]):
        d = xls.parse(L)
        plt.errorbar(d["Nt"], d["extinction"], yerr=d["extinction error bar"], fmt="o-", label=L, color=color)
    plt.xlabel("Nt")
    plt.ylabel("Extinction")
    plt.title("Fig.3: Extinction vs Nt")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_dir / "fig3_extinction.png", dpi=150)
    plt.close()


def plot_fig4(base: Path, out_dir: Path):
    xls = pd.ExcelFile(base / "Figure4_data.xlsx")
    for sh in ["L=2", "Gaussian"]:
        d = xls.parse(sh)
        time = d["time bins"]
        plt.figure(figsize=(7, 4))
        for col, color in zip(["Ng=0", "Ng=0.3", "Ng=1.0", "Ng=3.0"], ["C0", "C1", "C2", "C3"]):
            plt.plot(time, d[col], label=col, color=color)
        plt.xlabel("Time bins")
        plt.ylabel("Signal")
        plt.title(f"Fig.4: {sh}")
        plt.legend()
        plt.tight_layout()
        plt.savefig(out_dir / f"fig4_{sh.replace('=', '').replace(' ', '_')}.png", dpi=150)
        plt.close()


def main():
    base = Path("Additional data")
    out_dir = base / "plots_excel"
    out_dir.mkdir(exist_ok=True, parents=True)

    plot_fig2(base, out_dir)
    plot_fig3(base, out_dir)
    plot_fig4(base, out_dir)
    print(f"Plots written to {out_dir}")


if __name__ == "__main__":
    main()
