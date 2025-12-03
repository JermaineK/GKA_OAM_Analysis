"""
Visuals for newly parsed datasets:
- Entanglement .dat (Figure2/3) parsed CSVs.
- SPP/GRIN/MPlex text grids (3968750_extracted).
Outputs PNGs under Additional data/plots_newdata/.
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


def ensure_dir(path: Path):
    path.mkdir(parents=True, exist_ok=True)


def plot_entanglement_fig2_states(base: Path, out_dir: Path):
    """Bar plots of state populations per target for Fig2 Panel d-f."""
    csv_path = base / "fig2_panel_d_f_states.csv"
    if not csv_path.exists():
        return
    df = pd.read_csv(csv_path)
    for target, sub in df.groupby("Target"):
        plt.figure(figsize=(10, 4))
        plt.bar(sub["state"], sub["population"])
        plt.xticks(rotation=90, fontsize=7)
        plt.ylabel("Population (arb.)")
        plt.title(f"Fig2 d-f states: {target} (λ={sub['lambda_nm'].iloc[0]:.1f} nm, I={sub['I_e14_Wcm2'].iloc[0]:.3f}e14)")
        plt.tight_layout()
        plt.savefig(out_dir / f"fig2_states_{target}.png", dpi=150)
        plt.close()


def plot_entanglement_fig3_panels(base: Path, out_dir: Path):
    """Heatmaps for Fig3 Panel_a..f CSV numeric grids (assumed columns: x, y, value)."""
    for csv in sorted(base.glob("Panel_*.csv")):
        df = pd.read_csv(csv, header=None)
        if df.shape[1] == 3:
            x = df.iloc[:, 0].values
            y = df.iloc[:, 1].values
            z = df.iloc[:, 2].values
            plt.figure(figsize=(5, 4))
            sc = plt.scatter(x, y, c=z, s=8, cmap="viridis")
            plt.colorbar(sc, label="value")
            plt.xlabel("x")
            plt.ylabel("y")
            plt.title(csv.stem)
            plt.tight_layout()
            plt.savefig(out_dir / f"{csv.stem}.png", dpi=150)
            plt.close()
        else:
            # fallback: matrix heatmap
            mat = df.apply(pd.to_numeric, errors="coerce").values
            plt.figure(figsize=(5, 4))
            plt.imshow(mat, aspect="auto", cmap="viridis")
            plt.colorbar(label="value")
            plt.title(csv.stem)
            plt.tight_layout()
            plt.savefig(out_dir / f"{csv.stem}.png", dpi=150)
            plt.close()


def plot_grids_spp(root: Path, out_dir: Path):
    """Plot grid-like TXT datasets from 3968750_extracted."""
    txt_files = [
        "StepSPP.txt",
        "StepSPP_lens.txt",
        "SmoothSPP.txt",
        "SmoothSPP_lens.txt",
        "GRINSPP.txt",
        "GRINSPP_lens.txt",
        "DeMPlex_Mode0.txt",
        "MPlex_Mode2.txt",
    ]
    for name in txt_files:
        path = root / name
        if not path.exists():
            continue
        df = pd.read_csv(path, sep=r"\\s+|\\t", engine="python", header=None)
        mat = df.apply(pd.to_numeric, errors="coerce").values
        plt.figure(figsize=(7, 4))
        plt.imshow(mat, aspect="auto", cmap="coolwarm")
        plt.colorbar(label="value")
        plt.title(name)
        plt.tight_layout()
        plt.savefig(out_dir / f"{name}.png", dpi=150)
        plt.close()


def main():
    out_base = Path("Additional data/plots_newdata")
    ensure_dir(out_base)

    # Entanglement parsed
    ent_base = Path("Additional data/entanglement_parsed")
    if ent_base.exists():
        plot_entanglement_fig2_states(ent_base, out_base)
        plot_entanglement_fig3_panels(ent_base, out_base)

    # SPP/GRIN grids
    spp_root = Path("3968750_extracted")
    if spp_root.exists():
        plot_grids_spp(spp_root, out_base)

    print(f"Plots written to {out_base}")


if __name__ == "__main__":
    main()
