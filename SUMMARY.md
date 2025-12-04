# Findings and Evidence (GKA_OAM_Analysis)

This file aggregates key observations from the processed datasets to keep reproducible evidence in one place.

## DM3 image analyses (GKA pipeline)
- Datasets: `Additional data/converted_png` (log-spiral), `Additional data/converted_png_5761427` (probe/diffraction/ZL carbon).
- Parity/odd channel:
  - Log-spiral set: many frames show strong parity_odd (e.g., >0.7 for sorter-like outputs), consistent with boundary–odd activation.
  - Carbon spectra: parity_odd ranges ~0.35–0.88; diffraction frames lower, spectra higher (e.g., `002_ZL_carbon-20s+16eV.png` parity_odd≈0.83).
- Eta slope/knee:
  - Knees span single-digit to few-hundred pixels; e.g., ZL spectra yield eta_knee ≈ 6–50 px, diffraction frames much larger.
  - Eta_slope varies by modality (diffraction near zero/NaN; spectra positive).
- Summaries:
  - log-spiral summary: `gka_oam_image_summary_v4.csv` (root).
  - 5761427 batch summary: `Additional data/converted_png_5761427/gka_oam_image_summary_v4_5761427.csv`.

## Excel (Figure2/3/4) findings
- Source: `Additional data/Figure2/3/4_data.xlsx` and CSV derivatives.
- L-dependence:
  - Fig.2(g): transmission rises with L (0→5), VRR/2π drops (17.4→5.8 MHz) — reproducible OAM scaling.
  - Fig.3: extinction decreases with L (mean: L0≈0.306 > L1≈0.222 > L2≈0.165; peak at Nt=8 shrinks with L).
  - Fig.4: L=2 retains more integrated energy than Gaussian across all Ng; consistent ordering.
- Visuals: `Additional data/plots_excel/` (detuning vs transmission, extinction vs Nt, time traces).

## Entanglement .dat (NSDI) findings
- Parsed tables: `Additional data/entanglement_parsed/`
  - Fig2 Panel d–f: state populations per target (Magnesium, Beryllium, Neon0) with lambda/I/Keldysh; see `fig2_panel_d_f_states.csv`.
  - Fig3 Panel a–f: momentum/log-negativity grids saved as CSVs plus summary (`fig3_matrices_summary.csv`).
- Visuals: `Additional data/plots_newdata/` (state population bars; Panel_a–f heatmaps).

## SPP/GRIN text grids (3968750_extracted)
- Files: Step/SPP, Smooth/SPP, GRIN, DeMPlex/MPlex grids; visual heatmaps in `Additional data/plots_newdata/`.
- MAT files (`T_*.mat`) remain opaque (MATLAB objects); open in MATLAB if needed.

## How to regenerate
- DM3 conversion: see README DM3 section; run `gka_oam_extended.py` on converted PNGs.
- Excel plots: `visualize_excel_data.py`.
- Parsed new data plots: `visualize_parsed_newdata.py`.
- Invariant table: `compute_invariant_table.py` (writes `Results/summaries/invariant_meta.csv`).
- Symmetry filters: see `symmetry_filters.py` (parity-odd, spiral-harmonic, fixed log-spiral kernel, slow-tick projection).

## Evidence snapshot (high level)
- Boundary–odd / parity channel is active in multiple independent datasets (log-spiral DM3, carbon spectra).
- OAM index L drives measurable scalings (transmission/VRR, extinction, energy retention).
- Additional physics datasets (entanglement, SPP/GRIN masks) show structured, non-random patterns; visuals available for inspection.
