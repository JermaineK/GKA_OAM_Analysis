# GKA_OAM_Analysis

Analysis scripts and files for paper reproducibility. Standalone copy of the OAM/Talbot analysis assets (images + Python scripts).

## Contents
- `gka_oam_extended.py`: headless analyzer for images and .bin files (parity, coherence, BORGT slope, etc.).
- `gka_oam_visualize.py`: roll-up plots (histograms/scatter) from the summary CSV.
- `gka_oam_image_summary_v3.csv`: latest metrics table from the last run.
- Raw images under `Generalized angle-OAM Talbot effect/` and `OAM sorting/` (plus prior outputs/plots).

## Usage
```bash
# activate your venv (if you have one)
# python -m venv .venv && .venv\Scripts\activate  # optional
pip install -r requirements.txt  # see below for minimal list

# Run analysis (images; add --plots to save per-image I(theta))
python gka_oam_extended.py --root . --out-csv gka_oam_image_summary_v3.csv

# Example for BINs (if present)
python gka_oam_extended.py --root . \
  --bin-pattern "*.bin" --bin-shape 1024,1024 --bin-dtype uint16 \
  --out-csv gka_oam_with_bins.csv

# Make roll-up plots
python gka_oam_visualize.py
```

## Minimal Python deps
```
numpy
pandas
matplotlib
opencv-python-headless
```
SciPy is only needed if you later add the MAT-reader script.

## Notes
- BORGT slope is guarded with intensity thresholds and outlier filtering; see `borgt_status` for fit quality.
- Parity/oddness is often the cleanest discriminator between petals/Talbot vs sorter outputs.
- This is a straight copy from `C:\Weather warning project\Data and code` for a dedicated repo.

## DM3 (Gatan) conversion
- DM3 files under `Additional data/` can be converted to PNG via `hyperspy` + OpenCV. Example inline script:
  ```
  import hyperspy.api as hs, cv2, numpy as np
  from pathlib import Path
  root = Path("Additional data"); out = root/"converted_png"; out.mkdir(exist_ok=True)
  for dm3 in root.glob("*.dm3"):
      data = np.asarray(hs.load(dm3).data)
      if data.ndim > 2: data = data.reshape((-1, data.shape[-2], data.shape[-1]))[0]
      dmin, dmax = data.min(), data.max(); span = dmax - dmin if dmax != dmin else 1
      img = ((data - dmin)/span * 65535).clip(0, 65535).astype(np.uint16)
      cv2.imwrite(str(out/(dm3.stem+".png")), img)
  ```
- After conversion, run `gka_oam_extended.py --root Additional data/converted_png --out-csv gka_oam_image_summary_v4_additional.csv`.
