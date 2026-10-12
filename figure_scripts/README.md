# Figure scripts

Python scripts that make the data figures of the paper (Figs 10–14, S1 and S2 Figs) in the PLOS ONE format:
`<name>.tif` (RGB, 300 dpi, LZW compression, 7.5 in or 3.5 in wide) and `<name>_editable.pdf` (vector, live text).
Figs 1–9 are diagrams and micrographs prepared in Adobe Illustrator and are not made here.

## Requirements

Python 3.9 or later with numpy, scipy, pandas, matplotlib, pillow and h5py
(`pip install -r requirements.txt`, or an Anaconda installation, which has all of them).
The text is set in Arial. On Windows it is found automatically; elsewhere, copy `arial.ttf`, `arialbd.ttf`,
`ariali.ttf` and `arialbi.ttf` into a folder `fonts/` next to the scripts.

## Inputs

| Figure | Input |
|---|---|
| Fig 10, Fig 11, Fig 12 | the background-corrected result files: the unzipped result zip (mat data) from the Texas Data Repository, i.e. the folder that contains `czi mat results/` |
| Fig 13 | the ground-truth data: `gt_embryos.csv`, `stacks/` (annotations and the masks of CellProfiler, Cellpose 3 and StarDist) and `pipeline_masks/` (our pipeline), as used by `ground_truth/FRAP_GroundTruth_segmentation_comparison/` |
| Fig 14, S1 Fig | tables in `figure_data/` written by the ground-truth package (`score_gt_all.py`, `frap_track.py`): `GT_AP_vs_IoU_threshold.csv`, `GT_matched_nuclei.csv`, `GT_summary.csv`, `GT_bleached_track_summary.csv` |
| S2 Fig | `figure_data/S2_frames/*.npz`: the 9 example frames (both channels, our pipeline's labels and the other tool's labels), taken from the outputs of the tool-comparison packages in `comparisons/` |

## Running

```
python make_figures.py --results "<folder with czi mat results>" --gt "<ground-truth data folder>" --out figures
```

`--only Fig10 Fig14 ...` makes only some of the figures. All seven figures take about 1–2 minutes.

Some Pillow builds (seen with Anaconda on Windows) crash when writing LZW-compressed TIFF files. In that case run
with `--tiff png`, which writes lossless PNG files at the final size, and convert them in MATLAB with
`png_to_tiff('figures')` (LZW TIFF, 300 dpi).

## Outputs

- `figures/Fig10.tif`, `Fig11.tif`, ..., `S2_Fig.tif` and the matching `*_editable.pdf` files.
- The numbers shown in each figure, written to `figure_data/`:

| File | Content |
|---|---|
| `Fig10_data.csv` | per side and time point: mean and SD of the 6 neighbouring nuclei and compartments, the bleached nucleus and compartment, and both after normalization |
| `Fig11_data.csv` | per movie: side, fitted `k_in`, `k_out` (1/s), their 95% confidence half-widths, R² |
| `Fig12ab_data.csv` | per side and time point: mean simulated nuclear intensity with its 95% band, mean normalized cytoplasmic intensity |
| `Fig12c_data.csv` | per movie: N/C ratio at the end of the simulation |
| `Fig13_panels.csv`, `S2_Fig_panels.csv` | the embryo, frame and failure shown in each panel |
| `Fig14c_data.csv` | per method: % of the 523 post-bleach frames in which the bleached nucleus was segmented and tracked |

## Files

| File | What it does |
|---|---|
| `make_figures.py` | runs all figure scripts |
| `fig10_intensities.py` | Fig 10 |
| `fig11_fig12_rates.py` | Figs 11 and 12 (the simulation solves the model with each movie's fitted rates) |
| `fig13_examples_vs_ground_truth.py` | Fig 13 |
| `fig14_comparison_quantitative.py` | Fig 14 |
| `s1_fig_iou_histograms.py` | S1 Fig |
| `s2_fig_examples_vs_pipeline.py` | S2 Fig |
| `gt_tools.py` | ground-truth reading and matching (the same functions as in the ground-truth package) |
| `tiles.py` | drawing of the 12 × 12 µm example crops |
| `plos_style.py`, `mat_io.py` | Arial and output format; reading of MATLAB files |
| `png_to_tiff.m` | MATLAB: converts the PNG output of `--tiff png` into LZW TIFF files |
