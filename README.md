# FRAP image analysis: measuring nuclear import and export rates

Code for the paper:

S. V. Hiremath, E. Goyal, G. T. Reeves, C. M. Williams. *FRAP analysis:
Measuring biophysical kinetic parameters using image analysis.* Manuscript
under review.

The pipeline segments and tracks nuclei in FRAP movies of Dorsal-GFP in
*Drosophila* embryos, follows the photobleached nucleus even when it
disappears, and fits nuclear import and export rates. The comparison
packages benchmark its segmentation and tracking against CellProfiler,
Cellpose 3 and StarDist, and against ground-truth annotations.

## Contents

| Folder | What it contains | Language |
|---|---|---|
| [`pipeline/`](pipeline/) | The FRAP image-analysis pipeline: segmentation, tracking, intensity extraction, normalization and fitting | MATLAB |
| [`comparisons/FRAP_CellProfiler_segmentation_comparison/`](comparisons/FRAP_CellProfiler_segmentation_comparison/) | Our segmentation compared with CellProfiler 4.2.8 on 10 movies | MATLAB, CellProfiler |
| [`comparisons/FRAP_Cellpose_segmentation_comparison/`](comparisons/FRAP_Cellpose_segmentation_comparison/) | Our segmentation compared with Cellpose 3 | MATLAB, Python |
| [`comparisons/FRAP_StarDist_segmentation_comparison/`](comparisons/FRAP_StarDist_segmentation_comparison/) | Our segmentation compared with StarDist | MATLAB, Python |
| [`ground_truth/FRAP_GroundTruth_segmentation_comparison/`](ground_truth/FRAP_GroundTruth_segmentation_comparison/) | All four methods scored against Labkit ground truth on 10 further movies, plus bleached-nucleus and tracking analysis | MATLAB, Python |
| [`figure_scripts/`](figure_scripts/) | The data figures of the paper (Figs 10–14, S1 and S2 Figs) and the numbers behind them (`figure_data/`) | Python |

Each folder has its own README with requirements and step-by-step
instructions. The comparison packages each contain a copy of the pipeline's
segmentation and tracking functions (`code/frap_pipeline/`) so that they run
on their own; these copies are the same as in `pipeline/` apart from
display settings in `imshowcat.m`.

## Data

- Texas Data Repository, https://doi.org/10.18738/T8/M3TGCY: the image data (`.czi` movies), the
  fitted results per movie (`*_results.mat`, background-corrected) and the ground-truth masks of the
  10 annotated movies.

The intensities are background-corrected: the detector's zero-light level is subtracted from the
Dl-GFP images before the intensities are measured (see `pipeline/README.md`).

Raw images, masks and results are not stored in this repository.

## Licence

The code is released under the MIT licence (see [`LICENSE`](LICENSE)),
copyright (c) 2026 Sharva V. Hiremath, Etika Goyal, Gregory T. Reeves and
Cranos M. Williams.

Exception: `fraccircshift.m` (in `pipeline/` and in the ground-truth package)
is from the MATLAB File Exchange (Bernhard, 2014) and keeps its own
BSD-3-Clause licence, given at the top of the file. Bio-Formats, CellProfiler,
Cellpose and StarDist are not included here and have their own licences.
