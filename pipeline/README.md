# FRAP image-analysis pipeline (MATLAB)

MATLAB code that segments and tracks nuclei in Dorsal-GFP FRAP movies of
*Drosophila* embryos, extracts nuclear and cytoplasmic intensities, corrects
for unintentional bleaching, and fits nuclear import and export rates
(`k_in`, `k_out`).

## Requirements

- MATLAB (run with R2024b; R2020a or later should work), with the Image
  Processing, Optimization, and Statistics and Machine Learning toolboxes.
- Bio-Formats MATLAB toolbox (`bfmatlab`) version 5.5.1, which provides
  `bfopen` and `bfGetReader` for reading `.czi` files. It is not included
  here; download `bfmatlab.zip` for 5.5.1 from
  https://downloads.openmicroscopy.org/bio-formats/5.5.1/, unzip it and add
  the folder to the MATLAB path. Newer versions will probably work, but
  5.5.1 is the one used for the paper.

Add this folder and its subfolders to the MATLAB path before running:

```matlab
addpath(genpath('path\to\pipeline'))
```

## Data layout

Put the `.czi` movies in a folder called `Data` inside the folder you run
from, with one subfolder per imaging date:

```
Data/
  2022-10-11/Embryo1_dorsal_5x_1.8sec.czi
  2022-10-12/...
```

(To use another location, change `basepth2` at the top of
`script_gatherbleachfiles.m`.) `List of Dorsal FRAP data.xlsx` lists the
movies used in the paper.

## Running the pipeline

1. **`script_gatherbleachfiles.m`** finds all `.czi` files under `Data` and
   saves the list to `bleaching_filenames.mat` (and `Bleaching_stats.xlsx`).
2. **`script_bleachanalysis2.m`** loops over the files and, for each movie,
   calls `run_analyze_bleach` → `analyze_bleach_czi` → `fit_bleach_dl`.
   Run the movies through these two scripts rather than calling
   `run_analyze_bleach` directly: the file list made in step 1 is what the
   later steps use to find each movie's date folder.
   - For every movie MATLAB asks `Please enter bleached nuc label:`. The
     labelled nuclei of each frame are saved as images in
     `Image results/Tracking/<date>/<movie>/`; look up the label of the
     bleached nucleus there and type it in.
   - When the bleached nucleus disappears for a few frames after bleaching,
     the projection step (`reconstruct_missing_nucleus.m`) is switched on by
     uncommenting the marked block in `analyze_bleach_czi.m` (used for three
     movies in the paper).
3. Results are written to `czi mat results/<date>/<movie>/<movie>_results.mat`
   (one structure, `soln`), and images of the masks and labels to
   `Image results/`.
4. **`generate_boxplots.m`** draws the box plots of import and export rates
   on the dorsal and ventral sides from all files in `czi mat results`.
5. **`gather_cyto_nuc_ints.m`** makes the simulated recovery curves and N/C
   ratio box plots. The figures in the paper (Figs. 10-14, S1 and S2 Figs)
   are made by the scripts in [`../figure_scripts/`](../figure_scripts/).

## Background subtraction

After the nuclear and cytoplasmic intensities are first extracted,
`analyze_bleach_czi.m` subtracts the detector's zero-light level (the value
recorded where no light reaches the detector) from every Dl-GFP image and
extracts the intensities again:

- Dorsal side (mean nuclear Dl-GFP below mean cytoplasmic Dl-GFP): many
  pixels receive no light, so the level is measured in every frame with
  `Functions/ftn_calcbg.m` (centre of the most populated peak of the
  intensity histogram, by a Gaussian fit).
- Ventral side: the Dl-GFP image has no unlit pixels, so the median of the
  dorsal per-movie values, 103.3 counts (`bg_ventral`), is used.
- The subtraction is done in double precision, without clipping at zero.

Segmentation, tracking, normalization and fitting are unchanged.

## What the code does

| Step | Main file(s) |
|---|---|
| Read `.czi`, time stamps | `openczi.m`, `analyze_bleach_czi.m` |
| Pre-processing, marker-controlled watershed, Otsu threshold per region, mask cleaning | `find_nucleizoom_new_comp_set2.m` |
| Tracking by point registration with lookback | `track_nuclei_labels_exact.m` |
| Projection of a bleached nucleus that disappears | `reconstruct_missing_nucleus.m` |
| Nuclear and cytoplasmic intensities | `dl_vs_t.m` |
| Background (zero-light level) subtraction | `analyze_bleach_czi.m`, `Functions/ftn_calcbg.m` |
| Nearest neighbours and normalization for unintentional bleaching | `persistent_closest_neighbors.m`, `plot_neighbors.m`, `fit_bleach_dl.m` |
| Model and least-squares fit of `k_in`, `k_out` with 95% confidence intervals | `ftn_FRAP.m`, `fit_bleach_dl.m` |
| Helper functions | `Functions/`, `fraccircshift.m`, `traceobject.m` |

## Output fields (`soln`)

| Field | Meaning |
|---|---|
| `side` | `'dorsal'` if the mean nuclear Dl-GFP intensity is below the mean cytoplasmic intensity, otherwise `'ventral'` (the grouping used for the box plots) |
| `t` | frame times (s), bleaching at `t = 0` |
| `y_nuc`, `y_cyt` | intensity of the bleached nucleus and of its cytoplasmic compartment |
| `y_nucavg`, `y_cytavg`, `y_nucall`, `y_cytall` | intensities of the other nuclei / cytoplasm |
| `idx_bleach`, `locs` | label of the bleached nucleus and its nearest neighbours |
| `k_in`, `k_out`, `c0`, `r2` | fitted import and export rate constants (1/s), initial value, R² |
| `delta1`, `delta2`, `covar` | 68% and 95% confidence half-widths of `[c0, k_in, k_out]`, covariance |

## Other folders

- `legacy/ftn_imread.m`: an older reader for `.lsm` files. It is not used by
  the pipeline and needs functions that are not included.
