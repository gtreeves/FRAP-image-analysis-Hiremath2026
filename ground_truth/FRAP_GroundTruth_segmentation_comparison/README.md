# Scoring our FRAP nuclear segmentation, CellProfiler, Cellpose and StarDist against ground truth

This folder has everything needed to repeat the ground-truth comparison in the
paper. Ten Dorsal FRAP movies that were not used to tune any tool (five imaged
from the ventral side and five from the dorsal side) were segmented by our
FRAP pipeline and by CellProfiler, Cellpose 3 and StarDist. Every method is
then scored against ground-truth nuclei made in Labkit (semi-automatic
pixel classification, checked and corrected by hand) in every frame. It also
measures how consistently each method follows nuclei, and the bleached
nucleus in particular, from frame to frame.

The three tools are run exactly as in our separate CellProfiler, Cellpose and
StarDist packages (same settings, same rule for which channel each one
segments). Those packages compare each tool with our pipeline; this one
compares all four methods with the ground truth.

The raw images and the ground-truth labels are not included because of their
size (see [Data](#data)). You also need to download the Bio-Formats MATLAB
toolbox yourself (see [What you need](#what-you-need)).

## What's in the folder

```
FRAP_GroundTruth_segmentation_comparison/
  README.md                            this file
  run_all.m                            the script you run
  data/gt_embryos.csv                  which movies are used, and how
  data/cellprofiler/*.cppipe           the two CellProfiler pipelines
  environment/setup_cellpose_env.bat   sets up Cellpose 3
  environment/setup_stardist_env.bat   sets up StarDist
  environment/setup_scoring_env.bat    sets up the Python packages used for scoring
  code/frap_pipeline/                  our segmentation, tracking and reconstruction code, as used in the paper
  code/groundtruth/                    MATLAB code for exporting, running the pipeline and the tools
  code/cellpose/run_cellpose.py        runs Cellpose on the exported images
  code/stardist/run_stardist.py        runs StarDist on the exported images
  code/scoring/                        Python code for scoring, tables, figures, overlays and tracking
```

`data/gt_embryos.csv` has one row per movie:

- `file`: path of the `.czi` file inside the data folder.
- `gt_folder`: folder of that movie's label images inside the ground-truth
  folder (the same path as `file`).
- `use_for_comparison`: rows marked `Yes` are analysed.
- `case`: `Good` or `Bad` image quality; `side`: `Ventral` or `Dorsal`. Used
  when pooling the results.
- `filenum`: the movie's number in our full FRAP data list. The pipeline uses
  it to decide how to prepare the image it segments (lines 38–66 of
  `find_nucleizoom_new_comp_set2.m`).
- `bleached_label`: the label our tracking gives the photobleached nucleus.
- `reconstruct`: `yes` for the movies in which our workflow reconstructs the
  bleached nucleus with `reconstruct_missing_nucleus` (filenums 2, 14 and 20).
- `segmented_channel`: which channel is segmented, `mRFP` or
  `inverted Dorsal` (see [Which channel each tool segments](#which-channel-each-tool-segments)).

## What you need

We ran this on Windows 11 with:

- **MATLAB R2024b**, with the Image Processing, Computer Vision and Statistics
  and Machine Learning toolboxes. Any release from R2020a on should work.
- **Bio-Formats MATLAB toolbox (bfmatlab) 5.5.1**, which reads the `.czi`
  files. Download `bfmatlab.zip` for 5.5.1 from
  https://downloads.openmicroscopy.org/bio-formats/5.5.1/ and unzip it. The
  folder should contain `bfopen.m` and `bioformats_package.jar`. Bio-Formats
  is distributed by the Open Microscopy Environment under the GPL v2.
- **CellProfiler 4.2.8** (https://cellprofiler.org/releases), installed in
  the default location or wherever `CELLPROFILER_EXE` in `run_all.m` points.
- **Anaconda or Miniconda**, for three Python environments. Double-click each
  `.bat` file in `environment` once; each writes its progress to a log file
  next to it.
  - `cellpose3`: Python 3.11, Cellpose 3.1.1.3, PyTorch 2.14.1 for CUDA 12.6
    (an NVIDIA GPU is used if present; without one, see the note in the
    `.bat` file). About 3 GB.
  - `stardist`: Python 3.11, TensorFlow 2.15.1, StarDist 0.9.2, CSBDeep 0.8.2,
    NumPy 1.26.4. TensorFlow is kept at 2.15 because newer versions moved to
    Keras 3, which StarDist doesn't support yet. About 1 GB.
  - `frapscore`: Python 3.13 with NumPy, SciPy, pandas, Matplotlib,
    tifffile, h5py, openpyxl and Pillow, used for scoring.
- Optionally, a LaTeX installation, if you want the results tables as PDFs.

A full run takes about one to two hours, most of it our own pipeline and the
figures, and needs about 6 GB of free disk space. The pipeline code writes
some folders using Windows-style paths, so we have only run it on Windows.

## Data

The ten movies and their ground-truth labels (`FRAP_ground_truth_masks.zip`) are
available in the Texas Data Repository, https://doi.org/10.18738/T8/M3TGCY. Each movie has two channels, Dorsal-GFP and an mRFP nuclear
marker, and 60 frames of 512 × 512 pixels (about 0.062 µm per pixel, read
from each file), about 30 s apart. Put the movies in one folder and the labels
in another, keeping the date subfolders:

```
Data/
  2022-10-11/Embryo2_dorsal_4.5x_2sec-1.czi
  2022-11-02/embryo1_spot1.czi
  2022-11-02/embryo2_spot1.czi
  2022-11-21/embryo1_spot 1.czi
  2022-11-24/embryo1_spot1.czi
  2022-11-24/embryo1_spot2.czi
  2022-12-02/embryo1.czi
  2022-12-31/Embryo1_spot1.czi
  2023-01-03/embryo5_spot1.czi
  2023-01-06/embryo1_spot1.czi

GroundTruth/
  2022-10-11/Embryo2_dorsal_4.5x_2sec-1.czi/Ground truth_..._0001.tif ... _0039.tif
  2022-11-02/embryo1_spot1.czi/Ground truth_..._0001.tif ... _0060.tif
  ...
```

Each label image is a binary Labkit export of one frame; the number at the end
of the file name is the movie frame. All 60 frames are labelled, except in
the 2022-10-11 movie (frames 1–39 only) and the 2022-12-31 movie (frames
1–14 only).

## Running it

1. Set up the three Python environments as described above, and install
   CellProfiler. You only need to do this once.
2. Unzip this folder somewhere with a short path, for example
   `C:\FRAP_GT\`. Windows limits file paths to 260 characters, and the
   output folders nest a few levels deep.
3. Open `run_all.m` and set `DATA_DIR`, `GT_DIR` and `BIOFORMATS_DIR`.
   `CELLPROFILER_EXE` and `CONDA_EXE` already point to the standard install
   locations; change them if yours are somewhere else.
4. Run it from MATLAB:
   ```matlab
   run('C:\FRAP_GT\FRAP_GroundTruth_segmentation_comparison\run_all.m')
   ```
   Progress is shown in the command window and saved to `output\run_log.txt`.

If the run stops partway, just run it again; it skips anything already
finished. Set `REDO = true` to redo everything.

## What happens during a run

Each movie is given a short ID made from its `filenum` and date, for example
`E02_2022-10-11`. For every movie:

1. It reads the `.czi` file and saves each channel of each frame as a 16-bit
   TIFF (`output\1_tiffs`), with the pixel values unchanged.
2. It runs our segmentation and tracking exactly as in the paper and saves the
   tracked nuclear masks (`output\2_pipeline_masks`). As in our other
   comparisons, we use the masks before the final 3-pixel erosion, which the
   pipeline applies only when it measures intensities. For the movies marked
   `reconstruct = yes`, it then reconstructs the bleached nucleus with
   `reconstruct_missing_nucleus`, which rebuilds it in frames 2 to N from its
   frame-1 shape, moved with its neighbours. N is chosen from the pipeline's
   own masks, never from the ground truth: it is the last frame in which the
   pipeline, without reconstruction, loses the bleached nucleus or its mask
   falls below 60% of its frame-1 area. The frames used are listed in
   `reconstruction_frames.csv`.
3. It runs CellProfiler (`output\3_cellprofiler`).
4. It runs Cellpose 3 (`output\3_cellpose`).
5. It runs StarDist (`output\3_stardist`).
6. MATLAB then calls `code\scoring\run_scoring.py`, which scores every method
   against the ground truth and makes the tables and figures
   (`output\4_scoring_inputs` to `output\7_failure_examples`).

`reconstruct_missing_nucleus_track.m` is our `reconstruct_missing_nucleus`
adapted for this comparison: it also returns the tracked labels before the
final erosion, which is what we compare.

### Which channel each tool segments

Every tool segments the same channel our pipeline segments, but none of them
gets any of our pipeline's own preprocessing.

- **Six of the movies:** our pipeline segments the mRFP nuclear marker
  directly, so each tool gets the raw mRFP image.
- **Four movies (filenums 16, 17, 19 and 21):** these are dorsal-side movies,
  and for them our pipeline segments the Dorsal-GFP channel. Dorsal-GFP is
  mostly outside the nucleus on the dorsal side, so the raw channel is
  inverted (65535 minus the pixel value) to make nuclei bright. CellProfiler
  then smooths it (Gaussian, 7-pixel artifact diameter) and Cellpose smooths
  it (Gaussian, sigma 3 pixels); StarDist gets it unsmoothed.

Settings, the same as in our separate packages:

- **CellProfiler:** mRFP with adaptive Otsu (window 151, correction 0.8);
  inverted Dorsal with global Otsu after smoothing; nucleus diameter 40–110 px.
- **Cellpose 3:** pretrained `nuclei` model, diameter fixed at 75 px (its own
  size estimate fails on these images).
- **StarDist:** pretrained `2D_versatile_fluo` model with its own thresholds,
  1st–99.8th percentile normalisation, images shrunk to 0.4× before
  prediction so that nuclei are near the size the model was trained on; the
  labels are returned at full size.

### How the masks are scored

The ground-truth labels are the reference. Holes are filled, specks smaller
than 200 pixels are removed (Labkit's pixel classifier leaves some), and each
connected region is one nucleus. The image polarity is detected automatically,
and each label file is checked against its movie frame; a warning is printed
if they do not match.

In each frame, nuclei are paired one to one by intersection over union (IoU),
taking the best-overlapping pairs first. A pair counts as a match if its IoU
is at least 0.5. Nuclei touching the image edge are left out of the counts.
We report precision, recall and F1; the IoU and Dice of matched nuclei; the
area difference and centroid shift; and average precision,
AP = TP / (TP + FP + FN), at IoU thresholds from 0.5 to 0.95. The bleached
nucleus in the ground truth is the labelled nucleus whose mRFP intensity drops
most from frame 1 to 2, followed through later frames by overlap. For each
method we count the post-bleach frames in which it is matched.

**Tracking.** Ground-truth nuclei are linked between consecutive frames (one to
one, IoU ≥ 0.5). For each link we check whether a method finds the nucleus in
both frames and keeps its identity. For our pipeline this uses its own tracked
labels. The other tools do not track, so their objects are linked by the same
overlap rule, which is how they would usually be used. The bleached nucleus is
followed from frame 1 in the same way.

**FRAP recovery (supplementary).** The bleached nucleus is also measured
with each method's masks, the way our pipeline measures it, and the recovery
curve and fitted import/export rates are compared with those from the
ground-truth masks.

## Output

Everything goes into `output`:

- `run_log.txt`, `cellpose_input.csv`, `stardist_input.csv`,
  `scoring_embryos.csv`: the log, and the lists handed to each step.
- `1_tiffs`: the exported images.
- `2_pipeline_masks`: our masks for each movie, with and without
  reconstruction, as `.mat` files; `reconstruction_frames.csv`; the saved
  pipeline state used for reconstruction; and the pipeline's own diagnostic
  images.
- `3_cellprofiler`, `3_cellpose`, `3_stardist`: each tool's label image for
  every frame. The Cellpose and StarDist versions, settings and timings are
  recorded in a `*_run_info.csv` file in their folders.
- `4_scoring_inputs`: the labels packed for scoring.
- `5_results`: the main results.
  - `GT_results.xlsx`, with the sheets Summary (Good, Bad, Ventral, Dorsal
    and All), Per_embryo, Embryos (checks), Per_frame, Bleached_nucleus,
    AP_vs_IoU_threshold and Matched_nuclei; and `GT_summary.csv` and
    `GT_per_embryo.csv`.
  - The results tables as LaTeX for `\input{}`, and a standalone version
    (run `pdflatex` on it for a PDF).
  - `GT_IoU_histograms` and `GT_AP_and_bleached` as PDF, PNG and SVG.
  - `overlays_paper`: one figure per movie, plus every panel on its own at
    2048 × 2048 pixels with a 10 µm scale bar.
- `6_tracking_and_frap`: `FRAP_and_tracking_results.xlsx`, CSV tables, and
  the figures `Tracking_consistency`, `FRAP_curves_per_embryo` and
  `FRAP_accuracy`.
- `7_failure_examples`: overlays of where each tool did badly, judged
  against the ground truth, with our pipeline on the same region for
  comparison, and an index spreadsheet.

In the overlays the ground truth is yellow, drawn just inside each nucleus.
Our pipeline is magenta, CellProfiler green, Cellpose blue and StarDist
orange, each drawn just outside, so all outlines stay visible where they
coincide. The background is the image the method segmented.

With the versions above, a re-run should reproduce the paper's numbers. For
all ten movies pooled:

| | Our pipeline | CellProfiler | Cellpose 3 | StarDist |
|---|---|---|---|---|
| F1 | 0.998 | 0.983 | 0.998 | 0.937 |
| Missed / extra nuclei | 25 / 15 | 176 / 133 | 34 / 10 | 181 / 1042 |
| Median IoU | 0.881 | 0.930 | 0.926 | 0.840 |
| Bleached nucleus found (frames) | 523/523 | 382/523 | 512/523 | 486/523 |
| Bleached nucleus followed, own tracking | 523/523 | 358/523 | 465/523 | 486/523 |
| Nucleus links followed between frames | 99.8% | 97.5% | 99.5% | 97.3% |

Cellpose runs on the GPU and StarDist through TensorFlow; on different
hardware their floating-point results can shift very slightly, so the last
digit may occasionally differ.

## Adding movies

Add a row to `data/gt_embryos.csv`, put the `.czi` file in your data folder and
its label images in your ground-truth folder, and run `run_all.m` again.
Movies that are already done are skipped. Use the movie's real `filenum` from
our data list, because the pipeline's preprocessing depends on it. Set
`segmented_channel` to match the channel our pipeline segments for that
`filenum`. That is `inverted Dorsal` for filenums 5–9, 16, 17, 19, 21, 23 and
27, and `mRFP` otherwise. The failure-example figures in
`7_failure_examples` are chosen for the ten movies of the paper and are
skipped for movies that are not in the list.

## Credits and licence

- CellProfiler: Stirling et al., *BMC Bioinformatics* 22, 433 (2021).
- Cellpose: Stringer et al., *Nature Methods* 18, 100–106 (2021); Stringer &
  Pachitariu, *Nature Methods* (2025).
- StarDist: Schmidt et al., MICCAI 2018; Weigert et al., WACV 2020.
- Labkit: Arzt et al., *Frontiers in Computer Science* 4, 777728 (2022).
- Bio-Formats: Linkert et al., *J Cell Biol* 189, 777–782 (2010).
- Our code: MIT licence (see `LICENSE` in the top folder of the repository),
  except `code/frap_pipeline/fraccircshift.m`, which keeps its own BSD-3-Clause
  licence (MATLAB File Exchange, Bernhard, 2014). If you use the code, please
  cite the paper given in the main README.
