# Comparing our FRAP nuclear segmentation with Cellpose

This folder has everything needed to repeat the comparison in the paper
between the nuclear segmentation of our FRAP pipeline and Cellpose. It uses
the same ten Dorsal FRAP movies (five "good" and five "bad" image-quality
cases), the same embryo list and the same statistics as our CellProfiler
comparison, which is a separate package.

We used the pretrained "nuclei" model from Cellpose 3, with the nucleus
diameter set to 75 pixels. We chose it because it was trained specifically on
nuclear images, and like the other tools in our comparison it is told the
approximate size of the nuclei.

We fixed the diameter because Cellpose 3's automatic size estimate failed on
these images. It guessed 9–23 pixels for nuclei that are really about 70–84
pixels across, and the resulting segmentations were unusable. 75 pixels is the
median nucleus diameter we measured on these movies.

The raw images are not included because of their size (see [Data](#data)).
You also need to download the Bio-Formats MATLAB toolbox yourself (see
[What you need](#what-you-need)).

## What's in the folder

```
FRAP_Cellpose_segmentation_comparison/
  README.md                             this file
  run_all.m                             the script you run
  data/embryo_list.csv                  which embryos are used, and how
  environment/setup_cellpose_envs.bat   sets up the Cellpose installation
  code/frap_pipeline/                   our segmentation and tracking code, as used in the paper
  code/comparison/                      code for exporting, matching, statistics and figures
  code/cellpose/run_cellpose.py         runs Cellpose on the exported images
```

`data/embryo_list.csv` is the same file used in the CellProfiler package. It
has one row per movie:

- `file`: path of the `.czi` file inside the data folder.
- `use_for_comparison`: rows marked `Yes` are analysed.
- `case`: `Good` or `Bad`, the image-quality group used when pooling the results.
- `filenum`: the movie's number in our full FRAP data list. The pipeline uses
  it to decide how to prepare the image it segments (lines 38–66 of
  `find_nucleizoom_new_comp_set2.m`).
- `bleached_label`: the label our tracking gives the photobleached nucleus.
- `cellprofiler_input`: which channel is segmented, `mRFP` or
  `inverted Dorsal`. The column name comes from the CellProfiler package, and
  Cellpose uses it in the same way (see [Which channel Cellpose
  segments](#which-channel-cellpose-segments)).

## What you need

We ran this on Windows 11 with:

- **MATLAB R2024b**, with the Image Processing, Computer Vision and Statistics
  and Machine Learning toolboxes. Any release from R2020a on should work.
- **Bio-Formats MATLAB toolbox (bfmatlab) 5.5.1**, which reads the `.czi`
  files. Download `bfmatlab.zip` for 5.5.1 from
  https://downloads.openmicroscopy.org/bio-formats/5.5.1/ and unzip it. The
  folder should contain `bfopen.m` and `bioformats_package.jar`. Bio-Formats
  is distributed by the Open Microscopy Environment under the GPL v2.
- **Anaconda or Miniconda**, for the Cellpose installation described below.
- **An NVIDIA GPU (optional).** We used an RTX 4060 laptop GPU with 8 GB of
  memory. Cellpose also runs on the CPU, just more slowly.
- Optionally, a LaTeX installation, if you want the results tables as PDFs.

Cellpose runs in its own conda environment called `cellpose3`. It contains
Cellpose 3.1.1.3, Python 3.11 and PyTorch 2.14.1 built for CUDA 12.6. We pin
Cellpose 3 because the "nuclei" model was dropped from Cellpose 4.

Double-click `environment\setup_cellpose_envs.bat` to create it. It
downloads a few GB and writes its progress to `environment\setup_log.txt`.
The last line of that file shows the installed versions and whether the GPU
was found. Without an NVIDIA GPU, edit the batch file as described at its
top. The model weights are downloaded automatically the first time Cellpose
runs.

The pipeline code writes some folders using Windows-style paths, so we have
only run it on Windows. A full run takes roughly 45 minutes with a GPU and
needs about 3 GB of free disk space.

## Data

The ten movies are available in the Texas Data Repository, https://doi.org/10.18738/T8/M3TGCY. Each one
has two channels: Dorsal-GFP and an mRFP nuclear marker. Each has 60 frames
of 512 × 512 pixels (roughly 0.06–0.07 µm per pixel, read from each file).
Put them in one folder, keeping the date subfolders:

```
Data/
  2022-10-11/Embryo1_dorsal_5x_1.8sec.czi
  2022-10-11/Embryo2_dorsal_4.5x_2sec.czi
  2022-10-11/Embryo2_dorsal_5x_2sec.czi
  2022-10-12/embryo1_dorsal_3.37sec bleached.czi
  2022-10-12/embryo1_spot 2 dorsal_3.47sec bleached.czi
  2022-10-19/embryo1_spot2_dorsal.czi
  2022-10-25/embryo1_spot1_5x_i thimk its lateral side..czi
  2022-11-02/embryo1_spot2.czi
  2022-11-21/embryo1_spot 2.czi
  2022-11-29/embryo1_spot1.czi
```

## Running it

1. Set up the Cellpose environment as described above. You only need to do
   this once.
2. Unzip this folder somewhere with a short path, for example
   `C:\FRAP_Cellpose\`. Windows limits file paths to 260 characters, and the
   output folders nest a few levels deep.
3. Open `run_all.m` and set `DATA_DIR` and `BIOFORMATS_DIR`. `CONDA_EXE`
   already points to a standard Anaconda install; change it if yours is
   somewhere else.
4. Run it from MATLAB:
   ```matlab
   run('C:\FRAP_Cellpose\FRAP_Cellpose_segmentation_comparison\run_all.m')
   ```
   Progress is shown in the command window and saved to `output\run_log.txt`.

If the run stops partway, just run it again; it skips anything already
finished. Set `REDO = true` to redo everything.

## What happens during a run

Each movie is given a short ID made from its `filenum` and date, for example
`E03_2022-10-11`. The run then does the following for every movie:

1. It reads the `.czi` file and saves each channel of each frame as a 16-bit
   TIFF (`output\1_tiffs`), with the pixel values unchanged.
2. It runs our segmentation and tracking exactly as in the paper and saves the
   tracked nuclear masks (`output\2_pipeline_masks`). As in the CellProfiler
   comparison, we use the masks before the final 3-pixel erosion, which the
   pipeline applies only when it measures intensities. The log also confirms
   that the bleached-nucleus label in `embryo_list.csv` belongs to the nucleus
   that dims most between frames 1 and 2.
3. MATLAB calls `run_cellpose.py`, which segments the
   channel named in `cellprofiler_input` (see below) and saves one label
   image per frame (`output\3_cellpose`). The Cellpose version, settings and
   timings are recorded in `cellpose_run_info.csv` in each model's folder.
4. It matches nuclei between our pipeline and Cellpose and computes the
   statistics (`output\4_comparison`).
5. It draws outline overlays for figures (`output\5_overlays`).

Finally it pools all movies into Good, Bad and All groups and makes the
histograms and tables (`output\6_paper_outputs`), separately for each
Cellpose model.

### Which channel Cellpose segments

Cellpose always segments the same channel our pipeline segments, but it never
gets any of our pipeline's own preprocessing.

- **Seven of the movies:** our pipeline segments the mRFP nuclear marker
  directly, so Cellpose gets the raw mRFP image.
- **Three movies (filenums 5, 6 and 9):** these are dorsal-side movies, and
  for them our pipeline segments the Dorsal-GFP channel instead, after
  scaling and inverting it. Dorsal-GFP is mostly outside the nucleus on the
  dorsal side, so once inverted, the nuclei look bright. For these movies
  `run_cellpose.py` takes the raw Dorsal-GFP channel, inverts it
  (65535 minus the pixel value) and smooths it with a Gaussian filter
  (sigma 3 pixels). That is the same preparation the CellProfiler package
  uses (ImageMath "Invert", then Smooth with a 7-pixel artifact diameter), so
  Cellpose and CellProfiler get the same image. None of our scaling or
  filtering is used.

We tested this on movie 5 before running the rest. Unlike CellProfiler,
Cellpose finds the nuclei in the inverted Dorsal channel even without
smoothing (F1 1.000). Smoothing mainly tightens the outlines: the median IoU
went from 0.910 to 0.937.

Cellpose normalises every image itself (1st–99th percentile), so we don't
rescale anything.

### How the masks are compared

Our pipeline's masks are the reference. In each frame, nuclei are paired one
to one by intersection over union (IoU), taking the best-overlapping pairs
first. A pair counts as a match if its IoU is at least 0.5. Nuclei touching
the image edge are left out of the counts.

From the matches we report:

- **Detection:** precision, recall and F1, pooled over all frames.
- **Outline agreement:** IoU and Dice of matched nuclei, and the pixel-level
  Dice of each frame.
- **Size and position:** the difference in nuclear area (Cellpose minus
  pipeline, in µm²) and the centroid shift.

For the bleached nucleus, we also record its best IoU with any Cellpose
object in every post-bleach frame, and how often that reaches 0.5.

## Output

Everything goes into `output`:

- `run_log.txt` and `cellpose_input.csv`: the log, and the image folders given
  to Cellpose.
- `1_tiffs`: the exported images.
- `2_pipeline_masks`: our masks for each movie as a `.mat` file and as a
  multi-page TIFF, plus the pipeline's own diagnostic images.
- `3_cellpose\Cellpose3_nuclei`: Cellpose's label image for each frame, and
  `cellpose_run_info.csv`.
- `4_comparison\Cellpose3_nuclei`: per-movie results (summary, per-frame and
  per-nucleus tables, the bleached nucleus, and quick-look figures).
- `5_overlays\Cellpose3_nuclei`: figure-ready overlays, described below.
- `6_paper_outputs\Cellpose3_nuclei`: the pooled results. These are the Good, Bad and
  All IoU histograms (PDF, PNG and MATLAB figure) and the results table
  (Excel, a LaTeX table for `\input{}`, and a standalone LaTeX version; run
  `pdflatex` on the standalone one for a PDF).

Each overlay frame is a 2048 × 2048 PNG of selected frames, with a 10 µm
scale bar and no text. There is also a four-panel PDF/PNG per movie. The
background is the image Cellpose segmented: mRFP, or the inverted Dorsal-GFP
channel for the three dorsal movies. Our outline is magenta, drawn just
inside each nucleus. Cellpose's is green, drawn just outside, so both stay
visible where they coincide.

With the versions above, a re-run should reproduce the paper's numbers. For
all cases pooled:

| | Good | Bad | All |
|---|---|---|---|
| Matched nuclei | 4990 | 5471 | 10461 |
| F1 | 0.986 | 0.994 | 0.990 |
| Median IoU | 0.880 | 0.901 | 0.890 |
| Bleached nucleus matched (frames) | 288/295 | 295/295 | 583/590 |

On a different GPU or driver, Cellpose's floating-point results can shift
very slightly, so the last digit may occasionally differ.

## Adding movies

Add a row to `data/embryo_list.csv`, put the `.czi` file in your data folder,
and run `run_all.m` again. Movies that are already done are skipped. Use the
movie's real `filenum` from our data list, because the pipeline's
preprocessing depends on it. Set `cellprofiler_input` to match the channel
our pipeline segments for that `filenum`. That is `inverted Dorsal` for
filenums 5–9, 16, 17, 19, 21, 23 and 27, and `mRFP` otherwise.

## Credits and licence

- Cellpose 3: Stringer and Pachitariu, *Nature Methods* 22, 592–599 (2025).
- Bio-Formats: Linkert et al., *J Cell Biol* 189, 777–782 (2010).
- Our code: MIT licence (see `LICENSE` in the top folder of the repository).
  If you use it, please cite the paper given in the main README.
