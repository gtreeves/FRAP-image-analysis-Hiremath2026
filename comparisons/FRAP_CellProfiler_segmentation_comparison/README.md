# Comparing our FRAP nuclear segmentation with CellProfiler

This folder has everything needed to repeat the comparison in the paper
between the nuclear segmentation of our FRAP pipeline and CellProfiler 4.2.8.
We ran both methods on ten Dorsal FRAP movies of *Drosophila* embryos,
grouped into five "good" and five "bad" image-quality cases. We then measured
how well the two sets of nuclear masks agree.

The raw images are not included because of their size (see
[Data](#data)). You also need to download the Bio-Formats MATLAB toolbox
yourself (see [What you need](#what-you-need)).

## What's in the folder

```
FRAP_CellProfiler_segmentation_comparison/
  README.md                 this file
  run_all.m                 the script you run
  data/embryo_list.csv      which embryos are used, and how
  cellprofiler/             the four CellProfiler pipelines (.cppipe)
  code/frap_pipeline/       our segmentation and tracking code, as used in the paper
  code/comparison/          code for exporting, matching, statistics and figures
```

`data/embryo_list.csv` has one row per movie:

- `file`: path of the `.czi` file inside the data folder.
- `use_for_comparison`: rows marked `Yes` are analysed.
- `case`: `Good` or `Bad`, the image-quality group used when pooling the results.
- `filenum`: the movie's number in our full FRAP data list. The pipeline uses
  it to decide how to prepare the image it segments (lines 38–66 of
  `find_nucleizoom_new_comp_set2.m`).
- `bleached_label`: the label our tracking gives the photobleached nucleus.
- `cellprofiler_input`: which channel CellProfiler segments, `mRFP` or
  `inverted Dorsal` (see [Which channel CellProfiler
  segments](#which-channel-cellprofiler-segments)).

## What you need

We ran this on Windows 11 with:

- **MATLAB R2024b**, with the Image Processing, Computer Vision and Statistics
  and Machine Learning toolboxes. Any release from R2020a on should work.
- **CellProfiler 4.2.8**, from https://cellprofiler.org/releases. MATLAB runs
  it in the background, so you never need to open it.
- **Bio-Formats MATLAB toolbox (bfmatlab) 5.5.1**, which reads the `.czi`
  files. Download `bfmatlab.zip` for 5.5.1 from
  https://downloads.openmicroscopy.org/bio-formats/5.5.1/ and unzip it. The
  folder should contain `bfopen.m` and `bioformats_package.jar`. Newer
  versions will probably work, but 5.5.1 is the one we used. Bio-Formats is
  distributed by the Open Microscopy Environment under the GPL v2.
- Optionally, a LaTeX installation, if you want the results table as a PDF.

The pipeline code writes some folders using Windows-style paths, so we have
only run it on Windows. A full run takes about 45 minutes on a desktop
computer and needs about 4 GB of free disk space.

## Data

The ten movies are available in the Texas Data Repository, https://doi.org/10.18738/T8/M3TGCY. Each one
has two channels: Dorsal-GFP and an mRFP nuclear marker. Each has 60 frames
of 512 × 512 pixels (roughly 0.06–0.07 µm per pixel, read from each file). Put them in one folder,
keeping the date subfolders:

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

1. Unzip this folder somewhere with a short path, for example `C:\FRAP_CP\`.
   CellProfiler can't save files whose full path is longer than 260
   characters, and the output folders nest a few levels deep.
2. Open `run_all.m` and set the three paths at the top: `DATA_DIR`,
   `BIOFORMATS_DIR`, and `CELLPROFILER_EXE` (only needed if CellProfiler isn't
   installed in the default place).
3. Run it from MATLAB:
   ```matlab
   run('C:\FRAP_CP\FRAP_CellProfiler_segmentation_comparison\run_all.m')
   ```
   Progress is shown in the command window and saved to `output\run_log.txt`.

If the run stops partway, just run it again; it skips anything already
finished. Set `REDO = true` to redo everything.

To get the table as a PDF, run
`pdflatex results_table_standalone.tex` inside
`output\6_paper_outputs\CellProfiler_main`.

## What happens during a run

Each movie is given a short ID made from its `filenum` and date, for example
`E03_2022-10-11`. The run then does the following for every movie:

1. It reads the `.czi` file and saves each channel of each frame as a 16-bit
   TIFF (`output\1_tiffs`), with the pixel values unchanged.
2. It runs our segmentation and tracking exactly as in the paper and saves the
   tracked nuclear masks (`output\2_pipeline_masks`). We compare the masks as
   they are before the final 3-pixel erosion. The pipeline applies that
   erosion only when it measures intensities, and leaving it out makes the
   comparison fairer to both methods. As a check, the log confirms that the
   bleached-nucleus label in `embryo_list.csv` belongs to the nucleus that
   dims most between frames 1 and 2.
3. It runs CellProfiler on the channel named in `cellprofiler_input`
   (`output\3_cellprofiler`).
4. It matches nuclei between the two methods and computes the statistics
   (`output\4_comparison`).
5. It draws outline overlays for figures (`output\5_overlays`).

Finally it pools all movies into Good, Bad and All groups and makes the
histograms and tables (`output\6_paper_outputs`).

### Which channel CellProfiler segments

CellProfiler always segments the same channel our pipeline segments, but it
never gets any of our pipeline's own preprocessing.

- **Seven of the movies:** our pipeline segments the mRFP nuclear marker
  directly, so CellProfiler does too.
- **Three movies (filenums 5, 6 and 9):** these are dorsal-side movies, and
  for them our pipeline segments the Dorsal-GFP channel instead, after
  scaling and inverting it. Dorsal-GFP is
  mostly outside the nucleus on the dorsal side, so inverted, the nuclei look
  bright. For these movies CellProfiler also starts from the raw Dorsal-GFP
  channel, but does its own preparation. It inverts the image with ImageMath
  ("Invert"), then smooths it with its Smooth module (Gaussian filter,
  artifact diameter 7 pixels). It uses none of our scaling or filtering.

The smoothing step is needed because the Dorsal-GFP channel is much noisier
than mRFP. Without it, CellProfiler merges nuclei with patches of cytoplasm.
We tested this on movie 1, where both channels are usable: with inversion
alone, CellProfiler matched almost none of our nuclei (F1 0.26). With
inversion and smoothing it reached F1 0.97.

### The CellProfiler settings

All four pipelines use the standard modules (Images, Metadata,
NamesAndTypes, Groups). Nuclei are found with IdentifyPrimaryObjects, and each
frame's nuclei are saved as a 16-bit label image. All use a typical nucleus
diameter of 40–110 pixels and Otsu two-class thresholding. Touching nuclei
are split by shape, and holes are filled. Nuclei touching the image border are
kept here and removed later, in the matching step.

We report two settings. Each has one pipeline for mRFP movies and one for
inverted Dorsal movies.

**`CellProfiler_main` (main result):**

- mRFP movies: `FRAP_nuclei_adaptive`, adaptive Otsu with a 151-pixel window
  and a correction factor of 0.8. This keeps the dimmed, bleached nucleus in
  more frames.
- Dorsal movies: `FRAP_nuclei_invertedDorsal_smooth_global`, which inverts,
  smooths, then applies global Otsu with a correction factor of 1.0.

**`CellProfiler_default` (supplementary, closer to CellProfiler's defaults):**

- mRFP movies: `FRAP_nuclei_global_default`, global Otsu with a correction
  factor of 1.0.
- Dorsal movies: `FRAP_nuclei_invertedDorsal_global_default`, which inverts
  and applies global Otsu, without smoothing.

### How the masks are compared

Our pipeline's masks are the reference. In each frame, nuclei are paired one
to one by intersection over union (IoU), taking the best-overlapping pairs
first. A pair counts as a match if its IoU is at least 0.5. Nuclei touching
the image edge are left out of the counts.

From the matches we report:

- **Detection:** precision, recall and F1, pooled over all frames.
- **Outline agreement:** IoU and Dice of matched nuclei, and the pixel-level
  Dice of each frame.
- **Size and position:** the difference in nuclear area (CellProfiler minus
  pipeline, in µm²) and the centroid shift.

For the bleached nucleus, we also record its best IoU with any CellProfiler
object in every post-bleach frame, and how often that reaches 0.5.

## Output

Everything goes into `output`, with a separate subfolder per setting
(`CellProfiler_main`, `CellProfiler_default`) from step 3 on:

- `run_log.txt` and `embryos_used.csv`: the log, and the list of IDs with
  their files and input channel.
- `1_tiffs`: the images given to CellProfiler.
- `2_pipeline_masks`: our masks for each movie as a `.mat` file and as a
  multi-page TIFF, plus the pipeline's own diagnostic images.
- `3_cellprofiler`: CellProfiler's label image for each frame.
- `4_comparison`: per-movie results (summary, per-frame and per-nucleus
  tables, the bleached nucleus, and quick-look figures).
- `5_overlays` (main setting only): figure-ready overlays, described below.
- `6_paper_outputs`: the pooled results. These are the Good, Bad and All IoU
  histograms (PDF, PNG and MATLAB figure) and the results table (Excel, a
  LaTeX table for `\input{}`, and a standalone LaTeX version).

Each overlay frame is a 2048 × 2048 PNG of selected frames, with a 10 µm
scale bar and no text. There is also a four-panel PDF/PNG per movie. The
background is the image CellProfiler segmented: mRFP, or the inverted
Dorsal-GFP channel for the three dorsal movies. Our outline is magenta, drawn
just inside each nucleus. CellProfiler's is green, drawn just outside, so both
stay visible where they coincide.

Every step is deterministic, so a re-run should give the paper's numbers. For
`CellProfiler_main` with all ten movies pooled:

| | Good | Bad | All |
|---|---|---|---|
| Matched nuclei | 4880 | 5439 | 10319 |
| F1 | 0.970 | 0.985 | 0.978 |
| Median IoU | 0.894 | 0.898 | 0.896 |
| Bleached nucleus matched (frames) | 190/295 | 274/295 | 464/590 |

## Adding movies

Add a row to `data/embryo_list.csv`, put the `.czi` file in your data folder,
and run `run_all.m` again. Movies that are already done are skipped. Use the
movie's real `filenum` from our data list, because the pipeline's
preprocessing depends on it. Set `cellprofiler_input` to match the channel
our pipeline segments for that `filenum`. That is `inverted Dorsal` for
filenums 5–9, 16, 17, 19, 21, 23 and 27, and `mRFP` otherwise.

## Credits and licence

- Bio-Formats: Linkert et al., *J Cell Biol* 189, 777–782 (2010).
- CellProfiler: Stirling et al., *BMC Bioinformatics* 22, 433 (2021).
- Our code: MIT licence (see `LICENSE` in the top folder of the repository).
  If you use it, please cite the paper given in the main README.
