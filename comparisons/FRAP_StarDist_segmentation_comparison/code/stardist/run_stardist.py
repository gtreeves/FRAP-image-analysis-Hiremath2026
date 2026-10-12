"""Segment the nuclei of every embryo with StarDist.

    python run_stardist.py --embryos embryos.csv --scale 0.4 --out <dir>

embryos.csv columns: short_id, tif_dir, input
  tif_dir holds the exported frames <...>_C1_T###.tif (Dorsal-GFP) and
  <...>_C2_T###.tif (mRFP). input says which image StarDist gets:
    mrfp               the mRFP nuclear marker (C2), raw (default if no column)
    dorsal_inv         the Dorsal-GFP channel (C1) inverted: 65535 - I.
                       Used for the dorsal embryos. No smoothing: in our
                       trial on E05, smoothing made StarDist call patches of
                       cytoplasm between nuclei as extra nuclei.
    dorsal_inv_smooth  inverted and Gaussian smoothed (sigma 3 px), as given
                       to CellProfiler and Cellpose (trial only)
Model: pretrained 2D_versatile_fluo with its own default thresholds.
Each image is normalised to its 1st-99.8th percentile (the StarDist default
recipe). --scale rescales the image before prediction; the labels are
returned at the original size. 1 = no rescaling.
Writes <out>/<short_id>/<short_id>_T###_stardist.tif (uint16 labels)
and <out>/stardist_run_info.csv.
"""
import argparse, csv, glob, os, re, time, platform
from importlib.metadata import version as pkg_version
import numpy as np
from scipy.ndimage import gaussian_filter
import tifffile
from csbdeep.utils import normalize
from stardist.models import StarDist2D

ap = argparse.ArgumentParser()
ap.add_argument('--embryos', required=True)
ap.add_argument('--scale', type=float, default=1.0)
ap.add_argument('--out', required=True)
ap.add_argument('--redo', action='store_true')
a = ap.parse_args()

model = StarDist2D.from_pretrained('2D_versatile_fluo')
thr = model.thresholds
print(f"StarDist {pkg_version('stardist')}, TensorFlow {pkg_version('tensorflow')}, "
      f"model 2D_versatile_fluo, prob {thr.prob:.3f}, nms {thr.nms:.2f}, scale {a.scale}", flush=True)

os.makedirs(a.out, exist_ok=True)
info = []
old_info = {}
info_csv = os.path.join(a.out, 'stardist_run_info.csv')
if os.path.exists(info_csv):
    with open(info_csv, newline='') as f:
        for r0 in csv.DictReader(f):
            if not r0.get('input'): r0['input'] = 'mrfp'
            old_info[r0['short_id']] = r0
with open(a.embryos, newline='') as f:
    rows = list(csv.DictReader(f))
for r in rows:
    sid, tdir = r['short_id'], r['tif_dir']
    inp = (r.get('input') or 'mrfp').strip()
    assert inp in ('mrfp', 'dorsal_inv', 'dorsal_inv_smooth'), f'unknown input {inp}'
    ch = 'C2' if inp == 'mrfp' else 'C1'
    files = sorted(glob.glob(os.path.join(tdir, f'*_{ch}_T*.tif')))
    assert files, f'no {ch} frames in {tdir}'
    od = os.path.join(a.out, sid); os.makedirs(od, exist_ok=True)
    t0 = time.time(); n_obj = []
    for fn in files:
        t = int(re.search(r'_T(\d+)\.tif$', fn).group(1))
        outf = os.path.join(od, f'{sid}_T{t:03d}_stardist.tif')
        if os.path.exists(outf) and not a.redo:
            continue
        img = tifffile.imread(fn).astype(np.float32)
        if inp != 'mrfp':
            img = 65535.0 - img
        if inp == 'dorsal_inv_smooth':
            img = gaussian_filter(img, 3.0)
        img = normalize(img, 1, 99.8, axis=(0, 1))
        labels, _ = model.predict_instances(img, scale=None if a.scale == 1 else a.scale)
        tifffile.imwrite(outf, labels.astype(np.uint16))
        n_obj.append(int(labels.max()))
    dt = time.time() - t0
    if not n_obj and sid in old_info:          # nothing new for this embryo: keep its earlier record
        info.append(old_info[sid]); print(f'{sid}: already done', flush=True); continue
    print(f'{sid}: {len(files)} frames, {dt:.0f} s, median objects {np.median(n_obj) if n_obj else "-"}', flush=True)
    info.append(dict(short_id=sid, input=inp, frames=len(files), seconds=round(dt, 1),
                     median_objects=float(np.median(n_obj)) if n_obj else '',
                     model='2D_versatile_fluo', scale=a.scale, prob_thresh=round(float(thr.prob), 4),
                     nms_thresh=round(float(thr.nms), 3), normalization='percentile 1-99.8',
                     stardist=pkg_version('stardist'), csbdeep=pkg_version('csbdeep'),
                     tensorflow=pkg_version('tensorflow'), python=platform.python_version()))
fields = list(dict.fromkeys(k for d in info for k in d))
with open(info_csv, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fields, restval=''); w.writeheader(); w.writerows(info)
print('STARDIST DONE', flush=True)
