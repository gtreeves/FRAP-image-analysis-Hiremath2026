"""Segment the nuclei of every embryo with Cellpose.

Usage (inside the conda env that holds the Cellpose version you want):
    python run_cellpose.py --embryos embryos.csv --model nuclei --diameter 75 --out <dir>
(Cellpose 3, "nuclei" model)

embryos.csv columns: short_id, tif_dir, input
  tif_dir holds the exported frames <...>_C1_T###.tif (Dorsal-GFP) and
  <...>_C2_T###.tif (mRFP). input says which image Cellpose gets:
    mrfp               the mRFP nuclear marker (C2), raw (default if no column)
    dorsal_inv_smooth  the Dorsal-GFP channel (C1) inverted (65535 - I) and
                       Gaussian smoothed with sigma 3 px, the same preparation
                       CellProfiler gets (ImageMath Invert + Smooth, 7 px).
    dorsal_inv         inverted only (used in the trial)
  Nothing else is done to the image; Cellpose normalises it itself.
Writes <out>/<short_id>/<short_id>_T###_cellpose.tif (uint16 labels, 0 = background)
and <out>/cellpose_run_info.csv (versions, settings, timings, estimated diameters).
"""
import argparse, csv, glob, os, re, time, platform
import numpy as np
from scipy.ndimage import gaussian_filter
import tifffile
import torch
from importlib.metadata import version as pkg_version
from cellpose import models
CP_VERSION = pkg_version('cellpose')

ap = argparse.ArgumentParser()
ap.add_argument('--embryos', required=True)
ap.add_argument('--model', default='nuclei', choices=['nuclei'])
ap.add_argument('--out', required=True)
ap.add_argument('--redo', action='store_true')
ap.add_argument('--diameter', type=float, default=None,
                help='nuclei model only: nucleus diameter in px (None = Cellpose size estimate)')
a = ap.parse_args()

gpu = torch.cuda.is_available()
print(f'Cellpose {CP_VERSION}, torch {torch.__version__}, CUDA {gpu}, model {a.model}', flush=True)
major = int(str(CP_VERSION).split('.')[0])

assert major == 3, 'the nuclei model run uses Cellpose 3'
model = models.Cellpose(gpu=gpu, model_type='nuclei')
def segment(img):
    masks, flows, styles, diams = model.eval(img, diameter=a.diameter, channels=[0, 0])
    return masks, float(np.atleast_1d(diams)[0])

os.makedirs(a.out, exist_ok=True)
info = []
old_info = {}
info_csv = os.path.join(a.out, 'cellpose_run_info.csv')
if os.path.exists(info_csv):
    with open(info_csv, newline='') as f:
        for r0 in csv.DictReader(f):
            r0.setdefault('input', 'mrfp')
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
    t0 = time.time(); diams = []; n_obj = []
    for fn in files:
        t = int(re.search(r'_T(\d+)\.tif$', fn).group(1))
        outf = os.path.join(od, f'{sid}_T{t:03d}_cellpose.tif')
        if os.path.exists(outf) and not a.redo:
            continue
        img = tifffile.imread(fn).astype(np.float32)
        if inp != 'mrfp':
            img = 65535.0 - img
        if inp == 'dorsal_inv_smooth':
            img = gaussian_filter(img, 3.0)
        masks, d = segment(img)
        tifffile.imwrite(outf, masks.astype(np.uint16))
        diams.append(d); n_obj.append(int(masks.max()))
    dt = time.time() - t0
    if not n_obj and sid in old_info:          # nothing new for this embryo: keep its earlier record
        info.append(old_info[sid]); print(f'{sid}: already done', flush=True); continue
    print(f'{sid}: {len(files)} frames, {dt:.0f} s, median objects {np.median(n_obj) if n_obj else "-"}', flush=True)
    info.append(dict(short_id=sid, input=inp, frames=len(files), seconds=round(dt, 1),
                     median_objects=float(np.median(n_obj)) if n_obj else '',
                     median_estimated_diameter_px=float(np.nanmedian(diams)) if diams and not all(np.isnan(diams)) else '',
                     model=a.model, diameter_setting=a.diameter if a.diameter else 'auto', cellpose_version=str(CP_VERSION), torch=torch.__version__,
                     cuda=gpu, python=platform.python_version()))
fields = list(dict.fromkeys(k for d in info for k in d))
with open(info_csv, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fields, restval=''); w.writeheader(); w.writerows(info)
print('CELLPOSE DONE', flush=True)
