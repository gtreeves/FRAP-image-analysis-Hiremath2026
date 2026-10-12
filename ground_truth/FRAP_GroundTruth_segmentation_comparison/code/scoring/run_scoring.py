"""Score every method against the hand-drawn ground truth (step 6 of run_all.m).

    python run_scoring.py --embryos <output>/scoring_embryos.csv --gt-dir <GroundTruth folder> --out <output>

Reads what the MATLAB steps wrote into <output> (1_tiffs, 2_pipeline_masks, 3_cellprofiler,
3_cellpose, 3_stardist) and the Labkit label images, and writes
  4_scoring_inputs\          label stacks for scoring (.npy) and links to the pipeline masks
  5_results\                 results workbook, tables, histograms, AP figure, overlays
  6_tracking_and_frap\       tracking consistency and FRAP recovery per method
  7_failure_examples\        failure examples of the other tools, judged against the ground truth
Finished parts are not redone unless --redo is given.
"""
import argparse, csv, glob, os, re, shutil, subprocess, sys
import numpy as np
import tifffile

HERE = os.path.dirname(os.path.abspath(__file__))
METHOD_DIRS = {'CP': ('3_cellprofiler', '*_cpNuclei*'), 'Cellpose3': ('3_cellpose', '*_cellpose.tif'),
               'StarDist': ('3_stardist', '*_stardist.tif')}

def frame_no(f, pat):
    m = re.search(pat, os.path.basename(f)); return int(m.group(1)) if m else -1

def prepare(embryos, gt_dir, out, redo):
    IN = os.path.join(out, '4_scoring_inputs'); st = os.path.join(IN, 'stacks'); pm = os.path.join(IN, 'pipeline_masks')
    os.makedirs(st, exist_ok=True); os.makedirs(pm, exist_ok=True)
    rows = list(csv.DictReader(open(embryos, newline='', encoding='utf-8-sig')))
    keep = []
    for r in rows:
        eid = r['id']
        gfs = sorted(glob.glob(os.path.join(gt_dir, r['gt_folder'], '*.tif')), key=lambda f: frame_no(f, r'_(\d+)\.tif$'))
        if not gfs: print(f'{eid}: no ground-truth label images in {os.path.join(gt_dir, r["gt_folder"])}, skipped'); continue
        fr = [frame_no(f, r'_(\d+)\.tif$') for f in gfs]
        if redo or not os.path.exists(os.path.join(st, f'{eid}_GT.npy')):
            gt = np.stack([tifffile.imread(f).astype(np.uint8) for f in gfs])
            np.save(os.path.join(st, f'{eid}_GT.npy'), gt); np.save(os.path.join(st, f'{eid}_GT_frames.npy'), np.array(fr))
            # check that label file k belongs to movie frame k (correlation with the exported movie frames)
            tif = os.path.join(out, '1_tiffs', eid); T = len(glob.glob(os.path.join(tif, '*_C2_T*.tif')))
            if T:
                for i in [0, len(fr) // 2, len(fr) - 1]:
                    g = gt[i].astype(np.float32).ravel(); c = []
                    for t in range(1, T + 1):
                        c.append(max(abs(np.corrcoef(g, tifffile.imread(os.path.join(tif, f'{eid}_{ch}_T{t:03d}.tif')).astype(np.float32).ravel())[0, 1])
                                     for ch in ('C1', 'C2')))
                    best = int(np.argmax(c)) + 1
                    if best != fr[i]: print(f'WARNING {eid}: label file {fr[i]} matches movie frame {best} best - check the label export')
        for suf, (sub, pat) in METHOD_DIRS.items():
            npy = os.path.join(st, f'{eid}_{suf}.npy')
            if os.path.exists(npy) and not redo: continue
            fs = sorted(glob.glob(os.path.join(out, sub, eid, pat)), key=lambda f: frame_no(f, r'_T(\d+)_'))
            if not fs: print(f'{eid}: no {sub} label images, this method is skipped for it'); continue
            np.save(npy, np.stack([tifffile.imread(f).astype(np.uint16) for f in fs]))
        for kind in ('noreconstruction', 'reconstructed'):
            src = os.path.join(out, '2_pipeline_masks', f'{eid}_pipeline_masks_{kind}.mat'); dst = os.path.join(pm, os.path.basename(src))
            if os.path.exists(src) and not os.path.exists(dst):
                try: os.link(src, dst)          # hard link: no extra disk space
                except OSError: shutil.copy2(src, dst)
        r['gt_frames'] = f'{fr[0]}-{fr[-1]}'; keep.append(r)
        print(f'{eid}: {len(fr)} labelled frames ({fr[0]}-{fr[-1]})', flush=True)
    cols = ['id', 'file', 'filenum', 'case', 'side', 'bleached_label', 'reconstruct', 'gt_frames', 'gt_folder', 'segmented_channel']
    with open(os.path.join(IN, 'gt_embryos.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore'); w.writeheader(); w.writerows(keep)
    return IN, [r['id'] for r in keep]

def run(script, *args):
    cmd = [sys.executable, os.path.join(HERE, script), *args]
    print('\n>>', script, *args, flush=True)
    subprocess.run(cmd, check=True)

ap = argparse.ArgumentParser()
ap.add_argument('--embryos', required=True); ap.add_argument('--gt-dir', required=True); ap.add_argument('--out', required=True)
ap.add_argument('--redo', action='store_true')
a = ap.parse_args()
IN, ids = prepare(a.embryos, a.gt_dir, a.out, a.redo)
R = os.path.join(a.out, '5_results'); F = os.path.join(a.out, '6_tracking_and_frap'); X = os.path.join(a.out, '7_failure_examples')
if a.redo or not os.path.exists(os.path.join(R, 'GT_summary.csv')):
    run('score_gt_all.py', IN, R, 'FRAP pipeline,CellProfiler,Cellpose 3,StarDist')
run('table_gt_all.py', R)
run('figs_gt_all.py', R)
for eid in ids:                                   # one embryo per process keeps memory use low
    if a.redo or not os.path.exists(os.path.join(R, 'overlays_paper', f'{eid}_GT_overlays.png')):
        run('overlays_gt_paper.py', IN, R, eid)
if a.redo or not os.path.exists(os.path.join(F, 'frap_fits.csv')):
    run('frap_track.py', IN, F)
run('frap_figs.py', F)
if a.redo or not os.path.exists(os.path.join(X, 'Failure_examples_GT_index.csv')):
    run('make_fail_gt.py', IN, X)
print('\nSCORING DONE')
