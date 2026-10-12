import logging; logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
import sys, os, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from PIL import Image
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); from h5read import read
import importlib.util
src = open(os.path.join(HERE, 'score_gt_all.py')).read(); fn_src = src[src.index('def gt_instances'):src.index('THR = ')]
MIN_AREA = 200
exec(fn_src)
plt.rcParams.update({'font.family': ['Liberation Sans', 'Arial', 'DejaVu Sans']})
IN, OUT = sys.argv[1], sys.argv[2]; os.makedirs(f'{OUT}/overlays_paper/panels', exist_ok=True)
E = pd.read_csv(f'{IN}/gt_embryos.csv')
cGT = np.array([1.0, 0.85, 0.0]); cM = {'FRAP pipeline': np.array([0.93, 0.10, 0.62]), 'CellProfiler': np.array([0.15, 0.85, 0.25]), 'Cellpose 3': np.array([0.0, 0.75, 1.0]), 'StarDist': np.array([1.0, 0.5, 0.0])}
def band(L, r, inner):
    st = ndi.iterate_structure(ndi.generate_binary_structure(2, 1), r); out = np.zeros(L.shape, bool)
    for sl, lab in zip(ndi.find_objects(L), range(1, L.max() + 1)):
        if sl is None: continue
        sl = tuple(slice(max(s.start - r - 1, 0), s.stop + r + 1) for s in sl); m = L[sl] == lab
        out[sl] |= (m & ~ndi.binary_erosion(m, st)) if inner else (ndi.binary_dilation(m, st) & ~m)
    return out
ONLY = sys.argv[3:]
for _, e in E.iterrows():
    if ONLY and e.id not in ONLY: continue
    eid = e.id; pm = f'{IN}/pipeline_masks/{eid}_pipeline_masks_reconstructed.mat'
    if not os.path.exists(pm): pm = f'{IN}/pipeline_masks/{eid}_pipeline_masks_noreconstruction.mat'
    Lp = read(pm, 'L_track').astype(np.int32); px = float(read(pm, 'scalings').ravel()[0])
    dors = 'Dorsal' in e.segmented_channel
    I = read(pm, 'I_dl' if dors else 'I_nuc').astype(float)
    if dors: I = 65535 - I
    lims = np.percentile(I[:, :, 0], [0.5, 99.7])
    gt = np.load(f'{IN}/stacks/{eid}_GT.npy'); FR = np.load(f'{IN}/stacks/{eid}_GT_frames.npy').tolist()
    rows = [('FRAP pipeline', lambda k: Lp[:, :, k - 1])]
    for m, suf in (('CellProfiler', 'CP'), ('Cellpose 3', 'Cellpose3'), ('StarDist', 'StarDist')):
        f = f'{IN}/stacks/{eid}_{suf}.npy'
        if os.path.exists(f): X = np.load(f); rows.append((m, lambda k, X=X: X[k - 1].astype(int)))
    show = sorted(set([FR[0], FR[1], FR[len(FR) // 2], FR[-1]]))
    fig, ax = plt.subplots(len(rows), len(show), figsize=(3.0 * len(show), 3.15 * len(rows)), squeeze=False)
    for j, k in enumerate(show):
        A = gt_instances(gt[FR.index(k)])[0]; U = 4
        g = np.clip((I[:, :, k - 1] - lims[0]) / (lims[1] - lims[0]), 0, 1); g = np.kron(g, np.ones((U, U)))
        Au = np.kron(A, np.ones((U, U), int))
        for i, (m, B) in enumerate(rows):
            B = B(k); Bu = np.kron(B, np.ones((U, U), int)); rgb = np.dstack([g, g, g])
            rgb[band(Bu, 6, False)] = cM[m]; rgb[band(Au, 6, True)] = cGT
            Lb = int(round(10 / px * U)); hb = 10; mg = 40; rgb[rgb.shape[0] - mg - hb:rgb.shape[0] - mg, rgb.shape[1] - mg - Lb:rgb.shape[1] - mg] = 1.0
            Image.fromarray((rgb * 255).astype(np.uint8)).save(f'{OUT}/overlays_paper/panels/{eid}_frame{k:03d}_{m.replace(" ", "_")}.png', dpi=(600, 600))
            ax[i, j].imshow(rgb); ax[i, j].set_axis_off()
            if i == 0: ax[i, j].set_title(f'frame {k}', fontsize=10)
            if j == 0: ax[i, j].text(-0.04, 0.5, m, transform=ax[i, j].transAxes, rotation=90, va='center', ha='right', color=cM[m], fontsize=10)
    fig.suptitle(f'{eid} ({e.case}, {e.side}; background: {"inverted Dorsal-GFP" if dors else "mRFP"}).  Yellow = ground truth (inner edge), coloured = method (outer edge).  Scale bar 10 µm', fontsize=9)
    fig.tight_layout(); fig.savefig(f'{OUT}/overlays_paper/{eid}_GT_overlays.png', dpi=400); fig.savefig(f'{OUT}/overlays_paper/{eid}_GT_overlays.pdf', dpi=400); plt.close(fig); print(eid, flush=True)
