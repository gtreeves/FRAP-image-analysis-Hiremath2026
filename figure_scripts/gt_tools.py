"""Ground-truth helpers for Fig 13: reading the annotations and masks, and one-to-one matching at IoU >= 0.5.
These functions are the same as in ground_truth/FRAP_GroundTruth_segmentation_comparison/code/scoring/
(score_gt_all.py: gt_instances, touches, candidates, match; make_fail_gt.py: load)."""
import os, numpy as np, pandas as pd
from scipy import ndimage as ndi
from mat_io import read_h5

MIN_AREA = 200
SUF = {'CellProfiler': 'CP', 'Cellpose 3': 'Cellpose3', 'StarDist': 'StarDist'}
_E = {}
cache = {}


def embryos(gt_dir):
    if gt_dir not in _E: _E[gt_dir] = pd.read_csv(f'{gt_dir}/gt_embryos.csv').set_index('id')
    return _E[gt_dir]


def stacks_dir(gt_dir):
    """<gt_dir>/stacks, or <gt_dir>/work/stacks (layout of the working folder)."""
    s = os.path.join(gt_dir, 'stacks')
    return s if os.path.isdir(s) else os.path.join(gt_dir, 'work', 'stacks')


def eid_of(gt_dir, e):
    return [i for i in embryos(gt_dir).index if i.startswith(e + '_')][0]


def gt_instances(a):
    fg = a > 0
    lab0 = ndi.label(fg)[0]
    if lab0.max() > 0 and np.bincount(lab0.ravel())[1:].max() > 0.2 * fg.size:
        fg = ~fg
    fg = ndi.binary_fill_holes(fg)
    L, n = ndi.label(fg); area = np.bincount(L.ravel())
    keep = np.where(area >= MIN_AREA)[0]; keep = keep[keep > 0]
    lut = np.zeros(n + 1, np.int32); lut[keep] = np.arange(1, len(keep) + 1)
    return lut[L], int(n - len(keep))

def touches(L):
    e = np.concatenate([L[0], L[-1], L[:, 0], L[:, -1]]); return set(np.unique(e[e > 0]).tolist())

def candidates(A, B):
    """All overlapping (GT, method) pairs with their IoU, computed once per frame."""
    aA = np.bincount(A.ravel()); aB = np.bincount(B.ravel())
    m = (A > 0) & (B > 0)
    if m.any():
        key = A[m].astype(np.int64) * 100000 + B[m]
        u, cnt = np.unique(key, return_counts=True); pa, pb = u // 100000, u % 100000
    else:
        pa = pb = cnt = np.zeros(0, np.int64)
    iou = cnt / (aA[pa] + aB[pb] - cnt) if len(cnt) else np.zeros(0)
    o = np.argsort(-iou)
    return dict(c=list(zip(iou[o], pa[o], pb[o], cnt[o])), aA=aA, aB=aB, ba=touches(A), bb=touches(B),
                la=set(np.unique(A[A > 0]).tolist()), lb=set(np.unique(B[B > 0]).tolist()))

def match_c(C, thr=0.5):
    ua, ub, M0 = set(), set(), []
    for iou, a, b, it in C['c']:
        if iou < thr: break
        if a in ua or b in ub: continue
        ua.add(a); ub.add(b); M0.append((int(a), int(b), float(iou), int(it), int(C['aA'][a]), int(C['aB'][b])))
    M = [x for x in M0 if x[0] not in C['ba'] and x[1] not in C['bb']]
    fn = sorted(C['la'] - C['ba'] - {x[0] for x in M0})
    fp = sorted(C['lb'] - C['bb'] - {x[1] for x in M0})
    return M, fn, fp

def match(A, B, thr=0.5):
    return match_c(candidates(A, B), thr)

def best_iou(A, lab, B):
    m = A == lab; ov = B[m]; ov = ov[ov > 0]
    if ov.size == 0: return 0.0
    aB = np.bincount(B.ravel()); u, it = np.unique(ov, return_counts=True)
    return float(np.max(it / (m.sum() + aB[u] - it)))


def load(gt_dir, eid):
    if eid in cache: return cache[eid]
    E = embryos(gt_dir)
    pm = f'{gt_dir}/pipeline_masks/{eid}_pipeline_masks_reconstructed.mat'
    if not os.path.exists(pm): pm = f'{gt_dir}/pipeline_masks/{eid}_pipeline_masks_noreconstruction.mat'
    e = E.loc[eid]; dors = 'Dorsal' in e.segmented_channel
    D = dict(Lp=read_h5(pm, 'L_track').astype(np.int32), I=read_h5(pm, 'I_dl' if dors else 'I_nuc').astype(float), Inuc=read_h5(pm, 'I_nuc').astype(float),
             px=float(read_h5(pm, 'scalings').ravel()[0]), t=read_h5(pm, 't').ravel(), dors=dors)
    if dors: D['I'] = 65535 - D['I']
    gt = np.load(f'{stacks_dir(gt_dir)}/{eid}_GT.npy'); FR = np.load(f'{stacks_dir(gt_dir)}/{eid}_GT_frames.npy').tolist()
    G = {k: gt_instances(gt[i])[0] for i, k in enumerate(FR)}
    g1, g2 = G[FR[0]], G[FR[1]]; I1, I2 = D['Inuc'][:, :, FR[0] - 1], D['Inuc'][:, :, FR[1] - 1]; ratios = {}
    for lab in np.unique(g1[g1 > 0]):
        if lab in touches(g1): continue
        m1 = g1 == lab; ov = g2[m1]; ov = ov[ov > 0]
        if ov.size == 0: continue
        l2 = int(np.bincount(ov).argmax()); ratios[(int(lab), l2)] = I2[g2 == l2].mean() / I1[m1].mean()
    (b1, b2), _ = min(ratios.items(), key=lambda x: x[1]); bgt = {FR[0]: b1, FR[1]: b2}
    for k in FR[2:]:
        prev = G[k - 1] == bgt.get(k - 1, -1); cur = G[k][prev]; cur = cur[cur > 0]
        bgt[k] = int(np.bincount(cur).argmax()) if cur.size else 0
    D.update(G=G, bgt=bgt, FR=FR, M={m: np.load(f'{stacks_dir(gt_dir)}/{eid}_{s}.npy') for m, s in SUF.items()})
    cache[eid] = D; return D

