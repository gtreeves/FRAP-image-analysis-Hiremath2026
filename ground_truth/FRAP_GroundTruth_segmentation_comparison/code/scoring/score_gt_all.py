"""Score segmentation methods against the Labkit ground truth, all GT embryos.

Inputs (staged from the "Ground truth comparison" folder):
  gt_embryos.csv                         embryo list (id, case, side, bleached label, ...)
  stacks/<id>_GT.npy                     Labkit labels, one plane per labelled frame (+ _GT_frames.npy)
  stacks/<id>_CP.npy                     CellProfiler label stack (all frames)
  pipeline_masks/<id>_pipeline_masks_{reconstructed,noreconstruction}.mat

Ground truth: binary Labkit export, file k = movie frame k. Polarity is
detected (nuclei = the value forming many separate blobs). Holes are filled
and specks < 200 px removed. Each nucleus = one connected region.
Matching: one-to-one, highest IoU first, IoU >= 0.5; a matched pair is dropped
if either object touches the image edge; unmatched edge objects are dropped.
GT is the reference. AP(t) = TP/(TP+FP+FN) at IoU threshold t.
Bleached nucleus in the GT: the labelled nucleus whose mRFP intensity drops
most from frame 1 to 2, followed through later frames by overlap.
"""
import sys, os, json, numpy as np, pandas as pd
from scipy import ndimage as ndi
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from h5read import read

IN = sys.argv[1] if len(sys.argv) > 1 else 'in'
OUT = sys.argv[2] if len(sys.argv) > 2 else 'out'
METHODS = sys.argv[3].split(',') if len(sys.argv) > 3 else ['FRAP pipeline', 'CellProfiler']
MIN_AREA = 200
os.makedirs(OUT, exist_ok=True)
E = pd.read_csv(f'{IN}/gt_embryos.csv')

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

THR = np.round(np.arange(0.5, 0.951, 0.05), 2)
rows_f, rows_p, rows_b, rows_ap, rows_e = [], [], [], [], []
for _, e in E.iterrows():
    eid = e.id
    gt = np.load(f'{IN}/stacks/{eid}_GT.npy'); FR = np.load(f'{IN}/stacks/{eid}_GT_frames.npy').tolist()
    pm = f'{IN}/pipeline_masks/{eid}_pipeline_masks_reconstructed.mat'
    recon = os.path.exists(pm)
    if not recon: pm = f'{IN}/pipeline_masks/{eid}_pipeline_masks_noreconstruction.mat'
    px = float(read(pm, 'scalings').ravel()[0]); tvec = read(pm, 't').ravel()
    Lp = read(pm, 'L_track').astype(np.int32); Inuc = read(pm, 'I_nuc')
    src = {'FRAP pipeline': lambda k: Lp[:, :, k - 1]}
    if os.path.exists(f'{IN}/stacks/{eid}_CP.npy'):
        CPs = np.load(f'{IN}/stacks/{eid}_CP.npy'); src['CellProfiler'] = lambda k, C=CPs: C[k - 1].astype(np.int32)
    for extra in ('Cellpose 3', 'StarDist'):
        f = f'{IN}/stacks/{eid}_{extra.replace(" ", "")}.npy'
        if os.path.exists(f): X = np.load(f); src[extra] = lambda k, X=X: X[k - 1].astype(np.int32)
    G = {}; nspk = {}
    for i, k in enumerate(FR): G[k], nspk[k] = gt_instances(gt[i])
    # bleached nucleus in GT
    g1, g2 = G[FR[0]], G[FR[1]]
    I1, I2 = Inuc[:, :, FR[0] - 1].astype(float), Inuc[:, :, FR[1] - 1].astype(float)
    ratios = {}
    for lab in np.unique(g1[g1 > 0]):
        if lab in touches(g1): continue
        m1 = g1 == lab; ov = g2[m1]; ov = ov[ov > 0]
        if ov.size == 0: continue
        l2 = int(np.bincount(ov).argmax()); ratios[(int(lab), l2)] = I2[g2 == l2].mean() / I1[m1].mean()
    (bl1, bl2), rmin = min(ratios.items(), key=lambda x: x[1])
    bgt = {FR[0]: bl1, FR[1]: bl2}
    for k in FR[2:]:
        prev = G[k - 1] == bgt.get(k - 1, -1); cur = G[k][prev]; cur = cur[cur > 0]
        bgt[k] = int(np.bincount(cur).argmax()) if cur.size else 0
    # check: does the GT bleached nucleus coincide with the pipeline's bleached label in frame 1?
    pl = int(e.bleached_label); chk = best_iou(G[FR[0]], bl1, (Lp[:, :, FR[0] - 1] == pl).astype(np.int32))
    rows_e.append(dict(embryo=eid, case=e.case, side=e.side, frames_scored=len(FR), pipeline_reconstructed=recon,
                       bleach_ratio_frame1to2=round(rmin, 2), GT_bleached_vs_pipeline_label_IoU_frame1=round(chk, 3),
                       GT_specks_removed_per_frame=round(np.mean(list(nspk.values())), 1)))
    for meth in METHODS:
        if meth not in src: continue
        for k in FR:
            A = G[k]; B = src[meth](k)
            C = candidates(A, B)
            M, fn, fp = match_c(C)
            TP, FN, FP = len(M), len(fn), len(fp)
            inter = np.sum((A > 0) & (B > 0)); pd_ = 2 * inter / max((A > 0).sum() + (B > 0).sum(), 1)
            rows_f.append(dict(method=meth, embryo=eid, case=e.case, side=e.side, frame=k, t_s=round(float(tvec[k - 1]), 1),
                               n_gt=TP + FN, n_method=TP + FP, TP=TP, FP=FP, FN=FN, precision=TP / max(TP + FP, 1),
                               recall=TP / max(TP + FN, 1), F1=2 * TP / max(2 * TP + FP + FN, 1),
                               mean_IoU=np.mean([x[2] for x in M]) if M else np.nan, pixel_Dice=pd_))
            if M:
                cA = ndi.center_of_mass(np.ones_like(A), A, [x[0] for x in M]); cB = ndi.center_of_mass(np.ones_like(B), B, [x[1] for x in M])
                for x, ca, cb in zip(M, cA, cB):
                    rows_p.append(dict(method=meth, embryo=eid, case=e.case, side=e.side, frame=k, gt_label=x[0], method_label=x[1],
                                       IoU=x[2], Dice=2 * x[3] / (x[4] + x[5]), area_gt_um2=x[4] * px ** 2, area_method_um2=x[5] * px ** 2,
                                       area_diff_um2=(x[5] - x[4]) * px ** 2, centroid_shift_um=float(np.hypot(ca[0] - cb[0], ca[1] - cb[1]) * px),
                                       bleached=(x[0] == bgt[k])))
            for th in THR:
                Mt, fnt, fpt = match_c(C, th)
                rows_ap.append(dict(method=meth, embryo=eid, case=e.case, side=e.side, frame=k, threshold=th, TP=len(Mt), FP=len(fpt), FN=len(fnt)))
            if k != FR[0] and bgt.get(k):
                bi = best_iou(A, bgt[k], B)
                rows_b.append(dict(method=meth, embryo=eid, case=e.case, side=e.side, frame=k, t_s=round(float(tvec[k - 1]), 1), best_IoU=bi, matched=int(bi >= 0.5)))
    print(eid, 'done', 'recon' if recon else '', 'bleach ratio %.2f, GT-vs-pipeline bleached IoU %.2f' % (rmin, chk), flush=True)

F = pd.DataFrame(rows_f); P = pd.DataFrame(rows_p); B = pd.DataFrame(rows_b); AP = pd.DataFrame(rows_ap); EM = pd.DataFrame(rows_e)
def block(f, p, b, a, name):
    TP, FP, FN = f.TP.sum(), f.FP.sum(), f.FN.sum()
    g = a.groupby('threshold')[['TP', 'FP', 'FN']].sum(); apv = g.TP / (g.TP + g.FP + g.FN)
    return {'group': name, 'embryos': f.embryo.nunique(), 'frames': len(f), 'nuclei_GT': int(TP + FN), 'nuclei_method': int(TP + FP),
            'TP': int(TP), 'FP': int(FP), 'FN': int(FN), 'precision': TP / max(TP + FP, 1), 'recall': TP / max(TP + FN, 1),
            'F1': 2 * TP / max(2 * TP + FP + FN, 1), 'IoU_mean': p.IoU.mean(), 'IoU_SD': p.IoU.std(), 'IoU_median': p.IoU.median(),
            'IoU_Q1': p.IoU.quantile(.25), 'IoU_Q3': p.IoU.quantile(.75), 'Dice_mean': p.Dice.mean(), 'pixel_Dice_mean': f.pixel_Dice.mean(),
            'area_diff_um2_mean': p.area_diff_um2.mean(), 'area_diff_um2_SD': p.area_diff_um2.std(),
            'centroid_shift_um_mean': p.centroid_shift_um.mean(), 'centroid_shift_um_SD': p.centroid_shift_um.std(),
            'AP@0.5': apv.get(0.5), 'AP@0.75': apv.get(0.75), 'AP@0.9': apv.get(0.9), 'mAP(0.5:0.95)': apv.mean(),
            'bleached_matched': int(b.matched.sum()), 'bleached_frames': len(b), 'bleached_IoU_mean': b.best_IoU.mean()}
S = []
for meth in F.method.unique():
    sel = lambda d, q: d[(d.method == meth) & q(d)]
    for name, q in [('Good', lambda d: d.case == 'Good'), ('Bad', lambda d: d.case == 'Bad'), ('Ventral', lambda d: d.side == 'Ventral'),
                    ('Dorsal', lambda d: d.side == 'Dorsal'), ('All', lambda d: d.case.notna())]:
        r = block(sel(F, q), sel(P, q), sel(B, q), sel(AP, q), name); r = {'method': meth, **r}; S.append(r)
S = pd.DataFrame(S)
PE = []
for (meth, eid), f in F.groupby(['method', 'embryo'], sort=False):
    r = block(f, P[(P.method == meth) & (P.embryo == eid)], B[(B.method == meth) & (B.embryo == eid)], AP[(AP.method == meth) & (AP.embryo == eid)], eid)
    PE.append({'method': meth, 'embryo': eid, 'case': f.case.iloc[0], 'side': f.side.iloc[0], **{k: v for k, v in r.items() if k not in ('group', 'embryos')}})
PE = pd.DataFrame(PE)
APc = AP.groupby(['method', 'threshold'])[['TP', 'FP', 'FN']].sum().reset_index(); APc['AP'] = APc.TP / (APc.TP + APc.FP + APc.FN)
with pd.ExcelWriter(f'{OUT}/GT_results.xlsx') as w:
    S.to_excel(w, sheet_name='Summary', index=False); PE.to_excel(w, sheet_name='Per_embryo', index=False)
    EM.to_excel(w, sheet_name='Embryos', index=False); F.to_excel(w, sheet_name='Per_frame', index=False)
    B.to_excel(w, sheet_name='Bleached_nucleus', index=False); APc.to_excel(w, sheet_name='AP_vs_IoU_threshold', index=False)
    P.to_excel(w, sheet_name='Matched_nuclei', index=False)
S.to_csv(f'{OUT}/GT_summary.csv', index=False); PE.to_csv(f'{OUT}/GT_per_embryo.csv', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40)
print(S[['method', 'group', 'embryos', 'TP', 'FP', 'FN', 'F1', 'IoU_median', 'area_diff_um2_mean', 'AP@0.75', 'bleached_matched', 'bleached_frames']].round(3).to_string(index=False))
print(EM.to_string(index=False))
