"""FRAP recovery and tracking consistency for each segmentation method, on the ground-truth embryos.

FRAP (step 1): for every method the bleached nucleus is measured the way the FRAP pipeline does it
(analyze_bleach_czi.m, dl_vs_t.m, fit_bleach_dl.m, ftn_FRAP.m):
  - nuclear and cytoplasmic Dorsal-GFP intensity = mean of the pixels between the 10th and 90th percentile
  - cytoplasm of a nucleus = non-nuclear pixels whose nearest nucleus it is (nuclei dilated by 1 px,
    at most 25 px from the nucleus), the same idea as the pipeline's watershed territories
  - normalised by the mean of the 6 nearest other nuclei (and their cytoplasm) in each frame:
    y = y(t) / mean_others(t) * mean_others(frame 1)
  - fitted from frame 2 on with y(t) = exp(-k_out t) (c0 + k_in * int_0^t ycyt(s) exp(k_out s) ds),
    ycyt and ynuc smoothed with a 20-point moving average, one outlier-removal pass (3 SD)
The bleached nucleus of each method:
  ground truth   the GT bleached nucleus (largest mRFP drop frame 1->2, then followed by overlap)
  FRAP pipeline  its own tracked label (L_nuc = final masks as used for FRAP; L_track also reported)
  other methods  no tracking of their own: the object that overlaps the pipeline's bleached nucleus in
                 frame 1, then in every frame the object with the largest IoU with the last one found
                 (IoU >= 0.3); otherwise that frame is missing
The fit window is the frames with ground truth (all 60 except E02: 1-39, E20: 1-14).

Tracking (step 2): ground-truth nuclei are linked between consecutive frames (one-to-one, IoU >= 0.5).
For every such step each method is checked: is the nucleus found in both frames (IoU >= 0.5 with GT),
does it keep its identity (pipeline: same label; other methods: linked by the same one-to-one IoU rule),
and how much does its area change compared with the ground truth's own change.
"""
import sys, os, numpy as np, pandas as pd
from scipy import ndimage as ndi
from scipy.optimize import least_squares
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h5read import read
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'score_gt_all.py')).read()
MIN_AREA = 200
exec(src[src.index('def gt_instances'):src.index('THR = ')])

IN, OUT = sys.argv[1], sys.argv[2]; os.makedirs(OUT, exist_ok=True)
E = pd.read_csv(f'{IN}/gt_embryos.csv')
OTHER = [('CellProfiler', 'CP'), ('Cellpose 3', 'Cellpose3'), ('StarDist', 'StarDist')]

def tmean(I, m):
    v = I[m].astype(float)
    if v.size < 5: return np.nan
    lo, hi = np.percentile(v, [10, 90]); v = v[(v > lo) & (v < hi)]
    return v.mean() if v.size else np.nan

def mavg(y, span=20):
    span = span - 1 if span % 2 == 0 else span; h = span // 2; n = len(y); out = np.empty(n)
    for i in range(n):
        k = min(h, i, n - 1 - i); out[i] = np.mean(y[i - k:i + k + 1])
    return out

def ftn(P, t, ycyt):
    c0, kin, kout = P
    yI = ycyt * np.exp(kout * t)
    I = np.concatenate([[0], np.cumsum(0.5 * (yI[1:] + yI[:-1]) * np.diff(t))])
    return np.exp(-kout * t) * (c0 + kin * I)

def fit(t, ynuc, ycyt):
    ok = np.isfinite(ynuc) & np.isfinite(ycyt); t, ynuc, ycyt = t[ok], ynuc[ok], ycyt[ok]
    if len(t) < 6: return dict(c0=np.nan, k_in=np.nan, k_out=np.nan, R2=np.nan, n_fit=len(t))
    ycf = mavg(ycyt)
    for jj in range(2):
        ys = mavg(ynuc)
        lb = [np.min(ynuc) / 2, 0, 0]; ub = [2 * np.max(ynuc), np.inf, np.inf]
        x0 = [min(max(ynuc[0], lb[0] * 1.0001), ub[0] * 0.9999), 0.01, 0.01]
        r = least_squares(lambda P: ftn(P, t, ycf) - ys, x0, bounds=(lb, ub), max_nfev=3000, xtol=1e-12, ftol=1e-12)
        res = r.fun
        if jj == 0:
            v = np.abs(res - res.mean()) > 3 * res.std()
            if v.any(): t, ynuc, ycf = t[~v], ynuc[~v], ycf[~v]
            else: break
    resnorm = np.sum(r.fun ** 2)
    R2 = 1 - resnorm / np.sum((ynuc - ynuc.mean()) ** 2)
    return dict(c0=r.x[0], k_in=r.x[1], k_out=r.x[2], R2=R2, n_fit=len(t))

def cyto(L, labs, maxd=25):
    nuc = ndi.binary_dilation(L > 0)
    d, (ii, jj) = ndi.distance_transform_edt(L == 0, return_indices=True)
    near = L[ii, jj]; bg = (~nuc) & (d <= maxd)
    return {l: bg & (near == l) for l in labs}

def centroids(L):
    labs = np.unique(L[L > 0]);
    if len(labs) == 0: return labs, np.zeros((0, 2))
    return labs, np.array(ndi.center_of_mass(np.ones_like(L), L, labs))

def measure(Lk, I, lab, nn=6):
    """nuclear/cyt intensity of label lab and mean of its nn nearest neighbours in frame k."""
    if not lab: return [np.nan] * 4 + [0]
    labs, C = centroids(Lk)
    if lab not in labs: return [np.nan] * 4 + [0]
    c = C[list(labs).index(lab)]; d = np.hypot(*(C - c).T)
    others = [int(l) for l in labs[np.argsort(d)] if l != lab][:nn]
    cy = cyto(Lk, [lab] + others)
    yn, yc = tmean(I, Lk == lab), tmean(I, cy[lab])
    on = np.nanmean([tmean(I, Lk == l) for l in others]); oc = np.nanmean([tmean(I, cy[l]) for l in others])
    return [yn, yc, on, oc, int((Lk == lab).sum())]

def link_track(stack_get, frames, start_mask, thr=0.3):
    """naive overlap tracking of one nucleus: returns {frame: label or 0}."""
    out = {}; prev = start_mask
    for k in frames:
        L = stack_get(k); ov = L[prev]; ov = ov[ov > 0]
        lab = 0
        if ov.size:
            u, it = np.unique(ov, return_counts=True); aB = np.bincount(L.ravel())
            iou = it / (prev.sum() + aB[u] - it); j = np.argmax(iou)
            if iou[j] >= thr: lab = int(u[j])
        out[k] = lab
        if lab: prev = L == lab
    return out

def match_all(A, B, thr=0.5):
    C = candidates(A, B); ua, ub, M = set(), set(), {}
    for iou, a, b, it in C['c']:
        if iou < thr: break
        if a in ua or b in ub: continue
        ua.add(a); ub.add(b); M[int(a)] = int(b)
    return M

rows_c, rows_fit, rows_bt, rows_tr = [], [], [], []
for _, e in E.iterrows():
    eid = e.id
    pm = f'{IN}/pipeline_masks/{eid}_pipeline_masks_reconstructed.mat'
    if not os.path.exists(pm): pm = f'{IN}/pipeline_masks/{eid}_pipeline_masks_noreconstruction.mat'
    tvec = read(pm, 't').ravel(); px = float(read(pm, 'scalings').ravel()[0])
    Idl = read(pm, 'I_dl'); Inuc = read(pm, 'I_nuc')
    Lt = read(pm, 'L_track').astype(np.int32); Ln = read(pm, 'L_nuc').astype(np.int32)
    gt = np.load(f'{IN}/stacks/{eid}_GT.npy'); FR = np.load(f'{IN}/stacks/{eid}_GT_frames.npy').tolist()
    G = {k: gt_instances(gt[i])[0] for i, k in enumerate(FR)}
    stacks = {'FRAP pipeline (L_nuc, as used)': lambda k: Ln[:, :, k - 1], 'FRAP pipeline (L_track)': lambda k: Lt[:, :, k - 1]}
    for m, suf in OTHER:
        X = np.load(f'{IN}/stacks/{eid}_{suf}.npy'); stacks[m] = (lambda k, X=X: X[k - 1].astype(np.int32))
    # GT bleached nucleus (same rule as the scoring)
    g1, g2 = G[FR[0]], G[FR[1]]; I1, I2 = Inuc[:, :, FR[0] - 1].astype(float), Inuc[:, :, FR[1] - 1].astype(float)
    ratios = {}
    for lab in np.unique(g1[g1 > 0]):
        if lab in touches(g1): continue
        m1 = g1 == lab; ov = g2[m1]; ov = ov[ov > 0]
        if ov.size == 0: continue
        l2 = int(np.bincount(ov).argmax()); ratios[(int(lab), l2)] = I2[g2 == l2].mean() / I1[m1].mean()
    (bl1, bl2), _ = min(ratios.items(), key=lambda x: x[1]); bgt = {FR[0]: bl1, FR[1]: bl2}
    for k in FR[2:]:
        prev = G[k - 1] == bgt.get(k - 1, -1); cur = G[k][prev]; cur = cur[cur > 0]
        bgt[k] = int(np.bincount(cur).argmax()) if cur.size else 0
    pl = int(e.bleached_label)
    # bleached-nucleus label per frame for each method
    BL = {'Ground truth': (lambda k: G[k], bgt)}
    BL['FRAP pipeline (L_nuc, as used)'] = (stacks['FRAP pipeline (L_nuc, as used)'], {k: pl if (Ln[:, :, k - 1] == pl).any() else 0 for k in FR})
    BL['FRAP pipeline (L_track)'] = (stacks['FRAP pipeline (L_track)'], {k: pl if (Lt[:, :, k - 1] == pl).any() else 0 for k in FR})
    start = Lt[:, :, FR[0] - 1] == pl
    NAIVE = {}
    for m, _ in OTHER:
        NAIVE[m] = link_track(stacks[m], FR, start)
        asst = {}
        for k in FR:                       # GT-assisted: the object that best matches the GT bleached nucleus
            L = stacks[m](k); gm = G[k] == bgt.get(k, -1); ov = L[gm]; ov = ov[ov > 0]; lab = 0
            if ov.size:
                u, it = np.unique(ov, return_counts=True); aB = np.bincount(L.ravel())
                iou = it / (gm.sum() + aB[u] - it); j = np.argmax(iou)
                if iou[j] >= 0.5: lab = int(u[j])
            asst[k] = lab
        BL[m] = (stacks[m], asst)
    curves = {}
    for meth, (get, bl) in BL.items():
        Y = np.array([measure(get(k), Idl[:, :, k - 1], bl[k]) for k in FR], float)
        yn, yc, on, oc, area = Y.T
        t = tvec[np.array(FR) - 1]
        rel = (yn / on) / (yn[0] / on[0])                          # fraction of pre-bleach (neighbour-normalised)
        ynuc = yn[1:] / on[1:] * on[0]; ycyt = yc[1:] / oc[1:] * oc[0]
        f = fit(t[1:] - t[1], ynuc, ycyt)
        curves[meth] = rel
        # bleached-track checks against GT (other methods: their own naive overlap tracking)
        tb = NAIVE.get(meth, bl)
        area_t = np.array([(get(k) == tb[k]).sum() if tb[k] else 0 for k in FR], float) if meth in NAIVE else area
        wrong = missing = 0; ious = []
        for k in FR[1:]:
            L = get(k); lab = tb[k]
            if not lab: missing += 1; continue
            if meth == 'Ground truth': continue
            gm = G[k] == bgt.get(k, -1); mm = L == lab
            iou = (gm & mm).sum() / max((gm | mm).sum(), 1); ious.append(iou)
            if iou < 0.5: wrong += 1
        a = area_t[1:] * px ** 2; a = np.where(a > 0, a, np.nan)
        da = np.abs(np.diff(a)) / a[:-1]
        gaps = np.diff(np.concatenate([[0], (np.array([tb[k] for k in FR[1:]]) == 0).astype(int), [0]]))
        rows_bt.append(dict(embryo=eid, case=e.case, side=e.side, method=meth, post_bleach_frames=len(FR) - 1,
                            tracking='own label' if meth not in NAIVE else 'simple overlap linking',
                            frames_segmented_GT_assisted=int(sum(1 for k in FR[1:] if bl[k])),
                            frames_found=len(FR) - 1 - missing, frames_missing=missing, gap_runs=int((gaps == 1).sum()),
                            frames_wrong_nucleus=wrong if meth != 'Ground truth' else np.nan,
                            mean_IoU_with_GT_bleached=np.mean(ious) if ious else np.nan,
                            area_mean_um2=np.nanmean(a), area_CV=np.nanstd(a) / np.nanmean(a),
                            max_area_jump_pct=100 * np.nanmax(da) if np.isfinite(da).any() else np.nan,
                            steps_area_jump_gt20pct=int(np.nansum(da > 0.2))))
        f['k_in_over_k_out'] = f['k_in'] / f['k_out'] if f['k_out'] else np.nan
        rows_fit.append(dict(embryo=eid, case=e.case, side=e.side, method=meth, fit_frames=f'{FR[0]}-{FR[-1]}', **f,
                             t_half_s=np.log(2) / f['k_out'] if f['k_out'] and np.isfinite(f['k_out']) else np.nan,
                             recovery_end=np.nanmean(rel[-5:])))
        for k, tt, r_, a_ in zip(FR, t, rel, area):
            rows_c.append(dict(embryo=eid, method=meth, frame=k, t_s=tt - t[1], rel_intensity=r_, area_um2=a_ * px ** 2 if a_ else np.nan))
    for meth in BL:
        if meth == 'Ground truth': continue
        d = curves[meth][1:] - curves['Ground truth'][1:]
        rows_fit[[i for i, r in enumerate(rows_fit) if r['embryo'] == eid and r['method'] == meth][0]]['curve_RMSE_vs_GT'] = np.sqrt(np.nanmean(d ** 2))
    # ---------- tracking consistency, all nuclei ----------
    gl = {}
    for k0, k1 in zip(FR[:-1], FR[1:]):
        if k1 != k0 + 1: continue
        bd = touches(G[k0]) | touches(G[k1])
        gl[k0] = {a: b for a, b in match_all(G[k0], G[k1]).items() if a not in touches(G[k0]) and b not in touches(G[k1])}
    meths = ['FRAP pipeline (L_track)'] + [m for m, _ in OTHER]
    for meth in meths:
        get = stacks[meth]
        mt = {k: match_all(G[k], get(k)) for k in FR}   # GT label -> method label (no edge filter; edge GT links are skipped)
        own = {}
        if meth not in ('FRAP pipeline (L_track)',):
            for k0 in gl:
                own[k0] = match_all(get(k0), get(k0 + 1))
        n = cov = keep = 0; dA_m, dA_g = [], []
        brk = {}
        for k0, links in gl.items():
            A0, A1 = get(k0), get(k0 + 1); aA0, aA1 = np.bincount(A0.ravel()), np.bincount(A1.ravel())
            g0, g1_ = np.bincount(G[k0].ravel()), np.bincount(G[k0 + 1].ravel())
            for ga, gb in links.items():
                n += 1
                ma, mb = mt[k0].get(ga), mt[k0 + 1].get(gb)
                ok = ma is not None and mb is not None
                same = ok and ((ma == mb) if meth == 'FRAP pipeline (L_track)' else own[k0].get(ma) == mb)
                cov += ok; keep += same
                if ok:
                    dA_m.append((aA1[mb] - aA0[ma]) / aA0[ma]); dA_g.append((g1_[gb] - g0[ga]) / g0[ga])
        dA_m, dA_g = np.array(dA_m), np.array(dA_g)
        rows_tr.append(dict(embryo=eid, case=e.case, side=e.side, method=meth, GT_links=n,
                            found_both_frames=cov, identity_kept=keep,
                            pct_found=100 * cov / n, pct_identity_kept_of_found=100 * keep / max(cov, 1),
                            median_abs_area_change_pct=100 * np.median(np.abs(dA_m)) if len(dA_m) else np.nan,
                            GT_median_abs_area_change_pct=100 * np.median(np.abs(dA_g)) if len(dA_g) else np.nan,
                            pct_area_jumps_gt20=100 * np.mean(np.abs(dA_m) > 0.2) if len(dA_m) else np.nan,
                            GT_pct_area_jumps_gt20=100 * np.mean(np.abs(dA_g) > 0.2) if len(dA_g) else np.nan,
                            median_abs_area_change_error_pct=100 * np.median(np.abs(dA_m - dA_g)) if len(dA_m) else np.nan))
    print(eid, 'done', flush=True)

C = pd.DataFrame(rows_c); F = pd.DataFrame(rows_fit); B = pd.DataFrame(rows_bt); T = pd.DataFrame(rows_tr)
g = F[F.method == 'Ground truth'].set_index('embryo')
for col in ('k_in', 'k_out', 't_half_s', 'k_in_over_k_out'):
    F[f'{col}_ratio_to_GT'] = F[col] / F.embryo.map(g[col])
C.to_csv(f'{OUT}/frap_curves.csv', index=False); F.to_csv(f'{OUT}/frap_fits.csv', index=False)
B.to_csv(f'{OUT}/bleached_track.csv', index=False); T.to_csv(f'{OUT}/tracking_all_nuclei.csv', index=False)
print('saved')
