import logging; logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
"""Failure examples against the ground truth.  Usage: python make_fail_gt.py IN OUT
IN: gt_embryos.csv, stacks/ (GT, CP, Cellpose3, StarDist), pipeline_masks/.
Each example: A = whole field (method vs ground truth), B = zoom (method vs ground truth),
C = the same zoom for the FRAP pipeline vs ground truth."""
import sys, os, json, numpy as np, pandas as pd, matplotlib
from scipy import ndimage as ndi
from PIL import Image
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from h5read import read
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'score_gt_all.py')).read()
MIN_AREA = 200; exec(src[src.index('def gt_instances'):src.index('THR = ')])
plt.rcParams.update({'font.family': ['Liberation Sans', 'Arial', 'DejaVu Sans'], 'pdf.fonttype': 42})
IN, OUT = sys.argv[1], sys.argv[2]
E = pd.read_csv(f'{IN}/gt_embryos.csv').set_index('id')
cGT = np.array([1.0, 0.85, 0.0])
COL = {'FRAP pipeline': np.array([0.85, 0.11, 0.55]), 'CellProfiler': np.array([0.17, 0.63, 0.17]),
       'Cellpose 3': np.array([0.12, 0.47, 0.71]), 'StarDist': np.array([1.0, 0.50, 0.05])}
SUF = {'CellProfiler': 'CP', 'Cellpose 3': 'Cellpose3', 'StarDist': 'StarDist'}

EX = {
 'CellProfiler': [
  dict(e='E02', f=10, focus='bleached', tag='bleached_nucleus_missed', title='bleached nucleus not detected'),
  dict(e='E20', f=2, focus='bleached', tag='bleached_nucleus_missed_after_bleach', title='bleached nucleus not detected right after bleaching'),
  dict(e='E13', f=30, focus='bleached', tag='bleached_nucleus_missed', title='bleached nucleus not detected'),
  dict(e='E13', f=49, focus='bleached', tag='bleached_nucleus_undersized', title='bleached nucleus outlined far too small'),
  dict(e='E02', f=39, focus='fnfp', tag='nuclei_missed', title='nuclei not detected'),
  dict(e='E16', f=1, focus='worst', tag='oversized_outlines_dorsal', title='outlines too wide (inverted Dorsal channel)')],
 'Cellpose 3': [
  dict(e='E11', f=25, focus='bleached', tag='bleached_and_small_nucleus_missed', title='bleached nucleus and another small, dim nucleus not detected'),
  dict(e='E20', f=2, focus='bleached', tag='bleached_nucleus_missed_after_bleach', title='bleached nucleus not detected right after bleaching'),
  dict(e='E22', f=11, focus='worst', tag='undersized_outline', title='nucleus outlined far too small'),
  dict(e='E16', f=1, focus='worst', tag='oversized_outlines_dorsal', title='outlines too wide (inverted Dorsal channel)')],
 'StarDist': [
  dict(e='E02', f=5, focus='bleached', tag='bleached_nucleus_missed_after_bleach', title='bleached nucleus not detected after bleaching'),
  dict(e='E20', f=4, focus='bleached', tag='bleached_nucleus_missed', title='bleached nucleus not detected'),
  dict(e='E19', f=4, focus='fp', tag='extra_objects_in_cytoplasm', title='extra objects in the cytoplasm between nuclei (inverted Dorsal channel)'),
  dict(e='E22', f=20, focus='fnfp', tag='oversized_outlines_and_missed', title='outlines far too wide, scored as missed + extra'),
  dict(e='E16', f=10, focus='worst', tag='oversized_outlines_dorsal', title='outlines too wide (inverted Dorsal channel)')]}

def band_edges(L, r, inner):
    st = ndi.iterate_structure(ndi.generate_binary_structure(2, 1), r); out = np.zeros(L.shape, bool)
    for sl, lab in zip(ndi.find_objects(L), range(1, int(L.max()) + 1)):
        if sl is None: continue
        sl = tuple(slice(max(s.start - r - 1, 0), s.stop + r + 1) for s in sl); m = L[sl] == lab
        out[sl] |= (m & ~ndi.binary_erosion(m, st)) if inner else (ndi.binary_dilation(m, st) & ~m)
    return out

def render(I, G, B, lims, U, lw, col, fnG=(), fpB=()):
    g = np.clip((I - lims[0]) / (lims[1] - lims[0]), 0, 1)
    g = np.array(Image.fromarray((g * 255).astype(np.uint8)).resize((I.shape[1] * U, I.shape[0] * U), Image.BILINEAR)) / 255.
    Gu = np.kron(G, np.ones((U, U), int)); Bu = np.kron(B, np.ones((U, U), int)); rgb = np.dstack([g, g, g])
    r = max(1, int(round(lw * U)))
    for l in fnG: m = Gu == l; rgb[m] = rgb[m] * 0.5 + cGT * 0.5
    for l in fpB: m = Bu == l; rgb[m] = rgb[m] * 0.5 + col * 0.5
    rgb[band_edges(Bu, r, False)] = col; rgb[band_edges(Gu, r, True)] = cGT
    return rgb

def scalebar(rgb, px, U, um, margin=0.04):
    H, W, _ = rgb.shape; L = int(round(um / px * U)); h = max(3, int(H * 0.013)); m = int(W * margin)
    rgb[H - m - h:H - m, W - m - L:W - m] = 1.0; return rgb

def relabel(L):
    return ndi.label(L > 0)[0] if L.max() == 0 else L

cache = {}
def load(eid):
    if eid in cache: return cache[eid]
    pm = f'{IN}/pipeline_masks/{eid}_pipeline_masks_reconstructed.mat'
    if not os.path.exists(pm): pm = f'{IN}/pipeline_masks/{eid}_pipeline_masks_noreconstruction.mat'
    e = E.loc[eid]; dors = 'Dorsal' in e.segmented_channel
    D = dict(Lp=read(pm, 'L_track').astype(np.int32), I=read(pm, 'I_dl' if dors else 'I_nuc').astype(float), Inuc=read(pm, 'I_nuc').astype(float),
             px=float(read(pm, 'scalings').ravel()[0]), t=read(pm, 't').ravel(), dors=dors)
    if dors: D['I'] = 65535 - D['I']
    gt = np.load(f'{IN}/stacks/{eid}_GT.npy'); FR = np.load(f'{IN}/stacks/{eid}_GT_frames.npy').tolist()
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
    D.update(G=G, bgt=bgt, FR=FR, M={m: np.load(f'{IN}/stacks/{eid}_{s}.npy') for m, s in SUF.items()})
    cache[eid] = D; return D

index = []
for meth, items in EX.items():
    od = f'{OUT}/{meth.replace(" ", "_")}'; os.makedirs(od, exist_ok=True); col = COL[meth]
    for it in items:
        ids = [i for i in E.index if i.startswith(it['e'] + '_')]
        if not ids: print('skipped (embryo not in the list):', it['e'], it['tag']); continue
        eid = ids[0]; e = E.loc[eid]; D = load(eid); f = it['f']
        G = D['G'][f]; B = D['M'][meth][f - 1].astype(np.int32); P = D['Lp'][:, :, f - 1]; I = D['I'][:, :, f - 1]
        lims = np.percentile(D['I'][:, :, 0], [0.5, 99.7])
        M, fn, fp = match(G, B); Mp, fnp, fpp = match(G, P)
        bl = D['bgt'][f]; bi = best_iou(G, bl, B) if bl else np.nan; bip = best_iou(G, bl, P) if bl else np.nan
        if it['focus'] == 'bleached': cy, cx = ndi.center_of_mass(G == bl); mark = G == bl
        elif it['focus'] == 'worst':
            w = min(M, key=lambda x: x[2]); cy, cx = ndi.center_of_mass(G == w[0]); mark = G == w[0]
        else:
            objs = ([G == l for l in fn] if it['focus'] != 'fp' else []) + [B == l for l in fp]
            cs = [ndi.center_of_mass(o) for o in objs]; cy, cx = np.median([c[0] for c in cs]), np.median([c[1] for c in cs]); mark = None
        half = it.get('half', 95); cy = int(np.clip(cy, half, 512 - half)); cx = int(np.clip(cx, half, 512 - half)); y0, x0 = cy - half, cx - half
        sl = (slice(y0, y0 + 2 * half), slice(x0, x0 + 2 * half))
        inz = lambda L, ls: [l for l in ls if (L[sl] == l).any()]
        full = scalebar(render(I, G, B, lims, 4, 2.5, col, fn, fp), D['px'], 4, 10)
        zoom = scalebar(render(I[sl], G[sl], B[sl], lims, 10, 1.6, col, inz(G, fn), inz(B, fp)), D['px'], 10, 5)
        zp = scalebar(render(I[sl], G[sl], P[sl], lims, 10, 1.6, COL['FRAP pipeline'], inz(G, fnp), inz(P, fpp)), D['px'], 10, 5)
        tag = f"{it['e']}_frame{f:03d}_{it['tag']}"
        for nm, im in (('full', full), ('zoom', zoom), ('zoom_pipeline', zp)):
            Image.fromarray((im * 255).astype(np.uint8)).save(f'{od}/{tag}_{nm}.png', dpi=(600, 600))
        fig, ax = plt.subplots(1, 3, figsize=(10.6, 4.25), gridspec_kw=dict(wspace=0.05))
        for a, im in zip(ax, (full, zoom, zp)): a.imshow(im, interpolation='none'); a.set_axis_off()
        ax[0].add_patch(Rectangle((x0 * 4, y0 * 4), 2 * half * 4, 2 * half * 4, fill=False, ec='white', lw=1.3, ls=(0, (4, 3))))
        if mark is not None and mark.any():
            ly, lx = ndi.center_of_mass(mark); zy, zx = (ly - y0) * 10, (lx - x0) * 10
            for a in ax[1:]:
                a.annotate('', xy=(zx - 0.27 * half * 10, zy + 0.27 * half * 10), xytext=(zx - 0.6 * half * 10, zy + 0.6 * half * 10),
                           arrowprops=dict(arrowstyle='-|>', color='white', lw=2, mutation_scale=15))
        for a, l, s in zip(ax, 'ABC', (meth, meth + ' (zoom)', 'FRAP pipeline (same zoom)')):
            a.text(0.0, 1.015, l, transform=a.transAxes, fontsize=11, fontweight='bold', va='bottom')
            a.text(0.07, 1.017, s, transform=a.transAxes, fontsize=8.5, va='bottom')
        t = D['t'][f - 1]
        fig.suptitle(f"{meth} vs ground truth: {it['e']} ({e.case}, {e.side}), frame {f}, t = {t:.0f} s  –  {it['title']}", fontsize=10, y=0.99)
        bg = 'inverted Dorsal-GFP' if D['dors'] else 'mRFP'
        arrow = 'Arrow: bleached nucleus.' if it['focus'] == 'bleached' else ('Arrow: the nucleus with the lowest IoU.' if it['focus'] == 'worst' else '')
        leg = (f"Yellow: ground truth (inner edge). Coloured line: segmentation (outer edge). Yellow fill: nucleus missed; coloured fill: extra object.\n"
               f"Background: {bg} (the channel segmented).  Scale bars: 10 µm (A), 5 µm (B, C).  {arrow}\n"
               f"Frame: {meth} {len(M)} matched, {len(fn)} missed, {len(fp)} extra; bleached nucleus IoU {bi:.2f}.   "
               f"FRAP pipeline {len(Mp)} matched, {len(fnp)} missed, {len(fpp)} extra; bleached nucleus IoU {bip:.2f}.")
        fig.text(0.5, 0.015, leg, ha='center', va='bottom', fontsize=7.4, linespacing=1.5)
        fig.subplots_adjust(left=0.01, right=0.99, top=0.86, bottom=0.16)
        for ext in ('png', 'pdf'): fig.savefig(f'{od}/{tag}_figure.{ext}', dpi=600 if ext == 'png' else 300)
        plt.close(fig)
        index.append({'Method': meth, 'Folder': meth.replace(' ', '_'), 'Embryo': eid, 'Case': e.case, 'Side': e.side, 'Frame': f,
                      'Time after bleach (s)': round(float(t), 1), 'Problem shown': it['title'],
                      'Method: matched / missed / extra': f'{len(M)} / {len(fn)} / {len(fp)}', 'Method: bleached nucleus IoU': round(float(bi), 3),
                      'Pipeline: matched / missed / extra': f'{len(Mp)} / {len(fnp)} / {len(fpp)}', 'Pipeline: bleached nucleus IoU': round(float(bip), 3),
                      'File prefix': tag})
        print(meth, tag, len(M), len(fn), len(fp), round(float(bi), 2), '| pipe', len(Mp), len(fnp), len(fpp), round(float(bip), 2), flush=True)
X = pd.DataFrame(index); X.to_csv(f'{OUT}/Failure_examples_GT_index.csv', index=False, encoding='utf-8-sig'); X.to_excel(f'{OUT}/Failure_examples_GT_index.xlsx', index=False)
