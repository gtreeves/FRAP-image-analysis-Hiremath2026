"""S2 Fig: the same kinds of failure in 10 further embryos without ground truth, compared with our pipeline.
Input: figure_data/S2_frames/*.npz, one file per example frame (both channels, our pipeline's labels P, the other
tool's labels M, the bleached label; extracted from the tool-comparison packages in comparisons/)."""
import os, glob, numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from plos_style import setup, save
from gt_tools import match
from tiles import *


def make(frames_dir, out_dir, data_dir):
    setup(**{'font.size': 8})
    EX = [('CellProfiler', 'E01', 2, 'bleached', 'Bleached nucleus not detected\nright after bleaching'),
          ('StarDist', 'E01', 2, 'bleached', 'Bleached nucleus broken into\nfragments right after bleaching'),
          ('CellProfiler', 'E10', 5, 'bleached', 'Bleached nucleus only\npartly segmented'),
          ('Cellpose 3', 'E03', 51, 'fn', 'Dim nucleus not detected'),
          ('Cellpose 3', 'E12', 45, 'split', 'Nucleus split into\ntwo objects'),
          ('Cellpose 3', 'E01', 53, 'worst', 'Ragged, undersized outline\nof a small nucleus'),
          ('StarDist', 'E01', 16, 'worst', 'Outlines far larger\nthan the nuclei'),
          ('StarDist', 'E09', 17, 'fp', 'Extra objects in the\ncytoplasm between nuclei'),
          ('CellProfiler', 'E04', 58, 'worst', 'Ragged, undersized outline\nof a dim nucleus')]
    def channels(Z):
        P1 = Z['P1'] > 0; r = []
        for c in ('C1f', 'C2f'):
            a = Z[c].astype(float); r.append(a[P1].mean() / a[~P1].mean())
        return r
    DPI = 300
    sn = lambda v: round(v * DPI) / DPI
    W = 6.9; NR = 3; hdr = 0.62; tw, Hh = layout(W, 3, 2, NR, hdr=hdr)
    tw = 1.05; Hh = sn(0.12 + NR * (hdr + tw) + 0.12)
    fig = plt.figure(figsize=(W, Hh)); info = []
    for i, (m, e, f, focus, title) in enumerate(EX):
        fn_ = glob.glob(f'{frames_dir}/{e}_*_frame{f:03d}_{m.replace(" ", "_")}.npz')[0]; Z = np.load(fn_, allow_pickle=True)
        P = Z['P'].astype(np.int32); B = Z['M'].astype(np.int32); bl = int(Z['bleach']); dors = 'Dorsal' in str(Z['cp_input'])
        rr = channels(Z); cn = int(np.argmin(rr)) if dors else int(np.argmax(rr))     # dorsal: show Dl-GFP (nuclei dark); else mRFP (nuclei bright)
        img = Z[('C1', 'C2')[cn]].astype(float); lims = np.percentile(Z[('C1f', 'C2f')[cn]].astype(float), [0.5, 99.7])
        M, fnl, fpl = match(P, B); lab = None
        if focus == 'bleached': lab = bl
        elif focus == 'worst': lab = min(M, key=lambda x: x[2])[0]
        elif focus == 'split':
            best = None
            for p in np.unique(P[P > 0]):
                ov = B[P == p]; ov = ov[ov > 0]; u, c = np.unique(ov, return_counts=True); k = (c > 0.2 * (P == p).sum()).sum()
                if k >= 2 and (best is None or c.sum() > best[1]): best = (p, c.sum())
            lab = best[0]
        elif focus == 'fn': lab = max(fnl, key=lambda l: (P == l).sum())
        if lab is not None: cy, cx = ndi.center_of_mass(P == lab)
        else:
            objs = [P == l for l in fnl] if focus == 'fn' else [B == l for l in fpl]
            cs = [ndi.center_of_mass(o) for o in objs]; cy, cx = np.median([c[0] for c in cs]), np.median([c[1] for c in cs])
        y0, x0 = crop_box(cy, cx); ar = ((cy - y0) * U, (cx - x0) * U) if lab is not None else None
        r, c = divmod(i, 3); x = sn((W - (6 * tw + 3 * 0.03 + 2 * 0.20)) / 2 + c * (2 * tw + 0.03 + 0.20)); ytop = sn(Hh - 0.12 - r * (hdr + tw))
        fig.text(x / W, ytop / Hh, 'abcdefghi'[i], fontsize=10, va='top', ha='left')
        fig.text((x + 0.17) / W, ytop / Hh, title + f'\n({m})', fontsize=8, va='top', ha='left', linespacing=1.05)
        tiles = [('Our pipeline', tile(img, lims, dors, y0, x0, seg=P)), (m, tile(img, lims, dors, y0, x0, seg=B))]
        for j, (lb, rgb) in enumerate(tiles):
            draw(fig, [sn(x + j * (tw + 0.03)) / W, sn(ytop - hdr + 0.08 - tw) / Hh, tw / W, tw / Hh], rgb, lb, ar, bar=(j == 0))
        info.append(dict(panel='abcdefghi'[i], method=m, embryo=e, frame=f, channel_shown='Dl-GFP' if dors else 'Histone-RFP', failure=title.replace('\n', ' ')))
    pd.DataFrame(info).to_csv(os.path.join(data_dir, 'S2_Fig_panels.csv'), index=False)
    return save(fig, out_dir, 'S2_Fig', 7.5)
