"""Fig 13: typical failures of CellProfiler, Cellpose 3 and StarDist against the ground truth, with our pipeline on
the same frame. Input: the ground-truth data (gt_embryos.csv, stacks/, pipeline_masks/)."""
import os, numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from plos_style import setup, save
from gt_tools import load, eid_of, embryos, match
from tiles import *


def make(gt_dir, out_dir, data_dir):
    setup(**{'font.size': 8})
    EX = [('CellProfiler', 'E02', 2, 'bleached', 'Bleached nucleus not detected\nright after bleaching'),
          ('CellProfiler', 'E11', 2, 'bleached', 'Bleached nucleus only\npartly segmented'),
          ('Cellpose 3', 'E11', 25, 'bleached', 'Dim bleached nucleus not\ndetected later in the recovery'),
          ('CellProfiler', 'E02', 39, 'fn', 'Dim nucleus not detected'),
          ('Cellpose 3', 'E22', 29, ('g', 14), 'Nucleus split into\ntwo objects'),
          ('CellProfiler', 'E13', 49, 'bleached', 'Ragged, far too small outline\nof the bleached nucleus'),
          ('StarDist', 'E16', 10, 'worst', 'Outlines far larger\nthan the nuclei'),
          ('StarDist', 'E19', 4, 'fp', 'Extra objects in the\ncytoplasm between nuclei')]
    DPI = 300
    sn = lambda v: round(v * DPI) / DPI
    W = 6.9; NR = 4; tw, Hh = layout(W, 2, 3, NR, hdr=0.52)
    tw = 1.05; Hh = sn(0.12 + NR * (0.52 + tw) + 0.12)
    fig = plt.figure(figsize=(W, Hh)); info = []
    for i, (m, e, f, focus, title) in enumerate(EX):
        eid = eid_of(gt_dir, e); D = load(gt_dir, eid); G = D['G'][f]; B = D['M'][m][f - 1].astype(np.int32); P = D['Lp'][:, :, f - 1]
        img = (65535 - D['I'][:, :, f - 1]) if D['dors'] else D['I'][:, :, f - 1]; img0 = (65535 - D['I'][:, :, 0]) if D['dors'] else D['I'][:, :, 0]
        lims = np.percentile(img0, [0.5, 99.7]); M, fn, fp = match(G, B); arrow_lab = None
        if focus == 'bleached': arrow_lab = D['bgt'][f]
        elif focus == 'worst': arrow_lab = min(M, key=lambda x: x[2])[0]
        elif isinstance(focus, tuple): arrow_lab = focus[1]
        elif focus == 'fn':
            sc = [l for l in fn if not (G == l)[:2].any()]; arrow_lab = max(fn, key=lambda l: (G == l).sum())
        if arrow_lab is not None: cy, cx = ndi.center_of_mass(G == arrow_lab)
        else:
            objs = [G == l for l in fn] if focus == 'fn' else [B == l for l in fp]
            cs = [ndi.center_of_mass(o) for o in objs]; cy, cx = np.median([c[0] for c in cs]), np.median([c[1] for c in cs])
        y0, x0 = crop_box(cy, cx); ar = ((cy - y0) * U, (cx - x0) * U) if arrow_lab is not None else None
        r, c = divmod(i, 2); x = sn((W - (6 * tw + 4 * 0.03 + 0.20)) / 2 + c * (3 * tw + 0.06 + 0.20)); ytop = sn(Hh - 0.12 - r * (0.52 + tw))
        fig.text(x / W, ytop / Hh, 'abcdefgh'[i], fontsize=10, va='top', ha='left')
        fig.text((x + 0.17) / W, ytop / Hh, title + f'  ({m})', fontsize=8, va='top', ha='left', linespacing=1.05)
        tiles = [('Ground truth', tile(img, lims, D['dors'], y0, x0, ref=G)), ('Our pipeline', tile(img, lims, D['dors'], y0, x0, ref=G, seg=P)),
                 (m, tile(img, lims, D['dors'], y0, x0, ref=G, seg=B))]
        for j, (lab, rgb) in enumerate(tiles):
            draw(fig, [sn(x + j * (tw + 0.03)) / W, sn(ytop - 0.44 - tw) / Hh, tw / W, tw / Hh], rgb, lab, ar, bar=(j == 0))
        info.append(dict(panel='abcdefgh'[i], method=m, embryo=eid, side=embryos(gt_dir).loc[eid, 'side'], frame=f, minutes_after_bleach=round(float(D['t'][f - 1] - D['t'][1]) / 60, 1), failure=title.replace('\n', ' ')))
    pd.DataFrame(info).to_csv(os.path.join(data_dir, 'Fig13_panels.csv'), index=False)
    return save(fig, out_dir, 'Fig13', 7.5)
