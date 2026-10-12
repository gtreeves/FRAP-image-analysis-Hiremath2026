"""Fig 14: comparison with the ground truth (10 embryos, 533 frames). (a) average precision vs IoU threshold,
(b) area difference of each matched nucleus, (c) bleached nucleus segmented and tracked.
Input tables (in figure_data/, written by the ground-truth package, ground_truth/.../code/scoring/score_gt_all.py
and frap_track.py): GT_AP_vs_IoU_threshold.csv, GT_matched_nuclei.csv, GT_summary.csv, GT_bleached_track_summary.csv."""
import os, numpy as np, pandas as pd, matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from plos_style import setup, save

STYLE = {'font.size': 8, 'axes.linewidth': 1.2, 'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.major.width': 1.0,
         'ytick.major.width': 1.0, 'axes.labelweight': 'bold', 'axes.labelsize': 9}


def make(data_dir, out_dir):
    setup(**STYLE)
    M = ['FRAP pipeline', 'CellProfiler', 'Cellpose 3', 'StarDist']; LAB = ['Our pipeline', 'CellProfiler', 'Cellpose 3', 'StarDist']
    TL = ['Our\npipeline', 'Cell-\nProfiler', 'Cellpose\n3', 'Star-\nDist']
    COL = ['#d81b8c', '#2ca02c', '#1f77b4', '#ff7f0e']
    AP = pd.read_csv(os.path.join(data_dir, 'GT_AP_vs_IoU_threshold.csv')); MN = pd.read_csv(os.path.join(data_dir, 'GT_matched_nuclei.csv')); S = pd.read_csv(os.path.join(data_dir, 'GT_summary.csv'))
    BT = pd.read_csv(os.path.join(data_dir, 'GT_bleached_track_summary.csv'))
    fig, ax = plt.subplots(1, 3, figsize=(6.9, 2.25), gridspec_kw=dict(width_ratios=[1, 1.15, 1.15], wspace=0.42, left=0.07, right=0.995, top=0.86, bottom=0.23))
    # a: AP vs IoU threshold
    a = ax[0]
    for m, l, c in zip(M, LAB, COL):
        d = AP[AP.method == m].sort_values('threshold'); a.plot(d.threshold, d.AP, '-o', color=c, lw=1.6, ms=3, label=l)
    a.set_xlabel('IoU threshold'); a.set_ylabel('Average precision'); a.set_xlim(0.48, 0.97); a.set_ylim(0, 1.03)
    a.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9]); a.legend(frameon=False, fontsize=7.4, loc='lower left', handlelength=1.2, borderaxespad=0.2)
    # b: area difference per matched nucleus
    b = ax[1]; data = [MN[MN.method == m].area_diff_um2.values for m in M]
    bp = b.boxplot(data, widths=0.55, showfliers=False, patch_artist=True, medianprops=dict(color='k', lw=1.2), whiskerprops=dict(lw=1), capprops=dict(lw=1))
    for p, c in zip(bp['boxes'], COL): p.set_facecolor(matplotlib.colors.to_rgba(c, 0.55)); p.set_edgecolor(c); p.set_linewidth(1.2)
    b.axhline(0, color='0.4', lw=0.8, ls='--'); b.set_xticks(range(1, 5)); b.set_xticklabels(TL, fontsize=7.4)
    b.set_ylabel('Area − GT area (µm²)', fontsize=9.5)
    # c: bleached nucleus, segmented vs followed
    c = ax[2]; seg = [S[(S.method == m) & (S.group == 'All')].bleached_matched.iloc[0] / 523 * 100 for m in M]
    row = {'FRAP pipeline': 'FRAP pipeline (L_track)', 'CellProfiler': 'CellProfiler', 'Cellpose 3': 'Cellpose 3', 'StarDist': 'StarDist'}
    fol = []
    for m in M:
        r = BT[BT.method == row[m]].iloc[0]; fol.append((r['bleached nucleus followed (frames)'] - r['frames on the wrong nucleus']) / 523 * 100)
    x = np.arange(4); w = 0.38
    c.bar(x - w / 2, seg, w, color=[matplotlib.colors.to_rgba(cc, 0.55) for cc in COL], edgecolor=COL, lw=1.0, label='segmented')
    c.bar(x + w / 2, fol, w, color=COL, edgecolor='k', lw=1.0, hatch='////', label='tracked')
    for xi, s_, f_ in zip(x, seg, fol): c.text(xi + w / 2, f_ + 1.5, f'{f_:.0f}', ha='center', va='bottom', fontsize=7.4)
    c.set_xticks(x); c.set_xticklabels(TL, fontsize=7.4); c.set_ylim(0, 112); c.set_yticks([0, 25, 50, 75, 100])
    c.set_ylabel('Bleached nucleus\n(% of frames)', linespacing=1.0)
    c.legend(handles=[Patch(facecolor='0.75', edgecolor='0.5', label='segmented'), Patch(facecolor='0.75', edgecolor='k', hatch='////', label='tracked')],
             frameon=False, fontsize=7.4, loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.2), handlelength=1.2, columnspacing=1)
    for a_, L in zip(ax, 'abc'):
        a_.text(-0.33 if L == 'c' else -0.27, 1.12, L, transform=a_.transAxes, fontsize=10, va='top')
        for s in ['top', 'right']: a_.spines[s].set_visible(True)
    pd.DataFrame({'method': LAB, 'bleached_segmented_pct': seg, 'bleached_tracked_pct': fol}).to_csv(os.path.join(data_dir, 'Fig14c_data.csv'), index=False)
    return save(fig, out_dir, 'Fig14', 7.5)
