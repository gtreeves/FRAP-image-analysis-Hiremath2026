import logging; logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
import sys, os, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': ['Liberation Sans', 'Arial', 'DejaVu Sans'], 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
OUT = sys.argv[1] if len(sys.argv) > 1 else 'out'
X = pd.read_excel(f'{OUT}/GT_results.xlsx', sheet_name=None)
P, APc, B, S = X['Matched_nuclei'], X['AP_vs_IoU_threshold'], X['Bleached_nucleus'], X['Summary']
COL = {'FRAP pipeline': '#d81b8c', 'CellProfiler': '#2ca02c', 'Cellpose 3': '#1f77b4', 'StarDist': '#ff7f0e'}
meths = [m for m in COL if m in set(P.method)]
bins = np.arange(0.5, 1.0001, 0.01)
# 1) IoU histograms: rows = methods, cols = Good / Bad / All
fig, axs = plt.subplots(len(meths), 3, figsize=(7.2, 2.3 * len(meths)), sharex=True, squeeze=False)
tops = np.zeros((len(meths), 3))
for i, m in enumerate(meths):
    for j, g in enumerate(['Good', 'Bad', 'All']):
        ax = axs[i, j]; v = P[(P.method == m) & ((P.case == g) | (g == 'All'))].IoU.to_numpy()
        h, _, _ = ax.hist(v, bins=bins, color='0.35', edgecolor='white', linewidth=0.4); tops[i, j] = h.max()
        med = np.median(v); ax.plot([med, med], [0, h.max() * 1.04], 'k--', lw=1.1)
        ax.text(0.04, 0.95, f'n = {len(v):,}\nmedian = {med:.3f}\nmean ± SD = {v.mean():.3f} ± {v.std(ddof=1):.3f}', transform=ax.transAxes, va='top', fontsize=7)
        ax.spines[['top', 'right']].set_visible(False); ax.tick_params(labelsize=7); ax.set_xlim(0.5, 1)
        if i == 0: ax.set_title(f'{g} cases', fontsize=10, fontweight='bold')
        if j == 0: ax.set_ylabel(f'{m}\n\nNumber of nuclei', fontsize=8.5)
        if i == len(meths) - 1: ax.set_xlabel('IoU with ground truth', fontsize=8.5)
for j in range(3):
    for i in range(len(meths)): axs[i, j].set_ylim(0, tops[:, j].max() * 1.6)
fig.tight_layout(); [fig.savefig(f'{OUT}/GT_IoU_histograms.{x}', dpi=600 if x == 'png' else None) for x in ('png', 'pdf', 'svg')]; plt.close(fig)
# 2) AP vs threshold, and bleached-nucleus detection per embryo
fig, ax = plt.subplots(1, 2, figsize=(7.6, 3.0))
for m in meths:
    a = APc[APc.method == m]; ax[0].plot(a.threshold, a.AP, 'o-', ms=4, color=COL[m], label=m)
ax[0].set(xlabel='IoU threshold', ylabel='AP = TP / (TP + FP + FN)', ylim=(0, 1.02)); ax[0].legend(frameon=False, fontsize=8)
emb = list(dict.fromkeys(B.embryo)); w = 0.8 / len(meths)
for i, m in enumerate(meths):
    r = [B[(B.method == m) & (B.embryo == e)].matched.mean() * 100 for e in emb]
    ax[1].bar(np.arange(len(emb)) + (i - (len(meths) - 1) / 2) * w, r, w, color=COL[m], label=m)
ax[1].set_xticks(range(len(emb))); ax[1].set_xticklabels([e[:3] for e in emb], fontsize=7.5)
ax[1].set(ylabel='Bleached nucleus matched (% of frames)', ylim=(0, 105), xlabel='Embryo')
for a in ax: a.spines[['top', 'right']].set_visible(False); a.tick_params(labelsize=7.5)
fig.tight_layout(); [fig.savefig(f'{OUT}/GT_AP_and_bleached.{x}', dpi=600 if x == 'png' else None) for x in ('png', 'pdf', 'svg')]; plt.close(fig)
print('figs ok')
