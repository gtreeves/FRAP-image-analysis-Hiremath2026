"""S1 Fig: distribution of the IoU between each matched nucleus and its annotated outline, for each method.
Input: figure_data/GT_matched_nuclei.csv (written by the ground-truth package)."""
import os, numpy as np, pandas as pd, matplotlib
import matplotlib.pyplot as plt
from plos_style import setup, save
from fig14_comparison_quantitative import STYLE


def make(data_dir, out_dir):
    setup(**STYLE)
    M = ['FRAP pipeline', 'CellProfiler', 'Cellpose 3', 'StarDist']; LAB = ['Our pipeline', 'CellProfiler', 'Cellpose 3', 'StarDist']
    COL = ['#d81b8c', '#2ca02c', '#1f77b4', '#ff7f0e']
    MN = pd.read_csv(os.path.join(data_dir, 'GT_matched_nuclei.csv'))
    bins = np.arange(0.5, 1.0001, 0.01)
    fig, ax = plt.subplots(1, 4, figsize=(6.9, 1.95), sharey=True, gridspec_kw=dict(wspace=0.12, left=0.075, right=0.995, top=0.86, bottom=0.22))
    for a, m, l, c, L in zip(ax, M, LAB, COL, 'abcd'):
        v = MN[MN.method == m].IoU.to_numpy(); med = np.median(v)
        a.hist(v, bins=bins, color=matplotlib.colors.to_rgba(c, 0.6), edgecolor=c, linewidth=0.5)
        a.axvline(med, color='k', ls='--', lw=1.0)
        a.text(0.04, 0.96, f'{l}\nn = {len(v):,}\nmedian = {med:.2f}', transform=a.transAxes, va='top', fontsize=7.4)
        a.set_xlim(0.5, 1.0); a.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0]); a.set_xticklabels(['0.5', '', '0.7', '', '0.9', ''])
        a.set_xlabel('IoU with ground truth'); a.text(-0.08 if L != 'a' else -0.30, 1.12, L, transform=a.transAxes, fontsize=10, va='top')
    ax[0].set_ylabel('Number of nuclei'); ax[0].set_ylim(0, 1450)
    return save(fig, out_dir, 'S1_Fig', 7.5)
