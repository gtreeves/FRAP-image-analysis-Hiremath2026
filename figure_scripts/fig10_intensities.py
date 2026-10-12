"""Fig 10: Dl-GFP intensities of one ventral-side and one dorsal-side movie.
(a, c) bleached nucleus and compartment, with the mean and 95% band of the 6 nearest neighbours;
(b, d) bleached nucleus and compartment after normalization by the neighbours.
Input: the background-corrected result files (result zip, mat data)."""
import os, numpy as np, pandas as pd
import matplotlib.pyplot as plt
from plos_style import setup, save
from mat_io import load_soln

MOVIES = {'ventral': '2022-11-21/embryo1_spot 1', 'dorsal': '2022-11-24/embryo1_spot1'}


def fill_movmedian(A, w=10):
    """MATLAB fillmissing(A, 'movmedian', w), column by column."""
    A = np.array(A, float)
    for j in range(A.shape[1]):
        c = A[:, j]
        for i in np.where(np.isnan(c))[0]:
            v = c[max(0, i - w // 2):min(len(c), i + w // 2)]; v = v[~np.isnan(v)]
            if len(v): A[i, j] = np.median(v)
    return A


def load(results_dir, key):
    d = load_soln(os.path.join(results_dir, 'czi mat results', key, key.split('/')[1] + '_results.mat'))
    t = np.array(d.t, float); locs = np.atleast_1d(d.locs) - 1
    N = fill_movmedian(np.array(d.y_nucall, float)[:, locs]); C = fill_movmedian(np.array(d.y_cytall, float)[:, locs])
    yn = pd.Series(np.array(d.y_nuc, float)).interpolate(limit_direction='both').values
    yc = pd.Series(np.array(d.y_cyt, float)).interpolate(limit_direction='both').values
    return t, N, C, yn, yc


def panels(a, b, t, N, C, yn, yc, show_legend):
    mn, sn = N.mean(1), N.std(1, ddof=1); mc, sc = C.mean(1), C.std(1, ddof=1)
    ynn = yn[1:] / mn[1:] * mn[0]; ycn = yc[1:] / mc[1:] * mc[0]
    a.fill_between(t, mn - 1.96 * sn, mn + 1.96 * sn, color=(0.98, 0.85, 0.85), lw=0, label='95% conf nuc')
    a.fill_between(t, mc - 1.96 * sc, mc + 1.96 * sc, color=(0.75, 0.90, 0.78), lw=0, label='95% conf cyt')
    a.plot(t, mn, color=(0.93, 0.13, 0.13), lw=1.35, label='mean nuc'); a.plot(t, yn, color=(0.75, 0.30, 0.70), lw=1.35, label='bleach nuc')
    a.plot(t, yc, color=(0.20, 0.30, 0.65), lw=1.35, label='bleach cyt'); a.plot(t, mc, color=(0.05, 0.50, 0.20), lw=1.35, label='mean cyt')
    b.plot(t[1:], ynn, color=(0.93, 0.13, 0.13), lw=1.35, label='corrected nuc int.'); b.plot(t[1:], ycn, color=(0.20, 0.75, 0.95), lw=1.35, label='corrected cyt int.')
    if show_legend:
        h, l = a.get_legend_handles_labels(); o = [2, 3, 0, 4, 5, 1]
        a.legend([h[i] for i in o], [l[i] for i in o], ncol=2, loc='center', bbox_to_anchor=(0.5, 0.47), frameon=False, fontsize=8,
                 handlelength=1.2, columnspacing=0.6, labelspacing=0.15, borderaxespad=0.3)
        b.legend(loc='center', bbox_to_anchor=(0.5, 0.5), frameon=False, fontsize=8, handlelength=1.8)
    top = max(np.nanmax(mn + 1.96 * sn), np.nanmax(mc + 1.96 * sc), np.nanmax(ynn), np.nanmax(ycn), np.nanmax(yn))
    step = 2000 if top > 4000 else 500; ymax = np.ceil(top * 1.04 / step) * step
    for ax in (a, b):
        ax.set_xlim(-30, 1800); ax.set_ylim(0, ymax); ax.set_xticks([0, 500, 1000, 1500]); ax.set_yticks(np.arange(0, ymax + 1, step))
        ax.tick_params(labelsize=8.5, width=1.05, length=3, direction='in', top=True, right=True)
        for s in ax.spines.values(): s.set_linewidth(1.2)
        ax.set_xlabel('time (in sec)', fontsize=9.5, fontweight='bold'); ax.set_ylabel('Intensity (AU)', fontsize=9.5, fontweight='bold')
    return pd.DataFrame({'t_s': t, 'neighbours_nuc_mean': mn, 'neighbours_nuc_sd': sn, 'neighbours_cyt_mean': mc, 'neighbours_cyt_sd': sc,
                         'bleached_nuc': yn, 'bleached_cyt': yc,
                         'bleached_nuc_normalized': np.r_[np.nan, ynn], 'bleached_cyt_normalized': np.r_[np.nan, ycn]})


def make(results_dir, out_dir, data_dir):
    setup()
    fig = plt.figure(figsize=(7.5, 3.525)); tables = []
    for r, side in enumerate(['ventral', 'dorsal']):
        y0 = 0.58 if r == 0 else 0.125
        a = fig.add_axes([0.075, y0, 0.405, 0.335]); b = fig.add_axes([0.585, y0, 0.395, 0.335])
        df = panels(a, b, *load(results_dir, MOVIES[side]), show_legend=(side == 'dorsal'))
        df.insert(0, 'movie', MOVIES[side]); df.insert(0, 'side', side); tables.append(df)
        for ax in (a, b):
            ax.set_title(side + ' side', fontsize=9, fontweight='bold', pad=3)
            if r == 0: ax.set_xlabel('')
    for x, y, lab in [(0.005, 0.945, 'a'), (0.515, 0.945, 'b'), (0.005, 0.49, 'c'), (0.515, 0.49, 'd')]:
        fig.text(x, y, lab, fontsize=11)
    pd.concat(tables).to_csv(os.path.join(data_dir, 'Fig10_data.csv'), index=False)
    return save(fig, out_dir, 'Fig10', 7.5)
