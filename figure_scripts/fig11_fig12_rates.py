"""Figs 11 and 12: import and export rate constants (Fig 11) and the simulated recovery with the N/C ratio (Fig 12),
13 movies per side. Input: the background-corrected result files (result zip, mat data); the rate constants are the
fitted values stored in each file (soln.k_in, soln.k_out)."""
import os, glob, warnings, numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from plos_style import setup, save, ylabel_superscript
from mat_io import load_soln


def smooth(y, span=20):
    """MATLAB smooth(y, 20): moving average over 19 points, shrinking at the ends."""
    y = np.asarray(y, float); n = len(y); h = (span - 1 - (span % 2 == 0)) // 2; o = np.empty(n)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)       # all-NaN windows give NaN, as in MATLAB
        for i in range(n): k = min(i, n - 1 - i, h); o[i] = np.nanmean(y[i - k:i + k + 1])
    return o


def fill_movmedian(A, w=10):
    """MATLAB fillmissing(A, 'movmedian', w), column by column."""
    A = np.array(A, float)
    for j in range(A.shape[1]):
        c = A[:, j]
        for i in np.where(np.isnan(c))[0]:
            v = c[max(0, i - w // 2):min(len(c), i + w // 2)]; v = v[~np.isnan(v)]
            if len(v): A[i, j] = np.median(v)
    return A


def load_all(results_dir):
    data = {'ventral': [], 'dorsal': []}; rows = []
    root = os.path.join(results_dir, 'czi mat results')
    for f in sorted(glob.glob(os.path.join(root, '*', '*', '*_results.mat'))):
        d = load_soln(f)
        key = os.path.relpath(os.path.dirname(f), root).replace(os.sep, '/')
        side = str(d.side)
        N = np.array(d.y_nucall, float); C = np.array(d.y_cytall, float)
        locs = np.atleast_1d(d.locs) - 1; ib = int(d.idx_bleach) - 1
        mc = fill_movmedian(C[:, locs]).mean(1)
        t = np.array(d.t, float)[1:]
        ycr = smooth(C[1:, ib] / mc[1:]); ynr = smooth(N[1:, ib] / mc[1:])
        data[side].append(dict(key=key, t=t, ycr=ycr, ynr=ynr, kin=float(d.k_in), kout=float(d.k_out)))
        dl = np.atleast_1d(d.delta2)
        rows.append(dict(movie=key, side=side, k_in=float(d.k_in), k_out=float(d.k_out),
                         k_in_CI95_halfwidth=float(dl[1]), k_out_CI95_halfwidth=float(dl[2]), r2=float(d.r2)))
    return data, pd.DataFrame(rows)


def make(results_dir, out_dir, data_dir):
    setup(**{'axes.linewidth': 1.0})
    data, rates = load_all(results_dir)
    rates.to_csv(os.path.join(data_dir, 'Fig11_data.csv'), index=False)
    curves, ncr_rows = [], []
    RED=(0.85,0.13,0.10); BLUE=(0.0,0.447,0.741); GREY=(0.30,0.30,0.30)
    # ---------------- Fig. 11 ----------------
    fig, axs = plt.subplots(1, 2, figsize=(3.5, 2.45))
    style = {'ventral': [((0.87,0.66,0.80), RED, RED, 2.2), ((1.00,0.78,0.50), BLUE, (0.93,0.55,0.15), 2.2)],
             'dorsal':  [((0.92,0.56,0.45), RED, RED, 1.0), ((0.50,0.76,0.58), BLUE, (0.05,0.55,0.50), 1.0)]}
    for ax, side in zip(axs, ['ventral', 'dorsal']):
        vals = [np.array([m['kin'] for m in data[side]]), np.array([m['kout'] for m in data[side]])]
        for b, (v, (fc, ec, mc_, ww)) in enumerate(zip(vals, style[side])):
            ax.boxplot(v, positions=[b+1], widths=0.5, patch_artist=True, whis=1.5,
                       boxprops=dict(facecolor=fc, edgecolor=ec, linewidth=1.3),
                       medianprops=dict(color=mc_, linewidth=1.3),
                       whiskerprops=dict(color=(0.30,0.30,0.30), linewidth=ww*0.73), capprops=dict(color=(0.30,0.30,0.30), linewidth=ww*0.73),
                       capwidths=0.22,
                       flierprops=dict(marker='o', markerfacecolor='none', markeredgecolor=RED, markersize=3.6, markeredgewidth=0.8))
        top = max(v.max() for v in vals)
        step = 0.01 if side == 'ventral' else 0.002
        ymax = np.ceil(top / step) * step
        ax.set_ylim(0, ymax); ax.set_xlim(0.35, 2.65)
        ticks = np.arange(0, ymax + step/2, step); ax.set_yticks(ticks)
        ax.set_yticklabels(['0.0'] + [('%.2f' if side=='ventral' else '%.3f') % x for x in ticks[1:]])
        ax.set_yticks(np.arange(0, ymax + step/8, step/4), minor=True)
        ax.grid(axis='y', which='minor', color=(0.92,0.92,0.92), linewidth=0.6)
        ax.grid(axis='y', which='major', color=(0.86,0.86,0.86), linewidth=0.9); ax.set_axisbelow(True)
        ax.set_xticks([1, 2]); ax.set_xticklabels([r'$\mathbf{k}_{\mathbf{in}}$', r'$\mathbf{k}_{\mathbf{out}}$'], fontsize=10.5)
        for s_ in ['top', 'right']: ax.spines[s_].set_visible(False)
        for s_ in ['left', 'bottom']: ax.spines[s_].set_linewidth(1.75)
        ax.tick_params(axis='y', which='major', labelsize=9, width=1.6, length=3.6, direction='in', pad=3)
        ax.tick_params(axis='y', which='minor', length=0)
        ax.tick_params(axis='x', width=1.6, length=3.6, direction='in', pad=3)
        for fl in ax.lines: fl.set_clip_on(False)
        ax.set_title(side + ' side', fontsize=9, fontweight='bold', color=(0.15,0.15,0.15))
    fig.subplots_adjust(left=0.175, right=0.985, bottom=0.14, top=0.90, wspace=0.42)
    ylabel_superscript(fig, axs[0], 'rate constant (s', '\u22121', ')', fontsize=10)
    W11 = save(fig, out_dir, 'Fig11', 3.5)
    # ---------------- Fig. 12 ----------------
    fig = plt.figure(figsize=(7.5, 2.55))
    axA = fig.add_axes([0.065, 0.20, 0.25, 0.64]); axB = fig.add_axes([0.40, 0.20, 0.25, 0.64])
    frame = fig.add_axes([0.715, 0.12, 0.27, 0.78]); frame.set_xticks([]); frame.set_yticks([])
    for s in frame.spines.values(): s.set_linewidth(0.7)
    frame.patch.set_alpha(0)
    axC = fig.add_axes([0.79, 0.27, 0.17, 0.56])
    ncr = {}
    summary = []
    for ax, side in [(axA, 'ventral'), (axB, 'dorsal')]:
        ms = data[side]; L = min(len(m['t']) for m in ms); tg = ms[0]['t'][:L]
        yc = np.nanmean([m['ycr'][:L] for m in ms], 0); y0 = np.nanmean([m['ynr'][0] for m in ms])
        Y = []
        for m in ms:
            sol = solve_ivp(lambda tq, yq: m['kin'] * np.interp(tq, tg, yc) - m['kout'] * yq, (tg[0], tg[-1]), [y0],
                            t_eval=tg, rtol=1e-8, atol=1e-10)
            Y.append(sol.y[0])
        Y = np.array(Y); ncr[side] = Y[:, -1] / yc[-1]
        mean, sd = Y.mean(0), Y.std(0, ddof=1); lo, hi = mean - 1.96*sd, mean + 1.96*sd
        curves.append(pd.DataFrame({'side': side, 't_s': tg, 'simulated_nuc_mean': mean, 'simulated_nuc_95_low': lo, 'simulated_nuc_95_high': hi, 'normalized_cyt_mean': yc}))
        ncr_rows += [dict(side=side, movie=m['key'], NC_ratio_end=v) for m, v in zip(ms, ncr[side])]
        ax.fill_between(tg, lo, hi, color=(0.97,0.78,0.78), linewidth=0, label='95% conf nuc')
        ax.plot(tg, mean, color=(0.93,0.13,0.13), linewidth=1.4, label='mean simulated nuc')
        ax.plot(tg, yc, color=(0.22,0.33,0.65), linewidth=1.4, label='norm. cyt')
        h, l = ax.get_legend_handles_labels(); order = [1, 0, 2]
        ax.legend([h[i] for i in order], [l[i] for i in order], loc='upper left', bbox_to_anchor=(0.0, 1.0 if side == 'ventral' else 0.80), frameon=False, fontsize=8,
                  handlelength=1.8, borderaxespad=0.2, labelspacing=0.25)
        ymax = 5 if side == 'ventral' else 1
        while hi.max() > ymax * 0.98 or yc.max() > ymax * 0.98: ymax += 1 if side == 'ventral' else 0.2
        ax.set_xlim(0, 2000); ax.set_ylim(0, ymax); ax.set_xticks([0, 500, 1000, 1500, 2000])
        if side == 'dorsal': ax.set_yticks(np.round(np.arange(0, ymax + 0.01, 0.2), 1)); ax.set_yticklabels(['0'] + ['%.1f' % x for x in np.arange(0.2, ymax + 0.01, 0.2)])
        ax.tick_params(labelsize=8, width=1.05, length=2.5, direction='in', top=True, right=True)
        for s in ax.spines.values(): s.set_linewidth(1.1)
        ax.set_xlabel('time (in sec)', fontsize=8.5, fontweight='bold'); ax.set_ylabel('Relative intensity', fontsize=8.5, fontweight='bold')
        ax.set_title(side + ' side', fontsize=8.5, fontweight='bold', pad=4)
        summary.append((side, ncr[side].min(), ncr[side].max(), hi.max()))
    for i, (side, fc, ec, w) in enumerate([('ventral', (0.97,0.80,0.80), RED, 1.2), ('dorsal', (0.75,0.88,0.98), (0.12,0.56,0.95), 0.8)]):
        axC.boxplot(ncr[side], positions=[i+1], widths=0.45, patch_artist=True, whis=1.5,
                    boxprops=dict(facecolor=fc, edgecolor=ec, linewidth=w), medianprops=dict(color=ec, linewidth=w),
                    whiskerprops=dict(color=(0.1,0.1,0.1), linewidth=w), capprops=dict(color=(0.1,0.1,0.1), linewidth=w),
                    flierprops=dict(marker='o', markerfacecolor='none', markeredgecolor=ec, markersize=3.75))
    cmax = max(6, np.ceil(max(ncr['ventral'].max(), 1) + 0.5))
    axC.set_ylim(0, cmax); axC.set_xlim(0.4, 2.6); axC.set_yticks(range(0, int(cmax) + 1))
    axC.set_xticks([1, 2]); axC.set_xticklabels(['ventral\nside', 'dorsal\nside'], fontsize=8, fontweight='bold')
    axC.set_ylabel('N/C ratio', fontsize=8.5, fontweight='bold')
    for s in ['top', 'right']: axC.spines[s].set_visible(False)
    for s in ['left', 'bottom']: axC.spines[s].set_linewidth(1.05)
    axC.tick_params(axis='y', labelsize=8, width=1.05, length=2.5, direction='in'); axC.tick_params(axis='x', length=0)
    for x, lab in [(0.008, 'a'), (0.335, 'b'), (0.700, 'c')]: fig.text(x, 0.93, lab, fontsize=10)
    W12 = save(fig, out_dir, 'Fig12', 7.5)
    pd.concat(curves).to_csv(os.path.join(data_dir, 'Fig12ab_data.csv'), index=False)
    pd.DataFrame(ncr_rows).to_csv(os.path.join(data_dir, 'Fig12c_data.csv'), index=False)
    return W11, W12
