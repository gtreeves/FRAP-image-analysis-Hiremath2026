import logging; logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
"""Summaries, workbook and figures for frap_track.py output.  Usage: python frap_figs.py OUT"""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
OUT = sys.argv[1]
plt.rcParams.update({'font.family': ['Liberation Sans', 'Arial', 'DejaVu Sans'], 'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.linewidth': 0.8, 'svg.fonttype': 'none', 'pdf.fonttype': 42})
C = pd.read_csv(f'{OUT}/frap_curves.csv'); F = pd.read_csv(f'{OUT}/frap_fits.csv')
B = pd.read_csv(f'{OUT}/bleached_track.csv'); T = pd.read_csv(f'{OUT}/tracking_all_nuclei.csv')
PIPE = 'FRAP pipeline (L_nuc, as used)'
METH = [PIPE, 'CellProfiler', 'Cellpose 3', 'StarDist']
COL = {PIPE: '#d81b8c', 'FRAP pipeline (L_track)': '#d81b8c', 'CellProfiler': '#2ca02c', 'Cellpose 3': '#1f77b4', 'StarDist': '#ff7f0e', 'Ground truth': '#222222'}
SHORT = {PIPE: 'FRAP pipeline', 'FRAP pipeline (L_track)': 'FRAP pipeline (L_track)', 'CellProfiler': 'CellProfiler', 'Cellpose 3': 'Cellpose 3', 'StarDist': 'StarDist'}

# ---------------- FRAP summary ----------------
emb = list(dict.fromkeys(F.embryo)); nfr = B.drop_duplicates('embryo').set_index('embryo').post_bleach_frames
F['post_bleach_frames'] = F.embryo.map(nfr)
F['usable_fit'] = (F.n_fit >= 0.8 * F.post_bleach_frames) & F.k_out.gt(1e-5) & F.k_out.lt(0.05)
F['ratio_err_pct'] = 100 * (F.k_in_over_k_out_ratio_to_GT - 1)
F['kout_err_pct'] = 100 * (F.k_out_ratio_to_GT - 1)
GTU = F[F.method == 'Ground truth'].set_index('embryo').usable_fit
rows = []
main = [e for e in emb if not e.startswith('E20')]          # E20 has only 13 labelled post-bleach frames
for m in [PIPE, 'FRAP pipeline (L_track)', 'CellProfiler', 'Cellpose 3', 'StarDist']:
    d = F[F.method == m]; dm = d[d.embryo.isin(main)]; ok = dm[dm.usable_fit & dm.embryo.map(GTU)]
    rmse = d.curve_RMSE_vs_GT.dropna()
    rows.append({'method': m, 'embryos': len(d), 'usable FRAP fits (of 9, E20 excluded)': int(dm.usable_fit.sum()),
                 'curve RMSE vs GT, median': rmse.median(), 'curve RMSE vs GT, max': rmse.max(),
                 'embryos with curve RMSE > 0.05': int((rmse > 0.05).sum()),
                 'k_in/k_out error vs GT, median |%|': ok.ratio_err_pct.abs().median(),
                 'k_out error vs GT, median |%|': ok.kout_err_pct.abs().median(),
                 'k_in/k_out within 20% of GT (embryos)': f'{int((ok.ratio_err_pct.abs() <= 20).sum())}/{len(ok)}'})
FS = pd.DataFrame(rows)

# ---------------- tracking summary ----------------
rows = []
for m in [PIPE, 'FRAP pipeline (L_track)', 'CellProfiler', 'Cellpose 3', 'StarDist']:
    b = B[B.method == m]
    rows.append({'method': m, 'tracking': b.tracking.iloc[0], 'post-bleach frames': int(b.post_bleach_frames.sum()),
                 'bleached nucleus followed (frames)': int(b.frames_found.sum()),
                 'frames missing': int(b.frames_missing.sum()), 'frames on the wrong nucleus': int(b.frames_wrong_nucleus.sum()),
                 'track breaks (gap runs)': int(b.gap_runs.sum()),
                 'embryos followed without any miss or wrong frame': f"{int(((b.frames_missing == 0) & (b.frames_wrong_nucleus == 0)).sum())}/{len(b)}",
                 'bleached area CV, median': b.area_CV.median(), 'bleached area jumps > 20% (steps)': int(b.steps_area_jump_gt20pct.sum())})
gb = B[B.method == 'Ground truth']
rows.append({'method': 'Ground truth', 'tracking': 'overlap linking of the labels', 'post-bleach frames': int(gb.post_bleach_frames.sum()),
             'bleached nucleus followed (frames)': int(gb.frames_found.sum()), 'frames missing': int(gb.frames_missing.sum()),
             'bleached area CV, median': gb.area_CV.median(), 'bleached area jumps > 20% (steps)': int(gb.steps_area_jump_gt20pct.sum())})
BS = pd.DataFrame(rows)
rows = []
for m in ['FRAP pipeline (L_track)', 'CellProfiler', 'Cellpose 3', 'StarDist']:
    t = T[T.method == m]
    rows.append({'method': SHORT[m] if m != 'FRAP pipeline (L_track)' else 'FRAP pipeline', 'GT nucleus links (frame to next frame)': int(t.GT_links.sum()),
                 'found in both frames (%)': 100 * t.found_both_frames.sum() / t.GT_links.sum(),
                 'identity kept (% of found)': 100 * t.identity_kept.sum() / t.found_both_frames.sum(),
                 'followed correctly (% of all links)': 100 * t.identity_kept.sum() / t.GT_links.sum(),
                 'area change error, median of embryo medians (%)': t.median_abs_area_change_error_pct.median(),
                 'area jumps > 20% (% of steps, mean over embryos)': t.pct_area_jumps_gt20.mean(),
                 'GT area jumps > 20% (%)': t.GT_pct_area_jumps_gt20.mean()})
TS = pd.DataFrame(rows)

with pd.ExcelWriter(f'{OUT}/FRAP_and_tracking_results.xlsx') as xw:
    FS.to_excel(xw, sheet_name='FRAP_summary', index=False)
    F.to_excel(xw, sheet_name='FRAP_fits_per_embryo', index=False)
    C.to_excel(xw, sheet_name='FRAP_curves', index=False)
    BS.to_excel(xw, sheet_name='Bleached_track_summary', index=False)
    B.to_excel(xw, sheet_name='Bleached_track_per_embryo', index=False)
    TS.to_excel(xw, sheet_name='Tracking_summary', index=False)
    T.to_excel(xw, sheet_name='Tracking_per_embryo', index=False)
    pd.DataFrame({'item': ['curve', 'rel_intensity', 'fit', 'k_in, k_out', 'k_in_over_k_out', 'curve_RMSE_vs_GT', 'usable_fit',
                           'bleached track', 'identity kept', 'area change error'],
                  'meaning': ['bleached-nucleus Dorsal intensity, normalised by the 6 nearest nuclei',
                              '(nucleus/neighbours) at time t divided by the same at frame 1 (pre-bleach)',
                              'ftn_FRAP model from frame 2: dy/dt = k_in*ycyt - k_out*y, smoothed as in fit_bleach_dl.m',
                              'import and export rate constants (1/s)', 'steady-state nuclear/cytoplasmic ratio of the model',
                              'root-mean-square difference of rel_intensity from the ground-truth curve (post-bleach frames)',
                              'at least 80% of post-bleach frames available and 1e-5 < k_out < 0.05 1/s',
                              'pipeline: its own tracked label; other methods: simple overlap linking from the frame-1 bleached nucleus',
                              'same pipeline label, or linked one-to-one at IoU >= 0.5 (other methods)',
                              '|relative area change of the method| minus |that of the ground truth| between consecutive frames']}).to_excel(xw, sheet_name='Definitions', index=False)
FS.round(3).to_csv(f'{OUT}/FRAP_summary.csv', index=False); BS.round(3).to_csv(f'{OUT}/bleached_track_summary.csv', index=False); TS.round(3).to_csv(f'{OUT}/tracking_summary.csv', index=False)

def save(fig, name):
    for ext in ('png', 'pdf', 'svg'): fig.savefig(f'{OUT}/{name}.{ext}', dpi=600 if ext == 'png' else None, bbox_inches='tight')
    plt.close(fig)

# ---------------- Figure 1: FRAP curves per embryo ----------------
fig, axs = plt.subplots(2, 5, figsize=(10, 4.4), sharex=False)
for ax, e in zip(axs.ravel(), emb):
    for m in ['Ground truth'] + METH:
        d = C[(C.embryo == e) & (C.method == m)].iloc[1:]
        ax.plot(d.t_s / 60, d.rel_intensity, color=COL[m], lw=2 if m != 'Ground truth' else 1.6, ls='--' if m == 'Ground truth' else '-',
                label='Ground truth' if m == 'Ground truth' else SHORT[m], zorder=5 if m == 'Ground truth' else 2)
    row = F[(F.embryo == e)].iloc[0]
    ax.set_title(f'{e[:3]} ({row.case}, {row.side})', fontsize=8)
    ax.set_ylim(0, 1.8 if C[(C.embryo == e)].rel_intensity.max() > 1.3 else 1.3)
    ax.grid(axis='y', color='#e6e6e6', lw=0.6); ax.set_axisbelow(True)
for ax in axs[1]: ax.set_xlabel('Time after bleach (min)')
for ax in axs[:, 0]: ax.set_ylabel('Bleached nucleus intensity\n(fraction of pre-bleach)')
h, l = axs[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc='upper center', ncol=5, frameon=False, bbox_to_anchor=(0.5, 1.03))
fig.tight_layout(rect=(0, 0, 1, 0.95)); save(fig, 'FRAP_curves_per_embryo')

# ---------------- Figure 2: FRAP accuracy ----------------
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
rng = np.random.default_rng(1)
for i, m in enumerate(METH):
    v = F[(F.method == m)].curve_RMSE_vs_GT.dropna().values
    ax[0].scatter(i + rng.uniform(-0.12, 0.12, len(v)), v, s=22, color=COL[m], edgecolor='white', lw=0.6, zorder=3)
    ax[0].hlines(np.median(v), i - 0.25, i + 0.25, color='#222222', lw=1.5, zorder=4)
    ok = F[(F.method == m) & F.embryo.isin(main) & F.usable_fit & F.embryo.map(GTU)]
    w = ok.k_in_over_k_out_ratio_to_GT.values
    ax[1].scatter(i + rng.uniform(-0.12, 0.12, len(w)), w, s=22, color=COL[m], edgecolor='white', lw=0.6, zorder=3)
    if len(w): ax[1].hlines(np.median(w), i - 0.25, i + 0.25, color='#222222', lw=1.5, zorder=4)
    nf = 9 - len(ok)
    if nf: ax[1].text(i, 0.36, f'{nf} failed', ha='center', fontsize=7, color='#555555')
for a in ax: a.set_xticks(range(len(METH))); a.set_xticklabels([SHORT[m] for m in METH], rotation=20, ha='right'); a.grid(axis='y', color='#e6e6e6', lw=0.6); a.set_axisbelow(True)
ax[0].set_ylabel('Curve error vs ground truth (RMSE)'); ax[0].set_title('A  Recovery curve', loc='left', fontsize=9)
ax[1].set_yscale('log'); ax[1].axhline(1, color='#999999', lw=0.8, ls='--')
ax[1].set_yticks([0.5, 1, 2, 3]); ax[1].set_yticklabels(['0.5', '1', '2', '3']); ax[1].yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax[1].set_ylim(0.33, 3.2)
ax[1].set_ylabel('k$_{in}$/k$_{out}$ relative to ground truth'); ax[1].set_title('B  Fitted k$_{in}$/k$_{out}$ (E20 excluded)', loc='left', fontsize=9)
fig.tight_layout(); save(fig, 'FRAP_accuracy')

# ---------------- Figure 3: tracking ----------------
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
x = np.arange(len(METH))
fol = [B[B.method == m].frames_found.sum() - np.nan_to_num(B[B.method == m].frames_wrong_nucleus).sum() for m in METH]
tot = B[B.method == PIPE].post_bleach_frames.sum()
bars = ax[0].bar(x, [100 * f / tot for f in fol], color=[COL[m] for m in METH], width=0.6, edgecolor='white', lw=2)
for xi, f in zip(x, fol): ax[0].text(xi, 100 * f / tot + 1.5, f'{int(f)}/{tot}', ha='center', fontsize=7)
ax[0].set_ylabel('Post-bleach frames on the\ncorrect bleached nucleus (%)'); ax[0].set_ylim(0, 112)
ax[0].set_title('A  Bleached nucleus, own tracking', loc='left', fontsize=9)
tm = ['FRAP pipeline (L_track)', 'CellProfiler', 'Cellpose 3', 'StarDist']
vals = [100 - 100 * T[T.method == m].identity_kept.sum() / T[T.method == m].GT_links.sum() for m in tm]
ax[1].bar(x, vals, color=[COL[m] for m in tm], width=0.6, edgecolor='white', lw=2)
for xi, v in zip(x, vals): ax[1].text(xi, v + 0.05, f'{v:.1f}%', ha='center', fontsize=7)
ax[1].set_ylim(0, 3.3); ax[1].set_ylabel('Nucleus links lost or broken (%)')
ax[1].set_title('B  All nuclei, frame to frame', loc='left', fontsize=9)
for a in ax: a.set_xticks(x); a.set_xticklabels([SHORT[m] for m in METH], rotation=20, ha='right'); a.grid(axis='y', color='#e6e6e6', lw=0.6); a.set_axisbelow(True)
fig.tight_layout(); save(fig, 'Tracking_consistency')
print(FS.round(3).to_string()); print(BS.round(3).to_string()); print(TS.round(2).to_string())
