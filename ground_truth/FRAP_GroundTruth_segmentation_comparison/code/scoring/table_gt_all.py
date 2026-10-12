"""Build the LaTeX tables (GT_results_table.tex, _standalone.tex) from GT_summary.csv.
Usage: python table_gt_all.py OUT   (OUT holds GT_summary.csv; writes the .tex files there)"""
import sys, pandas as pd
OUT = sys.argv[1]
S = pd.read_csv(f'{OUT}/GT_summary.csv')
METHODS = list(dict.fromkeys(S['method']))
words = {1: 'one', 2: 'two', 3: 'three', 4: 'four', 5: 'five', 6: 'six', 7: 'seven', 8: 'eight', 9: 'nine', 10: 'ten'}
nall = int(S[(S.method == METHODS[0]) & (S.group == 'All')].embryos.iloc[0])

def row(r):
    f3 = lambda x: f'{x:.3f}'
    return {
        'Embryos': f'{int(r.embryos)}', 'Frames scored': f'{int(r.frames)}',
        'Nuclei, ground truth': f'{int(r.nuclei_GT)}', 'Nuclei, method': f'{int(r.nuclei_method)}',
        'Matched (TP)': f'{int(r.TP)}', 'Missed (FN)': f'{int(r.FN)}', 'Extra (FP)': f'{int(r.FP)}',
        '-1': None,
        'Precision': f3(r.precision), 'Recall': f3(r.recall), 'F1 score': f3(r.F1),
        '-2': None,
        'IoU, median [IQR]': f'{r.IoU_median:.3f} [{r.IoU_Q1:.3f}--{r.IoU_Q3:.3f}]',
        'IoU, mean $\\pm$ SD': f'{r.IoU_mean:.3f} $\\pm$ {r.IoU_SD:.3f}',
        'Dice, mean': f3(r.Dice_mean), 'Pixel Dice, mean': f3(r.pixel_Dice_mean),
        'Area difference ($\\mu$m$^2$), mean $\\pm$ SD': f'{r.area_diff_um2_mean:.2f} $\\pm$ {r.area_diff_um2_SD:.2f}',
        'Centroid shift ($\\mu$m), mean': f3(r.centroid_shift_um_mean),
        '-3': None,
        'AP at IoU 0.50 / 0.75 / 0.90': f"{r['AP@0.5']:.2f} / {r['AP@0.75']:.2f} / {r['AP@0.9']:.2f}",
        'mAP (IoU 0.50--0.95)': f"{r['mAP(0.5:0.95)']:.3f}",
        '-4': None,
        'Bleached nucleus matched (frames)': f'{int(r.bleached_matched)}/{int(r.bleached_frames)} ({100 * r.bleached_matched / r.bleached_frames:.0f}\\%)',
        'Bleached nucleus IoU, mean': f3(r.bleached_IoU_mean),
    }

def table(groups, label, what):
    cols = [(m, g) for m in METHODS for g in groups]
    R = {c: row(S[(S.method == c[0]) & (S.group == c[1])].iloc[0]) for c in cols}
    ng = len(groups)
    cap = (f'Segmentation accuracy against manual (Labkit) ground truth for {words.get(nall, nall)} Dorsal FRAP embryos, grouped by {what}. '
           'Nuclei are matched one-to-one at IoU $\\geq$ 0.5; nuclei touching the image edge are excluded. AP = TP/(TP+FP+FN). '
           'Area difference is method minus ground truth. Our pipeline uses the tracked masks before the final erosion, '
           'with the bleached nucleus reconstructed for embryos 2 and 20.')
    L = ['\\begin{table}[htbp]', '\\centering\\footnotesize\\setlength{\\tabcolsep}{4pt}', f'\\caption{{{cap}}}', f'\\label{{{label}}}',
         '\\resizebox{\\linewidth}{!}{%', '\\begin{tabular}{l' + 'c' * len(cols) + '}', '\\toprule',
         ' & ' + ' & '.join(f'\\multicolumn{{{ng}}}{{c}}{{{m}}}' for m in METHODS) + ' \\\\',
         ''.join(f'\\cmidrule(lr){{{2 + i * ng}-{1 + (i + 1) * ng}}}' for i in range(len(METHODS))),
         'Metric & ' + ' & '.join(g for _, g in cols) + ' \\\\', '\\midrule']
    for k in R[cols[0]]:
        if k.startswith('-'): L.append('\\midrule'); continue
        L.append(k + ' & ' + ' & '.join(R[c][k] for c in cols) + ' \\\\')
    L += ['\\bottomrule', '\\end{tabular}}', '\\end{table}']
    return '\n'.join(L)

body = table(['Good', 'Bad', 'All'], 'tab:gt_quality', 'image quality') + '\n\\clearpage\n' + table(['Ventral', 'Dorsal'], 'tab:gt_side', 'side')
open(f'{OUT}/GT_results_table.tex', 'w').write('% needs \\usepackage{booktabs} and \\usepackage{graphicx}\n' + body + '\n')
pre = ('\\documentclass[11pt]{article}\n\\usepackage[margin=1.5cm,landscape]{geometry}\n\\usepackage{booktabs}\n\\usepackage{graphicx}\n'
       '\\usepackage[T1]{fontenc}\n\\usepackage{helvet}\n\\renewcommand{\\familydefault}{\\sfdefault}\n\\pagestyle{empty}\n\\begin{document}\n')
open(f'{OUT}/GT_results_table_standalone.tex', 'w').write(pre + body + '\n\\end{document}\n')
print('tables written for', METHODS)
