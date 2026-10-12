"""Make the data figures of the paper: Figs 10-14 and S1-S2 Figs (PLOS ONE format).

Usage (from this folder):
    python make_figures.py --results <folder with "czi mat results"> --gt <ground-truth data folder> [--out figures] [--only Fig10 Fig11 ...]

  --results  the unzipped result zip (mat data): the folder that contains "czi mat results"   (Figs 10-12)
  --gt       the ground-truth data: gt_embryos.csv, stacks/, pipeline_masks/                    (Fig 13)
Figs 14, S1 and S2 use only the tables and frames in figure_data/.
Each figure is written as <name>.tif (RGB, 300 dpi, LZW, final width) and <name>_editable.pdf (vector, live text).
The numbers plotted in each figure are written to figure_data/ as CSV files."""
import argparse, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
DATA = os.path.join(HERE, 'figure_data')
ALL = ['Fig10', 'Fig11', 'Fig12', 'Fig13', 'Fig14', 'S1_Fig', 'S2_Fig']


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results'); ap.add_argument('--gt'); ap.add_argument('--out', default=os.path.join(HERE, 'figures'))
    ap.add_argument('--only', nargs='*', default=ALL, choices=ALL)
    ap.add_argument('--tiff', choices=['pil', 'png'], default='pil',
                    help="'png': write PNG instead of TIFF (for Pillow builds whose LZW TIFF writer crashes); convert with png_to_tiff.m")
    a = ap.parse_args(); os.environ['FIG_TIFF_WRITER'] = a.tiff; os.makedirs(a.out, exist_ok=True); todo = set(a.only)
    def run(name, f, *args):
        t0 = time.time(); w, h = f(*args); print(f'{name}: {w} x {h} px ({w / 300:.2f} x {h / 300:.2f} in), {time.time() - t0:.0f} s')
    if 'Fig10' in todo:
        import fig10_intensities as m; run('Fig10', m.make, a.results, a.out, DATA)
    if todo & {'Fig11', 'Fig12'}:
        import fig11_fig12_rates as m; t0 = time.time(); w11, w12 = m.make(a.results, a.out, DATA)
        print(f'Fig11: {w11[0]} x {w11[1]} px; Fig12: {w12[0]} x {w12[1]} px, {time.time() - t0:.0f} s')
    if 'Fig13' in todo:
        import fig13_examples_vs_ground_truth as m; run('Fig13', m.make, a.gt, a.out, DATA)
    if 'Fig14' in todo:
        import fig14_comparison_quantitative as m; run('Fig14', m.make, DATA, a.out)
    if 'S1_Fig' in todo:
        import s1_fig_iou_histograms as m; run('S1_Fig', m.make, DATA, a.out)
    if 'S2_Fig' in todo:
        import s2_fig_examples_vs_pipeline as m; run('S2_Fig', m.make, os.path.join(DATA, 'S2_frames'), a.out, DATA)


if __name__ == '__main__':
    main()
