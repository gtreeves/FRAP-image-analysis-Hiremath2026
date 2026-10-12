"""Shared settings for the paper figures: Arial text, PLOS ONE output (300 dpi TIFF at the final width)."""
import os, io, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
import matplotlib._mathtext as _mathtext
from PIL import Image

TIFF_WRITER = os.environ.get('FIG_TIFF_WRITER', 'pil')   # 'pil' (LZW TIFF via Pillow) or 'png' (see png_to_tiff.m)
FONT_DIRS = [os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fonts'), r'C:\Windows\Fonts', '/Library/Fonts']


def setup(**rc):
    """Register Arial (from a 'fonts' folder next to this file or the system fonts) and set the common style."""
    for d in FONT_DIRS:
        for f in ('arial.ttf', 'arialbd.ttf', 'ariali.ttf', 'arialbi.ttf', 'Arial.ttf', 'Arial Bold.ttf'):
            p = os.path.join(d, f)
            if os.path.exists(p):
                fm.fontManager.addfont(p)
    if 'Arial' not in {f.name for f in fm.fontManager.ttflist}:
        print('WARNING: Arial was not found; text will use another font. Copy arial*.ttf into figure_scripts/fonts.')
    _mathtext.SHRINK_FACTOR = 0.8          # sub- and superscripts at 80% so they stay >= 8 pt
    plt.rcParams.update({'font.family': 'Arial', 'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial',
                         'mathtext.it': 'Arial:italic', 'mathtext.bf': 'Arial:bold',
                         'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none'})
    plt.rcParams.update(rc)


def ylabel_superscript(fig, ax, before, sup, after, fontsize=10, fontweight='bold', raise_frac=0.40, sup_scale=0.75):
    """Y-axis label with a superscript, e.g. ylabel_superscript(fig, ax, 'rate constant (s', '\u22121', ')').
    Drawn from plain text pieces, because some matplotlib versions place mathtext superscripts too high with Arial."""
    ax.set_ylabel(before + sup + after, fontsize=fontsize, fontweight=fontweight)   # reserve the space and get the position
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    bb = ax.yaxis.label.get_window_extent(r); ax.set_ylabel('')
    px = fig.dpi / 72.0
    x_base = bb.x1 - 0.22 * fontsize * px          # baseline of the rotated text (glyphs extend to the left of it)
    kw = dict(rotation=90, rotation_mode='anchor', ha='left', va='baseline', fontweight=fontweight)
    parts = [(before, fontsize, 0.0), (sup, fontsize * sup_scale, raise_frac * fontsize), (after, fontsize, 0.0)]
    texts = [fig.text(0, 0, s, fontsize=fs, **kw) for s, fs, _ in parts]
    fig.canvas.draw()
    hs = [t_.get_window_extent(r).height for t_ in texts]
    y = (bb.y0 + bb.y1) / 2 - sum(hs) / 2; inv = fig.transFigure.inverted()
    for t_, (s, fs, rise), h in zip(texts, parts, hs):
        t_.set_position(inv.transform((x_base - rise * px, y))); y += h
    return texts


def save(fig, out_dir, name, width_in, dpi=300, margin_in=0.04):
    """Write <name>_editable.pdf (vector, live text) and <name>.tif (RGB, LZW, `dpi`, exactly `width_in` wide).

    The figure is cropped to its artwork with a small margin and scaled to the final width, as required by
    PLOS ONE (2.63-7.5 in wide, 300-600 dpi)."""
    os.makedirs(out_dir, exist_ok=True)
    fig_w = fig.get_size_inches()[0]
    pad = margin_in * fig_w / width_in
    fig.savefig(os.path.join(out_dir, name + '_editable.pdf'), bbox_inches='tight', pad_inches=pad)
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=4 * dpi * width_in / fig_w, bbox_inches='tight', pad_inches=pad)
    buf.seek(0)
    im = Image.open(buf).convert('RGB')
    W = int(round(width_in * dpi)); H = int(round(im.size[1] * W / im.size[0]))
    im = im.resize((W, H), Image.LANCZOS)                     # supersampled, then downsampled: smooth lines and tiles
    if TIFF_WRITER == 'png':
        # some Pillow builds crash when writing LZW TIFF; then write a lossless PNG at the final size and convert it
        # with png_to_tiff.m (MATLAB imwrite, LZW, 300 dpi)
        im.save(os.path.join(out_dir, name + '.png'), dpi=(dpi, dpi))
    else:
        im.save(os.path.join(out_dir, name + '.tif'), compression='tiff_lzw', dpi=(dpi, dpi))
    plt.close(fig)
    return W, H
