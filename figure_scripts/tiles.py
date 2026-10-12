"""Drawing helpers for the example figures (Fig 13, S2 Fig): 12 x 12 um crops shown 4x enlarged, with
outlines drawn on the inside (ground truth, yellow) or the outside (method, cyan) of each boundary."""
import numpy as np
from scipy import ndimage as ndi
from matplotlib.patches import FancyArrowPatch
from PIL import Image

YEL = np.array([1.0, 0.9, 0.0]); CYAN = np.array([0.0, 1.0, 1.0]); U = 4; R = 7; H = 96   # crop 192 px = 12 um
PX = 0.062392578125
def bandL(L, r, inner):
    st = ndi.iterate_structure(ndi.generate_binary_structure(2, 1), r); out = np.zeros(L.shape, bool)
    for sl, lab in zip(ndi.find_objects(L), range(1, int(L.max()) + 1)):
        if sl is None: continue
        sl = tuple(slice(max(s.start - r - 1, 0), s.stop + r + 1) for s in sl); m = L[sl] == lab
        out[sl] |= (m & ~ndi.binary_erosion(m, st)) if inner else (ndi.binary_dilation(m, st) & ~m)
    return out
def crop_box(cy, cx):
    cy = int(np.clip(round(cy), H, 512 - H)); cx = int(np.clip(round(cx), H, 512 - H)); return cy - H, cx - H
def tile(img, lims, green, y0, x0, ref=None, seg=None):
    sl = (slice(y0, y0 + 2 * H), slice(x0, x0 + 2 * H)); su = (slice(y0 * U, (y0 + 2 * H) * U), slice(x0 * U, (x0 + 2 * H) * U))
    g = np.clip((img[sl] - lims[0]) / (lims[1] - lims[0]), 0, 1)
    g = np.array(Image.fromarray((g * 255).astype(np.uint8)).resize((2 * H * U, 2 * H * U), Image.BILINEAR)) / 255.
    rgb = np.zeros(g.shape + (3,)); rgb[..., 1 if green else 0] = g
    if seg is not None: rgb[bandL(np.kron(seg, np.ones((U, U), int)), R, False)[su]] = CYAN
    if ref is not None: rgb[bandL(np.kron(ref, np.ones((U, U), int)), R, True)[su]] = YEL
    return rgb
def draw(fig, rect, rgb, label, arrow=None, bar=False):
    a = fig.add_axes(rect); a.imshow(rgb, interpolation='none'); a.set_xticks([]); a.set_yticks([])   # 'none': tile embedded at native 768 px
    for s in a.spines.values(): s.set_visible(False)
    a.set_title(label, fontsize=7.5, pad=2)
    if arrow is not None:
        by, bx = arrow; a.add_patch(FancyArrowPatch((bx - 230, by + 230), (bx - 80, by + 80), arrowstyle='-|>', mutation_scale=8, lw=1.3, color='white'))
    if bar:
        Lb = 5 / PX * U; yb = rgb.shape[0] - 55; xb = 45
        a.plot([xb, xb + Lb], [yb, yb], color='white', lw=1.5, solid_capstyle='butt')
        for xx in (xb, xb + Lb): a.plot([xx, xx], [yb - 18, yb + 18], color='white', lw=1.1)
        a.text(xb + Lb / 2, yb - 32, '5 µm', color='white', ha='center', va='bottom', fontsize=7.4)
    return a
def layout(fig_w, n_ex_row, n_tiles, n_rows, top=0.12, hdr=0.30, gap=0.20, lm=0.05):
    tw = (fig_w - 2 * lm - gap * (n_ex_row - 1) - 0.03 * n_ex_row * (n_tiles - 1)) / (n_ex_row * n_tiles)
    fig_h = top + n_rows * (hdr + tw) + 0.12
    return tw, fig_h
METH_SHORT = {'CellProfiler': 'CellProfiler', 'Cellpose 3': 'Cellpose 3', 'StarDist': 'StarDist'}
