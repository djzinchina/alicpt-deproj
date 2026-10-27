"""Small helpers used by the example notebooks of the AliCPT beam-systematics
and deprojection dataset.

This module is released under the GNU General Public License v3.0; see LICENSE.
"""
import numpy as np
import healpy as hp
import matplotlib as mpl
import matplotlib.pyplot as plt


# ------------------------------------------------------------------ colours --
def _planck_cmap(tone, rgb, name):
    """Build a 256-entry matplotlib colormap from PSM/Planck colour anchors."""
    x = np.arange(256)
    return mpl.colors.ListedColormap(
        np.column_stack([np.interp(x, tone, c) for c in rgb]) / 255., name=name)


# Diverging map for a Gaussian signal.
GAUSSIAN = _planck_cmap(
    np.array([0, 42, 85, 127, 170, 212, 255]),
    np.array([[0, 0, 0, 255, 255, 255, 100],
              [0, 112, 221, 237, 180, 75, 0],
              [255, 255, 255, 217, 0, 0, 0]]), 'gaussian')

# High-contrast map for foreground-dominated fields.
FOREGROUND = _planck_cmap(
    np.array([0, 13, 26, 39, 52, 65, 76, 77, 88, 101, 114, 127, 140, 153, 166,
              179, 192, 205, 218, 231, 255]),
    np.array([[0, 10, 30, 80, 191, 228, 241, 241, 245, 248, 249.9, 242.25, 204,
               165, 114, 127.5, 178.5, 204, 229.5, 242.25, 252.45],
              [0, 20, 184, 235, 239, 240, 241, 241, 240, 235, 204, 153, 76.5, 32,
               0, 127.5, 178.5, 204, 229.5, 242.25, 252.45],
              [255, 255, 255, 255, 250, 245, 212, 212, 175, 130, 38.25, 12.75, 0,
               32, 32, 153, 204, 229.5, 242.25, 249.9, 255]]), 'foreground')

# healpy applies its ``bgcolor``/``badcolor`` arguments only when the colormap
# is passed *by name*; for a Colormap object it keeps the object's own
# "under"/"bad" colours.  Set them explicitly, otherwise the area outside the
# projection (and hp.UNSEEN) is painted with the lowest colour, i.e. blue.
GAUSSIAN.set_under('white')
GAUSSIAN.set_bad('white')
FOREGROUND.set_under('white')
FOREGROUND.set_bad('white')


# -------------------------------------------------------------------- plots --
def plot_maps(maps, titles=None, proj='orth', unit='', vmin=None, vmax=None,
              rot=(175., 50., 0.), cmap=GAUSSIAN, outfile=None, dpi=300):
    """Draw one or several HEALPix maps side by side in Mollweide (``moll``)
    or half-sky orthographic (``orth``) projection and return the figure."""
    if isinstance(maps, np.ndarray) and maps.ndim == 1:
        maps = [maps]
    titles = [''] * len(maps) if titles is None else titles
    fig = plt.figure(figsize=(3.0 * len(maps), 2.6), dpi=dpi)
    view = hp.mollview if proj == 'moll' else hp.orthview
    extra = {} if proj == 'moll' else {'half_sky': True, 'rot': list(rot)}
    for i, m in enumerate(maps):
        view(np.asarray(m, dtype=np.float64), cmap=cmap, min=vmin, max=vmax,
             title=titles[i], unit=unit, fig=fig.number,
             sub=(1, len(maps), i + 1), **extra)
        hp.graticule(dpar=30., dmer=30., ls='-', lw=0.1)
    if outfile is not None:
        fig.savefig(outfile, bbox_inches='tight', pad_inches=0.1)
    return fig


# --------------------------------------------------------------- harmonics ---
def iqu2teb(iqu, mask=None, teb='te', lmax=None, return_alm=False):
    """Convert an IQU map into the requested T/E/B maps (or a_lm).

    ``mask`` is applied to the map before the transform and to zero the output
    maps.  The transforms use ``use_weights=False``, so no HEALPix weight file
    is required.
    """
    iqu = np.asarray(iqu, dtype=np.float64)
    nside = hp.get_nside(iqu[0])
    mask = (np.ones(hp.nside2npix(nside)) if mask is None
            else np.asarray(mask, dtype=np.float64))
    alms = hp.map2alm(iqu * np.array([mask, mask, mask]), lmax=lmax,
                      use_weights=False)
    teb = teb.lower()
    return np.array([alms[i] if return_alm
                     else hp.alm2map(alms[i], nside, lmax=lmax, pol=False)
                     * (mask > 0)
                     for i, pol in enumerate('teb') if pol in teb])


def clean_bmap(iqu, mask, lmax=None):
    """E-mode recycling: fit and subtract the T->B leakage template.

    The template ``B_t`` is the B map of the T+E part of the input, and the
    cleaned map is ``B_c - b0 - b1 * B_t``, evaluated on the mask and fitted by
    a linear regression over the unmasked pixels.
    """
    iqu = np.asarray(iqu, dtype=np.float64)
    mask = np.asarray(mask, dtype=np.float64)
    nside = hp.get_nside(iqu[0])
    lmax = 3 * nside - 1 if lmax is None else lmax
    m3 = np.array([mask, mask, mask])

    alm = hp.map2alm(iqu * m3, lmax=lmax, use_weights=False)
    b_c = hp.alm2map(alm[2], nside, lmax=lmax, pol=False) * mask

    alm_te = np.zeros_like(alm)
    alm_te[:2] = alm[:2]
    iqu_te = hp.alm2map(alm_te, nside, lmax=lmax, pol=True)
    b_t = hp.alm2map(hp.map2alm(iqu_te * m3, lmax=lmax, use_weights=False)[2],
                     nside, lmax=lmax, pol=False) * mask

    good = mask > 0.9
    x, y = b_t[good], b_c[good]
    slope = np.cov(x, y)[0, 1] / np.var(x, dtype=np.float64)
    intercept = np.mean(y, dtype=np.float64) - slope * np.mean(x, dtype=np.float64)
    return (b_c - intercept - slope * b_t) * mask
