"""Visualization procedures."""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as grd
from shed.experiment import find_index


def make_extent(t1s, t2s, scale=1.0):
    """Return extents for imshow or contour(f).

    Parameters
    ----------
    t1s: ndarray
        Time/frequency vertical axis (pump, w1).
    t2s: ndarray
        Time/frequency horizontal axis (probe, w3).

    Returns
    -------
    list of float
    """
    t1s = t1s*scale
    t2s = t2s*scale
    dt1 = (t1s[1]-t1s[0])/2
    dt2 = (t2s[1]-t2s[0])/2

    return [t2s[0]-dt2, t2s[-1]+dt2, t1s[0]-dt1, t1s[-1]+dt1]


def plot2d_im(freqs, spec2d, spec_linear=None):
    """2D imshow plot with decent defaults."""
    extent = make_extent(freqs[0], freqs[1])

    fig = plt.figure()
    if spec_linear is not None:
        gs = grd.GridSpec(nrows=2, ncols=2, height_ratios=[1, 6], width_ratios=[20, 1], figure=fig)
        ax1d = plt.subplot(gs[0, 0])
        ax2d = plt.subplot(gs[1, 0])
        axcbar = plt.subplot(gs[1, 1])
    else:
        gs = grd.GridSpec(nrows=1, ncols=2, width_ratios=[20, 1], figure=fig)
        ax2d = plt.subplot(gs[0])
        axcbar = plt.subplot(gs[1])

    absmax = np.max(np.abs(spec2d))
    linthresh = absmax/100.0
    cset = ax2d.imshow(spec2d, cmap='seismic', aspect='auto', extent=extent,
                       clim=(-absmax, absmax),
                       norm=colors.SymLogNorm(linthresh=linthresh, vmax=absmax, vmin=-absmax))
    ax2d.set(xlabel=r'Probe (cm$^{-1}$)', ylabel=r'Pump (cm$^{-1}$)')
    ax2d.axline((ax2d.get_xlim()[0], ax2d.get_ylim()[0]),
                (ax2d.get_xlim()[1], ax2d.get_ylim()[1]),
                color='gray')
    axcbar = fig.colorbar(cset, ax=ax2d, cax=axcbar)

    if spec_linear is not None:
        ax1d.plot(spec_linear[0], spec_linear[1])
        ax1d.set(
            xlim=(spec_linear[0][0], spec_linear[0][-1]),
            xticklabels=[],
        )

    fig.set_constrained_layout_pads(wspace=0.01, hspace=0.01, h_pad=0.01, w_pad=0.01)

    if spec_linear is not None:
        return {'ax2d': ax2d, 'axcbar': axcbar, 'ax1d': ax1d, 'fig': fig}
    else:
        return {'ax2d': ax2d, 'axcbar': axcbar, 'ax1d': None, 'fig': fig}


def plot1d_probe(pump_wn, freqs, spec2d, fig_dict=None, **plot_kwargs):
    """Plot probe spectrum at `pump_wn`."""
    ipump = find_index(freqs[0], pump_wn)

    if fig_dict is None:
        fig, ax = plt.subplots()
    else:
        fig, ax = fig_dict['fig'], fig_dict['ax']
    ax.plot(freqs[1], spec2d[ipump], **plot_kwargs)

    fig.set_constrained_layout_pads(wspace=0.01, hspace=0.01, h_pad=0.01, w_pad=0.01)

    return dict(fig=fig, ax=ax)
