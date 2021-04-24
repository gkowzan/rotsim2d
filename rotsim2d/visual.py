"""Visualization procedures."""
# * Imports
from typing import Optional, Tuple, Dict
from pathlib import Path
import subprocess as subp
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as grd
import matplotlib.cm as cm
from knickknacks.experiment import find_index

# * Matplotlib visualization
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


def plot2d_im(freqs, spec2d, spec_linear=None, scale='symlog', line=True, pthresh=100.0, absmax=None, fig_kwargs={}):
    """2D imshow plot with decent defaults."""
    # pylint: disable=too-many-locals,too-many-arguments,dangerous-default-value
    extent = make_extent(freqs[0], freqs[1])
    cmap = cm.get_cmap('RdBu').reversed()

    fig = plt.figure(**fig_kwargs)
    if spec_linear is not None:
        gs = grd.GridSpec(nrows=2, ncols=2, height_ratios=[1, 6], width_ratios=[20, 1], figure=fig)
        ax1d = plt.subplot(gs[0, 0])
        ax2d = plt.subplot(gs[1, 0])
        axcbar = plt.subplot(gs[1, 1])
    else:
        gs = grd.GridSpec(nrows=1, ncols=2, width_ratios=[20, 1], figure=fig)
        ax2d = plt.subplot(gs[0])
        axcbar = plt.subplot(gs[1])

    if absmax is None:
        absmax = np.max(np.abs(spec2d))
    linthresh = absmax/pthresh
    if scale == 'symlog':
        cset = ax2d.imshow(spec2d, cmap=cmap, aspect='auto', extent=extent,
                           clim=(-absmax, absmax), origin='lower',
                           norm=colors.SymLogNorm(linthresh=linthresh, vmax=absmax, vmin=-absmax))
    elif scale == 'linear':
        cset = ax2d.imshow(spec2d, cmap=cmap, aspect='auto', extent=extent,
                           clim=(-absmax, absmax), origin='lower')
    ax2d.set(xlabel=r'Probe (cm$^{-1}$)', ylabel=r'Pump (cm$^{-1}$)')

    if line:
        p = (ax2d.get_xlim()[1]+ax2d.get_xlim()[0])/2
        ax2d.axline((p, p), slope=1.0, color='gray', alpha=0.5, zorder=1)
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


def plot2d_animation(freqs: Tuple[np.ndarray], spec3d: np.ndarray, absmax: Optional[float]=None,
                     fig_kwargs: Dict={}):
    """Prepare Figure and callbacks for `matplotlib.animation.FuncAnimation`.

    Parameters
    ----------
    freqs : tuple of ndarray
        Array if pump frequencies, array of waiting times and array of
        probe frequencies.
    spec3d : ndarray
        3D array of data to plot, second dimension is the animated time.
    absmax : float
        Data limits for colorbar, taken from data otherwise.
    fig_kwargs : dict
        Keyword arguments for `plt.figure`.

    Returns
    -------
    fig : matplotlib.figure.Figure
    init : function
        Initialization closure function for FuncAnimation.
    update : function
        Update closure function for FuncAnimation.
    frames : ndarray
        NumPy array of frame numbers. 
    """
    extent = vis.make_extent(freqs[0], freqs[2])
    ts2 = freqs[1]
    cmap = cm.get_cmap('RdBu').reversed()
    fig = plt.figure(**fig_kwargs)
    gs = fig.add_gridspec(nrows=1, ncols=2, width_ratios=[20, 1], figure=fig)
    ax2d = fig.add_subplot(gs[0])
    axcbar = fig.add_subplot(gs[1])
    if absmax is None:
        absmax = np.max(np.abs(spec2d))
    cset = ax2d.imshow(np.zeros((spec3d.shape[0], spec3d.shape[2])),
                       cmap=cmap, aspect='auto', extent=extent, clim=(-absmax, absmax),
                       origin='lower')
    ax2d.set(xlabel=r'Probe (cm$^{-1}$)', ylabel=r'Pump (cm$^{-1}$)')
    title = ax2d.text(0.5, 0.99, "$t_2$", ha='center', va='top', transform=ax2d.transAxes)
    axcbar = fig.colorbar(cset, ax=ax2d, cax=axcbar)
    fig.set_constrained_layout_pads(wspace=0.01, hspace=0.01, h_pad=0.01, w_pad=0.01)

    def init():
        return [cset, title]

    def update(i):
        cset.set_data(spec3d[:, i, :])
        title.set_text('$t_2={:.2f}$ ps'.format(ts2[i]*1e12))
        return [cset, title]

    return fig, init, update, np.arange(spec2d.shape[1])


def plot2d_scatter(pl, fig_dict=None, line=True, vminmax=None, fig_kwargs={}, scatter_kwargs={}):
    # pylint: disable=dangerous-default-value,too-many-arguments
    if fig_dict:
        fig = fig_dict['fig']
        ax = fig_dict['ax']
        axcbar = fig_dict['axcbar']
    else:
        fig = plt.figure(**fig_kwargs)
        gs = grd.GridSpec(1, 2, width_ratios=[20, 1], figure=fig)
        ax = plt.subplot(gs[0])
        axcbar = plt.subplot(gs[1])
    if vminmax:
        vmin = -np.abs(vminmax)
    else:
        vmin = -np.max(np.abs(pl.sigs))

    sc = ax.scatter(pl.probes, pl.pumps, s=10.0, c=pl.sigs, cmap='seismic',
                    vmin=vmin, vmax=-vmin, **scatter_kwargs)
    if not fig_dict:
        fig.colorbar(sc, ax=ax, cax=axcbar)
    if line:
        ax.axline((pl.probes[0], pl.probes[0]),
                  slope=1.0, color='gray', alpha=0.5, zorder=0)

    fig.set_constrained_layout_pads(wspace=0.01, hspace=0.01, h_pad=0.01, w_pad=0.01)

    return {'ax': ax, 'axcbar': axcbar, 'fig': fig, 'sc': sc}


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


# * LaTeX double-sided Feynmann diagrams
LATEX_PRE = r"""\documentclass[tikz,border=3.14mm]{standalone}
\usetikzlibrary{matrix,fit,positioning}
\medmuskip=0mu
\tikzset{
mymat/.style={
    matrix of math nodes,
    left delimiter=|,right delimiter=|,
    align=center,
    column sep=-\pgflinewidth
},
mymats/.style={
    mymat,
    nodes={draw,fill=#1}
}
}

\begin{document}
\begin{tikzpicture}[font=\sffamily,
    every left delimiter/.style={xshift=.4em},
    every right delimiter/.style={xshift=-.4em}]
"""

LATEX_POST = r"""\end{tikzpicture}
\end{document}"""

dj_to_letter = {-2: "O", -1: "P", 0: "Q", 1: "R", 2: "S"}


def tikz_abstract_format(dnu: int, dj: int):
    if dnu==0 and dj==0:
        return "0"
    return str(dnu)+dj_to_letter[dj]


def tikz_diagram(leaf, index=1, first=True, direction=False, abstract=None, hspace="2cm"):
    """Return TikZ double-sided Feynmann diagram.

    With `abstract` ketbras are labeled relative to some (nu, j) ground state.

    Parameters
    ----------
    index : str
        Internal Tikz index of the current diagram.
    first : bool
        First diagram in the document.
    direction : bool
        Add phase-matching direction label.
    abstract : tuple of int, optional
        Tuple of ground state quantum numbers, (nu, j).
    hspace : str
        Horizontal space between pathways in LaTeX dimensions.
    """
    diag_rows = []
    for i, kb in enumerate(leaf.ketbras()):
        if not abstract:
            diag_rows.append(r"$|{:d},{:d}\rangle\langle {:d},{:d}|$\\".format(kb.knu, kb.kj, kb.bnu, kb.bj))
        else:
            taf = tikz_abstract_format
            gnu, gj = abstract
            knu, kj, bnu, bj = kb.knu-gnu, kb.kj-gj, kb.bnu-gnu, kb.bj-gj
            if i==0:
                diag_rows.append(r"$|0\rangle\langle 0|$\\")
            else:
                diag_rows.append(r"$|{:s}\rangle\langle {:s}|$\\".format(
                    taf(knu, kj), taf(bnu, bj)))
    diag_rows.reverse()

    diag_arrows = []
    ints = leaf.interactions()
    nints = len(ints)
    for i, li in enumerate(ints):
        print(li)
        if li.side == -1:
            # label_side = "right"
            arr_side = "east"
            xshift = "4mm"
        else:
            # label_side = "left"
            arr_side = "west"
            xshift = "-4mm"

        if li.side == -1 and li.sign == -1:
            arr = "stealth-"
            yshift = "-4mm"
        elif li.side == -1 and li.sign == 1:
            arr = "-stealth"
            yshift = "4mm"
        elif li.side == 1 and li.sign == 1:
            arr = "stealth-"
            yshift = "-4mm"
        elif li.side == 1 and li.sign == -1:
            arr = "-stealth"
            yshift = "4mm"

        if li.readout:
            arr = arr + ",dashed"

        # if li.readout:
        #     klabel = r"$\vec{\mu}$"
        # else:
        #     klabel = r"$\vec{{k}}_{:d}$".format(i+1)

        arrow = r"\draw[{arr:s}] (mat{index:d}-{ri:d}-1.south -| mat{index:d}.{arr_side:s}) -- ++({xshift:s},{yshift:s});".format(arr=arr, ri=nints-i, arr_side=arr_side, xshift=xshift, yshift=yshift, index=index)
        diag_arrows.append(arrow)

    note = ''
    if direction:
        note = r"\node[below] (mat{index:d}label) at (mat{index:d}.south)".format(index=index)
        if leaf.is_SI():
            note += r"{$S_{\mathrm{I}}$};"
        elif leaf.is_SII():
            note += r"{$S_{\mathrm{II}}$};"
        elif leaf.is_SIII():
            note += r"{$S_{\mathrm{III}}$};"

    if first:
        header = r"""\matrix[mymat] at (0,0) (mat{index:d})""".format(index=index)
    else:
        header = r"""\matrix[mymat,right={hspace:s} of mat{prev_index:d}] (mat{index:d})""".format(prev_index=index-1,
                                                                                                   index=index,
                                                                                                   hspace=hspace)
    ret = header + """ {{ {rows:s} }};
{arrows:s}
{note:s}
""".format(rows="\n".join(diag_rows), arrows="\n".join(diag_arrows), note=note)

    return ret


def tikz_diagrams(tree, direction=False, abstract=None, hspace="2cm"):
    """Return tikz code for series of double-sided Feynmann diagrams."""
    leaves = tree.leaves
    diagrams = [tikz_diagram(leaves[0], direction=direction, abstract=abstract)]
    for i, l in enumerate(leaves[1:], start=2):
        diagrams.append(tikz_diagram(l, index=i, first=False, direction=direction, abstract=abstract, hspace=hspace))

    return "\n".join(diagrams)


def latex_compile(path, doc):
    """Save `doc` to `path` and compile it."""
    doc_path = Path(path)
    doc_dir = doc_path.parent
    doc_dir.mkdir(parents=True, exist_ok=True)
    doc_path.write_text(doc)

    subp.run(['pdflatex', doc_path.name], cwd=doc_dir, check=True)
