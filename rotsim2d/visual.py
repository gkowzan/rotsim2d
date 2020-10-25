"""Visualization procedures."""
# * Imports
from pathlib import Path
import subprocess as subp
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as grd
from shed.experiment import find_index

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


def plot2d_im(freqs, spec2d, spec_linear=None, scale='symlog'):
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
    if scale == 'symlog':
        cset = ax2d.imshow(spec2d, cmap='seismic', aspect='auto', extent=extent,
                           clim=(-absmax, absmax), origin='lower',
                           norm=colors.SymLogNorm(linthresh=linthresh, vmax=absmax, vmin=-absmax))
    elif scale == 'linear':
        cset = ax2d.imshow(spec2d, cmap='seismic', aspect='auto', extent=extent,
                           clim=(-absmax, absmax), origin='lower')
    ax2d.set(xlabel=r'Probe (cm$^{-1}$)', ylabel=r'Pump (cm$^{-1}$)')
    ax2d.axline((ax2d.get_xlim()[0], ax2d.get_ylim()[0]),
                (ax2d.get_xlim()[1], ax2d.get_ylim()[1]),
                color='gray', alpha=0.5)
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


# * LaTeX double-sided Feynmann diagrams
LATEX_PRE = r"""\documentclass[tikz,border=3.14mm]{standalone}
\usetikzlibrary{matrix,fit,positioning}
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
Bazinga!\\
\begin{tikzpicture}[font=\sffamily,
    every left delimiter/.style={xshift=.4em},
    every right delimiter/.style={xshift=-.4em}]
"""

LATEX_POST = r"""\end{tikzpicture}
\end{document}"""

def tikz_diagram(leaf, index=1, first=True):
    """Return TikZ double-sided Feynmann diagram."""
    diag_rows = []
    for kb in reversed(leaf.ketbras()):
        diag_rows.append(r"$|{:d},{:d}\rangle\langle {:d},{:d}|$\\".format(kb.knu, kb.kj, kb.bnu, kb.bj))

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

    if first:
        header = r"""\matrix[mymat] at (0,0) (mat{index:d})""".format(index=index)
    else:
        header = r"""\matrix[mymat,right=2cm of mat{prev_index:d}] (mat{index:d})""".format(prev_index=index-1,
                                                                                            index=index)
    ret = header + """ {{ {rows:s} }};
{arrows:s}
""".format(rows="\n".join(diag_rows), arrows="\n".join(diag_arrows))

    return ret


def tikz_diagrams(tree):
    """Return tikz code for series of double-sided Feynmann diagrams."""
    leaves = tree.leaves
    diagrams = [tikz_diagram(leaves[0])]
    for i, l in enumerate(leaves[1:], start=2):
        diagrams.append(tikz_diagram(l, index=i, first=False))

    return "\n".join(diagrams)


def latex_compile(path, doc):
    """Save `doc` to `path` and compile it."""
    doc_path = Path(path)
    doc_dir = doc_path.parent
    doc_dir.mkdir(parents=True, exist_ok=True)
    doc_path.write_text(doc)

    subp.run(['pdflatex', doc_path.name], cwd=doc_dir, check=True)
