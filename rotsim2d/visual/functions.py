"""Visualization procedures."""
# * Imports
from typing import Optional, Tuple, Dict, Sequence, List, Mapping, Callable, Union
from copy import deepcopy
import re
from pathlib import Path
import subprocess as subp
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as grd
import matplotlib.cm as cm
from rotsim2d.utils import find_index
import rotsim2d.dressedleaf as dl
import rotsim2d.symbolic.functions as sym
from rotsim2d.visual.latexprinter import latex
import molspecutils.molecule as mol

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


def plot2d_im(freqs, spec2d, spec_linear=None, scale='symlog', line=True,
              pthresh=100.0, absmax=None, fig_kwargs={}):
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


def plot1d_animation(coords: Tuple[np.ndarray, np.ndarray],
                     spec2d: np.ndarray,
                     ylim: Optional[Tuple[float, float]]=None,
                     subplots_kwargs: Mapping={},
                     update_title: Callable=None,
                     update_data_limits: bool=False):
    """Prepare Figure and callbacks for `FuncAnimation` for 1D spectra.

    The returned objects can be used as::

    >> fanim = FuncAnimation(fig, update, frames, init, blit=True)

    Parameters
    ----------
    freqs : tuple of ndarray
        Array of animation coordinates and array of frequencies.
    spec2d : ndarray
        2D array of data to plot, first dimension is the animated coordinate.
    ylim : tuple of float
        Y limits, taken from data otherwise.
    fig_kwargs : dict
        Keyword arguments for `plt.figure`.
    update_title : Callable
        Function called to update the title.
    update_data_limits : bool
        Recalculate data limits (`absmax`) for each frame.

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
    fig, ax = plt.subplots(**subplots_kwargs)
    line, = ax.plot([], [])
    title = ax.text(0.5, 0.99, "", ha="center", va="top", transform=ax.transAxes)
    ax.set_xlim(np.min(coords[1]), np.max(coords[1]))
    if ylim:
        ymin, ymax = ylim
    else:
        ymin, ymax = np.min(spec2d), np.max(spec2d)
    ax.set_ylim(ymin, ymax)
    fig.set_constrained_layout_pads(wspace=0.01, hspace=0.01, h_pad=0.01, w_pad=0.01)

    def init():
        return [line, title]

    def update(i):
        line.set_data(coords[1], spec2d[i])
        if update_title is not None:
            update_title(title, i)
        if update_data_limits:
            ymin, ymax = np.min(spec2d[i]), np.max(spec2d[i])
            ax.set_ylim(ymin, ymax)

        return [line, title]

    return fig, init, update, np.arange(len(coords[0]))


def plot2d_animation(freqs: Tuple[np.ndarray, ...], spec3d: np.ndarray,
                     absmax: Optional[float]=None, fig_kwargs: Mapping={},
                     update_title: Callable=None,
                     update_data_limits: bool=False):
    """Prepare Figure and callbacks for `FuncAnimation` for 2D spectra.

    The returned objects can be used as::

    >> fanim = FuncAnimation(fig, update, frames, init, blit=True)

    Parameters
    ----------
    freqs : tuple of ndarray
        Array of pump frequencies, array of waiting times and array of
        probe frequencies.
    spec3d : ndarray
        3D array of data to plot, second dimension is the animated time.
    absmax : float
        Data limits for colorbar, taken from data otherwise.
    fig_kwargs : dict
        Keyword arguments for `plt.figure`.
    update_title : Callable
        Function called to update the title.
    update_data_limits : bool
        Recalculate data limits (`absmax`) for each frame.

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
    extent = make_extent(freqs[0], freqs[2])
    ts2 = freqs[1]
    cmap = cm.get_cmap('RdBu').reversed()
    fig = plt.figure(**fig_kwargs)
    gs = fig.add_gridspec(nrows=1, ncols=2, width_ratios=[20, 1], figure=fig)
    ax2d = fig.add_subplot(gs[0])
    axcbar = fig.add_subplot(gs[1])
    if absmax is None:
        absmax = np.max(np.abs(spec3d))
    cset = ax2d.imshow(np.zeros((spec3d.shape[0], spec3d.shape[2])),
                       cmap=cmap, aspect='auto', extent=extent, clim=(-absmax, absmax),
                       origin='lower')
    ax2d.set(xlabel=r'Probe (cm$^{-1}$)', ylabel=r'Pump (cm$^{-1}$)')
    title = ax2d.text(0.5, 0.99, "title", ha='center', va='top', transform=ax2d.transAxes)
    axcbar = fig.colorbar(cset, ax=ax2d, cax=axcbar)
    fig.set_constrained_layout_pads(wspace=0.01, hspace=0.01, h_pad=0.01, w_pad=0.01)

    def init():
        return [cset, title]

    def update(i):
        cset.set_data(spec3d[:, i, :])
        if update_title is not None:
            update_title(title, i)
        if update_data_limits:
            absmax = np.max(np.abs(spec3d[:, i, :]))
            cset.set_clim(-absmax, absmax)
        return [cset, title]

    return fig, init, update, np.arange(spec3d.shape[1])


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
        vmin = -np.max(np.abs(pl.intensities))

    skwargs = dict(s=10.0, cmap='seismic')
    skwargs.update(scatter_kwargs)
    sc = ax.scatter(pl.probes, pl.pumps, c=pl.intensities, vmin=vmin, vmax=-vmin,
                    **skwargs)
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
\newcommand{\textoverline}[1]{$\overline{\mbox{#1}}$}
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


def trans_label2latex(label: str) -> str:
    """Write ``label`` with overline instead of parentheses."""
    return re.sub(r'\(([A-Z0-9]{1,2})\)', r'\\ensuremath{{}\\mkern1mu\\overline{\\mkern-1mu\\mbox{\1}}}', label)


def R_label2latex(label: str) -> str:
    lat_label = r'\ensuremath{\mbox{R}_' + label[1]
    if len(label) == 3:
        lat_label += r'^\ast'

    return  lat_label + '}'


def tikz_diagram(leaf, index=1, first=True, direction=False, trans_label=None,
                 R_label=False, abstract=None, hspace="2cm") -> str:
    """Return TikZ double-sided Feynmann diagram.

    With `abstract` ketbras are labeled relative to some (nu, j) ground state.

    Parameters
    ----------
    leaf : rotsim2d.pathways.KetBra
        KetBra leaf instance to render.
    index : str, optional
        Internal Tikz index of the current diagram.
    first : bool, optional
        First diagram in the document.
    direction : bool, optional
        Add phase-matching direction label.
    trans_label : {'proper', 'degenerate'}, optional
        Add transition label, either unambiguous or degenerate one.
    R_label: bool, optional
        Add R_i, i=1,...,8, label.
    abstract : tuple of int, optional
        Tuple of ground state quantum numbers, (nu, j).
    hspace : str
        Horizontal space between pathways in LaTeX dimensions.
    """
    diag_rows = []
    for i, kb in enumerate(leaf.ketbras()):
        if not abstract:
            diag_rows.append(r"$|{:d},{:d}\rangle\langle {:d},{:d}|$\\".format(kb.ket.nu, kb.ket.j,
                                                                               kb.bra.nu, kb.bra.j))
        else:
            asl = dl.abstract_state_label
            gstate = mol.DiatomState(*abstract)
            if i==0:
                diag_rows.append(r"$|0\rangle\langle 0|$\\")
            else:
                diag_rows.append(r"$|{:s}\rangle\langle {:s}|$\\".format(
                    asl(kb.ket, gstate), asl(kb.bra, gstate)))
    diag_rows.reverse()

    diag_arrows = []
    ints = leaf.interactions()
    nints = len(ints)
    for i, li in enumerate(ints):
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

    notes = []
    if direction:
        if leaf.is_SI():
            notes.append(r"$S_{\mathrm{I}}$")
        elif leaf.is_SII():
            notes.append(r"$S_{\mathrm{II}}$")
        elif leaf.is_SIII():
            notes.append(r"$S_{\mathrm{III}}$")
    if trans_label:
        if trans_label == "proper":
            notes.append(trans_label2latex(dl.Pathway(leaf).trans_label))
        elif trans_label == "degenerate":
            notes.append(dl.Pathway(leaf).trans_label_deg)
    if R_label:
        notes.append(R_label2latex(leaf.R_label()))
    if notes:
        notes = ' {' + ', '.join(notes) + '};'
        note = r"\node[below] (mat{index:d}label) at (mat{index:d}.south)".format(index=index) +\
            notes
    else:
        note = ''

    if first:
        header = r"""\matrix[mymat] at (0,0) (mat{index:d})""".format(index=index)
    else:
        header = r"""\matrix[mymat,right={hspace:s} of mat{prev_index:d}] (mat{index:d})""".format(
            prev_index=index-1, index=index, hspace=hspace)
    ret = header + """ {{ {rows:s} }};
{arrows:s}
{note:s}
""".format(rows="\n".join(diag_rows), arrows="\n".join(diag_arrows), note=note)

    return ret


def tikz_diagrams(tree, direction=False, trans_label=None, R_label=False,
                  abstract=None, hspace="2cm"):
    """Return tikz code for series of double-sided Feynmann diagrams.

    Parameters
    ----------
    tree : rotsim2d.pathways.KetBra
        Root of KetBra excitation tree. Generate pathway diagrams for all leaves.
    direction : bool, optional
        Add phase-matching direction label.
    trans_label : {'proper', 'degenerate'}, optional
        Add transition label, either unambiguous or degenerate one.
    R_label: bool, optional
        Add R_i, i=1,...,8, label.
    abstract : tuple of int, optional
        Tuple of ground state quantum numbers, (nu, j).
    hspace : str
        Horizontal space between pathways in LaTeX dimensions.

    Returns
    -------
    str
        Tikz code to for diagrams. Needs to be placed in a LaTeX document to
        compile.

    Examples
    --------
    Create KetBras, filter and render diagrams:

    >>> import rotsim2d.pathways as pw
    >>> import rotsim2d.visual as vis
    >>> kbs = pw.gen_pathways([5], meths=[pw.only_SI], rotor='linear')
    >>> vis.latex_compile('diagrams.tex',
                          vis.latex_document(
                              vis.tikz_diagrams(kbs, abstract=(0, 5))))
    """
    leaves = tree.leaves
    diagrams = [tikz_diagram(leaves[0], direction=direction, trans_label=trans_label, R_label=R_label,
                             abstract=abstract)]
    for i, l in enumerate(leaves[1:], start=2):
        diagrams.append(tikz_diagram(l, index=i, first=False,
                                     direction=direction, trans_label=trans_label, R_label=R_label,
                                     abstract=abstract, hspace=hspace))

    return "\n".join(diagrams)


def latex_document(fragment: str) -> str:
    """Wrap provided fragment in a LaTeX document."""
    return LATEX_PRE + fragment + LATEX_POST


def latex_compile(path: Union[str, Path], doc: str):
    """Save `doc` to `path` and compile it."""
    doc_path = Path(path)
    doc_dir = doc_path.parent
    doc_dir.mkdir(parents=True, exist_ok=True)
    doc_path.write_text(doc)

    subp.run(['pdflatex', doc_path.name], cwd=doc_dir, check=True)


# * LaTeX table of coefficients
def classified_table_prep(classified: Dict[sym.RFactor, dl.Pathway],
                          highj: bool=False) -> List[List]:
    """Convert dict of R-factors and pathways to a table of coefficients.

    Returns a list of rows. First element of a row is a list of `trans_label`s,
    remaining elements of a row are `c` coefficients.
    """
    def labels(l: Sequence[dl.Pathway], highj: bool=False):
        l1 = [pw.trans_label for pw in l]
        if highj:
            l1 = list(set([re.sub(r'[)(0-9]', '', x) for x in l1]))
            l1.sort()
            l1 = [r'\textbf{{{:s}}}'.format(v) for v in l1]
        else:
            l1.sort(key=lambda x: re.sub(r'[)(0-9]', '', x))
            l1 = [trans_label2latex(v) for v in l1]
        return l1

    # represent pathways with transition labels
    classified = {k: labels(v, highj) for k, v in classified.items()}
    table = []
    for k, v in classified.items():
        if highj:
            table.append([v] + list(k.tuple[1:]))
        else:
            table.append([v] + list(k.tuple))

    return table


def merge_tables(dtable: Mapping[str, List[List]]) -> List[List]:
    """Merge different tables of coefficients into one.

    Duplicate `trans_label`s are distinguished by subscripted `dtable` keys.
    """
    dtable = deepcopy(dtable)
    table = []
    keys = list(dtable.keys())
    # uniquify labels
    uniquify = set()
    for row in dtable[keys[0]]:
        labels, coeffs = row[0], row[1:]
        for il in range(len(labels)):
            # loop over other tables, rows, find labels to uniquify
            for fkey in keys[1:]:
                for frow in dtable[fkey]:
                    flabels, fcoeffs = frow[0], frow[1:]
                    for ifl in range(len(flabels)):
                        if labels[il] == flabels[ifl] and coeffs != fcoeffs:
                            uniquify.add(labels[il])
    # actually uniquify the labels
    for label in uniquify:
        for key in keys:
            for row in dtable[key]:
                labels = row[0]
                for il in range(len(labels)):
                    if labels[il] == label:
                        row[0][il] = labels[il]+r'\textsubscript{{{:s}}}'.format(key)
    # now that labels are either unique or mean the same thing, merge the tables
    while len(dtable[keys[0]]):
        row = dtable[keys[0]].pop()
        labels, coeffs = set(row[0]), row[1:]
        for fkey in keys[1:]:
            for ifrow in range(len(dtable[fkey])):
                if coeffs == dtable[fkey][ifrow][1:]:
                    frow = dtable[fkey].pop(ifrow)
                    labels.update(frow[0])
                    break
        labels = list(labels)
        labels.sort()
        table.append([labels] + coeffs)

    for key in keys:
        assert len(dtable[key]) == 0

    return table


def classified_table_render(table: List[List]) -> str:
    """Render table of coefficients as LaTeX table body."""
    str_table = [' & '.join([', '.join(x[0])] +
                        [latex(xx, mode='inline') for xx in x[1:]])
             for x in table]
    str_table.sort()

    return '\\\\\n'.join(str_table)


def classified_table(classified: Dict[sym.RFactor, dl.Pathway],
                     highj: bool=False) -> str:
    """Format dict of R-factors and pathways as a table of coefficients.

    Uses custom LaTeX printer from `.latexprinter` module. With `highj` True
    don't print the first coefficient.
    """
    table = classified_table_prep(classified, highj)

    return classified_table_render(table)
