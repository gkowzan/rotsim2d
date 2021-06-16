"""Propagate time evolution of coherences based on pathway trees.

**TODO:**

- refactor in terms of :class:`rotsim2d.dressedleaf.DressedPathway` and get rid
  of `pop` attribute in `KetBra`.
"""
import logging
from typing import List, Union, Tuple, Optional, Sequence
import numpy as np
import scipy.constants as C
try:
    import pyfftw
    pyfftw.config.PLANNER_EFFORT = 'FFTW_ESTIMATE'
    import pyfftw.interfaces.scipy_fftpack as fftp
except ModuleNotFoundError:
    import scipy.fftpack as fftp
from rotsim2d import pathways as pw
from rotsim2d.couple import four_couple
import rotsim2d.dressedleaf as dl

e2i = C.epsilon_0*C.c/2         #: Electric field to intensity
xs2cm = 1e2/C.c                 #: Cross section in m^2 Hz to cm^2 cm^{-1}
ignore_missing = False          #: silently ignore dynamics of missing lines

log = logging.getLogger(__name__)


def aid(x):
    return x.__array_interface__['data'][0]


def pws_autospan(pws, margin=5.0*30e9):
    """Return min/max pump/probe frequencies from `pws`."""
    pump_freqs, probe_freqs = [pw.nu(0) for pw in pws], [pw.nu(2) for pw in pws]
    pump_min, pump_max = min(pump_freqs), max(pump_freqs)
    probe_min, probe_max = min(probe_freqs), max(probe_freqs)

    return ((pump_min-margin, pump_max+margin, 2*margin+pump_max-pump_min),
            (probe_min-margin, probe_max+margin, 2*margin+probe_max-probe_min))


def aligned_fs(fsmin: float, fsmax: float, df: float):
    def align(f: float):
        return np.ceil(f/df).astype(np.int) if f < 0 else np.floor(f/df).astype(np.int)
    return np.arange(align(fsmin), align(fsmax)+1)*df


def leaf_term(nu: float, gam: float, coord: np.ndarray, domain: str):
    """Return either time or frequency domain response."""
    if domain == 'f':
        return 1.0/(gam - 1.0j*(coord-nu))
    elif domain == 't':
        return np.exp(-2.0*np.pi*coord*(1.0j*nu+gam))


def dressed_leaf_response(dl: dl.DressedLeaf, coords: Sequence[Optional[np.ndarray]],
                          domains: Sequence[str], freq_shifts: Optional[Sequence[float]]) -> np.ndarray:
    """Calculate response for a single DressedLeaf."""
    # validate inputs
    if len(coords) != len(domains):
        raise ValueError("len(coords) != len(domains)")
    for d in domains:
        if d not in ('t', 'f'):
            raise ValueError("domain can either be 't' or 'f'")

    freq_shifts = freq_shifts if freq_shifts else [0.0]*len(coords)
    resps = []
    for i, coord in enumerate(coords):
        if coord is None:
            continue

        if dl.nu(i) > 0.0:
            nu = dl.nu(i)-freq_shifts[i]
        elif dl.nu(i) < 0.0:
            nu = dl.nu(i)+freq_shifts[i]
        else:
            nu = dl.nu(i)

        resps.append(leaf_term(nu, dl.gamma(i), coord, domains[i]))

    resp = resps[0]*dl.intensity()
    if len(resps) > 1:
        for next_resp in resps[1:]:
            resp = resp*next_resp
    return resp


class CrossSectionMixin:
    def leaf_cross_section(self, leaf: pw.KetBra, times: List[float]) -> np.complex:
        """Calculate interaction cross section for a single leaf."""
        resp = self.leaf_response(leaf, times)
        if leaf.parent.readout:
            leaf = leaf.parent.parent
        pair = ((leaf.bnu, leaf.bj), (leaf.knu, leaf.kj))
        try:
            nu = self.sys_params['line_params'][pair]['nu']
        except KeyError:
            try:
                # this should always work 
                nu = self.sys_params['line_params'][(pair[1], pair[0])]['nu']
            except KeyError as e:
                print(list(leaf.ancestors) + [leaf])
                raise e

        return resp*2*np.pi*nu/C.c/C.epsilon_0/6.0

    def cross_section(self, kb: pw.KetBra, times: List[float]) -> np.complex:
        r"""Calculate interaction cross section.

        The returned quantity has units of `m^2 Hz`. Mulitply by `10e2/C.c` to
        get a HITRAN line intensity (cross section integrated over frequency
        range) in `cm^2 cm^{-1}`.
        """
        return sum(self.leaf_cross_section(leaf, times) for leaf in kb.diagonals()
                   if self.filter(leaf))

    def response(self, kb: pw.KetBra, times: List[float]) -> np.complex:
        """Calculate the total response.

        Parameters
        ----------
        kb: pw.KetBra
            Root of the excitation tree.
        times: list of float
            Delays between pulses (and readout)

        Returns
        -------
        resp: complex
            Response function.
        """
        return sum(self.leaf_response(leaf, times) for leaf in kb.diagonals()
                   if self.filter(leaf))

    def leaf_response(self, leaf: pw.KetBra, freqs: List[Union[np.ndarray, float]]) -> Union[np.ndarray, np.complex]:
        r"""Return single pathway response spectrum."""
        freq_shift = [0.0]*len(freqs) if self.freq_shift is None else self.freq_shift
        resps, const = [], np.complex(1.0)

        kb_series = [x for x in leaf.ancestors if isinstance(x, pw.KetBra)] + [leaf]
        # generate four-fold dipole operator
        wkets, wbras = [], []
        for i in range(1, len(kb_series)):
            kb, kbp = kb_series[i], kb_series[i-1]

            # dipole interaction reduced matrix element
            if kb.parent.side is pw.Side.KET:
                pair = ((kbp.knu, kbp.kj), (kb.knu, kb.kj))
                wkets.insert(0, (kb.kj,  kb.parent.angle))
            else:
                pair = ((kbp.bnu, kbp.bj), (kb.bnu, kb.bj))
                wbras.append((kb.parent.parent.bj, kb.parent.angle))
            try:
                mu = self.sys_params['line_params'][pair]['mu']
            except KeyError:
                mu = self.sys_params['line_params'][(pair[1], pair[0])]['mu']

            if kb.parent.readout:
                const *= mu*kb.root.pop
                wbras.extend(wkets)
                const *= four_couple([x[0] for x in wbras], [x[1] for x in wbras])
                # print("Js:", [x[0] for x in wbras], "Angles:", [x[1] for x in wbras],
                #       'Four-couple:', four_couple([x[0] for x in wbras], [x[1] for x in wbras]))
                break
            const *= 1.0j/C.hbar*kb.parent.side*mu
            if freqs[i-1] is None:
                continue

            # nu and gam; pairs in line_params have lower nu first
            pair = ((kb.bnu, kb.bj), (kb.knu, kb.kj))
            try:
                nu = self.sys_params['line_params'][pair]['nu'] - freq_shift[i-1]
                gam = self.sys_params['line_params'][pair]['gam']
            except KeyError:
                try:
                    nu = -self.sys_params['line_params'][(pair[1], pair[0])]['nu'] + freq_shift[i-1]
                    gam = self.sys_params['line_params'][(pair[1], pair[0])]['gam']
                except KeyError as e:
                    if ignore_missing:
                        resps.append(0.0)
                        log.debug("Missing `nu` and `gam` for %s, ignoring pathway", pair)
                        break
                    else:
                        raise e

            resps.append(self.leaf_term(nu, gam, freqs[i-1]))

        if len(resps) == 1:
            return resps[0]*const
        else:
            resp = resps[0]*const
            for i in range(len(resps)-1):
                resp = resp*resps[i+1]
            return resp

    
class Spectrum(CrossSectionMixin):
    """Calculate spectrum based on excitation tree.

    Parameters
    ----------
    sys_params: dict of dict
        {'elevels': {(nu, j): elevel}, 'populations': {(nu, j): pop},
        'line_params': {((nupp, jpp), (nup, jp)): {'mu': dipole element, 'gam':  dephasing,
        'nu': position}}}.
    """
    def __init__(self, sys_params: dict, freq_shift: List[float]=None, filter=None):
        self.sys_params = sys_params
        if filter is None:
            self.filter = lambda x: True
        else:
            self.filter = filter
        self.freq_shift = freq_shift

    def leaf_term(self, nu, gam, freqs):
        return 1.0/(gam - 1.0j*(freqs-nu))


class Propagator(CrossSectionMixin):
    """Calculate useful physical quantities based on excitation tree.

    Parameters
    ----------
    sys_params: dict of dict
        {'elevels': {(nu, j): elevel}, 'populations': {(nu, j): pop},
        'line_params': {((nupp, jpp), (nup, jp)): {'mu': dipole element, 'gam':  dephasing,
        'nu': position}}}.
    """
    def __init__(self, sys_params: dict, freq_shift: List[float]=None, filter=None):
        self.sys_params = sys_params
        if filter is None:
            self.filter = lambda x: True
        else:
            self.filter = filter
        self.freq_shift = freq_shift

    def cross_section_spectrum(self, kb: pw.KetBra, times: List[float]) -> np.complex:
        r"""Return cross section spectrum."""
        dt = times[-1][1]-times[-1][0]
        xs = self.cross_section(kb, times)
        if isinstance(xs, (list, tuple, np.ndarray)):
            xs_spectrum = fftp.fft(xs)*dt
        else:
            xs_spectrum = 0.0

        return xs_spectrum

    def cross_section_freqs(self, times: List[float]) -> np.real:
        r"""Return optical frequencies of the spectrum."""
        dt = times[-1][1]-times[-1][0]

        return fftp.fftfreq(len(times[-1]), d=dt)

    def leaf_term(self, nu, gam, times):
        return np.exp(-2.0*np.pi*times*(1.0j*nu+gam))


class MultiPropagator:
    """Execute :class:`Propagator` methods on multiple pathway trees.

    Parameters
    ----------
    pathways: list of :class:`pw.KetBra`
        Roots of pathways.
    """
    def __init__(self, pathways: list, prop: Propagator):
        self.pathways = pathways
        self.prop = prop

    def response(self, times: List[float]) -> np.complex:
        """Calculate the total response.

        Parameters
        ----------
        times: list of float
            Delays between pulses (and readout)

        Returns
        -------
        resp: complex
            Response function.
        """
        return sum(self.prop.response(kb, times) for kb in self.pathways)

    def cross_section(self, times: List[float]) -> np.complex:
        r"""Calculate interaction cross section.

        The returned quantity has units of `m^2 Hz`. Mulitply by `10e2/C.c` to
        get a HITRAN line intensity (cross section integrated over frequency
        range) in `cm^2 cm^{-1}`.
        """
        return sum(self.prop.cross_section(kb, times) for kb in self.pathways)

    def cross_section_spectrum(self, times: List[float]) -> np.complex:
        r"""Return cross section spectrum."""
        dt = times[-1][1]-times[-1][0]
        xs = self.cross_section(times)
        if isinstance(xs, (list, tuple, np.ndarray)):
            xs_spectrum = fftp.fft(xs)*dt
        else:
            xs_spectrum = 0.0

        return xs_spectrum


def times_dims(times: List[Union[float, np.ndarray]]) -> Tuple[int]:
    dims = []
    for d in times:
        try:
            dims.append(d.size)
        except AttributeError:
            pass

    return tuple(dims)


def power2electric_field(P: float, r: float):
    """Convert optical power (W) to electric field (N/C)."""
    return np.sqrt(P/np.pi/r**2/e2i)


def electric_field2power(E: float, r: float):
    """Convert electric field (N/C) to optical power (W)."""
    return np.abs(E)**2*np.pi*r**2*e2i


def electric_field(E0: float, xs: float, cd: float):
    """Return electric field generated by molecular sample.

    Assumes Dirac-delta electric field.

    Parameters
    ----------
    E0: float
        Source time-integrated electric field, V s/m.
    xs: float
        Cross section, m^2/s.
    cd: float
        Column density (density*length), molecules/m^2.

    Returns
    -------
    float
    """
    return 1.0j*E0*xs*cd


def polarization(E0: float, resp: float, n: float):
    """Return polarization density generated by molecular sample.

    Assumes Dirac-delta electric field.

    Parameters
    ----------
    E0: float
        Time-integrated electric field, N/C*s.
    resp: float
        [\mu]^2/J*Hz, result of :method:`Propagator.cross_section`.
    n: float
        Number density (density), molecules/m^3.

    Returns
    -------
    float
    """
    return resp*n*E0


def absorptive(spec2d):
    """Create absorptive spectrum from 2D reph/non-reph spectrum.

    Works for ket-only and full pathways.
    """
    return np.imag(spec2d - spec2d[:,::-1])


def dispersive(spec2d):
    """Create dispersive spectrum from 2D reph/non-reph spectrum.

    Works for ket-only and full pathways.
    """
    return np.real(spec2d + spec2d[:,::-1])
