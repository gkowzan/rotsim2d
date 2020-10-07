"""Propagate time evolution of coherences based on pathway trees.

General idea is to dress up the tree with dipole couplings, pulse properties,
dephasing and evolution times and then collapse the tree into a scalar.

TODO:
- add selection of k
- generate 2D spectra
"""
from typing import List
import numpy as np
import scipy.constants as C
import pyfftw.interfaces.scipy_fftpack as fftp
from rotsim2d import pathways as pw
from spectroscopy import happier
import shed.units as u

e2i = C.epsilon_0*C.c/2         #: Electric field to intensity
xs2cm = 1e2/C.c                 #: Cross section in m^2 Hz to cm^2 cm^{-1}


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

    
class Spectrum(CrossSectionMixin):
    """Calculate spectrum based on excitation tree.

    Parameters
    ----------
    sys_params: dict of dict
        {'elevels': {(nu, j): elevel}, 'populations': {(nu, j): pop},
        'line_params': {((nupp, jpp), (nup, jp)): {'mu': dipole element, 'gam':  dephasing,
        'nu': position}}}.
    """
    def __init__(self, sys_params: dict, freq_shift: float=0.0, filter=None):
        self.sys_params = sys_params
        if filter is None:
            self.filter = lambda x: True
        else:
            self.filter = filter
        self.freq_shift = freq_shift

    def leaf_response(self, leaf: pw.KetBra, freqs: List[float]) -> np.complex:
        r"""Return single pathway response spectrum."""
        resp = np.complex(1.0)
        kb_series = [x for x in leaf.ancestors if isinstance(x, pw.KetBra)] + [leaf]
        for i in range(1, len(kb_series)):
            kb, kbp = kb_series[i], kb_series[i-1]

            # dipole matrix element
            if kb.parent.side is pw.Side.KET:
                pair = ((kbp.knu, kbp.kj), (kb.knu, kb.kj))
            else:
                pair = ((kbp.bnu, kbp.bj), (kb.bnu, kb.bj))
            try:
                mu = self.sys_params['line_params'][pair]['mu']
            except KeyError:
                mu = self.sys_params['line_params'][(pair[1], pair[0])]['mu']

            if kb.parent.readout:
                resp *= mu*kb.root.pop
                break

            # nu and gam
            # pairs in line_params have lower nu first
            pair = ((kb.bnu, kb.bj), (kb.knu, kb.kj))
            try:
                nu = self.sys_params['line_params'][pair]['nu'] - self.freq_shift
                gam = self.sys_params['line_params'][pair]['gam']
            except KeyError:
                try:
                    nu = -self.sys_params['line_params'][(pair[1], pair[0])]['nu'] + self.freq_shift
                    gam = self.sys_params['line_params'][(pair[1], pair[0])]['gam']
                except KeyError:
                    resp *= 1.0j/C.hbar*kb.parent.side*mu
                    continue

            resp = resp*1.0j/C.hbar*kb.parent.side*mu/(freqs[i-1]-nu - 1.0j*gam)

        return resp

            
class Propagator(CrossSectionMixin):
    """Calculate useful physical quantities based on excitation tree.

    Parameters
    ----------
    sys_params: dict of dict
        {'elevels': {(nu, j): elevel}, 'populations': {(nu, j): pop},
        'line_params': {((nupp, jpp), (nup, jp)): {'mu': dipole element, 'gam':  dephasing,
        'nu': position}}}.
    """
    def __init__(self, sys_params: dict, freq_shift: float=0.0, filter=None):
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
            xs_spectrum = fftp.fft(xs, planner_effort='FFTW_ESTIMATE')*dt
        else:
            xs_spectrum = 0.0

        return xs_spectrum

    def cross_section_freqs(self, times: List[float]) -> np.real:
        r"""Return optical frequencies of the spectrum."""
        dt = times[-1][1]-times[-1][0]

        return fftp.fftfreq(len(times[-1]), d=dt)

    # @profile
    def leaf_response(self, leaf: pw.KetBra, times: List[float]) -> np.complex:
        r"""Calculate response from a single leaf.

        The returned quantity has units of :math:`[\mu]^n[\hbar]^{-n}` for
        `n`-th order excitation or of :math:`[\mu]^{n+1}[\hbar]^{-n}` for `n`-th
        order excitaion followed by readout. Multiplication by electric field
        gives n-th order excited fraction in the first case, and n-th order
        polarization in the second case.
        """
        resp = np.complex(1.0)
        kb_series = [x for x in leaf.ancestors if isinstance(x, pw.KetBra)] + [leaf]
        for i in range(1, len(kb_series)):
            kb, kbp = kb_series[i], kb_series[i-1]

            # dipole matrix element
            if kb.parent.side is pw.Side.KET:
                pair = ((kbp.knu, kbp.kj), (kb.knu, kb.kj))
            else:
                pair = ((kbp.bnu, kbp.bj), (kb.bnu, kb.bj))
            try:
                mu = self.sys_params['line_params'][pair]['mu']
            except KeyError:
                mu = self.sys_params['line_params'][(pair[1], pair[0])]['mu']

            if kb.parent.readout:
                resp *= mu*kb.root.pop
                break

            # nu and gam
            # pairs in line_params have lower nu first
            pair = ((kb.bnu, kb.bj), (kb.knu, kb.kj))
            try:
                nu = self.sys_params['line_params'][pair]['nu'] - self.freq_shift
                gam = self.sys_params['line_params'][pair]['gam']
            except KeyError:
                try:
                    nu = -self.sys_params['line_params'][(pair[1], pair[0])]['nu'] + self.freq_shift
                    gam = self.sys_params['line_params'][(pair[1], pair[0])]['gam']
                except KeyError:
                    nu = 0.0
                    gam = 0.0

            resp *= 1.0j/C.hbar*kb.parent.side*mu*np.exp(-2.0*np.pi*times[i-1]*(1.0j*nu+gam))

        return resp


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
            xs_spectrum = fftp.fft(xs, planner_effort='FFTW_ESTIMATE')*dt
        else:
            xs_spectrum = 0.0

        return xs_spectrum


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
    """Create absorptive spectrum from 2D reph/non-reph spectrum."""
    pass
