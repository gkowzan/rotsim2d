# pylint: disable=redefined-outer-name
from functools import partial
from pathlib import Path
from pprint import pprint

import pytest
import rotsim2d.couple as cp
import rotsim2d.dressedleaf as dl
import rotsim2d.pathways as pw
import rotsim2d.symbolic.functions as sym
from molspecutils.alchemy.meta import hitran_cache
from molspecutils.molecule import CH3ClAlchemyMode, DiatomState, SymTopState
from sqlalchemy import create_engine


@pytest.fixture
def partial_pws():
    return partial(pw.gen_pathways, [5], rotor='symmetric',
                   kiter_func=lambda x: [1], pump_overlap=False)

def test_onecolor(partial_pws):
    kbs = partial_pws(meths=[pw.only_dfwm])
    assert all(x.is_dfwm() for x in kbs[0].leaves)
    pws = dl.Pathway.from_kb_list(kbs)
    assert all(x.leaf.color_tier() == 1 for x in pws)


def test_twocolor(partial_pws):
    kbs = partial_pws(meths=[pw.only_twocolor])
    assert all(x.is_twocolor() for x in kbs[0].leaves)
    pws = dl.Pathway.from_kb_list(kbs)
    assert all(x.leaf.color_tier() == 2 for x in pws)


def test_threecolor(partial_pws):
    kbs = partial_pws(meths=[pw.only_threecolor])
    assert all(x.is_threecolor() for x in kbs[0].leaves)
    pws = dl.Pathway.from_kb_list(kbs)
    assert all(x.leaf.color_tier() == 3 for x in pws)


@pytest.fixture
def pws():
    """Generate a list of pathways for further testing."""
    kbs = pw.gen_pathways([5], rotor='symmetric',
                          kiter_func=lambda x: [1], pump_overlap=False)
    pws = dl.Pathway.from_kb_list(kbs)

    return pws


def test_pathway_equality(pws):
    assert pws[0] == pws[0]
    for pw2 in pws[1:]:
        assert pws[0] != pw2


def test_nineteen_geo_labels(pws):
    geo_labels = {pw.geo_label for pw in pws}

    assert len(geo_labels) == 19


def test_nineteen_trans_labels(pws):
    trans_labels = {pw.trans_label for pw in pws}

    assert len(trans_labels)


def pw_rme_sign(p):
    sides = [li.side for li in p.leaf.interactions()]
    sign = 1.0
    for pair, side in zip(p.transitions, sides):
        if side == pw.Side.BRA:
            pair = pair[::-1]
        if pair[0].j > pair[1].j:
            sign *= -1

    return sign


def test_pathway_sign(pws):
    """Check if pathway has the correct sign."""
    for p in pws:
        rfactor = p.geometric_factor()
        brute_rfactor = cp.brute_geometric_factor(p.js, 0)
        assert pytest.approx(rfactor) == brute_rfactor
        rme_sign = pw_rme_sign(p)
        brute_rme_sign = cp.brute_rme(p.js, 1)
        assert rme_sign*brute_rme_sign > 0 # same sign
        assert rfactor*rme_sign > 0
        assert brute_rfactor*brute_rme_sign > 0

@pytest.fixture(scope='module')
def ch3cl_mode():
    ch3cl_mode = CH3ClAlchemyMode()

    return ch3cl_mode

@pytest.fixture
def dressed_pws(ch3cl_mode):
    kbs = pw.gen_pathways([5], rotor='symmetric',
                          kiter_func=lambda x: [1], pump_overlap=False)
    dpws = dl.DressedPathway.from_kb_list(kbs, ch3cl_mode, 296.0)

    return dpws

def test_dressedpathway_rme_sign(dressed_pws):
    for p in dressed_pws:
        rme_sign = pw_rme_sign(p)
        rme = 1.0
        sides = [li.side for li in p.leaf.interactions()]
        for pair, side in zip(p.transitions, sides):
            if side == pw.Side.BRA:
                pair = pair[::-1]
            rme *= p.vib_mode.mu(pair)

        assert rme_sign*rme > 0

def test_dressedpathway_rme_sign(dressed_pws):
    for p in dressed_pws:
        rme_sign = pw_rme_sign(p)
        rme = 1.0
        sides = [li.side for li in p.leaf.interactions()]
        for pair, side in zip(p.transitions, sides):
            if side == pw.Side.BRA:
                pair = pair[::-1]
            rme *= p.vib_mode.mu(pair)

        assert rme_sign*rme > 0

def test_dressedpathway_equality(dressed_pws):
    assert dressed_pws[0] == dressed_pws[0]
    for pw2 in dressed_pws[1:]:
        assert dressed_pws[0] != pw2

def test_gammas_are_positive(dressed_pws):
    for pw in dressed_pws:
        for i in range(3):
            assert pw.gamma(i) >= 0.0

def test_intensities_are_nonzero(dressed_pws):
    for pw in dressed_pws:
        assert abs(pw.intensity()) >= 0.0


def test_coherence_frequencies_are_real(dressed_pws):
    for pw in dressed_pws:
        for i in range(3):
            assert pw.nu(i).imag == 0

def test_missing_level_handling(ch3cl_mode):
    kbs = pw.gen_pathways([20], rotor='symmetric',
                          kiter_func=lambda x: [20])
    dpws = dl.DressedPathway.from_kb_list(kbs, ch3cl_mode, 296.0)

@pytest.fixture
def peak_list(dressed_pws):
    return dl.peak_list(dressed_pws)


def test_peak_list_frequencies_are_real(peak_list):
    assert all(pump.imag == 0 for pump in peak_list.pumps)
    assert all(probe.imag == 0 for probe in peak_list.probes)


def test_peak_list_contains_dfwm(dressed_pws):
    sl = [(DiatomState(nu=0, j=0), DiatomState(nu=0, j=0)),
          (DiatomState(nu=1, j=-1), DiatomState(nu=0, j=0)),
          (DiatomState(nu=0, j=0), DiatomState(nu=0, j=0)),
          (DiatomState(nu=1, j=-1), DiatomState(nu=0, j=0)),
          (DiatomState(nu=0, j=0), DiatomState(nu=0, j=0))]
    pl, dll = dl.peak_list(dressed_pws, return_dls=True)

    assert len(dl.equiv_peaks(sl, pl, dll)) > 0
