import pytest
from pathlib import Path
from numbers import Complex
from sqlalchemy import create_engine
import rotsim2d.pathways as pw
import rotsim2d.dressedleaf as dl
from spectroscopy.alchemy.meta import hitran_cache
from spectroscopy.molecule import CH3ClAlchemyMode, DiatomState


@pytest.fixture
def pws():
    """Generate a list of pathways for further testing."""
    kbs = pw.gen_pathways([5], [0, 1, 2, 3], rotor='symmetric',
                          kiter_func=lambda x: [1], pump_overlap=False)
    pws = dl.Pathway.from_kb_list(kbs)

    return pws


def test_nineteen_geo_labels(pws):
    geo_labels = {pw.geo_label for pw in pws}

    assert len(geo_labels) == 19


def test_nineteen_trans_labels(pws):
    trans_labels = {pw.trans_label for pw in pws}

    assert len(trans_labels)


@pytest.fixture(scope='module')
def ch3cl_mode():
    sql_path = Path(hitran_cache) / 'CH3Cl.sqlite3'
    engine = create_engine("sqlite:///" + str(sql_path))
    ch3cl_mode = CH3ClAlchemyMode(engine)

    return ch3cl_mode


@pytest.fixture
def dressed_pws(ch3cl_mode):
    kbs = pw.gen_pathways([5], [0, 1, 2, 3], rotor='symmetric',
                          kiter_func=lambda x: [1], pump_overlap=False)
    dpws = dl.DressedPathway.from_kb_list(kbs, ch3cl_mode, 296.0)

    return dpws


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
