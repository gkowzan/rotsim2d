from rotsim2d.dressedleaf import geometric_labels
import rotsim2d.symbolic.results as symr
import rotsim2d.symbolic.functions as sym

def test_gfactors():
    for js, label in geometric_labels.items():
        args = [symr.J_i+j for j in js]
        gfacs = tuple(
            (symr.factor(symr.powdenest(sym.gfactor_expr(*(args + [i])), force=True), deep=True)
             for i in range(3)))
        for i in range(3):
            assert symr.factor(symr.powdenest(symr.nsimplify(
                gfacs[i]-symr.gfactors[label][i]), force=True), deep=True) == symr.S(0)


def test_gfactors_highj():
    gfactors_highj = {k: tuple(symr.nsimplify(symr.limit(x*(2*symr.J_i+1)**(3/2), symr.J_i, symr.oo)) for x in v)
                      for k, v in symr.gfactors.items()}
    for k in gfactors_highj.keys():
        for i in range(3):
            assert symr.factor(symr.nsimplify(
                symr.gfactors_highj[k][i]-gfactors_highj[k][i]), deep=True) == symr.S(0)

def test_rfactors():
    rfactors = {k: sym.rfactorize(symr.gfactors[k], symr.T00_exprs)
                for k in symr.gfactors.keys()}
    for k in rfactors.keys():
        assert symr.factor(symr.nsimplify(
            symr.rfactors[k]-rfactors[k]), deep=True) == symr.S(0)

def test_rfactors_dict():
    rfactors_dict = {k: sym.rfactorize(symr.gfactors[k], symr.T00_exprs, cfac=False, coeffs=True)
                     for k in symr.gfactors.keys()}
    for label in rfactors_dict.keys():
        for coeff_label in rfactors_dict[label].keys():
            assert symr.factor(symr.nsimplify(
                symr.rfactors_dict[label][coeff_label]-rfactors_dict[label][coeff_label]
            ), deep=True) == symr.S(0)

def test_rfactors_xxxx():
    rfactors_xxxx = {k: symr.factor(symr.powdenest(
        v.subs({symr.phi: symr.S(0), symr.phj: symr.S(0), symr.phk: symr.S(0), symr.phl: symr.S(0)}), force=True), deep=True)
                     for k, v in symr.rfactors.items()}
    for k in rfactors_xxxx.keys():
        assert symr.factor(symr.nsimplify(
            symr.rfactors_xxxx[k]-rfactors_xxxx[k]), deep=True) == symr.S(0)

def test_rfactors_relative():
    for k in symr.gfactors:
        actual = sym.rfactorize(symr.gfactors[k], symr.T00_exprs, cfac=False, relative=True)
        assert symr.factor(symr.nsimplify(actual-symr.rfactors_relative[k]), deep=True) == symr.S(0)

def test_rfactors_highj():
    for k in symr.gfactors_highj:
        actual = sym.rfactorize(symr.gfactors_highj[k], symr.T00_exprs, cfac=False)
        assert symr.factor(symr.nsimplify(actual-symr.rfactors_highj[k]), deep=True) == symr.S(0)

def test_rfactors_highj_dict():
    for k in symr.gfactors_highj:
        expected = sym.rfactorize(symr.gfactors_highj[k], symr.T00_exprs, cfac=False, coeffs=True)
        for kk in symr.rfactors_highj_dict[k]:
            assert expected[kk]-symr.rfactors_highj_dict[k][kk] == symr.S(0)


def test_vaccaro_angles():
    expected = [-sym.atan(sym.Rational(3, 4)), sym.atan(sym.Rational(1, 2)),
                -sym.atan(2), -sym.atan(sym.Rational(1, 3))]
    actual = sym.detection_angles([sym.pi/2, sym.pi/4, sym.pi/2])

    assert expected == list(actual.keys())


def test_gk_angles():
    expected = [sym.atan(sym.Rational(3, 2)), sym.atan(4), -sym.pi/4,
                sym.atan(sym.Rational(1, 4)), sym.atan(sym.Rational(2, 3))]
    actual = sym.detection_angles([0, sym.pi/4, sym.pi/2])

    assert expected == list(actual.keys())
