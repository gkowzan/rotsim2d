"""Derived symbolic expressions for angular momentum dependence and polarization
dependence of four-fold dipole interaction operator.

All the contents can obtained by using :mod:`rotsim2d.symbolic.functions` functions.  See
:download:`examples/generate_angular.symbolic_results.py <../../examples/generate_angular.symbolic_results.py>`.
"""
# * Imports
from rotsim2d.symbolic.common import *
from molspecutils.molecule import DiatomState
from numpy import array

# * G-factors
# ** General expressions
#: Isotropic version of angular momentum-dependent combinations of three
#: Wigner-6j factors involved in evaluating four-fold dipole interaction
#: operator.
gfactors = {
    'PPR': (0, 0, sqrt(5) / (5 * sqrt(2 * J_i + 1) * Abs(2 * J_i - 1))),
    'RRP': (0, 0, sqrt(5) / (5 * sqrt(2 * J_i + 1) * (2 * J_i + 3))),
    'PRP': (1 / (3 * (2 * J_i + 1)**Rational(3, 2)),
            -sqrt(3) * (J_i + 1) / (6 * J_i * (2 * J_i + 1)**Rational(3, 2)),
            sqrt(5) * (J_i + 1) * (2 * J_i + 3) /
            (30 * J_i * (2 * J_i - 1) * (2 * J_i + 1)**Rational(3, 2))),
    'RPR': (1 / (3 * (2 * J_i + 1)**Rational(3, 2)),
            -sqrt(3) * J_i / (6 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
            sqrt(5) * J_i * (2 * J_i - 1) /
            (30 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2) * (2 * J_i + 3))),
    'RPP': (1 / (3 * (2 * J_i + 1)**Rational(3, 2)),
            sqrt(3) / (6 * (2 * J_i + 1)**Rational(3, 2)),
            sqrt(5) / (30 * (2 * J_i + 1)**Rational(3, 2))),
    'PRR': (1 / (3 * (2 * J_i + 1)**Rational(3, 2)),
            sqrt(3) / (6 * (2 * J_i + 1)**Rational(3, 2)),
            sqrt(5) / (30 * (2 * J_i + 1)**Rational(3, 2))),
    'PQQ':
    (0, sqrt(3) * (J_i - 1) / (6 * J_i * (2 * J_i - 1) * sqrt(2 * J_i + 1)),
     -sqrt(5) * (J_i + 1) / (10 * J_i * sqrt(2 * J_i + 1) * Abs(2 * J_i - 1))),
    'QPQ':
    (0, -sqrt(3) * sqrt(J_i - 1) * sqrt(J_i + 1) /
     (6 * J_i * sqrt(2 * J_i - 1) * (2 * J_i + 1)), -sqrt(5) * sqrt(J_i - 1) *
     sqrt(J_i + 1) / (10 * J_i * sqrt(2 * J_i - 1) * (2 * J_i + 1))),
    'QPR':
    (0, sqrt(3) * (J_i + 1) / (6 * J_i * (2 * J_i + 1)**Rational(3, 2)),
     -sqrt(5) * Abs(J_i - 1) / (10 * J_i * (2 * J_i + 1)**Rational(3, 2))),
    'QQP':
    (-1 / (3 * (2 * J_i + 1)**Rational(3, 2)),
     sqrt(3) / (6 * J_i * (2 * J_i + 1)**Rational(3, 2)),
     sqrt(5) * (2 * J_i + 3) / (30 * J_i * (2 * J_i + 1)**Rational(3, 2))),
    'QQQ': (1 / (3 * (2 * J_i + 1)**Rational(3, 2)),
            -sqrt(3) / (6 * J_i * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
            sqrt(5) * (2 * J_i - 1) * (2 * J_i + 3) /
            (30 * J_i * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2))),
    'QRP':
    (0, sqrt(3) * J_i / (6 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
     -sqrt(5) * (J_i + 2) / (10 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2))),
    'RPQ': (-1 / (3 * (2 * J_i + 1)**Rational(3, 2)),
            -sqrt(3) / (6 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
            sqrt(5) * Abs(2 * J_i - 1) / (30 * (J_i + 1) *
                                          (2 * J_i + 1)**Rational(3, 2))),
    'RQP': (0, -sqrt(3) * sqrt(J_i) * sqrt(J_i + 2) /
            (6 * (J_i + 1) * (2 * J_i + 1) * sqrt(2 * J_i + 3)),
            -sqrt(5) * sqrt(J_i) * sqrt(J_i + 2) /
            (10 * (J_i + 1) * (2 * J_i + 1) * sqrt(2 * J_i + 3))),
    'RQQ':
    (0,
     sqrt(3) * (J_i + 2) / (6 * (J_i + 1) * sqrt(2 * J_i + 1) * (2 * J_i + 3)),
     -sqrt(5) * J_i / (10 * (J_i + 1) * sqrt(2 * J_i + 1) * (2 * J_i + 3))),
    'PQR': (0, -sqrt(3) * sqrt((J_i - 1) * (2 * J_i - 1)) * sqrt(J_i + 1) /
            (6 * J_i * (2 * J_i - 1) * (2 * J_i + 1)),
            -sqrt(5) * sqrt(J_i - 1) * sqrt(J_i + 1) /
            (10 * J_i * sqrt(2 * J_i - 1) * (2 * J_i + 1))),
    'PRQ': (-1 / (3 * (2 * J_i + 1)**Rational(3, 2)),
            sqrt(3) / (6 * J_i * (2 * J_i + 1)**Rational(3, 2)), sqrt(5) *
            (2 * J_i + 3) / (30 * J_i * (2 * J_i + 1)**Rational(3, 2))),
    'QQR': (-1 / (3 * (2 * J_i + 1)**Rational(3, 2)),
            -sqrt(3) / (6 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
            sqrt(5) * Abs(2 * J_i - 1) / (30 * (J_i + 1) *
                                          (2 * J_i + 1)**Rational(3, 2))),
    'QRQ': (0, -sqrt(3) * sqrt(J_i) * sqrt(J_i + 2) /
            (6 * (J_i + 1) * (2 * J_i + 1) * sqrt(2 * J_i + 3)),
            -sqrt(5) * sqrt(J_i) * sqrt(J_i + 2) /
            (10 * (J_i + 1) * (2 * J_i + 1) * sqrt(2 * J_i + 3)))
}

# ** High-J limit
# gfactors_highj = {k: tuple(nsimplify(limit(x*(2*J_i+1)**(3/2), J_i, oo)) for x in v)
#                   for k, v in gfactors.items()}
#: G-factors multiplied by :math:`(2J_i+1)^{3/2}` and with :math:`J_i\to\infty`.
gfactors_highj = {
    'PPR': (0, 0, sqrt(5) / 5),
    'RRP': (0, 0, sqrt(5) / 5),
    'PRP': (Rational(1, 3), -sqrt(3) / 6, sqrt(5) / 30),
    'RPR': (Rational(1, 3), -sqrt(3) / 6, sqrt(5) / 30),
    'RPP': (Rational(1, 3), sqrt(3) / 6, sqrt(5) / 30),
    'PRR': (Rational(1, 3), sqrt(3) / 6, sqrt(5) / 30),
    'PQQ': (0, sqrt(3) / 6, -sqrt(5) / 10),
    'QPQ': (0, -sqrt(3) / 6, -sqrt(5) / 10),
    'QPR': (0, sqrt(3) / 6, -sqrt(5) / 10),
    'QQP': (Rational(-1, 3), 0, sqrt(5) / 15),
    'QQQ': (Rational(1, 3), 0, 2 * sqrt(5) / 15),
    'QRP': (0, sqrt(3) / 6, -sqrt(5) / 10),
    'RPQ': (Rational(-1, 3), 0, sqrt(5) / 15),
    'RQP': (0, -sqrt(3) / 6, -sqrt(5) / 10),
    'RQQ': (0, sqrt(3) / 6, -sqrt(5) / 10),
    'PQR': (0, -sqrt(3) / 6, -sqrt(5) / 10),
    'PRQ': (Rational(-1, 3), 0, sqrt(5) / 15),
    'QQR': (Rational(-1, 3), 0, sqrt(5) / 15),
    'QRQ': (0, -sqrt(3) / 6, -sqrt(5) / 10)
}

#: Numerical values of high-J G-factors.
gfactors_highj_numeric = {'PPR': array([ 0.       ,  0.       ,  0.4472136]),
 'PQQ': array([ 0.        ,  0.28867513, -0.2236068 ]),
 'PQR': array([ 0.        , -0.28867513, -0.2236068 ]),
 'PRP': array([ 0.33333333, -0.28867513,  0.0745356 ]),
 'PRQ': array([-0.33333333,  0.        ,  0.1490712 ]),
 'PRR': array([ 0.33333333,  0.28867513,  0.0745356 ]),
 'QPQ': array([ 0.        , -0.28867513, -0.2236068 ]),
 'QPR': array([ 0.        ,  0.28867513, -0.2236068 ]),
 'QQP': array([-0.33333333,  0.        ,  0.1490712 ]),
 'QQQ': array([ 0.33333333,  0.        ,  0.2981424 ]),
 'QQR': array([-0.33333333,  0.        ,  0.1490712 ]),
 'QRP': array([ 0.        ,  0.28867513, -0.2236068 ]),
 'QRQ': array([ 0.        , -0.28867513, -0.2236068 ]),
 'RPP': array([ 0.33333333,  0.28867513,  0.0745356 ]),
 'RPQ': array([-0.33333333,  0.        ,  0.1490712 ]),
 'RPR': array([ 0.33333333, -0.28867513,  0.0745356 ]),
 'RQP': array([ 0.        , -0.28867513, -0.2236068 ]),
 'RQQ': array([ 0.        ,  0.28867513, -0.2236068 ]),
 'RRP': array([ 0.       ,  0.       ,  0.4472136])}

# * Linear polarization
#: Expressions for isotropic four-fold polarization tensor components (dummy angles)
T00_exprs = [
    cos(phi - phj)*cos(phk - phl)/3,
    sqrt(3)*sin(phi - phj)*sin(phk - phl)/6,
    sqrt(5)*(cos(phi - phj - phk + phl) +
             cos(phi - phj + phk - phl) +
             6*cos(phi + phj - phk - phl))/60
]

#: Expressions for isotropic four-fold polarization tensor components (experimental angles)
T00_theta_exprs = [
    cos(theta_i - theta_j)*cos(theta_k - theta_l)/3,
    sqrt(3)*sin(theta_i - theta_j)*sin(theta_k - theta_l)/6,
    sqrt(5)*(cos(theta_i - theta_j - theta_k + theta_l) +
             cos(theta_i - theta_j + theta_k - theta_l) +
             6*cos(theta_i + theta_j - theta_k - theta_l))/60
]
# * R-factors
# ** General expressions
#: R-factors as dictionaries of coefficients
rfactors_dict = {
    'PRR': {
        'c0': 1 / (60 * (2 * J_i - 1) * sqrt(2 * J_i + 1)),
        'c12': 6,
        'c13': 1,
        'c14': 1
    },
    'PR(P)': {
        'c0': 1 / (60 * (2 * J_i - 1) * sqrt(2 * J_i + 1)),
        'c12': 6,
        'c13': 1,
        'c14': 1
    },
    'PQQ': {
        'c0': -1 / (60 * J_i * (2 * J_i - 1) * sqrt(2 * J_i + 1)),
        'c12': 3 * J_i + 3,
        'c13': 3 * J_i - 2,
        'c14': 3 - 2 * J_i
    },
    'PQR': {
        'c0':
        -sqrt(J_i - 1) * sqrt(J_i + 1) / (60 * J_i * sqrt(2 * J_i - 1) *
                                          (2 * J_i + 1)),
        'c12':
        3,
        'c13':
        -2,
        'c14':
        3
    },
    'PQ(P)': {
        'c0': -1 / (60 * J_i * (2 * J_i - 1) * sqrt(2 * J_i + 1)),
        'c12': 3 * J_i + 3,
        'c13': 3 - 2 * J_i,
        'c14': 3 * J_i - 2
    },
    'PQ(Q)': {
        'c0':
        -sqrt(J_i - 1) * sqrt(J_i + 1) / (60 * J_i * sqrt(2 * J_i - 1) *
                                          (2 * J_i + 1)),
        'c12':
        3,
        'c13':
        3,
        'c14':
        -2
    },
    'PPP': {
        'c0': 1 / (60 * J_i * (2 * J_i - 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i**2 + 5 * J_i + 3,
        'c13': 12 * J_i**2 - 2,
        'c14': 2 * J_i**2 - 5 * J_i + 3
    },
    'PPQ': {
        'c0': 1 / (60 * J_i * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i + 3,
        'c13': -3 * J_i - 2,
        'c14': 3 - 3 * J_i
    },
    'PPR': {
        'c0': 1 / (60 * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 1,
        'c13': 1,
        'c14': 6
    },
    'PP(P)': {
        'c0': 1 / (60 * J_i * (2 * J_i - 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i**2 + 5 * J_i + 3,
        'c13': 2 * J_i**2 - 5 * J_i + 3,
        'c14': 12 * J_i**2 - 2
    },
    'PP(Q)': {
        'c0': 1 / (60 * J_i * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i + 3,
        'c13': 3 - 3 * J_i,
        'c14': -3 * J_i - 2
    },
    'PP(R)': {
        'c0': 1 / (60 * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 1,
        'c13': 6,
        'c14': 1
    },
    'P2P2P': {
        'c0': 1 / (60 * (2 * J_i - 1) * sqrt(2 * J_i + 1)),
        'c12': 6,
        'c13': 1,
        'c14': 1
    },
    'P2P(P)': {
        'c0': 1 / (60 * (2 * J_i - 1) * sqrt(2 * J_i + 1)),
        'c12': 6,
        'c13': 1,
        'c14': 1
    },
    'P2Q2Q': {
        'c0': -1 / (60 * J_i * (2 * J_i - 1) * sqrt(2 * J_i + 1)),
        'c12': 3 * J_i + 3,
        'c13': 3 * J_i - 2,
        'c14': 3 - 2 * J_i
    },
    'P2Q2P': {
        'c0':
        -sqrt(J_i - 1) * sqrt(J_i + 1) / (60 * J_i * sqrt(2 * J_i - 1) *
                                          (2 * J_i + 1)),
        'c12':
        3,
        'c13':
        -2,
        'c14':
        3
    },
    'P2Q(P)': {
        'c0': -1 / (60 * J_i * (2 * J_i - 1) * sqrt(2 * J_i + 1)),
        'c12': 3 * J_i + 3,
        'c13': 3 - 2 * J_i,
        'c14': 3 * J_i - 2
    },
    'P2Q(Q)': {
        'c0':
        -sqrt(J_i - 1) * sqrt(J_i + 1) / (60 * J_i * sqrt(2 * J_i - 1) *
                                          (2 * J_i + 1)),
        'c12':
        3,
        'c13':
        3,
        'c14':
        -2
    },
    'P2R2R': {
        'c0': 1 / (60 * J_i * (2 * J_i - 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i**2 + 5 * J_i + 3,
        'c13': 12 * J_i**2 - 2,
        'c14': 2 * J_i**2 - 5 * J_i + 3
    },
    'P2R2Q': {
        'c0': 1 / (60 * J_i * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i + 3,
        'c13': -3 * J_i - 2,
        'c14': 3 - 3 * J_i
    },
    'P2R2P': {
        'c0': 1 / (60 * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 1,
        'c13': 1,
        'c14': 6
    },
    'P2R(P)': {
        'c0': 1 / (60 * J_i * (2 * J_i - 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i**2 + 5 * J_i + 3,
        'c13': 2 * J_i**2 - 5 * J_i + 3,
        'c14': 12 * J_i**2 - 2
    },
    'P2R(Q)': {
        'c0': 1 / (60 * J_i * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i + 3,
        'c13': 3 - 3 * J_i,
        'c14': -3 * J_i - 2
    },
    'P2R(R)': {
        'c0': 1 / (60 * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 1,
        'c13': 6,
        'c14': 1
    },
    'P(P)R': {
        'c0': 1 / (60 * (2 * J_i - 1) * sqrt(2 * J_i + 1)),
        'c12': 1,
        'c13': 6,
        'c14': 1
    },
    'P(P)Q': {
        'c0': 1 / (60 * J_i * (2 * J_i - 1) * sqrt(2 * J_i + 1)),
        'c12': 2 * J_i - 3,
        'c13': -3 * J_i - 3,
        'c14': 2 - 3 * J_i
    },
    'P(P)P': {
        'c0': 1 / (60 * J_i * (2 * J_i - 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i**2 - 5 * J_i + 3,
        'c13': 2 * J_i**2 + 5 * J_i + 3,
        'c14': 12 * J_i**2 - 2
    },
    'P(P)2P': {
        'c0': 1 / (60 * (2 * J_i - 1) * sqrt(2 * J_i + 1)),
        'c12': 1,
        'c13': 6,
        'c14': 1
    },
    'P(P)2Q': {
        'c0': 1 / (60 * J_i * (2 * J_i - 1) * sqrt(2 * J_i + 1)),
        'c12': 2 * J_i - 3,
        'c13': -3 * J_i - 3,
        'c14': 2 - 3 * J_i
    },
    'P(P)2R': {
        'c0': 1 / (60 * J_i * (2 * J_i - 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i**2 - 5 * J_i + 3,
        'c13': 2 * J_i**2 + 5 * J_i + 3,
        'c14': 12 * J_i**2 - 2
    },
    'P(P)(R)': {
        'c0': 1 / (60 * (2 * J_i - 1) * sqrt(2 * J_i + 1)),
        'c12': 1,
        'c13': 1,
        'c14': 6
    },
    'P(P)(Q)': {
        'c0': 1 / (60 * J_i * (2 * J_i - 1) * sqrt(2 * J_i + 1)),
        'c12': 2 * J_i - 3,
        'c13': 2 - 3 * J_i,
        'c14': -3 * J_i - 3
    },
    'P(P)(P)': {
        'c0': 1 / (60 * J_i * (2 * J_i - 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i**2 - 5 * J_i + 3,
        'c13': 12 * J_i**2 - 2,
        'c14': 2 * J_i**2 + 5 * J_i + 3
    },
    'P(P)(2P)': {
        'c0': 1 / (60 * (2 * J_i - 1) * sqrt(2 * J_i + 1)),
        'c12': 1,
        'c13': 1,
        'c14': 6
    },
    'P(P)(2Q)': {
        'c0': 1 / (60 * J_i * (2 * J_i - 1) * sqrt(2 * J_i + 1)),
        'c12': 2 * J_i - 3,
        'c13': 2 - 3 * J_i,
        'c14': -3 * J_i - 3
    },
    'P(P)(2R)': {
        'c0': 1 / (60 * J_i * (2 * J_i - 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i**2 - 5 * J_i + 3,
        'c13': 12 * J_i**2 - 2,
        'c14': 2 * J_i**2 + 5 * J_i + 3
    },
    'P(Q)Q': {
        'c0':
        -sqrt(J_i - 1) * sqrt(J_i + 1) / (60 * J_i * sqrt(2 * J_i - 1) *
                                          (2 * J_i + 1)),
        'c12':
        3,
        'c13':
        3,
        'c14':
        -2
    },
    'P(Q)P': {
        'c0': -1 / (60 * J_i * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 3 * J_i - 3,
        'c13': -2 * J_i - 3,
        'c14': 3 * J_i + 2
    },
    'P(Q)2Q': {
        'c0':
        -sqrt(J_i - 1) * sqrt(J_i + 1) / (60 * J_i * sqrt(2 * J_i - 1) *
                                          (2 * J_i + 1)),
        'c12':
        3,
        'c13':
        3,
        'c14':
        -2
    },
    'P(Q)2R': {
        'c0': -1 / (60 * J_i * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 3 * J_i - 3,
        'c13': -2 * J_i - 3,
        'c14': 3 * J_i + 2
    },
    'P(Q)(R)': {
        'c0':
        -sqrt(J_i - 1) * sqrt(J_i + 1) / (60 * J_i * sqrt(2 * J_i - 1) *
                                          (2 * J_i + 1)),
        'c12':
        3,
        'c13':
        -2,
        'c14':
        3
    },
    'P(Q)(Q)': {
        'c0': -1 / (60 * J_i * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 3 * J_i - 3,
        'c13': 3 * J_i + 2,
        'c14': -2 * J_i - 3
    },
    'P(Q)(2P)': {
        'c0':
        -sqrt(J_i - 1) * sqrt(J_i + 1) / (60 * J_i * sqrt(2 * J_i - 1) *
                                          (2 * J_i + 1)),
        'c12':
        3,
        'c13':
        -2,
        'c14':
        3
    },
    'P(Q)(2Q)': {
        'c0': -1 / (60 * J_i * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 3 * J_i - 3,
        'c13': 3 * J_i + 2,
        'c14': -2 * J_i - 3
    },
    'P(R)P': {
        'c0': 1 / (60 * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 6,
        'c13': 1,
        'c14': 1
    },
    'P(R)2R': {
        'c0': 1 / (60 * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 6,
        'c13': 1,
        'c14': 1
    },
    'P(R)(R)': {
        'c0': 1 / (60 * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 6,
        'c13': 1,
        'c14': 1
    },
    'P(R)(2P)': {
        'c0': 1 / (60 * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 6,
        'c13': 1,
        'c14': 1
    },
    'QRQ': {
        'c0':
        -sqrt(J_i - 1) * sqrt(J_i + 1) / (60 * J_i * sqrt(2 * J_i - 1) *
                                          (2 * J_i + 1)),
        'c12':
        3,
        'c13':
        -2,
        'c14':
        3
    },
    'QRR': {
        'c0': -1 / (60 * J_i * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 3 * J_i - 3,
        'c13': 3 * J_i + 2,
        'c14': -2 * J_i - 3
    },
    'QR(P)': {
        'c0':
        -sqrt(J_i - 1) * sqrt(J_i + 1) / (60 * J_i * sqrt(2 * J_i - 1) *
                                          (2 * J_i + 1)),
        'c12':
        3,
        'c13':
        3,
        'c14':
        -2
    },
    'QR(Q)': {
        'c0': -1 / (60 * J_i * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 3 * J_i - 3,
        'c13': -2 * J_i - 3,
        'c14': 3 * J_i + 2
    },
    'QQP': {
        'c0': 1 / (60 * J_i * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i + 3,
        'c13': -3 * J_i - 2,
        'c14': 3 - 3 * J_i
    },
    'QQQ': {
        'c0': 1 / (60 * J_i * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 4 * J_i**2 + 4 * J_i - 3,
        'c13': 4 * J_i**2 + 4 * J_i + 2,
        'c14': 4 * J_i**2 + 4 * J_i - 3
    },
    'QQR': {
        'c0': 1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i - 1,
        'c13': -3 * J_i - 1,
        'c14': -3 * J_i - 6
    },
    'QQ(P)': {
        'c0': 1 / (60 * J_i * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i + 3,
        'c13': 3 - 3 * J_i,
        'c14': -3 * J_i - 2
    },
    'QQ(Q)': {
        'c0': 1 / (60 * J_i * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 4 * J_i**2 + 4 * J_i - 3,
        'c13': 4 * J_i**2 + 4 * J_i - 3,
        'c14': 4 * J_i**2 + 4 * J_i + 2
    },
    'QQ(R)': {
        'c0': 1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i - 1,
        'c13': -3 * J_i - 6,
        'c14': -3 * J_i - 1
    },
    'QPP': {
        'c0': -1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 3 * J_i + 6,
        'c13': 3 * J_i + 1,
        'c14': 1 - 2 * J_i
    },
    'QPQ': {
        'c0':
        -sqrt(J_i) * sqrt(J_i + 2) / (60 * (J_i + 1) *
                                      (2 * J_i + 1) * sqrt(2 * J_i + 3)),
        'c12':
        3,
        'c13':
        -2,
        'c14':
        3
    },
    'QP(Q)': {
        'c0': -1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 3 * J_i + 6,
        'c13': 1 - 2 * J_i,
        'c14': 3 * J_i + 1
    },
    'QP(R)': {
        'c0':
        -sqrt(J_i) * sqrt(J_i + 2) / (60 * (J_i + 1) *
                                      (2 * J_i + 1) * sqrt(2 * J_i + 3)),
        'c12':
        3,
        'c13':
        3,
        'c14':
        -2
    },
    'Q2P2Q': {
        'c0':
        -sqrt(J_i - 1) * sqrt(J_i + 1) / (60 * J_i * sqrt(2 * J_i - 1) *
                                          (2 * J_i + 1)),
        'c12':
        3,
        'c13':
        -2,
        'c14':
        3
    },
    'Q2P2P': {
        'c0': -1 / (60 * J_i * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 3 * J_i - 3,
        'c13': 3 * J_i + 2,
        'c14': -2 * J_i - 3
    },
    'Q2P(P)': {
        'c0':
        -sqrt(J_i - 1) * sqrt(J_i + 1) / (60 * J_i * sqrt(2 * J_i - 1) *
                                          (2 * J_i + 1)),
        'c12':
        3,
        'c13':
        3,
        'c14':
        -2
    },
    'Q2P(Q)': {
        'c0': -1 / (60 * J_i * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 3 * J_i - 3,
        'c13': -2 * J_i - 3,
        'c14': 3 * J_i + 2
    },
    'Q2Q2R': {
        'c0': 1 / (60 * J_i * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i + 3,
        'c13': -3 * J_i - 2,
        'c14': 3 - 3 * J_i
    },
    'Q2Q2Q': {
        'c0': 1 / (60 * J_i * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 4 * J_i**2 + 4 * J_i - 3,
        'c13': 4 * J_i**2 + 4 * J_i + 2,
        'c14': 4 * J_i**2 + 4 * J_i - 3
    },
    'Q2Q2P': {
        'c0': 1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i - 1,
        'c13': -3 * J_i - 1,
        'c14': -3 * J_i - 6
    },
    'Q2Q(P)': {
        'c0': 1 / (60 * J_i * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i + 3,
        'c13': 3 - 3 * J_i,
        'c14': -3 * J_i - 2
    },
    'Q2Q(Q)': {
        'c0': 1 / (60 * J_i * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 4 * J_i**2 + 4 * J_i - 3,
        'c13': 4 * J_i**2 + 4 * J_i - 3,
        'c14': 4 * J_i**2 + 4 * J_i + 2
    },
    'Q2Q(R)': {
        'c0': 1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i - 1,
        'c13': -3 * J_i - 6,
        'c14': -3 * J_i - 1
    },
    'Q2R2R': {
        'c0': -1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 3 * J_i + 6,
        'c13': 3 * J_i + 1,
        'c14': 1 - 2 * J_i
    },
    'Q2R2Q': {
        'c0':
        -sqrt(J_i) * sqrt(J_i + 2) / (60 * (J_i + 1) *
                                      (2 * J_i + 1) * sqrt(2 * J_i + 3)),
        'c12':
        3,
        'c13':
        -2,
        'c14':
        3
    },
    'Q2R(Q)': {
        'c0': -1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 3 * J_i + 6,
        'c13': 1 - 2 * J_i,
        'c14': 3 * J_i + 1
    },
    'Q2R(R)': {
        'c0':
        -sqrt(J_i) * sqrt(J_i + 2) / (60 * (J_i + 1) *
                                      (2 * J_i + 1) * sqrt(2 * J_i + 3)),
        'c12':
        3,
        'c13':
        3,
        'c14':
        -2
    },
    'Q(P)R': {
        'c0':
        -sqrt(J_i - 1) * sqrt(J_i + 1) / (60 * J_i * sqrt(2 * J_i - 1) *
                                          (2 * J_i + 1)),
        'c12':
        3,
        'c13':
        3,
        'c14':
        -2
    },
    'Q(P)Q': {
        'c0': -1 / (60 * J_i * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 3 * J_i - 3,
        'c13': -2 * J_i - 3,
        'c14': 3 * J_i + 2
    },
    'Q(P)2P': {
        'c0':
        -sqrt(J_i - 1) * sqrt(J_i + 1) / (60 * J_i * sqrt(2 * J_i - 1) *
                                          (2 * J_i + 1)),
        'c12':
        3,
        'c13':
        3,
        'c14':
        -2
    },
    'Q(P)2Q': {
        'c0': -1 / (60 * J_i * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 3 * J_i - 3,
        'c13': -2 * J_i - 3,
        'c14': 3 * J_i + 2
    },
    'Q(P)(Q)': {
        'c0':
        -sqrt(J_i - 1) * sqrt(J_i + 1) / (60 * J_i * sqrt(2 * J_i - 1) *
                                          (2 * J_i + 1)),
        'c12':
        3,
        'c13':
        -2,
        'c14':
        3
    },
    'Q(P)(P)': {
        'c0': -1 / (60 * J_i * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 3 * J_i - 3,
        'c13': 3 * J_i + 2,
        'c14': -2 * J_i - 3
    },
    'Q(P)(2Q)': {
        'c0':
        -sqrt(J_i - 1) * sqrt(J_i + 1) / (60 * J_i * sqrt(2 * J_i - 1) *
                                          (2 * J_i + 1)),
        'c12':
        3,
        'c13':
        -2,
        'c14':
        3
    },
    'Q(P)(2R)': {
        'c0': -1 / (60 * J_i * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 3 * J_i - 3,
        'c13': 3 * J_i + 2,
        'c14': -2 * J_i - 3
    },
    'Q(Q)R': {
        'c0': 1 / (60 * J_i * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i + 3,
        'c13': 3 - 3 * J_i,
        'c14': -3 * J_i - 2
    },
    'Q(Q)Q': {
        'c0': 1 / (60 * J_i * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 4 * J_i**2 + 4 * J_i - 3,
        'c13': 4 * J_i**2 + 4 * J_i - 3,
        'c14': 4 * J_i**2 + 4 * J_i + 2
    },
    'Q(Q)P': {
        'c0': 1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i - 1,
        'c13': -3 * J_i - 6,
        'c14': -3 * J_i - 1
    },
    'Q(Q)2P': {
        'c0': 1 / (60 * J_i * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i + 3,
        'c13': 3 - 3 * J_i,
        'c14': -3 * J_i - 2
    },
    'Q(Q)2Q': {
        'c0': 1 / (60 * J_i * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 4 * J_i**2 + 4 * J_i - 3,
        'c13': 4 * J_i**2 + 4 * J_i - 3,
        'c14': 4 * J_i**2 + 4 * J_i + 2
    },
    'Q(Q)2R': {
        'c0': 1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i - 1,
        'c13': -3 * J_i - 6,
        'c14': -3 * J_i - 1
    },
    'Q(Q)(R)': {
        'c0': 1 / (60 * J_i * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i + 3,
        'c13': -3 * J_i - 2,
        'c14': 3 - 3 * J_i
    },
    'Q(Q)(Q)': {
        'c0': 1 / (60 * J_i * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 4 * J_i**2 + 4 * J_i - 3,
        'c13': 4 * J_i**2 + 4 * J_i + 2,
        'c14': 4 * J_i**2 + 4 * J_i - 3
    },
    'Q(Q)(P)': {
        'c0': 1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i - 1,
        'c13': -3 * J_i - 1,
        'c14': -3 * J_i - 6
    },
    'Q(Q)(2P)': {
        'c0': 1 / (60 * J_i * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i + 3,
        'c13': -3 * J_i - 2,
        'c14': 3 - 3 * J_i
    },
    'Q(Q)(2Q)': {
        'c0': 1 / (60 * J_i * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 4 * J_i**2 + 4 * J_i - 3,
        'c13': 4 * J_i**2 + 4 * J_i + 2,
        'c14': 4 * J_i**2 + 4 * J_i - 3
    },
    'Q(Q)(2R)': {
        'c0': 1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i - 1,
        'c13': -3 * J_i - 1,
        'c14': -3 * J_i - 6
    },
    'Q(R)Q': {
        'c0': -1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 3 * J_i + 6,
        'c13': 1 - 2 * J_i,
        'c14': 3 * J_i + 1
    },
    'Q(R)P': {
        'c0':
        -sqrt(J_i) * sqrt(J_i + 2) / (60 * (J_i + 1) *
                                      (2 * J_i + 1) * sqrt(2 * J_i + 3)),
        'c12':
        3,
        'c13':
        3,
        'c14':
        -2
    },
    'Q(R)2Q': {
        'c0': -1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 3 * J_i + 6,
        'c13': 1 - 2 * J_i,
        'c14': 3 * J_i + 1
    },
    'Q(R)2R': {
        'c0':
        -sqrt(J_i) * sqrt(J_i + 2) / (60 * (J_i + 1) *
                                      (2 * J_i + 1) * sqrt(2 * J_i + 3)),
        'c12':
        3,
        'c13':
        3,
        'c14':
        -2
    },
    'Q(R)(R)': {
        'c0': -1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 3 * J_i + 6,
        'c13': 3 * J_i + 1,
        'c14': 1 - 2 * J_i
    },
    'Q(R)(Q)': {
        'c0':
        -sqrt(J_i) * sqrt(J_i + 2) / (60 * (J_i + 1) *
                                      (2 * J_i + 1) * sqrt(2 * J_i + 3)),
        'c12':
        3,
        'c13':
        -2,
        'c14':
        3
    },
    'Q(R)(2P)': {
        'c0': -1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 3 * J_i + 6,
        'c13': 3 * J_i + 1,
        'c14': 1 - 2 * J_i
    },
    'Q(R)(2Q)': {
        'c0':
        -sqrt(J_i) * sqrt(J_i + 2) / (60 * (J_i + 1) *
                                      (2 * J_i + 1) * sqrt(2 * J_i + 3)),
        'c12':
        3,
        'c13':
        -2,
        'c14':
        3
    },
    'RRP': {
        'c0': 1 / (60 * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 1,
        'c13': 1,
        'c14': 6
    },
    'RRQ': {
        'c0': 1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i - 1,
        'c13': -3 * J_i - 1,
        'c14': -3 * J_i - 6
    },
    'RRR': {
        'c0':
        1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2) * (2 * J_i + 3)),
        'c12': 2 * J_i**2 - J_i,
        'c13': 12 * J_i**2 + 24 * J_i + 10,
        'c14': 2 * J_i**2 + 9 * J_i + 10
    },
    'RR(P)': {
        'c0': 1 / (60 * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 1,
        'c13': 6,
        'c14': 1
    },
    'RR(Q)': {
        'c0': 1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i - 1,
        'c13': -3 * J_i - 6,
        'c14': -3 * J_i - 1
    },
    'RR(R)': {
        'c0':
        1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2) * (2 * J_i + 3)),
        'c12': 2 * J_i**2 - J_i,
        'c13': 2 * J_i**2 + 9 * J_i + 10,
        'c14': 12 * J_i**2 + 24 * J_i + 10
    },
    'RQP': {
        'c0':
        -sqrt(J_i) * sqrt(J_i + 2) / (60 * (J_i + 1) *
                                      (2 * J_i + 1) * sqrt(2 * J_i + 3)),
        'c12':
        3,
        'c13':
        -2,
        'c14':
        3
    },
    'RQQ': {
        'c0': -1 / (60 * (J_i + 1) * sqrt(2 * J_i + 1) * (2 * J_i + 3)),
        'c12': 3 * J_i,
        'c13': 3 * J_i + 5,
        'c14': -2 * J_i - 5
    },
    'RQ(Q)': {
        'c0':
        -sqrt(J_i) * sqrt(J_i + 2) / (60 * (J_i + 1) *
                                      (2 * J_i + 1) * sqrt(2 * J_i + 3)),
        'c12':
        3,
        'c13':
        3,
        'c14':
        -2
    },
    'RQ(R)': {
        'c0': -1 / (60 * (J_i + 1) * sqrt(2 * J_i + 1) * (2 * J_i + 3)),
        'c12': 3 * J_i,
        'c13': -2 * J_i - 5,
        'c14': 3 * J_i + 5
    },
    'RPP': {
        'c0': 1 / (60 * sqrt(2 * J_i + 1) * (2 * J_i + 3)),
        'c12': 6,
        'c13': 1,
        'c14': 1
    },
    'RP(R)': {
        'c0': 1 / (60 * sqrt(2 * J_i + 1) * (2 * J_i + 3)),
        'c12': 6,
        'c13': 1,
        'c14': 1
    },
    'R2P2R': {
        'c0': 1 / (60 * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 1,
        'c13': 1,
        'c14': 6
    },
    'R2P2Q': {
        'c0': 1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i - 1,
        'c13': -3 * J_i - 1,
        'c14': -3 * J_i - 6
    },
    'R2P2P': {
        'c0':
        1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2) * (2 * J_i + 3)),
        'c12': 2 * J_i**2 - J_i,
        'c13': 12 * J_i**2 + 24 * J_i + 10,
        'c14': 2 * J_i**2 + 9 * J_i + 10
    },
    'R2P(P)': {
        'c0': 1 / (60 * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 1,
        'c13': 6,
        'c14': 1
    },
    'R2P(Q)': {
        'c0': 1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 2 * J_i - 1,
        'c13': -3 * J_i - 6,
        'c14': -3 * J_i - 1
    },
    'R2P(R)': {
        'c0':
        1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2) * (2 * J_i + 3)),
        'c12': 2 * J_i**2 - J_i,
        'c13': 2 * J_i**2 + 9 * J_i + 10,
        'c14': 12 * J_i**2 + 24 * J_i + 10
    },
    'R2Q2R': {
        'c0':
        -sqrt(J_i) * sqrt(J_i + 2) / (60 * (J_i + 1) *
                                      (2 * J_i + 1) * sqrt(2 * J_i + 3)),
        'c12':
        3,
        'c13':
        -2,
        'c14':
        3
    },
    'R2Q2Q': {
        'c0': -1 / (60 * (J_i + 1) * sqrt(2 * J_i + 1) * (2 * J_i + 3)),
        'c12': 3 * J_i,
        'c13': 3 * J_i + 5,
        'c14': -2 * J_i - 5
    },
    'R2Q(Q)': {
        'c0':
        -sqrt(J_i) * sqrt(J_i + 2) / (60 * (J_i + 1) *
                                      (2 * J_i + 1) * sqrt(2 * J_i + 3)),
        'c12':
        3,
        'c13':
        3,
        'c14':
        -2
    },
    'R2Q(R)': {
        'c0': -1 / (60 * (J_i + 1) * sqrt(2 * J_i + 1) * (2 * J_i + 3)),
        'c12': 3 * J_i,
        'c13': -2 * J_i - 5,
        'c14': 3 * J_i + 5
    },
    'R2R2R': {
        'c0': 1 / (60 * sqrt(2 * J_i + 1) * (2 * J_i + 3)),
        'c12': 6,
        'c13': 1,
        'c14': 1
    },
    'R2R(R)': {
        'c0': 1 / (60 * sqrt(2 * J_i + 1) * (2 * J_i + 3)),
        'c12': 6,
        'c13': 1,
        'c14': 1
    },
    'R(P)R': {
        'c0': 1 / (60 * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 6,
        'c13': 1,
        'c14': 1
    },
    'R(P)2P': {
        'c0': 1 / (60 * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 6,
        'c13': 1,
        'c14': 1
    },
    'R(P)(P)': {
        'c0': 1 / (60 * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 6,
        'c13': 1,
        'c14': 1
    },
    'R(P)(2R)': {
        'c0': 1 / (60 * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 6,
        'c13': 1,
        'c14': 1
    },
    'R(Q)R': {
        'c0': -1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 3 * J_i + 6,
        'c13': 1 - 2 * J_i,
        'c14': 3 * J_i + 1
    },
    'R(Q)Q': {
        'c0':
        -sqrt(J_i) * sqrt(J_i + 2) / (60 * (J_i + 1) *
                                      (2 * J_i + 1) * sqrt(2 * J_i + 3)),
        'c12':
        3,
        'c13':
        3,
        'c14':
        -2
    },
    'R(Q)2P': {
        'c0': -1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 3 * J_i + 6,
        'c13': 1 - 2 * J_i,
        'c14': 3 * J_i + 1
    },
    'R(Q)2Q': {
        'c0':
        -sqrt(J_i) * sqrt(J_i + 2) / (60 * (J_i + 1) *
                                      (2 * J_i + 1) * sqrt(2 * J_i + 3)),
        'c12':
        3,
        'c13':
        3,
        'c14':
        -2
    },
    'R(Q)(Q)': {
        'c0': -1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 3 * J_i + 6,
        'c13': 3 * J_i + 1,
        'c14': 1 - 2 * J_i
    },
    'R(Q)(P)': {
        'c0':
        -sqrt(J_i) * sqrt(J_i + 2) / (60 * (J_i + 1) *
                                      (2 * J_i + 1) * sqrt(2 * J_i + 3)),
        'c12':
        3,
        'c13':
        -2,
        'c14':
        3
    },
    'R(Q)(2Q)': {
        'c0': -1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2)),
        'c12': 3 * J_i + 6,
        'c13': 3 * J_i + 1,
        'c14': 1 - 2 * J_i
    },
    'R(Q)(2R)': {
        'c0':
        -sqrt(J_i) * sqrt(J_i + 2) / (60 * (J_i + 1) *
                                      (2 * J_i + 1) * sqrt(2 * J_i + 3)),
        'c12':
        3,
        'c13':
        -2,
        'c14':
        3
    },
    'R(R)R': {
        'c0':
        1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2) * (2 * J_i + 3)),
        'c12': 2 * J_i**2 + 9 * J_i + 10,
        'c13': 2 * J_i**2 - J_i,
        'c14': 12 * J_i**2 + 24 * J_i + 10
    },
    'R(R)Q': {
        'c0': 1 / (60 * (J_i + 1) * sqrt(2 * J_i + 1) * (2 * J_i + 3)),
        'c12': 2 * J_i + 5,
        'c13': -3 * J_i,
        'c14': -3 * J_i - 5
    },
    'R(R)P': {
        'c0': 1 / (60 * sqrt(2 * J_i + 1) * (2 * J_i + 3)),
        'c12': 1,
        'c13': 6,
        'c14': 1
    },
    'R(R)2P': {
        'c0':
        1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2) * (2 * J_i + 3)),
        'c12': 2 * J_i**2 + 9 * J_i + 10,
        'c13': 2 * J_i**2 - J_i,
        'c14': 12 * J_i**2 + 24 * J_i + 10
    },
    'R(R)2Q': {
        'c0': 1 / (60 * (J_i + 1) * sqrt(2 * J_i + 1) * (2 * J_i + 3)),
        'c12': 2 * J_i + 5,
        'c13': -3 * J_i,
        'c14': -3 * J_i - 5
    },
    'R(R)2R': {
        'c0': 1 / (60 * sqrt(2 * J_i + 1) * (2 * J_i + 3)),
        'c12': 1,
        'c13': 6,
        'c14': 1
    },
    'R(R)(R)': {
        'c0':
        1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2) * (2 * J_i + 3)),
        'c12': 2 * J_i**2 + 9 * J_i + 10,
        'c13': 12 * J_i**2 + 24 * J_i + 10,
        'c14': 2 * J_i**2 - J_i
    },
    'R(R)(Q)': {
        'c0': 1 / (60 * (J_i + 1) * sqrt(2 * J_i + 1) * (2 * J_i + 3)),
        'c12': 2 * J_i + 5,
        'c13': -3 * J_i - 5,
        'c14': -3 * J_i
    },
    'R(R)(P)': {
        'c0': 1 / (60 * sqrt(2 * J_i + 1) * (2 * J_i + 3)),
        'c12': 1,
        'c13': 1,
        'c14': 6
    },
    'R(R)(2P)': {
        'c0':
        1 / (60 * (J_i + 1) * (2 * J_i + 1)**Rational(3, 2) * (2 * J_i + 3)),
        'c12': 2 * J_i**2 + 9 * J_i + 10,
        'c13': 12 * J_i**2 + 24 * J_i + 10,
        'c14': 2 * J_i**2 - J_i
    },
    'R(R)(2Q)': {
        'c0': 1 / (60 * (J_i + 1) * sqrt(2 * J_i + 1) * (2 * J_i + 3)),
        'c12': 2 * J_i + 5,
        'c13': -3 * J_i - 5,
        'c14': -3 * J_i
    },
    'R(R)(2R)': {
        'c0': 1 / (60 * sqrt(2 * J_i + 1) * (2 * J_i + 3)),
        'c12': 1,
        'c13': 1,
        'c14': 6
    }
}

# ** High J limit
#: R-factors in high-J limit as dictionaries of coefficients
rfactors_highj_dict = {
    'PRR': {
        'c0': 1,
        'c12': Rational(1, 10),
        'c13': Rational(1, 60),
        'c14': Rational(1, 60)
    },
    'PR(P)': {
        'c0': 1,
        'c12': Rational(1, 10),
        'c13': Rational(1, 60),
        'c14': Rational(1, 60)
    },
    'PQQ': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'PQR': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'PQ(P)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'PQ(Q)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'PPP': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 10),
        'c14': Rational(1, 60)
    },
    'PPQ': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'PPR': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 60),
        'c14': Rational(1, 10)
    },
    'PP(P)': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 60),
        'c14': Rational(1, 10)
    },
    'PP(Q)': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'PP(R)': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 10),
        'c14': Rational(1, 60)
    },
    'P2P2P': {
        'c0': 1,
        'c12': Rational(1, 10),
        'c13': Rational(1, 60),
        'c14': Rational(1, 60)
    },
    'P2P(P)': {
        'c0': 1,
        'c12': Rational(1, 10),
        'c13': Rational(1, 60),
        'c14': Rational(1, 60)
    },
    'P2Q2Q': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'P2Q2P': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'P2Q(P)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'P2Q(Q)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'P2R2R': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 10),
        'c14': Rational(1, 60)
    },
    'P2R2Q': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'P2R2P': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 60),
        'c14': Rational(1, 10)
    },
    'P2R(P)': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 60),
        'c14': Rational(1, 10)
    },
    'P2R(Q)': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'P2R(R)': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 10),
        'c14': Rational(1, 60)
    },
    'P(P)R': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 10),
        'c14': Rational(1, 60)
    },
    'P(P)Q': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'P(P)P': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 60),
        'c14': Rational(1, 10)
    },
    'P(P)2P': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 10),
        'c14': Rational(1, 60)
    },
    'P(P)2Q': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'P(P)2R': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 60),
        'c14': Rational(1, 10)
    },
    'P(P)(R)': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 60),
        'c14': Rational(1, 10)
    },
    'P(P)(Q)': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'P(P)(P)': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 10),
        'c14': Rational(1, 60)
    },
    'P(P)(2P)': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 60),
        'c14': Rational(1, 10)
    },
    'P(P)(2Q)': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'P(P)(2R)': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 10),
        'c14': Rational(1, 60)
    },
    'P(Q)Q': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'P(Q)P': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'P(Q)2Q': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'P(Q)2R': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'P(Q)(R)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'P(Q)(Q)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'P(Q)(2P)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'P(Q)(2Q)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'P(R)P': {
        'c0': 1,
        'c12': Rational(1, 10),
        'c13': Rational(1, 60),
        'c14': Rational(1, 60)
    },
    'P(R)2R': {
        'c0': 1,
        'c12': Rational(1, 10),
        'c13': Rational(1, 60),
        'c14': Rational(1, 60)
    },
    'P(R)(R)': {
        'c0': 1,
        'c12': Rational(1, 10),
        'c13': Rational(1, 60),
        'c14': Rational(1, 60)
    },
    'P(R)(2P)': {
        'c0': 1,
        'c12': Rational(1, 10),
        'c13': Rational(1, 60),
        'c14': Rational(1, 60)
    },
    'QRQ': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'QRR': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'QR(P)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'QR(Q)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'QQP': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'QQQ': {
        'c0': 1,
        'c12': Rational(1, 15),
        'c13': Rational(1, 15),
        'c14': Rational(1, 15)
    },
    'QQR': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'QQ(P)': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'QQ(Q)': {
        'c0': 1,
        'c12': Rational(1, 15),
        'c13': Rational(1, 15),
        'c14': Rational(1, 15)
    },
    'QQ(R)': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'QPP': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'QPQ': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'QP(Q)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'QP(R)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'Q2P2Q': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'Q2P2P': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'Q2P(P)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'Q2P(Q)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'Q2Q2R': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'Q2Q2Q': {
        'c0': 1,
        'c12': Rational(1, 15),
        'c13': Rational(1, 15),
        'c14': Rational(1, 15)
    },
    'Q2Q2P': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'Q2Q(P)': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'Q2Q(Q)': {
        'c0': 1,
        'c12': Rational(1, 15),
        'c13': Rational(1, 15),
        'c14': Rational(1, 15)
    },
    'Q2Q(R)': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'Q2R2R': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'Q2R2Q': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'Q2R(Q)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'Q2R(R)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'Q(P)R': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'Q(P)Q': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'Q(P)2P': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'Q(P)2Q': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'Q(P)(Q)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'Q(P)(P)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'Q(P)(2Q)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'Q(P)(2R)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'Q(Q)R': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'Q(Q)Q': {
        'c0': 1,
        'c12': Rational(1, 15),
        'c13': Rational(1, 15),
        'c14': Rational(1, 15)
    },
    'Q(Q)P': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'Q(Q)2P': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'Q(Q)2Q': {
        'c0': 1,
        'c12': Rational(1, 15),
        'c13': Rational(1, 15),
        'c14': Rational(1, 15)
    },
    'Q(Q)2R': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'Q(Q)(R)': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'Q(Q)(Q)': {
        'c0': 1,
        'c12': Rational(1, 15),
        'c13': Rational(1, 15),
        'c14': Rational(1, 15)
    },
    'Q(Q)(P)': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'Q(Q)(2P)': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'Q(Q)(2Q)': {
        'c0': 1,
        'c12': Rational(1, 15),
        'c13': Rational(1, 15),
        'c14': Rational(1, 15)
    },
    'Q(Q)(2R)': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'Q(R)Q': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'Q(R)P': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'Q(R)2Q': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'Q(R)2R': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'Q(R)(R)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'Q(R)(Q)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'Q(R)(2P)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'Q(R)(2Q)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'RRP': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 60),
        'c14': Rational(1, 10)
    },
    'RRQ': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'RRR': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 10),
        'c14': Rational(1, 60)
    },
    'RR(P)': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 10),
        'c14': Rational(1, 60)
    },
    'RR(Q)': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'RR(R)': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 60),
        'c14': Rational(1, 10)
    },
    'RQP': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'RQQ': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'RQ(Q)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'RQ(R)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'RPP': {
        'c0': 1,
        'c12': Rational(1, 10),
        'c13': Rational(1, 60),
        'c14': Rational(1, 60)
    },
    'RP(R)': {
        'c0': 1,
        'c12': Rational(1, 10),
        'c13': Rational(1, 60),
        'c14': Rational(1, 60)
    },
    'R2P2R': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 60),
        'c14': Rational(1, 10)
    },
    'R2P2Q': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'R2P2P': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 10),
        'c14': Rational(1, 60)
    },
    'R2P(P)': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 10),
        'c14': Rational(1, 60)
    },
    'R2P(Q)': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'R2P(R)': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 60),
        'c14': Rational(1, 10)
    },
    'R2Q2R': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'R2Q2Q': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'R2Q(Q)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'R2Q(R)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'R2R2R': {
        'c0': 1,
        'c12': Rational(1, 10),
        'c13': Rational(1, 60),
        'c14': Rational(1, 60)
    },
    'R2R(R)': {
        'c0': 1,
        'c12': Rational(1, 10),
        'c13': Rational(1, 60),
        'c14': Rational(1, 60)
    },
    'R(P)R': {
        'c0': 1,
        'c12': Rational(1, 10),
        'c13': Rational(1, 60),
        'c14': Rational(1, 60)
    },
    'R(P)2P': {
        'c0': 1,
        'c12': Rational(1, 10),
        'c13': Rational(1, 60),
        'c14': Rational(1, 60)
    },
    'R(P)(P)': {
        'c0': 1,
        'c12': Rational(1, 10),
        'c13': Rational(1, 60),
        'c14': Rational(1, 60)
    },
    'R(P)(2R)': {
        'c0': 1,
        'c12': Rational(1, 10),
        'c13': Rational(1, 60),
        'c14': Rational(1, 60)
    },
    'R(Q)R': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'R(Q)Q': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'R(Q)2P': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'R(Q)2Q': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'R(Q)(Q)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'R(Q)(P)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'R(Q)(2Q)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(-1, 20),
        'c14': Rational(1, 30)
    },
    'R(Q)(2R)': {
        'c0': 1,
        'c12': Rational(-1, 20),
        'c13': Rational(1, 30),
        'c14': Rational(-1, 20)
    },
    'R(R)R': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 60),
        'c14': Rational(1, 10)
    },
    'R(R)Q': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'R(R)P': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 10),
        'c14': Rational(1, 60)
    },
    'R(R)2P': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 60),
        'c14': Rational(1, 10)
    },
    'R(R)2Q': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'R(R)2R': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 10),
        'c14': Rational(1, 60)
    },
    'R(R)(R)': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 10),
        'c14': Rational(1, 60)
    },
    'R(R)(Q)': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'R(R)(P)': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 60),
        'c14': Rational(1, 10)
    },
    'R(R)(2P)': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 10),
        'c14': Rational(1, 60)
    },
    'R(R)(2Q)': {
        'c0': 1,
        'c12': Rational(1, 30),
        'c13': Rational(-1, 20),
        'c14': Rational(-1, 20)
    },
    'R(R)(2R)': {
        'c0': 1,
        'c12': Rational(1, 60),
        'c13': Rational(1, 60),
        'c14': Rational(1, 10)
    }
}

# * R-factor classes
#: Map between high-J R-factor coefficient tuples and :math:`\Theta_\alpha` classes
theta_labels = {
    (Rational(1, 10), Rational(1, 60), Rational(1, 60)): '1',
    (Rational(1, 20), Rational(1, 20), Rational(1, -30)): '2',
    (Rational(1, 20), Rational(-1, 30), Rational(1, 20)): '3',
    (Rational(1, 15), Rational(1, 15), Rational(1, 15)): '4',
    (Rational(1, 60), Rational(1, 10), Rational(1, 60)): '5',
    (Rational(1, 60), Rational(1, 60), Rational(1, 10)): '6',
    (Rational(1, 30), Rational(-1, 20), Rational(-1, 20)): '7'
}

# * Polarization conditions
#: Dictionary of polarization conditions and angle substitutions
magic_angle = atan(sqrt(2))
pac_angle = asin(2*sqrt(7)/7)
angle_subs = {
    'MA': dict(zip(thetas, (0, 0, magic_angle, magic_angle))),
    'AMA': dict(zip(thetas, (0, magic_angle, 0, magic_angle))),
    'Vac2': dict(zip(thetas, (pi/2, pi/4, pi/2, atan(Rational(1,2))))),
    'MMA': dict(zip(thetas, (0, magic_angle, magic_angle, 0))),
    'PAC': dict(zip(thetas, (0, 0, pac_angle, -pac_angle))),
    'MPAC': dict(zip(thetas, (0, pac_angle, -pac_angle, 0))),
    'APAC': dict(zip(thetas, (0, pac_angle, 0, -pac_angle))),
    'Vac1': dict(zip(thetas, (pi/2, pi/4, pi/2, -atan(Rational(3,4))))),
}
