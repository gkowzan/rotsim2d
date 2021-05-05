"""Derived symbolic expressions for angular momentum dependence and polarization
dependence of four-fold dipole interaction operator.

All the contents can obtained by using :mod:`rotsim2d.angular.symbolic` functions.  See
:download:`examples/generate_angular.symbolic_results.py <../../examples/generate_angular.symbolic_results.py>`.
"""
# * Imports
from rotsim2d.symbolic.common import *
from molspecutils.molecule import DiatomState

# * G-factors
# ** General expressions
#: :meta hide-value:
#: Isotropic version of angular momentum-dependent combinations of three
#: Wigner-6j factors involved in evaluating four-fold dipole interaction
#: operator.
gfactors = {'PPR': (0, 0, sqrt(5)/(5*sqrt(2*J_i + 1)*Abs(2*J_i - 1))),
 'RRP': (0, 0, sqrt(5)/(5*sqrt(2*J_i + 1)*(2*J_i + 3))),
 'PRP': (1/(3*(2*J_i + 1)**(Rational(3, 2))),
  -sqrt(3)*(J_i + 1)/(6*J_i*(2*J_i + 1)**(Rational(3, 2))),
  sqrt(5)*(J_i + 1)*(2*J_i + 3)/(30*J_i*(2*J_i - 1)*(2*J_i + 1)**(Rational(3, 2)))),
 'RPR': (1/(3*(2*J_i + 1)**(Rational(3, 2))),
  -sqrt(3)*J_i/(6*(J_i + 1)*(2*J_i + 1)**(Rational(3, 2))),
  sqrt(5)*J_i*(2*J_i - 1)/(30*(J_i + 1)*(2*J_i + 1)**(Rational(3, 2))*(2*J_i + 3))),
 'RPP': (1/(3*(2*J_i + 1)**(Rational(3, 2))),
  sqrt(3)/(6*(2*J_i + 1)**(Rational(3, 2))),
  sqrt(5)/(30*(2*J_i + 1)**(Rational(3, 2)))),
 'PRR': (1/(3*(2*J_i + 1)**(Rational(3, 2))),
  sqrt(3)/(6*(2*J_i + 1)**(Rational(3, 2))),
  sqrt(5)/(30*(2*J_i + 1)**(Rational(3, 2)))),
 'PQQ': (0,
  -sqrt(3)*(J_i - 1)/(6*J_i*(2*J_i - 1)*sqrt(2*J_i + 1)),
  sqrt(5)*(J_i + 1)/(10*J_i*sqrt(2*J_i + 1)*Abs(2*J_i - 1))),
 'QPQ': (0,
  -sqrt(3)*sqrt(J_i - 1)*sqrt(J_i + 1)/(6*J_i*sqrt(2*J_i - 1)*(2*J_i + 1)),
  -sqrt(5)*sqrt(J_i - 1)*sqrt(J_i + 1)/(10*J_i*sqrt(2*J_i - 1)*(2*J_i + 1))),
 'QPR': (0,
  -sqrt(3)*(J_i + 1)/(6*J_i*(2*J_i + 1)**(Rational(3, 2))),
  sqrt(5)*Abs(J_i - 1)/(10*J_i*(2*J_i + 1)**(Rational(3, 2)))),
 'QQP': (1/(3*(2*J_i + 1)**(Rational(3, 2))),
  -sqrt(3)/(6*J_i*(2*J_i + 1)**(Rational(3, 2))),
  -sqrt(5)*(2*J_i + 3)/(30*J_i*(2*J_i + 1)**(Rational(3, 2)))),
 'QQQ': (1/(3*(2*J_i + 1)**(Rational(3, 2))),
  -sqrt(3)/(6*J_i*(J_i + 1)*(2*J_i + 1)**(Rational(3, 2))),
  sqrt(5)*(2*J_i - 1)*(2*J_i + 3)/(30*J_i*(J_i + 1)*(2*J_i + 1)**(Rational(3, 2)))),
 'QRP': (0,
  -sqrt(3)*J_i/(6*(J_i + 1)*(2*J_i + 1)**(Rational(3, 2))),
  sqrt(5)*(J_i + 2)/(10*(J_i + 1)*(2*J_i + 1)**(Rational(3, 2)))),
 'RPQ': (1/(3*(2*J_i + 1)**(Rational(3, 2))),
  sqrt(3)/(6*(J_i + 1)*(2*J_i + 1)**(Rational(3, 2))),
  -sqrt(5)*Abs(2*J_i - 1)/(30*(J_i + 1)*(2*J_i + 1)**(Rational(3, 2)))),
 'RQP': (0,
  -sqrt(3)*sqrt(J_i)*sqrt(J_i + 2)/(6*(J_i + 1)*(2*J_i + 1)*sqrt(2*J_i + 3)),
  -sqrt(5)*sqrt(J_i)*sqrt(J_i + 2)/(10*(J_i + 1)*(2*J_i + 1)*sqrt(2*J_i + 3))),
 'RQQ': (0,
  -sqrt(3)*(J_i + 2)/(6*(J_i + 1)*sqrt(2*J_i + 1)*(2*J_i + 3)),
  sqrt(5)*J_i/(10*(J_i + 1)*sqrt(2*J_i + 1)*(2*J_i + 3))),
 'PQR': (0,
  -sqrt(3)*sqrt((J_i - 1)*(2*J_i - 1))*sqrt(J_i + 1)/(6*J_i*(2*J_i - 1)*(2*J_i + 1)),
  -sqrt(5)*sqrt(J_i - 1)*sqrt(J_i + 1)/(10*J_i*sqrt(2*J_i - 1)*(2*J_i + 1))),
 'PRQ': (1/(3*(2*J_i + 1)**(Rational(3, 2))),
  -sqrt(3)/(6*J_i*(2*J_i + 1)**(Rational(3, 2))),
  -sqrt(5)*(2*J_i + 3)/(30*J_i*(2*J_i + 1)**(Rational(3, 2)))),
 'QQR': (1/(3*(2*J_i + 1)**(Rational(3, 2))),
  sqrt(3)/(6*(J_i + 1)*(2*J_i + 1)**(Rational(3, 2))),
  -sqrt(5)*Abs(2*J_i - 1)/(30*(J_i + 1)*(2*J_i + 1)**(Rational(3, 2)))),
 'QRQ': (0,
  -sqrt(3)*sqrt(J_i)*sqrt(J_i + 2)/(6*(J_i + 1)*(2*J_i + 1)*sqrt(2*J_i + 3)),
  -sqrt(5)*sqrt(J_i)*sqrt(J_i + 2)/(10*(J_i + 1)*(2*J_i + 1)*sqrt(2*J_i + 3)))}

# ** High-J limit
# gfactors_highj = {k: tuple(nsimplify(limit(x*(2*J_i+1)**(3/2), J_i, oo)) for x in v)
#                   for k, v in gfactors.items()}
#: :meta hide-value:
#: G-factors multiplied by `(2*J_i+1)**(3/2)` and with :math:`J_i\to\infty`.
gfactors_highj = {'PPR': (0, 0, sqrt(5)/5),
 'RRP': (0, 0, sqrt(5)/5),
 'PRP': (Rational(1, 3), -sqrt(3)/6, sqrt(5)/30),
 'RPR': (Rational(1, 3), -sqrt(3)/6, sqrt(5)/30),
 'RPP': (Rational(1, 3), sqrt(3)/6, sqrt(5)/30),
 'PRR': (Rational(1, 3), sqrt(3)/6, sqrt(5)/30),
 'PQQ': (0, -sqrt(3)/6, sqrt(5)/10),
 'QPQ': (0, -sqrt(3)/6, -sqrt(5)/10),
 'QPR': (0, -sqrt(3)/6, sqrt(5)/10),
 'QQP': (Rational(1, 3), 0, -sqrt(5)/15),
 'QQQ': (Rational(1, 3), 0, 2*sqrt(5)/15),
 'QRP': (0, -sqrt(3)/6, sqrt(5)/10),
 'RPQ': (Rational(1, 3), 0, -sqrt(5)/15),
 'RQP': (0, -sqrt(3)/6, -sqrt(5)/10),
 'RQQ': (0, -sqrt(3)/6, sqrt(5)/10),
 'PQR': (0, -sqrt(3)/6, -sqrt(5)/10),
 'PRQ': (Rational(1, 3), 0, -sqrt(5)/15),
 'QQR': (Rational(1, 3), 0, -sqrt(5)/15),
 'QRQ': (0, -sqrt(3)/6, -sqrt(5)/10)}

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
# rfactors = {k: rfactorize(gfactors[k], T00_exprs) for k in gfactors.keys()}
#: :meta hide-value:
#: R-factors, 14 expressions are unique, 12 if one ignores the common factor
rfactors = {'PPR': (cos(phi - phj - phk + phl) + cos(phi - phj + phk - phl) + 6*cos(phi + phj - phk - phl))/(60*(2*J_i - 1)*sqrt(2*J_i + 1)),
 'RRP': (cos(phi - phj - phk + phl) + cos(phi - phj + phk - phl) + 6*cos(phi + phj - phk - phl))/(60*sqrt(2*J_i + 1)*(2*J_i + 3)),
 'PRP': ((12*J_i**2 - 2)*cos(phi - phj + phk - phl) + (2*J_i**2 - 5*J_i + 3)*cos(phi - phj - phk + phl) + (2*J_i**2 + 5*J_i + 3)*cos(phi + phj - phk - phl))/(60*J_i*(2*J_i - 1)*(2*J_i + 1)**(Rational(3, 2))),
 'RPR': ((2*J_i**2 - J_i)*cos(phi + phj - phk - phl) + (2*J_i**2 + 9*J_i + 10)*cos(phi - phj - phk + phl) + (12*J_i**2 + 24*J_i + 10)*cos(phi - phj + phk - phl))/(60*(J_i + 1)*(2*J_i + 1)**(Rational(3, 2))*(2*J_i + 3)),
 'RPP': (6*cos(phi - phj - phk + phl) + cos(phi - phj + phk - phl) + cos(phi + phj - phk - phl))/(60*(2*J_i + 1)**(Rational(3, 2))),
 'PRR': (6*cos(phi - phj - phk + phl) + cos(phi - phj + phk - phl) + cos(phi + phj - phk - phl))/(60*(2*J_i + 1)**(Rational(3, 2))),
 'PQQ': ((3 - 2*J_i)*cos(phi - phj - phk + phl) + (3*J_i - 2)*cos(phi - phj + phk - phl) + (3*J_i + 3)*cos(phi + phj - phk - phl))/(60*J_i*(2*J_i - 1)*sqrt(2*J_i + 1)),
 'QPQ': -sqrt(J_i - 1)*sqrt(J_i + 1)*(3*cos(phi - phj - phk + phl) - 2*cos(phi - phj + phk - phl) + 3*cos(phi + phj - phk - phl))/(60*J_i*sqrt(2*J_i - 1)*(2*J_i + 1)),
 'QPR': ((-2*J_i - 3)*cos(phi - phj - phk + phl) + (3*J_i - 3)*cos(phi + phj - phk - phl) + (3*J_i + 2)*cos(phi - phj + phk - phl))/(60*J_i*(2*J_i + 1)**(Rational(3, 2))),
 'QQP': -((3 - 3*J_i)*cos(phi - phj - phk + phl) + (-3*J_i - 2)*cos(phi - phj + phk - phl) + (2*J_i + 3)*cos(phi + phj - phk - phl))/(60*J_i*(2*J_i + 1)**(Rational(3, 2))),
 'QQQ': ((4*J_i**2 + 4*J_i - 3)*cos(phi - phj - phk + phl) + (4*J_i**2 + 4*J_i - 3)*cos(phi + phj - phk - phl) + (4*J_i**2 + 4*J_i + 2)*cos(phi - phj + phk - phl))/(60*J_i*(J_i + 1)*(2*J_i + 1)**(Rational(3, 2))),
 'QRP': ((1 - 2*J_i)*cos(phi - phj - phk + phl) + (3*J_i + 1)*cos(phi - phj + phk - phl) + (3*J_i + 6)*cos(phi + phj - phk - phl))/(60*(J_i + 1)*(2*J_i + 1)**(Rational(3, 2))),
 'RPQ': -((-3*J_i - 6)*cos(phi - phj - phk + phl) + (-3*J_i - 1)*cos(phi - phj + phk - phl) + (2*J_i - 1)*cos(phi + phj - phk - phl))/(60*(J_i + 1)*(2*J_i + 1)**(Rational(3, 2))),
 'RQP': -sqrt(J_i)*sqrt(J_i + 2)*(3*cos(phi - phj - phk + phl) - 2*cos(phi - phj + phk - phl) + 3*cos(phi + phj - phk - phl))/(60*(J_i + 1)*(2*J_i + 1)*sqrt(2*J_i + 3)),
 'RQQ': (3*J_i*cos(phi + phj - phk - phl) + (-2*J_i - 5)*cos(phi - phj - phk + phl) + (3*J_i + 5)*cos(phi - phj + phk - phl))/(60*(J_i + 1)*sqrt(2*J_i + 1)*(2*J_i + 3)),
 'PQR': -(3*sqrt(J_i - 1)*sqrt(2*J_i**2 + 3*J_i + 1)*cos(phi - phj - phk + phl) - 2*sqrt(J_i - 1)*sqrt(2*J_i**2 + 3*J_i + 1)*cos(phi - phj + phk - phl) + 3*sqrt(2*J_i**3 + J_i**2 - 2*J_i - 1)*cos(phi + phj - phk - phl))/(60*J_i*sqrt(2*J_i - 1)*(2*J_i + 1)**(Rational(3, 2))),
 'PRQ': -((3 - 3*J_i)*cos(phi - phj - phk + phl) + (-3*J_i - 2)*cos(phi - phj + phk - phl) + (2*J_i + 3)*cos(phi + phj - phk - phl))/(60*J_i*(2*J_i + 1)**(Rational(3, 2))),
 'QQR': -((-3*J_i - 6)*cos(phi - phj - phk + phl) + (-3*J_i - 1)*cos(phi - phj + phk - phl) + (2*J_i - 1)*cos(phi + phj - phk - phl))/(60*(J_i + 1)*(2*J_i + 1)**(Rational(3, 2))),
 'QRQ': -sqrt(J_i)*sqrt(J_i + 2)*(3*cos(phi - phj - phk + phl) - 2*cos(phi - phj + phk - phl) + 3*cos(phi + phj - phk - phl))/(60*(J_i + 1)*(2*J_i + 1)*sqrt(2*J_i + 3))}

# rfactors_dict = {k: rfactorize(gfactors[k], T00_exprs, coeffs=True) for k in gfactors.keys()}
rfactors_dict = {'PPR': {'c0': 1/(60*(2*J_i - 1)*sqrt(2*J_i + 1)),
  'c12': 6,
  'c13': 1,
  'c14': 1},
 'RRP': {'c0': 1/(60*sqrt(2*J_i + 1)*(2*J_i + 3)),
  'c12': 6,
  'c13': 1,
  'c14': 1},
 'PRP': {'c0': 1/(60*J_i*(2*J_i - 1)*(2*J_i + 1)**(Rational(3, 2))),
  'c12': 2*J_i**2 + 5*J_i + 3,
  'c13': 12*J_i**2 - 2,
  'c14': 2*J_i**2 - 5*J_i + 3},
 'RPR': {'c0': 1/(60*(J_i + 1)*(2*J_i + 1)**(Rational(3, 2))*(2*J_i + 3)),
  'c12': 2*J_i**2 - J_i,
  'c13': 12*J_i**2 + 24*J_i + 10,
  'c14': 2*J_i**2 + 9*J_i + 10},
 'RPP': {'c0': 1/(60*(2*J_i + 1)**(Rational(3, 2))), 'c12': 1, 'c13': 1, 'c14': 6},
 'PRR': {'c0': 1/(60*(2*J_i + 1)**(Rational(3, 2))), 'c12': 1, 'c13': 1, 'c14': 6},
 'PQQ': {'c0': 1/(60*J_i*(2*J_i - 1)*sqrt(2*J_i + 1)),
  'c12': 3*J_i + 3,
  'c13': 3*J_i - 2,
  'c14': 3 - 2*J_i},
 'QPQ': {'c0': -sqrt(J_i - 1)*sqrt(J_i + 1)/(60*J_i*sqrt(2*J_i - 1)*(2*J_i + 1)),
  'c12': 3,
  'c13': -2,
  'c14': 3},
 'QPR': {'c0': 1/(60*J_i*(2*J_i + 1)**(Rational(3, 2))),
  'c12': 3*J_i - 3,
  'c13': 3*J_i + 2,
  'c14': -2*J_i - 3},
 'QQP': {'c0': -1/(60*J_i*(2*J_i + 1)**(Rational(3, 2))),
  'c12': 2*J_i + 3,
  'c13': -3*J_i - 2,
  'c14': 3 - 3*J_i},
 'QQQ': {'c0': 1/(60*J_i*(J_i + 1)*(2*J_i + 1)**(Rational(3, 2))),
  'c12': 4*J_i**2 + 4*J_i - 3,
  'c13': 4*J_i**2 + 4*J_i + 2,
  'c14': 4*J_i**2 + 4*J_i - 3},
 'QRP': {'c0': 1/(60*(J_i + 1)*(2*J_i + 1)**(Rational(3, 2))),
  'c12': 3*J_i + 6,
  'c13': 3*J_i + 1,
  'c14': 1 - 2*J_i},
 'RPQ': {'c0': -1/(60*(J_i + 1)*(2*J_i + 1)**(Rational(3, 2))),
  'c12': 2*J_i - 1,
  'c13': -3*J_i - 1,
  'c14': -3*J_i - 6},
 'RQP': {'c0': -sqrt(J_i)*sqrt(J_i + 2)/(60*(J_i + 1)*(2*J_i + 1)*sqrt(2*J_i + 3)),
  'c12': 3,
  'c13': -2,
  'c14': 3},
 'RQQ': {'c0': 1/(60*(J_i + 1)*sqrt(2*J_i + 1)*(2*J_i + 3)),
  'c12': 3*J_i,
  'c13': 3*J_i + 5,
  'c14': -2*J_i - 5},
 'PQR': {'c0': -1/(60*J_i*sqrt(2*J_i - 1)*(2*J_i + 1)),
  'c12': 3*sqrt(J_i**2 - 1),
  'c13': -2*sqrt(J_i - 1)*sqrt(J_i + 1),
  'c14': 3*sqrt(J_i - 1)*sqrt(J_i + 1)},
 'PRQ': {'c0': -1/(60*J_i*(2*J_i + 1)**(Rational(3, 2))),
  'c12': 2*J_i + 3,
  'c13': -3*J_i - 2,
  'c14': 3 - 3*J_i},
 'QQR': {'c0': -1/(60*(J_i + 1)*(2*J_i + 1)**(Rational(3, 2))),
  'c12': 2*J_i - 1,
  'c13': -3*J_i - 1,
  'c14': -3*J_i - 6},
 'QRQ': {'c0': -sqrt(J_i)*sqrt(J_i + 2)/(60*(J_i + 1)*(2*J_i + 1)*sqrt(2*J_i + 3)),
  'c12': 3,
  'c13': -2,
  'c14': 3}}

# ** XXXX-polarization
# rfactors_xxxx = {k: factor(powdenest(v.subs({phi: S(0), phj: S(0), phk: S(0), phl: S(0)}), force=True), deep=True)
#                  for k, v in rfactors.items()}
rfactors_xxxx = {'PPR': 2/(15*(2*J_i - 1)*sqrt(2*J_i + 1)),
 'RRP': 2/(15*sqrt(2*J_i + 1)*(2*J_i + 3)),
 'PRP': (4*J_i**2 + 1)/(15*J_i*(2*J_i - 1)*(2*J_i + 1)**(Rational(3, 2))),
 'RPR': (4*J_i**2 + 8*J_i + 5)/(15*(J_i + 1)*(2*J_i + 1)**(Rational(3, 2))*(2*J_i + 3)),
 'RPP': 2/(15*(2*J_i + 1)**(Rational(3, 2))),
 'PRR': 2/(15*(2*J_i + 1)**(Rational(3, 2))),
 'PQQ': (J_i + 1)/(15*J_i*(2*J_i - 1)*sqrt(2*J_i + 1)),
 'QPQ': -sqrt(J_i - 1)*sqrt(J_i + 1)/(15*J_i*sqrt(2*J_i - 1)*(2*J_i + 1)),
 'QPR': (J_i - 1)/(15*J_i*(2*J_i + 1)**(Rational(3, 2))),
 'QQP': (J_i - 1)/(15*J_i*(2*J_i + 1)**(Rational(3, 2))),
 'QQQ': (3*J_i**2 + 3*J_i - 1)/(15*J_i*(J_i + 1)*(2*J_i + 1)**(Rational(3, 2))),
 'QRP': (J_i + 2)/(15*(J_i + 1)*(2*J_i + 1)**(Rational(3, 2))),
 'RPQ': (J_i + 2)/(15*(J_i + 1)*(2*J_i + 1)**(Rational(3, 2))),
 'RQP': -sqrt(J_i)*sqrt(J_i + 2)/(15*(J_i + 1)*(2*J_i + 1)*sqrt(2*J_i + 3)),
 'RQQ': J_i/(15*(J_i + 1)*sqrt(2*J_i + 1)*(2*J_i + 3)),
 'PQR': -sqrt(J_i - 1)*sqrt(J_i + 1)/(15*J_i*sqrt(2*J_i - 1)*(2*J_i + 1)),
 'PRQ': (J_i - 1)/(15*J_i*(2*J_i + 1)**(Rational(3, 2))),
 'QQR': (J_i + 2)/(15*(J_i + 1)*(2*J_i + 1)**(Rational(3, 2))),
 'QRQ': -sqrt(J_i)*sqrt(J_i + 2)/(15*(J_i + 1)*(2*J_i + 1)*sqrt(2*J_i + 3))}

# ** R-factors relative to XXXX
# rfactors_relative = {k: rfactorize(gfactors[k], T00_exprs, relative=True) for k in gfactors.keys()}
rfactors_relative = {'PPR': cos(phi - phj - phk + phl)/8 + cos(phi - phj + phk - phl)/8 + 3*cos(phi + phj - phk - phl)/4,
 'RRP': cos(phi - phj - phk + phl)/8 + cos(phi - phj + phk - phl)/8 + 3*cos(phi + phj - phk - phl)/4,
 'PRP': ((12*J_i**2 - 2)*cos(phi - phj + phk - phl) + (2*J_i**2 - 5*J_i + 3)*cos(phi - phj - phk + phl) + (2*J_i**2 + 5*J_i + 3)*cos(phi + phj - phk - phl))/(4*(4*J_i**2 + 1)),
 'RPR': ((2*J_i**2 - J_i)*cos(phi + phj - phk - phl) + (2*J_i**2 + 9*J_i + 10)*cos(phi - phj - phk + phl) + (12*J_i**2 + 24*J_i + 10)*cos(phi - phj + phk - phl))/(4*(4*J_i**2 + 8*J_i + 5)),
 'RPP': 3*cos(phi - phj - phk + phl)/4 + cos(phi - phj + phk - phl)/8 + cos(phi + phj - phk - phl)/8,
 'PRR': 3*cos(phi - phj - phk + phl)/4 + cos(phi - phj + phk - phl)/8 + cos(phi + phj - phk - phl)/8,
 'PQQ': ((3 - 2*J_i)*cos(phi - phj - phk + phl) + (3*J_i - 2)*cos(phi - phj + phk - phl) + (3*J_i + 3)*cos(phi + phj - phk - phl))/(4*(J_i + 1)),
 'QPQ': 3*cos(phi - phj - phk + phl)/4 - cos(phi - phj + phk - phl)/2 + 3*cos(phi + phj - phk - phl)/4,
 'QPR': ((-2*J_i - 3)*cos(phi - phj - phk + phl) + (3*J_i - 3)*cos(phi + phj - phk - phl) + (3*J_i + 2)*cos(phi - phj + phk - phl))/(4*(J_i - 1)),
 'QQP': -((3 - 3*J_i)*cos(phi - phj - phk + phl) + (-3*J_i - 2)*cos(phi - phj + phk - phl) + (2*J_i + 3)*cos(phi + phj - phk - phl))/(4*(J_i - 1)),
 'QQQ': ((4*J_i**2 + 4*J_i - 3)*cos(phi - phj - phk + phl) + (4*J_i**2 + 4*J_i - 3)*cos(phi + phj - phk - phl) + (4*J_i**2 + 4*J_i + 2)*cos(phi - phj + phk - phl))/(4*(3*J_i**2 + 3*J_i - 1)),
 'QRP': ((1 - 2*J_i)*cos(phi - phj - phk + phl) + (3*J_i + 1)*cos(phi - phj + phk - phl) + (3*J_i + 6)*cos(phi + phj - phk - phl))/(4*(J_i + 2)),
 'RPQ': -((-3*J_i - 6)*cos(phi - phj - phk + phl) + (-3*J_i - 1)*cos(phi - phj + phk - phl) + (2*J_i - 1)*cos(phi + phj - phk - phl))/(4*(J_i + 2)),
 'RQP': 3*cos(phi - phj - phk + phl)/4 - cos(phi - phj + phk - phl)/2 + 3*cos(phi + phj - phk - phl)/4,
 'RQQ': (3*J_i*cos(phi + phj - phk - phl) + (-2*J_i - 5)*cos(phi - phj - phk + phl) + (3*J_i + 5)*cos(phi - phj + phk - phl))/(4*J_i),
 'PQR': (6*sqrt(J_i - 1)*sqrt(2*J_i - 1)*cos(phi + phj - phk - phl) + (-5*sqrt(J_i - 1)*sqrt(2*J_i - 1) + sqrt(2*J_i**2 - 3*J_i + 1))*cos(phi - phj + phk - phl) + (5*sqrt(J_i - 1)*sqrt(2*J_i - 1) + sqrt(2*J_i**2 - 3*J_i + 1))*cos(phi - phj - phk + phl))/(8*sqrt(J_i - 1)*sqrt(2*J_i - 1)),
 'PRQ': -((3 - 3*J_i)*cos(phi - phj - phk + phl) + (-3*J_i - 2)*cos(phi - phj + phk - phl) + (2*J_i + 3)*cos(phi + phj - phk - phl))/(4*(J_i - 1)),
 'QQR': -((-3*J_i - 6)*cos(phi - phj - phk + phl) + (-3*J_i - 1)*cos(phi - phj + phk - phl) + (2*J_i - 1)*cos(phi + phj - phk - phl))/(4*(J_i + 2)),
 'QRQ': 3*cos(phi - phj - phk + phl)/4 - cos(phi - phj + phk - phl)/2 + 3*cos(phi + phj - phk - phl)/4}

# ** High-J R-factors
# These are already multiplied by common factor of (2*J_i+1)^(1.5) and there is
# no common factor between all cosine terms (c0=1)
# rfactors_highj = {k: rfactorize(gfactors_highj[k], T00_exprs) for k in gfactors_highj.keys()}
rfactors_highj = {'PPR': cos(phi - phj - phk + phl)/60 + cos(phi - phj + phk - phl)/60 + cos(phi + phj - phk - phl)/10,
 'RRP': cos(phi - phj - phk + phl)/60 + cos(phi - phj + phk - phl)/60 + cos(phi + phj - phk - phl)/10,
 'PRP': cos(phi - phj - phk + phl)/60 + cos(phi - phj + phk - phl)/10 + cos(phi + phj - phk - phl)/60,
 'RPR': cos(phi - phj - phk + phl)/60 + cos(phi - phj + phk - phl)/10 + cos(phi + phj - phk - phl)/60,
 'RPP': cos(phi - phj - phk + phl)/10 + cos(phi - phj + phk - phl)/60 + cos(phi + phj - phk - phl)/60,
 'PRR': cos(phi - phj - phk + phl)/10 + cos(phi - phj + phk - phl)/60 + cos(phi + phj - phk - phl)/60,
 'PQQ': -cos(phi - phj - phk + phl)/30 + cos(phi - phj + phk - phl)/20 + cos(phi + phj - phk - phl)/20,
 'QPQ': -cos(phi - phj - phk + phl)/20 + cos(phi - phj + phk - phl)/30 - cos(phi + phj - phk - phl)/20,
 'QPR': -cos(phi - phj - phk + phl)/30 + cos(phi - phj + phk - phl)/20 + cos(phi + phj - phk - phl)/20,
 'QQP': cos(phi - phj - phk + phl)/20 + cos(phi - phj + phk - phl)/20 - cos(phi + phj - phk - phl)/30,
 'QQQ': cos(phi - phj - phk + phl)/15 + cos(phi - phj + phk - phl)/15 + cos(phi + phj - phk - phl)/15,
 'QRP': -cos(phi - phj - phk + phl)/30 + cos(phi - phj + phk - phl)/20 + cos(phi + phj - phk - phl)/20,
 'RPQ': cos(phi - phj - phk + phl)/20 + cos(phi - phj + phk - phl)/20 - cos(phi + phj - phk - phl)/30,
 'RQP': -cos(phi - phj - phk + phl)/20 + cos(phi - phj + phk - phl)/30 - cos(phi + phj - phk - phl)/20,
 'RQQ': -cos(phi - phj - phk + phl)/30 + cos(phi - phj + phk - phl)/20 + cos(phi + phj - phk - phl)/20,
 'PQR': -cos(phi - phj - phk + phl)/20 + cos(phi - phj + phk - phl)/30 - cos(phi + phj - phk - phl)/20,
 'PRQ': cos(phi - phj - phk + phl)/20 + cos(phi - phj + phk - phl)/20 - cos(phi + phj - phk - phl)/30,
 'QQR': cos(phi - phj - phk + phl)/20 + cos(phi - phj + phk - phl)/20 - cos(phi + phj - phk - phl)/30,
 'QRQ': -cos(phi - phj - phk + phl)/20 + cos(phi - phj + phk - phl)/30 - cos(phi + phj - phk - phl)/20}

rfactors_highj_dict = {'PPR': {'c0': 1, 'c12': Rational(1, 10), 'c13': Rational(1, 60), 'c14': Rational(1, 60)},
 'RRP': {'c0': 1, 'c12': Rational(1, 10), 'c13': Rational(1, 60), 'c14': Rational(1, 60)},
 'PRP': {'c0': 1, 'c12': Rational(1, 60), 'c13': Rational(1, 10), 'c14': Rational(1, 60)},
 'RPR': {'c0': 1, 'c12': Rational(1, 60), 'c13': Rational(1, 10), 'c14': Rational(1, 60)},
 'RPP': {'c0': 1, 'c12': Rational(1, 60), 'c13': Rational(1, 60), 'c14': Rational(1, 10)},
 'PRR': {'c0': 1, 'c12': Rational(1, 60), 'c13': Rational(1, 60), 'c14': Rational(1, 10)},
 'PQQ': {'c0': 1, 'c12': Rational(1, 20), 'c13': Rational(1, 20), 'c14': -Rational(1, 30)},
 'QPQ': {'c0': 1, 'c12': -Rational(1, 20), 'c13': Rational(1, 30), 'c14': -Rational(1, 20)},
 'QPR': {'c0': 1, 'c12': Rational(1, 20), 'c13': Rational(1, 20), 'c14': -Rational(1, 30)},
 'QQP': {'c0': 1, 'c12': -Rational(1, 30), 'c13': Rational(1, 20), 'c14': Rational(1, 20)},
 'QQQ': {'c0': 1, 'c12': Rational(1, 15), 'c13': Rational(1, 15), 'c14': Rational(1, 15)},
 'QRP': {'c0': 1, 'c12': Rational(1, 20), 'c13': Rational(1, 20), 'c14': -Rational(1, 30)},
 'RPQ': {'c0': 1, 'c12': -Rational(1, 30), 'c13': Rational(1, 20), 'c14': Rational(1, 20)},
 'RQP': {'c0': 1, 'c12': -Rational(1, 20), 'c13': Rational(1, 30), 'c14': -Rational(1, 20)},
 'RQQ': {'c0': 1, 'c12': Rational(1, 20), 'c13': Rational(1, 20), 'c14': -Rational(1, 30)},
 'PQR': {'c0': 1, 'c12': -Rational(1, 20), 'c13': Rational(1, 30), 'c14': -Rational(1, 20)},
 'PRQ': {'c0': 1, 'c12': -Rational(1, 30), 'c13': Rational(1, 20), 'c14': Rational(1, 20)},
 'QQR': {'c0': 1, 'c12': -Rational(1, 30), 'c13': Rational(1, 20), 'c14': Rational(1, 20)},
 'QRQ': {'c0': 1, 'c12': -Rational(1, 20), 'c13': Rational(1, 30), 'c14': -Rational(1, 20)}}

#: Cosine coefficients in R-factors associated with pathways being in a coherent state during the waiting time
magic_classes = [(Rational(1, 10), Rational(1, 60), Rational(1, 60)),
                 (Rational(1, 20), Rational(1, 20), Rational(1, -30)),
                 (Rational(1, 20), Rational(-1, 30), Rational(1, 20))]
#: Cosine coefficients in R-factors associated with pathways being in a population state during the waiting time
muggle_classes = [(Rational(1, 15), Rational(1, 15), Rational(1, 15)),
                  (Rational(1, 60), Rational(1, 10), Rational(1, 60)),
                  (Rational(1, 60), Rational(1, 60), Rational(1, 10)),
                  (Rational(1, 30), Rational(-1, 20), Rational(-1, 20))]
