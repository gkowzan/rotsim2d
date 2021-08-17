"""Common definitions for :mod:`rotsim2d.symbolic.functions` and
:mod:`rotsim2d.symbolic.results`."""
from sympy import *
#: Initial angular momentum
J_i = symbols("J_i", integer=True, nonnegative=True)
#: Dummy angles
phi, phj, phk, phl = symbols(r"\phi_i \phi_j \phi_k \phi_l", real=True)
thetas = symbols(r"\theta_i \theta_j \theta_k \theta_l", real=True)
thetas_dict = dict(zip(('omg1', 'omg2', 'omg3', 'mu'),
                       thetas))
#: Pulse and detection angles
theta_i, theta_j, theta_k, theta_l = thetas
#: Dummy variables for derivations
x0, x1, x2, x3, x4, x1x2, x3x4 = symbols("x0 x1 x2 x3 x4 x1x2 x3x4", real=True)
#: Rotational constants
nu0, nu1, nu2, B0, B1, B2 = symbols(r"\nu_0 \nu_1 \nu_2 B_0 B_1 B_2", real=True, nonnegative=True)
