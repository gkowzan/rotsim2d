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
x0, x1, x2, x3 = symbols("x0 x1 x2 x3", real=True)
