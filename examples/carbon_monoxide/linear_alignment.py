"""Show evolution of molecular alignment after second-order interaction, i.e.
during waiting time or between pump and probe in pump-probe scheme.

TODO:
- Calculate Zare's orientation and alignment parameters.
- Calculate these parameters for all DM elements after second-order excitation.
  Compare with some known results.
"""
import numpy as np
import rotsim2d.couple as cp
from gkpywigxjpf import table_init, temp_init, table_free, temp_free
from sympy import symbols, integrate, cos, pi, sin, simplify
from sympy.functions.special.spherical_harmonics import Ynm, Ynm_c

table_init(800, 3)
temp_init(800)
# * Test cos2_jm
j = 40
cos_val = 0.0
for m in range(-j, j+1):
    cos_val += cp.cos2_jm(j, m, j, m)

# * Tear down
temp_free()
table_free()

# * Sympy
r, phi, theta = symbols("r, phi, theta", real=True)
val = integrate(cos(theta)**2*sin(theta), (theta, 0, pi))

val = integrate(integrate(simplify(Ynm_c(2, 0, theta, phi).expand(func=True)*
                                   Ynm(0, 0, theta, phi).expand(func=True)*
                                   sin(theta)*
                                   cos(theta)**2), (theta, 0, pi)), (phi, 0, 2*pi))
print(np.arccos(np.sqrt(np.float(val.evalf())))*180.0/np.pi)
