"""Verify that my ESA line strengths are correct.

Line strengths check out.
"""
# * Imports
import numpy as np
import scipy.constants as C
import spectroscopy.happier as hap
import rotsim2d.propagate as prop
import shed.units as u

# * HITRAN line strengths
line01_0_1 = 9.48e-20
line01_5_6 = 4.34e-19
line12_0_1 = 5.586e-24
line12_5_6 = 2.563e-23

# * dipole matrix elements
co_levels = hap.energy_levels([0,1], 5, 1)
co_params = hap.line_params([(0,1), (1,2)], 5, 1)
co_pops = hap.equilibrium_pops(co_levels, 296.0, 5, 1)

def resp(mu):
    return 1.0j/C.hbar*mu**2

def xs(resp, nu):
    return resp*np.pi*nu/C.c/C.epsilon_0/6.0


my_line01_0_1 = xs(resp(co_params[((0,0), (1,1))]['mu']), u.wn2nu(co_params[((0,0), (1,1))]['nu']))*co_pops[(0,0)]
my_line01_5_6 = xs(resp(co_params[((0,5), (1,6))]['mu']), u.wn2nu(co_params[((0,5), (1,6))]['nu']))*co_pops[(0,5)]/(2*5+1)
my_line12_0_1 = xs(resp(co_params[((1,0), (2,1))]['mu']), u.wn2nu(co_params[((1,0), (2,1))]['nu']))*co_pops[(1,0)]
my_line12_5_6 = xs(resp(co_params[((1,5), (2,6))]['mu']), u.wn2nu(co_params[((1,5), (2,6))]['nu']))*co_pops[(1,5)]/(2*5+1)
