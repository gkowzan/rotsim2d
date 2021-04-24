"""Time-domain evolution examples.

They look wrong.
"""
# * Imports and functions
import numpy as np
import scipy.constants as C
import matplotlib.pyplot as plt
import molspecutils.foreign.hapi3 as h3
from molspecutils import happier
import knickknacks.units as u
import rotsim2d.propagate as prop
import rotsim2d.pathways as pw
plt.ion()

# * System properties
pressure, Tgas, length = 5.0, 296.0, 3.0
conc = 1e-7*happier.volumeConcentration(pressure, Tgas)*1e6
co_levels = happier.energy_levels([0,1,2], 5, 1)
co_pops = happier.equilibrium_pops(co_levels, Tgas, 5, 1)
co_params = happier.line_params([(0, 1)], 5, 1)
co_env = {'T': Tgas, 'p': pressure}
co_params = happier.apply_env(co_params, co_env)
co_params = {k: {'mu': v['mu'], 'gam': u.wn2nu(v['gam']), 'nu': u.wn2nu(v['nu'])}
             for k, v in co_params.items()}

# * Response function parameters
dt = 1/u.wn2nu(2200.0)/10.0
T = 1/co_params[((0,0), (1,1))]['gam']*9         # relaxation time constant
# dt = 1e-15
# T = 1e-12
N = np.round(T/dt).astype(np.int)
ts = np.arange(N)*dt

# * Calculate molecular coherence evolution
root_prop = prop.Propagator({'elevels': co_levels, 'populations': co_pops,
                             'line_params': co_params},
                            # filter=lambda kb: kb.parent.parent==pw.KetBra(0, 0, 1, 1)\
                            # or kb.parent.parent==pw.KetBra(1, 1, 0, 0),
                            filter=lambda kb: True)

resp, resp_xs = np.zeros(ts.size, dtype=np.complex), np.zeros(ts.size, dtype=np.complex)
for j in range(50):
    root = pw.KetBra(0, j, 0, j)
    root = pw.excite(root, 'omg1', 'both')
    root = pw.readout(root)
    pop_factor = co_pops[(0, j)]
    resp += pop_factor*root_prop.response(root, [ts])
    resp_xs += pop_factor*root_prop.cross_section(root, [ts])

# ** Cross section (Hz)
fig, ax = plt.subplots()
ax.plot(ts, resp_xs)
ax.set(xlabel='Time (s)', ylabel=r'Optical density (m$^2$ Hz)')

# ** Cross section (cm-1)
fig, ax = plt.subplots()
ax.plot(ts, resp_xs*prop.xs2cm)
ax.set(xlabel='Time (s)', ylabel=r'Optical density (cm$^2$ cm$^{-1}$)')

# ** Optical density
fig, ax = plt.subplots()
ax.plot(ts, resp_xs*conc*length)
ax.set(xlabel='Time (s)', ylabel=r'Optical density (1/s)')

# ** Polarization density
# For frequency comb pulses, frep=100 MHz, beam radius=50 um, avg. power=100 mW.
# Pulses are Dirac deltas.
Pt, r = 100e-3/100e6, 50e-6
efield0 = prop.power2electric_field(Pt, r) #: source electric field

fig, ax = plt.subplots()
ax.plot(ts, prop.polarization(efield0, resp, conc))
ax.set(xlabel='Time (s)', ylabel=r'Polarization density (C/m$^2$)')

# ** Generated electric field
resp_efield = prop.electric_field(efield0, resp_xs, conc *length)

fig, ax = plt.subplots()
ax.plot(ts, np.imag(resp_efield)/efield0)
ax.set(xlabel='Time (s)', ylabel=r'Electric field (N/C)')

# ** Generated beam power
fig, ax = plt.subplots()
ax.plot(ts, prop.electric_field2power(resp_efield, r))
ax.set(xlabel='Time (s)', ylabel=r'Power (W)')
