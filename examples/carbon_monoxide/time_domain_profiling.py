"""Time-domain evolution examples."""
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

mb_per_float = 8/2**20
mb_per_complex = 16/2**20

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

# resp_xs_spec = np.zeros(ts.size, dtype=np.complex)
# for j in range(50):
#     root = pw.KetBra(0, j, 0, j)
#     root = pw.excite(root, 'omg1', 'both')
#     root = pw.readout(root)
#     pop_factor = co_pops[(0, j)]/(2*j+1)
#     resp_xs_spec += pop_factor*root_prop.cross_section_spectrum(root, [ts])
#     # resp += pop_factor*root_prop.response(root, [ts])
#     # resp_xs += pop_factor*root_prop.cross_section(root, [ts])

# ** MultiPropagator
pathways = [pw.readout(pw.excite(pw.KetBra(0, j, 0, j, pop=co_pops[(0, j)]/(2*j+1)), 'omg1', 'both')) for j in range(50)]
multi_prop = prop.MultiPropagator(pathways, root_prop)
resp_xs_spec = multi_prop.cross_section_spectrum([ts])

# import timeit
# print(timeit.timeit('multi_prop.cross_section_spectrum([ts])', globals=globals(), number=5))
