# * Imports and functions
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import rotsim2d.pathways as pw
import rotsim2d.propagate as prop
import rotsim2d.visual as vis
import knickknacks.units as u
from molspecutils import happier

prop.ignore_missing = False
plt.ion()

OUTPUT = Path('/mnt/d/DFCS/Stony Brook/rotsim2d_results/carbon_monoxide')
# * System properties
pressure, Tgas, length = 20.0, 296.0, 3.0
conc = 1e-7*happier.volumeConcentration(pressure, Tgas)*1e6
co_levels = happier.energy_levels([0,1,2], 5, 1)
co_pops = happier.equilibrium_pops(co_levels, Tgas, 5, 1)
co_params = happier.line_params([(0, 1), (1, 2), (2, 3), (0, 2)], 5, 1)
co_params_suppl = happier.generate_line_params(
    [(0, 0, 0), (1, 1, 0), (2, 2, 0), (0, 0, 2), (1, 1, 2), (0, 2, 2), (0, 2, 0)], co_levels, co_params)
# co_params_suppl = happier.generate_line_params(
#     [(0, 0, 0), (1, 1, 0), (2, 2, 0)], co_levels, co_params)
co_params.update(co_params_suppl)
co_env = {'T': Tgas, 'p': pressure}
co_params = happier.apply_env(co_params, co_env, to_nu=True)

# * Response function parameters
df = 10e9
F = u.wn2nu(700.0)
N = np.round(F/df).astype(np.int)
fs = np.arange(-N, N+1)*df
fs_cm = u.nu2wn(fs)

# * Calculate molecular coherence spectrum
root_prop = prop.Spectrum({'elevels': co_levels, 'populations': co_pops,
                             'line_params': co_params},
                            # filter=lambda kb: kb.parent.parent==pw.KetBra(0, 0, 1, 1)\
                            # or kb.parent.parent==pw.KetBra(1, 1, 0, 0),
                            filter=lambda kb: True,
                            freq_shift=[u.wn2nu(1800.0), 0.0, u.wn2nu(1800.0)])

# ** Third order --- full
pws = []
for j in range(0, 16):
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'])
    root = pw.readout(root)
    root = pw.remove_nondiagonal(root)
    root.pop = co_pops[(0, j)]
    pws.append(root)
multi_prop = prop.MultiPropagator(pws, root_prop)
resp_xs_full = multi_prop.response([fs[:, np.newaxis], None, fs[np.newaxis, :]])

# *** Plot
fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=prop.absorptive(resp_xs_full),
                         scale='symlog')
fig_dict['fig'].canvas.set_window_title('Frequency domain -- absorptive -- all pathways')

# vis.latex_compile(OUTPUT / 'diagrams/test_diagrams.tex',
#                   vis.LATEX_PRE + vis.tikz_diagrams(pws[5]) + vis.LATEX_POST)

# ** Third order --- no interstates
pws = []
for j in range(0, 16):
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'])
    root = pw.remove_interstates(root)
    root = pw.readout(root)
    root = pw.remove_nondiagonal(root)
    root.pop = co_pops[(0, j)]
    pws.append(root)
multi_prop = prop.MultiPropagator(pws, root_prop)
resp_xs_noints = multi_prop.response([fs[:, np.newaxis], None, fs[np.newaxis, :]])

# *** Plot
fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=prop.absorptive(resp_xs_noints),
                         scale='linear')
fig_dict['fig'].canvas.set_window_title('Frequency domain -- absorptive -- no interstates')


# ** Third order --- no overtones
pws = []
for j in range(0, 16):
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'])
    root = pw.remove_overtones(root)
    root = pw.readout(root)
    root = pw.remove_nondiagonal(root)
    root.pop = co_pops[(0, j)]
    pws.append(root)
multi_prop = prop.MultiPropagator(pws, root_prop)
resp_xs_noover = multi_prop.response([fs[:, np.newaxis], None, fs[np.newaxis, :]])

# *** Plot
fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=prop.absorptive(resp_xs_noover),
                         scale='linear')
fig_dict['fig'].canvas.set_window_title('Frequency domain -- absorptive -- no overtones')
