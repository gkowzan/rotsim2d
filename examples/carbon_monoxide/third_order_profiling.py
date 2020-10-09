# * Imports and functions
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from spectroscopy import happier
import shed.units as u
import rotsim2d.propagate as prop
import rotsim2d.pathways as pw
import rotsim2d.visual as vis
plt.ion()

OUTPUT = Path('/mnt/d/DFCS/Stony Brook/rotsim2d_results/carbon_monoxide')
# * System properties
pressure, Tgas, length = 5.0, 296.0, 3.0
conc = 1e-7*happier.volumeConcentration(pressure, Tgas)*1e6
co_levels = happier.energy_levels([0,1,2], 5, 1)
co_pops = happier.equilibrium_pops(co_levels, Tgas, 5, 1)
co_params = happier.line_params([(0, 1), (1, 2), (2, 3), (0, 2)], 5, 1)
co_params_suppl = happier.generate_line_params(
    [(0, 0, 0), (1, 1, 0), (2, 2, 0), (0, 0, 2), (1, 1, 2), (0, 2, 2), (0, 2, 0)], co_levels, co_params)
co_params.update(co_params_suppl)
co_env = {'T': Tgas, 'p': pressure}
co_params = happier.apply_env(co_params, co_env, to_nu=True)

# * Response function parameters
df = 10e9
F = u.wn2nu(700.0)
N = np.round(F/df).astype(np.int)
fs = np.arange(-N, N)*df
fs_cm = u.nu2wn(fs)

# * Calculate molecular coherence spectrum
root_prop = prop.Spectrum({'elevels': co_levels, 'populations': co_pops,
                             'line_params': co_params},
                            # filter=lambda kb: kb.parent.parent==pw.KetBra(0, 0, 1, 1)\
                            # or kb.parent.parent==pw.KetBra(1, 1, 0, 0),
                            filter=lambda kb: True,
                            freq_shift=[u.wn2nu(1800.0), 0.0, u.wn2nu(1800.0)])

# ** With interstate coherences
# *** Third order -- nonrephasing
import time
ttt = time.time()
# resp_xs_nonreph = np.zeros((fs.size, fs.size), dtype=np.complex)
for j in range(5, 6):
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'])
    root = pw.remove_rephasing(root)
    root = pw.only_between(root, pw.KetBra(1, 4, 0, 5), pw.KetBra(1, 4, 0, 5))
    root = pw.readout(root)
    root = pw.remove_nondiagonal(root)
    resp_xs_nonreph = co_pops[(0, j)]*root_prop.response(root, [fs[:, np.newaxis], 10.0, fs[np.newaxis, :]])
print(time.time() - ttt)
