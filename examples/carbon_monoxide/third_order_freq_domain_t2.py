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

df2 = 100e9
F2 = u.wn2nu(5000.0)
N2 = np.round(F2/df2).astype(np.int)
fs2 = np.arange(-N2, N2)*df2
fs2_cm = u.nu2wn(fs2)

# * Calculate molecular coherence spectrum
root_prop = prop.Spectrum({'elevels': co_levels, 'populations': co_pops,
                             'line_params': co_params},
                            # filter=lambda kb: kb.parent.parent==pw.KetBra(0, 0, 1, 1)\
                            # or kb.parent.parent==pw.KetBra(1, 1, 0, 0),
                            filter=lambda kb: True,
                            freq_shift=[u.wn2nu(1800.0), 0.0, 0.0])

# ** With interstate coherences
# *** Third order -- nonrephasing
import time
ttt = time.time()
resp_xs_nonreph = np.zeros((fs.size, fs2.size), dtype=np.complex)
for j in range(5, 6):
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'])
    root = pw.remove_rephasing(root)
    root = pw.only_between(root, pw.KetBra(1, 4, 0, 5), pw.KetBra(1, 4, 0, 5))
    root = pw.readout(root)
    root = pw.remove_nondiagonal(root)
    # resp_xs_nonreph += co_pops[(0, j)]*root_prop.response(root, [fs[:, np.newaxis], 10.0, fs[np.newaxis, :]])
    resp_xs_nonreph += co_pops[(0, j)]*root_prop.response(root, [fs[:, np.newaxis], fs2[np.newaxis, :], None])
print(time.time()-ttt)

# *** Third order -- rephasing
resp_xs_reph = np.zeros((fs.size, fs2.size), dtype=np.complex)
for j in range(5, 6):
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'])
    root = pw.remove_nonrephasing(root)
    root = pw.only_between(root, pw.KetBra(1, 4, 0, 5), pw.KetBra(1, 4, 0, 5))
    root = pw.readout(root)
    root = pw.remove_nondiagonal(root)
    # resp_xs_reph += co_pops[(0, j)]*root_prop.response(root, [fs[:, np.newaxis], 10.0, fs[np.newaxis, :]])
    resp_xs_reph += co_pops[(0, j)]*root_prop.response(root, [fs[:, np.newaxis], fs2[np.newaxis, :], None])

# *** Plot
fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs2_cm),
                         spec2d=np.real(resp_xs_reph))
fig_dict['fig'].suptitle(r'$\vec{k}_1 = -\vec{k}_3$')
# fig_dict['ax2d'].set(xlim=(2000, 2250), ylim=(2000, 2250))
# fig_dict['fig'].savefig(str(OUTPUT / 'P5P5_rephasing.png'))

fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs2_cm),
                         spec2d=np.imag(resp_xs_nonreph[::-1]))
fig_dict['fig'].suptitle(r'$\vec{k}_1 = \vec{k}_3$')
# fig_dict['ax2d'].set(xlim=(2000, 2250), ylim=(2000, 2250))
# fig_dict['fig'].savefig(str(OUTPUT / 'P5P5_nonrephasing.png'))

fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs2_cm),
                         spec2d=prop.absorptive(resp_xs_nonreph+resp_xs_reph))
fig_dict['fig'].suptitle('purely absorptive')
# fig_dict['ax2d'].set(xlim=(2000, 2250), ylim=(2000, 2250))
# fig_dict['fig'].savefig(str(OUTPUT / 'P5P5_absorptive.png'))

fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs2_cm),
                         spec2d=prop.dispersive(resp_xs_nonreph+resp_xs_reph))
fig_dict['fig'].suptitle('purely dispersive')
# fig_dict['ax2d'].set(xlim=(2000, 2250), ylim=(2000, 2250))
# fig_dict['fig'].savefig(str(OUTPUT / 'P5P5_dispersive.png'))

fig_1d = vis.plot1d_probe(u.nu2wn(co_params[((0,5), (1,4))]['nu'])-1800.0,
                          (fs_cm, fs_cm), np.imag(resp_xs_reph), label='rephasing')
vis.plot1d_probe(u.nu2wn(co_params[((0,5), (1,4))]['nu'])-1800.0,
                 (fs_cm, fs_cm), np.imag(resp_xs_nonreph), fig_dict=fig_1d, label='non-rephasing')
vis.plot1d_probe(u.nu2wn(co_params[((0,5), (1,4))]['nu'])-1800.0,
                 (fs_cm, fs_cm), np.imag(resp_xs_reph+resp_xs_nonreph), fig_dict=fig_1d, label='absorptive')
fig_1d['ax'].legend(loc='best')
