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
co_env = {'T': Tgas, 'p': pressure}
co_params = happier.apply_env(co_params, co_env)
co_params = {k: {'mu': v['mu'], 'gam': u.wn2nu(v['gam']), 'nu': u.wn2nu(v['nu'])}
             for k, v in co_params.items()}

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
                            freq_shift=u.wn2nu(1800.0))

# ** Without interstate coherences
# *** Third order -- nonrephasing
resp_xs_nonreph = np.zeros((fs.size, fs.size), dtype=np.complex)
for j in range(0, 40):
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['both', 'both', 'both'])
    root = pw.remove_rephasing(root)
    root = pw.remove_interstates(root)
    root = pw.readout(root)
    root = pw.remove_nondiagonal(root)
    resp_xs_nonreph += co_pops[(0, j)]*root_prop.response(root, [fs[:, np.newaxis], 0.0, fs[np.newaxis, :]])

# *** Third order -- rephasing
resp_xs_reph = np.zeros((fs.size, fs.size), dtype=np.complex)
for j in range(0, 40):
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['both', 'both', 'both'])
    root = pw.remove_nonrephasing(root)
    root = pw.remove_interstates(root)
    root = pw.readout(root)
    root = pw.remove_nondiagonal(root)
    resp_xs_reph += co_pops[(0, j)]*root_prop.response(root, [fs[:, np.newaxis], 0.0, fs[np.newaxis, :]])

# *** Plot
fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=np.imag(resp_xs_reph))
fig_dict['fig'].suptitle(r'$\vec{k}_1 = -\vec{k}_3$')
fig_dict['ax2d'].set(xlim=(2000, 2250), ylim=(2000, 2250))
fig_dict['fig'].savefig(str(OUTPUT / 'P5P5_nointerstate_rephasing.png'))

fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=np.imag(resp_xs_nonreph[::-1]))
fig_dict['fig'].suptitle(r'$\vec{k}_1 = \vec{k}_3$')
fig_dict['ax2d'].set(xlim=(2000, 2250), ylim=(2000, 2250))
fig_dict['fig'].savefig(str(OUTPUT / 'P5P5_nointerstate_nonrephasing.png'))

fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=np.imag(resp_xs_reph[1:] + resp_xs_nonreph[1:][::-1]))
fig_dict['fig'].suptitle('purely absorptive')
fig_dict['ax2d'].set(xlim=(2000, 2250), ylim=(2000, 2250))
fig_dict['fig'].savefig(str(OUTPUT / 'P5P5_nointerstate_absorptive.png'))

fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=np.imag(resp_xs_reph[1:] - resp_xs_nonreph[1:][::-1]))
fig_dict['fig'].suptitle('purely dispersive')
fig_dict['ax2d'].set(xlim=(2000, 2250), ylim=(2000, 2250))
fig_dict['fig'].savefig(str(OUTPUT / 'P5P5_nointerstate_dispersive.png'))

fig_1d = vis.plot1d_probe(u.nu2wn(co_params[((0,5), (1,4))]['nu'])-1800.0,
                          (fs_cm, fs_cm), np.imag(resp_xs_reph), label='rephasing')
vis.plot1d_probe(u.nu2wn(co_params[((0,5), (1,4))]['nu'])-1800.0,
                 (fs_cm, fs_cm), np.imag(resp_xs_nonreph), fig_dict=fig_1d, label='non-rephasing')
vis.plot1d_probe(u.nu2wn(co_params[((0,5), (1,4))]['nu'])-1800.0,
                 (fs_cm, fs_cm), np.imag(resp_xs_reph+resp_xs_nonreph), fig_dict=fig_1d, label='absorptive')
fig_1d['ax'].legend(loc='best')

# ** With interstate coherences
# *** Third order -- nonrephasing
resp_xs_nonreph = np.zeros((fs.size, fs.size), dtype=np.complex)
for j in range(5, 6):
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['both', 'both', 'both'])
    root = pw.remove_rephasing(root)
    root = pw.only_between(root, pw.KetBra(1, 4, 0, 5), pw.KetBra(1, 4, 0, 5))
    root = pw.readout(root)
    root = pw.remove_nondiagonal(root)
    resp_xs_nonreph += co_pops[(0, j)]*root_prop.response(root, [fs[:, np.newaxis], 0.0, fs[np.newaxis, :]])

# *** Third order -- rephasing
resp_xs_reph = np.zeros((fs.size, fs.size), dtype=np.complex)
for j in range(5, 6):
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['both', 'both', 'both'])
    root = pw.remove_nonrephasing(root)
    root = pw.only_between(root, pw.KetBra(1, 4, 0, 5), pw.KetBra(1, 4, 0, 5))
    root = pw.readout(root)
    root = pw.remove_nondiagonal(root)
    resp_xs_reph += co_pops[(0, j)]*root_prop.response(root, [fs[:, np.newaxis], 0.0, fs[np.newaxis, :]])

# *** Plot
fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=np.imag(resp_xs_reph))
fig_dict['fig'].suptitle(r'$\vec{k}_1 = -\vec{k}_3$')
fig_dict['ax2d'].set(xlim=(2000, 2250), ylim=(2000, 2250))
fig_dict['fig'].savefig(str(OUTPUT / 'P5P5_rephasing.png'))

fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=np.imag(resp_xs_nonreph[::-1]))
fig_dict['fig'].suptitle(r'$\vec{k}_1 = \vec{k}_3$')
fig_dict['ax2d'].set(xlim=(2000, 2250), ylim=(2000, 2250))
fig_dict['fig'].savefig(str(OUTPUT / 'P5P5_nonrephasing.png'))

fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=np.imag(resp_xs_reph[1:] + resp_xs_nonreph[1:][::-1]))
fig_dict['fig'].suptitle('purely absorptive')
fig_dict['ax2d'].set(xlim=(2000, 2250), ylim=(2000, 2250))
fig_dict['fig'].savefig(str(OUTPUT / 'P5P5_absorptive.png'))

fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=np.imag(resp_xs_reph[1:] - resp_xs_nonreph[1:][::-1]))
fig_dict['fig'].suptitle('purely dispersive')
fig_dict['ax2d'].set(xlim=(2000, 2250), ylim=(2000, 2250))
fig_dict['fig'].savefig(str(OUTPUT / 'P5P5_dispersive.png'))

fig_1d = vis.plot1d_probe(u.nu2wn(co_params[((0,5), (1,4))]['nu'])-1800.0,
                          (fs_cm, fs_cm), np.imag(resp_xs_reph), label='rephasing')
vis.plot1d_probe(u.nu2wn(co_params[((0,5), (1,4))]['nu'])-1800.0,
                 (fs_cm, fs_cm), np.imag(resp_xs_nonreph), fig_dict=fig_1d, label='non-rephasing')
vis.plot1d_probe(u.nu2wn(co_params[((0,5), (1,4))]['nu'])-1800.0,
                 (fs_cm, fs_cm), np.imag(resp_xs_reph+resp_xs_nonreph), fig_dict=fig_1d, label='absorptive')
fig_1d['ax'].legend(loc='best')
