# * Imports and constants
from pathlib import Path
from pprint import pprint

import matplotlib.gridspec as grd
import matplotlib.pyplot as plt
import numpy as np
import rotsim2d.couple as cp
import rotsim2d.pathways as pw
import rotsim2d.dressedleaf as dl
import rotsim2d.visual as vis
import knickknacks.units as u
from gkpywigxjpf import wigxjpf
from molspecutils import happier
plt.ion()

OUTPUT = Path('/mnt/d/DFCS/Stony Brook/rotsim2d_results/carbon_monoxide')
j = 5
ma = 54.7356103172453*np.pi/180.0 # magic angle

# * System parameters
pressure, Tgas, length = 15.0, 296.0, 3.0
conc = 1e-7*happier.volumeConcentration(pressure, Tgas)*1e6
co_levels = happier.energy_levels([0,1,2], 5, 1)
co_pops = happier.equilibrium_pops(co_levels, Tgas, 5, 1)
co_params = happier.line_params([(0, 1), (1, 2), (2, 3), (0, 2), (0, 0), (1, 1), (2, 2)], 5, 1)
co_params_suppl = happier.generate_line_params(
    [(0, 0, 0), (1, 1, 0), (2, 2, 0), (0, 0, 2), (1, 1, 2), (0, 2, 2), (0, 2, 0)], co_levels, co_params)
co_params_suppl2 = happier.generate_line_params([(0, 1, 0), (1, 2, 0)], co_levels, co_params)
co_params.update(co_params_suppl)
co_params.update(co_params_suppl2)
co_env = {'T': Tgas, 'p': pressure}
co_params = happier.apply_env(co_params, co_env, to_nu=True)
sys_params = dict(elevels=co_levels, populations=co_pops, line_params=co_params)

# * Single non-rephasing
# ** Linear
j = 5
nonreph = pw.KetBra(0,j, 0,j)
nonreph = pw.multi_excite(nonreph, ['omg1', 'omg2', 'omg3'],
                          parts=['ket', 'both', 'both'],
                          light_angles=[1, 2, 3])
nonreph = pw.remove_rephasing(nonreph)
nonreph = pw.readout(nonreph, 4)
dls_nonreph = [dl.DressedLeaf(l, sys_params) for l in nonreph.leaves]
dls_by_js = dl.split_by_js(dls_nonreph)
linear_deg_js = list(dls_by_js.keys())
linear_js = dl.undegenerate_js(linear_deg_js)

# ** Symmetric top
j = 5
nonreph = pw.KetBra(0,j, 0,j)
nonreph = pw.multi_excite(nonreph, ['omg1', 'omg2', 'omg3'],
                          parts=['ket', 'both', 'both'],
                          light_angles=[1, 2, 3], rotor='symmetric')
nonreph = pw.remove_rephasing(nonreph)
nonreph = pw.readout(nonreph, 4, rotor='symmetric')
dls_nonreph = [dl.DressedLeaf(l, sys_params) for l in nonreph.leaves]
dls_by_js = dl.split_by_js(dls_nonreph)
symmetric_deg_js = list(dls_by_js.keys())
symmetric_js = dl.undegenerate_js(symmetric_deg_js)


n = len(dls_nonreph)
prs, pus, sigs = np.empty(n), np.empty(n), np.empty(n)
with wigxjpf(300, 6):
    for i in range(n):
        prs[i] = u.nu2wn(dls_nonreph[i].nus[2])
        pus[i] = u.nu2wn(dls_nonreph[i].nus[0])
        # sigs[i] = np.imag(dls_nonreph[i].const)
        sigs[i] = np.imag(dls_nonreph[i].intensity())

fig, ax = plt.subplots(figsize=(7, 3))
ax.scatter(prs, pus, s=20.0, c=-sigs, cmap='seismic', vmin=-np.abs(np.max(sigs)), vmax=np.abs(np.max(sigs)))
ax.set(xlabel=r'Probe (cm$^{-1}$)', ylabel=r'Pump (cm$^{-1}$)')

for dl in dls_nonreph:
    # ax.text(u.nu2wn(dl.nus[2]), u.nu2wn(dl.nus[0]), '\n'.join(dl.peak).replace("><", "-").replace("|", ""), ha='center', va='top', fontsize=8)
    if '|1,4>' in dl.peak[0]:
        align = 'bottom'
    else:
        align = 'top'
    ax.text(u.nu2wn(dl.nus[2]), u.nu2wn(dl.nus[0]), '\n'.join(dl.peak), ha='center', va=align, fontsize=10)
