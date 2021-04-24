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
# co_params_suppl = happier.generate_line_params(
#     [(0, 0, 0), (1, 1, 0), (2, 2, 0)], co_levels, co_params)
co_params.update(co_params_suppl)
co_params.update(co_params_suppl2)
co_env = {'T': Tgas, 'p': pressure}
co_params = happier.apply_env(co_params, co_env, to_nu=True)
sys_params = dict(elevels=co_levels, populations=co_pops, line_params=co_params)

# * Pathways
# ** Plot multiple non-rephasing
sopr = -36.87/180.0*np.pi
pprr = 26.57/180.0*np.pi
# pws = pw.gen_pathways(20, [np.pi/2, pprr, np.pi/2, np.pi/4], co_pops, [pw.remove_rephasing, pw.remove_interstates])
# pws = pw.gen_pathways(range(1, 20), [np.pi/2, sopr, np.pi/2, np.pi/4],
#                       co_pops, [pw.remove_rephasing, pw.remove_interstates])
# pws = pw.gen_pathways(20, [0.0, np.arctan(2/np.tan(np.pi/4)), 0.0, np.pi/4], co_pops, [pw.remove_rephasing, pw.remove_interstates])
pws = pw.gen_pathways(range(20), [0.0]*4, co_pops, [pw.remove_rephasing])
ll = dl.dress_pws(pws, sys_params)
pl = dl.peak_list(ll)

# ** Plot
vis.plot2d_scatter(pl, line=True)

# ** Plot single non-rephasing
# [pi+phj, ma, 0, 0] <- [0, pi+ma, 0, ma]
j = 5
nonreph = pw.KetBra(0,j, 0,j)
# nonreph = pw.multi_excite(nonreph, ['omg1', 'omg2', 'omg3'],
#                           parts=['ket', 'both', 'both'],
#                           light_angles=[0.0, np.pi+ma, 0.0])
nonreph = pw.multi_excite(nonreph, ['omg1', 'omg2', 'omg3'],
                          parts=['ket', 'both', 'both'],
                          light_angles=[1, 2, 1])
nonreph = pw.remove_rephasing(nonreph)
# nonreph = pw.remove_interstates(nonreph)
# nonreph = pw.readout(nonreph, ma)
nonreph = pw.readout(nonreph, 4)
dls_nonreph = [dl.DressedLeaf(l, sys_params) for l in nonreph.leaves]
dls_by_js = dl.split_by_js(dls_nonreph)
dls_by_peaks = pw.split_by_peaks(dls_nonreph)
dls_by_pols = pw.split_by_pols(dls_nonreph)
dls_by_pols_js = pw.split_by_pols_js(dls_nonreph)
dls_by_highjs_pols = pw.split_by_pols_highjs(dls_nonreph)

jsets = list(dls_by_js.keys())
jsets.sort()
jsets_nondeg = []
for jset in jsets:
    if jset not in jsets_nondeg and dl.perm_js(jset) not in jsets_nondeg:
        jsets_nondeg.append(jset)

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
