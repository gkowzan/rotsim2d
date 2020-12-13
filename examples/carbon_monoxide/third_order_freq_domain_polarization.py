"""Investigate polarization dependence."""
# * Imports
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from spectroscopy import happier
import shed.units as u
import rotsim2d.propagate as prop
prop.ignore_missing = False
import rotsim2d.pathways as pw
import rotsim2d.visual as vis
import pywigxjpf as wig
plt.ion()

OUTPUT = Path('/mnt/d/DFCS/Stony Brook/rotsim2d_results/carbon_monoxide/polarization')
# * System properties
pressure, Tgas, length = 15.0, 296.0, 3.0
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

# * Calculate molecular coherence spectrum -- full
root_prop = prop.Spectrum({'elevels': co_levels, 'populations': co_pops,
                             'line_params': co_params},
                            # filter=lambda kb: kb.parent.parent==pw.KetBra(0, 0, 1, 1)\
                            # or kb.parent.parent==pw.KetBra(1, 1, 0, 0),
                            filter=lambda kb: True,
                            freq_shift=[u.wn2nu(1800.0), 0.0, u.wn2nu(1800.0)])
jmax = 20

# ** XXXX
pws = []
for j in range(0, jmax):
    root = pw.KetBra(0, j, 0, j)

    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'],
                           # light_angles=[54.7356/180.0*np.pi, 0.0, 54.7356/180.0*np.pi]
                           # light_angles=[np.pi/2, 0.0, np.pi/2]
                           )
    root = pw.readout(root)
    root = pw.remove_nondiagonal(root)
    root.pop = co_pops[(0, j)]/np.sqrt(2*j+1)
    pws.append(root)
multi_prop = prop.MultiPropagator(pws, root_prop)
with wig.wigxjpf(300, 6):
    resp_xs_xxxx = multi_prop.response([fs[:, np.newaxis], None, fs[np.newaxis, :]])

# ** XYXY
pws = []
for j in range(0, jmax):
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'],
                           # light_angles=[54.7356/180.0*np.pi, 0.0, 54.7356/180.0*np.pi]
                           light_angles=[np.pi/2, 0.0, np.pi/2]
                           )
    root = pw.readout(root)
    root = pw.remove_nondiagonal(root)
    root.pop = co_pops[(0, j)]/np.sqrt(2*j+1)
    pws.append(root)
multi_prop = prop.MultiPropagator(pws, root_prop)
with wig.wigxjpf(300, 6):
    resp_xs_xyxy = multi_prop.response([fs[:, np.newaxis], None, fs[np.newaxis, :]])

# ** X(MA)X(MA)
pws = []
for j in range(0, jmax):
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'],
                           light_angles=[54.7356/180.0*np.pi, 0.0, 54.7356/180.0*np.pi]
                           # light_angles=[np.pi/2, 0.0, np.pi/2]
                           )
    root = pw.readout(root)
    root = pw.remove_nondiagonal(root)
    root.pop = co_pops[(0, j)]/np.sqrt(2*j+1)
    pws.append(root)
multi_prop = prop.MultiPropagator(pws, root_prop)
with wig.wigxjpf(300, 6):
    resp_xs_xmaxma = multi_prop.response([fs[:, np.newaxis], None, fs[np.newaxis, :]])

# ** Bracamonte
pws = []
for j in range(0, jmax):
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'],
                           light_angles=[np.pi/2, np.pi/4, np.pi/2]
                           # light_angles=[np.pi/2, 0.0, np.pi/2]
                           )
    root = pw.readout(root, -18.43/180.0*np.pi)
    root = pw.remove_nondiagonal(root)
    root.pop = co_pops[(0, j)]/np.sqrt(2*j+1)
    pws.append(root)
multi_prop = prop.MultiPropagator(pws, root_prop)
with wig.wigxjpf(300, 6):
    resp_xs_bracamonte = multi_prop.response([fs[:, np.newaxis], None, fs[np.newaxis, :]])

# ** Bracamonte 2
pws = []
for j in range(0, jmax):
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'],
                           light_angles=[np.pi/2, np.pi/4, np.pi/2]
                           # light_angles=[np.pi/2, 0.0, np.pi/2]
                           )
    root = pw.readout(root, 26.57/180.0*np.pi)
    root = pw.remove_nondiagonal(root)
    root.pop = co_pops[(0, j)]/np.sqrt(2*j+1)
    pws.append(root)
multi_prop = prop.MultiPropagator(pws, root_prop)
with wig.wigxjpf(300, 6):
    resp_xs_bracamonte2 = multi_prop.response([fs[:, np.newaxis], None, fs[np.newaxis, :]])

# ** Mine
pws = []
for j in range(0, jmax):
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'],
                           light_angles=[np.pi/2, np.pi/4, np.pi/2]
                           # light_angles=[np.pi/2, 0.0, np.pi/2]
                           )
    root = pw.readout(root, np.pi/4)
    root = pw.remove_nondiagonal(root)
    root.pop = co_pops[(0, j)]/np.sqrt(2*j+1)
    pws.append(root)
multi_prop = prop.MultiPropagator(pws, root_prop)
with wig.wigxjpf(300, 6):
    resp_xs_mine = multi_prop.response([fs[:, np.newaxis], None, fs[np.newaxis, :]])

fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=prop.absorptive(resp_xs_bracamonte2),
                         scale='linear')
fig_dict['fig'].suptitle('Bracamonte 2')
fig_dict['ax2d'].set(xlim=(2000, 2250), ylim=(2000, 2250))

fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=prop.absorptive(resp_xs_xxxx),
                         scale='linear', line=False,
                         fig_kwargs=dict(figsize=(4,3),
                                         constrained_layout=True))
fig_dict['fig'].suptitle('XXXX')
fig_dict['ax2d'].set(xlim=(2030, 2220), ylim=(2030, 2220))
fig_dict['fig'].savefig(str(OUTPUT / 'XXXX.png'))
fig_dict['fig'].savefig(str(OUTPUT / 'XXXX.pdf'))

fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=prop.absorptive(resp_xs_bracamonte),
                         scale='linear')
fig_dict['fig'].suptitle('Bracamonte')
fig_dict['ax2d'].set(xlim=(2000, 2250), ylim=(2000, 2250))

fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=prop.absorptive(resp_xs_xyxy),
                         scale='linear', line=False,
                         fig_kwargs=dict(figsize=(4,3),
                                         constrained_layout=True))
fig_dict['fig'].suptitle('XYXY')
fig_dict['ax2d'].set(xlim=(2030, 2220), ylim=(2030, 2220))
fig_dict['fig'].savefig(str(OUTPUT / 'XYXY.png'))
fig_dict['fig'].savefig(str(OUTPUT / 'XYXY.pdf'))

fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=prop.absorptive(resp_xs_xmaxma),
                         scale='linear', line=False,
                         fig_kwargs=dict(figsize=(4,3),
                                         constrained_layout=True))
fig_dict['fig'].suptitle('X(MA)X(MA)')
fig_dict['ax2d'].set(xlim=(2030, 2220), ylim=(2030, 2220))
fig_dict['fig'].savefig(str(OUTPUT / 'XMAXMA.png'))
fig_dict['fig'].savefig(str(OUTPUT / 'XMAXMA.pdf'))

fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=prop.absorptive(resp_xs_xmaxma)-prop.absorptive(resp_xs_xxxx),
                         scale='linear', line=False,
                         fig_kwargs=dict(figsize=(4,3),
                                         constrained_layout=True))
fig_dict['fig'].suptitle('XXXX+X(MA)X(MA)')
fig_dict['ax2d'].set(xlim=(2030, 2220), ylim=(2030, 2220))
fig_dict['fig'].savefig(str(OUTPUT / 'XXXX_XMAXMA.png'))
fig_dict['fig'].savefig(str(OUTPUT / 'XXXX_XMAXMA.pdf'))

# * Calculate molecular coherence spectrum -- no interstates
root_prop = prop.Spectrum({'elevels': co_levels, 'populations': co_pops,
                             'line_params': co_params},
                            # filter=lambda kb: kb.parent.parent==pw.KetBra(0, 0, 1, 1)\
                            # or kb.parent.parent==pw.KetBra(1, 1, 0, 0),
                            filter=lambda kb: True,
                            freq_shift=[u.wn2nu(1800.0), 0.0, u.wn2nu(1800.0)])
jmax = 20

# ** XXXX
pws = []
for j in range(0, jmax):
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'],
                           # light_angles=[54.7356/180.0*np.pi, 0.0, 54.7356/180.0*np.pi]
                           # light_angles=[np.pi/2, 0.0, np.pi/2]
                           )
    # root = pw.remove_interstates(root)
    root = pw.remove_overtones(root)
    root = pw.readout(root)
    root = pw.remove_nondiagonal(root)
    root.pop = co_pops[(0, j)]/np.sqrt(2*j+1)
    pws.append(root)
multi_prop = prop.MultiPropagator(pws, root_prop)
with wig.wigxjpf(300, 6):
    resp_xs_xxxx = multi_prop.response([fs[:, np.newaxis], None, fs[np.newaxis, :]])

# ** XYXY
pws = []
for j in range(0, jmax):
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'],
                           # light_angles=[54.7356/180.0*np.pi, 0.0, 54.7356/180.0*np.pi]
                           light_angles=[np.pi/2, 0.0, np.pi/2]
                           )
    # root = pw.remove_interstates(root)
    root = pw.remove_overtones(root)
    root = pw.readout(root)
    root = pw.remove_nondiagonal(root)
    root.pop = co_pops[(0, j)]/np.sqrt(2*j+1)
    pws.append(root)
multi_prop = prop.MultiPropagator(pws, root_prop)
with wig.wigxjpf(300, 6):
    resp_xs_xyxy = multi_prop.response([fs[:, np.newaxis], None, fs[np.newaxis, :]])

# ** X(MA)X(MA)
pws = []
for j in range(0, jmax):
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'],
                           light_angles=[54.7356/180.0*np.pi, 0.0, 54.7356/180.0*np.pi]
                           # light_angles=[np.pi/2, 0.0, np.pi/2]
                           )
    # root = pw.remove_interstates(root)
    root = pw.remove_overtones(root)
    root = pw.readout(root)
    root = pw.remove_nondiagonal(root)
    root.pop = co_pops[(0, j)]/np.sqrt(2*j+1)
    pws.append(root)
multi_prop = prop.MultiPropagator(pws, root_prop)
with wig.wigxjpf(300, 6):
    resp_xs_xmaxma = multi_prop.response([fs[:, np.newaxis], None, fs[np.newaxis, :]])

# ** Bracamonte
pws = []
for j in range(0, jmax):
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'],
                           light_angles=[0.0, np.pi/4, -np.pi/4]
                           # light_angles=[np.pi/2, 0.0, np.pi/2]
                           )
    root = pw.remove_overtones(root)
    root = pw.readout(root, 0.0)
    root = pw.remove_nondiagonal(root)
    root.pop = co_pops[(0, j)]/np.sqrt(2*j+1)
    pws.append(root)
multi_prop = prop.MultiPropagator(pws, root_prop)
with wig.wigxjpf(300, 6):
    resp_xs_bracamonte = multi_prop.response([fs[:, np.newaxis], None, fs[np.newaxis, :]])

# ** Bracamonte 2
pws = []
for j in range(0, jmax):
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'],
                           light_angles=[54.7356/180.0*np.pi, 54.7356/180.0*np.pi, 0.0]
                           # light_angles=[np.pi/2, 0.0, np.pi/2]
                           )
    root = pw.remove_overtones(root)
    root = pw.readout(root, 0.0)
    root = pw.remove_nondiagonal(root)
    root.pop = co_pops[(0, j)]/np.sqrt(2*j+1)
    pws.append(root)
multi_prop = prop.MultiPropagator(pws, root_prop)
with wig.wigxjpf(300, 6):
    resp_xs_bracamonte2 = multi_prop.response([fs[:, np.newaxis], None, fs[np.newaxis, :]])

# ** Mine
pws = []
for j in range(0, jmax):
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'],
                           light_angles=[np.pi/2, np.pi/4, np.pi/2]
                           # light_angles=[np.pi/2, 0.0, np.pi/2]
                           )
    root = pw.remove_interstates(root)
    root = pw.readout(root, np.pi/4)
    root = pw.remove_nondiagonal(root)
    root.pop = co_pops[(0, j)]/np.sqrt(2*j+1)
    pws.append(root)
multi_prop = prop.MultiPropagator(pws, root_prop)
with wig.wigxjpf(300, 6):
    resp_xs_mine = multi_prop.response([fs[:, np.newaxis], None, fs[np.newaxis, :]])

fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=prop.absorptive(resp_xs_bracamonte2),
                         scale='linear', line=False, fig_kwargs=dict(figsize=(4,3)))
fig_dict['fig'].suptitle('Bracamonte 2')
fig_dict['ax2d'].set(xlim=(2000, 2250), ylim=(2000, 2250))

fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=prop.absorptive(resp_xs_xxxx),
                         scale='linear', line=False,
                         fig_kwargs=dict(figsize=(4,3),
                                         constrained_layout=True))
fig_dict['fig'].suptitle('XXXX (no overtones)')
fig_dict['ax2d'].set(xlim=(2030, 2220), ylim=(2030, 2220))
fig_dict['fig'].savefig(str(OUTPUT / 'XXXX_noovertones.png'))
fig_dict['fig'].savefig(str(OUTPUT / 'XXXX_noovertones.pdf'))

fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=prop.absorptive(resp_xs_bracamonte),
                         scale='linear', line=False,
                         fig_kwargs=dict(figsize=(4,3),
                                         constrained_layout=True))
fig_dict['fig'].suptitle(r'X($\pi/4$)(-$\pi/4$)X (no overtones)')
fig_dict['ax2d'].set(xlim=(2030, 2220), ylim=(2030, 2220))
fig_dict['fig'].savefig(str(OUTPUT / 'XPI4PI4X_noovertones.png'))
fig_dict['fig'].savefig(str(OUTPUT / 'XPI4PI4X_noovertones.pdf'))

fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=prop.absorptive(resp_xs_xyxy),
                         scale='linear', line=False,
                         fig_kwargs=dict(figsize=(4,3),
                                         constrained_layout=True))
fig_dict['fig'].suptitle('XYXY (no overtones)')
fig_dict['ax2d'].set(xlim=(2030, 2220), ylim=(2030, 2220))
fig_dict['fig'].savefig(str(OUTPUT / 'XYXY_noovertones.png'))
fig_dict['fig'].savefig(str(OUTPUT / 'XYXY_noovertones.pdf'))

fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=prop.absorptive(resp_xs_xmaxma),
                         scale='linear', line=False,
                         fig_kwargs=dict(figsize=(4,3),
                                         constrained_layout=True))
fig_dict['fig'].suptitle('X(MA)X(MA) (no overtones)')
fig_dict['ax2d'].set(xlim=(2030, 2220), ylim=(2030, 2220))
fig_dict['fig'].savefig(str(OUTPUT / 'XMAXMA_noovertones.png'))
fig_dict['fig'].savefig(str(OUTPUT / 'XMAXMA_noovertones.pdf'))

fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=prop.absorptive(resp_xs_bracamonte2)-prop.absorptive(resp_xs_xxxx),
                         scale='linear')
fig_dict['fig'].suptitle('XXXX+X(MA)X(MA)')
fig_dict['ax2d'].set(xlim=(2000, 2250), ylim=(2000, 2250))
