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

OUTPUT = Path('/mnt/d/DFCS/Stony Brook/rotsim2d_results/carbon_monoxide/summer_internship_ad')

# * System properties
def co_props(co_env, I):
    co_levels = happier.energy_levels([0,1,2], 5, I)
    co_pops = happier.equilibrium_pops(co_levels, co_env['T'], 5, I)
    co_params = happier.line_params([(0, 1), (1, 2), (2, 3), (0, 2)], 5, I)
    co_params_suppl = happier.generate_line_params(
        [(0, 0, 0), (1, 1, 0), (2, 2, 0), (0, 0, 2), (1, 1, 2), (0, 2, 2), (0, 2, 0)], co_levels, co_params)
    co_params.update(co_params_suppl)
    co_params = happier.apply_env(co_params, co_env, to_nu=True)

    return {'elevels': co_levels, 'populations': co_pops, 'line_params': co_params}


iso1_props = co_props({'p': 5.0, 'T': 296.0}, 1)
iso2_props = co_props({'p': 5.0, 'T': 296.0}, 2)
iso3_props = co_props({'p': 5.0, 'T': 296.0}, 3)

# * Response function parameters
df = 10e9
F = u.wn2nu(700.0)
N = np.round(F/df).astype(np.int)
fs = np.arange(-N, N+1)*df
fs_cm = u.nu2wn(fs)
jmax = 20

def response_calc(props, pols, resp_xs, scale=1.0):
    root_prop = prop.Spectrum(props, filter=lambda kb: True,
                              freq_shift=[u.wn2nu(1800.0), 0.0, u.wn2nu(1800.0)])

    pws = []
    for j in range(0, jmax):
        root = pw.KetBra(0, j, 0, j)
        root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                               parts=['ket', 'both', 'both'],
                               light_angles=pols[:3])
        # root = pw.remove_rephasing(root)
        root = pw.remove_interstates(root)
        root = pw.readout(root, pols[3])
        root.pop = props['populations'][(0, j)]/np.sqrt(2*j+1)
        pws.append(root)
    multi_prop = prop.MultiPropagator(pws, root_prop)
    with wig.wigxjpf(300, 6):
        resp_xs += multi_prop.response([fs[:, np.newaxis], None, fs[np.newaxis, :]])*scale

    return resp_xs


# * new
sopr = -36.87*np.pi/180
pprr = 26.57*np.pi/180
resp_xs_new = np.zeros((fs.size, fs.size), dtype=np.complex)
# ** iso1
resp_xs_new = response_calc(iso1_props, [np.pi/2, sopr, np.pi/2, np.pi/4], resp_xs_new)
# ** Plot
fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=prop.absorptive(resp_xs_new),
                         scale='symlog', line=True,
                         fig_kwargs=dict(figsize=(3.5, 2.7), constrained_layout=True))
# fig_dict['fig'].suptitle('LLXX (no overtones)')
fig_dict['ax2d'].set(xlim=(1990, 2220), ylim=(1990, 2220))
fig_dict['axcbar'].set_ticks([])
fig_dict['axcbar'].set_label(r'$\Delta$OD')
# fig_dict['fig'].savefig(str(OUTPUT / 'new_freq_abs.pdf'))
# fig_dict['fig'].savefig(str(OUTPUT / 'new_freq_abs.png'), dpi=300.0)


# * XXXX
resp_xs_xxxx = np.zeros((fs.size, fs.size), dtype=np.complex)

# ** iso1
resp_xs_xxxx = response_calc(iso1_props, [0.0, 0.0, 0.0, 0.0], resp_xs_xxxx)

# ** iso2
resp_xs_xxxx = response_calc(iso2_props, [0.0, 0.0, 0.0, 0.0], resp_xs_xxxx)

# ** iso3
resp_xs_xxxx = response_calc(iso3_props, [0.0, 0.0, 0.0, 0.0], resp_xs_xxxx)

# ** Plot
fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=prop.absorptive(resp_xs_xxxx),
                         scale='symlog', line=True,
                         fig_kwargs=dict(figsize=(3.5, 2.7), constrained_layout=True))
# fig_dict['fig'].suptitle('LLXX (no overtones)')
fig_dict['ax2d'].set(xlim=(1990, 2220), ylim=(1990, 2220))
fig_dict['axcbar'].set_ticks([])
fig_dict['axcbar'].set_label(r'$\Delta$OD')
fig_dict['fig'].savefig(str(OUTPUT / 'XXXX_freq_abs.pdf'))
fig_dict['fig'].savefig(str(OUTPUT / 'XXXX_freq_abs.png'), dpi=300.0)

# * XXXX 2
resp_xs_xxxx2 = np.zeros((fs.size, fs.size), dtype=np.complex)

# ** iso2
resp_xs_xxxx2 = response_calc(iso2_props, [0.0, 0.0, 0.0, 0.0], resp_xs_xxxx2)

# ** iso3
resp_xs_xxxx2 = response_calc(iso3_props, [0.0, 0.0, 0.0, 0.0], resp_xs_xxxx2)

# ** Plot
fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=prop.absorptive(resp_xs_xxxx2),
                         scale='symlog', line=True, pthresh=10.0,
                         fig_kwargs=dict(figsize=(3.5, 2.7), constrained_layout=True))
# fig_dict['fig'].suptitle('LLXX (no overtones)')
fig_dict['ax2d'].set(xlim=(2050, 2150), ylim=(2050, 2150))
fig_dict['axcbar'].set_ticks([])
fig_dict['axcbar'].set_label(r'$\Delta$OD')
fig_dict['fig'].savefig(str(OUTPUT / 'XXXX_freq_abs2.pdf'))
fig_dict['fig'].savefig(str(OUTPUT / 'XXXX_freq_abs2.png'), dpi=300.0)

# * X(\pi+MA)X(MA)
ma = 54.7356103172453*np.pi/180.0 # magic angle
# ** iso1
root_prop = prop.Spectrum(iso1_props, filter=lambda kb: True,
                          freq_shift=[u.wn2nu(1800.0), 0.0, u.wn2nu(1800.0)])

pws = []
for j in range(0, jmax):
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'],
                           light_angles=[0.0, np.pi+ma, 0.0])
    root = pw.remove_rephasing(root)
    root = pw.readout(root, ma)
    root.pop = iso1_props['populations'][(0, j)]/np.sqrt(2*j+1)
    pws.append(root)
multi_prop = prop.MultiPropagator(pws, root_prop)
with wig.wigxjpf(300, 6):
    resp_xs_twozero = multi_prop.response([fs[:, np.newaxis], None, fs[np.newaxis, :]])

# ** iso2
root_prop = prop.Spectrum(iso2_props, filter=lambda kb: True,
                          freq_shift=[u.wn2nu(1800.0), 0.0, u.wn2nu(1800.0)])

pws = []
for j in range(0, jmax):
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'],
                           light_angles=[0.0, np.pi+ma, 0.0])
    root = pw.remove_rephasing(root)
    root = pw.readout(root, ma)
    root.pop = iso2_props['populations'][(0, j)]/np.sqrt(2*j+1)
    pws.append(root)
multi_prop = prop.MultiPropagator(pws, root_prop)
with wig.wigxjpf(300, 6):
    resp_xs_twozero += multi_prop.response([fs[:, np.newaxis], None, fs[np.newaxis, :]])

# ** iso3
root_prop = prop.Spectrum(iso3_props, filter=lambda kb: True,
                          freq_shift=[u.wn2nu(1800.0), 0.0, u.wn2nu(1800.0)])

pws = []
for j in range(0, jmax):
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'],
                           light_angles=[0.0, np.pi+ma, 0.0])
    root = pw.remove_rephasing(root)
    root = pw.readout(root, ma)
    root.pop = iso3_props['populations'][(0, j)]/np.sqrt(2*j+1)
    pws.append(root)
multi_prop = prop.MultiPropagator(pws, root_prop)
with wig.wigxjpf(300, 6):
    resp_xs_twozero += multi_prop.response([fs[:, np.newaxis], None, fs[np.newaxis, :]])

# ** Plot
fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=prop.absorptive(resp_xs_twozero),
                         scale='linear', line=False,
                         fig_kwargs=dict(figsize=(4,3), constrained_layout=True))
# fig_dict['fig'].suptitle('LLXX (no overtones)')
fig_dict['ax2d'].set(xlim=(1990, 2220), ylim=(1990, 2220))

# * X(\pi+MA)(MA)X
ma = 54.7356103172453*np.pi/180.0 # magic angle
# ** iso1
root_prop = prop.Spectrum(iso1_props, filter=lambda kb: True,
                          freq_shift=[u.wn2nu(1800.0), 0.0, u.wn2nu(1800.0)])

pws = []
for j in range(0, jmax):
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'],
                           light_angles=[0.0, np.pi+ma, ma])
    root = pw.remove_rephasing(root)
    root = pw.readout(root, 0.0)
    root.pop = iso1_props['populations'][(0, j)]/np.sqrt(2*j+1)
    pws.append(root)
multi_prop = prop.MultiPropagator(pws, root_prop)
with wig.wigxjpf(300, 6):
    resp_xs_twozero2 = multi_prop.response([fs[:, np.newaxis], None, fs[np.newaxis, :]])

# ** iso2
root_prop = prop.Spectrum(iso2_props, filter=lambda kb: True,
                          freq_shift=[u.wn2nu(1800.0), 0.0, u.wn2nu(1800.0)])

pws = []
for j in range(0, jmax):
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'],
                           light_angles=[0.0, np.pi+ma, ma])
    root = pw.remove_rephasing(root)
    root = pw.readout(root, 0.0)
    root.pop = iso2_props['populations'][(0, j)]/np.sqrt(2*j+1)
    pws.append(root)
multi_prop = prop.MultiPropagator(pws, root_prop)
with wig.wigxjpf(300, 6):
    resp_xs_twozero2 += multi_prop.response([fs[:, np.newaxis], None, fs[np.newaxis, :]])

# ** iso3
root_prop = prop.Spectrum(iso3_props, filter=lambda kb: True,
                          freq_shift=[u.wn2nu(1800.0), 0.0, u.wn2nu(1800.0)])

pws = []
for j in range(0, jmax):
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'],
                           light_angles=[0.0, np.pi+ma, ma])
    root = pw.remove_rephasing(root)
    root = pw.readout(root, 0.0)
    root.pop = iso3_props['populations'][(0, j)]/np.sqrt(2*j+1)
    pws.append(root)
multi_prop = prop.MultiPropagator(pws, root_prop)
with wig.wigxjpf(300, 6):
    resp_xs_twozero2 += multi_prop.response([fs[:, np.newaxis], None, fs[np.newaxis, :]])

# ** Plot
fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=prop.absorptive(resp_xs_twozero2),
                         scale='linear', line=False,
                         fig_kwargs=dict(figsize=(4,3), constrained_layout=True))
# fig_dict['fig'].suptitle('LLXX (no overtones)')
fig_dict['ax2d'].set(xlim=(1990, 2220), ylim=(1990, 2220))
