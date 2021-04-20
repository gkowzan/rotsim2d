# * Imports and functions
from typing import Tuple
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.cm as cm
from spectroscopy import happier
import shed.units as u
import rotsim2d.propagate as prop
prop.ignore_missing = False
import rotsim2d.pathways as pw
import rotsim2d.dressedleaf as dl
import rotsim2d.visual as vis
import pywigxjpf as wig
plt.ion()
PAPER_OUTPUT = Path('/mnt/d/DFCS/ARTYKULY/2021 rotsim PRL/figures')


def co_props(co_env, I):
    co_levels = happier.energy_levels([0,1,2], 5, I)
    co_pops = happier.equilibrium_pops(co_levels, co_env['T'], 5, I)
    co_params = happier.line_params([(0, 1), (1, 2), (2, 3), (0, 2), (0, 0), (1, 1), (2, 2)], 5, 1)
    co_params_suppl = happier.generate_line_params(
        [(0, 0, 0), (1, 1, 0), (2, 2, 0), (0, 0, 2), (1, 1, 2), (0, 2, 2), (0, 2, 0)], co_levels, co_params)
    co_params_suppl2 = happier.generate_line_params([(0, 1, 0), (1, 2, 0)], co_levels, co_params)
    co_params.update(co_params_suppl)
    co_params.update(co_params_suppl2)
    co_params = happier.apply_env(co_params, co_env, to_nu=True)

    return {'elevels': co_levels, 'populations': co_pops, 'line_params': co_params}


def response_calc(pws, props, resp_xs, scale=1.0):
    root_prop = prop.Spectrum(props, filter=lambda kb: True,
                              freq_shift=[u.wn2nu(1800.0), 0.0, u.wn2nu(1800.0)])
    multi_prop = prop.MultiPropagator(pws, root_prop)
    with wig.wigxjpf(300, 6):
        resp_xs += multi_prop.response([fs[:, np.newaxis], None, fs[np.newaxis, :]])*scale

    return resp_xs


def load_npz(fpath, name):
    temp = np.load(fpath)
    arr = temp[name]
    temp.close()
    
    return arr

iso1_props = co_props({'p': 15.0,  'T': 296.0}, 1)

# * Response function parameters
df = 10e9
F = u.wn2nu(700.0)
N = np.round(F/df).astype(np.int)
fs = np.arange(-N, N+1)*df
fs_cm = u.nu2wn(fs)
jmax = 30

# * Base case
# ** XXXX (SI)
# pws = pw.gen_pathways(range(jmax), [0, 0, 0, 0], iso1_props['populations'], [pw.only_SI])
# resp_xs_xxxx_SI = np.zeros((fs.size, fs.size), dtype=np.complex)
# resp_xs_xxxx_SI = response_calc(pws, iso1_props, resp_xs_xxxx_SI)
# np.savez(PAPER_OUTPUT / 'resp_xs_xxxx_SI.npz', resp_xs_xxxx_SI=resp_xs_xxxx_SI)
resp_xs_xxxx_SI = load_npz(PAPER_OUTPUT / 'resp_xs_xxxx_SI.npz', 'resp_xs_xxxx_SI')

# ** XXXX (SII)
# pws = pw.gen_pathways(range(jmax), [0, 0, 0, 0], iso1_props['populations'], [pw.only_SII])
# resp_xs_xxxx = np.zeros((fs.size, fs.size), dtype=np.complex)
# resp_xs_xxxx = response_calc(pws, iso1_props, resp_xs_xxxx)
# np.savez(PAPER_OUTPUT / 'resp_xs_xxxx.npz', resp_xs_xxxx=resp_xs_xxxx)
resp_xs_xxxx = load_npz(PAPER_OUTPUT / 'resp_xs_xxxx.npz', 'resp_xs_xxxx')

# ** Combine
fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=prop.absorptive(resp_xs_xxxx+resp_xs_xxxx_SI),
                         scale='linear', line=True,
                         fig_kwargs=dict(figsize=(3.5, 2.7), constrained_layout=True))
# fig_dict['fig'].suptitle('LLXX (no overtones)')
fig_dict['ax2d'].set(xlim=(2030, 2230), ylim=(2030, 2230))
fig_dict['axcbar'].set_ticks([])
fig_dict['axcbar'].set_label(r'$\Delta$OD')
# fig_dict['fig'].savefig(NIST_OUTPUT / 'full_absorptive_spectrum.png', dpi=300.0)
# None

# * Absorptive diagonal
# ** Zero diagonal, set waiting time (SI) - t2 movie
# *** Prepare time frequency axes
def aligned_fs(fsmin: float, fsmax: float, df: float):
    def align(f: float):
        return np.ceil(f/df).astype(np.int) if f < 0 else np.floor(f/df).astype(np.int)
    return np.arange(align(fsmin), align(fsmax)+1)*df
    
df = 10e9
fs1, fs2 = aligned_fs(0, u.wn2nu(110), df), aligned_fs(0, u.wn2nu(150), df)
fs1_cm = u.nu2wn(fs1)
fs2_cm = u.nu2wn(fs2)
B = 57898343910.56784
TRCS = 1/4/B
ts2 = np.linspace(0, TRCS, 200)

# *** Calculate
pws = pw.gen_pathways(range(10), [0, np.pi/4, np.pi/2, -np.pi/4], iso1_props['populations'], [pw.only_SII])
# removes antidiagonal constant peaks
resp_xs_inter_SI = np.zeros((fs1.size, ts2.size, fs2.size), dtype=np.complex)
dl_list = dl.dress_pws(pws, iso1_props)
for i, t2 in enumerate(ts2):
    print(i)
    with wig.wigxjpf(300, 6):
        for dli in dl_list:
            resp_xs_inter_SI[:, i, :] += prop.dressed_leaf_response(
                dli, [fs1[:, None], t2, fs2[None, :]],
                ['f', 't', 'f'], [u.wn2nu(2090.0), 0.0, u.wn2nu(2050.0)])
# np.savez(PAPER_OUTPUT / 'resp_xs_inter_SI.npz', resp_xs_inter_SI=resp_xs_inter_SI)
# resp_xs_inter_SI = load_npz(PAPER_OUTPUT / 'resp_xs_inter_SI.npz', 'resp_xs_inter_SI')

# *** Set up animation
def plot2d_animation(freqs, spec3d, absmax=None, fig_kwargs={}):
    extent = vis.make_extent(freqs[0], freqs[2])
    ts2 = freqs[1]
    cmap = cm.get_cmap('RdBu').reversed()
    fig = plt.figure(**fig_kwargs)
    gs = fig.add_gridspec(nrows=1, ncols=2, width_ratios=[20, 1], figure=fig)
    ax2d = fig.add_subplot(gs[0])
    axcbar = fig.add_subplot(gs[1])
    if absmax is None:
        absmax = np.max(np.abs(spec2d))
    cset = ax2d.imshow(np.zeros((spec3d.shape[0], spec3d.shape[2])),
                       cmap=cmap, aspect='auto', extent=extent, clim=(-absmax, absmax),
                       origin='lower')
    ax2d.set(xlabel=r'Probe (cm$^{-1}$)', ylabel=r'Pump (cm$^{-1}$)')
    title = ax2d.text(0.5, 0.99, "$t_2$", ha='center', va='top', transform=ax2d.transAxes)
    axcbar = fig.colorbar(cset, ax=ax2d, cax=axcbar)
    fig.set_constrained_layout_pads(wspace=0.01, hspace=0.01, h_pad=0.01, w_pad=0.01)

    def init():
        return [cset, title]

    def update(i):
        cset.set_data(spec3d[:, i, :])
        title.set_text('$t_2={:.2f}$ ps'.format(ts2[i]*1e12))
        return [cset, title]

    return fig, init, update, np.arange(spec2d.shape[1])

spec2d = np.real(resp_xs_inter_SI)-np.imag(resp_xs_inter_SI)
fig, init, update, frames = plot2d_animation((fs1_cm, ts2, fs2_cm), spec2d)
ani = FuncAnimation(fig, update, frames, init, blit=True, repeat=False, interval=250)
# ani.save('/mnt/d/DFCS/Stony Brook/rotsim2d_results/carbon_monoxide/waiting_time/diagonal.mp4')

# *** Compare with base case
absmax = np.max(np.abs(prop.absorptive(resp_xs_inter_SI)))
fig_dict = vis.plot2d_im(freqs=(fs1_cm, fs2_cm),
                         spec2d=np.imag(resp_xs_inter_SI),
                         scale='linear', line=False, absmax=absmax,
                         fig_kwargs=dict(figsize=(3.5, 2.7), constrained_layout=True))
#fig_dict['ax2d'].set(xlim=(-2030, -2230), ylim=(2030, 2230))

# ** Zero diagonal, set waiting time (SI)
pws = pw.gen_pathways(range(jmax), [0, np.pi/4, np.pi/2, -np.pi/4], iso1_props['populations'], [pw.only_SI])
resp_xs_inter_SI = np.zeros((fs.size, fs.size), dtype=np.complex)
dl_list = dl.dress_pws(pws, iso1_props)
with wig.wigxjpf(300, 6):
    for dli in dl_list:
        resp_xs_inter_SI[...] += prop.dressed_leaf_response(
            dli, [fs[:, None], 0.0, fs[None, :]],
            ['f', 't', 'f'], [u.wn2nu(1800.0), 0.0, u.wn2nu(1800.0)])
# np.savez(PAPER_OUTPUT / 'resp_xs_inter_SI.npz', resp_xs_inter_SI=resp_xs_inter_SI)
# resp_xs_inter_SI = load_npz(PAPER_OUTPUT / 'resp_xs_inter_SI.npz', 'resp_xs_inter_SI')

# *** Compare with base case
absmax = np.max(np.abs(prop.absorptive(resp_xs_inter_SI)))
fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=prop.absorptive(resp_xs_inter_SI),
                         scale='linear', line=False, absmax=absmax,
                         fig_kwargs=dict(figsize=(3.5, 2.7), constrained_layout=True))
fig_dict['ax2d'].set(xlim=(2030, 2230), ylim=(2030, 2230))

fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=prop.absorptive(resp_xs_xxxx_SI),
                         scale='linear', line=False, absmax=absmax,
                         fig_kwargs=dict(figsize=(3.5, 2.7), constrained_layout=True))
fig_dict['ax2d'].set(xlim=(2030, 2230), ylim=(2030, 2230))

# ** Zero antidiagonal and interstates (SII)
pws = pw.gen_pathways(range(jmax), [0, np.pi/4, 0, -np.arctan(4/3)], iso1_props['populations'], [pw.only_SII])
resp_xs_anti_inter = np.zeros((fs.size, fs.size), dtype=np.complex)
resp_xs_anti_inter = response_calc(pws, iso1_props, resp_xs_anti_inter)
# np.savez(PAPER_OUTPUT / 'resp_xs_anti_inter.npz', resp_xs_anti_inter=resp_xs_anti_inter)
resp_xs_anti_inter = load_npz(PAPER_OUTPUT / 'resp_xs_anti_inter.npz', 'resp_xs_anti_inter')

# ** Combine
fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=prop.absorptive(resp_xs_anti_inter-3*resp_xs_inter_SI),
                         scale='linear', line=True,
                         fig_kwargs=dict(figsize=(3.5, 2.7), constrained_layout=True))
# fig_dict['fig'].suptitle('LLXX (no overtones)')
fig_dict['ax2d'].set(xlim=(2030, 2230), ylim=(2030, 2230))
fig_dict['axcbar'].set_ticks([])
fig_dict['axcbar'].set_label(r'$\Delta$OD')
# fig_dict['fig'].savefig(PAPER_OUTPUT / 'absorptive_diagonal.pdf')
# fig_dict['fig'].savefig(PAPER_OUTPUT / 'absorptive_diagonal.png', dpi=300.0)

# * Absorptive antidiagonal
# ** Zero diagonal and interstates (SI)
# pws = pw.gen_pathways(range(20), [0, np.pi/4, np.pi/2, np.arctan(2/3)], iso1_props['populations'], [pw.only_SI])
# resp_xs_anti_inter_SI = np.zeros((fs.size, fs.size), dtype=np.complex)
# resp_xs_anti_inter_SI = response_calc(pws, iso1_props, resp_xs_anti_inter_SI)
# np.savez(PAPER_OUTPUT / 'resp_xs_anti_inter_SI.npz', resp_xs_anti_inter_SI=resp_xs_anti_inter_SI)
resp_xs_anti_inter_SI = load_npz(PAPER_OUTPUT / 'resp_xs_anti_inter_SI.npz', 'resp_xs_anti_inter_SI')

# ** Zero diagonal and set waiting time (SII)
# pws = pw.gen_pathways(range(20), [0, np.pi/4, np.pi/2, -np.pi/4], iso1_props['populations'], [pw.only_SII])
# resp_xs_diag = np.zeros((fs.size, fs.size), dtype=np.complex)
# # resp_xs_diag = response_calc(pws, iso1_props, resp_xs_diag)
# dl_list = dl.dress_pws(pws, iso1_props)
# with wig.wigxjpf(300, 6):
#     for dli in dl_list:
#         resp_xs_diag[...] += prop.dressed_leaf_response(dli, [fs[:, None], 2.12e-12, fs[None, :]],
#                                                         ['f', 't', 'f'],
#                                                         [u.wn2nu(1800.0), 0.0, u.wn2nu(1800.0)])
# np.savez(PAPER_OUTPUT / 'resp_xs_diag.npz', resp_xs_diag=resp_xs_diag)
resp_xs_diag = load_npz(PAPER_OUTPUT / 'resp_xs_diag.npz', 'resp_xs_diag')

# ** Combine
fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs_cm+1800.0),
                         spec2d=prop.absorptive(resp_xs_anti_inter_SI-1.5*resp_xs_diag),
                         scale='linear', line=True,
                         fig_kwargs=dict(figsize=(3.5, 2.7), constrained_layout=True))
# fig_dict['fig'].suptitle('LLXX (no overtones)')
fig_dict['ax2d'].set(xlim=(2030, 2230), ylim=(2030, 2230))
fig_dict['axcbar'].set_ticks([])
fig_dict['axcbar'].set_label(r'$\Delta$OD')
# fig_dict['fig'].savefig(PAPER_OUTPUT / 'absorptive_antidiagonal.pdf')
# fig_dict['fig'].savefig(PAPER_OUTPUT / 'absorptive_antidiagonal.png', dpi=300.0)
# None

ax1_data = prop.absorptive(resp_xs_anti_inter-3*resp_xs_inter_SI)
ax2_data = prop.absorptive(resp_xs_anti_inter_SI-1.5*resp_xs_diag)
ax0_data = prop.absorptive(resp_xs_xxxx+resp_xs_xxxx_SI)

# * Paper figure
extent = vis.make_extent(fs_cm+1800.0, fs_cm+1800.0)

FIGW = 2*(3+3/8)
from shed.matplotlib import jqsrt_rc, compact
import matplotlib as mpl
import matplotlib.ticker as tck
mpl.rcParams.update(jqsrt_rc)
mpl.rcParams.update(compact)
mpl.rcParams['text.usetex'] = True

fig = plt.figure(figsize=(FIGW, 2.2))
gs = fig.add_gridspec(nrows=1, ncols=4, width_ratios=[20, 20, 20, 1])

absmax = np.max(np.abs(ax0_data))
ax0 = fig.add_subplot(gs[0])
ax0_cset = ax0.imshow(ax0_data, cmap='seismic', aspect='auto', extent=extent,
                      clim=(-absmax, absmax), origin='lower')
ax0.set(xlabel=r'$\tilde{\nu}_3$ (cm$^{-1}$)', ylabel=r'$\tilde{\nu}_1$ (cm$^{-1}$)',
        xlim=(2030, 2230), ylim=(2030, 2230))
ax0.text(-0.02, 1.0, r'(a)', transform=ax0.transAxes, ha='right', va='top')
ax0.yaxis.set_major_locator(tck.MultipleLocator(50.0))

ax1 = fig.add_subplot(gs[1])
ax1_cset = ax1.imshow(ax1_data, cmap='seismic', aspect='auto', extent=extent,
                      clim=(-absmax, absmax), origin='lower')
ax1.set(xlabel=r'$\tilde{\nu}_3$ (cm$^{-1}$)', yticklabels=[],
        xlim=(2030, 2230), ylim=(2030, 2230))
ax1.text(-0.02, 1.0, r'(b)', transform=ax1.transAxes, ha='right', va='top')
ax1.yaxis.set_major_locator(tck.MultipleLocator(50.0))

ax2 = fig.add_subplot(gs[2])
ax2_cset = ax2.imshow(ax2_data, cmap='seismic', aspect='auto', extent=extent,
                      clim=(-absmax, absmax), origin='lower')
ax2.set(xlabel=r'$\tilde{\nu}_3$ (cm$^{-1}$)', yticklabels=[],
        xlim=(2030, 2230), ylim=(2030, 2230))
ax2.text(-0.02, 1.0, r'(c)', transform=ax2.transAxes, ha='right', va='top')
ax2.yaxis.set_major_locator(tck.MultipleLocator(50.0))

axc = fig.add_subplot(gs[3])
axcbar = fig.colorbar(ax0_cset, cax=axc)
axcbar.set_ticks([])
axcbar.set_label(r'$\Delta$OD')

fig.set_constrained_layout_pads(wspace=0.03, hspace=0.01, h_pad=0.01, w_pad=0.01)
fig.savefig(PAPER_OUTPUT / 'suppress_line_shapes.pdf')
fig.savefig(PAPER_OUTPUT / 'suppress_line_shapes.png')
