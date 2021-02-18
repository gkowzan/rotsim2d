# * Imports
import time
from pathlib import Path
import numpy as np
import pyfftw.interfaces.scipy_fftpack as fftp
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


iso1_props = co_props({'p': 15.0, 'T': 296.0}, 1)
iso2_props = co_props({'p': 15.0, 'T': 296.0}, 2)
iso3_props = co_props({'p': 15.0, 'T': 296.0}, 3)

# * Isotopologue-dependent rotational coherence
jmax = 20

def rcs_calc(props):
    dt = 10e-15
    T = 5e-11
    N = np.round(T/dt).astype(np.int)
    ts = np.arange(N)*dt

    root_prop = prop.Propagator(props,
                                filter=lambda kb: kb.total_side() == -1,
                                freq_shift=[0.0, 0.0, 0.0])

    pws = []
    for j in range(0, jmax):
        root = pw.KetBra(0, j, 0, j)
        root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                               parts=['ket', 'both', 'both'])
        root = pw.only_SII(root)
        root = pw.only_interstates(root)
        root = pw.readout(root)
        root.pop = props['populations'][(0, j)]/np.sqrt(2*j+1)
        pws.append(root)
    multi_prop = prop.MultiPropagator(pws, root_prop)
    with wig.wigxjpf(300, 6):
        resp_xs_time = multi_prop.response([None, ts, None])

    return ts, resp_xs_time


# ** Calculate
ts, iso1_rcs = rcs_calc(iso1_props)
# iso2_rcs = rcs_calc(iso2_props)[1]
# iso3_rcs = rcs_calc(iso3_props)[1]

# ** Plot
fig, ax = plt.subplots(figsize=(4.0, 4.0))
ax.plot(ts*1e12, np.imag(iso1_rcs), label='iso1')
# ax.plot(ts*1e12, np.imag(iso2_rcs), label='iso2')
# ax.plot(ts*1e12, np.imag(iso3_rcs), label='iso3')
ax.legend(loc='best')
ax.set(xlabel=r'Waiting time (ps)', ylabel='Molecular coherence',
       xlim=(0, 20)
       )

# ** Single-leaf response
def rcs_single_plot(props):
    dt = 10e-15
    T = 5e-11
    N = np.round(T/dt).astype(np.int)
    ts = np.arange(N)*dt

    root_prop = prop.Propagator(props,
                                filter=lambda kb: kb.total_side() == 1,
                                freq_shift=[0.0, 0.0, 0.0])

    fig_re, ax_re = plt.subplots(figsize=(4.0, 3.0))
    fig_im, ax_im = plt.subplots(figsize=(4.0, 3.0))
    for j in range(0, 25):
        root = pw.KetBra(0, j, 0, j)
        root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                               parts=['ket', 'both', 'both'])
        root = pw.only_SI(root)
        root = pw.only_interstates(root)
        root = pw.readout(root)
        root.pop = props['populations'][(0, j)]/np.sqrt(2*j+1)
        for l in root.leaves:
            with wig.wigxjpf(300, 6):
                leaf_resp = root_prop.leaf_response(l, [None, ts, None])
            ax_re.plot(ts*1e12, np.real(leaf_resp))
            ax_im.plot(ts*1e12, np.imag(leaf_resp))
    ax_re.set(xlabel=r'Waiting time (ps)', ylabel='Molecular coherence',
              xlim=(0, 4.7))
    ax_re.set_title('Real part', fontsize=9.0)
    ax_im.set(xlabel=r'Waiting time (ps)', ylabel='Molecular coherence',
              xlim=(0, 4.7))
    ax_im.set_title('Imaginary part', fontsize=9.0)


rcs_single_plot(iso2_props)


# * Response function parameters
dt = 1.0/(2.0*u.wn2nu(343.0))/2.0
T = 1.0/iso1_props['line_params'][((0,0), (1,1))]['gam']*5.0         # relaxation time constant
N = np.round(T/dt).astype(np.int)
# t2 = 4375e-15
t2 = 2190e-15
# for other isotopologues
# t2 = 2290e-15
ts = np.arange(N)*dt
jmax = 20

def response_calc(props, pols, resp_xs, t2=0.0):
    root_prop = prop.Propagator(props, filter=lambda kb: True,
                                freq_shift=[u.wn2nu(1800.0), 0.0, u.wn2nu(1800.0)])

    for j in range(0, jmax):
        print(j, time.asctime())
        root = pw.KetBra(0, j, 0, j)
        root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                               parts=['ket', 'both', 'both'],
                               light_angles=pols[:3])
        root = pw.remove_rephasing(root)
        root = pw.readout(root, pols[3])
        pop_factor = props['populations'][(0, j)]/np.sqrt(2*j+1)
        for i, t1 in enumerate(ts):
            with wig.wigxjpf(300, 6):
                resp_xs[i, :] += pop_factor*root_prop.response(root, [t1, t2, ts])

    return resp_xs


def iso_response_calc(pols, t2=0.0):
    resp_xs = np.zeros((ts.size, ts.size), dtype=np.complex)
    response_calc(iso1_props, pols, resp_xs, t2)
    response_calc(iso2_props, pols, resp_xs, t2)
    response_calc(iso3_props, pols, resp_xs, t2)

    # Fourier transform
    root_prop = prop.Propagator(iso1_props, filter=lambda kb: True,
                                freq_shift=[u.wn2nu(1800.0), 0.0, u.wn2nu(1800.0)])
    resp_xs_freqs_wn = fftp.fftshift(u.nu2wn(root_prop.cross_section_freqs([ts])))
    resp_xs_spec_fft = fftp.fftshift(fftp.fft2(resp_xs, planner_effort='FFTW_ESTIMATE'))
    if resp_xs_freqs_wn.size % 2 == 0:
        # Make positive and negative frequencies symmetric around DC.
        resp_xs_freqs_wn = resp_xs_freqs_wn[1:]
        resp_xs_spec_fft = resp_xs_spec_fft[1:, 1:]
    # Invert the order elements to agree with Mukamel's FT convention.
    resp_xs_spec_fft = resp_xs_spec_fft[::-1, ::-1]

    return resp_xs_freqs_wn, resp_xs_spec_fft

# * XXXX
xxxx_wns, xxxx_spec = iso_response_calc([0.0, 0.0, 0.0, 0.0], t2=0.0)

# ** Plot
fig_dict = vis.plot2d_im(freqs=(xxxx_wns+1800.0, xxxx_wns+1800.0),
                         spec2d=prop.absorptive(xxxx_spec),
                         scale='linear', line=False,
                         fig_kwargs=dict(figsize=(3.5, 2.7), constrained_layout=True))
# fig_dict['fig'].suptitle('LLXX (no overtones)')
fig_dict['ax2d'].set(xlim=(1990, 2220), ylim=(1990, 2220))
fig_dict['axcbar'].set_ticks([])
fig_dict['axcbar'].set_label(r'$\Delta$OD')
fig_dict['fig'].savefig(str(OUTPUT / 'XXXX_t.pdf'))
fig_dict['fig'].savefig(str(OUTPUT / 'XXXX_t.png'), dpi=300.0)

# * XXXX (t2)
xxxxt2_wns, xxxxt2_spec = iso_response_calc([0.0, 0.0, 0.0, 0.0], t2=2190e-15)

# ** Plot
fig_dict = vis.plot2d_im(freqs=(xxxxt2_wns+1800.0, xxxxt2_wns+1800.0),
                         spec2d=prop.absorptive(xxxxt2_spec),
                         scale='linear', line=False,
                         fig_kwargs=dict(figsize=(3.5, 2.7), constrained_layout=True))
# fig_dict['fig'].suptitle('LLXX (no overtones)')
fig_dict['ax2d'].set(xlim=(1990, 2220), ylim=(1990, 2220))
fig_dict['axcbar'].set_ticks([])
fig_dict['axcbar'].set_label(r'$\Delta$OD')
fig_dict['fig'].savefig(str(OUTPUT / 'XXXXt2_t.pdf'))
fig_dict['fig'].savefig(str(OUTPUT / 'XXXXt2_t.png'), dpi=300.0)

# * X(\pi+MA)X(MA)
ma = 54.7356103172453*np.pi/180.0 # magic angle
xmaxma_wns, xmaxma_spec = iso_response_calc([0.0, np.pi+ma, 0.0, ma], t2=2190e-15)

# ** Plot
fig_dict = vis.plot2d_im(freqs=(xmaxma_wns+1800.0, xmaxma_wns+1800.0),
                         spec2d=prop.absorptive(xmaxma_spec),
                         scale='linear', line=False,
                         fig_kwargs=dict(figsize=(3.5, 2.7), constrained_layout=True))
# fig_dict['fig'].suptitle('LLXX (no overtones)')
fig_dict['ax2d'].set(xlim=(1990, 2220), ylim=(1990, 2220))
fig_dict['axcbar'].set_ticks([])
fig_dict['axcbar'].set_label(r'$\Delta$OD')
fig_dict['fig'].savefig(str(OUTPUT / 'XMAXMA_t.pdf'))
fig_dict['fig'].savefig(str(OUTPUT / 'XMAXMA_t.png'), dpi=300.0)

# * X(\pi+MA)(MA)X
ma = 54.7356103172453*np.pi/180.0 # magic angle
xmaxma2_wns, xmaxma2_spec = iso_response_calc([0.0, np.pi+ma, ma, 0.0], t2=2190e-15)

# ** Plot
fig_dict = vis.plot2d_im(freqs=(xmaxma2_wns+1800.0, xmaxma2_wns+1800.0),
                         spec2d=prop.absorptive(xmaxma2_spec),
                         scale='linear', line=False,
                         fig_kwargs=dict(figsize=(3.5, 2.7), constrained_layout=True))
# fig_dict['fig'].suptitle('LLXX (no overtones)')
fig_dict['ax2d'].set(xlim=(1990, 2220), ylim=(1990, 2220))
fig_dict['ax2d'].set(xlim=(1990, 2220), ylim=(1990, 2220))
fig_dict['axcbar'].set_ticks([])
fig_dict['axcbar'].set_label(r'$\Delta$OD')
fig_dict['fig'].savefig(str(OUTPUT / 'XMAMAX_t.pdf'))
fig_dict['fig'].savefig(str(OUTPUT / 'XMAMAX_t.png'), dpi=300.0)
