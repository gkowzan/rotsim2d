# * Imports and functions
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyfftw.interfaces.scipy_fftpack as fftp
import rotsim2d.pathways as pw
import rotsim2d.propagate as prop
import rotsim2d.visual as vis
import shed.units as u
import pywigxjpf as wig
from spectroscopy import happier

prop.ignore_missing = False
plt.ion()

OUTPUT = Path('/mnt/d/DFCS/Stony Brook/rotsim2d_results/carbon_monoxide/polarization/waiting_time')
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

df2 = 50e9
F2 = u.wn2nu(5000.0)
N2 = np.round(F2/df2).astype(np.int)
fs2 = np.arange(-N2, N2)*df2
fs2_cm = u.nu2wn(fs2)

# * Calculate molecular coherence spectrum
root_prop = prop.Spectrum({'elevels': co_levels, 'populations': co_pops, 'line_params': co_params},
                            # filter=lambda kb: kb.parent.parent==pw.KetBra(0, 0, 1, 1)\
                            # or kb.parent.parent==pw.KetBra(1, 1, 0, 0),
                          filter=lambda kb: kb.total_side() == -1,
                          # filter=lambda kb: True,
                          freq_shift=[u.wn2nu(1800.0), 0.0, 0.0])

# ** Third order --- full -- waiting frequency
pws = []
for j in range(0, 16):
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'])
    root = pw.readout(root)
    root = pw.remove_nondiagonal(root)
    root.pop = co_pops[(0, j)]/np.sqrt(2*j+1)
    pws.append(root)
multi_prop = prop.MultiPropagator(pws, root_prop)
with wig.wigxjpf(300, 6):
    resp_xs_full = multi_prop.response([fs[:, np.newaxis], fs2[np.newaxis, :], None])

# *** Plot
fig_dict = vis.plot2d_im(freqs=(fs_cm+1800.0, fs2_cm),
                         spec2d=prop.absorptive(resp_xs_full),
                         scale='linear')
fig_dict['fig'].canvas.set_window_title('Frequency domain -- absorptive -- all pathways -- waiting time')
fig_dict['ax2d'].set_xlabel(r'Waiting axis (cm$^{-1}$)')

# * Calculate molecular coherence time domain signal
# ** Time domain response function parameters
dt = 10e-15
T = 5e-11
N = np.round(T/dt).astype(np.int)
ts = np.arange(N)*dt

root_prop = prop.Propagator({'elevels': co_levels, 'populations': co_pops, 'line_params': co_params},
                            # filter=lambda kb: kb.parent.parent==pw.KetBra(0, 0, 1, 1)\
                            # or kb.parent.parent==pw.KetBra(1, 1, 0, 0),
                            filter=lambda kb: kb.total_side() == -1,
                            freq_shift=[0.0, 0.0, 0.0])

# ** Third order --- full -- waiting time
pws = []
for j in range(0, 25):
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'])
    root = pw.remove_overtones(root)
    root = pw.only_interstates(root)
    root = pw.readout(root)
    root = pw.remove_nondiagonal(root)
    root.pop = co_pops[(0, j)]
    # resp_xs_time = root_prop.response(root, [None, ts, None])
    pws.append(root)
multi_prop = prop.MultiPropagator(pws, root_prop)
with wig.wigxjpf(300, 6):
    resp_xs_time = multi_prop.response([None, ts, None])

# ** Check single leaf responses
fig_re, ax_re = plt.subplots(figsize=(4.0, 3.0))
fig_im, ax_im = plt.subplots(figsize=(4.0, 3.0))
for j in range(0, 25):
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'])
    root = pw.remove_overtones(root)
    root = pw.only_interstates(root)
    root = pw.readout(root)
    root = pw.remove_nondiagonal(root)
    root.pop = co_pops[(0, j)]
    for l in root.leaves:
        with wig.wigxjpf(300, 6):
            leaf_resp = root_prop.leaf_response(l, [None, ts, None])
        ax_re.plot(ts*1e12, np.real(leaf_resp))
        ax_im.plot(ts*1e12, np.imag(leaf_resp))
ax_re.set(xlabel=r'Waiting time (ps)', ylabel='Molecular coherence',
          xlim=(4.1, 4.7))
ax_re.set_title('Real part', fontsize=9.0)
ax_im.set(xlabel=r'Waiting time (ps)', ylabel='Molecular coherence',
          xlim=(1.9, 2.6))
ax_im.set_title('Imaginary part', fontsize=9.0)
# fig_re.savefig(OUTPUT/ 'waiting_time_individual_real.png', dpi=300.0)
# fig_re.savefig(OUTPUT/ 'waiting_time_individual_real.pdf', dpi=300.0)
# fig_im.savefig(OUTPUT/ 'waiting_time_individual_imag.png', dpi=300.0)
# fig_im.savefig(OUTPUT/ 'waiting_time_individual_imag.pdf', dpi=300.0)

# *** plot
fig, ax = plt.subplots(figsize=(4.0, 3.0))
ax.plot(ts*1e12, np.real(resp_xs_time), label='real')
ax.plot(ts*1e12, np.imag(resp_xs_time), label='imag')
ax.legend(loc='best')
ax.set(xlabel=r'Waiting time (ps)', ylabel='Molecular coherence',
       xlim=(0, 20)
       )
# fig.canvas.set_window_title('Time-domain -- waiting time dependence')
# fig.savefig(OUTPUT / 'waiting_time_revivals.png', dpi=300.0)
# fig.savefig(OUTPUT / 'waiting_time_revivals.pdf', dpi=300.0)

# * Calculate molecular coherence time domain signal -- different waiting times
# ** Time domain response function parameters
dt = 1.0/(2.0*u.wn2nu(343.0))/2.0
T = 1.0/co_params[((0,0), (1,1))]['gam']*4.0/5.0         # relaxation time constant
N = np.round(T/dt).astype(np.int)
# t2 = 4375e-15
t2 = 2190e-15
# t2= 0.0
ts = np.arange(N)*dt

root_prop = prop.Propagator({'elevels': co_levels, 'populations': co_pops, 'line_params': co_params},
                            # filter=lambda kb: kb.parent.parent==pw.KetBra(0, 0, 1, 1)\
                            # or kb.parent.parent==pw.KetBra(1, 1, 0, 0),
                            filter=lambda kb: True,
                            freq_shift=[u.wn2nu(1800.0), 0.0, u.wn2nu(1800.0)])

# ** Third order --- full -- pump and detection time
resp_xs = np.zeros((ts.size, ts.size), dtype=np.complex)
for j in range(0, 15):
    # third order
    print(j, time.asctime())
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['ket', 'both', 'both'],
                           # light_angles=[54.7356/180.0*np.pi, 0.0, 54.7356/180.0*np.pi]
                           light_angles=[np.pi/2, np.pi/4, np.pi/2])
    root = pw.remove_rephasing(root)
    # root = pw.remove_overtones(root)
    root = pw.readout(root, 26.57/180.0*np.pi)
    pop_factor = co_pops[(0, j)]/np.sqrt(2*j+1)
    for i, t1 in enumerate(ts):
        with wig.wigxjpf(300, 6):
            resp_xs[i, :] += pop_factor*root_prop.response(root, [t1, t2, ts])

# *** Fourier transform in pump axis
resp_xs_freqs_wn = fftp.fftshift(u.nu2wn(root_prop.cross_section_freqs([ts])))
resp_xs_spec_fft = fftp.fftshift(fftp.fft2(resp_xs, planner_effort='FFTW_ESTIMATE'))
if resp_xs_freqs_wn.size % 2 == 0:
    # Make positive and negative frequencies symmetric around DC.
    resp_xs_freqs_wn = resp_xs_freqs_wn[1:]
    resp_xs_spec_fft = resp_xs_spec_fft[1:, 1:]
# Invert the order elements to agree with Mukamel's FT convention.
resp_xs_spec_fft = resp_xs_spec_fft[::-1, ::-1]

# *** Plot
fig_dict = vis.plot2d_im(freqs=(resp_xs_freqs_wn+1800.0, resp_xs_freqs_wn+1800.0),
                         spec2d=prop.absorptive(resp_xs_spec_fft),
                         scale='linear')
fig_dict['fig'].set_size_inches(4.0, 3.0)
# fig_dict['fig'].canvas.set_window_title('Time domain -- absorptive -- no overtones')
fig_dict['ax2d'].set_title(r'$t_2 = {:.3f}$ ps'.format(t2*1e12), fontsize=9.0)
fig_dict['ax2d'].set(xlim=(2030, 2230), ylim=(2030, 2230))
fig_dict['fig'].savefig(OUTPUT / '2d_{:s}.png'.format("{:.3f}".format(t2*1e12).replace('.', '_')), dpi=300.0)
fig_dict['fig'].savefig(OUTPUT / '2d_{:s}.pdf'.format("{:.3f}".format(t2*1e12).replace('.', '_')), dpi=300.0)
