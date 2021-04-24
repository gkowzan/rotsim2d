# * Imports and functions
import numpy as np
import scipy.constants as C
import matplotlib.pyplot as plt
import molspecutils.foreign.hapi3 as h3
import pyfftw.interfaces.scipy_fftpack as fftp
from molspecutils import happier
import knickknacks.units as u
from knickknacks.experiment import find_index
import rotsim2d.propagate as prop
import rotsim2d.pathways as pw
import rotsim2d.visual as vis
plt.ion()

# * System properties
pressure, Tgas, length = 5.0, 296.0, 3.0
conc = 1e-7*happier.volumeConcentration(pressure, Tgas)*1e6
co_levels = happier.energy_levels([0,1,2], 5, 1)
co_pops = happier.equilibrium_pops(co_levels, Tgas, 5, 1)
co_params = happier.line_params([(0, 1), (1, 2), (2, 3)], 5, 1)
co_env = {'T': Tgas, 'p': pressure}
co_params = happier.apply_env(co_params, co_env)
co_params = {k: {'mu': v['mu'], 'gam': u.wn2nu(v['gam']), 'nu': u.wn2nu(v['nu'])}
             for k, v in co_params.items()}

# * Response function parameters
dt = 1.0/(2.0*u.wn2nu(343.0))/2.0
T = 1.0/co_params[((0,0), (1,1))]['gam']*3.0/10         # relaxation time constant
# dt = 1e-15
# T = 1e-12
N = np.round(T/dt).astype(np.int)
#t1 = 1/co_params[((0,0), (1,1))]['gam']/20.0
t2 = 0.0
ts = np.arange(N)*dt

# * Calculate molecular coherence evolution
root_prop = prop.Propagator({'elevels': co_levels, 'populations': co_pops,
                             'line_params': co_params},
                            # filter=lambda kb: kb.parent.parent==pw.KetBra(0, 0, 1, 1)\
                            # or kb.parent.parent==pw.KetBra(1, 1, 0, 0),
                            filter=lambda kb: True,
                            freq_shift=u.wn2nu(1800.0))
resp_xs_freqs = root_prop.cross_section_freqs([ts])
resp_xs_freqs_wn = u.nu2wn(resp_xs_freqs)
wnmin, wnmax = 240, 430
imin, imax = find_index(resp_xs_freqs_wn, wnmin), find_index(resp_xs_freqs_wn, wnmax)
Nsmall = imax-imin

# ** First order
resp_xs_spec_o1 = np.zeros(ts.size, dtype=np.complex)
for j in range(5, 6):
    # first order
    root = pw.KetBra(0, j, 0, j)
    root = pw.excite(root, 'omg1', 'both')
    root = pw.readout(root)
    pop_factor = co_pops[(0, j)]/(2*j+1)
    resp_xs_spec_o1 += pop_factor*root_prop.cross_section_spectrum(root, [ts])

# ** Third order
import time
resp_xs_spec = np.zeros((ts.size, Nsmall), dtype=np.complex)
for j in range(5, 6):
    # third order
    print(j, time.asctime())
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['both', 'both', 'both'])
    # root = pw.excite(root, 'omg1', 'both')
    root = pw.readout(root)
    pop_factor = co_pops[(0, j)]/(2*j+1)
    for i, t1 in enumerate(ts):
        resp_xs_spec[i, :] += pop_factor*root_prop.cross_section_spectrum(root, [t1, t2, ts])[imin:imax]

# *** Fourier transform in pump axis
resp_xs_spec_fft = np.zeros((ts.size, Nsmall), dtype=np.complex)
for i in range(Nsmall):
    resp_xs_spec_fft[:, i] = fftp.fft(resp_xs_spec[:, i],
                                      planner_effort='FFTW_ESTIMATE')

# *** Plot
fig_dict = vis.plot2d_im(freqs=(resp_xs_freqs_wn[imin:imax],
                                resp_xs_freqs_wn[imin:imax]),
                         spec2d=np.imag(resp_xs_spec_fft[imin:imax, :][::-1]),
                         spec_linear=(resp_xs_freqs_wn[imin:imax],
                                      np.imag(resp_xs_spec_o1[imin:imax])))

fig_dict = vis.plot2d_im(freqs=(resp_xs_freqs_wn[imin:imax],
                                resp_xs_freqs_wn[imin:imax]),
                         spec2d=np.real(resp_xs_spec_fft + resp_xs_spec_fft[::-1]),
                         spec_linear=(resp_xs_freqs_wn[imin:imax],
                                      np.imag(resp_xs_spec_o1[imin:imax])))

# ** Third order -- no ESA
import time
resp_xs_spec_no_esa = np.zeros((ts.size, Nsmall), dtype=np.complex)
for j in range(5, 6):
    # third order
    print(j, time.asctime())
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['both', 'both', 'both'])
    # root = pw.excite(root, 'omg1', 'both')
    root = pw.remove_esa(root)
    root = pw.readout(root)
    pop_factor = co_pops[(0, j)]/(2*j+1)
    for i, t1 in enumerate(ts):
        resp_xs_spec_no_esa[i, :] += pop_factor*root_prop.cross_section_spectrum(root, [t1, t2, ts])[imin:imax]

# *** Fourier transform in pump axis
resp_xs_spec_no_esa_fft = np.zeros((ts.size, Nsmall), dtype=np.complex)
for i in range(Nsmall):
    resp_xs_spec_no_esa_fft[:, i] = fftp.fft(resp_xs_spec_no_esa[:, i],
                                             planner_effort='FFTW_ESTIMATE')

# *** Plot
fig_dict = vis.plot2d_im(freqs=(resp_xs_freqs_wn[imin:imax],
                                resp_xs_freqs_wn[imin:imax]),
                         spec2d=np.imag(resp_xs_spec_no_esa_fft[imin:imax, :][::-1]),
                         spec_linear=(resp_xs_freqs_wn[imin:imax],
                                      np.imag(resp_xs_spec_o1[imin:imax])))


# ** Third order -- nonrephasing
import time
resp_xs_spec_nonreph = np.zeros((ts.size, Nsmall), dtype=np.complex)
for j in range(5, 6):
    # third order
    print(j, time.asctime())
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['both', 'both', 'both'])
    root = pw.remove_rephasing(root)
    root = pw.readout(root)
    pop_factor = co_pops[(0, j)]/(2*j+1)
    for i, t1 in enumerate(ts):
        resp_xs_spec_nonreph[i, :] += pop_factor*root_prop.cross_section_spectrum(root, [t1, t2, ts])[imin:imax]

# *** Fourier transform in pump axis
resp_xs_spec_nonreph_fft = np.zeros((ts.size, Nsmall), dtype=np.complex)
for i in range(Nsmall):
    resp_xs_spec_nonreph_fft[:, i] = fftp.fft(resp_xs_spec_nonreph[:, i],
                                              planner_effort='FFTW_ESTIMATE')

# *** Plot
fig_dict = vis.plot2d_im(freqs=(resp_xs_freqs_wn[imin:imax],
                                resp_xs_freqs_wn[imin:imax]),
                         spec2d=np.real(resp_xs_spec_nonreph_fft[::-1]),
                         spec_linear=(resp_xs_freqs_wn[imin:imax],
                                      np.imag(resp_xs_spec_o1[imin:imax])))

# ** Third order -- rephasing
import time
resp_xs_spec_reph = np.zeros((ts.size, Nsmall), dtype=np.complex)
for j in range(5, 6):
    # third order
    print(j, time.asctime())
    root = pw.KetBra(0, j, 0, j)
    root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'],
                           parts=['both', 'both', 'both'])
    root = pw.remove_nonrephasing(root)
    root = pw.readout(root)
    pop_factor = co_pops[(0, j)]/(2*j+1)
    for i, t1 in enumerate(ts):
        resp_xs_spec_reph[i, :] += pop_factor*root_prop.cross_section_spectrum(root, [t1, t2, ts])[imin:imax]

# *** Fourier transform in pump axis
resp_xs_spec_reph_fft = np.zeros((ts.size, Nsmall), dtype=np.complex)
for i in range(Nsmall):
    resp_xs_spec_reph_fft[:, i] = fftp.fft(resp_xs_spec_reph[:, i],
                                           planner_effort='FFTW_ESTIMATE')

# *** Plot
fig_dict = vis.plot2d_im(freqs=(resp_xs_freqs_wn[imin:imax],
                                resp_xs_freqs_wn[imin:imax]),
                         spec2d=np.imag(resp_xs_spec_reph_fft[::-1]),
                         spec_linear=(resp_xs_freqs_wn[imin:imax],
                                      np.imag(resp_xs_spec_o1[imin:imax])))

# ** Third order -- rephasing-nonrephasing
# *** Plot
fig_dict = vis.plot2d_im(freqs=(resp_xs_freqs_wn[imin:imax],
                                resp_xs_freqs_wn[imin:imax]),
                         spec2d=np.imag(resp_xs_spec_reph_fft[::-1]+
                                        resp_xs_spec_nonreph_fft),
                         spec_linear=(resp_xs_freqs_wn[imin:imax],
                                      np.imag(resp_xs_spec_o1[imin:imax])))

# * Graph
root = pw.KetBra(0, 5, 0, 5)
root = pw.multi_excite(root, ['omg1', 'omg2', 'omg3'])
root = pw.readout(root)
root = pw.remove_nondiagonal(root)
root.savepng('/mnt/c/Users/heroy/Desktop/third_order_readout.png')

root = pw.remove_rephasing(root)
root.savepng('/mnt/c/Users/heroy/Desktop/third_order_no_esa.png')
