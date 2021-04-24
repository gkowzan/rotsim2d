"""Compare HAPI and rotsim2d calculations of linear absorption."""
# * Imports and functions
import numpy as np
import scipy.constants as C
import matplotlib.pyplot as plt
import molspecutils.foreign.hapi3 as h3
from molspecutils import happier
from knickknacks.experiment import find_index
import knickknacks.units as u
import rotsim2d.propagate as prop
import rotsim2d.pathways as pw


def mps(M, T):
    """Calculate the most probable speed.

    Parameters
    ----------
    M : float
        Mass in g/mol.
    T : float
        Temperature in Kelvin.

    Returns
    -------
    float
        The most probable speed in m/s.
    """
    return np.sqrt(2*C.R*T/M*1e3)


# * Physical conditions
Tgas, Tref, pressure = 296.0, 296.0, 5.0
wnmin, wnmax = 2000, 2280
conc = 1e-7*happier.volumeConcentration(pressure, Tgas)*1e6
length = 1.0            # mm

# * Calculate spectrum with HITRAN
# ** Load data
h3.db_begin('data')
if 'carbon_monoxide' not in h3.getTableList():
    h3.fetch('carbon_monoxide', 5, 1, 2000, 5000)

# ** Select the fundamental band of carbon monoxide
co_conds = ('AND', ('=', 'molec_id', 5), ('=', 'local_iso_id', 1),
            ('MATCH', ('STR', '\s+0$'), 'global_lower_quanta'),
            ('MATCH', ('STR', '\s+1$'), 'global_upper_quanta'))
h3.select('carbon_monoxide', Conditions=co_conds, DestinationTableName='CO_band01')

# ** Line-shape function
def lorentz(nu, nu0, gam, trans_xs):
    return trans_xs/np.pi/(nu-nu0-1.0j*gam)

N = np.round((wnmax-wnmin)/u.nu2wn(50e6)).astype(np.int)
wns = np.linspace(wnmin, wnmax, N)
spec = np.zeros(wns.size, dtype=np.complex)

for trans_xs, nu0, gam, shift in zip(*h3.getColumns('CO_band01', ['sw', 'nu', 'gamma_air', 'delta_air'])):
    spec += lorentz(wns, nu0+shift*pressure, gam*pressure, trans_xs)

# * Calculate spectrum with rotsim2d
# ** System properties
co_levels = happier.energy_levels([0,1,2], 5, 1)
co_pops = happier.equilibrium_pops(co_levels, 296.0, 5, 1)
co_params = happier.line_params([(0, 1)], 5, 1)
co_env = {'T': Tgas, 'p': pressure}
co_params = happier.apply_env(co_params, co_env)
co_params = {k: {'mu': v['mu'], 'gam': u.wn2nu(v['gam']), 'nu': u.wn2nu(v['nu'])}
             for k, v in co_params.items()}

# ** Response function parameters
dt = 1/u.wn2nu(2200.0)/4
T = 1/co_params[((0,0), (1,1))]['gam']*9         # relaxation time constant
# dt = 1e-15
# T = 1e-12
N = np.round(T/dt).astype(np.int)
ts = np.arange(N)*dt

# ** Calculate cross section spectrum
root_prop = prop.Propagator({'elevels': co_levels, 'populations': co_pops,
                             'line_params': co_params},
                            # filter=lambda kb: kb.parent.parent==pw.KetBra(0, 0, 1, 1)\
                            # or kb.parent.parent==pw.KetBra(1, 1, 0, 0),
                            filter=lambda kb: True)
resp_xs_spec = np.zeros(ts.size, dtype=np.complex)
for j in range(50):
    root = pw.KetBra(0, j, 0, j)
    root = pw.excite(root, 'omg1', 'both')
    root = pw.readout(root)
    pop_factor = co_pops[(0, j)]
    resp_xs_spec += pop_factor*root_prop.cross_section_spectrum(root, [ts])
resp_xs_freqs = root_prop.cross_section_freqs([ts])

# * Plot
# mpl.rcParams['figure.dpi'] = 192.0
plt.ion()
fig, ax = plt.subplots()
ax.plot(wns, np.imag(spec), label='HITRAN')
ax.plot(u.nu2wn(resp_xs_freqs), -2*np.imag(resp_xs_spec)*1e4, '--', label='rotsim2d')
ax.set(
    # xlim=(wnmin, wnmax),
    ylim=(0.0, 1e-18),
    xlabel=r'Wavenumber (cm$^{-1}$)',
    ylabel=r'Absorption cross section (cm$^{2}$)')
ax.legend(loc='best')

# * Calculate spectrum with rotsim2d -- frequency shift
# ** System properties
co_levels = happier.energy_levels([0,1,2], 5, 1)
co_pops = happier.equilibrium_pops(co_levels, 296.0, 5, 1)
co_params = happier.line_params([(0, 1)], 5, 1)
co_env = {'T': Tgas, 'p': pressure}
co_params = happier.apply_env(co_params, co_env)
co_params = {k: {'mu': v['mu'], 'gam': u.wn2nu(v['gam']), 'nu': u.wn2nu(v['nu'])}
             for k, v in co_params.items()}

# ** Response function parameters
dt = 1/(2*u.wn2nu(343.0))/2
T = 1/co_params[((0,0), (1,1))]['gam']*5         # relaxation time constant
# dt = 1e-15
# T = 1e-12
N = np.round(T/dt).astype(np.int)
ts = np.arange(N)*dt

# ** Calculate cross section spectrum
root_prop = prop.Propagator({'elevels': co_levels, 'populations': co_pops,
                             'line_params': co_params},
                            # filter=lambda kb: kb.parent.parent==pw.KetBra(0, 0, 1, 1)\
                            # or kb.parent.parent==pw.KetBra(1, 1, 0, 0),
                            filter=lambda kb: True,
                            freq_shift=u.wn2nu(1800.0))
resp_xs_spec = np.zeros(ts.size, dtype=np.complex)
for j in range(50):
    root = pw.KetBra(0, j, 0, j)
    root = pw.excite(root, 'omg1', 'both')
    root = pw.readout(root)
    pop_factor = co_pops[(0, j)]
    resp_xs_spec += pop_factor*root_prop.cross_section_spectrum(root, [ts])
resp_xs_freqs = root_prop.cross_section_freqs([ts])

# * Plot
# mpl.rcParams['figure.dpi'] = 192.0
plt.ion()
fig, ax = plt.subplots()
ax.plot(wns-1800.0, np.imag(spec), label='HITRAN')
ax.plot(u.nu2wn(resp_xs_freqs), -2*np.imag(resp_xs_spec)*1e4, '--', label='rotsim2d')
ax.set(
    # xlim=(wnmin, wnmax),
    ylim=(0.0, 1e-18),
    xlabel=r'Wavenumber (cm$^{-1}$)',
    ylabel=r'Absorption cross section (cm$^{2}$)')
ax.legend(loc='best')
