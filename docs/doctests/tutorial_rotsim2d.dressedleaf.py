import rotsim2d.pathways as pw
import rotsim2d.dressedleaf as dl
import molspecutils.molecule as mol
kb1 = pw.KetBra(mol.DiatomState(nu=0, j=1), mol.DiatomState(nu=0, j=1))
kb1 = pw.multi_excite(kb1, light_names=['omg1', 'omg2', 'omg3'],
                      parts=['ket', 'both', 'both'],
                      light_angles=[0]*3)
kb1 = pw.readout(kb1)
pw1 = dl.Pathway(kb1.leaves[0])
print()
pw1.pprint()

pws1 = dl.Pathway.from_kb_tree(kb1)
print(pws1[:2])

js = (0, 1, 2, 3)
kbs = pw.gen_pathways(js, pols=[0]*4, meths=[pw.only_SII, pw.only_twocolor],
                      rotor='symmetric', kiter_func=lambda j: range(j+1))
pws2 = dl.Pathway.from_kb_list(kbs)

pws1_by_peaks = dl.split_by_peaks(pws1)
dl.print_dl_dict(pws1_by_peaks, fields=['trans_label', 'tw_coherence', 'js'])

from molspecutils.molecule import CH3ClAlchemyMode
import rotsim2d.visual as vis

ch3cl_mode = CH3ClAlchemyMode()
kbs = pw.gen_pathways(range(1, 10), [0]*4, meths=[pw.only_SII, pw.only_twocolor],
                      rotor='symmetric', kiter_func=lambda x: range(x if x<10 else 10))
dressed_pws = dl.DressedPathway.from_kb_list(kbs, ch3cl_mode, 296.0)
peaks, dls = dl.peak_list(dressed_pws, return_dls=True)
fig_dict = vis.plot2d_scatter(peaks, fig_kwargs=dict(figsize=(4.5, 3.2)))
fig_dict['ax'].set(xlabel='Probe (cm$^{-1}$)', ylabel='Pump (cm$^{-1}$)')
fig_dict['fig'].tight_layout()
fig_dict['fig'].savefig('../images/ch3cl_2d_plot.png')
