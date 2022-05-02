.. currentmodule:: rotsim2d.dressedleaf

###########################################
Manipulating individual excitation pathways
###########################################

Individual abstract pathways
============================

:class:`Pathway` represents an individual excitation pathway (a double-sided
Feynman diagram) not associated with any specific molecular vibrational
mode. It can be instantiated with a leaf node of
:class:`rotsim2d.pathways.KetBra` excitation tree:

.. ipython:: python

   import rotsim2d.pathways as pw
   import rotsim2d.dressedleaf as dl
   import molspecutils.molecule as mol
   from pprint import pprint
   kb1 = pw.KetBra(mol.DiatomState(nu=0, j=1), mol.DiatomState(nu=0, j=1))
   kb1 = pw.multi_excite(kb1, light_names=['omg1', 'omg2', 'omg3'],
                         parts=['ket', 'both', 'both'])
   kb1 = pw.readout(kb1)
   pw1 = dl.Pathway(kb1.leaves[0])
   pw1.pprint()

Transition chain label classifies the transitions induced by each light
interaction by a sequence of three letters: P, Q or R. This is the same notation
as the one used in [1]_.

A list of all pathways in a given excitation tree can be obtained by using the
class method :meth:`Pathway.from_kb_tree`:

.. ipython:: python

   pws1 = dl.Pathway.from_kb_tree(kb1)
   pprint(pws1[:2])

With :meth:`Pathway.from_kb_list` we can obtain a list pathways from a list of
excitation trees:

.. ipython:: python

   js = (0, 1, 2, 3)
   kbs = pw.gen_pathways(js, meths=[pw.only_SII, pw.only_twocolor],
                         rotor='symmetric', kiter_func='range(j+1)')

Such potentially long list of pathways can be classified by 2D peaks they
generate and printed out for easier inspection:

.. ipython:: python

   pws1_by_peaks = dl.split_by_peaks(pws1)
   dl.print_dl_dict(pws1_by_peaks, fields=['trans_label', 'tw_coherence', 'js'])

`tw_coherence` tells whether the system is in a coherence state or population
state between second and third pulse. `fields` argument determines what
information is printed for each pathway. Valid field names are defined in
:attr:`Pathway.fields`.

Pathways associated with vibrational modes
==========================================

To associate a pathway with some concrete vibrational mode, we need an instance
of a class implementing the :class:`molspecutils.molecule.VibrationalMode`
interface.  Two such classes are :class:`molspecutils.molecule.COAlchemyMode`
and :class:`molspecutils.molecule.CH3ClAlchemyMode`. Let's use the latter:

.. ipython:: python

   from molspecutils.molecule import CH3ClAlchemyMode
   import rotsim2d.visual as vis

   ch3cl_mode = CH3ClAlchemyMode()
   kbs = pw.gen_pathways(range(1, 10), meths=[pw.only_SII, pw.only_twocolor],
                         rotor='symmetric', kiter_func='range(j if j<10 else 10)')
   dressed_pws = dl.DressedPathway.from_kb_list(kbs, ch3cl_mode, 296.0)                    
   peaks = dl.Peak2DList.from_dp_list(dressed_pws, tw=1.0e-12, angles=[0.0]*4)
   fig_dict = vis.plot2d_scatter(peaks)                                                    
   fig_dict['fig'].savefig('images/ch3cl_2d_plot.png')

.. image:: /images/ch3cl_2d_plot.png
   :width: 80%
   :align: center

As before, we can filter the excitation tree with functions from
:mod:`rotsim2d.pathways`. We can also investigate the polarization dependence
and waiting time dependence by changing `angles` and `tw` arguments of
:func:`Peak2DList.from_dp_list`. :func:`~rotsim2d.visual.plot2d_scatter` is a
convenience functions for plotting a list of 2D peaks.
	   
.. [1] D. Murdock, L. A. Burns, P. H. Vaccaro. J. Phys. Chem. A 113, 13184-13198 (2009). doi:10.1021/jp903970d.
