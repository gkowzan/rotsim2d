.. rotsim2d documentation master file, created by
   sphinx-quickstart on Fri Sep 18 13:46:05 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Simulate rotationally-resolved 2D infrared spectra
==================================================

`rotsim2d` is a library to model two-dimensional infrared (2DIR) rovibrational
spectra. It uses Mukamelian density-matrix perturbation theory treatment of
nonlinear spectroscopy. It assumes explicit time-ordering of interactions,
i.e. it is suitable for time-domain spectroscopy with time-separated optical
pulses.
	    
.. toctree::
   :maxdepth: 2

   install
   intro
   api

Getting started
---------------

**Manipulating excitation tree**

Create :class:`KetBra` in ground vibrational state

.. testcode::
   
   import rotsim2d.pathways as pw
   import molspecutils.molecule as mol
   kb1 = pw.KetBra(mol.DiatomState(nu=0, j=1), mol.DiatomState(nu=0, j=1))
   print(kb1)

.. testoutput::

   |0,1><0,1|

and excite it,

.. testcode::

   pw.excite(kb1, light_name='omg1')
   kb1.print_tree()
   
.. testoutput::

   |0,1><0,1|
   └── omg1
       ├── |1,0><0,1|
       └── |1,2><0,1|
   
`kb1` is modified in-place. Use `readout` to simulate measurement and see which
branches would survive taking the trace of the density matrix.

.. testcode::

   pw.readout(kb1)
   kb1.print_tree()

.. testoutput::

   |0,1><0,1|
   └── omg1
       ├── |1,0><0,1|
       │   └── mu
       │       └── |0,1><0,1|
       └── |1,2><0,1|
           └── mu
               └── |0,1><0,1|
