Installation
============

**Dependencies**

- attrs
- asteval
- h5py
- numpy
- scipy
- anytree
- sympy
- matplotlib
- :external:ref:`molspecutils <molspecutils>`

**Install**

.. highlight:: sh

The package is available on PyPI and can be most easily installed with :program:`pip`::

  pip install rotsim2d

This will install the library.
To install GUI applications install the following package::

  pip install rotsim2d_apps

Molecular data
++++++++++++++

`rotsim2d` needs data on rovibrational states and transitions to produce 2D peak plots, simulate 2D lineshapes or time-dependent signals.
For simulations of transitions within an isolated vibrational mode, :class:`rotsim2d.dressedleaf.DressedPathway` needs to be provided with an object implementing `VibrationalMode` interface defined in :mod:`molspecutils.molecule`.
The `molspecutils <https://gitlab.com/allisonlab/mdcs/spectroscopy>`_ package currently contains two classes implementing the interface: `COAlchemyMode` and `CH3ClAlchemyMode`, where the latter is for the :math:`\nu_3` mode.
The molecular data is not bundled with the package and during the first run :mod:`molspecutils` will use HAPI [#f1]_ to fetch transitions line lists from HITRAN [#f2]_ and cache them for future use.

.. rubric:: Footnotes

.. [#f1] |HAPI|
.. [#f2] |HITRAN|
