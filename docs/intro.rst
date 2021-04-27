Introduction
============

`rotsim2d` is a library to model two-dimensional infrared (2DIR) rovibrational
spectra. It uses Mukamelian [1]_ density-matrix perturbation theory treatment of
nonlinear spectroscopy. It assumes explicit time-ordering of interactions,
i.e. it is suitable for time-domain spectroscopy with time-separated optical
pulses.

Module :mod:`rotsim2d.pathways` contains classes and functions to produce and
manipulate trees of states that represent all possible n\ :sup:`th` order dipole
excitations starting from some molecular ground state. A specific pathway (a
double-sided Feynmann diagram) can be represented by :class:`Pathway` from
:mod:`rotsim2d.dressedleaf`, which also allows to calculate the
polarization-angular momentum coupling factor for the pathway. An abstract
pattern of excitations can be attached to some specific molecular vibrational
mode by a :class:`DressedPathway` instance.  At this level, one can produce a 2D
map of resonance peaks (see :class:`Peak2D`, :class:`Peak2DList`,
:func:`peak_list`) without calculating line shapes or explicit time dependence.

Module :mod:`rotsim2d.propagate` provides classes and functions for generating
actual time-dependent or frequency-domain experimental
signals. :mod:`rotsim2d.visual` provides convenience functions to visualize 2D
peak maps, 2D spectra and double-sided Feynmann diagrams with LaTeX/TikZ.

Modules :mod:`rotsim2d.symbolic.functions` and :mod:`rotsim2d.symbolic.results`
provide SymPy expressions and functions used to derive polarization and angular
momentum dependence of four-fold dipole interaction operator.

Actual data on vibrational modes of real molecules is provided by
:mod:`molspecutils` package, which (so far) wraps data from the HITRAN database
[2]_.

.. [1] S. Mukamel. "Principles of nonlinear spectroscopy". Oxford University Press, New York, 1995.
.. [2] I.E. Gordon *et al.* "The HITRAN2016 molecular spectroscopic database". J. Quant. Spectrosc. Radiat. Transf. *203*, 3-69 (2017).
