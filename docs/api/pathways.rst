.. currentmodule:: rotsim2d.pathways

#################		   
rotsim2d.pathways
#################

.. automodule:: rotsim2d.pathways

Enums
=====

.. autosummary::
   :toctree: ../generated/

   Side
   KSign

LightInteraction
================
.. autosummary::
   :toctree: ../generated/

   LightInteraction

Attributes
----------

.. autosummary::
   :toctree: ../generated/

   LightInteraction.name
   LightInteraction.side
   LightInteraction.sign
   LightInteraction.angle
   LightInteraction.readout
   LightInteraction.tier

Methods
-------

.. autosummary::
   :toctree: ../generated/

   LightInteraction.is_absorption
   LightInteraction.evolve
   LightInteraction.conj

KetBra
======

.. autosummary::
   :toctree: ../generated/

   KetBra

Attributes
----------

.. autosummary::
   :toctree: ../generated/

   KetBra.ket
   KetBra.bra
   KetBra.name

   
Simple predicates
-----------------
These methods are primarily useful for leaf nodes.

.. autosummary::
   :toctree: ../generated/

   KetBra.is_diagonal
   KetBra.is_rephasing
   KetBra.is_SI
   KetBra.is_SII
   KetBra.is_SIII
   KetBra.is_Pinitial
   KetBra.is_Rinitial
   KetBra.is_Qinitial
   KetBra.is_Qbranch
   KetBra.is_esa
   KetBra.is_sep
   KetBra.is_gshb
   KetBra.is_doublequantum
   KetBra.is_interstate
   KetBra.is_dfwm
   KetBra.is_twocolor
   KetBra.is_threecolor

Comparison
----------

.. autosummary::
   :toctree: ../generated/

   KetBra.is_between
   KetBra.is_same_branch
   KetBra.is_equiv_pathway
   KetBra.is_pathway
   KetBra.is_some_pathway
   
Miscellaneous
-------------

.. autosummary::
   :toctree: ../generated/

   KetBra.copy
   KetBra.evolve
   KetBra.get
   KetBra.conj
   KetBra.savepng
   KetBra.normalized
   KetBra.print_tree
   KetBra.kb_ancestor
   KetBra.diagonals
   KetBra.ketbras
   KetBra.interactions
   KetBra.interaction
   KetBra.ksigns
   KetBra.total_ksign
   KetBra.sides
   KetBra.total_side
   KetBra.transitions
   KetBra.color_tier
   KetBra.to_statelist
   

Top-level functions
===================

.. _tree-filtering-functions:

Tree-filtering functions
------------------------

.. autosummary::
   :toctree: ../generated/

   only_between
   only_pathway
   only_some_pathway
   make_remove
   make_only

The tree-filtering functions below were generated with :func:`make_only` and
:func:`make_remove`. Their names should be self-explanatory. For example::

  only_SI = make_only(lambda kb: kb.is_SI())

.. autosummary::
   :toctree: ../generated/

   only_SI
   only_SII
   only_SIII
   only_nonrephasing
   only_rephasing
   remove_nondiagonal
   remove_overtones
   only_esa
   remove_esa
   only_sep
   only_gshb
   only_dfwm
   only_twocolor
   remove_threecolor
   remove_interstates
   only_interstates
   only_ketside
   remove_ketside
   only_Pinitial
   only_Rinitial
   only_Qinitial
   only_Qbranch
   remove_Qbranch


Tree manipulation and generation
--------------------------------

.. autosummary::
   :toctree: ../generated/

   prune
   readout
   flip_readout
   copy_chain
   conjugate_chain
   zannify_chain
   zannify_tree
   excited_states_symtop
   excited_states_diatom
   excite
   multi_excite
   gen_roots
   gen_excitations
   gen_pathways
   geometric_factor
