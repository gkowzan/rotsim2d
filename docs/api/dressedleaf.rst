.. currentmodule:: rotsim2d.dressedleaf

####################		   
rotsim2d.dressedleaf
####################

.. automodule:: rotsim2d.dressedleaf

Pathway
=======

.. autosummary::
   :toctree: ../generated/

   Pathway

Attributes
----------

.. autosummary::
   :toctree: ../generated/

   Pathway.leaf
   Pathway.colors
   Pathway.experimental_label
   Pathway.geo_label
   Pathway.peak_label
   Pathway.trans_label
   Pathway.trans_label_deg
   Pathway.coherences
   Pathway.transitions
   Pathway.js
   Pathway.light_inds
   Pathway.const
   Pathway.tw_coherence
   Pathway.peak
   Pathway.abstract_peak

Methods
-------

.. autosummary::
   :toctree: ../generated/

   Pathway.from_kb_tree
   Pathway.from_kb_list
   Pathway.T00s
   Pathway.gfactors
   Pathway.geometric_factor
   Pathway.custom_str
   Pathway.print_diagram
   Pathway.pprint


DressedPathway
==============

.. autosummary::
   :toctree: ../generated/

   DressedPathway


Attributes
----------

.. autosummary::
   :toctree: ../generated/

   DressedPathway.vib_mode
   DressedPathway.T
   DressedPathway.const


Methods
-------

.. autosummary::
   :toctree: ../generated/

   DressedPathway.from_kb_tree
   DressedPathway.from_kb_list
   DressedPathway.from_params_dict
   DressedPathway.gamma
   DressedPathway.nu
   DressedPathway.amplitude
   DressedPathway.pump_fraction
   DressedPathway.pprint


Peak2D
======

.. autosummary::
   :toctree: ../generated/

   Peak2D

Attributes
----------

.. autosummary::
   :toctree: ../generated/

   Peak2D.pump_wl
   Peak2D.probe_wl
   Peak2D.peak
   Peak2D.amplitude
   Peak2D.intensity
   Peak2D.max_intensity
   Peak2D.dp_list

Methods
-------

.. autosummary::
   :toctree: ../generated/

   Peak2D.max_abs_coeff
   Peak2D.from_dp_list


Peak2DList
==========

.. autosummary::
   :toctree: ../generated/

   Peak2DList

Attributes
----------

.. autosummary::
   :toctree: ../generated/

   Peak2DList.pumps
   Peak2DList.probes
   Peak2DList.peaks
   Peak2DList.amplitudes
   Peak2DList.intensities
   Peak2DList.max_intensities

Methods
-------

.. autosummary::
   :toctree: ../generated/

   Peak2DList.from_dp_list
   Peak2DList.from_file
   Peak2DList.from_params_dict
   Peak2DList.to_file
   Peak2DList.get_by_peak
   Peak2DList.sort_by_amplitudes


Top-level functions
===================

Labels
------

.. autosummary::
   :toctree: ../generated/

   abstract_format
   abstract_state_label
   abstract_pair_label
   abstract_line_label

Pathways
--------

.. autosummary::
   :toctree: ../generated/

   split_by_js
   split_by_peaks
   pprint_dllist
   print_dl_dict
   print_dl_tuple_dict

Peaks
-----

.. autosummary::
   :toctree: ../generated/

   equiv_peaks

Miscellaneous
-------------

.. autosummary::
   :toctree: ../generated/

   perm_pols
   perm_js
   undegenerate_js
