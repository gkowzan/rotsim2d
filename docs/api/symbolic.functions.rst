.. currentmodule:: rotsim2d.symbolic.functions

###########################
rotsim2d.symbolic.functions
###########################

.. automodule:: rotsim2d.symbolic.functions


Miscellaneous functions
=======================

.. autosummary::
   :toctree: ../generated/

    gfactor_expr
    T1q
    T00_phis
    dl_gfactors
    dl_T00s


RFactor
========

.. autosummary::
   :toctree: ../generated/

   RFactor

Attributes
----------

.. autosummary::
   :toctree: ../generated/

    RFactor.expr
    RFactor.dict
    RFactor.tuple
    RFactor.trigs
    RFactor.angles
    RFactor.angles_type


Methods
-------

.. autosummary::
   :toctree: ../generated/

    RFactor.from_gterms
    RFactor.from_expr
    RFactor.from_pathway
    RFactor.dict_to_expr
    RFactor.expr_xxxx
    RFactor.expr_relative
    RFactor.numeric
    RFactor.numeric_rel
    RFactor.det_angle


Top-level functions
-------------------

.. autosummary::
   :toctree: ../generated/

    classify_dls
    solve_det_angle
    suppression_angles
    common_angles
    classify_suppression
    pathway_angles
    detection_angles


RFactorPathways
===============

.. autosummary::
   :toctree: ../generated/

    RFactorPathways

Attributes
----------

.. autosummary::
   :toctree: ../generated/

    RFactorPathways.rfactor
    RFactorPathways.pws
    RFactorPathways.props
    RFactorPathways.det_angle
    RFactorPathways.trans_labels
    RFactorPathways.trans_labels_deg
    RFactorPathways.peak_labels

Methods
-------

.. autosummary::
   :toctree: ../generated/

    RFactorPathways.from_pwlist
    RFactorPathways.calc_det_angle
    RFactorPathways.add_property

Top-level functions
-------------------

.. autosummary::
   :toctree: ../generated/

    calc_rfactors
    calc_angle_funcs
    calc_angle_amps
    calc_relative_exprs
    calc_amplitude_exprs
    optimize_contrast


Rotational coherence functions
==============================
