.. barotropy documentation master file, created by
   sphinx-quickstart on Thu Sep 28 17:55:22 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to barotropy's documentation!
=====================================

.. module:: barotropy
   :synopsis: A Python package to generate barotropic fluid property models.

``barotropy`` is a Python package to generate barotropic fluid property models. 
A barotropic model is a simplified mathematical representation that assumes 
the fluid's properties are only a function of pressure. This simplification 
can be useful to model flows in which the fluid approximately follows a known 
thermodynamic process, such as an isentropic process within a turbine or compressor.

Using a barotropic model in the context of a CFD simulation does not imply 
that the flow is frictionless since the momentum equations can have dissipation 
terms due to viscous friction. Instead, the implication of using the barotropic 
model is that the variation of fluid properties caused by viscous dissipation 
(i.e., entropy generation and heating) is ignored.

The thermodynamic properties are computed using the CoolProp library. Once computed, 
the fluid properties are automatically fitted with piece-wise polynomials and exported 
as `Fluent Expressions <https://ansyshelp.ansys.com/account/secured?returnurl=/Views/Secured/corp/v232/en/flu_ug/flu_ug_expressions_intro.html>`_ 
to be used in conjunction with Fluent's pressure-based solver.

The fluid properties between the saturation line and the spinodal line can be 
computed according to phase-equilibrium or by extrapolating the Helmholtz-energy 
equation of state into the metastable region. Check the 
`documentation <./documentation/thermodynamic_properties.rst>`_ for more information 
about extrapolating the equation of state beyond the saturation line and its limitations.

Use the panel to the left or the table of contents below to navigate the documentation.

.. toctree::
   :maxdepth: 2
   :caption: Contents:


   source/installation
   source/tutorials
   source/barotropic_model
   source/bibliography
   source/api/modules

