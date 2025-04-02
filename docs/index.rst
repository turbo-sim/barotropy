
Welcome to barotropy's documentation!
=====================================

.. module:: barotropy
   :synopsis: A Python package to generate barotropic fluid property models.

Enhance your two-phase turbomachinery CFD simulations with **barotropy**,  
a Python package designed to generate **barotropic fluid property models** for use in CFD simlations.


A **barotropic model** assumes that thermophysical properties such as **density**, **viscosity**, and **speed of sound** depend **only on pressure**.  
This simplification is highly accurate for flows in which the fluid approximately follows a known 
thermodynamic process, such as a polytropic processes in turbines and compressors. This simplification enables:

- Efficient and stable CFD convergence
- Physically plausible extrapolation near critical or two-phase conditions
- Seamless integration with pressure-based CFD solvers such as Fluent or CFX

Using a barotropic model in the context of a CFD simulation does not imply 
that the flow is frictionless. This is because momentum equations still have dissipation 
terms due to viscous friction. Instead, the implication of using the barotropic 
model is that the variation of fluid properties caused by viscous dissipation described by the polytropic process.

The thermodynamic properties are computed using the CoolProp library. Once computed, 
the fluid properties are automatically fitted with piece-wise polynomials and exported 
as Expressions to be used in conjunction with Fluent's or CFX's pressure-based solvers.

Two-phase fluid properties between the saturation line and the spinodal line can be 
computed according to phase-equilibrium or by extrapolating the Helmholtz-energy 
equation of state into the metastable region. Check the 
`documentation <./documentation/thermodynamic_properties.rst>`_ for more information 
about extrapolating the equation of state beyond the saturation line and its limitations.


.. raw:: html

    <p align="center">
      <img src="_static/two-phase_turbine.jpg"  width="85%" />
    </p>

    <p align="center">
      <img src="_static/sCO2_compressor.jpg"  width="85%" />
    </p>





Table of contents
------------------
Use the panel to the left or the table of contents below to navigate the documentation.



.. toctree::
   :maxdepth: 2

   source/introduction
   source/tutorials
   source/barotropic_model
   source/developer_guide
   source/bibliography
   source/api/barotropy

