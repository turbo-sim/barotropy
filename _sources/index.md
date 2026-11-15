# Welcome to barotropy's documentation!

Conduct your two-phase turbomachinery CFD simulations with **barotropy**, a Python package to generate **barotropic fluid property models** for CFD simulations.


A **barotropic model** assumes that thermophysical properties such as **density**, **viscosity**, and **speed of sound** depend **only on pressure**.  This simplification is highly accurate for flows in which the fluid approximately follows a known thermodynamic process, such as a polytropic processes in turbines and compressors. This simplification enables:

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
equation of state into the metastable region. See the [barotropic model](barotropic_model.md) for the modeling overview.



```{image} ../_static/two-phase_turbine.jpg
:width: 85%
:align: center
:alt: Two-phase turbine
```

```{image} ../_static/sCO2_compressor.jpg
:width: 95%
:align: center
:alt: Supercritical CO2 compressor
```

## Table of contents

```{toctree}
:maxdepth: 2

introduction
tutorials
barotropic_model
developer_guide
references/bibliography
api_reference
```
