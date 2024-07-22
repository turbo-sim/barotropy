# Notes on condensation models



**Classical nucleation theory and its application to condensing steam flow calculations**


The paper discusses the classical theory of the homogeneous nucleation of water droplets from supersaturated vapour and its application in predicting condensation in steam nozzles. The first part consists of a review of classical nucleation theory, focusing on the many modifications made to the original Becker –Do ̈ring theory and providing some new insights into recent developments. **It is concluded that the predictive accuracy required for engineering calculations is not yet attainable with a theory derived from first principles.**

**Semi-analytical model for the prediction of the Wilson point for homogeneously condensing steam flows**

It was established that the location of the nucleation onset is particularly sensitive to the steam heat capacity ratio γ (Bakhtar et al., 2005). Especially when approaching high-pressures, small variations of this parameter can increase the discrepancy between the theoretical solution and the measurements. Moreover, the surface tension is usually affected by considerable uncertainties, and existing correlations do not take into account any droplet curvature effects (Lai and Kadambi, 1990). As a consequence, due to the exponential dependence of Js on σs, the theoretical Wilson pressure Pw and the droplet properties are far from the experimental data. Therefore, empirical coefficients are customarily introduced (Bakhtar et al., 2005; Lai and Kadambi, 1990) to correct the parameters σs, Js, Gs in order to reach a better accuracy. Following Azzini and Pini (2017), σs, Js, Gs are then multiplied by an empirical factor, yielding to

**Numerical study of heterogeneous condensation in the de Laval nozzle to guide the compressor performance optimization in a compressed air energy storage system**

- Condensation of water in humid air
  - Homogeneous
  - Heterogeneous
- One-fluid mixture model
  - Mass, momentum equations for mixture
  - Energy equation with source term for condensation enthalpy
  - Turbulence equations
- 3 additional scalar transport equations:
  1. Number of homogeneous condensation droplets per unit of mass
  2. Mass fraction of homogenous condensation phase
  3. Mass fraction of heterogeneous condensation phase
- The additional equations have source terms that have to be modeled:
  1. Directly given by the nucleation rate $J$
  2. Sum of two terms, one involving the nucleation rate and the critical radius and another involving the droplet growth
  3. One term involving the droplet growth
- Additional fundamental physical models required:
  - Nucleation rate
  - Critical radius
  - Droplet growth rate
- **Critical radius** is modeled according to the Kelvin (?) equation
- **Nucleation rate** is obtained from molecular theory of gases for an ideal gas, and there are corrections for non-isothermal cases. It depends on the surface ternsion.
- **Droplet growth rate** is based on semi-empirical considerations. Depends on wether the interaction between gas and bubbles hold the continuum hypothesis or the interaction should be described using molecular dynamics considerations.
  - Variety of models for continuum: Gyamathy, Fuchs, Young, Gill.
  - Knudsen-Hertz model for rarefied cases
  - Blending model in-between
- Still not clear to me how they compute the droplet radius at the CFD level.




## My thoughts
- Condensation models require several sub-models that are a combination of theoretical considerations based on the ideal-gas assumption and empirical models
- The accuracy of the analytical foundation close to the critical point is uncertain as the fluid does not behave as an ideal gas
- There are several additional equations that depend on several fitting constants to achieve good matching against experimental data.
- The solution of the CFD problem is sensitive to the choice of condensation model and to the fitting constants of each submodel
- Authors often do not mention anything about the fitting of the constants.
- At the end, the model can only be accurate if the fitting constants are tuned against experimental data
- Tuning these constants is not trivial because they often convey very little physical significance and it is difficult to undertand how the affect the simulation results onless one is a subject matter expert.
- As of now, accurate predictions cannot be achieved unless the model fitting parameters are fitted against experimental data


## Alternative approach
- Condensation models not suited for flashing
- Phase change adds new modeling equations, higher computational cost, more difficult to converge.
- Considering this landscape, we deem it sensible to develop a different modeling strategy based a simplied description of the phase change.
- Barotropic model where metastable effects are considered until the fluid reaches the wilson line. Then the metastable and equilibrium properties are blended for a certain metastable penetration depth.
- The use has two tuning parameters, the location of the onset of condensation (the location of the Wilson line) in terms of degree of supersaturation, or as a certain vapor quality. And penetration depth where the metastable effects fade out (specified as the width of a smooth-step function)
- This approach relies on only two intuitive fitting parameters and was able to replicate the behavior of the fluid expansion to the right and to the left of the saturation line
- Suitable for turbomachinery applications, where the fluid process can be approximated by a polytropic process and where the focus is on the overall performance of the machine, not on the behavior of each phase and their interactions. 
