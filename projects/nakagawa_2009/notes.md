

## Friday 22-11-2024 nozzle simulations

- We are running simulations on two cases:
  - Nozzle 2, expansion from 91 bar. Expansion from close to the critical enrtropy. It is not expected to exhibit non-equilibrium phase change.
  - Nozzle 2, expansion from 11 bar. Expansion farther away from critical entropy. It is expected to exhibit non-equilibrium phase change.
- The roughness of the nozzle walls is uncertaint. We conducted a sensitivity analysis on the roughness height for the 91 bar expansion to determine a suitable value that results in pressure distribution predictions that agree with the experimental data.

## Saturday 23-11-2024 nozzle simulations
- I ran a inviscid simulation to assess the losses due to friction effects. The losses in this nozzle are extremely high due to the small width and large influence of surface roughness. The efficiency of the nozzle expansions is about 10 %, where 0 % is an isenthalpic expansion (valve) while 100 % is a perfect isentropic process.
- I conducted a sensitivity analysis on the efficiency of the polytropic process but this factor only has a modest effect on the density of the fluid and does not lead to significant differences in results. Perhaps the reason is that the density of supercritical CO2 is very high?
  - I did a preliminary sensitivity analysis of the influence of the onset of phase change (cases 400s):
  - I obtained unexpected results that were difficult to explain. The reason why I was getting weird results is that I inadvertedly set the surface roughness to zero, rather than the value 5 micrometers.
  - I did not realize about the problem with the surface roughness until hours later. Before realizing, this I thought the problem could be about the settings of the metastable and that widening the q_transition variable could help
  - I designed cases 500s and 600s where Q_onset is small and Q_width increases to from 0.02 to 0.08. I also did not get good results in these cases, but I realized about the surface rouhness problems


## Sunday 24-11-2024 nozzle simulations
- I organized the notes with my thoughts about work so far
- I have to create a new output_old folder with all these results, and make nice graphs to organize my mind. Then do the new simulations


## Todo
- Rerun:
  - Roughness sensitivity for both cases
  - Mesh convergence for the best cases
- [x] Prepare boundary conditions for Jiawei
- [x] Change SNOPT settings for Srinivas (it did not have an impact)
- Paper Simone finish
- Ppaer Bartosz finish intro
- Email to Amit
- Cleanup barotropy new features and push


## Improvements to barotropy
I improved my python code to:
- Handle inviscid and viscous cases
- Print the residuals for any number of residual variables
- Improved the barotropy calculations to make them more robust:
  - Add the custom solver as fallback option for the equilibrium property calculations.
  - Added a recalculation of the inlet state as initial guess just before the recalculation of the fluid states during ODE postprocessing.
  - Increased the tolerance of the ODE solver from 1e-9 to 1e-12 to force tinier steps and a higher resolution of the polytropic process. Having a larger number of datapoints is advantagenous for the polynomial fitting, specially in the sigmoid region, which needs to be fitted with 3-to-6 order polynomial and requires multiple datapoints to avoid ill-conditioned polyfitting.
  - Addded the option to provide the degrees of the polynomials as a single scalar or as a list of scalarts to the BarotropicModel class. In addition, I added a smart strategy to automatically reduce the order of the polynomials when the number of datapoints in a fitting segment is smaller than the degree of the polynomial to be fitted.



## To do in the future
I have to check if the inviscid simulation agrees well the the IHE prediction calculated in the following way:
- Maximize the mass flux function. This is the product of density and velocity, with velocity calculated from the inlet stagnation enthalpy and entropy and the local density
- Determine the Mach number at the thoat for this maximum mass function value
- Solve an algebraic equation for the mass balance between the throat and the exit. Again, the velocity is computed at the exit density using the inlet stagnation enthalpy and entropy.
- The Mach number and velocity at the exit of the nozzle can be computed

## Other things I did
- Power cycle club on Monday
- Meeting with Andrea to discuss the need for rounding the nozzles for 1D simulations
- Coordinating with Alberto for mesh generation
- Helped Srinivas with the turboflow simulation settings and design of experiments
- Helped Jiawei with the KAIST compressor simulations, plus sent email to KAIST
- Discussion with Lasse about PhD thesis
- Meeting with Ivan about bachelor project
- Prepare little example for Jiawei to to the assignments
- MAN activities beyond simulations
  - Had a look at the blade modeler
  - Tour around the manufacturing facility
  - Running with Sebastiano
  - Afterwork on Wednesday

