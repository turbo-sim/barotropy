
Add better polynomial fit (custom optimization?)
Provide order of each polynomial as list / check number of elements is correct /prevent negative slope.

## Software package
- Translate the function for the two-component barotropic model
- Translate the functions for plotting the two-component barotropic model
- Verify the polynomials for vapor quality and void fraction
- Create function to compute the pseudocritical point/line
- Function for s-Q function calls (fzero or fsolve?)
- Function for sonic point calculation working on two-phase?
- Delete matlab code once it is all traslated
- Make demos scripts about:
  - Single fluid model (single-phase region, two-phase equilibrium, two-phase blending)
    - Well posed (CO2). Ill-posed (nitrogen)
  - Two-component model (easy modelig, no metastable)


## Testing
- Create tests for basic fluid properties
- Create tests intermediate-complexity functions
  - rhoT HEOS evaluation
  - custom solver (check it agrees with coolprop)
  - spinodal point calculation (regression test)
  - spinodal line calculation (regression)
  - Subcooling/superheating
  - Supersaturation
  - States to dict, states to dict 2D
- Create tests for barotropy
  - Isentropic case, verify entropy conserved (both 1-fluid and 2-fluid)
  - Zero efficiency case, verify enthalpy conserved (both 1-fluid and 2-fluid)
  - Test that polynomials have error below threshold for several properties
  - Test that polynomials are continuous at the endpoints for several properties
  - Regression test on CFX and Fluent expressions for equilibrium, blending, and two-component cases

## Documentation
- Make installation instructions:
  - Basic (pip install)
  - Recommended (conda environment)
  - Developer
- Add developer instructions (installation, CI/CD, testing)
- Add tutorials
- Add theory about barotropic model
- Add API documentation about the standalone blending.
- Add API documentation about two-component model
- Add API documentation of fluent automation? (clean it?)




## Simulations

- LAES cases
- Lettieri CO2 data (important!)
- Zhang steam nozzle
- Simoneau and Hendriks nitrogen nozzle data (critical!)



## Critical
- [x] Move the calculation of the blended states as a function of q_onset and q_transition to the properties calculation. Not within the ODE
- [x] Generalize the polytropic ODE behavior to dividing in several steps with linspace entropy and logspace pressure.
- [x] Make a test script showing the generalized quality calculation
- [x] Make a test script showing the degree of superheating-subcooling calculation
- [x] Define polynomial order
- [x] Check matlab function create_barotropic_model()
- [x] Check that both single-component and two-component functions work
- [x] Postprocess the integrated function according to the type of multiphase model:
	- [x] If metastable (just one polynomial)
	- [x] If equilibrium 2 polynomials
	- [x] If blending 3 polynomials
- [x] Fit polynomials for each branch and ensure continuity at the endpoints
- [x] Check matlab function plot_barotropic model, do I have to translate?
- [x] Translate functions for
	- [x] Fluent expressions
	- [x] CFX expressions
- [x] Functions to calculate the spinodal line in the rho-T diagram
	- [x] Optimization problem formulation?



## Medium critical
- Compute function to calculate pseudocritical point and line and make script illustrating calculations
- Make script illustrating phase diagram
- Check if the Q-s solver works well as 2 variables and compare to the function I wrote in MATLAB and the one I wrote in Python in KAIST

## Not critical
- [ ] Functions to plot T-s and p-s diagrams not specific for sCO2
- [ ] Move compute_inlet_state() to sCO2 specific project
- [ ] figure out where is the function compute_outlet_state_isentropic() used
- [ ] utilities to convert from actual-to-equivalent turbomachinery inlet conditions
- [ ] Isentropic nozzle calculations. figure out where they are used



## Long term
- Redo Hendriks cases with the smoothed model
- Fix the automation of the Elliot 1982
- Carry out simulations for Elliot 1968
- Carry out simulations for the Nakakawa/Lettieri cases
