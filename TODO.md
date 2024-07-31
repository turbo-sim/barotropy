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