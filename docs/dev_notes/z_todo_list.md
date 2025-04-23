
# Changelog


## Sunday 24-11-2024
I improved my python code to:
- Handle inviscid and viscous cases
- Print the residuals for any number of residual variables
- Improved the barotropy calculations to make them more robust:
  - Add the custom solver as fallback option for the equilibrium property calculations.
  - Added a recalculation of the inlet state as initial guess just before the recalculation of the fluid states during ODE postprocessing.
  - Increased the tolerance of the ODE solver from 1e-9 to 1e-12 to force tinier steps and a higher resolution of the polytropic process. Having a larger number of datapoints is advantagenous for the polynomial fitting, specially in the sigmoid region, which needs to be fitted with 3-to-6 order polynomial and requires multiple datapoints to avoid ill-conditioned polyfitting.
  - Addded the option to provide the degrees of the polynomials as a single scalar or as a list of scalarts to the BarotropicModel class. In addition, I added a smart strategy to automatically reduce the order of the polynomials when the number of datapoints in a fitting segment is smaller than the degree of the polynomial to be fitted.


# Planned changes

- [ ] Add better polynomial fit. Perhaps not needed now that the ODE solver tolerance has been improved.
- [ ] Prevent negative slope. Perhaps not needed now that the ODE solver tolerance has been improved.
- [ ] Provide order of each polynomial as list / check number of elements is correct



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




# (OLD MATLAB NOTES)
## Progress made
- I created a MATLAB function to generate polynomials for the barotropic model
	- Single phase region
	- Metastable region
	- Two-phase region (equilibrium)
- Compute single-phase or metastable properties using the Helmholtz energy (rho, T)
	- Compute single-phase or metastable properties using T-s function calls (1D root)
	- Compute single-phase or metastable properties using p-s function calls (2D root)
- Compute liquid spinodal point as a function of temperature (1D root)
- Compute vapor spinodal point as a function of temperature (1D root)
- Compute spinodal point as a function of entropy (2D root - robust initial guess)
- Plotted the spinodal lines and isotherm lines in thermodynamic diagrams
- Plotted thermodynamic properties along an isentropic expansion
	- 3 regions: single-phase, metastable, two-phase
	- The setup with 3 regions is necessary for the case when extrapolating to metastable states, otherwise the accuracy of the polynomial fit is poor.
	- Regions separated by saturation point and spinodal point
	- Polynomial fitted to the thermodynamic variables in the regions. The polynomials are fitted using the reduced pressure as independent variable to have a good scaling for the problem and avoid ill-conditioning.
	- Plot comparing the polynomial and the calculation and show relative error
	- Evaluation of the polynomial in two different ways:
		- Piecewise function with (possibly) discontinuities (Matlab uses the Horner's rule for polynomial evaluation: N sums + N multiplications)
		- Blended function with a sigmoid activation function (hyperbolic tangent). This leads to a single continuous and smooth function (perhaps better properties for CFD simulation)
	- Export the polynomial coefficients to Excel and Fluent
	- I learned that there is a discontinuity in temperature (pressure) when the pressure (temperature) that sweeps the isentrope crosses the spinodal line. This is because the properties are calculated according to equilibrium beyond the spinodal line (rather than from the Helmholtz energy equation). I believe that it makes more sense to sweep the isentrope using pressure as variable and have a discontinuity in temperature (the temperature field is anyway not computed in the CFD barotropic model). Perhaps it makes sense to use the blended model in these cases even if it is not physical.
- Created function to compute the pseudocritical line
	- Locii of points with maximum isobaric heat capacity at a certain temperatrue
	- Computation based on solving a 1D optimization problem (fminunc)
	- Re-use the previous pseudocritical point as initial guess to speed up computations and improve convergence
	- The pseudocritical line is close the the line of critical density, but it is not exactly the same

## Work to do
- Finish documentation
- Refactor code to simplify main computation function

# Nitrogen nozzle simulation
## Progress made
- Automatic isentrope calculation
- Automatic polynomial generation
- Automatic expression generation
- Extrapolation to negative pressures to ensure convergence
- Convergence nitrogen case with Fluent
- I checked the flow is choked as it should. I decreased the pressure at the exit from 2.5 bar to 2 bar and got the same mass flow rate, good indication that the flow is choked :)
- I checked the speed of sound definition is consistent
	- Yes, the speed of sound computed with the expression agrees well with the speed of sound computed internally by fluent
	- I tested with meld and there is only a few decimals of difference between the two speeds of sound, probably because Ansys is computing the speed of sound internally by finite differences.
	- Smoothing does not give improvements and sometimes might crash if the isentropic bulk modulus is negative (imaginary speed of sound)
- Finding the correct settings for nice convergence
- The simulation converges for:
	- Constant density
	- Linear with small slope
	- Linear with big slope
	- High compressibility with negative pressure safeguard
- Can we run the density-based solver without energy equation? Probably not, I get SIGSEGV error


## Work to do:
- Run new simulation closer to the critical point (good convergence for some cases)
- Check sensitivity with respect to:
	- y+ value
	- turbulence intensity boundary condition
	- turbulence viscosity ratio boundary condition
	- roughness of the call
- Improve mesh?
- Better post-processing
- Use pyfluent to automate parametric analyses?







# Fluent automation

## 23.09.2023
Improved the functionality of my Fluent python functions:
- Print transcript to python console in real time
- Plot flow variables for 1D nozzle
- Export many integrated quantities to .out file and postprocess this file with pandas to retrieve dataframe or excel file 
- Create function to run fluent simulations
- Implement function to gradually change the iteration settings


## 26.09.2023

- Refactored the automate simulations script to be more compact:
  - Use Excel file with index and tag to schedule all simulations
  - The iteration is now based on the excel file content (instead of being a simple iteration loop over some indices).
  - Great flexibility for sensitivity analyses of any type (e.g., changing outlet pressure)
  - Add residual tolerance termination criterion to automation
  - Improve the export of .xy files to be configurable from the automation script
- The generation/rendering of the Fluent journal was vastly improved:
  - Output journals preserve spacings and semi-colons for readibility
  - The rendering function checks if all the "template_map" keys were used and raises and error if there is any key that is un-used.
  - New function added to validate the journal and check that there are no un-rendered variables
  - New function implemented to add the barotropic expression definition commands ot the journal
  - New function implemented to add the solution strategy commands to the journal
  - New function implemented to add the export .xy file commands to the journal
  - New function added to encapsulate all the functionality related to the journal creation
  - Validations and checks added to minimize the risk of unnoticed human errors
- Refactored functions to parse fluent files and added docstrings:
  - parse_fluent_out handles the header in a better way
  - parse_fluent_xy elies less on the parse module for faster speed. 
  - parse_fluent_xy execution logic was refactored to make it easier to understand
  - added docstrings to read_residual_file
- Added functionality to plot the nozzle profiles
  - No code redundancy for different plots thanks to helper functions
  - Easily customizable with several lines per plot
  - Easily extensibel to other plots (like y+ distribution)
- Ran simulations on the mesh with inflation layers:
  - No siginificant difference with respect to the case with no inflation layers
  - Adding roughness to the simulations changes the results a bit, but does not explain the big deviation
  - Adding a lot of roughness decreases the blade speed



## To-do
- Use pyfluent for
  - Meshing
  - Scheduling simulations
- Create meshes for
  - Nozzle (internal flow)
  - Nozzle (exhaust jet to atmosphere)
  - Linear cascade
  - Turbomachine (turbogrid)