## List of functions

### Fluid properties
- `compute_properties_metastable_Td`
  - Compute all thermodynamic properties from the Helmholtz energy equation
  - This function bypasses the equilibrium relations from CoolProp
  - Explicit evaluation without equation solving
- `compute_properties_metastable_ps`
  - This function bypasses the equilibrium relations from CoolProp
  - This function must solver 2 nonlinear equations to find the temperature-density inputs that matches the pressure-entropy arguments.
- `FluidCoolProp_2Phase`
  - Create object to interact with CoolProp
  - Adds functionality to compute HEM properties in the 2-phase region
  - Robust error handling in MATLAB
- `compute_spinodal_point`
  - Compute spinodal point using temperature as input variable
  - The function solves 1 nonlinear equation to get the point where the isothermal bulk modulus is zero (or a minimum with the robust formulation)
  - Functionality to find a good initial guess along the liquid or vapor lines
- `compute_spinodal_point_entropy`
  - Compute the spinodal line using entropy as input variable
  - The function solvers 2 nonlinear equations to get the point where the isothermal bulk modulus is zero and the input entropy is satisfied.
  - Initial guess based on pre-computing the entire spinodal line and using the closest point to the input entropy
- `compute_spinodal_line`
  - Sweep the spinodal line from the critical temperature to the triple temperature
  - Use previous point as initial guess to ensure quick convergence
  - The resulting spinodal in the p-s plane does not range from the triple pressure to the critical pressure. The liquid line goes to negative values and the vapor line ends before the triple pressure. This is a limitation of the equation of state
- `compute_saturation_line`
  - Compute the fluid properties along the vapor and liquid saturation lines
- `compute_saturation_point_entropy`
  - Compute the state along the saturation line that has the entropy specified as parameter
  - This function is necessary because the CoolProp Q-Smass function calls do not work
- `compute_pseudocritical_point`
  - Compute the point in which the isobaric heat capacity is a local maximum as a function of temperature
  - Only relevant for temperatures above the critical temperature
- `compute_pseudocritical_line`
  - Compute the locii of points in which the isobaric heat capacity is a local maximum

### Barotropic model
- `compute_barotropic_model_segment`
  - Compute properties along isentrope using p-s function calls
  - Fit a polynomial to data
- `create_barotropic_model`
  - Create barotropic model by defining several isentrope segments
  - The number of segments (1-phase, metastable, 2-phase) depends on the input entropy level
  - The function considers different cases with if-statements
- `evaluate_isentrope_polynomials`
  - Concatenates the polynomials segments at the connecting points
  - Uses Horners method for polynomial evaluation
  - Has the option to use a sigmoid blending to smooth the breakpoints
- `export_fluent_expressions`
  - Creates fluent expressions for the piecewise polynomial

### Plotting

- `set_plot_options`
  - Default settings for pretty plots
- `plot_phase_diagram`
  - Plot the phase envelope for any thermodynamic input pair
  - Include saturation line or not
  - Include spinodal line or not
  - Include vapor quality isolines or not
  - Include pseudo-critical line or not
- `plot_barotropic_process_ps`
  - Plot the isentrope segments in the p-s diagram
  - Include all options from `plot_phase_diagram`
- `plot_barotropic_polynomials`
  - Plot the piwcewise polynomial of the fluid properties as a function of pressure
  - Plot the fitting error as a function of pressure


