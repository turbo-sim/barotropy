
## Work progress

- I ran simulations from cases 88 to 117
- Some cases did not converge
- Problems with the cases that converged
  - Turbulence viscosity warning
    - Fluent displayed a warning indicating that the turbulence viscosity was hitting the upper limit
    - This may be the case for a situation with a shear flow as it happens after the oblique shock wave.
    - The simulation converges well without warnings after lifting the limit
    - I should do a sensitivity study comparing the three turbulence models to see if we get the same result
      - $k$-$\epsilon$
      - $k$-$\omega$ SST
      - Spalart-Allmaras
  - Unphysical fluctuations near saturation point:
    - Not happening for all cases? Very noticeable in some
    - The pressure distribution takes a strange shape close to the saturation point
    - Perhaps it is caused by the slope non-differentiability
      - Value discontinuity -> not the culprit (change in polynomial coefficients)
      - Slope discontinuity -> probably the culprit

- The results obtained so far are very promising. Further refinement/automation is required

## Functionality implemented

### Scripts
- **Automate Fluent execution**
  - Specify the cases to run by case index
  - Loop over folder structure systematically
  - Create unique results folders with timestamps
  - Read boundary conditions from excel file and impose them on the template
  - Specify solver settings (threads, iterations and timeout) on template file
  - Specify output filenames in the template file
  - Specify barotropic model expression in template file
  - Render template file with the correct values
  - Use pexpect to run Fluent as a child subprocess in the background
  - Send the journal to Fluent to configure and run the simulation
  - Finish the pexpect process when the execution is complete or the timeout is reached
  - Run the fluent cleanup file if the timeout is reached
- **Plot the pressure validation**
  - Loop over folder structure systematically
  - Read the .xy files with simulation data
  - Read the .csv files with experimental data
  - Plot the results with custom plot settings
- **Plot mass flow rate validation**
  - Loop over folde rstructure systematically
  - Read the .out files with the simulation data
  - Read the excel file to get all the experimental mass flows
  - Create two lists of experimental/simulation mass flows
  - Print them as a markdown table comparing relative/absolute deviations
  - Make a y=x plot to compare the experimental/simulation mass flows

### Functions
- Make function to render a fluent journal template in which some elements are specified between curly brackets {}
- Make function to read an expression file with the barotropic model
- Make function to read xy files as pandas dataframe
- Make function to read report .out files as pandas dataframe
- Make function to read transcript .trn files and extract residuals as pandas dataframe
  - Improve functionality to read just the lines from the last line
  - Live reading of a transcript file that is being updated
  - Live plotting of a transcript file that is being updated


## To do list
- [ ] Fix the problem with the wiggles
- [ ] Make extra plots for density, speed of sound, and mach number
- [ ] Sentitivity analysis turbulence models
- [ ] Sensitivity analysis roughness
- [ ] Sensitivity analysis number of cores
- [ ] Try PyFluent to streamline automation
- [ ] Try CFD post to make publication quality plots
- [ ] Translate the barotropic model to python (barotroPY?)
- [ ] Improve the barotropic model with blending