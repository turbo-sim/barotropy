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