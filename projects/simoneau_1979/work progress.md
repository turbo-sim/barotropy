  

### Work Progress
- [x] Create Barotropic Fluid model
	- [x] MATLAB script
	- [x] Equilibrium/metastable properties
	- [x] Robust spinodal computations
	- [x] Export Fluent expressions
- [x] Run simulation case away from critical point
- [x] Investigate which are the best settings for the solver
- [x] Harmonize the name of the experimental datafiles (python script)
- [x] Create MATLAB script to plot all the experimental cases on a p-s diagram (read excel)
- [x] Create MATLAB script to automate the generation of the barotropic model for each case
	- [x] User input is just the case numbers. Example: \[92, 98, 101\]
	- [x] Save polynomial expressions and plots in a folder for each case. For example `simulations/case_092_nitrogen_T130_00K_P42_50bar/`
	- [x] Plot all the cases considered in a p-s diagram and save in the `simulations` directory
- [ ] Create python function to plot simulation vs experimental results
	- [ ] Read simulation data from `simulations/case_092_nitrogen_T130_00K_P42_50bar/`
	- [ ] Read experimental data from `experiments/case_092_nitrogen_T130_00K_P42_50bar/`

I updated the harminization script to write a neew execel file with the new .csv data files filenames

I created functions to read text and replace the characters in brackets with a rendering mapping. This is similar to what jinja engine doe sin the context of Flask HTML rendering

I added commands to add the boundary conditions
- pressure in
- pressure out
- turbulence instensity and viscosity ratio 
To the simulations.
The boundary conditions are read from the excel file, not from the filename (more precision)

Improve script to delete existing datafiles when running a new simulation, good practice?


The simulation folders are created with a unique datetime identifyer so we do not have to worry about deleting old simulation folders. All previous results are saved.



I created script to read
- [ ] 
- [ ] Run simulation close to critical point
	- [ ] Plot the simulation vs experimental results


 TODO

 Create a loop structure for:

  - [ ] Read case name
  - [ ] Find expression directory
  - [ ] process expression file to get name and value for density/viscosity
  - [ ] change directory to the correct place (location of the expression)
  - [ ] save fluent simulation data in a "results" directory

 Add wildcard characters to the template that can be replaced with:
  - [ ] expression name and value
  - [ ] assigning expression to fluid
  - [ ] setting pressure boundary conditions
  - [ ] run simulation

 Add functionality for further reading/plotting
 - [ ] Plot the residual history to ensure convergence
 - [ ] Plot comparison with experimental data
 - [ ] Gather all the mass flow rates from file and make a comparison against experimental data for all cases

