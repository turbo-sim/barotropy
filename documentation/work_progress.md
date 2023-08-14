
# Barotropic model
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
- Run new simulation closer to the critical point
- Check sensitivity with respect to:
	- y+ value
	- turbulence intensity boundary condition
	- turbulence viscosity ratio boundary condition
	- roughness of the call
- Improve mesh?
- Better post-processing
- Use pyfluent to automate parametric analyses?



