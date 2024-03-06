# nitrogen_P6_00bar_T95_00K case

## Lessons learned to achieve good convergence
- Use a very low under-relaxation factor for density
	- Not needed for simulations with constant density
	- Not needed for simulations with linear density variation
	- Needed(?) for simulations with smooth density changes at the saturation point
	- Needed for simulations with non-smooth density changes at the saturation point
- Use high-order under-relaxation when using second order discretization.
- Use the differentiable gradient limiter for more monotonous behavior (not critical). High-order under-relaxation is essential even when using the differentiable limiter.
- SST k-omega model only:
	- Switch off the "Production Limiter" setting
	- Switch off the "Production Kato-Launder" setting
	- Switch off the "Low-Reynolds correction" setting
- Do not use the PRESTO! pressure interpolation scheme

## Lessons learned about things that are not important
- The time scale factor and length scale method do not have a big impact on the convergence rate and the final value of the residuals. However, setting a big time-scale factor can cause instability, so it is sensible to use a conservative value.
- The multigrid settings did not have a strong impact on the final residual level
- Changing the explicit under-relaxation factors of pressure, momentum and turbulence quantities does not have a strong impact on the final residuals. Perhaps the convergence rate can be accelerated tuning these factors. I observed that the residual cycles around the linear downwards trend depend on these explicit under-relaxation terms.
- I could not get convergence when using the "standard steady-state solver" rather than the "pseudo-time steady-state solver"


## Spalart-Allmaras turbulence model
### First order discretization
- I observed good convergence when using first order discretization
- The solution quickly reaches the point where the normalized & scaled residuals become flat at a small value. The magnitude the the residuals depends on two factors:
	- Weak dependence on the time factor or the length scale of the problem. Changing any of these two variables changes the time step of the simulation and the final value of the residuals. The final residual value does not have a strong dependence on this factors.
	- Strong dependence on the under-relaxation factor for density. It turns out that the final residual level can be reduced by reducing the under-relaxation factor for density to very low values (e.g., 1e-3, 1e-6, 1e-9). I believe that this is happening because the the flow field (pressure along the axis) crosses the saturation point several times leading to abrupt changes in the density field. As a result, it is necessary to update the density very slowly to ensure a smooth convergence to the flow solution. **This should be verified running a simulation with a smooth variation of density**
- The key to right convergence is a low under-relaxation factor for density, much more than the time or length factors used to determine the time step
- The figure below shows the convergence of the residuals when using the Spalart-Allmaras turbulence model. The solution was stopped several times to tighten the under-relaxation factor for density. As a result, we can see several levels/steps of convergence. ![[Pasted image 20230729151135.png]]
- It is not effective to reduce the "pseudo-time explicit under-relaxation factors" for the other quantities, including pressure, momentum and turbulence variables:
	- The solution is not stable when the relaxation factors are too high
	- The solution convergence rate is slow when the relaxation factors are too small. The solution may even be unstable (or wiggly) when the under-relaxation of the turbulence quantities are small (perhaps because the flow and turbulence model equations are not solved together, but in a segregated way)
- I did not get convergence with the segregated pressure-based solvers (SIMPLE, PISO)

### Second order discretization
- I also observed good convergence when using second order discretization schemes.
- I did not need to initialize the solution from a converged first-order simulation, or use high-order under-relaxation factors to achieve good convergence-
- The final residual level was reduced as the under-relaxation factor for density was tightened ![[Pasted image 20230729153531.png]]

## SST k-omega turbulence model
### First order discretization
- I initially struggled to get tight convergence when using the SST turbulence model, even when using first-order discretization for the flow and the turbulence model.
- Then I tried to switch-off the "low-Reynolds" correction and the "production limiter" options and I was able to get good convergence. This suggests that one of these options (or perhaps both) are introducing numerical difficulties to solve the equations. (further investigation below)![[Pasted image 20230729162211.png]]

### Second order discretization 
- Low-Reynolds correction and Production options were switched off.
- **Hot initialization**: I switched to second order from a converged first-order simulation and I was able to achieve very low residuals (no screenshot).
- **Cold initialization**: I tried to re-run the simulation starting from a "hybrid initialization" flow-field, but the residuals were not reduced to very low levels (the residuals stalled and showed fluctuations). ![[Pasted image 20230730012938.png]]
- **High-order term relaxation**: I tried to run the case using "high-order" under-relaxation with a factor of 0.25 and the simulations converged reliably:
	- Second order for flow and first order for turbulence model ![[Pasted image 20230729204046.png]]
	- Second order both for the flow and turbulence model![[Pasted image 20230730011716.png]]
	- The residuals decrease following a linear downwards trend with a fluctuating component. I also confirmed that reducing the pseudo-time explicit under-relaxation factor for density leads to lower final residuals values (tighter convergence)
- **Differentiable limiter (part 1)**: I tried to activate the differentiable gradient limiter and obtained more reliable convergence to very low residuals. Note that the solver starts the linear downwards trend towards full convergence at iteration ~200 (instead of iteration ~600 for the standard limiter). This is probably because the solver struggles with the non-differentiable points of the standard limiter ![[Pasted image 20230730143350.png]]
- **Differentiable limiter (part 2)**: I also tried to use the differentiable gradient limiter without using the high-order term under-relaxation factors. The solution did not converge reliably to very low residuals. This indicates that under-relaxing the high order terms is essential (the solution did not converge to very low residuals after enabling higher order term relaxation) ![[Pasted image 20230730144345.png]]
- **Differentiable limiter (part 3)**: I tried to use the differentiable limiter with the "cell-to-node limiting" setting. In this problem the convergence was similar to the case with "cell-to-face limiting". Perhaps the results would be different for a stretched mesh. ![[Pasted image 20230730151053.png]]
- **PRESTO! pressure interpolation**: I tried to run the simulations using the "PRESTO!" pressure scheme instead of the "Second Order" pressure scheme. The solution did not converge reliably to very low residuals. This suggest that the "Second Order" pressure scheme is preferable to achieve very tight convergence. The behavior of the residuals is still erratic even after switching to the "Second Order" pressure scheme. ![[Pasted image 20230730152144.png]]
- **Time factor sensitivity**: I tried to initialize the solution with a high time factor, but the solution diverges. I tried to increase the time-step factor in the middle of the simulation, but the convergence rate did not improve significantly. Based on these observations it seems that choosing the conservative length scale and setting a time factor of 0.50-0.75 leads to reliable convergence (piano piano si va lontano)
- **Production Limiter**: I tried to run the simulation activating the "Production Limiter" setting of the SST k-omega turbulence model. It turns out that the solution cannot converge to very low residuals when this option is active. This is probably because this setting introduces a non-differentiable limiter (minimum function) into the turbulence model. The figure below illustrates how the solution stalls when the production limiter is active, and how it converges reliably to very low residuals when this setting is switched off. ![[Pasted image 20230730124528.png]]
- **Production Kato-Launder**: I tried to run the simulation activating the "Production Kato-Launder" setting of the SST k-omega model. It turns out that the solution cannot converge to very low residuals when this option is active. I am not sure about the reason for this behavior. The figure below illustrates how the solution stalls when the Kato-Launder production is active. The solution remained stalled even after switching off the Kato-Launder setting (iteration ~1600) ![[Pasted image 20230730131944.png]]
- **Low-Reynolds correction**: I tried to run the simulation activating the "Low-Reynolds correction" setting of the SST k-omega model. It turns out that the solution cannot converge to very low residuals when this option is active. I am not sure about the reason for this behavior. The figure below illustrates how the solution stalls when the Low-Reynolds correction is active. The solution converges reliable to very low residuals when this setting is switched off![[Pasted image 20230730140315.png]]
- **Constant density case**: I ran a case setting a constant density (liquid nitrogen) for the simulation. In this case, the density discretization method disappeared from the "Spatial Discretization" box of the "Solution Methods" task page (as expected). The convergence of the flow towards very low residuals was smooth. This suggests that the change in density slope when the fluid crosses the saturation point is a challenge for the solver. ![[Pasted image 20230730175049.png]]
- **Linear density case**: I ran a case setting a linear density variation for the simulation. In this case, I set the density discretization method to second-order upwind and the explicit under-relaxation factor for density to 1.00.  The convergence of the flow towards very low residuals was smooth. This confirms that setting a very low under-relaxation for density is not necessary under normal circumstances, but it is an effective way to help the solver converge in cases when the density is not a continuous/smooth function (such as in the barotropic model with two-phase flow).![[Pasted image 20230730181534.png]]
- **Smoothed density polynomials**: I ran a case using the isentropic model for density. Instead of using a piece-wise polynomial function I blended the single-phase and two-phase polynomials with a sigmoid function (hyperbolic tangent) to avoid having an abrupt change in density slope at the saturation point.
	- The first simulation with the smoothed polynomial diverged. I suspect that this is because of the increase of density with pressure (imaginary speed of sound).
	- Then I created another expression with a more aggressive blending to limit the overshoot (positive slope) at the saturation point. The flow converged smoothly when using the an under-relaxation factor of 1e-9 for the density.
	- Characteristics of the residuals:
		- Faster convergence rate than with the if-statement
		- Similar final residual level as with the if-statement
		- The amplitude of the residual wiggles was smaller.
	- This suggest that it is possible to achieve better convergence when using a density function with continuous derivatives at the saturation point. ![[Pasted image 20230730184756.png]]
	- Then I changed the density-under-relaxation factor to 1.00 and re-ran the simulations. The residual convergence was still good, although the final residual levels were not as low as the the previous case (it would take more iterations to reach similar residual levels)![[Pasted image 20230730190459.png]]
	- In addition, I checked that the speed of sound is going bananas, even though the simulation converges (which is the important thing). I am not sure how to interpret that the internal speed of sound is wrong, but it does not look good.![[Pasted image 20230730190743.png]]
- **p-norm density polynomials**: I ran the simulations using a p-norm approximation to create a smooth density function at the saturation point. I used a under-relaxation factor of 1.00 for density. In this case the residuals almost stalled and the convergence rate was very small with high-amplitude wiggles.![[Pasted image 20230730192344.png]]
- Then I tightened the density under-relaxation factor to 1e-9 and obtained a better convergence trend. The residuals decreased rapidly to very low values and remained constant. The speed of sound for this case was normal (as expected)
![[Pasted image 20230730195011.png]] ![[Pasted image 20230730195038.png]]
- Using an moderate value of under-relaxation factor, like 0.25, helps the convergence, but the residuals eventually stall. This suggests that a very small under-relaxation factor is required for tight convergence. ![[Pasted image 20230730221316.png]]
- These results suggest that:
	- If-statement requires a very small density under-relaxation factor (not a problem)
	- Blending with hyperbolic tangent is associated with an overshoot of density that leads to unphysical sound speeds. Despite this limitation the simulation converges reliably to the correct solution (even with no density under-relaxation)
	- The p-norm approximation leads to a smooth density function with no overshoots and physical values for the speed of sound. However, a tight density under-relaxation factor is required for convergence?


If-model, 1e-9 under-relaxation

![[Pasted image 20230730225952.png]]
![[Pasted image 20230730230043.png]]

## Related information from theory guide
- **Discretization methods**
	- The pressure-based solver has the following discretization options to calculate the flow variables at the faces of the control volumes
	- Mass-flux discretization (product density normal velocity):
		- Distance-weighted high-order velocity interpolation
		- Momentum-weighted high-order velocity interpolation
		- Both methods compute a face mass flux and apply a pressure correction similar to that proposed by Rhie and Chow to prevent checkerboarding
	- Pressure interpolation
		- Linear
		- Standard
		- Body-force weighted
		- Second order (central differencing with gradients)
		- PRESTO! (pressure staggering option)
	- Density interpolation
		- Arithmetic average for incompressible
		- Any upwind scheme for compressible
	- Scalar transport (including u-v-w, temperature, and turbulence quantities)
		- First order upwind (does not need gradient)
		- Second order upwind (does need gradient)
		- QUICK (blend of second order upwind and central interpolation)
		- MUSCL (blend of second order upwind and central differencing)
- **Evaluation of gradients**
	- Green-Gauss cell-based (cheap and innacurate)
	- Green-Gauss node-based (expensive and accurate)
	- Least-squares (moderate cost and accuracy)
- **Gradient limiters**
	- Section 24.3.4. of theory guide and Section 35.8 of the user's guide
	- Gradient Limiters, also known as slope limiters, are used on the second-order upwind scheme to prevent spurious oscillations, which would otherwise appear in the solution flow field near shocks, discontinuities, or near rapid local changes in the flow field.
		- Standard limiter (non-differentiable, minmod function))
		- Multidimensional limiter (non-differentiable)
		- Differentiable limiter
	- One disadvantage with non-differentiable limiters is that they tend to stall the apparent residual’s convergence after a few orders of reduction in residual magnitude. This annoying behavior can be directly traced to the non-differentiable nature of the limiting functions. Therefore, the differentiable limiter uses a smooth function to impose the monotonicity condition while allowing the residuals to converge.
	- Limiter directions
		- Cell-to-face limiting: limited value evaluated at the cell face center
		- Cell-to-node limiting: limited value evaluated at the vertices of a cell (more conservative and improved monotonicity-robustness when computing flows with strong gradients on unstructured grids)
- **High Order Term Relaxation (HOTR)**
	- **Section 24.3.1.9 of the theory guide:** Higher order schemes can be written as a first-order scheme plus additional terms for the higher order scheme. The higher-order relaxation can be applied to these additional terms. The under-relaxation of high order terms follows the standard formulation for any generic property $\phi$	$$\phi_{new} = \phi_{old} + f(\phi_{high\,order} - \phi_{old}) $$where $\phi_{high\,order}$ is the new value of the property computed by the high-order scheme (e.g., second order upwind method) and $f$ is the under-relaxation factor between zero and one. The default value of for steady-state cases is 0.25 and for transient cases is 0.75. The same factor is applied to all equations solved.
	- **Section 35.2.5. of the user's guide:** The purpose of the relaxation of high order terms is to improve the startup and the general solution behavior of flow simulations when higher order spatial discretizations are used. It has also shown to **prevent convergence stalling** in some cases. Such high-order terms can be of significant importance in certain cases and lead to numerical instabilities. This is particularly true at aggressive solution settings. In such cases, high order relaxation is a useful strategy to minimize your interaction during the solution. This can be an effective alternative to starting the solution first order, then switching to second order spatial discretization at a later stage.
- **Solver under-relaxation factors**
	- The steady-state pressure-based solver uses two types of under-relaxation when the "Pseudo time method" is set to "Global time step"
		- Implicit relaxation of the equations. Pseudo-time integration is an advanced for of implicit under-relaxation of the equations where the time-step is dynamically computed depending on the flow solution
		- In addition, it is possible to specify explicit under-relaxation of flow variables to help the convergence of the problem (including, pressure, density, momentum, and turbulence quantities)
	- Define the explicit under-relaxation parameters of the flow variables in the "Solution Controls" task page (the default values are usually okay)
	- Set the "Pseudo time settings" in the "Run Calculation" task page
		- Number of iterations
		- Time scale factor
		- Length scale method
- **Pseudo-time settings**
	- The coupled pressure-based solver tries to converge to a steady-state solution using pseudo-time integration with a global time step.
	- The time step for the integration is selected as the minimum of several time scales
		- Convection
		- Dynamic
		- Buoyancy
		- Gravitational
		- Rotational
		- Acoustic
	- Most time scales are defined as a length scale divided by a representative velocity scale.
	- The length scale of the problem can be defined as:
		- Aggressive time scale: maximum of the volumetric and domain lengths
		- Conservative time scale: minimum of the volumetric and domain lengths
		- User-defined time scale: arbitrary length scale (useful for far-field domains)
	- Volumetric length scale: Cubic root of the volume of the domain (3D simulations) or square root of the area of the domain (2D simulations)
	- Domain length scale: Maximum of the x, y, or z extents of the flow domain.
- **Production limiters for two-equation turbulence models**
	- Section 4.20 of the theory guide.
	- A disadvantage of standard two-equation turbulence models (both $k-\epsilon$ and $k-\omega$) models is the excessive generation of the turbulence energy, in the vicinity of stagnation points.
	- In order to avoid the buildup of turbulent kinetic energy in the stagnation regions, the production term in the turbulence equations can be limited by:
		- Menter approach (min function as a limiter)
		- Kato and Launder approach (using vorticity instead of strain rate)
	- I experienced that both production limiters prevent a tight convergence of the problem to very low residual levels. Since I do not know what is the impact of the production limiters on a turbomachinery case (I could assess later on) it seems reasonable to switch these settings of to achieve good convergence properties.
- **Low-Reynolds Number Correction**
	- Section 4.4.1.3.1 of the theory guide
	- This correction introduces a less-than-one coefficient that damps the turbulence viscosity when the turbulence Reynolds number is low.
	- EI experienced that this setting prevents a tight convergence of the residuals.

