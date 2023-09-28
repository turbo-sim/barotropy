# Notes on "Elliot 1979" Report

## **Summary of Report Elliot 1979**
- Title: *Theory and Tests of Two-Phase Turbines*

## **Introduction**
- **Two-Phase Turbines Benefits**:
  - Steam-water mixtures
  - Partial-evaporation ORC
  - Immiscible mixtures of liquid and gas
- **Potential Applications**:
  - Geothermal power
  - Waste-heat recovery
  - Refrigerant expansion
  - Solar power
  - Bottoming cycles
- **Test Setup**:
  - Equipment for one- and two-stage turbines
  - Fluids: water-nitrogen mixture, Refrigerant 22
  - Challenges with water-steam mixture:
    - Drop sizes might be too large for homogeneous fluid assumption.
    - Liquid drops potentially deviating and impacting blades.
    - Possible formation of a liquid film.

## **Apparatus**
- **Acceleration Factors**:
  - Gas: Pressure gradient
  - Liquid: Pressure gradient and drag
- **Nozzle Dimensions**:
  - Inlet angle: converging at 10 degrees.
  - Exit angle: diverging at 2.5 degrees.
  - Nozzle throat diameter: 13.1mm.
  - Turbine test cut: 20 degrees (70-degree inlet angle).
  - Upstream nozzle diameter before ellipse: 27.6 mm.
  - Throat area changes gradually.

## **Nozzle Calculations**
- Velocity derived from the uncorrected thrust equation.
- Absence of details on thrust and flow rate measurements.
- Estimation of max droplet diameter: Weber number criterion, considering surface tension and slip velocity. Lower surface tension results in smaller droplets/bubbles.

## **Two-Phase Flows Data**
- **Steam and Water**:
  - Source: Lawrence Livermore Laboratory (LLL)
  - Reference: Alger 1975, ref 11.
- **Water and Nitrogen**:
  - Report contains experimental data.
  - Prior data from Elliot on nozzles is more relevant.
- **R22 (Refrigerant 22)**:
  - Limited data available for validating nozzle performance with this refrigerant.



# Barotropic Model as UDF
- Developed barotropic model as a C++ UDF.
  - Supports polynomials of all orders.
  - Utilizes Horner's rule for evaluation.
  - Function added to export model as user-defined-scalars for validation.
- Ran simulations comparing UDF and expressions (Elliot1979 case).
- Identical outcomes from UDF and expression-based simulations (when there are no density deviations).
- UDF allows for explicit speed of sound definition.
  - Only affects equation convergence (user manual says it "stabilizes" the solution)
  - No impact on the final solution: pressure and density fields align, despite internal speed variances between UDF and expressions.
- Expressions are simpler and more time-efficient -> they are the preferred approach

# Density under-relaxation
- Probem encountered:
  - Starting with overly high density under-relaxation prevents correct convergence.
  - Mismatch between Fluent's internal density field and barotropic model's pressure-related density prediction.
  - Confirmed mismatch not due to polynomial evaluation but density under-relaxation.
- Suggested strategy:
  - Start with high density under-relaxation (0.5-1.0).
  - If residuals stall or fluctuate, decrease density under-relaxation to enhance convergence.
  - Residual stalling often linked to shock waves or rapid gradients.
- Results/conclusion:
  - Implemented a step-wise solution approach, tightening density under-relaxation from 1 to 1e-9.
  - Confirmed that reducing density under-relaxation slows density field updates.
  - Even with small density under-relaxation, the density field still changes gradually when nearing convergence (the density field is not completely frozen)


# Nozzle thrust equation
- The paper's average velocity is based on the nozzle thrust equation: 
$$F = \dot{m}\,v $$
- This equation doesn't account for thrust due to differences in exit and backpressure when nozzle backpressure isn't adapted. Corrected equation:
  $$ F = \dot{m}\,v + A_{e}(p_e - p_b) $$
- Where $p_e$ is nozzle exit pressure and $p_b$ is backpressure (ambient pressure boundary condition).
- Infuence of inlet on thrust equation:
  - Question raised: Is thrust influenced by inlet velocity and pressure?
  - Possible answer: Likely not if inlet ports connect to flexible hoses that don't transmit reaction force. In such cases, the reaction force might only be felt in the thrust meter's stiff fittings.
- Simulation results:
  - Observed both overexpanded and underexpanded nozzle cases depending on mixture ratio:
  - To predict exit pressure-velocity (and thrust) accurately for overexpanded cases, it might be necessary to extend computational domain to ambient.
  - Some overexpanded nozzle cases, with shock waves at the boundary condition, didn't converge without reducing density under-relaxation.
  - Tests without shock waves converged well even without under-relaxation.
- Advised to adjust the computational domain to include the atmosphere downstream of the nozzle. This captures effects downstream, ensuring boundary conditions don't impact the domain's last cell and allowing for accurate modeling of shocks and expansion waves (would results change?)

## To-do
- [x] Run simulations on mesh with inflation layers
- [x] Analyze sensitivity of wall roughness model
- [ ] Create mesh with Pyfluent and add downstream domain for under-expanded case
- [ ] Run validation of the other water-nitrogen nozzle case:
	- More data including pressure distribution, mass flow, thrust for exit velocity
	- Bigger nozzle: less effects of roughness or scale
	- Crate mesh for this case
- [ ] Read the report and understand what data they present and how does the experimental apparatus works