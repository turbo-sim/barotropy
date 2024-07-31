

## LAES Turbine notes

- The spinodal line of nitrogen approximately follows the line of 10% vapor quality
- I computed the barotropic model properties for different values of phase-change t onset ranging from 0% to 8% (before spinodal) as detailed in the Excel file
- I generated fitted polynomials for these properties
  - The degree of the polynomial for the transition region should be relatively small (maybe 4).
  - Otherwise there could be numerical roundoff problems in CFX/Fluent (is it because of single-precision?)
  - The degree of the polynomials for the single-phase and equilibrium regions can be high to achieve good agreement with CoolProp data
- Both Fluent and CFX expressions are working (Wednesday 31-07-2024)
