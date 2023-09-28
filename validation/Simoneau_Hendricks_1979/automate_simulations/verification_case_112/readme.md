
### Case 112 analysis
- The solution for this case has unphysical wiggles close to the saturation point
- The wiggles in the pressure distribution remain even if the residuals are very low and flat
- I originally throught that this problem could be caused by:
  1. The small discontinuity of the density function at the saturation point
  2. The discontinuity of the density slope at the saturation point
- I corrected the first problem adjustiing the first coefficient of the polynomial fitting, but the unphysical wiggles persisten. 
- However, I discovered an strange behavior. If I change the fluent expressions (even without any meaningful change like adding zero or multiplying by one) and run the simulation again the residuals are resetted and the simulation can converge to a different solution.
- I suspect that this strange behavior is because changing the expression triggers a recompilation of the expression code and perhaps clearing some internal variables in Fluent (but I have not idea about what is really going on).
- Interestingly, when this happens the solver converges to a new solution with less wiggles.
- I repeated this process a couple of times until I obtained a solution where the unphysical wiggles dissapeared
- The take-away messages would be:
  - The wiggles are numerical artifacts most likely caused by the slope discontinuity
  - The approach I used to update the expressions is not practical, but it shows that Fluent can converge to a solution without wiggles

