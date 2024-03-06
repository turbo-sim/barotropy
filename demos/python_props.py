import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
import copy

from barotropy import Fluid, NonlinearSystemProblem, PropertyRoot, NonlinearSystemSolver, set_plot_options

set_plot_options(grid=False)



fluid = Fluid('Water', backend="HEOS")

# # 
# props = compute_properties_metastable_rhoT(10, 500, fluid.abstractstate)
# print("Metastable properties of water")
# print(f"{'Property':35} {'value':6}")
# for key, value in props.items():
#     print(f"{key:35} {value:.6e}")

# Check that the metastable property calculations match in the single-phase region
p, T =  101325, 400
props_stable = fluid.set_state(CP.PT_INPUTS, p, T)
print()
print(f"Properties of water at p={p} Pa and T={T} K")
print(f"{'Property':35} {'Equilibrium':>15} {'Metastable':>15} {'Deviation':>15}")
props_metastable = fluid.set_state_metastable_rhoT(props_stable["rho"], props_stable["T"])
for key in props_stable.keys():
    value_stable = props_stable[key]
    value_metastable = props_metastable[key]
    print(f"{key:35} {value_stable:+15.6e} {value_metastable:+15.6e} {(value_stable - value_metastable)/value_stable:+15.6e}")

    
# # print(CP.DmassT_INPUTS)
# props_metastable = fluid.set_state_metastable



p = props_stable["p"]
h = props_stable["h"]
problem = PropertyRoot("p", p, "h", h, fluid)

print(p, h)

x0 = [props_stable["rho"]*1.2, props_stable["T"]+50]
res = problem.get_values(x0)

print(res)


solver = NonlinearSystemSolver(problem, x0, display=True, plot=True)
solution = solver.solve(method="lm", x0=x0)

print(props_stable["rho"], props_stable["T"])


plt.show()


# TODO clearn metastable state calculation example for general inputs
# TODO make test calculations of metastable with different inputs p,h, T,s, h,s
# TODO make the tests from different initial points not very far away from the actual solution
# TODO clean nonlinear system and optimizer plotting:
    # - Plot each function call that is not comming from the gradient evaluation (line searches)
    # - The hybrid algorithm fro equation solving does line serches in the same iteration, important to display these points, or things do not make sense
    # - Important to have ticks only on the integer values for the plotting of the solver
    # - Add option to get_values(print_progress=True) and set explicitly to false when the evaluation is from the get_jacobian finite diffeence function
    # - Add counter for F-count and gradF-count in both optimization and equation solver. Option do not print on gradient finite difference