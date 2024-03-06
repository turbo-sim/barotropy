import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
import copy

import barotropy as bpy
import barotropy.fluids.core as props
import barotropy.sCO2_utilities as sco2

from barotropy.PySolverView import NonlinearSystemSolver, OptimizationSolver


bpy.set_plot_options(grid=False)


fluid = props.Fluid("CO2", backend="HEOS")

T = 1.01*fluid.critical_point.T
s = 0.95*fluid.critical_point.s


# Solve as root finder
state = fluid.set_state(props.SmassT_INPUTS, s, T)
sonic_problem = props.SonicStateProblem(fluid, props.SmassT_INPUTS, s, T)
solver = NonlinearSystemSolver(sonic_problem, sonic_problem.initial_guess, display=True, plot=True)
solution = solver.solve(method='hybr')
density = solution.x[0]
temperature = solution.x[1]
sonicState = fluid.set_state(props.DmassT_INPUTS, density, temperature)

# Solve as optimization
state = fluid.set_state(props.SmassT_INPUTS, s, T)
sonic_problem = props.SonicStateProblem2(fluid, props.SmassT_INPUTS, s, T)
solver = OptimizationSolver(sonic_problem, sonic_problem.initial_guess, display=True, plot=True)
solution = solver.solve(method="slsqp")
density = solution.x[0]
temperature = solution.x[1]
sonicState2 = fluid.set_state(props.DmassT_INPUTS, density, temperature)



# Create figure
fig, ax = plt.subplots(figsize=(6.8, 4.8))
ax.set_xlabel("Entropy [J/kg/K]")
ax.set_ylabel("Temperature [K]")

props.plot_phase_diagram(
    "s",
    "T",
    fluid,
    axes=ax,
    plot_critical_point=True,
    plot_quality_isolines=True,
    plot_saturation_line=True,
    plot_pseudocritical_line=False,
)

ax.plot(state.s, state.T,
        marker='o',
        markersize=4.5,
        markerfacecolor="w",
        linewidth=0.75,
        linestyle="none",
        markeredgewidth=1.00,
        # color='b',
        # label=label,
    )

ax.plot(sonicState.s, sonicState.T,
        marker='^',
        markersize=4.5,
        markerfacecolor="w",
        linewidth=0.75,
        linestyle="none",
        markeredgewidth=1.00,
        # color='b',
        # label=label,
    )

ax.plot(sonicState2.s, sonicState2.T,
        marker='v',
        markersize=4.5,
        markerfacecolor="w",
        linewidth=0.75,
        linestyle="none",
        markeredgewidth=1.00,
        # color='b',
        # label=label,
    )


plt.show()
