import os
import datetime
import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial import Polynomial

from scipy.integrate import solve_ivp
from scipy.integrate._ivp.ivp import METHODS as ODE_METHODS

from . import fluid_properties as props
from . import pysolver_view as psv
from .plot_options import COLORS_MATLAB

PROCESS_TYPES = ["polytropic", "adiabatic"]
PROPERTY_CALCULATION_TYPES = ["blending", "equilibrium", "metastable"]

from functools import wraps


def ensure_data_available(method):
    """
    Decorator to ensure that poly_handles and poly_breakpoints are available before executing a method.
    """

    @wraps(method)
    def wrapper(self, *args, **kwargs):
        if not self.poly_handles or not self.poly_breakpoints:
            raise ValueError(
                "Polynomials or breakpoints are not initialized. Please fit the model first."
            )
        return method(self, *args, **kwargs)

    return wrapper


def barotropic_model_one_fluid(
    fluid_name,
    T_in,
    p_in,
    p_out,
    efficiency,
    calculation_type=None,
    blending_onset=None,
    blending_width=None,
    ODE_solver="lsoda",
    ODE_tolerance=1e-8,
    HEOS_solver="hybr",
    HEOS_tolerance=1e-6,
    HEOS_max_iter=100,
    HEOS_print_convergence=False,
    n_eval=None,
    polynomial_vars=[
        "density",
        "viscosity",
        "speed_sound",
        "void_fraction",
        "vapor_quality",
    ],
    polynomial_order=10,
):
    r"""
    Simulates a polytropic process for a given fluid from an inlet state to a specified outlet pressure.

    This function integrates the differential equation describing the polytropic process, which may involve phase changes, using the relationship:

    .. math::

        \frac{\text{d}h}{\text{d}p} = \frac{\eta_p}{\rho}

    where :math:`h` is the specific enthalpy, :math:`p` is the pressure, :math:`\eta_p` is the process efficiency, and :math:`\rho` is the density of the fluid.

    The evaluation of fluid properties varies based on the ``calculation_type`` parameter:

    - ``calculation_type='equilibrium'``: Computes properties assuming thermodynamic equilibrium using the CoolProp property solver.
    - ``calculation_type='metastable'``: Computes properties based on the Helmholtz energy equations of state, iterating on density and temperature.
    - ``calculation_type='blending'``: Blends equilibrium and metastable properties according to the fluid vapor quality.

    In single-phase regions, all calculation types yield identical results. The differences arise in the evaluation within the two-phase region, where the handling of equilibrium and metastable fluid properties differs.

    .. note::

        The choice of ODE solver can significantly impact the accuracy and stability of the simulation. The ``BDF`` and ``LSODA`` solvers are recommended for general use, particularly for stiff problems involving equilibrium phase change or narrow blending widths, while `RK45` may be suitable for non-stiff problems with no phase change or smooth phase transitions.

    Parameters
    ----------
    fluid_name : str
        The name of the fluid to be used in the simulation (e.g., 'Water').
    T_in : float
        Inlet temperature of the fluid in Kelvin.
    p_in : float
        Inlet pressure of the fluid in Pascals.
    p_out : float
        Outlet pressure of the fluid in Pascals.
    efficiency : float
        The efficiency of the polytropic process, dimensionless.
    calculation_type : str, optional
        The type of calculation to perform. Options include:
        - 'equilibrium': Computes equilibrium properties only.
        - 'metastable': Computes metastable properties only.
        - 'blending': Computes both equilibrium and metastable properties, blends them, and returns the blended properties.
        Default is None, which will raise an error.
    blending_onset : float, optional
        The onset of blending in the process, typically a value between 0 and 1. Required if `calculation_type` is 'blending'.
    blending_width : float, optional
        The width of the blending region, typically a value between 0 and 1. Required if `calculation_type` is 'blending'.
    ODE_solver : str, optional
        The solver to use for the ODE integration. Valid options include:
        - 'RK23', 'RK45', 'DOP853', 'Radau', 'BDF', 'LSODA'.
        Default is 'LSODA'. Recommended solvers are 'BDF', 'LSODA', or 'Radau' for stiff problems or 'RK45' for non-stiff problems with smooth blending.
    ODE_tolerance : float, optional
        The relative and absolute tolerance for the ODE solver. Default is 1e-8.
    HEOS_solver : str, optional
        The solver algorithm for the metastable state calculations. Default is 'hybr'.
    HEOS_tolerance : float, optional
        The tolerance for the HEOS solver. Default is 1e-6.
    HEOS_max_iter : int, optional
        The maximum number of iterations for the HEOS solver. Default is 100.
    HEOS_print_convergence : bool, optional
        If True, prints convergence information for the HEOS solver. Default is False.
    n_eval : int, optional
        The number of evaluation points for the solution. Default is None, which means the solver chooses the points.

    Returns
    -------
    states : list
        A list of states (dictionaries) representing the properties of the fluid at each evaluation point.
    solution : scipy.integrate.OdeResult
        The result of the ODE integration containing information about the solver process.

    Examples
    --------
    >>> states, solution = barotropic_model_one_fluid(
    ...     fluid_name='Water',
    ...     T_in=300.0,
    ...     p_in=101325.0,
    ...     p_out=200000.0,
    ...     efficiency=0.85,
    ...     calculation_type='blending',
    ...     blending_onset=0.1,
    ...     blending_width=0.05,
    ...     ODE_solver='LSODA',
    ...     ODE_tolerance=1e-8,
    ...     HEOS_solver='hybr',
    ...     HEOS_tolerance=1e-6,
    ...     HEOS_max_iter=100,
    ...     HEOS_print_convergence=False,
    ...     n_eval=100
    ... )
    """

    # Initialize fluid and compute inlet state
    fluid = props.Fluid(name=fluid_name, backend="HEOS", exceptions=True)
    state_in = fluid.get_state(
        props.PT_INPUTS, p_in, T_in, supersaturation=True, generalize_quality=True
    )

    # Determine the type of phase change process
    phase_change_type = (
        "flashing" if state_in.s < fluid.critical_point.s else "condensation"
    )

    # Initial guess for the metastable state (could be updated later)
    rhoT_guess_metastable = [state_in.rho, state_in.T]

    # Differential equation defining the polytropic process
    def odefun(t, y):
        nonlocal rhoT_guess_metastable  # Allows modification within the function scope

        # Rename input variables
        p = t
        (h,) = y

        # Do calculations according to specified type
        if calculation_type == "equilibrium":
            # Compute equilibrium state using CoolProp solver
            state_eq = fluid.get_state(
                input_type=props.HmassP_INPUTS,
                prop_1=h,
                prop_2=p,
                generalize_quality=True,
                supersaturation=True,
            )
            dhdp = np.atleast_1d(efficiency / state_eq["rho"])
            return dhdp, state_eq

        elif calculation_type == "metastable":
            # Compute metastable state using custom solver
            state_meta = fluid.get_state_metastable(
                prop_1="h",
                prop_1_value=h,
                prop_2="p",
                prop_2_value=p,
                rhoT_guess=rhoT_guess_metastable,
                generalize_quality=False,
                supersaturation=True,
                solver_algorithm=HEOS_solver,
                solver_tolerance=HEOS_tolerance,
                solver_max_iterations=HEOS_max_iter,
                print_convergence=HEOS_print_convergence,
            )

            # Update guess for the next integration step
            rhoT_guess_metastable = [state_meta.rho, state_meta.T]
            dhdp = np.atleast_1d(efficiency / state_meta["rho"])
            return dhdp, state_meta

        elif calculation_type == "blending":
            # Compute equilibrium state using CoolProp solver
            state_eq = fluid.get_state(
                input_type=props.HmassP_INPUTS,
                prop_1=h,
                prop_2=p,
                generalize_quality=True,
                supersaturation=True,
            )

            # Compute metastable state using custom solver
            state_meta = fluid.get_state_metastable(
                prop_1="h",
                prop_1_value=h,
                prop_2="p",
                prop_2_value=p,
                rhoT_guess=rhoT_guess_metastable,
                generalize_quality=False,
                supersaturation=True,
                solver_algorithm=HEOS_solver,
                solver_tolerance=HEOS_tolerance,
                solver_max_iterations=HEOS_max_iter,
                print_convergence=HEOS_print_convergence,
            )

            # Update guess for the next integration step
            rhoT_guess_metastable = [state_meta.rho, state_meta.T]

            # Compute blended state
            state_blend = props.blend_properties(
                phase_change=phase_change_type,
                props_equilibrium=state_eq,
                props_metastable=state_meta,
                blending_variable="Q",
                blending_onset=blending_onset,
                blending_width=blending_width,
            )

            # Compute right-hand-side of the ODE
            dhdp = np.atleast_1d(efficiency / state_blend["rho"])
            return dhdp, props.FluidState(state_blend, fluid.name)

        else:
            raise ValueError(
                f"Invalid calculation_type='{calculation_type}'. "
                f"Valid options are: {', '.join(PROPERTY_CALCULATION_TYPES)}."
            )

    # Check if the provided ODE_solver is valid
    valid_solvers_heos = psv.nonlinear_system.SOLVER_OPTIONS
    if HEOS_solver not in valid_solvers_heos:
        error_message = (
            f"Invalid HEOS solver '{HEOS_solver}' provided. "
            f"Valid solvers are: {', '.join(valid_solvers_heos)}. "
        )
        raise ValueError(error_message)

    # Check if the provided ODE_solver is valid
    valid_solvers_ode = list(ODE_METHODS.keys())
    if ODE_solver not in valid_solvers_ode:
        error_message = (
            f"Invalid ODE solver '{ODE_solver}' provided. "
            f"Valid solver are: {', '.join(valid_solvers_ode)}. "
            "Recommended solvers: 'BDF', 'LSODA' or 'Radau' stiff problems involving equilibrium property calculations or blending calculations with a narrow blending width. 'RK45' can be used for non-stiff problems with a wide blending width."
        )
        raise ValueError(error_message)

    # Solve polytropic expansion differential equation
    ode_sol = solve_ivp(
        lambda p, h: odefun(p, h)[0],  # Get only first output
        [p_in, p_out],
        [state_in.h],
        t_eval=np.linspace(p_in, p_out, n_eval) if n_eval else None,
        method=ODE_solver,
        rtol=ODE_tolerance,
        atol=ODE_tolerance,
    )
    if not ode_sol.success:
        raise Exception(ode_sol.message)

    # Postprocess solution
    states = postprocess_ode(ode_sol.t, ode_sol.y, odefun)

    # Fit polynomials to data
    polyfit = PolynomialFitter(
        fluid,
        states,
        polynomial_vars,
        polynomial_order,
        calculation_type=calculation_type,
    )

    return polyfit


def postprocess_ode(t, y, ode_handle):
    """
    Post-processes the output of an ordinary differential equation (ODE) solver.

    This function takes the time points and corresponding ODE solution matrix,
    and for each time point, it calls a user-defined ODE handling function to
    process the state of the ODE system. It collects the results into a
    dictionary where each key corresponds to a state variable and the values
    are numpy arrays of that state variable at each integration step

    Parameters
    ----------
    t : array_like
        Integration points at which the ODE was solved, as a 1D numpy array.
    y : array_like
        The solution of the ODE system, as a 2D numpy array with shape (n,m) where
        n is the number of points and m is the number of state variables.
    ode_handle : callable
        A function that takes in a integration point and state vector and returns a tuple,
        where the first element is ignored (can be None) and the second element
        is a dictionary representing the processed state of the system.

    Returns
    -------
    ode_out : dict
        A dictionary where each key corresponds to a state variable and each value
        is a numpy array containing the values of that state variable at each integration step.
    """
    # Initialize ode_out as a dictionary
    ode_out = {}
    for t_i, y_i in zip(t, y.T):
        _, out = ode_handle(t_i, y_i)

        for key, value in out.items():
            # Initialize with an empty list
            if key not in ode_out:
                ode_out[key] = []
            # Append the value to list of current key
            ode_out[key].append(value)

    # Convert lists to numpy arrays
    for key in ode_out:
        ode_out[key] = np.array(ode_out[key])

    return ode_out


class PolynomialFitter:

    def __init__(self, fluid, states, variables, order, calculation_type):
        self.fluid = fluid
        self.states = states
        self.variables = variables
        self.order = order
        self.calculation_type = calculation_type
        self.poly_handles = {}
        self.poly_breakpoints = None

        self.variable_plot_labels = {
            "density": "Density (kg/m$^3$)",
            "viscosity": "Viscosity (Pa·s)",
            "speed_sound": "Speed of sound (m/s)",
            "void_fraction": "Void fraction",
            "vapor_quality": "Vapor quality",
        }

        if self.calculation_type == "equilibrium":
            self.fit_equilibrium()
        elif self.calculation_type == "blending":
            self.fit_blending()
        else:
            raise ValueError("Unsupported calculation type")

    def fit_equilibrium(self):

        # Scale pressure by inlet pressure to improve polynomial conditioning
        p = self.states["p"] / self.states["p"][0]

        # Determine points within the two-phase region
        eps = 1e-9
        mask_1phase = np.abs(self.states["supersaturation_degree"]) > eps
        mask_2phase = np.abs(self.states["supersaturation_degree"]) <= eps

        # Determine phase change pressure
        index = np.where(mask_2phase)[0][0] if mask_2phase.any() else len(p)
        self.poly_breakpoints = [p[index]] if mask_2phase.any() else [p[-1]]

        # Fit polynomials to data
        for var in self.variables:

            self.poly_handles[var] = []

            if mask_1phase.any():  # Single-phase data
                y_1p = np.array(self.states[var])[mask_1phase]
                p_1p = p[mask_1phase]
                poly_1p = Polynomial.fit(p_1p, y_1p, deg=self.order).convert()
                self.poly_handles[var].append(poly_1p)

            if mask_2phase.any():  # Two-phase data
                y_2p = np.array(self.states[var])[mask_2phase]
                p_2p = p[mask_2phase]
                poly_2p = Polynomial.fit(p_2p, y_2p, deg=self.order).convert()

                # Ensure continuity by adjusting the first coefficient
                y_1p_breakpoint = poly_1p(self.poly_breakpoints[0])
                y_2p_breakpoint = poly_2p(self.poly_breakpoints[0])
                coeffs_2p = poly_2p.coef
                coeffs_2p[0] -= y_2p_breakpoint - y_1p_breakpoint
                poly_2p = Polynomial(coeffs_2p)
                self.poly_handles[var].append(poly_2p)

    def fit_blending(self):
        # Get the scaled pressure values
        p = self.states["p"] / self.states["p"][0]

        # Define masks for different regions based on states["x"]
        mask_region_1 = self.states["x"] > 1
        mask_region_2 = (self.states["x"] >= 0) & (self.states["x"] <= 1)
        mask_region_3 = self.states["x"] < 0

        # Determine breakpoints in terms of scaled pressure
        p_region_2_start = p[mask_region_2][0]
        p_region_2_end = p[mask_region_2][-1]
        self.poly_breakpoints = [p_region_2_start, p_region_2_end]

        # Fit polynomials to data
        for var in self.variables:
            self.poly_handles[var] = []

            # Fit polynomial for Region 1 (x > 1)
            y_1 = np.array(self.states[var])[mask_region_1]
            p_1 = p[mask_region_1]
            poly_1 = Polynomial.fit(p_1, y_1, deg=self.order).convert()
            self.poly_handles[var].append(poly_1)

            # Fit polynomial for Region 2 (0 <= x <= 1)
            y_2 = np.array(self.states[var])[mask_region_2]
            p_2 = p[mask_region_2]
            poly_2 = Polynomial.fit(p_2, y_2, deg=self.order).convert()

            # Ensure continuity by adjusting the first coefficient
            y_1_end = poly_1(p_region_2_start)
            y_2_start = poly_2(p_region_2_start)
            coeffs_2 = poly_2.coef
            coeffs_2[0] -= y_2_start - y_1_end
            poly_2 = Polynomial(coeffs_2)
            self.poly_handles[var].append(poly_2)

            # Fit polynomial for Region 3 (x < 0)
            y_3 = np.array(self.states[var])[mask_region_3]
            p_3 = p[mask_region_3]
            poly_3 = Polynomial.fit(p_3, y_3, deg=self.order).convert()

            # Ensure continuity by adjusting the first coefficient
            y_2_end = poly_2(p_region_2_end)
            y_3_start = poly_3(p_region_2_end)
            coeffs_3 = poly_3.coef
            coeffs_3[0] -= y_3_start - y_2_end
            poly_3 = Polynomial(coeffs_3)
            self.poly_handles[var].append(poly_3)

    @ensure_data_available
    def evaluate_polynomials(self, p_scaled, variable):
        """
        Evaluates the polynomials at a given pressure with respect to the breakpoints.
        """
        # Initialize an array to store the evaluated values for the specified variable
        p_scaled = np.atleast_1d(p_scaled)
        evaluated_values = np.zeros_like(p_scaled, dtype=float)

        # Ensure poly_breakpoints is array (decreasing pressure)
        breakpoints = np.atleast_1d(self.poly_breakpoints)
        breakpoints = np.concatenate(([+np.inf], breakpoints, [-np.inf]))

        # Handle multiple breakpoints and corresponding polynomials
        if variable in self.poly_handles:
            for i, poly in enumerate(self.poly_handles.get(variable, [])):
                mask = (p_scaled < breakpoints[i]) & (p_scaled >= breakpoints[i + 1])
                if np.any(mask):  # Check if there are any True values in the mask
                    evaluated_values[mask] = poly(p_scaled[mask])
        else:
            raise ValueError(f"No polynomials found for variable '{variable}'")

        return evaluated_values

    @ensure_data_available
    def plot_polynomial(self, var, n_points=100):
        """
        Plots the evaluated polynomials for the specified variable from start to end pressure,
        including additional markers for the breakpoints

        Parameters
        ----------
        var : str
            The variable for which to plot the polynomial.
        n_points : int, optional
            The number of points to generate between breakpoints for plotting. Default is 100.
        """

        # Create a figure with two vertically stacked subplots
        fig, ax = plt.subplots(figsize=(6, 5))

        # Mapping variable names to axis labels
        ax.set_xlabel("Pressure (Pa)")
        ax.set_ylabel(self.variable_plot_labels.get(var, var))

        # Generate pressure values for the current segment
        p_eval = self._generate_scaled_pressure_intervals(n_points)
        prop_eval = self.evaluate_polynomials(p_eval, var)
        ax.plot(
            p_eval * self.states["p"][0],  # Original pressure values
            prop_eval,
            color=COLORS_MATLAB[0],
        )

        # Plot points just before and after the breakpoints to check continuity
        eps = 1e-12
        for bp in self.poly_breakpoints:

            # Just before the breakpoint
            x_minus = bp - eps
            y_minus = self.evaluate_polynomials(x_minus, var)
            ax.plot(
                [x_minus * self.states["p"][0]],
                [y_minus],
                marker="o",
                color=COLORS_MATLAB[0],
            )

            # Just after the breakpoint
            x_plus = bp + eps
            y_plus = self.evaluate_polynomials(x_plus, var)
            ax.plot(
                [x_plus * self.states["p"][0]],
                [y_plus],
                marker="o",
                color=COLORS_MATLAB[0],
            )

        return fig

    @ensure_data_available
    def plot_polynomial_and_error(self, var):
        """
        Plots the variable from polynomials and states, and the relative error between them.

        Parameters
        ----------
        var : str
            The variable to plot (e.g., 'density', 'viscosity', etc.).
        """

        # Create a figure with two vertically stacked subplots
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(6, 5))

        # Mapping variable names to axis labels
        x_label = "Pressure (Pa)"
        y_label = self.variable_plot_labels.get(
            var, var
        )  # Default to variable name if not found

        # Generate pressure values and evaluate the variable
        p_eval_original = self.states["p"]
        p_eval_scaled = p_eval_original / self.states["p"][0]
        poly_eval = self.evaluate_polynomials(p_eval_scaled, var)
        state_values = self.states[var]

        # Compute the relative error
        relative_error = 100 * np.abs(poly_eval - state_values) / np.abs(state_values)

        # Plot the variable from polynomials and states in the upper subplot
        ax1.plot(
            p_eval_original,
            poly_eval,
            label=f"Polynomial fit",
            linestyle="-",
            color=COLORS_MATLAB[1],
        )
        ax1.plot(
            p_eval_original,
            state_values,
            label=f"Original values",
            marker="o",
            markersize=2.5,
            linestyle="none",
            color=COLORS_MATLAB[0],
        )
        ax1.set_ylabel(y_label)
        ax1.legend(loc="lower right")

        # Plot the relative error in the lower subplot
        ax2.plot(
            p_eval_original,
            relative_error,
            label="Relative Error",
            linestyle="-",
            markersize=2.5,
            color=COLORS_MATLAB[1],
        )
        ax2.set_xlabel(x_label)
        ax2.set_ylabel("Relative error (%)")
        plt.tight_layout(pad=1)

        return fig

    @ensure_data_available
    def _generate_scaled_pressure_intervals(self, n_points=100):
        """
        Generate scaled pressure intervals between breakpoints with a specified number of points per interval.

        Parameters
        ----------
        n_points : int, optional
            The number of points to generate per interval. Default is 100.

        Returns
        -------
        numpy.ndarray
            A numpy array of scaled pressure values covering all intervals, concatenated into a single array.
        """

        # Ensure breakpoints are array-like and create a list with start, breakpoints, and end
        p_start = 1.00
        p_end = self.states["p"][-1] / self.states["p"][0]
        all_points = [p_start] + list(np.atleast_1d(self.poly_breakpoints)) + [p_end]

        # Iterate over the intervals defined by start, breakpoints, and end
        p_scaled_intervals = []
        for i in range(len(all_points) - 1):
            p_eval = np.linspace(all_points[i], all_points[i + 1], n_points)
            p_scaled = p_eval / p_start  # Scale with respect to the initial pressure
            p_scaled_intervals.append(p_scaled)

        return np.concatenate(p_scaled_intervals)

    def export_fluent_expressions(self, case_name=None, output_dir=None):
        """
        Exports the polynomial expressions for use in Fluent CFD software.

        Parameters
        ----------
        fluid_name : str, optional
            The name of the fluid. Default is None, which uses the class's fluid name.
        case_name : str, optional
            The case name for the export file. Default is based on the fluid name.
        output_dir : str, optional
            The directory where the export file will be saved. Default is current directory.
        """
        fluid_name = self.fluid.name
        case_name = case_name or f"barotropic_model_{fluid_name}"
        output_dir = output_dir or f"barotropic_model_{fluid_name}"

        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        file_path = os.path.join(output_dir, 'fluent_expressions.txt')
        with open(file_path, 'wt') as fid:
            # File header
            fid.write(f'Fluent expressions of {fluid_name} properties along polytropic process\n')
            fid.write(f'Creation datetime: {datetime.datetime.now()}\n\n')

            # Define the units of the fluid properties
            units = {
                'density': "kg/m^3",
                'viscosity': "Pa·s",
                'speed_of_sound': "m/s",
                'void_fraction': "none",
                'vapor_quality': "none",
            }

            # Generate expressions for each property
            for var in self.variables:
                if var in units:
                    unit = units[var]
                    expressions = self.generate_piecewise_expression(var, unit, 'horner')
                    fid.write(f'{case_name}_{var}\n')
                    fid.write(f'{expressions}\n\n')

    def generate_piecewise_expression(self, var, unit, polynomial_form):
        """
        Generates a piecewise polynomial expression for a given variable.

        Parameters
        ----------
        var : str
            The variable name for which to generate the expression.
        unit : str
            The unit of the variable.
        polynomial_form : str
            The form of the polynomial ('horner' or 'standard').

        Returns
        -------
        str
            The piecewise polynomial expression.
        """
        expression = self.polynomial_expression(self.poly_handles[var][0].coef, polynomial_form)
        for i in range(1, len(self.poly_handles[var])):
            new_expression = self.polynomial_expression(self.poly_handles[var][i].coef, polynomial_form)
            expression = self.if_expression(expression, new_expression, self.poly_breakpoints[i - 1])

        # Add fluid property unit
        if unit != "none":
            expression = f'{expression} * 1 [{unit}]'
        return expression

    @staticmethod
    def polynomial_expression(coefficients, polynomial_form):
        """
        Generates a polynomial expression string.

        Parameters
        ----------
        coefficients : array-like
            Coefficients of the polynomial.
        polynomial_form : str
            The form of the polynomial ('horner' or 'standard').

        Returns
        -------
        str
            The polynomial expression string.
        """
        if polynomial_form == "horner":
            polynomial_string = f'{coefficients[0]:+0.6e}'
            for i in range(1, len(coefficients)):
                polynomial_string = f'{coefficients[i]:+0.6e} + (AbsolutePressure/{1.0} [Pa]) *\n({polynomial_string})'
        elif polynomial_form == "standard":
            terms = [f'{coefficients[i]:+0.6e} * (AbsolutePressure/{1.0} [Pa])^{i}' for i in range(len(coefficients))]
            polynomial_string = ' \n'.join(terms)
        else:
            raise ValueError("The polynomial form must be 'horner' or 'standard'")
        return polynomial_string

    @staticmethod
    def if_expression(expression_1, expression_2, transition_pressure):
        """
        Generates an IF expression string for piecewise functions.

        Parameters
        ----------
        expression_1 : str
            The expression for the condition being true.
        expression_2 : str
            The expression for the condition being false.
        transition_pressure : float
            The transition pressure where the expressions change.

        Returns
        -------
        str
            The IF expression string.
        """
        return f'IF(AbsolutePressure>={transition_pressure} [Pa], \n{expression_1}, \n{expression_2})'
