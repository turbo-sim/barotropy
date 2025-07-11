import os
import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pysolver_view as psv
import coolpropx as props

from functools import wraps
from numpy.polynomial import Polynomial
from scipy.integrate import solve_ivp
from scipy.integrate._ivp.ivp import METHODS as ODE_METHODS

from . import graphics
from . import utilities as utils

COLORS_MATLAB = graphics.COLORS_MATLAB
PROCESS_TYPES = ["polytropic", "adiabatic"]
CALCULATION_TYPES = ["blending", "equilibrium", "metastable"]



class BarotropicModel:
    """
    Coordinates the simulation and polynomial fitting for a barotropic process.

    Parameters
    ----------
    fluid_name : str or list of str
        The name(s) of the fluid(s) for the barotropic model.

        - Specify a single string (e.g., 'co2') to use the single-component model.
        - Specify a list of two strings (e.g., ['water', 'nitrogen']) to use the two-component model.

    T_in : float
        Inlet temperature of the fluid in Kelvin.

    p_in : float
        Inlet pressure of the fluid in Pascals.

    rho_in : float
        Inlet density of the fluid in kilogram per cubic meter.

        .. note::
            Applicable only to the one-component model.

    p_min: float
        Minimum pressure for the barotropic process integration in Pascals.

    p_max : float
        Maximum pressure for the barotropic process integration in Pascals.

    efficiency : float
        The efficiency of the polytropic process, dimensionless.

    mixture_ratio : float
        Mass ratio of the first to the second fluid in the mixture.

        .. note::
            Applicable only to the two-component model.

    process_type : str, optional
        The type of polytropic process that the fluid experiences. Must be 'expansion' or 'compression'.

    calculation_type : str, optional
        The type of fluid property calculation for the one-component model. Options include:

        - ``equilibrium``: Computes equilibrium properties only.
        - ``metastable``: Computes metastable properties only.
        - ``blending``: Computes both equilibrium and metastable properties, blends them, and returns the blended properties.

        .. note::
            Applicable only to the one-component model.

    blending_onset : float, optional
        The onset of blending in the process, typically a value between 0 and 1. Required when `calculation_type` is ``blending``.

        .. note::
            Applicable only to the one-component model.

    blending_width : float, optional
        The width of the blending region, typically a value between 0 and 1. Required when `calculation_type` is ``blending``.

        .. note::
            Applicable only to the one-component model.

    HEOS_solver : str, optional
        The solver algorithm used to compute the metastable states. Valid options:

        .. list-table::
            :widths: 20 50
            :header-rows: 1

            * - Solver name
              - Description
            * - ``hybr``
              -  Powell's hybrid trust-region algorithm
            * - ``lm``
              - Levenberg–Marquardt algorithm

        See `Scipy root() <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root.html>`_ for more info.
        Recommended solvers: Both ``hybr`` and ``lm`` work well for the tested cases.

        .. note::
            Applicable only to the one-component model.

    HEOS_tolerance : float, optional
        The tolerance for the HEOS solver.

        .. note::
            Applicable only to the one-component model.

    HEOS_max_iter : int, optional
        The maximum number of iterations for the HEOS solver.

        .. note::
            Applicable only to the one-component model.

    HEOS_print_convergence : bool, optional
        If True, prints convergence information for the HEOS solver.

        .. note::
            Applicable only to the one-component model.

    ODE_solver : str, optional
        The solver to use for the ODE integration. Valid options:

        .. list-table::
            :widths: 20 50
            :header-rows: 1

            * - Solver name
              - Description
            * - ``RK23``
              - Explicit Runge-Kutta method of order 3(2)
            * - ``RK45``
              - Explicit Runge-Kutta method of order 5(4)
            * - ``DOP853``
              - Explicit Runge-Kutta method of order 8
            * - ``Radau``
              - Implicit Runge-Kutta method of the Radau IIA family of order 5
            * - ``BDF``
              - Implicit multi-step variable-order (1 to 5) method based on a backward differentiation formula for the derivative approximation
            * - ``LSODA``
              - Adams/BDF method with automatic stiffness detection and switching

        See `Scipy solver_ivp() <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html>`_  for more info.
        Recommended solvers: ``BDF``, ``LSODA``, or ``Radau`` for stiff problems or ``RK45`` for non-stiff problems with smooth blending.

    ODE_tolerance : float, optional
        The relative and absolute tolerance for the ODE solver.

    polynomial_degree : int
        Degree of the polynomials to fit.

        .. note::

            When `calculation_type` is ``blending`` the degree of the polynomial in the blending region is set to 4 to achieve sufficient accuracy while preventing numerical round-off errors associated with single-precision arithmetic in CFD solvers.

    polynomial_format : str, optional
        Type of polynomial representation (``horner`` or ``standard``).

    polynomial_variables : list of str
        A list of variable names to fit polynomials to, such as 'density', 'viscosity', 'speed_sound', 'void_fraction', 'vapor_quality'.

    output_dir : str
        The directory where output will be saved.

    """

    def __init__(
        self,
        fluid_name,
        p_in: float,
        T_in: float = None,
        rho_in: float = None,
        p_out: float = None,
        p_min: float = None,
        p_max: float = None,
        efficiency: float = 1.00,
        mixture_ratio: float = None,
        process_type: str = None,
        calculation_type: str = None,
        blending_onset: float = None,
        blending_width: float = None,
        HEOS_solver: str = "hybr",
        HEOS_tolerance: float = 1e-6,
        HEOS_max_iter: int = 100,
        HEOS_print_convergence: bool = False,
        ODE_solver: str = "LSODA",
        ODE_tolerance: float = 1e-8,
        polynomial_degree: int = 8,
        polynomial_format: str = "horner",
        polynomial_variables: list = [
            "density",
            "viscosity",
            "speed_sound",
            "void_fraction",
            "vapor_quality",
        ],
        output_dir: str = "barotropic_model",
    ):

        # Validate inputs and find the correct model type
        self.model_type = self._resolve_model_type(
            fluid_name,
            mixture_ratio,
            calculation_type,
            blending_onset,
            blending_width,
        )

        # Define inlet state
        self.fluid_name = fluid_name

        if self.model_type == "one-component":
            state_in = self._resolve_inlet_state(fluid_name, p_in, T_in, rho_in)
            self.p_in = state_in.p
            self.T_in = state_in.T
            self.rho_in = state_in.rho
        else:
            self.p_in = p_in
            self.T_in = T_in
            self.rho_in = None  # not used for two-component models

        # Define pressure range
        self.p_min, self.p_max, self.process_type = self._resolve_pressure_range(
            p_in=self.p_in,
            p_out=p_out,
            p_min=p_min,
            p_max=p_max,
            process_type=process_type,
        )

        # Assign optional arguments with defaults
        self.mixture_ratio = mixture_ratio
        self.efficiency = efficiency
        self.calculation_type = calculation_type
        self.blending_onset = blending_onset
        self.blending_width = blending_width
        self.ODE_solver = ODE_solver
        self.ODE_tolerance = ODE_tolerance
        self.HEOS_solver = HEOS_solver
        self.HEOS_tolerance = HEOS_tolerance
        self.HEOS_max_iter = HEOS_max_iter
        self.HEOS_print_convergence = HEOS_print_convergence
        self.polynomial_degree = polynomial_degree
        self.polynomial_format = polynomial_format
        self.polynomial_variables = polynomial_variables
        self.output_dir = output_dir

        # Initialize variables
        self.states = None
        self.ode_solution = None
        self.poly_fitter = None


    def _resolve_model_type(
        self,
        fluid_name,
        mixture_ratio,
        calculation_type,
        blending_onset,
        blending_width,
    ):
        # Validate inputs for one-component and two-component models
        if isinstance(fluid_name, str):
            # One-component model
            if calculation_type not in CALCULATION_TYPES:
                allowed_types = ", ".join(f"'{ptype}'" for ptype in CALCULATION_TYPES)
                msg = (
                    f"Invalid parameter calculation_type={calculation_type}. "
                    f"It must be one of {allowed_types} for the one-component model."
                )
                raise ValueError(msg)
            elif mixture_ratio is not None:
                msg = f"Parameter mixture_ratio={mixture_ratio} must be None for the one-component model."
                raise ValueError(msg)
            else:
                return "one-component"

        elif (
            isinstance(fluid_name, list)
            and len(fluid_name) == 2
            and all(isinstance(f, str) for f in fluid_name)
        ):
            # Two-component model
            if mixture_ratio is None or mixture_ratio <= 0:
                msg = f"Parameter mixture_ratio={mixture_ratio} must be a positive scalar for the two-component model."
                raise ValueError(msg)
            elif (
                calculation_type is not None
                or blending_onset is not None
                or blending_width is not None
            ):
                msg = (
                    "The following parameters must be None for the two-component model:\n"
                    f"  - calculation_type: {calculation_type}\n"
                    f"  - blending_onset: {blending_onset}\n"
                    f"  - blending_width: {blending_width}\n"
                )
                raise ValueError(msg)
            return "two-component"

        else:
            raise ValueError(
                f"Input parameter 'fluid_name'={fluid_name} must be either:\n"
                f"  1. A single string for the one-component barotropic model\n"
                f"  2. A list of two strings for the two-component barotropic model."
            )


    def _resolve_inlet_state(self, fluid_name, p_in=None, T_in=None, rho_in=None):
        """
        Resolves the inlet thermodynamic state from any valid combination of (p_in, T_in, rho_in).

        Returns:
            state_in
        Raises:
            ValueError if input combination is invalid
        """

        fluid = props.Fluid(name=fluid_name, backend="HEOS")
        num_given = sum(x is not None for x in [p_in, T_in, rho_in])
        if num_given < 2:
            raise ValueError("At least two of (p_in, T_in, rho_in) must be provided.")
        if num_given > 2:
            raise ValueError("Only two of (p_in, T_in, rho_in) should be provided. The third will be calculated.")

        # Determine which input pair was provided and compute the missing property
        if p_in is not None and T_in is not None:
            state = fluid.get_state(props.PT_INPUTS, p_in, T_in)

        elif p_in is not None and rho_in is not None:
            state = fluid.get_state(props.DmassP_INPUTS, rho_in, p_in)

        elif T_in is not None and rho_in is not None:
            state = fluid.get_state(props.DmassT_INPUTS, rho_in, T_in)

        else:
            raise ValueError("Unexpected error during state resolution.")

        return state


    def _resolve_pressure_range(self, p_in, p_out, p_min, p_max, process_type):
        """
        Resolve p_min, p_max, and process_type from user input.

        Args:
            p_in (float): Inlet pressure [Pa]
            p_out (float or None): Outlet pressure [Pa]
            p_min (float or None): Minimum pressure [Pa]
            p_max (float or None): Maximum pressure [Pa]
            process_type (str or None): "expansion" or "compression"

        Returns:
            Tuple (p_min, p_max, process_type)

        Raises:
            ValueError: If input combinations are invalid
        """

        # Case 1: process_type is explicitly given → p_min and p_max are mandatory
        if process_type is not None:
            if process_type not in ("expansion", "compression"):
                raise ValueError(
                    f"Invalid process_type='{process_type}'. Must be 'expansion' or 'compression'."
                )
            if p_min is None or p_max is None:
                raise ValueError(
                    f"When process_type is specified, both p_min and p_max must be provided. "
                    f"Got p_min={p_min}, p_max={p_max}."
                )
            if p_in < p_min or p_in > p_max:
                raise ValueError(
                    f"Inlet pressure p_in={p_in:.2e} must lie within p_min={p_min:.2e} and p_max={p_max:.2e}."
                )
            return p_min, p_max, process_type

        # Case 2: Infer everything from p_in and p_out
        if p_out is None:
            raise ValueError("If process_type is not specified, both p_in and p_out must be provided.")

        if p_min is not None or p_max is not None:
            raise ValueError(
                "When using p_in and p_out to infer the process type, both p_min and p_max must be None."
            )

        process_type = "expansion" if p_out < p_in else "compression"
        p_min = min(p_in, p_out)
        p_max = max(p_in, p_out)

        return p_min, p_max, process_type


    def solve(self):
        """
        Solves the equations for the one-component or two-component barotropic model and stores the fluid properties.
        The type of calculation performed is selected automatically depending on the number of fluid names defined when initializing the class.

        See Also
        --------
        barotropic_model_one_component :
            Calculation of fluid properties for the one-component model.
        barotropic_model_two_component :
            Calculation of fluid properties for the two-component model.
        """
        # Manually pass all parameters to the ODE solver function
        if self.model_type == "one-component":

            # Integrate upward from p_in to p_max:
            states_up, ode_solution_up = barotropic_model_one_component(
                fluid_name=self.fluid_name,
                p_in=self.p_in,
                rho_in=self.rho_in,
                p_out=self.p_max,
                efficiency=self.efficiency,
                process_type=self.process_type,
                calculation_type=self.calculation_type,
                blending_onset=self.blending_onset,
                blending_width=self.blending_width,
                ODE_solver=self.ODE_solver,
                ODE_tolerance=self.ODE_tolerance,
                HEOS_solver=self.HEOS_solver,
                HEOS_tolerance=self.HEOS_tolerance,
                HEOS_max_iter=self.HEOS_max_iter,
                HEOS_print_convergence=self.HEOS_print_convergence,
            )

            # Integrate downward from p_in to p_min:
            states_down, ode_solution_down = barotropic_model_one_component(
                fluid_name=self.fluid_name,
                p_in=self.p_in,
                rho_in=self.rho_in,
                p_out=self.p_min,
                efficiency=self.efficiency,
                process_type=self.process_type,
                calculation_type=self.calculation_type,
                blending_onset=self.blending_onset,
                blending_width=self.blending_width,
                ODE_solver=self.ODE_solver,
                ODE_tolerance=self.ODE_tolerance,
                HEOS_solver=self.HEOS_solver,
                HEOS_tolerance=self.HEOS_tolerance,
                HEOS_max_iter=self.HEOS_max_iter,
                HEOS_print_convergence=self.HEOS_print_convergence,
            )

            # Store solution
            self.state_in = {key: value[0] for key, value in states_up.items()}
            self.states = {
                key: np.concatenate((np.flipud(states_up[key][:-1]), states_down[key]))
                for key in states_up
            }
            self.ode_solution = (ode_solution_up, ode_solution_down)


        elif self.model_type == "two-component":
            states_up, ode_solution_up = barotropic_model_two_component(
                fluid_name_1=self.fluid_name[0],
                fluid_name_2=self.fluid_name[1],
                mixture_ratio=self.mixture_ratio,
                T_in=self.T_in,
                p_in=self.p_in,
                p_out=self.p_max,
                efficiency=self.efficiency,
                process_type=self.process_type,
                ODE_solver=self.ODE_solver,
                ODE_tolerance=self.ODE_tolerance,
            )

            states_down, ode_solution_down = barotropic_model_two_component(
                fluid_name_1=self.fluid_name[0],
                fluid_name_2=self.fluid_name[1],
                mixture_ratio=self.mixture_ratio,
                T_in=self.T_in,
                p_in=self.p_in,
                p_out=self.p_min,
                efficiency=self.efficiency,
                process_type=self.process_type,
                ODE_solver=self.ODE_solver,
                ODE_tolerance=self.ODE_tolerance,
            )

            # Store solution
            self.state_in = {key: value[0] for key, value in states_up.items()}
            self.states = {
                key: np.concatenate((np.flipud(states_up[key][:-1]), states_down[key]))
                for key in states_up
            }
            self.ode_solution = (ode_solution_up, ode_solution_down)

        else:
            msg = f"Invalid value for model_type={self.model_type}. Something went wrong during input validation."
            raise ValueError(msg)

    def fit_polynomials(self):
        """
        Fits polynomials to the states using the PolynomialFitter class.

        See Also
        --------
        PolynomialFitter :
            Class used for generating fitting polynomials.
        """

        if not self.states and not self.ode_solution:
            self.solve()
            # msg = "Fluid properties are not computed. Please call the 'solve()' method first."
            # raise ValueError(msg)

        # Fit the polynomials
        self.poly_fitter = PolynomialFitter(
            states=self.states,
            state_in=self.state_in,
            variables=self.polynomial_variables,
            degree=self.polynomial_degree,
            model_type=self.model_type,
            calculation_type=self.calculation_type,
            output_dir=self.output_dir,
        )
        self.poly_fitter.fit_polynomials()

        # Automatically initialize the ExpressionExporter
        self.exporter = ExpressionExporter(
            poly_fitter=self.poly_fitter,
            poly_format=self.polynomial_format,
        )

    def export_expressions_fluent(self, output_dir=None):
        """
        Exports the polynomial expressions in a format suitable for Ansys Fluent.

        Parameters
        ----------
        output_dir : str, optional
            The directory where the expressions will be saved. It uses a default directory if not provided.

        See Also
        --------
        ExpressionExporter :
            Class for exporting polynomial expressions for use in CFD software.
        """
        if not self.exporter:
            msg = "ExpressionExporter not initialized. Please call the method 'fit_polynomials()' first."
            raise ValueError(msg)

        output_dir = output_dir or self.output_dir
        self.exporter.export_expressions_fluent(output_dir=output_dir)

    def export_expressions_cfx(self, output_dir=None):
        """
        Exports the polynomial expressions in a format suitable for Ansys CFX.

        Parameters
        ----------
        output_dir : str, optional
            The directory where the expressions will be saved. It uses a default directory if not provided.

        See Also
        --------
        ExpressionExporter :
            Class for exporting polynomial expressions for use in CFD software.
        """
        if not self.exporter:
            msg = "ExpressionExporter not initialized. Please call the method 'fit_polynomials()' first."
            raise ValueError(msg)

        output_dir = output_dir or self.output_dir
        self.exporter.export_expressions_cfx(output_dir=output_dir)


def barotropic_model_one_component(
    fluid_name,
    p_in,
    rho_in,
    p_out,
    efficiency,
    process_type=None,
    calculation_type=None,
    blending_onset=None,
    blending_width=None,
    HEOS_solver="hybr",
    HEOS_tolerance=1e-6,
    HEOS_max_iter=100,
    HEOS_print_convergence=False,
    ODE_solver="lsoda",
    ODE_tolerance=1e-8,
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

    p_in : float
        Inlet pressure of the fluid in Pascal.

    rho_in : float
        Inlet density of the fluid in kilogram per cubic meter.

    p_out : float
        Outlet pressure of the fluid in Pascals.

    efficiency : float
        The efficiency of the polytropic process, dimensionless.

    calculation_type : str
        The type of calculation to perform. Options include:

        - ``equilibrium``: Computes equilibrium properties only.
        - ``metastable``: Computes metastable properties only.
        - ``blending``: Computes both equilibrium and metastable properties, blends them, and returns the blended properties.

    blending_onset : float, optional
        The onset of blending in the process, typically a value between 0 and 1. Required when `calculation_type` is ``blending``.
    blending_width : float, optional
        The width of the blending region, typically a value between 0 and 1. Required when `calculation_type` is ``blending``.

    HEOS_solver : str, optional
        The solver algorithm used to compute the metastable states. Valid options:

        .. list-table::
            :widths: 20 50
            :header-rows: 1

            * - Solver name
              - Description
            * - ``hybr``
              -  Powell's hybrid trust-region algorithm
            * - ``lm``
              - Levenberg–Marquardt algorithm

        See `Scipy root() <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root.html>`_ for more info.
        Recommended solvers: Both ``hybr`` and ``lm`` work well for the tested cases.

    HEOS_tolerance : float, optional
        The tolerance for the HEOS solver.

    HEOS_max_iter : int, optional
        The maximum number of iterations for the HEOS solver.

    HEOS_print_convergence : bool, optional
        If True, prints convergence information for the HEOS solver.

    ODE_solver : str, optional
        The solver to use for the ODE integration. Valid options:

        .. list-table::
            :widths: 20 50
            :header-rows: 1

            * - Solver name
              - Description
            * - ``RK23``
              - Explicit Runge-Kutta method of order 3(2)
            * - ``RK45``
              - Explicit Runge-Kutta method of order 5(4)
            * - ``DOP853``
              - Explicit Runge-Kutta method of order 8
            * - ``Radau``
              - Implicit Runge-Kutta method of the Radau IIA family of order 5
            * - ``BDF``
              - Implicit multi-step variable-order (1 to 5) method based on a backward differentiation formula for the derivative approximation
            * - ``LSODA``
              - Adams/BDF method with automatic stiffness detection and switching

        See `Scipy solver_ivp() <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html>`_  for more info.
        Recommended solvers: ``BDF``, ``LSODA``, or ``Radau`` for stiff problems or ``RK45`` for non-stiff problems with smooth blending.

    ODE_tolerance : float, optional
        The relative and absolute tolerance for the ODE solver.

    Returns
    -------
    states : dictionary of arrays
        A dictionary of Numpy arrays representing the properties of the fluid at each evaluation point.
    solution : scipy.integrate.OdeResult
        The result of the ODE integration containing information about the solver process.
    """

    # Check if the provided HEOS solver is valid
    valid_solvers_heos = psv.nonlinear_system.SOLVER_OPTIONS
    if HEOS_solver not in valid_solvers_heos:
        error_message = (
            f"Invalid HEOS solver '{HEOS_solver}' provided. "
            f"Valid solvers are: {', '.join(valid_solvers_heos)}. "
        )
        raise ValueError(error_message)

    # Check if the provided ODE solver is valid
    valid_solvers_ode = list(ODE_METHODS.keys())
    if ODE_solver not in valid_solvers_ode:
        error_message = (
            f"Invalid ODE solver '{ODE_solver}' provided. "
            f"Valid solver are: {', '.join(valid_solvers_ode)}. "
            "Recommended solvers: 'BDF', 'LSODA' or 'Radau' stiff problems involving equilibrium property calculations or blending calculations with a narrow blending width. 'RK45' can be used for non-stiff problems with a wide blending width."
        )
        raise ValueError(error_message)
    
    # Pre-process efficiency value
    if not (0.0 <= efficiency <= 1.0):
        raise ValueError(f"Efficiency must be between 0 and 1. Provided: {efficiency:.3f}")

    if process_type == "compression":
        if efficiency == 0.0:
            raise ValueError("Efficiency cannot be zero for compression (division by zero).")
        efficiency = 1 / efficiency
    elif process_type != "expansion":
        raise ValueError(f"Invalid process_type='{process_type}'. Must be 'expansion' or 'compression'.")


    # Initialize fluid and compute inlet state
    fluid = props.Fluid(name=fluid_name, backend="HEOS", exceptions=True)
    state_in = fluid.get_state(
        props.DmassP_INPUTS, rho_in, p_in, supersaturation=True, generalize_quality=True
    )

    # Determine the type of phase change process
    phase_change_type = (
        "flashing" if state_in.s < fluid.critical_point.s else "condensation"
    )

    # Initial guess for the metastable state (is updated at each iteration)
    rhoT_guess_metastable = [state_in.rho, state_in.T]

    # Differential equation defining the polytropic process
    def odefun(t, y):
        nonlocal rhoT_guess_metastable  # Allows modification within the function scope

        # Rename arguments
        p = t
        (h,) = y

        # Do calculations according to specified type
        if calculation_type == "equilibrium":
            
            # Compute equilibrium thermodynamic state
            try:  # Compute equilibrium state using CoolProp solver
                state_eq = fluid.get_state(
                    input_type=props.HmassP_INPUTS,
                    prop_1=h,
                    prop_2=p,
                    generalize_quality=False,
                    supersaturation=True,
                )

            except:  # Switch custom solver if something goes wrong
                print("Warning: CoolProp solver failed, switching to custom solver...")
                state_eq = fluid.get_state_equilibrium(
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
            rhoT_guess_metastable = [state_eq.rho, state_eq.T]
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

            if not utils.is_float(blending_onset) or not utils.is_float(blending_width):
                msg = f"The variables blending_onset={blending_onset} and blending_width={blending_width} must be floats when calculation_type='blending'."
                raise ValueError(msg)

            # Compute equilibrium thermodynamic state
            try:  # Compute equilibrium state using CoolProp solver
                state_eq = fluid.get_state(
                    input_type=props.HmassP_INPUTS,
                    prop_1=h,
                    prop_2=p,
                    generalize_quality=False,
                    supersaturation=True,
                )

            except:  # Switch custom solver if something goes wrong
                print("Warning: CoolProp solver failed, switching to custom solver...")
                state_eq = fluid.get_state_equilibrium(
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

            # Determine if metastable state should be used (avoid metastable computations too deep into the 2-phase region)
            q = state_eq["vapor_quality"]
            use_metastable = (
                (phase_change_type == "flashing" and q <= 1.0 * (blending_onset + blending_width)) or
                (phase_change_type == "condensation" and q >= (blending_onset - blending_width))
            )

            # Compute metastable state or use equilibrium state
            # Do the calculation only within the blending region to avoid using the HEOS solver too deep into the two-phase region
            if use_metastable:
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
                # Use only previous metastable calculation as initial guess (never equilibrium calculation)
                # Some ODE solvers do corrective steps at return to higher pressures after doing equilibrium calculations
                # The custom HEOS solver may fail if it uses the equilibrium state as initial guess because the density changes significantly
                rhoT_guess_metastable = [state_meta.rho, state_meta.T]

            else:
                # Take the equilibrium state outside the blending region
                state_meta = state_eq

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
                f"Valid options are: {', '.join(CALCULATION_TYPES)}."
            )

    # Solve polytropic expansion differential equation
    ode_sol = solve_ivp(
        fun=lambda p, h: odefun(p, h)[0],  # Get only first output
        t_span=[p_in, p_out],
        y0=[state_in.h],
        method=ODE_solver,
        rtol=ODE_tolerance,
        atol=ODE_tolerance,
    )
    if not ode_sol.success:
        raise Exception(ode_sol.message)

    # Postprocess solution
    # Start postprocessing from inlet state, otherwise initial guess will be too far away
    rhoT_guess_metastable = [state_in.rho, state_in.T] 
    states = utils.postprocess_ode(ode_sol.t, ode_sol.y, odefun)

    return states, ode_sol


def barotropic_model_two_component(
    fluid_name_1,
    fluid_name_2,
    mixture_ratio,
    T_in,
    p_in,
    p_out,
    efficiency,
    process_type=None,
    ODE_solver="lsoda",
    ODE_tolerance=1e-8,
):
    """
    Simulates a polytropic process for a mixture of two different fluids.

    TODO: add model equations and explanation

    Parameters
    ----------
    fluid_name_1 : str
        The name of the first component of the mixture.
    fluid_name_2 : str
        The name of the second component of the mixture.
    mixture_ratio : float
        Mass ratio of the first to the second fluid in the mixture.
    T_in : float
        Inlet temperature of the mixture in Kelvin.
    p_in : float
        Inlet pressure of the mixture in Pascals.
    p_out : float
        Outlet pressure of the mixture in Pascals.
    efficiency : float
        The efficiency of the polytropic process, (between zero and one).
    ODE_solver : str, optional
        The solver to use for the ODE integration. Valid options:

        .. list-table::
            :widths: 20 50
            :header-rows: 1

            * - Solver name
              - Description
            * - ``RK23``
              - Explicit Runge-Kutta method of order 3(2)
            * - ``RK45``
              - Explicit Runge-Kutta method of order 5(4)
            * - ``DOP853``
              - Explicit Runge-Kutta method of order 8
            * - ``Radau``
              - Implicit Runge-Kutta method of the Radau IIA family of order 5
            * - ``BDF``
              - Implicit multi-step variable-order (1 to 5) method based on a backward differentiation formula for the derivative approximation
            * - ``LSODA``
              - Adams/BDF method with automatic stiffness detection and switching

        See `Scipy solver_ivp() <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html>`_  for more info.
        Recommended solvers: ``BDF``, ``LSODA``, or ``Radau`` for stiff problems or ``RK45`` for non-stiff problems with smooth blending.

    ODE_tolerance : float, optional
        The relative and absolute tolerance for the ODE solver.

    Returns
    -------
    states : dictionary of arrays
        A dictionary of Numpy arrays representing the properties of the fluid at each evaluation point.
    solution : scipy.integrate.OdeResult
        The result of the ODE integration containing information about the solver process.
    """

    # Check if the provided ODE_solver is valid
    valid_solvers_ode = list(ODE_METHODS.keys())
    if ODE_solver not in valid_solvers_ode:
        error_message = (
            f"Invalid ODE solver '{ODE_solver}' provided. "
            f"Valid solver are: {', '.join(valid_solvers_ode)}."
        )
        raise ValueError(error_message)

    # Pre-process efficiency value
    if not (0.0 <= efficiency <= 1.0):
        raise ValueError(f"Efficiency must be between 0 and 1. Provided: {efficiency:.3f}")

    if process_type == "compression":
        if efficiency == 0.0:
            raise ValueError("Efficiency cannot be zero for compression (division by zero).")
        efficiency = 1 / efficiency
    elif process_type != "expansion":
        raise ValueError(f"Invalid process_type='{process_type}'. Must be 'expansion' or 'compression'.")

    # Calculate mass fractions of each component (constant values)
    y_1 = mixture_ratio / (1 + mixture_ratio)
    y_2 = 1 / (1 + mixture_ratio)

    # Initialize fluid and compute inlet state
    fluid_1 = props.Fluid(name=fluid_name_1, backend="HEOS", exceptions=True)
    fluid_2 = props.Fluid(name=fluid_name_2, backend="HEOS", exceptions=True)

    # Compute the inlet enthalpy of the mixture (ODE initial value)
    props_in_1 = fluid_1.get_state(props.PT_INPUTS, p_in, T_in)
    props_in_2 = fluid_2.get_state(props.PT_INPUTS, p_in, T_in)
    h_in = y_1 * props_in_1.h + y_2 * props_in_2.h

    # Define the ODE system
    def odefun(t, y):

        # Rename arguments
        p = t
        h, T = y

        # Compute fluid states
        state_1 = fluid_1.get_state(props.PT_INPUTS, p, T)
        state_2 = fluid_2.get_state(props.PT_INPUTS, p, T)

        # Compute mixture thermodynamic properties
        state = props.calculate_mixture_properties(state_1, state_2, y_1, y_2)

        # Add individual phases to the mixture properties
        for key, value in state_1.items():
            state[f"{key}_1"] = value
        for key, value in state_2.items():
            state[f"{key}_2"] = value

        # Compute right-hand-side of the ODE
        dhdp = efficiency / state["rho"]
        dTdp = (dhdp - state["dhdp_T"]) / state["cp"]

        return [dhdp, dTdp], state

    # Solve polytropic expansion differential equation
    ode_sol = solve_ivp(
        fun=lambda p, h: odefun(p, h)[0],  # Get only first output
        t_span=[p_in, p_out],
        y0=[h_in, T_in],
        method=ODE_solver,
        rtol=ODE_tolerance,
        atol=ODE_tolerance,
    )
    if not ode_sol.success:
        raise Exception(ode_sol.message)

    # Postprocess solution
    states = utils.postprocess_ode(ode_sol.t, ode_sol.y, odefun)

    return states, ode_sol


def _ensure_data_available(method):
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


class PolynomialFitter:
    """
    # TODO update docstring
    Fits polynomials to the thermodynamic properties of a fluid across various states.

    polynomial_degree : int
        Degree of the polynomials to fit.

        .. note::

            When `calculation_type` is ``blending`` the degree of the polynomial in the blending region is set to 4 to achieve sufficient accuracy while preventing numerical round-off errors associated with single-precision arithmetic in CFD solvers.

    polynomial_format : str, optional
        Type of polynomial representation (``horner`` or ``standard``). Default is 'horner'.

    polynomial_variables : list of str
        A list of variable names to fit polynomials to, such as 'density', 'viscosity', 'speed_sound', 'void_fraction', 'vapor_quality'.

    output_dir : str
        The directory where output will be saved.

    """

    def __init__(
        self,
        states,
        state_in,
        variables,
        degree,
        calculation_type,
        model_type,
        output_dir="barotropic_model",
    ):

        # Rename arguments
        self.states = states
        self.state_in = state_in
        self.variables = variables
        self.poly_degree = degree
        self.model_type = model_type
        self.calculation_type = calculation_type
        self.output_dir_default = output_dir

        # Initialize variables
        self.poly_handles = {}
        self.poly_breakpoints = None
        self.p_in = self.states["p"][0]
        self.p_out = self.states["p"][-1]
        self.p_in_scaled = self.p_in / self.p_in
        self.p_out_scaled = self.p_out / self.p_in

    def fit_polynomials(self):
        """
        Fits polynomials to the data obtained from the ODE solution based on the specified calculation type.
        """
        if self.calculation_type == "equilibrium":
            self._fit_equilibrium()
        elif self.calculation_type == "blending":
            self._fit_blending()
        elif self.calculation_type in ["metastable"] or self.model_type == "two-component":
            self._fit_single_segment()
        else:
            raise ValueError(
                f"Invalid calculation_type='{self.calculation_type}'. "
                f"Valid options are: {', '.join(CALCULATION_TYPES)}."
            )
        
    def _fit_single_segment(self):
        """Fits a single polynomial segment when 'calculation_type' is not specified'"""
        # Scale pressure by inlet pressure to improve polynomial conditioning
        p_scaled = self.states["p"] / self.p_in

        # Determine the polynomial limits
        self.poly_breakpoints = [self.p_in_scaled, self.p_out_scaled]
        self.poly_degree = 6

        # Fit polynomials to data
        for var in self.variables:
            y = np.array(self.states[var])
            poly = Polynomial.fit(p_scaled, y, deg=self.poly_degree).convert()
            self.poly_handles[var] = [poly]

    def _fit_equilibrium(self):
        """
        Fits polynomials for calculation_type="equilibrium".

        This method scales pressure, identifies the phase transition, and fits polynomials for the specified variables based on single-phase and two-phase regions.
        """
        # Scale pressure by inlet pressure to improve polynomial conditioning
        p_scaled = self.states["p"] / self.p_in

        # Determine points within the two-phase region
        # eps = 1e-9
        # mask_1phase = np.abs(self.states["supersaturation_degree"]) > eps
        # mask_2phase = np.abs(self.states["supersaturation_degree"]) <= eps
        mask_1phase = ~self.states["is_two_phase"]
        mask_2phase = self.states["is_two_phase"]

        # Determine the phase change pressure if it exists and the breakpoints for the different polynomial segments
        if mask_2phase.any():
            p_transition = p_scaled[np.where(mask_2phase)[0][0]]
            self.poly_breakpoints = [self.p_in_scaled, p_transition, self.p_out_scaled]
        else:
            self.poly_breakpoints = [self.p_in_scaled, self.p_out_scaled]

        # TODO Clean
        if isinstance(self.poly_degree, list) and len(self.poly_degree) == 2:
            degree_1, degree_2 = self.poly_degree
        elif isinstance(self.poly_degree, (int, float)):
            degree_1 = degree_2 = self.poly_degree
        else:
            print("poly_degree must be a scalar (int or float) or a list with 2 items.")
            print("Switching to default values")
            degree_1 = 4
            degree_2 = 4
            # raise ValueError("poly_degree must be a scalar (int or float) or a list with 2 items.")

        # Fit polynomials to data
        for var in self.variables:
            self.poly_handles[var] = []

            if mask_1phase.any():  # Single-phase data
                y_1p = np.array(self.states[var])[mask_1phase]
                p_1p = p_scaled[mask_1phase]
                poly_1p = Polynomial.fit(p_1p, y_1p, deg=degree_1).convert()
                self.poly_handles[var].append(poly_1p)

            if mask_2phase.any():  # Two-phase data
                y_2p = np.array(self.states[var])[mask_2phase]
                p_2p = p_scaled[mask_2phase]
                poly_2p = Polynomial.fit(p_2p, y_2p, deg=degree_2).convert()

                # # Ensure continuity by adjusting the first coefficient
                # y_1p_breakpoint = poly_1p(p_transition)
                # y_2p_breakpoint = poly_2p(p_transition)
                # coeffs_2p = poly_2p.coef
                # coeffs_2p[0] -= y_2p_breakpoint - y_1p_breakpoint
                # poly_2p = Polynomial(coeffs_2p)
                self.poly_handles[var].append(poly_2p)
                

    def _fit_blending(self):
        """
        Fits polynomials for calculation_type="blending".

        This method identifies different regions based on the blending parameter "x" and fits polynomials to ensure continuity across regions.
        """

        # Get the scaled pressure values
        p = self.states["p"] / self.p_in

        # Define masks for different regions based on states["x"]
        # limit = 0.96
        # limit = 0.8
        limit = 1
        mask_region_1 = self.states["x"] >= limit
        mask_region_2 = (self.states["x"] > 0) & (self.states["x"] < limit)
        mask_region_3 = self.states["x"] <= 0

        # Determine breakpoints in terms of scaled pressure
        p_region_2_start = p[mask_region_2][0]
        p_region_2_end = p[mask_region_2][-1]
        self.poly_breakpoints = [
            self.p_in_scaled,
            p_region_2_start,
            p_region_2_end,
            self.p_out_scaled,
        ]

        # Fit polynomials to data
        for var in self.variables:
            self.poly_handles[var] = []

            # TODO Clean
            if isinstance(self.poly_degree, list) and len(self.poly_degree) == 3:
                degree_1, degree_2, degree_3 = self.poly_degree
            elif isinstance(self.poly_degree, (int, float)):
                degree_1 = degree_2 = degree_3 = self.poly_degree
            else:
                raise ValueError("The variable poly_degree must be a scalar or a list with 3 items.")

            # Fit polynomial for Region 1 (x > 1)
            y_1 = np.array(self.states[var])[mask_region_1]
            p_1 = p[mask_region_1]
            degree_1 = min(len(y_1) - 1, degree_1)
            poly_1 = Polynomial.fit(p_1, y_1, deg=degree_1).convert()
            self.poly_handles[var].append(poly_1)

            # Fit polynomial for Region 2 (0 <= x <= 1)
            y_2 = np.array(self.states[var])[mask_region_2]
            p_2 = p[mask_region_2]
            degree_2 = min(len(y_2) - 1, degree_2)
            poly_2 = Polynomial.fit(p_2, y_2, deg=degree_2).convert()

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
            degree_3 = min(len(y_3) - 1, degree_3)
            poly_3 = Polynomial.fit(p_3, y_3, deg=degree_3).convert()

            # Ensure continuity by adjusting the first coefficient
            y_2_end = poly_2(p_region_2_end)
            y_3_start = poly_3(p_region_2_end)
            coeffs_3 = poly_3.coef
            coeffs_3[0] -= y_3_start - y_2_end
            poly_3 = Polynomial(coeffs_3)
            self.poly_handles[var].append(poly_3)

    @_ensure_data_available
    def evaluate_polynomial(self, p, variable):
        r"""
        Evaluates the polynomials for a given variable at a specified pressure values.

        The function evaluates the polynomial values at a given physical pressure `p`
        for a variable (e.g., density). It automatically determines the correct branch
        of the piecewise polynomial based on the breakpoints, and applies safeguards for 
        pressures outside the limits to ensure sensible values and continuity.

        The evaluation follows a piecewise definition depending on the range of :math:`\hat{p}=p/p_{\text{in}}`:

        .. math::

            \phi(\hat{p}) = 
            \begin{cases}
            a_1 \, e^\left(\frac{\hat{p} - \hat{p}_{\text{out}}}{a_2}\right) & \text{ for } & \hat{p} < \hat{p}_{\text{out}} \\
            \sum_{i=0}^{d} b_{i,1} \, \hat{p}^i & \text{ for }  & \hat{p}_{\text{ out}} \leq \hat{p} \leq \hat{p}_{1} \\
            & \;\; \vdots  & \\
            \sum_{i=0}^{d} b_{i,n} \, \hat{p}^i & \text{ for }  & \hat{p}_{n} \leq \hat{p} \leq \hat{p}_{\text{in}} \\
            c_1 + c_2 \left(\hat{p} - \hat{p}_{\text{in}}\right) & \text{ for }  & \hat{p} > \hat{p}_{\text{in}}
            \end{cases}

        where:

        - :math:`\phi(p)` is the variable being fitted (e.g., density, viscosity).
        - :math:`a_1` and :math:`a_2` are constants of the exponential function used to extrapolate the properties below the outlet pressure. The numerical values are determined to match the polynomial value and its first derivative at the endpoint :math:`\hat{p}=\hat{p}_{\text{out}}`.
        - :math:`b_{i,j}` is the :math:`i`-th polynomial coefficient of the :math:`j`-th polynomial segment, where :math:`d` is the degree of the polynomials, :math:`n` is the number of breakpoints, and :math:`n+1` is the number of polynomial segments. Additionally, :math:`[\hat{p}_{\text{out}},\, \hat{p}_{1},\,\ldots,\, \hat{p}_{n},\,\hat{p}_{\text{in}}]` are normalized pressures at the breakpoints.
        - :math:`c_1` and :math:`c_2` are constants of the linear function used to extrapolate the properties above the inlet pressure. The numerical values are determined to match the polynomial value and its first derivative at the endpoint :math:`\hat{p}=\hat{p}_{\text{in}}`.

        .. note::

            The safeguards for pressures outside the specified limits help to avoid numerical issues by applying 
            exponential decay for pressures below :math:`p_{\text{out}}` and linear extrapolation for pressures above 
            :math:`p_{\text{in}}`. These measures prevent properties from becoming negative or excessively large, which 
            can occur with very high or low (even negative) pressures during internal iterations of CFD solvers.


        Parameters
        ----------
        p : array-like
            The physical pressure values at which to evaluate the polynomials. 
            Given in Pascals (Pa).
        variable : str
            The variable for which to evaluate the polynomial, such as 'density', 'viscosity', 
            'speed_of_sound', or 'void_fraction'.

        Returns
        -------
        var_values : ndarray
            The evaluated values of the variable at the given physical pressures.

        """

        # Note:
        # The argument to the polynomials in self.poly_handles() should be the normalized pressures (not physical pressures)
        # The values in self.poly_breakpoints are normalized pressures (not physical pressures)

        # Convert input to array and scale
        p = np.atleast_1d(p)
        p_scaled = p / self.p_in

        # Reverse order of breakpoints if pressure is ascending (compression vs expansion)
        bps = np.array(self.poly_breakpoints)
        bps = bps[::-1] if self.p_in < self.p_out else bps

        # Loop over polynomial segments
        var_values = np.zeros_like(p, dtype=float)
        if variable in self.poly_handles.keys():
            for i, poly in enumerate(self.poly_handles.get(variable, [])):
                mask = (p_scaled <= bps[i]) & (p_scaled > bps[i + 1])
                if np.any(mask):
                    var_values[mask] = poly(p_scaled[mask])
        else:
            raise ValueError(f"No polynomials found for variable '{variable}'")

        # Safeguard for pressures lower than the lowest breakpoint (exponential decay)
        p_min_scaled = min(self.p_in_scaled, self.p_out_scaled)
        mask_lower = p_scaled <= p_min_scaled
        if np.any(mask_lower):
            poly_last = self.poly_handles[variable][-1]
            a1 = poly_last(p_min_scaled)
            a2 = a1 / poly_last.deriv(m=1)(p_min_scaled)
            exp_values = a1 * np.exp((p_scaled[mask_lower] - p_min_scaled) / a2)
            var_values[mask_lower] = exp_values

        # Safeguard for pressures higher than the highest breakpoint (linear extrapolation)
        p_max_scaled = max(self.p_in_scaled, self.p_out_scaled)
        mask_upper = p_scaled > p_max_scaled 
        if np.any(mask_upper):
            poly_first = self.poly_handles[variable][0]
            b1 = poly_first(p_max_scaled )
            b2 = poly_first.deriv(m=1)(p_max_scaled )
            extrap = b1 + b2 * (p_scaled[mask_upper] - p_max_scaled)
            var_values[mask_upper] = extrap

        return var_values

    @_ensure_data_available
    def plot_polynomial(
        self, var, p_eval=None, showfig=True, savefig=False, output_dir=None
    ):
        """
        Plots the polynomials for the specified variable. If the pressure values for plotting are not provided, a suitable range is generated internally.

        Parameters
        ----------
        var : str
            The variable for which to plot the polynomial.

        p_eval : array-like, optional
            The pressure values at which to plot the polynomial. If not provided, automatic values are generated.

        showfig : bool, optional
            If True, displays the plot.

        savefig : bool, optional
            If True, saves the plot to a file.

        output_dir : str, optional
            The directory where the plot will be saved if `savefig` is True. If not specified, uses the default output directory of the class.

        Returns
        -------
        fig : matplotlib.figure.Figure
            The figure object containing the plotted data.
        """

        # If output_dir is not given, set it to the default directory
        output_dir = output_dir or self.output_dir_default
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        # Create a figure and axis
        fig, ax = plt.subplots(figsize=(6, 5))

        # Set axis labels
        ax.set_xlabel("Pressure (Pa)")
        ax.set_ylabel(props.LABEL_MAPPING.get(var, var))

        # Plot specified values or segments between breakpoints
        if p_eval is not None:

            # Plot polynomial at specified pressures
            prop_eval = self.evaluate_polynomial(p_eval, var)
            ax.plot(p_eval, prop_eval, color=COLORS_MATLAB[0])

        else:

            # Reverse order of breakpoints if pressure is ascending (compression vs expansion)
            bps = np.array(self.poly_breakpoints)
            bps = bps[::-1] if self.p_in < self.p_out else bps

            # Concatenate to form the complete list of breakpoints
            breakpoints = [0.8 * bps[0]] + self.poly_breakpoints + [1.2 * bps[-1]]

            # Iterate over the breakpoints and plot the segments
            for i in range(len(breakpoints) - 1):

                # Generate segments between breakpoints
                bp1 = breakpoints[i]
                bp2 = breakpoints[i + 1]
                p_segment = np.linspace(bp1, bp2, 100) * self.p_in

                # Evaluate the polynomial over the segment
                prop_eval = self.evaluate_polynomial(p_segment, var)
                ax.plot(p_segment, prop_eval, color=COLORS_MATLAB[i])

                # Plot points just before breakpoints to check continuity
                eps = 1e-12
                x_minus = (bp2 - eps) * self.p_in
                y_minus = self.evaluate_polynomial(x_minus, var)
                ax.plot(x_minus, y_minus, marker="o", color=COLORS_MATLAB[i])

                # Plot points just after breakpoints to check continuity
                x_plus = (bp2 + eps) * self.p_in
                y_plus = self.evaluate_polynomial(x_plus, var)
                ax.plot(x_plus, y_plus, marker="o", color=COLORS_MATLAB[i])

        plt.tight_layout(pad=1)

        if savefig:
            file_path = os.path.join(output_dir, f"barotropic_model_{var}")
            graphics.savefig_in_formats(fig, file_path)

        if not showfig:
            plt.close(fig)

        return fig

    @_ensure_data_available
    def plot_polynomial_and_error(
        self, var, num_points=500, showfig=True, savefig=False, output_dir=None
    ):
        r"""
        Plots the barotropic polynomials and original data points from the ODE solution to illustrate
        the quality of the fit. Additionally, plots the relative error between the polynomial and the data.

        Parameters
        ----------
        var : str
            The variable to plot (e.g., 'density', 'viscosity', etc.).

        num_points : int, optional
            Number of plot points for each of the polynomial segments.

        showfig : bool, optional
            If True, displays the plot.

        savefig : bool, optional
            If True, saves the plot to a file.

        output_dir : str, optional
            The directory where the plot will be saved if `savefig` is True. If not specified, uses the default output directory of the class.

        Returns
        -------
        fig : matplotlib.figure.Figure
            The figure object containing the plotted data and error.
        """
        # If output_dir is not given, set it to the default directory
        output_dir = output_dir or self.output_dir_default
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        # Create a figure with two vertically stacked subplots
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(6, 5))

        # Mapping variable names to axis labels (variable name if not found)
        x_label = "Pressure (Pa)"
        y_label = props.LABEL_MAPPING.get(var, var)

        # Generate pressure values for the state points and calculate relative error
        p_states = self.states["p"]
        prop_states = self.states[var]
        relative_error = (
            100 * (self.evaluate_polynomial(p_states, var) - prop_states) / (np.abs(prop_states) + 1e-6)
        )

        # Plot polynomial and data segments for each breakpoint
        extrap_lower = [0.0 * self.poly_breakpoints[-1]]
        extrap_upper = [1.2 * self.poly_breakpoints[0]]
        breakpoints = extrap_upper + self.poly_breakpoints + extrap_lower

        for i in range(len(breakpoints) - 1):
            # Define segment range between breakpoints
            bp1, bp2 = breakpoints[i], breakpoints[i + 1]
            p_segment = np.linspace(bp1, bp2, num_points) * self.p_in

            # Evaluate polynomial on segment
            prop_poly_segment = self.evaluate_polynomial(p_segment, var)
            ax1.plot(
                p_segment,
                prop_poly_segment,
                linestyle="-",
                # marker="+",
                color=COLORS_MATLAB[1],
            )

        # Plot original data points
        ax1.plot(
            p_states,
            prop_states,
            label="Original values",
            marker="o",
            markersize=2.5,
            linestyle="none",
            color=COLORS_MATLAB[0],
        )
        ax1.plot(
            self.state_in["p"],
            self.state_in[var],
            color="k",
            linewidth=1.25,
            marker="o",
            markersize=2.5,
            markerfacecolor="w",
        )
        ax1.set_ylabel(y_label)
        ax1.legend(loc="lower right")

        # Automatically adjust y-limits if data is nearly constant
        ymin, ymax = np.min(prop_states), np.max(prop_states)
        yrange = ymax - ymin
        if yrange < 1e-6 or yrange / np.abs(np.mean(prop_states)) < 0.01:
            ycenter = np.mean(prop_states)
            ymargin = 0.1 * np.abs(ycenter)
            ax1.set_ylim(ycenter - ymargin, ycenter + ymargin)


        # Plot relative error
        ax2.plot(
            p_states,
            relative_error,
            label="Relative Error",
            linestyle="-",
            marker="none",
            markersize=2.5,
            color=COLORS_MATLAB[1],
        )
        ax2.plot(
            p_states,
            relative_error,
            linestyle="none",
            marker="o",
            markersize=2.5,
            color=COLORS_MATLAB[0],
        )
        ax2.set_xlabel(x_label)
        ax2.set_ylabel("Relative error (%)")
        plt.tight_layout(pad=1)

        if savefig:
            file_path = os.path.join(output_dir, f"barotropic_model_error_{var}")
            graphics.savefig_in_formats(fig, file_path)

        if not showfig:
            plt.close(fig)

        return fig

    @_ensure_data_available
    def plot_phase_diagram(self, fluid, var_x, var_y, savefig=True, showfig=True, output_dir=None, plot_spinodal_line=True, plot_quality_isolines=True):
        """
        Plots the barotropic process in the phase diagram of the specified fluid.

        .. note::
            This function is applicable to single-component systems

        Parameters
        ----------
        fluid : Fluid
            Fluid object used to plot the phase diagram.
        var_x : str
            Variable for the x-axis (e.g., "s" for entropy).
        var_y : str
            Variable for the y-axis (e.g., "T" for temperature).
        savefig : bool, optional
            If True, saves the figure to the specified directory. Default is True.
        showfig : bool, optional
            If True, displays the figure after plotting. Default is True.
        output_dir : str, optional
            Directory where the plot will be saved if `savefig` is True. If not specified, 
            uses the default output directory of the class.

        Returns
        -------
        fig : matplotlib.figure.Figure
            The figure object containing the plotted phase diagram.


        """
        # If output_dir is not given, set it to the default directory
        output_dir = output_dir or self.output_dir_default
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        # Create a figure and axis
        fig, ax = plt.subplots(figsize=(6.0, 5.0))
        fig.suptitle(
            f"Barotropic process for {fluid.name}",
            fontsize=14,
            y=0.95,  # Adjust this value to control the distance
        )

        # Apply the mapping if it exists; otherwise, use the variable name
        ax.set_xlabel(props.LABEL_MAPPING.get(var_x, var_x))
        ax.set_ylabel(props.LABEL_MAPPING.get(var_y, var_y))

        # Plot phase diagram for the first subplot
        fig, ax = fluid.plot_phase_diagram(
            var_x,
            var_y,
            axes=ax,
            plot_critical_point=True,
            plot_saturation_line=True,
            plot_spinodal_line=plot_spinodal_line,
            plot_quality_isolines=plot_quality_isolines,
            N=50,
        )

        # Plots contour of void fraction
        if plot_quality_isolines:
            var_z = "vapor_quality"
            # var_z = "void_fraction"
            prop_dict = props.compute_quality_grid(fluid, num_points=100, quality_levels=np.linspace(0.0, 1.0, 100))
            range_z = np.linspace(0.00, 1, 11)
            colors = plt.cm.Blues(np.linspace(0.25, 0.75, 256)[::-1])
            colormap = mcolors.LinearSegmentedColormap.from_list("", colors)
            contour = ax.contourf(
                prop_dict[var_x],
                prop_dict[var_y],
                prop_dict[var_z],
                range_z,
                cmap=colormap,
            )

            # Add contour lines with black edges
            contour_lines = ax.contour(
                prop_dict[var_x],
                prop_dict[var_y],
                prop_dict[var_z],
                range_z,
                colors='black',
                linewidths=0.5,
            )

        # Plot the calculated states
        ax.plot(
            self.states[var_x],
            self.states[var_y],
            color=graphics.COLORS_MATLAB[1],
            linewidth=1.25,
        )
        ax.plot(
            self.states[var_x][0],
            self.states[var_y][0],
            color=graphics.COLORS_MATLAB[1],
            linewidth=1.25,
            marker="o",
            markerfacecolor="w",
        )
        ax.plot(
            self.states[var_x][-1],
            self.states[var_y][-1],
            color=graphics.COLORS_MATLAB[1],
            linewidth=1.25,
            marker="o",
            markerfacecolor="w",
        )
        ax.plot(
            self.state_in[var_x],
            self.state_in[var_y],
            color="k",
            linewidth=1.25,
            marker="o",
            markerfacecolor="w",
        )

        plt.tight_layout(pad=1)

        # graphics.scale_graphics_x(fig, 1/fluid.critical_point.s, mode="multiply")
        # ax.set_xlim([0.9, 1.3])
        # graphics.scale_graphics_y(fig, 1/fluid.critical_point.T, mode="multiply")
        # ax.set_ylim([0.85, 1.05])
        if savefig:
            file_path = os.path.join(output_dir, f"barotropic_process_{fluid.name}")
            graphics.savefig_in_formats(fig, file_path)

        if not showfig:
            plt.close(fig)

        return fig


class ExpressionExporter:
    """
    Exports expressions describing the barotropic model for use within CFD solvers.

    Parameters
    ----------
    poly_fitter : PolynomialFitter
        An instance of PolynomialFitter containing the polynomial coefficients and breakpoints.

    poly_format : str, optional
        Type of polynomial representation ('horner' or 'standard').
    """

    def __init__(self, poly_fitter, poly_format="horner"):

        # Rename arguments
        self.poly_fitter = poly_fitter
        self.poly_format = poly_format

        # Define the units of the fluid properties
        self.units = {
            "density": "kg/m^3",
            "viscosity": "Pa*s",
            "speed_sound": "m/s",
            "void_fraction": "none",
            "vapor_quality": "none",
        }

        # Define syntax specifics for different solvers
        self.solver_syntax = {
            "fluent": {"if": "IF", "pressure": "AbsolutePressure"},
            "cfx": {"if": "if", "pressure": "Absolute Pressure"},
        }

        # Format of the numbers in the exported expressions
        self.num_format = "0.16e"

        # Define default output directory
        self.output_dir_default = poly_fitter.output_dir_default

    def export_expressions_fluent(self, output_dir=None):
        """
        Exports the barotropic model polynomials as expressions for use in Fluent CFD software.

        Parameters
        ----------
        output_dir : str, optional
            The directory where the expressions will be saved. It uses a default directory if not provided.
        """
        self._export_expressions(solver="fluent", output_dir=output_dir)

    def export_expressions_cfx(self, output_dir=None):
        """
        Exports the barotropic model polynomials as expressions for use in CFX CFD software.

        Parameters
        ----------
        output_dir : str, optional
            The directory where the expressions will be saved. It uses a default directory if not provided.
        """
        self._export_expressions(solver="cfx", output_dir=output_dir)

    def _export_expressions(self, solver, output_dir=None):
        """
        A generic method to generate expressions for a CFD solver, given specific keywords.

        Parameters
        ----------
        solver : str
            The solver for which to generate expressions ('fluent' or 'cfx').
        output_dir : str, optional
            The directory where the expressions will be saved. It uses a default directory if not provided.
        """
        # If output_dir is not given, set it to the default directory
        output_dir = output_dir or self.output_dir_default
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        # Define syntax based on solver
        self.syntax = self.solver_syntax[solver]

        # Create expressions
        file_path = os.path.join(output_dir, f"{solver}_expressions.txt")
        with open(file_path, "wt") as file:
            # File header
            file.write(f"{solver.upper()} expressions for barotropic properties\n")
            file.write(f"Creation datetime: {datetime.datetime.now()}\n\n")

            # Generate expressions for each property
            for var in self.poly_fitter.variables:
                if var in self.units:
                    unit = self.units[var]
                    expressions = self.generate_expressions(var, unit)
                    file.write(f"barotropic_{var}\n")
                    file.write(f"{expressions}\n\n")
                else:
                    raise ValueError(f"Missing unit for variable: {var}")

    def generate_expressions(self, var, unit):
        """
        Generates a piecewise polynomial expression with extrapolation safeguards for a given variable.

        Parameters
        ----------
        var : str
            The variable name for which to generate the expression.
        unit : str
            The unit of the variable.

        Returns
        -------
        str
            The piecewise polynomial expression with extrapolation safeguards.
        """

        # Retrieve polynomial coefficients
        p_in = self.poly_fitter.p_in
        p_out = self.poly_fitter.p_out
        p_in_scaled = self.poly_fitter.p_in_scaled
        p_out_scaled = self.poly_fitter.p_out_scaled
        polynomials = self.poly_fitter.poly_handles[var]
        breakpoints = [p_in * bp for bp in self.poly_fitter.poly_breakpoints]

        # Reverse the order of pressures if the process is a compression rather than an expansion
        if p_out > p_in:
            p_in_scaled, p_out_scaled = p_out_scaled, p_in_scaled
            breakpoints = breakpoints[::-1]
            # polynomials = polynomials[::-1]  # TODO Not sure if needed
            
        # Initialize expressions
        expression = self.polynomial_expression(polynomials[0].coef)
        for i in range(1, len(polynomials)):
            new_expression = self.polynomial_expression(polynomials[i].coef)
            expression = self.if_expression(expression, new_expression, breakpoints[i])

        # Safeguard for pressures lower than p_out (Exponential decay)
        poly_last = polynomials[-1]
        a1 = poly_last(p_out_scaled)
        a2 = a1 / poly_last.deriv(m=1)(p_out_scaled)
        exp_safeguard = f"{a1:{self.num_format}} * exp(({self.syntax['pressure']} / {p_in:{self.num_format}} [Pa] - {p_out_scaled}) / ({a2:{self.num_format}}))"
        expression = self.if_expression(expression, exp_safeguard, breakpoints[-1])

        # Safeguard for pressures higher than p_in (Linear extrapolation)
        poly_first = polynomials[0]
        c1 = poly_first(p_in_scaled)
        c2 = poly_first.deriv(m=1)(p_in_scaled)
        lin_extrap = f"{c1:{self.num_format}} + {c2:{self.num_format}} * ({self.syntax['pressure']} / {p_in:{self.num_format}} [Pa] - {p_in_scaled:{self.num_format}})"
        expression = self.if_expression(lin_extrap, expression, breakpoints[0])

        # Add fluid property unit
        if unit != "none":
            expression = f"({expression}) * 1 [{unit}]"
        return expression

    def polynomial_expression(self, coefficients):
        """
        Generates a polynomial expression string.

        Parameters
        ----------
        coefficients : array-like
            Coefficients of the polynomial.

        Returns
        -------
        str
            The polynomial expression string.
        """
        p_scale = self.poly_fitter.p_in
        if self.poly_format == "horner":
            coefficients = list(reversed(coefficients))
            polynomial_string = f"{coefficients[0]:{self.num_format}}"
            for i in range(1, len(coefficients)):
                polynomial_string = f"{coefficients[i]:{self.num_format}} + ({self.syntax['pressure']} / {p_scale:{self.num_format}} [Pa]) *\n({polynomial_string})"
        elif self.poly_format == "standard":
            terms = []
            for i in range(len(coefficients)):
                summand = f"{coefficients[i]:{self.num_format}} * ({self.syntax['pressure']} / {p_scale:{self.num_format}} [Pa])^{i}"
                if i < len(coefficients) - 1:
                    summand += " +"
                terms.append(summand)
            polynomial_string = " \n".join(terms)
        else:
            raise ValueError("The polynomial format must be 'horner' or 'standard'")
        return polynomial_string

    def if_expression(self, expression_1, expression_2, transition_pressure):
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
        return f"{self.syntax['if']}({self.syntax['pressure']} >= {transition_pressure:{self.num_format}} [Pa], \n{expression_1}, \n{expression_2})"
