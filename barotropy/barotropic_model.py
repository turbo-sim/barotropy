import os
import datetime
import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial import Polynomial

from scipy.integrate import solve_ivp
from scipy.integrate._ivp.ivp import METHODS as ODE_METHODS

from . import graphics
from . import pysolver_view as psv
from . import fluid_properties as props


COLORS_MATLAB = graphics.COLORS_MATLAB
PROCESS_TYPES = ["polytropic", "adiabatic"]
PROPERTY_CALCULATION_TYPES = ["blending", "equilibrium", "metastable"]

from functools import wraps


class BarotropicModelOneComponent:
    """
    Coordinates the simulation and polynomial fitting for a barotropic process.

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

        - ``equilibrium``: Computes equilibrium properties only.
        - ``metastable``: Computes metastable properties only.
        - ``blending``: Computes both equilibrium and metastable properties, blends them, and returns the blended properties.

        Default is None, which will raise an error.
    blending_onset : float, optional
        The onset of blending in the process, typically a value between 0 and 1. Required when `calculation_type` is ``blending``.
    blending_width : float, optional
        The width of the blending region, typically a value between 0 and 1. Required when `calculation_type` is ``blending``.
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
        Default is ``LSODA``.
        
        Recommended solvers: ``BDF``, ``LSODA``, or ``Radau`` for stiff problems or ``RK45`` for non-stiff problems with smooth blending.

    ODE_tolerance : float, optional
        The relative and absolute tolerance for the ODE solver. Default is 1e-8.
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

        See `Scipy root() <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root.html>`_ for more info. Default is ``hybr``.
        
        Recommended solvers: Both ``hybr`` and ``lm`` work well for the tested cases.
    HEOS_tolerance : float, optional
        The tolerance for the HEOS solver. Default is 1e-6.
    HEOS_max_iter : int, optional
        The maximum number of iterations for the HEOS solver. Default is 100.
    HEOS_print_convergence : bool, optional
        If True, prints convergence information for the HEOS solver. Default is False.
    polynomial_degree : int
        Degree of the polynomials to fit. Default is 8.

        .. note::

            When `calculation_type` is ``blending`` the degree of the polynomial in the blending region is set to 4 to achieve sufficient accuracy while preventing numerical round-off errors associated with single-precision arithmetic in CFD solvers.

    polynomial_format : str, optional
        Type of polynomial representation (``horner`` or ``standard``). Default is 'horner'.
    polynomial_variables : list of str
        A list of variable names to fit polynomials to, such as 'density', 'viscosity',
        'speed_sound', 'void_fraction', 'vapor_quality'.
    output_dir : str
        The directory where output will be saved.

    """

    def __init__(
        self,
        fluid_name: str,
        T_in: float,
        p_in: float,
        p_out: float,
        efficiency: float = 1.00,
        calculation_type: str = None,
        blending_onset: float = None,
        blending_width: float = None,
        ODE_solver: str = "lsoda",
        ODE_tolerance: float = 1e-8,
        HEOS_solver: str = "hybr",
        HEOS_tolerance: float = 1e-6,
        HEOS_max_iter: int = 100,
        HEOS_print_convergence: bool = False,
        polynomial_degree: int = 8,
        polynomial_format: str = "horner",
        polynomial_variables: list = [
            "density",
            "viscosity",
            "speed_sound",
            # "void_fraction",
            # "vapor_quality",
        ],
        output_dir: str = "barotropic_model",
    ):

        # Rename mandatory arguments
        self.fluid_name = fluid_name
        self.T_in = T_in
        self.p_in = p_in
        self.p_out = p_out

        # Assign optional arguments with defaults
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

    def compute_properties(self):
        """
        Solves the ODE for the barotropic process and stores the states.

        See Also
        --------
        barotropic_model_one_fluid : 
            Function used to compute the fluid properties.
        """
        # Manually pass all parameters to the ODE solver function
        self.states, self.ode_solution = barotropic_model_one_fluid(
            fluid_name=self.fluid_name,
            T_in=self.T_in,
            p_in=self.p_in,
            p_out=self.p_out,
            efficiency=self.efficiency,
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

    def fit_polynomials(self):
        """
        Fits polynomials to the states using the PolynomialFitter class.
        
        See Also
        --------
        PolynomialFitter : 
            Class used for generating fitting polynomials.
        """

        if not self.states and not self.ode_solution:
            msg = "Fluid properties are not computed. Please call the method 'compute_properties()' first."
            raise ValueError(msg)

        # Fit the polynomials
        self.poly_fitter = PolynomialFitter(
            states=self.states,
            variables=self.polynomial_variables,
            degree=self.polynomial_degree,
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

        - ``equilibrium``: Computes equilibrium properties only.
        - ``metastable``: Computes metastable properties only.
        - ``blending``: Computes both equilibrium and metastable properties, blends them, and returns the blended properties.

        Default is None, which will raise an error.
    blending_onset : float, optional
        The onset of blending in the process, typically a value between 0 and 1. Required when `calculation_type` is ``blending``.
    blending_width : float, optional
        The width of the blending region, typically a value between 0 and 1. Required when `calculation_type` is ``blending``.
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
        Default is ``LSODA``.
        
        Recommended solvers: ``BDF``, ``LSODA``, or ``Radau`` for stiff problems or ``RK45`` for non-stiff problems with smooth blending.

    ODE_tolerance : float, optional
        The relative and absolute tolerance for the ODE solver. Default is 1e-8.
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

        See `Scipy root() <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root.html>`_ for more info. Default is ``hybr``.
        
        Recommended solvers: Both ``hybr`` and ``lm`` work well for the tested cases.
    HEOS_tolerance : float, optional
        The tolerance for the HEOS solver. Default is 1e-6.
    HEOS_max_iter : int, optional
        The maximum number of iterations for the HEOS solver. Default is 100.
    HEOS_print_convergence : bool, optional
        If True, prints convergence information for the HEOS solver. Default is False.

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
        method=ODE_solver,
        rtol=ODE_tolerance,
        atol=ODE_tolerance,
    )
    if not ode_sol.success:
        raise Exception(ode_sol.message)

    # Postprocess solution
    states = postprocess_ode(ode_sol.t, ode_sol.y, odefun)

    return states, ode_sol


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
    Fits polynomials to the thermodynamic properties of a fluid across various states.

    polynomial_degree : int
        Degree of the polynomials to fit. Default is 8.

        .. note::

            When `calculation_type` is ``blending`` the degree of the polynomial in the blending region is set to 4 to achieve sufficient accuracy while preventing numerical round-off errors associated with single-precision arithmetic in CFD solvers.

    polynomial_format : str, optional
        Type of polynomial representation (``horner`` or ``standard``). Default is 'horner'.
    polynomial_variables : list of str
        A list of variable names to fit polynomials to, such as 'density', 'viscosity',
        'speed_sound', 'void_fraction', 'vapor_quality'.
    output_dir : str
        The directory where output will be saved.

    """

    def __init__(
        self,
        states,
        variables,
        degree,
        calculation_type,
        output_dir="barotropic_model",
    ):

        # Rename arguments
        self.states = states
        self.variables = variables
        self.degree = degree
        self.calculation_type = calculation_type
        self.output_dir_default = output_dir

        # Initialize variables
        self.poly_handles = {}
        self.poly_breakpoints = None
        self.p_in = self.states["p"][0]
        self.p_out = self.states["p"][-1]
        self.p_in_scaled = self.p_in / self.p_in
        self.p_out_scaled = self.p_out / self.p_in
        self.plot_labels = {
            "density": "Density (kg/m$^3$)",
            "viscosity": "Viscosity (Pa·s)",
            "speed_sound": "Speed of sound (m/s)",
            "void_fraction": "Void fraction",
            "vapor_quality": "Vapor quality",
        }

    def fit_polynomials(self):
        """
        Fits polynomials to the data obtained from the ODE solution based on the specified calculation type.
        """
        if self.calculation_type == "equilibrium":
            self._fit_equilibrium()
        elif self.calculation_type == "blending":
            self._fit_blending()
        else:
            raise ValueError("Unsupported calculation type")

    def _fit_equilibrium(self):
        """
        Fits polynomials for calculation_type="equilibrium".

        This method scales pressure, identifies the phase transition, and fits polynomials
        for the specified variables based on single-phase and two-phase regions.
        """
        # Scale pressure by inlet pressure to improve polynomial conditioning
        p_scaled = self.states["p"] / self.p_in

        # Determine points within the two-phase region
        eps = 1e-9
        mask_1phase = np.abs(self.states["supersaturation_degree"]) > eps
        mask_2phase = np.abs(self.states["supersaturation_degree"]) <= eps

        # Determine the phase change pressure if it exists and the breakpoints for the different polynomial segments
        if mask_2phase.any():
            p_transition = p_scaled[np.where(mask_2phase)[0][0]]
            self.poly_breakpoints = [self.p_in_scaled, p_transition, self.p_out_scaled]
        else:
            self.poly_breakpoints = [self.p_in_scaled, self.p_out_scaled]

        # Fit polynomials to data
        for var in self.variables:
            self.poly_handles[var] = []

            if mask_1phase.any():  # Single-phase data
                y_1p = np.array(self.states[var])[mask_1phase]
                p_1p = p_scaled[mask_1phase]
                poly_1p = Polynomial.fit(p_1p, y_1p, deg=self.degree).convert()
                self.poly_handles[var].append(poly_1p)

            if mask_2phase.any():  # Two-phase data
                y_2p = np.array(self.states[var])[mask_2phase]
                p_2p = p_scaled[mask_2phase]
                poly_2p = Polynomial.fit(p_2p, y_2p, deg=self.degree).convert()

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

        This method identifies different regions based on the blending
        parameter "x" and fits polynomials to ensure continuity across regions.
        """

        # Get the scaled pressure values
        p = self.states["p"] / self.p_in

        # Define masks for different regions based on states["x"]
        mask_region_1 = self.states["x"] > 1
        mask_region_2 = (self.states["x"] >= 0) & (self.states["x"] <= 1)
        mask_region_3 = self.states["x"] < 0

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

            # Fit polynomial for Region 1 (x > 1)
            y_1 = np.array(self.states[var])[mask_region_1]
            p_1 = p[mask_region_1]
            poly_1 = Polynomial.fit(p_1, y_1, deg=self.degree).convert()
            self.poly_handles[var].append(poly_1)

            # Fit polynomial for Region 2 (0 <= x <= 1)
            # TODO Hardcoded degree 4 to prevent numerical noise
            y_2 = np.array(self.states[var])[mask_region_2]
            p_2 = p[mask_region_2]
            poly_2 = Polynomial.fit(p_2, y_2, deg=4).convert()  

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
            poly_3 = Polynomial.fit(p_3, y_3, deg=self.degree).convert()

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
        # The argument to the polynomials in self.poly_handles() should be normalized pressures (not physical pressures)
        # The values in self.poly_breakpoints are normalized pressures (not physical pressures)

        # Convert input to array if not already
        p = np.atleast_1d(p)
        var_values = np.zeros_like(p, dtype=float)

        # Define scaled pressures for internal calculations
        p_scaled = p / self.p_in
        bps = np.array(self.poly_breakpoints)
        if variable in self.poly_handles.keys():
            for i, poly in enumerate(self.poly_handles.get(variable, [])):
                mask = (p_scaled <= bps[i]) & (p_scaled > bps[i + 1])
                if np.any(mask):
                    var_values[mask] = poly(p_scaled[mask])
        else:
            raise ValueError(f"No polynomials found for variable '{variable}'")

        # Safeguard for pressures lower than the lowest breakpoint (exponential decay)
        mask_lower = p_scaled <= self.p_out_scaled
        if np.any(mask_lower):
            poly_last = self.poly_handles[variable][-1]
            a1 = poly_last(self.p_out_scaled)
            a2 = a1 / poly_last.deriv(m=1)(self.p_out_scaled)
            exp_values = a1 * np.exp((p_scaled[mask_lower] - self.p_out_scaled) / a2)
            var_values[mask_lower] = exp_values

        # Safeguard for pressures higher than the highest breakpoint (linear extrapolation)
        mask_upper = p_scaled > self.p_in_scaled
        if np.any(mask_upper):
            poly_first = self.poly_handles[variable][0]
            b1 = poly_first(self.p_in_scaled)
            b2 = poly_first.deriv(m=1)(self.p_in_scaled)
            extrap = b1 + b2 * (p_scaled[mask_upper] - self.p_in_scaled)
            var_values[mask_upper] = extrap

        return var_values

    @_ensure_data_available
    def plot_polynomial(self, var, p_eval=None, showfig=True, savefig=False, output_dir=None):
        """
        Plots the polynomials for the specified variable. If the pressure values for plotting are not provided, a suitable range is generated internally.

        Parameters
        ----------
        var : str
            The variable for which to plot the polynomial.
        p_eval : array-like, optional
            The pressure values at which to plot the polynomial. If not provided, automatic values are generated.
        showfig : bool, optional
            If True, displays the plot. Default is True.
        savefig : bool, optional
            If True, saves the plot to a file. Default is False.
        output_dir : str, optional
            The directory where the plot will be saved if `savefig` is True.
            If not specified, uses the default output directory of the class.

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
        ax.set_ylabel(self.plot_labels.get(var, var))

        # Plot specified values or segments between breakpoints
        if p_eval is not None:

            # Plot polynomial at specified pressures
            prop_eval = self.evaluate_polynomial(p_eval, var)
            ax.plot(p_eval, prop_eval, color=COLORS_MATLAB[0])

        else:

            # Iterate over the breakpoints and plot the segments
            extrap_lower = [0.8 * self.poly_breakpoints[-1]]
            extrap_upper = [1.2 * self.poly_breakpoints[0]]
            breakpoints = extrap_upper + self.poly_breakpoints + extrap_lower
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
    def plot_polynomial_and_error(self, var, showfig=True, savefig=False, output_dir=None):
        r"""
        Plots the barotropic polynomials and original data points from the ODE solution to illustrate
        the quality of the fit. Additionally, plots the relative error between the polynomial and the data.

        Parameters
        ----------
        var : str
            The variable to plot (e.g., 'density', 'viscosity', etc.).
        showfig : bool, optional
            If True, displays the plot. Default is True.
        savefig : bool, optional
            If True, saves the plot to a file. Default is False.
        output_dir : str, optional
            The directory where the plot will be saved if `savefig` is True.
            If not specified, uses the default output directory of the class.

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
        y_label = self.plot_labels.get(var, var)

        # Generate pressure values and evaluate the variable
        pressures = self.states["p"]
        prop_poly = self.evaluate_polynomial(pressures, var)
        prop_states = self.states[var]

        # Compute the relative error
        relative_error = 100 * np.abs(prop_poly - prop_states) / np.abs(prop_states)

        # Plot the variable from polynomials and states in the upper subplot
        ax1.plot(
            pressures,
            prop_poly,
            label=f"Polynomial fit",
            linestyle="-",
            color=COLORS_MATLAB[1],
        )
        ax1.plot(
            pressures,
            prop_states,
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
            pressures,
            relative_error,
            label="Relative Error",
            linestyle="-",
            markersize=2.5,
            color=COLORS_MATLAB[1],
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


class ExpressionExporter:
    """
    Exports expressions describing the barotropic model for use within CFD solvers.

    Parameters
    ----------
    poly_fitter : PolynomialFitter
        An instance of PolynomialFitter containing the polynomial coefficients and breakpoints.
    poly_format : str, optional
        Type of polynomial representation ('horner' or 'standard'). Default is 'horner'.
    """

    def __init__(self, poly_fitter, poly_format="horner"):
        # Rename arguments
        self.poly_fitter = poly_fitter
        self.poly_format = poly_format

        # Define the units of the fluid properties
        self.units = {
            "density": "kg/m^3",
            "viscosity": "Pa*s",
            "speed_of_sound": "m/s",
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
                    file.write(f"barotropic_model_{var}\n")
                    file.write(f"{expressions}\n\n")


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
        p_in_scaled = self.poly_fitter.p_in_scaled
        p_out_scaled = self.poly_fitter.p_out_scaled
        polynomials = self.poly_fitter.poly_handles[var]
        breakpoints = [p_in * bp for bp in self.poly_fitter.poly_breakpoints]

        # Initialize expressions
        expression = self.polynomial_expression(polynomials[0].coef)
        for i in range(1, len(polynomials)):
            new_expression = self.polynomial_expression(polynomials[i].coef)
            expression = self.if_expression(expression, new_expression, breakpoints[i])

        # Safeguard for pressures lower than p_out (Exponential decay)
        poly_last = polynomials[-1]
        a1 = poly_last(p_out_scaled)
        a2 = a1 / poly_last.deriv(m=1)(p_out_scaled)
        exp_safeguard = f"{a1:{self.num_format}} * exp(({self.syntax['pressure']} / {p_in:{self.num_format}} [Pa] - {p_out_scaled}) / {a2:{self.num_format}})"
        expression = self.if_expression(expression, exp_safeguard, breakpoints[-1])

        # Safeguard for pressures higher than p_in (Linear extrapolation)
        poly_first = polynomials[0]
        c1 = poly_first(p_in_scaled)
        c2 = poly_first.deriv(m=1)(p_in_scaled)
        lin_extrap = f"{c1:{self.num_format}} + {c2:{self.num_format}} * ({self.syntax['pressure']} / {p_in:{self.num_format}} [Pa] - {p_in_scaled:{self.num_format}})"
        expression = self.if_expression(lin_extrap, expression, breakpoints[0])

        # Add fluid property unit
        if unit != "none":
            expression = f"{expression} * 1 [{unit}]"
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


