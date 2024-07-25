import numpy as np

from scipy.integrate import solve_ivp

from . import fluid_properties as props

PROCESS_TYPES = ["polytropic", "adiabatic"]
MULTIPHASE_MODELS = ["blended", "blending", "equilibrium", "metastable"]

def barotropic_model_one_fluid(
        fluid_name=None,
        T_in=None,
        p_in=None,
        p_out=None,
        efficiency=None,
        number_of_points=None,
        multiphase_model="blended",
        calculation_type="polytropic",
        tolerance=1e-6,
    ):

    # Compute process as an adiabatic or polytropic expansion
    if calculation_type == "adiabatic":
        states = calculate_adiabatic_process(
            fluid_name=fluid_name,
            T_in=T_in,
            p_in=p_in,
            p_out=p_out,
            efficiency=efficiency,
            # multiphase_model=multiphase_model,
            number_of_points=number_of_points,
        )

    elif calculation_type == "polytropic":
        states, _ = polytropic_process_one_fluid(
            fluid_name=fluid_name,
            T_in=T_in,
            p_in=p_in,
            p_out=p_out,
            efficiency=efficiency,
            multiphase_model=multiphase_model,
            number_of_points=number_of_points,
            tolerance=tolerance,
        )

    else:
        options = ", ".join(f"'{k}'" for k in PROCESS_TYPES)
        raise ValueError(
            f"Invalid process type '{calculation_type}'. Available options: {options}"
        )

    return states


def polytropic_process_one_fluid(
    fluid_name,
    T_in,
    p_in,
    p_out,
    efficiency,
    multiphase_model="blended",
    q_onset=0.00,
    q_transition=0.01,
    tolerance=1e-6,
    number_of_points=None,
):
    """ """

    # Differential equation defining the polytropic process
    def odefun(t, y):
        # Keep variables in memory to use previous state as initial guess
        nonlocal rho_guess, T_guess

        # Rename input variables
        p = t
        (h,) = y

        # Compute equilibrium state
        state_eq = fluid.get_state(props.HmassP_INPUTS, h, p, generalize_quality=True)

        # Compute metastable state
        state_meta = fluid.get_state_metastable("p", p, "hmass", h, rho_guess, T_guess)
        rho_guess = state_meta.rho  # Update guess for the next iteration
        T_guess = state_meta.T  # Update guess for the next iteration
        
        # Determine the actual thermodynamic state
        if multiphase_model == "equilibrium":
            state = state_eq
        elif multiphase_model == "metastable":
            state = state_meta
        elif multiphase_model == "blended" or multiphase_model == "blending":
            x = 1 - (state_eq["Q"] - q_onset) / q_transition
            # sigma = sigmoid_smoothstep(x)
            sigma = sigmoid_smootherstep(x)
            state = {}
            for key in state_eq.keys():
                state[key] = state_eq[key]*(1-sigma) + state_meta[key]*sigma
        else:
            options = ", ".join(f"'{k}'" for k in MULTIPHASE_MODELS)
            raise ValueError(
                f"Invalid multiphase model {multiphase_model}'. Available options: {options}"
        )

        # Compute right-hand-side of the ODE
        dhdp = efficiency / state["rho"]

        return dhdp, state
    
    # Compute inlet state
    fluid = props.Fluid(name=fluid_name, backend="HEOS", exceptions=True)
    state_in = fluid.get_state(props.PT_INPUTS, p_in, T_in)
    rho_guess = state_in.rho
    T_guess = state_in.T

    # Solve polytropic expansion differential equation
    solution = solve_ivp(
        lambda p, h: odefun(p, h)[0],
        [state_in.p, p_out],
        [state_in.h],
        t_eval=np.linspace(state_in.p, p_out, number_of_points)
        if number_of_points
        else None,
        method="LSODA",  # "RK45", "Radau", "LSODA"
        rtol=tolerance,
        atol=tolerance,
    )
    if not solution.success:
        raise Exception(solution.message)

    # Postprocess solution
    states = postprocess_ode(solution.t, solution.y, odefun)

    return states, solution


def calculate_adiabatic_process(
    fluid_name,
    T_in,
    p_in,
    p_out,
    efficiency,
    number_of_points=None,
):
    # Compute inlet state
    fluid = props.Fluid(name=fluid_name, backend="HEOS", exceptions=True)
    state_in = fluid.get_state(props.PT_INPUTS, p_in, T_in)

    # Compute outlet state
    state_out_s = fluid.get_state(props.PSmass_INPUTS, p_out, state_in.s)
    h_out = state_in.h - efficiency * (state_in.h - state_out_s.h)
    state_out = fluid.get_state(props.HmassP_INPUTS, h_out, p_out)

    # Prepare list of states
    if number_of_points is not None:
        p_array = np.logspace(np.log10(state_in.p), np.log10(state_out.p), number_of_points)
        s_array = np.linspace(state_in.s, state_out.s, number_of_points)
        states = []
        for p, s in zip(p_array, s_array):
            states.append(fluid.get_state(props.PSmass_INPUTS, p, s))
    else:
        states = [state_in, state_out]

    return props.states_to_dict(states)


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


# def sigmoid_hyperbolic(x, x0=0.5, alpha=1):
#     r"""
#     Compute the sigmoid hyperbolic function.

#     This function calculates a sigmoid function based on the hyperbolic tangent.
#     The formula is given by:

#     .. math::
#         \sigma(x) = \frac{1 + \tanh\left(\frac{x - x_0}{\alpha}\right)}{2}

#     Parameters
#     ----------
#     x : array_like
#         Input data.
#     x0 : float, optional
#         Center of the sigmoid function. Default is 0.5.
#     alpha : float, optional
#         Scale parameter. Default is 0.1.

#     Returns
#     -------
#     array_like
#         Output of the sigmoid hyperbolic function.

#     """
#     x = np.array(x)  # Ensure x is a NumPy array for vectorized operations
#     sigma = (1 + np.tanh((x - x0) / alpha)) / 2
#     return sigma


# def sigmoid_rational(x, n, m):
#     r"""
#     Compute the sigmoid rational function.

#     This function calculates a sigmoid function using an algebraic approach based on a rational function.
#     The formula is given by:

#     .. math::
#         \sigma(x) = \frac{x^n}{x^n + (1-x)^m}

#     Parameters
#     ----------
#     x : array_like
#         Input data.
#     n : int
#         Power of the numerator.
#     m : int
#         Power of the denominator.

#     Returns
#     -------
#     array_like
#         Output of the sigmoid algebraic function.

#     """
#     x = np.array(x)  # Ensure x is a NumPy array for vectorized operations
#     x = np.where(x < 0, 0, x)  # Set x to 0 where x < 0
#     x = np.where(x > 1, 1, x)  # Set x to 1 where x > 1
#     sigma = x**n / (x**n + (1 - x) ** m)
#     return sigma


# def sigmoid_smoothstep(x):
#     """
#     Compute the smooth step function.

#     This function calculates a smooth step function with zero first-order
#     derivatives at the endpoints. More information available at:
#     https://resources.wolframcloud.com/FunctionRepository/resources/SmoothStep/

#     Parameters
#     ----------
#     x : array_like
#         Input data.

#     Returns
#     -------
#     array_like
#         Output of the sigmoid smoothstep function.

#     """
#     x = np.array(x)  # Ensure x is a NumPy array for vectorized operations
#     x = np.where(x < 0, 0, x)  # Set x to 0 where x < 0
#     x = np.where(x > 1, 1, x)  # Set x to 1 where x > 1
#     return x**2 * (3 - 2 * x)


# def sigmoid_smootherstep(x):
#     """
#     Compute the smoother step function.

#     This function calculates a smoother step function with zero second-order
#     derivatives at the endpoints. It is a modification of the smoothstep
#     function to provide smoother transitions.

#     Parameters
#     ----------
#     x : array_like
#         Input data.

#     Returns
#     -------
#     array_like
#         Output of the sigmoid smootherstep function.

#     """
#     x = np.array(x)  # Ensure x is a NumPy array for vectorized operations
#     x = np.where(x < 0, 0, x)  # Set x to 0 where x < 0
#     x = np.where(x > 1, 1, x)  # Set x to 1 where x > 1
#     return 6 * x**5 - 15 * x**4 + 10 * x**3
