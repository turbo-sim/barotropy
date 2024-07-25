import scipy
import numpy as np
import CoolProp.CoolProp as CP

from .. import math as math
from .. import pysolver_view as psv



# Define property aliases
PROPERTY_ALIAS = {
    "P": "p",
    "rho": "rhomass",
    "d": "rhomass",
    "u": "umass",
    "h": "hmass",
    "s": "smass",
    "cv": "cvmass",
    "cp": "cpmass",
    "a": "speed_sound",
    "Z": "compressibility_factor",
    "mu": "viscosity",
    "k": "conductivity",
}



# ------------------------------------------------------------------------------------ #
# Equilibrum property calculations through CoolProp
# ------------------------------------------------------------------------------------ #


def compute_properties_1phase(
    abstract_state,
    generalize_quality=False,
    compute_subcooling=False,
    compute_superheating=False,
):
    """Get single-phase properties from CoolProp low level interface

    Direct call to coolprop

    TODO To be completed
    """

    # Fluid properties available from CoolProp
    props = {}
    AS = abstract_state
    props["T"] = AS.T()
    props["p"] = AS.p()
    props["rhomass"] = AS.rhomass()
    props["umass"] = AS.umass()
    props["hmass"] = AS.hmass()
    props["smass"] = AS.smass()
    props["gibbsmass"] = AS.gibbsmass()
    props["cvmass"] = AS.cvmass()
    props["cpmass"] = AS.cpmass()
    props["gamma"] = props["cpmass"] / props["cvmass"]
    props["compressibility_factor"] = AS.compressibility_factor()
    props["speed_sound"] = AS.speed_sound()
    props["isentropic_bulk_modulus"] = props["rhomass"] * props["speed_sound"] ** 2
    props["isentropic_compressibility"] = 1 / props["isentropic_bulk_modulus"]
    props["isothermal_bulk_modulus"] = 1 / AS.isothermal_compressibility()
    props["isothermal_compressibility"] = AS.isothermal_compressibility()
    isobaric_expansion_coefficient = AS.isobaric_expansion_coefficient()
    props["isobaric_expansion_coefficient"] = isobaric_expansion_coefficient
    props["viscosity"] = AS.viscosity()
    props["conductivity"] = AS.conductivity()

    # Generalized quality outside the two-phase region
    if generalize_quality:
        quality = calculate_generalized_quality(AS)
    else:
        quality = np.nan

    props["Q"] = quality
    props["quality_mass"] = np.nan
    props["quality_volume"] = np.nan

    # Calculate departure from saturation properties
    if compute_subcooling:
        props["subcooling"] = calculate_subcooling(AS)

    if compute_superheating:
        props["superheating"] = calculate_superheating(AS)

    # Add properties as aliases
    for key, value in PROPERTY_ALIAS.items():
        props[key] = props[value]

    return props


def compute_properties_2phase(
    abstract_state,
    compute_subcooling=True,
    compute_superheating=True,
):
    """Get two-phase properties from mixing rules and single-phase CoolProp properties

    Homogeneous equilibrium model

    State formulas for T=T, p=p, mfrac/vfrac(rho), h-s-g-u-cp-cv, mu-k, a
    TODO To be completed

    """

    # Instantiate new AbstractState to compute saturation properties without changing the state of the class
    AS = abstract_state
    cloned_AS = CP.AbstractState(AS.backend_name(), AS.name())

    # Basic properties of the two-phase mixture
    T_mix = AS.T()
    p_mix = AS.p()
    rho_mix = AS.rhomass()
    u_mix = AS.umass()
    h_mix = AS.hmass()
    s_mix = AS.smass()
    gibbs_mix = AS.gibbsmass()

    # Saturated liquid properties
    cloned_AS.update(CP.QT_INPUTS, 0.00, T_mix)
    rho_L = cloned_AS.rhomass()
    cp_L = cloned_AS.cpmass()
    cv_L = cloned_AS.cvmass()
    k_L = cloned_AS.conductivity()
    mu_L = cloned_AS.viscosity()
    a_L = cloned_AS.speed_sound()
    dsdp_L = cloned_AS.first_saturation_deriv(CP.iSmass, CP.iP)

    # Saturated vapor properties
    cloned_AS.update(CP.QT_INPUTS, 1.00, T_mix)
    rho_V = cloned_AS.rhomass()
    cp_V = cloned_AS.cpmass()
    cv_V = cloned_AS.cvmass()
    k_V = cloned_AS.conductivity()
    mu_V = cloned_AS.viscosity()
    a_V = cloned_AS.speed_sound()
    dsdp_V = cloned_AS.first_saturation_deriv(CP.iSmass, CP.iP)

    # Volume fractions of vapor and liquid
    vfrac_V = (rho_mix - rho_L) / (rho_V - rho_L)
    vfrac_L = 1.00 - vfrac_V

    # Mass fractions of vapor and liquid
    mfrac_V = (1 / rho_mix - 1 / rho_L) / (1 / rho_V - 1 / rho_L)
    mfrac_L = 1.00 - mfrac_V

    # Heat capacities of the two-phase mixture
    cp_mix = mfrac_L * cp_L + mfrac_V * cp_V
    cv_mix = mfrac_L * cv_L + mfrac_V * cv_V

    # Transport properties of the two-phase mixture
    k_mix = vfrac_L * k_L + vfrac_V * k_V
    mu_mix = vfrac_L * mu_L + vfrac_V * mu_V

    # Compressibility factor of the two-phase mixture
    M = AS.molar_mass()
    R = AS.gas_constant()
    Z_mix = p_mix / (rho_mix * (R / M) * T_mix)

    # Speed of sound of the two-phase mixture
    B1 = vfrac_L / (rho_L * a_L**2) + vfrac_V / (rho_V * a_V**2)
    B2 = vfrac_L * rho_L / cp_L * dsdp_L**2 + vfrac_V * rho_V / cp_V * dsdp_V**2
    compressibility_HEM = B1 + T_mix * B2
    if mfrac_V < 1e-6:  # Avoid discontinuity when Q_v=0
        a_HEM = a_L
    elif mfrac_V > 1.0 - 1e-6:  # Avoid discontinuity when Q_v=1
        a_HEM = a_V
    else:
        a_HEM = (1 / rho_mix / compressibility_HEM) ** 0.5

    # Store properties in dictionary
    props = {}
    props["T"] = T_mix
    props["p"] = p_mix
    props["rhomass"] = rho_mix
    props["umass"] = u_mix
    props["hmass"] = h_mix
    props["smass"] = s_mix
    props["gibbsmass"] = gibbs_mix
    props["cvmass"] = cv_mix
    props["cpmass"] = cp_mix
    props["gamma"] = props["cpmass"] / props["cvmass"]
    props["compressibility_factor"] = Z_mix
    props["speed_sound"] = a_HEM
    props["isentropic_bulk_modulus"] = rho_mix * a_HEM**2
    props["isentropic_compressibility"] = (rho_mix * a_HEM**2) ** -1
    props["isothermal_bulk_modulus"] = np.nan
    props["isothermal_compressibility"] = np.nan
    props["isobaric_expansion_coefficient"] = np.nan
    props["viscosity"] = mu_mix
    props["conductivity"] = k_mix
    props["Q"] = mfrac_V
    props["quality_mass"] = mfrac_V
    props["quality_volume"] = vfrac_V

    # Calculate departure from saturation properties
    if compute_subcooling:
        props["subcooling"] = calculate_subcooling(AS)

    if compute_superheating:
        props["superheating"] = calculate_superheating(AS)

    # Add properties as aliases
    for key, value in PROPERTY_ALIAS.items():
        props[key] = props[value]

    return props


def calculate_generalized_quality(abstract_state, alpha=10):
    r"""
    Calculate the generalized quality of a fluid, extending the quality calculation beyond the two-phase region if necessary.

    Below the critical temperature, the quality is calculated from the specific volume of the saturated liquid and vapor states.
    Above the critical temperature, an artificial two-phase region is defined around the critical density line to allow for a finite-width quality computation.

    The quality, :math:`Q`, is given by:

    .. math::

        Q = \frac{v - v(T, Q=0)}{v(T, Q=1) - v(T, Q=0)}

    where :math:`v=1/\rho` is the specific volume and :math:`T` is the temperature.

    Additionally, this function applies smoothing limiters to ensure the quality is bounded between -1 and 2.
    These limiters prevent the quality value from taking arbitrarily large values, enhancing stability and robustness of downstream calculation.
    The limiters use the `logsumexp` method for smooth transitions, with a specified alpha value controlling the smoothness.

    Parameters
    ----------
    abstract_state : CoolProp.AbstractState
        CoolProp abstract state of the fluid.

    alpha : float
        Smoothing factor of the quality-calculation limiters.

    Returns
    -------
    float
        The calculated quality of the fluid.
    """
    # Instantiate new abstract state to compute saturation properties without changing the state of the class
    state_bis = CP.AbstractState(abstract_state.backend_name(), abstract_state.name())

    # Extend quality calculation beyond the two-phase region
    # Checking if subcritical using temperature works better than with pressure
    if abstract_state.T() < abstract_state.T_critical():
        # Saturated liquid
        state_bis.update(CP.QT_INPUTS, 0.00, abstract_state.T())
        rho_liq = state_bis.rhomass()

        # Saturated vapor
        state_bis.update(CP.QT_INPUTS, 1.00, abstract_state.T())
        rho_vap = state_bis.rhomass()

    else:
        # For states at or above the critical temperature, the concept of saturation states is not applicable
        # Instead, an artificial two-phase region is created around the pseudo-critical density line (line of critical density)
        # The width of the artificial two-phase region is assumed to increase linearly with temperature

        # Rename properties
        T = abstract_state.T()
        T_crit = abstract_state.T_critical()
        rho_crit = abstract_state.rhomass_critical()

        # Define pseudocritical region
        T_hat = 1.5 * T_crit
        rho_hat_liq = 1.1 * rho_crit
        rho_hat_vap = 0.9 * rho_crit
        rho_liq = rho_crit + (rho_hat_liq - rho_crit) * (T - T_crit) / (T_hat - T_crit)
        rho_vap = rho_crit + (rho_hat_vap - rho_crit) * (T - T_crit) / (T_hat - T_crit)

    # Compute quality according to definition
    rho = abstract_state.rhomass()
    quality = (1 / rho - 1 / rho_liq) / (1 / rho_vap - 1 / rho_liq + 1e-6)

    # Apply smoothing limiters so that the quality is bounded between [-1, 2]
    # The margin is defined as delta_Q=1 to the left of Q=0 and to the right of Q=1
    quality = math.smooth_minimum(quality, +2, method="logsumexp", alpha=alpha)
    quality = math.smooth_maximum(quality, -1, method="logsumexp", alpha=alpha)

    return quality


def calculate_superheating(abstract_state):
    r"""
    Calculate the degree of superheating for a given CoolProp abstract state.

    This function computes the superheating of a fluid using the CoolProp library's
    abstract state class. It handles both subcritical and supercritical conditions
    to provide a measure of superheating for any thermodynamic state. This results in
    a continuous variation of superheating in the two-phase region, which is necessary
    to achieve in reliable convergence of systems of equations and optimization problems
    involving the degree of superheating.

    Calculation cases:
        - Physically meaningful superheating for subcritical states in the vapor region:

        .. math::

          \Delta T = T - T(p, Q=1) \quad \text{for} \quad h > h(p, Q=1)

        - Artificial superheating for subcritical states in the liquid and two-phase regions:

        .. math::

          \Delta T = \frac{h - h(p, Q=1)}{c_p(p, Q=1)}

        - Artificial superheating for supercritical states defined using the critical density line:

        .. math::

          \Delta T = T - T(p, \rho_{\text{crit}})

    Parameters
    ----------
    abstract_state : CoolProp.AbstractState
        The abstract state of the fluid for which the superheating is to be calculated.

    Returns
    -------
    float
        The degree of superheating in Kelvin.

    Examples
    --------
    >>> import CoolProp as CP
    >>> abstract_state = CP.AbstractState("HEOS", "water")
    >>> abstract_state.update(CP.PT_INPUTS, 101325, 120 + 273.15)
    >>> superheating = calculate_superheating(abstract_state)
    >>> print(f"Superheating is {superheating:+0.3f} K")
    Superheating is +20.026 K

    >>> abstract_state = CP.AbstractState("HEOS", "water")
    >>> abstract_state.update(CP.PQ_INPUTS, 101325, 0.95)
    >>> superheating = calculate_superheating(abstract_state)
    >>> print(f"Superheating is {superheating:+0.3f} K")
    Superheating is -54.244 K
    """
    # Instantiate new abstract state to compute saturation properties without changing the state of the class
    AS = abstract_state
    sat_state = CP.AbstractState(AS.backend_name(), AS.name())

    # Compute triple pressure
    sat_state.update(CP.QT_INPUTS, 1.00, AS.Ttriple())
    p_triple = sat_state.p()

    # Check if the pressure is below the critical pressure of the fluid
    if AS.p() < AS.p_critical():

        # Compute triple pressure (needed to avoid error at low pressure)
        sat_state.update(CP.QT_INPUTS, 1.00, AS.Ttriple())
        p_triple = sat_state.p()

        # Set the saturation state of the fluid at the given pressure
        sat_state.update(CP.PQ_INPUTS, max(p_triple, AS.p()), 1.00)

        # Check if the fluid is in the two-phase or liquid regions
        if AS.hmass() < sat_state.hmass():
            # Below the vapor saturation enthalpy, define superheating as the normalized difference in enthalpy
            # The normalization is done using the specific heat capacity at saturation (cp)
            # This provides a continuous measure of superheating, even in the two-phase region
            superheating = (AS.hmass() - sat_state.hmass()) / sat_state.cpmass()
        else:
            # Outside the two-phase region, superheating is the difference in temperature
            # from the saturation temperature at the same pressure
            superheating = AS.T() - sat_state.T()
    else:
        # For states at or above the critical pressure, the concept of saturation temperature is not applicable
        # Instead, use a 'pseudo-critical' state for comparison, where the density is set to the critical density
        # but the pressure is the same as the state of interest
        rho_crit = AS.rhomass_critical()
        sat_state.update(CP.DmassP_INPUTS, rho_crit, AS.p())

        # Define superheating as the difference in enthalpy from this 'pseudo-critical' state
        # This approach extends the definition of superheating to conditions above the critical pressure
        superheating = AS.T() - sat_state.T()

    return superheating


def calculate_subcooling(abstract_state):
    r"""
    Calculate the degree of subcooling for a given CoolProp abstract state.

    This function computes the subcooling of a fluid using the CoolProp library's
    abstract state class. It handles both subcritical and supercritical conditions
    to provide a measure of subcooling for any thermodynamic state. This results in
    a continuous variation of subcooling in the two-phase region, which is necessary
    to achieve reliable convergence of systems of equations and optimization problems
    involving the degree of subcooling.

    Calculation cases:
        - Physically meaningful subcooling for subcritical states in the liquid region:

        .. math::

          \Delta T = T(p, Q=0) - T \quad \text{for} \quad h < h(p, Q=0)

        - Artificial subcooling for subcritical states in the vapor and two-phase regions:

        .. math::

          \Delta T = \frac{h(p, Q=0) - h}{c_p(p, Q=0)}

        - Artificial subcooling for supercritical states defined using the critical density line:

        .. math::

          \Delta T = T(p, \rho_{\text{crit}}) - T

    Parameters
    ----------
    abstract_state : CoolProp.AbstractState
        The abstract state of the fluid for which the subcooling is to be calculated.

    Returns
    -------
    float
        The degree of subcooling in Kelvin.

    Examples
    --------
    >>> import CoolProp as CP
    >>> abstract_state = CP.AbstractState("HEOS", "water")
    >>> abstract_state.update(CP.PT_INPUTS, 101325, 25+273.15)
    >>> subcooling = bpy.calculate_subcooling(abstract_state)
    >>> print(f"Subcooling is {subcooling:+0.3f} K")
    Subcooling is +74.974 K

    >>> abstract_state = CP.AbstractState("HEOS", "water")
    >>> abstract_state.update(CP.PQ_INPUTS, 101325, 0.05)
    >>> subcooling = bpy.calculate_subcooling(abstract_state)
    >>> print(f"Subcooling is {subcooling:+0.3f} K")
    Subcooling is -26.763 K
    """

    # Instantiate new abstract state to compute saturation properties without changing the state of the class
    AS = abstract_state
    sat_state = CP.AbstractState(AS.backend_name(), AS.name())

    # Check if the pressure is below the critical pressure of the fluid
    if AS.p() < AS.p_critical():

        # Compute triple pressure (needed to avoid error at low pressure)
        sat_state.update(CP.QT_INPUTS, 0.00, AS.Ttriple())
        p_triple = sat_state.p()

        # Set the saturation state of the fluid at the given pressure
        sat_state.update(CP.PQ_INPUTS, max(p_triple, AS.p()), 0.00)

        # Check if the fluid is in the two-phase or vapor regions
        if AS.hmass() > sat_state.hmass():
            # Above the liquid saturation enthalpy, define superheating as the normalized difference in enthalpy
            # The normalization is done using the specific heat capacity at saturation (cp)
            # This provides a continuous measure of superheating, even in the two-phase region
            subcooling = (sat_state.hmass() - AS.hmass()) / sat_state.cpmass()
        else:
            # Outside the two-phase region, superheating is the difference in temperature
            # from the saturation temperature at the same pressure
            subcooling = sat_state.T() - AS.T()
    else:
        # For states at or above the critical pressure, the concept of saturation temperature is not applicable
        # Instead, use a 'pseudo-critical' state for comparison, where the density is set to the critical density
        # but the pressure is the same as the state of interest
        rho_crit = AS.rhomass_critical()
        sat_state.update(CP.DmassP_INPUTS, rho_crit, AS.p())

        # Define superheating as the difference in enthalpy from this 'pseudo-critical' state
        # This approach extends the definition of superheating to conditions above the critical pressure
        subcooling = sat_state.T() - AS.T()

    return subcooling


# ------------------------------------------------------------------------------------ #
# Metastable property calculations using Helmholtz equations of state
# ------------------------------------------------------------------------------------ #
def compute_properties_metastable(
    abstract_state,
    prop_1,
    prop_1_value,
    prop_2,
    prop_2_value,
    rho_guess,
    T_guess,
    solver_algorithm="hybr",
    tolerance=1e-6,
    print_convergence=False,
):
    r"""
    Determine the thermodynamic state for the given fluid property pair by iterating on the 
    density-temperature native inputs of the Helmholtz energy equations of state.

    This function uses a non-linear root finder to determine the solution of the nonlinear system given by:

    .. math::

        R_1(\rho,\, T) = y_1 - y_1(\rho,\, T) = 0 \\
        R_2(\rho,\, T) = y_2 - y_2(\rho,\, T) = 0

    where :math:`(y_1,\, y_2)` is the given thermodynamic property pair (e.g., enthalpy and pressure),
    while density and temperature :math:`(\rho,\, T)` are the independent variables that the solver 
    iterates until the residual of the problem is driven to zero. The calculations :math:`y_1(\rho,\, T)` and
    :math:`y_1(\rho,\, T)` are performed by evaluating the Helmholtz energy equation of state.

    See Also
    --------
    :func:`compute_properties_metastable_rhoT` : Evaluation of Helmholtz energy equation of state.

    Parameters
    ----------
    prop_1 : str
        The the type of the first thermodynamic property.
    prop_1_value : float
        The the numerical value of the first thermodynamic property.
    prop_2 : str
        The the type of the second thermodynamic property.
    prop_2_value : float
        The the numerical value of the second thermodynamic property.
    rho_guess : float
        Initial guess for density
    T_guess : float
        Initial guess for temperature
    method : str, optional
        Method to be used for solving the nonlinear system. Available solvers are:

        - :code:`hybr`:  Uses MINPACK's 'hybrd' method, which is is a modification of Powell's hybrid method (default).
        - :code:`lm`: The Levenberg-Marquardt method, which blends the steepest descent and the Gauss-Newton methods.

    tolerance : float, optional
        Tolerance for the solver termination. Defaults to 1e-6.
    print_convergence : bool, optional
        If True, displays the convergence progress. Defaults to False.

    Returns
    -------
    dict       
        Thermodynamic properties corresponding to the given input pair.
    """

    # Define problem (find root of temperature-entropy residual)
    AS = abstract_state
    problem = _HelmholtzFlashCalculation(prop_1, prop_1_value, prop_2, prop_2_value, AS)

    # Define root-finding solver object
    solver = psv.NonlinearSystemSolver(
        problem,
        method=solver_algorithm,
        tolerance=tolerance,
        max_iterations=100,
        print_convergence=print_convergence,
        update_on="function",
    )

    # Define initial guess and solve the problem
    x0 = np.asarray([rho_guess, T_guess])
    rho, T = solver.solve(x0)
    props = compute_properties_metastable_rhoT(AS, rho, T)

    # Check if solver converged
    if not solver.success:
        raise ValueError(
            f"Metastable property calculation did not converge.\n{solver.message}"
        )

    return props


def compute_properties_metastable_rhoT(abstract_state, rho, T):
    r"""
    Compute the thermodynamic properties of a fluid using temperature-density calls to the Helmholtz energy equation of state (HEOS).

    Parameters
    ----------
    abstract_state : CoolProp.AbstractState
        The abstract state of the fluid for which the properties are to be calculated.
    rho : float
        Density of the fluid (kg/m³).
    T : float
        Temperature of the fluid (K).

    Returns
    -------
    dict
        Thermodynamic properties computed at the given density and temperature.

    Notes
    -----
    The Helmholtz energy equation of state (HEOS) expresses the Helmholtz energy as an explicit function
    of temperature and density:

    .. math::
        \Phi = \Phi(\rho, T)

    In dimensionless form, the Helmholtz energy is given by:

    .. math::
        \phi(\delta, \tau) = \frac{\Phi(\delta, \tau)}{R T}

    where:

    - :math:`\phi` is the dimensionless Helmholtz energy
    - :math:`R` is the gas constant of the fluid
    - :math:`\delta = \rho / \rho_c` is the reduced density
    - :math:`\tau = T_c / T` is the inverse of the reduced temperature

    Thermodynamic properties can be derived from the Helmholtz energy and its partial derivatives.
    The following table summarizes the expressions for various properties:

    .. list-table:: Helmholtz energy thermodynamic relations
        :widths: 30 70
        :header-rows: 1

        * - Property name
          - Mathematical relation
        * - Pressure
          - .. math:: Z = \frac{p}{\rho R T} = \delta \cdot \phi_{\delta}
        * - Entropy
          - .. math:: \frac{s}{R} = \tau \cdot \phi_{\tau} - \phi
        * - Internal energy
          - .. math:: \frac{u}{R T} = \tau \cdot \phi_{\tau}
        * - Enthalpy
          - .. math:: \frac{h}{R T} = \tau \cdot \phi_{\tau} + \delta \cdot \phi_{\delta}
        * - Isochoric heat capacity
          - .. math:: \frac{c_v}{R} = -\tau^2 \cdot \phi_{\tau \tau}
        * - Isobaric heat capacity
          - .. math:: \frac{c_p}{R} = -\tau^2 \cdot \phi_{\tau \tau} + \frac{(\delta \cdot \phi_{\delta} - \tau \cdot \delta \cdot \phi_{\tau \delta})^2}{2 \cdot \delta \cdot \phi_{\delta} + \delta^2 \cdot \phi_{\delta \delta}}
        * - Isobaric expansivity
          - .. math:: \alpha \cdot T = \frac{\delta \cdot \phi_{\delta} - \tau \cdot \delta \cdot \phi_{\tau \delta}}{2 \cdot \delta \cdot \phi_{\delta} + \delta^2 \cdot \phi_{\delta \delta}}
        * - Isothermal compressibility
          - .. math:: \rho R T \beta = \left(2 \cdot \delta \cdot \phi_{\delta} + \delta^2 \cdot \phi_{\delta \delta} \right)^{-1}
        * - Isothermal bulk modulus
          - .. math:: \frac{K_T}{\rho R T} = 2 \cdot \delta \cdot \phi_{\delta} + \delta^2 \cdot \phi_{\delta \delta}
        * - Isentropic bulk modulus
          - .. math:: \frac{K_s}{\rho R T} = 2 \cdot \delta \cdot \phi_{\delta} + \delta^2 \ \cdot \phi_{\delta \delta} - \frac{(\delta \cdot \phi_{\delta} - \tau \cdot \delta \cdot \phi_{\tau \delta})^2}{\tau^2 \cdot \phi_{\tau \tau}}
        * - Joule-Thompson coefficient
          - .. math:: \rho R \mu_{\mathrm{JT}} = - \frac{\delta \cdot \phi_{\delta} + \tau \cdot \delta \cdot \phi_{\tau \delta} + \delta^2 \cdot \phi_{\delta \delta}}{(\delta \cdot \phi_{\delta} - \tau \cdot \delta \cdot \phi_{\tau \delta})^2 - \tau^2 \cdot \phi_{\tau \tau} \cdot (2 \cdot \delta \cdot \phi_{\delta} + \delta^2 \cdot \phi_{\delta \delta})}

    Where the following abbreviations are used:

    - :math:`\phi_{\delta} = \left(\frac{\partial \phi}{\partial \delta}\right)_{\tau}`
    - :math:`\phi_{\tau} = \left(\frac{\partial \phi}{\partial \tau}\right)_{\delta}`
    - :math:`\phi_{\delta \delta} = \left(\frac{\partial^2 \phi}{\partial \delta^2}\right)_{\tau, \tau}`
    - :math:`\phi_{\tau \tau} = \left(\frac{\partial^2 \phi}{\partial \tau^2}\right)_{\delta, \delta}`
    - :math:`\phi_{\tau \delta} = \left(\frac{\partial^2 \phi}{\partial \tau \delta}\right)_{\delta, \tau}`

    This function can be used to estimate metastable properties using the equation of state beyond the saturation lines.
    """

    # Update thermodynamic state
    AS = abstract_state
    AS = CP.AbstractState(AS.backend_name(), AS.name())
    AS.update(CP.DmassT_INPUTS, rho, T)

    # Get fluid constant properties
    R = AS.gas_constant()
    M = AS.molar_mass()
    T_crit = AS.T_critical()
    rho_crit = AS.rhomass_critical()

    # Compute reduced variables
    tau = T_crit / T
    delta = rho / rho_crit

    # Compute from the Helmholtz energy derivatives
    alpha = AS.alpha0() + AS.alphar()
    dalpha_dTau = AS.dalpha0_dTau() + AS.dalphar_dTau()
    dalpha_dDelta = AS.dalpha0_dDelta() + AS.dalphar_dDelta()
    d2alpha_dTau2 = AS.d2alpha0_dTau2() + AS.d2alphar_dTau2()
    d2alpha_dDelta2 = AS.d2alpha0_dDelta2() + AS.d2alphar_dDelta2()
    d2alpha_dDelta_dTau = AS.d2alpha0_dDelta_dTau() + AS.d2alphar_dDelta_dTau()

    # Compute thermodynamic properties from Helmholtz energy EOS
    props = {}
    props["T"] = T
    props["p"] = (R / M) * T * rho * delta * dalpha_dDelta
    props["rhomass"] = rho
    props["umass"] = (R / M) * T * (tau * dalpha_dTau)
    props["hmass"] = (R / M) * T * (tau * dalpha_dTau + delta * dalpha_dDelta)
    props["smass"] = (R / M) * (tau * dalpha_dTau - alpha)
    props["gibbsmass"] = (R / M) * T * (alpha + delta * dalpha_dDelta)
    props["cvmass"] = (R / M) * (-(tau**2) * d2alpha_dTau2)
    props["cpmass"] = (R / M) * (
        -(tau**2) * d2alpha_dTau2
        + (delta * dalpha_dDelta - delta * tau * d2alpha_dDelta_dTau) ** 2
        / (2 * delta * dalpha_dDelta + delta**2 * d2alpha_dDelta2)
    )
    props["gamma"] = props["cpmass"] / props["cvmass"]
    props["compressibility_factor"] = delta * dalpha_dDelta
    a_square = (R / M * T) * (
        (2 * delta * dalpha_dDelta + delta**2 * d2alpha_dDelta2)
        - (delta * dalpha_dDelta - delta * tau * d2alpha_dDelta_dTau) ** 2
        / (tau**2 * d2alpha_dTau2)
    )
    props["speed_sound"] = np.sqrt(a_square) if a_square > 0 else np.nan
    props["isentropic_bulk_modulus"] = (rho * R / M * T) * (
        (2 * delta * dalpha_dDelta + delta**2 * d2alpha_dDelta2)
        - (delta * dalpha_dDelta - delta * tau * d2alpha_dDelta_dTau) ** 2
        / (tau**2 * d2alpha_dTau2)
    )
    props["isentropic_compressibility"] = 1 / props["isentropic_bulk_modulus"]
    props["isothermal_bulk_modulus"] = (
        R / M * T * rho * (2 * delta * dalpha_dDelta + delta**2 * d2alpha_dDelta2)
    )
    props["isothermal_compressibility"] = 1 / (
        R / M * T * rho * (2 * delta * dalpha_dDelta + delta**2 * d2alpha_dDelta2)
    )
    props["isobaric_expansion_coefficient"] = (
        (1 / T)
        * (delta * dalpha_dDelta - delta * tau * d2alpha_dDelta_dTau)
        / (2 * delta * dalpha_dDelta + delta**2 * d2alpha_dDelta2)
    )
    props["viscosity"] = AS.viscosity()
    props["conductivity"] = AS.conductivity()
    props["Q"] = np.nan
    props["quality_mass"] = np.nan
    props["quality_volume"] = np.nan

    # Add properties as aliases
    for key, value in PROPERTY_ALIAS.items():
        props[key] = props[value]

    return props


class _HelmholtzFlashCalculation(psv.NonlinearSystemProblem):
    """Auxiliary class for property evaluations using Helmholtz energy EoS"""

    def __init__(self, prop_1, prop_1_value, prop_2, prop_2_value, AS):
        self.prop_1 = prop_1
        self.prop_2 = prop_2
        self.prop_1_value = prop_1_value
        self.prop_2_value = prop_2_value
        self.AS = AS

    def residual(self, x):
        """
        Compute the residuals for the given density and temperature.

        Parameters
        ----------
        x : list
            List containing the values for density and temperature.

        Returns
        -------
        np.ndarray
            Array containing residuals (difference) for the two properties.
        """
        # Ensure x can be indexed and contains exactly two elements
        if not hasattr(x, "__getitem__") or len(x) != 2:
            msg = f"Input x={x} must be a list, tuple or numpy array containing exactly two elements: density and temperature."
            raise ValueError(msg)
        rho, T = x
        state = compute_properties_metastable_rhoT(self.AS, rho, T)
        res_1 = 1 - state[self.prop_1] / self.prop_1_value
        res_2 = 1 - state[self.prop_2] / self.prop_2_value
        return np.asarray([res_1, res_2])


# ------------------------------------------------------------------------------------ #
# Spinodal point calculations
# ------------------------------------------------------------------------------------ #


def compute_spinodal_point(
    T,
    abstract_state,
    branch,
    rho_guess=None,
    N_trial=50,
    tolerance=1e-6,
    print_convergence=False,
):
    r"""
    Compute the vapor or liquid spinodal point of a fluid at a given temperature.

    Parameters
    ----------
    T : float
        Temperature of the fluid (K).
    abstract_state : CoolProp.AbstractState
        The abstract state of the fluid for which the spinodal point is to be calculated.
    branch : str
        Branch of the spinodal line used to determine the density initial guess.
        Options: 'liquid' or 'vapor'.
    rho_guess : float, optional
        Initial guess for the density. If provided, this value will be used directly.
        If not provided, the density initial guess will be generated based on a number of trial points.
    N_trial : int, optional
        Number of trial points to generate the density initial guess. Default is 50.
    tolerance : float, optional
        Tolerance for the solver termination. Defaults to 1e-6.
    print_convergence : bool, optional
        If True, displays the convergence progress. Defaults to False.

    Returns
    -------
    dict
        Thermodynamic properties at the spinodal point.

    Raises
    ------
    ValueError
        If the input temperature is higher than the critical temperature or lower than
        the triple temperature.

    Notes
    ------
    When a single-phase fluid undergoes a thermodynamic process and enters the two-phase region it
    can exist in a single-phase state that is different from the equilibrium two-phase state.
    Such states are know as metastable states and they are only possible in the thermodynamic
    region between the saturation lines and the spinodal lines. If the thermodynamic process
    continues and crosses the spinodal lines metastable states become unstable and the transition
    to two-distinct phase occurs rapidly. Therefore, the spinodal line represents the locus of points
    that separates the region where a mixture is thermodynamically unstable and prone to phase separation
    from the region where metastable states are physically possible.

    In mathematical terms, the spinodal line is defined as the loci of thermodynamic states in which the isothermal bulk modulus of the fluid is zero:

    .. math::

        K_T = \rho \left( \frac{\partial p}{\partial \rho} \right)_T = 0

    More precisely, a vapor spinodal point is the first local maximum of a isotherm line in a pressure-density diagram as the density increases.
    Conversely, a liquid spinodal point is the first local minimum of a isotherm line in a pressure-density diagram as the density decreases.
    The spinodal lines and phase envelope of carbon dioxide according to the HEOS developed by :cite:`span_new_1996` are illustrated in the figure below

    .. image:: /_static/spinodal_points_CO2.svg
        :alt: Pressure-density diagram and spinodal points for carbon dioxide.


    Some equations of state are not well-posed and do not satisfy the condition :math:`K_T=0` within the two-phase region.
    This is exemplified by the nitrogen HEOS developed by :cite:`span_reference_2000`.

    .. image:: /_static/spinodal_points_nitrogen.svg
        :alt: Pressure-density diagram and "pseudo" spinodal points for carbon dioxide.

    As seen in the figure, this HEOS is not well-posed because there are isotherms that do not have a local minimum/maximum corresponding to a state with zero isothermal bulk modulus.
    In such cases, this function returns the inflection point of the isotherms (corresponding to the point closest to zero bulk modulus) as the spinodal point.

    """

    # Instantiate new abstract state to compute saturation properties without changing the state of the class
    AS = abstract_state
    AS = CP.AbstractState(AS.backend_name(), AS.name())

    # Check that the inlet temperature is lower than the critical temperature
    T_critical = AS.T_critical()
    if T >= T_critical:
        msg = f"T={T:.3f}K must be less than T_critical={AS.T_critical:.3f}K"
        raise ValueError(msg)

    # Check that the inlet temperature is greater than the triple temperature
    T_triple = AS.Ttriple()
    if T < T_triple:
        msg = f"T={T:.3f}K must be greater than T_triple={T_triple:.3f}K"
        raise ValueError(msg)

    # Create spinodal point optimization problem
    problem = _SpinodalPointProblem(T, abstract_state, branch)

    # Create solver object
    solver = psv.OptimizationSolver(
        problem=problem,
        library="scipy",  
        method="bfgs", # "l-bfgs-b", "slsqp" "bfgs"
        tolerance=tolerance,
        print_convergence=print_convergence,
    )

    # Generate initial guess if not provided
    if rho_guess is None:
        rho_guess = problem.generate_density_guess(N_trial)

    # Solve the problem using the provided or generated initial guess
    rho_opt = solver.solve(rho_guess).item()
    props = compute_properties_metastable_rhoT(AS, rho_opt, T)

    # I tried different combinations of solvers and initial guesses
    # SLSQP is very fast, but also very aggresive and can lead to unpredictable results
    # BFGS is reliable when the initial guess is good
    # Achieving a good initial guess is possible by:
    #  1. Using previous point in the spinodal line and having high resolution on the spinodal line
    #  2. Generating an initial guess for each point with the initial guess strategy
    # Option 2. seems the most reliable, and even if the computational cost can be a bit hight
    # it seems to be the most effective way to get accurate spinodal lines.

    return props


class _SpinodalPointProblem(psv.OptimizationProblem):
    """Auxiliary class for the determination of the liquid and vapor spinodal points"""

    def __init__(self, T, abstract_state, branch):

        # Declare class attributes
        self.T = T
        self.branch = branch
        self.abstract_state = abstract_state

        # Compute saturation liquid density (used to determine initial guess)
        self.abstract_state.update(CP.QT_INPUTS, 0.00, self.T)
        self.rho_liq = self.abstract_state.rhomass()

        # Calculate saturation vapor density (used to determine initial guess)
        self.abstract_state.update(CP.QT_INPUTS, 1.00, self.T)
        self.rho_vap = self.abstract_state.rhomass()

        # Caclulate the critical density (used to determine initial guess)
        self.rho_crit = self.abstract_state.rhomass_critical()

    def generate_density_guess(self, N):
        """Generate a density initial guess that is close to the first local minima of the absolute value of the bulk modulus"""

        # Generate candidate densities between saturation and the critical value
        if self.branch == "liquid":
            rho_array = np.linspace(self.rho_liq, self.rho_crit, N)
        elif self.branch == "vapor":
            rho_array = np.linspace(self.rho_vap, self.rho_crit, N)
        else:
            msg = f"Invalid value for parameter branch={self.branch}. Options: 'liquid' or 'vapor'"
            raise ValueError(msg)

        # Evaluate residual vector at trial densities
        residual = np.abs([self.fitness(rho) for rho in rho_array])

        # Return the first local minima in the residual vector
        for i in range(1, N - 1):
            if residual[i - 1] > residual[i] < residual[i + 1]:
                self.rho_guess = rho_array[i + 1]
                return self.rho_guess

    def fitness(self, rho):
        """
        Compute the objective function of the optimization problem: the absolute value of the isothermal bulk modulus.

        This function uses the absolute value of residual to solve an optimization problem
        rather than a non-linear equation having the bulk modulus as residual. This approach
        is adopted because some fluids (e.g., nitrogen) have ill-posed EoS that not have a well-defined
        spinodal points where the isothermal bulk modulus is zero.

        Not a good idea to scale the bulk modulus by pressure because it can take negative values or zero when evaluated with the HEOS.
        """
        props = compute_properties_metastable_rhoT(self.abstract_state, rho, self.T)
        return np.atleast_1d(np.abs(props["isothermal_bulk_modulus"]))

    def get_bounds(self):
        """Compute the bounds of the optimization problem."""
        # if self.branch == "liquid":
        #     return [(self.rho_crit,), (self.rho_liq,)]
        # elif self.branch == "vapor":
        #     return [(self.rho_vap,), (self.rho_crit,)]
        # else:
        #     raise ValueError(f"Invalid value for branch={self.branch}")
        return None
    
    def get_nec(self):
        return 0

    def get_nic(self):
        return 0



