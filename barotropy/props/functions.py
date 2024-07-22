import copy
import scipy
import numpy as np
import CoolProp.CoolProp as CP

from .. import pysolver_view as solver

# ------------------------------------------------------------------------------------ #
# Equilibrum property calculations through CoolProp
# ------------------------------------------------------------------------------------ #


def compute_properties_1phase(AS, generalize_quality=True):
    """Get single-phase properties from CoolProp low level interface"""

    # Fluid properties available from CoolProp
    props = {}
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

    return props


def compute_properties_2phase(AS):
    """Get two-phase properties from mixing rules and single-phase CoolProp properties"""

    # Instantiate new AbstractState to compute saturation properties without changing the state of the class
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
    speed_sound_L = cloned_AS.speed_sound()
    dsdp_L = cloned_AS.first_saturation_deriv(CP.iSmass, CP.iP)

    # Saturated vapor properties
    cloned_AS.update(CP.QT_INPUTS, 1.00, T_mix)
    rho_V = cloned_AS.rhomass()
    cp_V = cloned_AS.cpmass()
    cv_V = cloned_AS.cvmass()
    k_V = cloned_AS.conductivity()
    mu_V = cloned_AS.viscosity()
    speed_sound_V = cloned_AS.speed_sound()
    dsdp_V = cloned_AS.first_saturation_deriv(CP.iSmass, CP.iP)

    # Volume fractions of vapor and liquid
    vol_frac_V = (rho_mix - rho_L) / (rho_V - rho_L)
    vol_frac_L = 1.00 - vol_frac_V

    # Mass fractions of vapor and liquid
    mass_frac_V = (1 / rho_mix - 1 / rho_L) / (1 / rho_V - 1 / rho_L)
    mass_frac_L = 1.00 - mass_frac_V

    # Heat capacities of the two-phase mixture
    cp_mix = mass_frac_L * cp_L + mass_frac_V * cp_V
    cv_mix = mass_frac_L * cv_L + mass_frac_V * cv_V

    # Transport properties of the two-phase mixture
    k_mix = vol_frac_L * k_L + vol_frac_V * k_V
    mu_mix = vol_frac_L * mu_L + vol_frac_V * mu_V

    # Compressibility factor of the two-phase mixture
    M = AS.molar_mass()
    R = AS.gas_constant()
    Z_mix = p_mix / (rho_mix * (R / M) * T_mix)

    # Speed of sound of the two-phase mixture
    mechanical_equilibrium = vol_frac_L / (rho_L * speed_sound_L**2) + vol_frac_V / (
        rho_V * speed_sound_V**2
    )
    thermal_equilibrium = T_mix * (
        vol_frac_L * rho_L / cp_L * dsdp_L**2
        + vol_frac_V * rho_V / cp_V * dsdp_V**2
    )
    compressibility_HEM = mechanical_equilibrium + thermal_equilibrium
    if mass_frac_V < 1e-6:  # Avoid discontinuity when Q_v=0
        a_HEM = speed_sound_L
    elif mass_frac_V > 1.0 - 1e-6:  # Avoid discontinuity when Q_v=1
        a_HEM = speed_sound_V
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
    props["Q"] = mass_frac_V
    props["quality_mass"] = mass_frac_V
    props["quality_volume"] = vol_frac_V
    return props


def calculate_generalized_quality(abstract_state):
    """
    Calculates the generalized quality of a fluid, extending the quality calculation
    beyond the two-phase region if necessary.

    Parameters:
    - backend: The backend used by CoolProp for property calculations.
    - name: The name of the fluid.
    - pressure: The pressure of the fluid state.
    - critical_point: A dictionary containing the critical point properties, including 'p' (pressure) and 'rho' (density).

    Returns:
    - quality: The calculated quality of the fluid. Returns np.nan if the quality is not generalized.
    """

    # Instantiate new fluid object to compute saturation properties without changing the state of the class
    cloned_state = CP.AbstractState(
        abstract_state.backend_name(), abstract_state.name()
    )

    # Extend quality calculation beyond the two-phase region
    if abstract_state.p() < abstract_state.p_critical():
        # Saturated liquid
        cloned_state.update(CP.PQ_INPUTS, abstract_state.p(), 0.00)
        h_liq = cloned_state.hmass()

        # Saturated vapor
        cloned_state.update(CP.PQ_INPUTS, abstract_state.p(), 1.00)
        h_vap = cloned_state.hmass()

        # Generalized quality
        quality = (abstract_state.hmass() - h_liq) / (h_vap - h_liq)
    else:
        # For states at or above the critical pressure, the concept of saturation states is not applicable
        # Instead, use a 'pseudo-critical' state for comparison, where the density is set to the critical density
        # but the pressure is the same as the state of interest
        # Use a band of a certain width to prevent a discontinuity
        cloned_state.update(
            CP.DmassP_INPUTS, abstract_state.rhomass_critical(), abstract_state.p()
        )
        quality = (abstract_state.hmass() - 0.99 * cloned_state.hmass()) / (
            1.01 * cloned_state.hmass() - 0.99 * cloned_state.hmass()
        )

    return quality


# ------------------------------------------------------------------------------------ #
# Metastable property calculations using Helmholtz equations of state
# ------------------------------------------------------------------------------------ #
def compute_properties_metastable_rhoT(rho, T, AS):
    """
    Compute the thermodynamic properties of a fluid using the Helmholtz
    energy equation of state. All properties thermodynamic properties can
    be derived as combinations of the Helmholtz energy and its
    derivatives with respect to density and pressure.

    This function can be used to estimate metastable properties using the
    equation of state beyond the saturation lines.
    """

    # Update thermodynamic state
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
    props["speed_sound"] = (
        (R / M * T)
        * (
            (2 * delta * dalpha_dDelta + delta**2 * d2alpha_dDelta2)
            - (delta * dalpha_dDelta - delta * tau * d2alpha_dDelta_dTau) ** 2
            / (tau**2 * d2alpha_dTau2)
        )
    ) ** 0.5
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

    return props


class HelmholtzResidual(solver.NonlinearSystemProblem):
    """
    Find the root for thermodynamic state by iterating on the density-temperature
    native inputs to the helmholtz energy equations of state.

    Attributes
    ----------
    prop_1 : str
        The first property to be compared.
    prop_1_value : float
        The value of the first property.
    prop_2 : str
        The second property to be compared.
    prop_2_value : float
        The value of the second property.

    Methods
    -------
    get_values(x)
        Calculates the residuals based on the given input values of density and temperature.
    """

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
            raise ValueError(
                f"Input x={x} must be a list, tuple or numpy array containing exactly two elements: density and temperature."
            )

        rho, T = x
        state = compute_properties_metastable_rhoT(rho, T, self.AS)

        residual = np.asarray(
            [
                (1 - state[self.prop_1] / self.prop_1_value),
                (1 - state[self.prop_2] / self.prop_2_value),
            ]
        )

        return residual
