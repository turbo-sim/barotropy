import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
import copy

from .PySolverView import NonlinearSystemSolver, NonlinearSystemProblem
      


class Fluid:
    """
    Represents a fluid with various thermodynamic properties computed via CoolProp.

    Attributes
    ----------
    fluid_name : str
        Name of the fluid.
    backend : str
        Backend used for CoolProp, default is 'HEOS'.
    throw_exceptions : bool
        Determines if exceptions should be raised during state calculations. Default is True.
    converged_flag : bool
        Flag indicating whether properties calculations converged.
    properties : dict
        Dictionary of various fluid properties.
    
    Methods
    -------
    compute_critical_properties():
        Calculate the critical state of the fluid.
    compute_triple_properties():
        Calculate the triple state (liquid) of the fluid.
    set_state(input_type, prop_1, prop_2):
        Set the thermodynamic state of the fluid.
    compute_properties_1phase():
        Get single-phase properties from CoolProp.
    compute_properties_2phase():
        Get two-phase properties from mixing rules.
    set_state_metastable():
        Set the metastable state of the fluid.
    set_state_metastable_rhoT(rho, T):
        Set the metastable state using rho and T.

    ...

    Class Attributes
    ----------------
    aliases : dict
        Dictionary mapping common property names to their CoolProp equivalent.
    """

    aliases = {"P": "p",
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
               "k": "conductivity"}
            
    
    def __init__(self, fluid_name, backend="HEOS", throw_exceptions=True):
        self.fluid_name = fluid_name
        self.backend = backend
        self.abstractstate = CP.AbstractState(backend, fluid_name)
        self.throw_exceptions = throw_exceptions
        self.converged_flag = False
        self.properties = {}
        self.compute_critical_properties()
        self.compute_triple_properties()


    def compute_critical_properties(self):
        """ Calculate the critical state of the fluid"""
        self.T_critical = self.abstractstate.T_critical()
        self.rho_critical = self.abstractstate.rhomass_critical()
        self.p_critical = self.abstractstate.p_critical()
        self.set_state(CP.DmassT_INPUTS, self.rho_critical, self.T_critical)
        self.critical_point = copy.deepcopy(self.properties)


    def compute_triple_properties(self):
        """Caclulate the tiple state (liquid) of the fluid"""
        self.T_triple = self.abstractstate.Ttriple()
        self.set_state(CP.QT_INPUTS, 0.00, self.T_triple)
        self.triple_point_liquid = copy.deepcopy(self.properties)
        self.set_state(CP.QT_INPUTS, 1.00, self.T_triple)
        self.triple_point_vapor = copy.deepcopy(self.properties)
        self.p_triple = self.abstractstate.p()

        
    def set_state(self, input_type, prop_1, prop_2):

        try:

            # Update Coolprop thermodynamic state
            self.abstractstate.update(input_type, prop_1, prop_2)

            # Retrieve single-phase properties
            if self.abstractstate.phase() != CP.iphase_twophase:
                self.properties = self.compute_properties_1phase()
            else:
                self.properties = self.compute_properties_2phase()

            # Add properties as aliases
            for key, value in self.aliases.items():
                self.properties[key] = self.properties[value]

            # No errors computing the properies
            self.converged_flag = True

        # Something went wrong while computing the properties
        except Exception as e:
            self.converged_flag = False
            if self.throw_exceptions:
                raise e
            
        return self.properties
    

    def compute_properties_1phase(self):
        """Get single-phase properties from CoolProp""" 
        
        properties = {}
        properties['T'] = self.abstractstate.T()
        properties['p'] = self.abstractstate.p()
        properties['rhomass'] = self.abstractstate.rhomass()
        properties['umass'] = self.abstractstate.umass()
        properties['hmass'] = self.abstractstate.hmass()
        properties['smass'] = self.abstractstate.smass()
        properties['gibbsmass'] = self.abstractstate.gibbsmass()
        properties['cvmass'] = self.abstractstate.cvmass()
        properties['cpmass'] = self.abstractstate.cpmass()
        properties['gamma'] = properties['cpmass']/properties['cvmass']
        properties['compressibility_factor'] = self.abstractstate.compressibility_factor()
        properties['speed_sound'] = self.abstractstate.speed_sound()
        properties['isentropic_bulk_modulus'] = self.abstractstate.rhomass() * self.abstractstate.speed_sound() ** 2
        properties['isentropic_compressibility'] = 1 / properties["isentropic_bulk_modulus"]
        properties['isothermal_bulk_modulus'] = 1 / self.abstractstate.isothermal_compressibility()
        properties['isothermal_compressibility'] = self.abstractstate.isothermal_compressibility()
        properties['isobaric_expansion_coefficient'] = self.abstractstate.isobaric_expansion_coefficient()
        properties['viscosity'] = self.abstractstate.viscosity()
        properties['conductivity'] = self.abstractstate.conductivity()
        properties['Q'] = np.nan
        properties['quality_mass'] = np.nan
        properties['quality_volume'] = np.nan

        return properties


    def compute_properties_2phase(self):
        """Get two-phase properties from mixing rules""" 
        
        # Basic properties of the two-phase mixture
        T_mix = self.abstractstate.T()
        p_mix = self.abstractstate.p()
        rho_mix = self.abstractstate.rhomass()
        u_mix = self.abstractstate.umass()
        h_mix = self.abstractstate.hmass()
        s_mix = self.abstractstate.smass()
        gibbs_mix = self.abstractstate.gibbsmass()

        # Instantiate new fluid object to compute saturation properties without changing the state of the class
        temp = CP.AbstractState(self.backend, self.fluid_name)

        # Saturated liquid properties
        temp.update(CP.QT_INPUTS, 0.00, T_mix)
        rho_L = temp.rhomass()
        cp_L = temp.cpmass()
        cv_L = temp.cvmass()
        k_L = temp.conductivity()
        mu_L = temp.viscosity()
        speed_sound_L = temp.speed_sound()
        dsdp_L = temp.first_saturation_deriv(CP.iSmass, CP.iP)

        # Saturated vapor properties
        temp.update(CP.QT_INPUTS, 1.00, T_mix)
        rho_V = temp.rhomass()
        cp_V = temp.cpmass()
        cv_V = temp.cvmass()
        k_V = temp.conductivity()
        mu_V = temp.viscosity()
        speed_sound_V = temp.speed_sound()
        dsdp_V = temp.first_saturation_deriv(CP.iSmass, CP.iP)

        # Volume fractions of vapor and liquid
        vol_frac_V = (rho_mix - rho_L) / (rho_V - rho_L)
        vol_frac_L = 1.00 - vol_frac_V

        # Mass fractions of vapor and liquid
        mass_frac_V = (1/rho_mix - 1/rho_L) / (1/rho_V - 1/rho_L)
        mass_frac_L = 1.00 - mass_frac_V

        # Heat capacities of the two-phase mixture
        cp_mix = mass_frac_L*cp_L + mass_frac_V*cp_V
        cv_mix = mass_frac_L*cv_L + mass_frac_V*cv_V

        # Transport properties of the two-phase mixture
        k_mix = vol_frac_L*k_L + vol_frac_V*k_V
        mu_mix = vol_frac_L*mu_L + vol_frac_V*mu_V

        # Compressibility factor of the two-phase mixture
        M = self.abstractstate.molar_mass()
        R = self.abstractstate.gas_constant()
        Z_mix = p_mix / (rho_mix * (R/M) * T_mix)

        # Speed of sound of the two-phase mixture
        mechanical_equilibrium = vol_frac_L/(rho_L*speed_sound_L**2) + vol_frac_V/(rho_V*speed_sound_V**2)
        thermal_equilibrium = T_mix*(vol_frac_L*rho_L/cp_L*dsdp_L**2 + vol_frac_V*rho_V/cp_V*dsdp_V**2)
        compressibility_HEM = mechanical_equilibrium + thermal_equilibrium
        if mass_frac_V < 1e-6:  # Avoid discontinuity when Q_v=0
            a_HEM = speed_sound_L
        elif mass_frac_V > 1.0 - 1e-6:  # Avoid discontinuity when Q_v=1
            a_HEM = speed_sound_V
        else:
            a_HEM = (1/rho_mix/compressibility_HEM)**0.5

        # Store properties in dictionary
        properties = {}
        properties['T'] = T_mix
        properties['p'] = p_mix
        properties['rhomass'] = rho_mix
        properties['umass'] = u_mix
        properties['hmass'] = h_mix
        properties['smass'] = s_mix
        properties['gibssmass'] = gibbs_mix
        properties['cvmass'] = cv_mix
        properties['cpmass'] = cp_mix
        properties['gamma'] = properties['cpmass']/properties['cvmass']
        properties['compressibility_factor'] = Z_mix
        properties['speed_sound'] = a_HEM
        properties['isentropic_bulk_modulus'] = rho_mix * a_HEM**2
        properties['isentropic_compressibility'] = 1 / properties["isentropic_bulk_modulus"]
        properties['isothermal_bulk_modulus'] = np.nan
        properties['isothermal_compressibility'] = np.nan
        properties['isobaric_expansion_coefficient'] = np.nan
        properties['viscosity'] = mu_mix
        properties['conductivity'] = k_mix
        properties['Q'] = mass_frac_V
        properties['quality_mass'] = mass_frac_V
        properties['quality_volume'] = vol_frac_V

        return properties
    


    def set_state_metastable(self, prop_1, prop_1_value, prop_2, prop_2_value, rho_guess, T_guess):




        # problem = PropertyRoot()



        return

    def set_state_metastable_rhoT(self, rho, T):

        # TODO: Add check to see if we are inside thespinodal and return two phase properties if yes
        # TODO: implement root finding functionality to accept p-h, T-s, p-s arguments [good initial guess required]
        # TODO: can it be generalized so that it uses equilibrium as initial guess with any inputs? (even if they are T-d)
        try:

            # Update Coolprop thermodynamic state
            self.properties = self.compute_properties_metastable_rhoT(rho, T, self.abstractstate)

            # Add properties as aliases
            for key, value in self.aliases.items():
                self.properties[key] = self.properties[value]

            # No errors computing the properies
            self.converged_flag = True

        # Something went wrong while computing the properties
        except Exception as e:
            self.converged_flag = False
            if self.throw_exceptions:
                raise e
            
        return self.properties
    





    @staticmethod
    def compute_properties_metastable_rhoT(rho, T, abstract_state):
        """
        Compute the thermodynamic properties of a fluid using the Helmholtz
        energy equation of state. All properties thermodynamic properties can
        be derived as combinations of the Helmholtz energy and its
        derivatives with respect to density and pressure.
        
        This function can be used to estimate metastable properties using the
        equation of state beyond the saturation lines.
        """
        
        # Update thermodynamic state
        abstract_state.update(CP.DmassT_INPUTS, rho, T)
        
        # Get fluid constant properties
        R = abstract_state.gas_constant()
        M = abstract_state.molar_mass()
        T_crit = abstract_state.T_critical()
        rho_crit = abstract_state.rhomass_critical()

        # Compute reduced variables
        tau = T_crit / T
        delta = rho / rho_crit

        # Compute from the Helmholtz energy derivatives
        alpha = abstract_state.alpha0() + abstract_state.alphar()
        dalpha_dTau = abstract_state.dalpha0_dTau() + abstract_state.dalphar_dTau()
        dalpha_dDelta = abstract_state.dalpha0_dDelta() + abstract_state.dalphar_dDelta()
        d2alpha_dTau2 = abstract_state.d2alpha0_dTau2() + abstract_state.d2alphar_dTau2()
        d2alpha_dDelta2 = abstract_state.d2alpha0_dDelta2() + abstract_state.d2alphar_dDelta2()
        d2alpha_dDelta_dTau = abstract_state.d2alpha0_dDelta_dTau() + abstract_state.d2alphar_dDelta_dTau()
        
        # Compute thermodynamic properties from Helmholtz energy EOS
        properties = {}
        properties['T'] = T
        properties['p'] = (R/M)*T*rho*delta*dalpha_dDelta
        properties['rhomass'] = rho
        properties['umass'] = (R/M)*T*(tau*dalpha_dTau)
        properties['hmass'] = (R/M)*T*(tau*dalpha_dTau+delta*dalpha_dDelta)
        properties['smass'] = (R/M)*(tau*dalpha_dTau - alpha)
        properties['gibbsmass'] = (R/M)*T*(alpha + delta*dalpha_dDelta)
        properties['cvmass'] = (R/M)*(-tau**2*d2alpha_dTau2)
        properties['cpmass'] = (R/M)*(-tau**2*d2alpha_dTau2 + (delta*dalpha_dDelta - delta*tau*d2alpha_dDelta_dTau)**2/(2*delta*dalpha_dDelta + delta**2*d2alpha_dDelta2))
        properties['gamma'] = properties['cpmass']/properties['cvmass']
        properties['compressibility_factor'] = delta*dalpha_dDelta
        properties['speed_sound'] = ((R/M)*T*(2*delta*dalpha_dDelta + delta**2*d2alpha_dDelta2 - (delta*dalpha_dDelta - delta*tau*d2alpha_dDelta_dTau)**2/(tau**2*d2alpha_dTau2)))**0.5
        properties['isentropic_bulk_modulus'] = rho*(R/M)*T*(2*delta*dalpha_dDelta + delta**2*d2alpha_dDelta2 - (delta*dalpha_dDelta - delta*tau*d2alpha_dDelta_dTau)**2/(tau**2*d2alpha_dTau2))
        properties['isentropic_compressibility'] = 1 / properties["isentropic_bulk_modulus"]
        properties['isothermal_bulk_modulus'] = (R/M)*T*rho*(2*delta*dalpha_dDelta + delta**2*d2alpha_dDelta2)
        properties['isothermal_compressibility'] = 1/((R/M)*T*rho*(2*delta*dalpha_dDelta + delta**2*d2alpha_dDelta2))
        properties['isobaric_expansion_coefficient'] = 1/T*(delta*dalpha_dDelta - delta*tau*d2alpha_dDelta_dTau)/(2*delta*dalpha_dDelta + delta**2*d2alpha_dDelta2)
        properties['viscosity'] = abstract_state.viscosity()
        properties['conductivity'] = abstract_state.conductivity()
        properties['Q'] = np.nan
        properties['quality_mass'] = np.nan
        properties['quality_volume'] = np.nan
        
        return properties



    def get_property(self, propname):
        """Get the value of a single property"""
        if propname in self.properties:
            return self.properties[propname]
        else:
            valid_options = "\n\t".join(self.properties.keys())
            raise ValueError(f"The requested property '{propname}' is not available. The valid options are:\n\t{valid_options}")

     
    def compute_properties_meanline(self, input_type, prop_1, prop_2):
        """Extract fluid properties for meanline model"""

        # Compute properties in the normal way
        self.set_state(input_type, prop_1, prop_2)
        
        # Store a subset of the properties in a dictionary
        fluid_properties={}
        property_subset = ["p", "T", "h", "s", "d", "Z", "a", "mu", "k", "cp", "cv", "gamma"]
        for item in property_subset:
            fluid_properties[item] = self.properties[item]

        return fluid_properties


def plot_phase_diagram(prop_x, prop_y, fluid, axes=None, N_points=200, plot_saturation_line=True, 
                       plot_critical_point=True, plot_triple_point=False, plot_spinodal_line=False, 
                       spinodal_line_method='standard', spinodal_line_color=0.5 * np.array([1, 1, 1]), 
                       spinodal_line_width=0.75, plot_quality_isolines=False, plot_pseudocritical_line=False, 
                       quality_levels=np.linspace(0.1, 1.0, 10), quality_labels=False, show_in_legend=False):
    
    if axes is None:
        axes = plt.gca()
    
    # Plot saturation line
    if plot_saturation_line:
        sat_liq, sat_vap = compute_saturation_line(fluid, N_points=N_points)
        x = sat_liq[prop_x] + sat_vap[prop_x]
        y = sat_liq[prop_y] + sat_vap[prop_y]
        label_value = "Saturation line" if show_in_legend else "_no_legend_"
        axes.plot(x, y, 'k', linewidth=1.25, label=label_value)
    
    # Plot spinodal line
    if plot_spinodal_line:
        spinodal_liq, spinodal_vap = compute_spinodal_line(fluid, N_points=N_points, method=spinodal_line_method)
        x = spinodal_liq[prop_x] + spinodal_vap[prop_x]
        y = spinodal_liq[prop_y] + spinodal_vap[prop_y]
        label_value = "Spinodal line" if show_in_legend else "_no_legend_"
        axes.plot(x, y, color=spinodal_line_color, linewidth=spinodal_line_width, label=label_value)
    
    # Compute vapor quality isocurves
    if plot_quality_isolines:
        t1 = np.logspace(np.log10(1-0.9999), np.log10(0.1), int(N_points/2))
        t2 = np.logspace(np.log10(0.1), np.log10(1-(fluid.T_triple)/fluid.T_critical), int(N_points/2))
        T_quality = (1-np.hstack((t1, t2)))*fluid.T_critical
        Q_quality = quality_levels
        x_quality, y_quality = np.zeros((len(Q_quality), len(T_quality))), np.zeros((len(Q_quality), len(T_quality)))
        for i, q in enumerate(Q_quality):
            for j, T in enumerate(T_quality):
                fluid.set_state(CP.QT_INPUTS, q, T)
                x_quality[i, j] = fluid.properties[prop_x]
                y_quality[i, j] = fluid.properties[prop_y]
        CS = axes.contour(x_quality, y_quality, np.tile(Q_quality, (len(T_quality), 1)).T, quality_levels, colors='black', linestyles=':', linewidths=0.75)
        if quality_labels:
            axes.clabel(CS, inline=True, fontsize=9, rightside_up=True)

    # Plot pseudocritical line
    if plot_pseudocritical_line:
        pseudo_critical_line = compute_pseudocritical_line(fluid)
        label_value = "Pseudocritical line" if show_in_legend else "_no_legend_"
        axes.plot(pseudo_critical_line[prop_x], pseudo_critical_line[prop_y], color='black', linestyle='-', linewidth=0.25, label=label_value)
    
    # Plot critical point
    if plot_critical_point:
        label_value = "Critical point" if show_in_legend else "_no_legend_"
        axes.plot(fluid.critical_point[prop_x], fluid.critical_point[prop_y], 'ko', markersize=4.5, markerfacecolor='w', label=label_value)

    # Plot triple point
    if plot_triple_point:
        label_value = "Triple point liquid" if show_in_legend else "_no_legend_"
        axes.plot(fluid.triple_point_liquid[prop_x], fluid.triple_point_liquid[prop_y], 'ko', markersize=4.5, markerfacecolor='w', label=label_value)
        label_value = "Triple point vapor" if show_in_legend else "_no_legend_"
        axes.plot(fluid.triple_point_vapor[prop_x], fluid.triple_point_vapor[prop_y], 'ko', markersize=4.5, markerfacecolor='w', label=label_value)
       

def compute_saturation_line(fluid, N_points=100):
    
    # Initialize objects to store properties
    prop_names = fluid.properties.keys()
    liquid_line = {name: [] for name in prop_names}
    vapor_line = {name: [] for name in prop_names}

    # Define temperature array with refinement close to the critical point
    ratio = 1 - fluid.T_triple / fluid.T_critical
    t1 = np.logspace(np.log10(1-0.9999), np.log10(ratio/10), int(np.ceil(N_points/2)))
    t2 = np.logspace(np.log10(ratio/10), np.log10(ratio), int(np.floor(N_points/2)))
    T_sat = (1 - np.concatenate([t1, t2])) * fluid.T_critical

    # Loop over temperatures and property names in an efficient way
    for T in T_sat:

        # Compute liquid saturation line
        for name in prop_names:
            fluid.set_state(CP.QT_INPUTS, 0.00, T)
            liquid_line[name].append(fluid.properties[name])

        # Compute vapor saturation line
        for name in prop_names:
            fluid.set_state(CP.QT_INPUTS, 1.00, T)
            vapor_line[name].append(fluid.properties[name])

    # Add critical point as part of the spinodal line
    for name in prop_names:
        liquid_line[name] = [fluid.critical_point[name]] + liquid_line[name]
        vapor_line[name] = [fluid.critical_point[name]] + vapor_line[name]

    # Re-format for easy concatenation
    for name in prop_names:
        liquid_line[name] = list(reversed(liquid_line[name]))

    return liquid_line, vapor_line



def compute_spinodal_line(fluid, N_points=100, method='standard'):
    raise NotImplementedError("The 'compute_spinodal_line' function has not been implemented yet.")

def compute_pseudocritical_line(fluid):
    raise NotImplementedError("The 'compute_pseudocritical_line' function has not been implemented yet.")



class PropertyRoot(NonlinearSystemProblem):
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
    fluid : object
        An instance of the fluid class which has the helmholtz energy equations of state.

    Methods
    -------
    get_values(x)
        Calculates the residuals based on the given input values of density and temperature.
    """

    def __init__(self, prop_1, prop_1_value, prop_2, prop_2_value, fluid):
        self.prop_1 = prop_1
        self.prop_2 = prop_2
        self.prop_1_value = prop_1_value
        self.prop_2_value = prop_2_value
        self.fluid = fluid

    def get_values(self, x):
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
            raise ValueError("Input x must be a list, tuple or numpy array containing exactly two elements: density and temperature.")

        rho, T = x
        props = self.fluid.set_state_metastable_rhoT(rho, T)

        residual = np.asarray([self.prop_1_value - props[self.prop_1],
                               self.prop_2_value - props[self.prop_2]])
                
        return residual




    
            
# def compute_properties_metastable_rhoT(rho, T, fluid):
#     """
#     Compute the thermodynamic properties of a fluid using the Helmholtz
#     energy equation of state. All properties thermodynamic properties can
#     be derived as combinations of the Helmholtz energy and its
#     derivatives with respect to density and pressure.
    
#     This function can be used to estimate metastable properties using the
#     equation of state beyond the saturation lines.
#     """
    
#     # Update thermodynamic state
#     fluid.update(CP.DmassT_INPUTS, rho, T)
    
#     # Get fluid constant properties
#     R = fluid.gas_constant()
#     M = fluid.molar_mass()
#     T_crit = fluid.T_critical()
#     rho_crit = fluid.rhomass_critical()

#     # Compute reduced variables
#     tau = T_crit / T
#     delta = rho / rho_crit

#     # Compute from the Helmholtz energy derivatives
#     alpha = fluid.alpha0() + fluid.alphar()
#     dalpha_dTau = fluid.dalpha0_dTau() + fluid.dalphar_dTau()
#     dalpha_dDelta = fluid.dalpha0_dDelta() + fluid.dalphar_dDelta()
#     d2alpha_dTau2 = fluid.d2alpha0_dTau2() + fluid.d2alphar_dTau2()
#     d2alpha_dDelta2 = fluid.d2alpha0_dDelta2() + fluid.d2alphar_dDelta2()
#     d2alpha_dDelta_dTau = fluid.d2alpha0_dDelta_dTau() + fluid.d2alphar_dDelta_dTau()
    
#     # Compute thermodynamic properties from Helmholtz energy EOS
#     properties = {}
#     properties['T'] = T
#     properties['p'] = (R/M)*T*rho*delta*dalpha_dDelta
#     properties['rhomass'] = rho
#     properties['umass'] = (R/M)*T*(tau*dalpha_dTau)
#     properties['hmass'] = (R/M)*T*(tau*dalpha_dTau+delta*dalpha_dDelta)
#     properties['smass'] = (R/M)*(tau*dalpha_dTau - alpha)
#     properties['gibbsmass'] = (R/M)*T*(alpha + delta*dalpha_dDelta)
#     properties['cvmass'] = (R/M)*(-tau**2*d2alpha_dTau2)
#     properties['cpmass'] = (R/M)*(-tau**2*d2alpha_dTau2 + (delta*dalpha_dDelta - delta*tau*d2alpha_dDelta_dTau)**2/(2*delta*dalpha_dDelta + delta**2*d2alpha_dDelta2))
#     properties['gamma'] = properties['cpmass']/properties['cvmass']
#     properties['compressibility_factor'] = delta*dalpha_dDelta
#     properties['speed_sound'] = ((R/M)*T*(2*delta*dalpha_dDelta + delta**2*d2alpha_dDelta2 - (delta*dalpha_dDelta - delta*tau*d2alpha_dDelta_dTau)**2/(tau**2*d2alpha_dTau2)))**0.5
#     properties['isentropic_bulk_modulus'] = rho*(R/M)*T*(2*delta*dalpha_dDelta + delta**2*d2alpha_dDelta2 - (delta*dalpha_dDelta - delta*tau*d2alpha_dDelta_dTau)**2/(tau**2*d2alpha_dTau2))
#     properties['isentropic_compressibility'] = 1 / properties["isentropic_bulk_modulus"]
#     properties['isothermal_bulk_modulus'] = (R/M)*T*rho*(2*delta*dalpha_dDelta + delta**2*d2alpha_dDelta2)
#     properties['isothermal_compressibility'] = 1/((R/M)*T*rho*(2*delta*dalpha_dDelta + delta**2*d2alpha_dDelta2))
#     properties['isobaric_expansion_coefficient'] = 1/T*(delta*dalpha_dDelta - delta*tau*d2alpha_dDelta_dTau)/(2*delta*dalpha_dDelta + delta**2*d2alpha_dDelta2)
#     properties['viscosity'] = fluid.viscosity()
#     properties['conductivity'] = fluid.conductivity()
#     properties['Q'] = np.nan
#     properties['quality_mass'] = np.nan
#     properties['quality_volume'] = np.nan
    
#     return properties




if __name__ == "__main__":
    fluid = Fluid('Water', backend="HEOS")

    # # Properties of liquid water
    # props_stable = fluid.set_state(CP.PT_INPUTS, 101325, 300)
    # print()
    # print("Properties of liquid water")
    # print(f"{'Property':35} {'value':6}")
    # for key, value in props_stable.items():
    #     print(f"{key:35} {value:.6e}")

    # # Properties of water/steam mixture
    # props = fluid.set_state(CP.QT_INPUTS, 0.5, 300)
    # print()
    # print("Properties of water/steam mixture")
    # print(f"{'Property':35} {'value':6}")
    # for key, value in props.items():
    #     print(f"{key:35} {value:.6e}")

    # # Get subset of properties for meanline code
    # props = fluid.compute_properties_meanline(CP.QT_INPUTS, 0.5, 300)
    # print()
    # print("Properties for the meanline code")
    # print(f"{'Property':15} {'value':6}")
    # for key, value in props.items():
    #     print(f"{key:15} {value:.6e}")

    # # 
    # props = compute_properties_metastable_rhoT(10, 500, fluid.abstractstate)
    # print("Metastable properties of water")
    # print(f"{'Property':35} {'value':6}")
    # for key, value in props.items():
    #     print(f"{key:35} {value:.6e}")

    # Check that the metastable property calculations match in the single-phase region
    p, T =  101325, 300
    props_stable = fluid.set_state(CP.PT_INPUTS, p, T)
    print()
    print(f"Properties of water at p={p} Pa and T={T} K")
    print(f"{'Property':35} {'Equilibrium':>15} {'Metastable':>15} {'Deviation':>15}")
    props_metastable = fluid.set_state_metastable_rhoT(props_stable["rho"], props_stable["T"])
    for key in props_stable.keys():
        value_stable = props_stable[key]
        value_metastable = props_metastable[key]
        print(f"{key:35} {value_stable:+15.6e} {value_metastable:+15.6e} {(value_stable - value_metastable)/value_stable:+15.6e}")

    


    
