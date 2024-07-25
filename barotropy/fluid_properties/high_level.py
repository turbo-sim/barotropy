import copy
import scipy
import numpy as np

import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP

from . import low_level as props
from . import utilities as utils

from ..pysolver_view import (
    NonlinearSystemSolver,
    NonlinearSystemProblem,
    OptimizationProblem,
    OptimizationSolver,
)


MEANLINE_PROPERTIES = [
    "p",
    "T",
    "h",
    "s",
    "d",
    "Z",
    "a",
    "mu",
    "k",
    "cp",
    "cv",
    "gamma",
]

# Dynamically add INPUTS fields to the module
# for attr in dir(CP):
#     if attr.endswith('_INPUTS'):
#         globals()[attr] = getattr(CP, attr)

# Statically add phase indices to the module (IDE autocomplete)
iphase_critical_point = CP.iphase_critical_point
iphase_gas = CP.iphase_gas
iphase_liquid = CP.iphase_liquid
iphase_not_imposed = CP.iphase_not_imposed
iphase_supercritical = CP.iphase_supercritical
iphase_supercritical_gas = CP.iphase_supercritical_gas
iphase_supercritical_liquid = CP.iphase_supercritical_liquid
iphase_twophase = CP.iphase_twophase
iphase_unknown = CP.iphase_unknown

# Statically add INPUT fields to the module (IDE autocomplete)
QT_INPUTS = CP.QT_INPUTS
PQ_INPUTS = CP.PQ_INPUTS
QSmolar_INPUTS = CP.QSmolar_INPUTS
QSmass_INPUTS = CP.QSmass_INPUTS
HmolarQ_INPUTS = CP.HmolarQ_INPUTS
HmassQ_INPUTS = CP.HmassQ_INPUTS
DmolarQ_INPUTS = CP.DmolarQ_INPUTS
DmassQ_INPUTS = CP.DmassQ_INPUTS
PT_INPUTS = CP.PT_INPUTS
DmassT_INPUTS = CP.DmassT_INPUTS
DmolarT_INPUTS = CP.DmolarT_INPUTS
HmolarT_INPUTS = CP.HmolarT_INPUTS
HmassT_INPUTS = CP.HmassT_INPUTS
SmolarT_INPUTS = CP.SmolarT_INPUTS
SmassT_INPUTS = CP.SmassT_INPUTS
TUmolar_INPUTS = CP.TUmolar_INPUTS
TUmass_INPUTS = CP.TUmass_INPUTS
DmassP_INPUTS = CP.DmassP_INPUTS
DmolarP_INPUTS = CP.DmolarP_INPUTS
HmassP_INPUTS = CP.HmassP_INPUTS
HmolarP_INPUTS = CP.HmolarP_INPUTS
PSmass_INPUTS = CP.PSmass_INPUTS
PSmolar_INPUTS = CP.PSmolar_INPUTS
PUmass_INPUTS = CP.PUmass_INPUTS
PUmolar_INPUTS = CP.PUmolar_INPUTS
HmassSmass_INPUTS = CP.HmassSmass_INPUTS
HmolarSmolar_INPUTS = CP.HmolarSmolar_INPUTS
SmassUmass_INPUTS = CP.SmassUmass_INPUTS
SmolarUmolar_INPUTS = CP.SmolarUmolar_INPUTS
DmassHmass_INPUTS = CP.DmassHmass_INPUTS
DmolarHmolar_INPUTS = CP.DmolarHmolar_INPUTS
DmassSmass_INPUTS = CP.DmassSmass_INPUTS
DmolarSmolar_INPUTS = CP.DmolarSmolar_INPUTS
DmassUmass_INPUTS = CP.DmassUmass_INPUTS
DmolarUmolar_INPUTS = CP.DmolarUmolar_INPUTS

# Define dictionary with dynamically generated fields
PHASE_INDEX = {attr: getattr(CP, attr) for attr in dir(CP) if attr.startswith("iphase")}
INPUT_PAIRS = {attr: getattr(CP, attr) for attr in dir(CP) if attr.endswith("_INPUTS")}
INPUT_PAIRS = sorted(INPUT_PAIRS.items(), key=lambda x: x[1])


class Fluid:
    """
    Represents a fluid with various thermodynamic properties computed via CoolProp.

    This class provides a convenient interface to CoolProp for various fluid property calculations.

    Critical and triple point properties are computed upon initialization and stored internally for convenience.

    Attributes
    ----------
    name : str
        Name of the fluid.
    backend : str
        Backend used for CoolProp, default is 'HEOS'.
    exceptions : bool
        Determines if exceptions should be raised during state calculations. Default is True.
    converged_flag : bool
        Flag indicating whether properties calculations converged.
    properties : dict
        Dictionary of various fluid properties. Accessible directly as attributes (e.g., `fluid.p` for pressure).
    critical_point : FluidState
        Properties at the fluid's critical point.
    triple_point_liquid : FluidState
        Properties at the fluid's triple point in the liquid state.
    triple_point_vapor : FluidState
        Properties at the fluid's triple point in the vapor state.

    Methods
    -------
    get_state(input_type, prop_1, prop_2):
        Set the thermodynamic state of the fluid using specified property inputs.

    Examples
    --------

    Calculating properties with Fluid.get_state()

    >>> fluid = bpy.Fluid(name="Water", backend="HEOS")
    >>> state = fluid.get_state(bpy.PT_INPUTS, 101325, 300)
    >>> print(f"Water density is {state.rho:0.2f} kg/m3 at p={state.p:0.2f} Pa and T={state.T:0.2f} K")
    Water density is 996.56 kg/m3 at p=101325.00 Pa and T=300.00 K

    >>> fluid = bpy.Fluid(name="Air", backend="HEOS")
    >>> state = fluid.get_state(bpy.PT_INPUTS, 101325, 300)
    >>> print(f"Air heat capacity ratio is {state.gamma:0.2f} at p={state.p:0.2f} Pa and T={state.T:0.2f} K")
    Air heat capacity ratio is 1.40 at p=101325.00 Pa and T=300.00 K


    Accessing critical point properties:

    >>> fluid.critical_point.p  # Retrieves critical pressure
    >>> fluid.critical_point['T']  # Retrieves critical temperature

    Accessing triple point properties:

    >>> fluid.triple_point_liquid.h  # Retrieves liquid enthalpy at the triple point
    >>> fluid.triple_point_vapor.s  # Retrieves vapor entropy at the triple point

    """

    def __init__(
        self,
        name,
        backend="HEOS",
        exceptions=True,
        generalize_quality=True,
        compute_superheating=True,
        compute_subcooling=True,
        # initialize_critical=True,
        # initialize_triple=True,
    ):
        self.name = name
        self.backend = backend
        self._AS = CP.AbstractState(backend, name)
        self.exceptions = exceptions
        self.converged_flag = False
        self._properties = {}

        # Define calculations peformend
        self.generalize_quality = generalize_quality
        self.compute_subcooling = compute_subcooling
        self.compute_superheating = compute_superheating

        # Initialize variables
        self.sat_liq = None
        self.sat_vap = None
        self.spdl_liq = None
        self.spdl_vap = None
        self.pseudo_critical_line = None
        # self.quality_grid = None
        self.q_mesh = None
        self.graphic_elements = {}

        # Get critical and triple point properties
        if self._AS.fluid_param_string("pure") == "true":
            self.critical_point = self._compute_critical_point()
            self.triple_point_liquid = self._compute_triple_point_liquid()
            self.triple_point_vapor = self._compute_triple_point_vapor()

        # # Assign critical point properties
        # if initialize_critical:
        #     self.critical_point = self._compute_critical_point()

        # # Assign triple point properties
        # if initialize_triple:
        #     self.triple_point_liquid = self._compute_triple_point_liquid()
        #     self.triple_point_vapor = self._compute_triple_point_vapor()

        # Pressure and temperature limits
        self.p_min = 1
        self.p_max = self._AS.pmax()
        self.T_min = self._AS.Tmin()
        self.T_max = self._AS.Tmax()

    def _compute_critical_point(self):
        """Calculate the properties at the critical point"""
        rho_crit, T_crit = self._AS.rhomass_critical(), self._AS.T_critical()
        return self.get_state(DmassT_INPUTS, rho_crit, T_crit)

    def _compute_triple_point_liquid(self):
        """Calculate the properties at the triple point (liquid state)"""
        return self.get_state(QT_INPUTS, 0.00, self._AS.Ttriple())

    def _compute_triple_point_vapor(self):
        """Calculate the properties at the triple point (vapor state)"""
        return self.get_state(QT_INPUTS, 1.00, self._AS.Ttriple())

    def get_state(
        self,
        input_type,
        prop_1,
        prop_2,
    ):
        """
        Set the thermodynamic state of the fluid based on input properties.

        This method updates the thermodynamic state of the fluid in the CoolProp ``abstractstate`` object
        using the given input properties. It then calculates either single-phase or two-phase
        properties based on the current phase of the fluid.

        If the calculation of properties fails, `converged_flag` is set to False, indicating an issue with
        the property calculation. Otherwise, it's set to True.

        Aliases of the properties are also added to the ``Fluid.properties`` dictionary for convenience.

        Parameters
        ----------
        input_type : int
            The variable pair used to define the thermodynamic state. This should be one of the
            predefined input pairs in CoolProp, such as ``PT_INPUTS`` for pressure and temperature.
            For all available input pairs, refer to :ref:`this list <module-input-pairs-table>`.
        prop_1 : float
            The first property value corresponding to the input type.
        prop_2 : float
            The second property value corresponding to the input type.

        Returns
        -------
        dict
            A dictionary of computed properties for the current state of the fluid. This includes both the
            raw properties from CoolProp and any additional alias properties.

        Raises
        ------
        Exception
            If `throw_exceptions` attribute is set to True and an error occurs during property calculation,
            the original exception is re-raised.


        """
        try:
            # Update Coolprop thermodynamic state
            self._AS.update(input_type, prop_1, prop_2)

            # Retrieve single-phase properties
            if self._AS.phase() != CP.iphase_twophase:
                self._properties = props.compute_properties_1phase(
                    self._AS,
                    self.generalize_quality,
                    self.compute_subcooling,
                    self.compute_superheating,
                )
            else:
                self._properties = props.compute_properties_2phase(
                    self._AS,
                    self.compute_subcooling,
                    self.compute_superheating,
                )

            # No errors computing the properies
            self.converged_flag = True

        # Something went wrong while computing the properties
        except Exception as e:
            self.converged_flag = False
            if self.exceptions:
                raise e

        # Return inmutable object
        return FluidState(self._properties, self.name)

    def get_state_metastable(
        self,
        prop_1,
        prop_1_value,
        prop_2,
        prop_2_value,
        rho_guess=None,
        T_guess=None,
        solver_algorithm="hybr",
        print_convergence=False,
    ):
        r"""
        Calculate the thermodynamic state based the Helmholtz HEOS (suitable to compute metastable states)
        
        Parameters
        ----------
        prop_1 : str
            The type of the first thermodynamic property (e.g., 'rho', 'T').
        prop_1_value : float
            The numerical value of the first thermodynamic property.
        prop_2 : str
            The type of the second thermodynamic property.
        prop_2_value : float
            The numerical value of the second thermodynamic property.
        rho_guess : float, optional
            Initial guess for density.
        T_guess : float, optional
            Initial guess for temperature.
        solver_algorithm : str, optional
            Method to be used for solving the nonlinear system. Defaults to 'hybr'.
        print_convergence : bool, optional
            If True, displays the convergence progress. Defaults to False.
        
        Returns
        -------
        dict
            Thermodynamic properties corresponding to the given input pair.
        
        Raises
        ------
        ValueError
            If the input property types are not density and temperature, and a valid initial guess is not provided.
        """
                

        # TODO: Add check to see if we are inside thespinodal and return two phase properties if yes

        # Check if state is outside spinodla line
        outside_spinodal = True

        # Ensure prop_1_value and prop_2_value are scalar numbers
        if not utils.is_float(prop_1_value) or not utils.is_float(prop_2_value):
            raise ValueError(f"Both prop_1_value and prop_2_value must be scalar numbers. Received: prop_1_value={prop_1_value}, prop_2_value={prop_2_value}")
        
        try:

            # Directly call the HEOS if the arguments are density-temperature
            if prop_1 == "rho" and prop_2 == "T":
                self._properties = props.compute_properties_metastable_rhoT(
                    self._AS, prop_1_value, prop_2_value
                )

            # Directly call the HEOS if the arguments are temperature-density
            elif prop_1 == "T" and prop_2 == "rho":
                self._properties = props.compute_properties_metastable_rhoT(
                    self._AS, prop_2_value, prop_1_value
                )

            # Use the solver for other input pairs
            else:

                # Validate initial guesses for rho and T if prop_1 and prop_2 are not 'rho' and 'T'
                if rho_guess is None or T_guess is None or not utils.is_float(rho_guess) or not utils.is_float(T_guess):
                    raise ValueError(f"If the input property types are not density and temperature, a valid initial guess must be provided for density and temperature. Received: rho_guess={rho_guess}, T_guess={T_guess}.")
                
                # Solve system of equations to determine the state
                self._properties = props.compute_properties_metastable(
                    self._AS,
                    prop_1,
                    prop_1_value,
                    prop_2,
                    prop_2_value,
                    rho_guess,
                    T_guess,
                    solver_algorithm=solver_algorithm,
                    print_convergence=print_convergence,
                )

        except Exception as e:
            raise RuntimeError(f"Failed to compute metastable properties: {str(e)}")





        #         # if T_guess is None:

        #     # else:
        #     #     self._properties = self.compute_properties_2phase()

        #     # No errors computing the properies
        #     self.converged_flag = True

        # # Something went wrong while computing the properties
        # except Exception as e:
        #     self.converged_flag = False
        #     if self.exceptions:
        #         raise e

        # Return inmutable object
        return FluidState(self._properties, self.name)

    def get_property(self, propname):
        """Get the value of a single property"""
        if propname in self._properties:
            return self._properties[propname]
        else:
            valid_options = "\n\t".join(self._properties.keys())
            raise ValueError(
                f"The requested property '{propname}' is not available. The valid options are:\n\t{valid_options}"
            )

    def compute_properties_meanline(self, input_type, prop_1, prop_2):
        """Extract fluid properties for meanline model"""

        # Compute properties in the normal way
        self.get_state(input_type, prop_1, prop_2)

        # Store a subset of the properties in a dictionary
        fluid_properties = {}
        for item in MEANLINE_PROPERTIES:
            fluid_properties[item] = self._properties[item]

        return fluid_properties

    def plot_phase_diagram(
        self,
        x_variable="s",
        y_variable="T",
        axes=None,
        N=100,
        plot_saturation_line=True,
        plot_critical_point=True,
        plot_triple_point_liquid=False,
        plot_triple_point_vapor=False,
        plot_spinodal_line=False,
        spinodal_line_color=0.5 * np.array([1, 1, 1]),
        spinodal_line_width=1.25,
        plot_quality_isolines=False,
        plot_pseudocritical_line=False,
        quality_levels=np.linspace(0.1, 1.0, 10),
        quality_labels=False,
        show_in_legend=False,
        # **kwargs,
    ):
        if axes is None:
            axes = plt.gca()

        # Saturation line
        if plot_saturation_line:
            if self.sat_liq is None or self.sat_vap is None:
                self.sat_liq, self.sat_vap = utils.compute_saturation_line(self, N)
            x = self.sat_liq[x_variable] + self.sat_vap[x_variable]
            y = self.sat_liq[y_variable] + self.sat_vap[y_variable]
            label = self._get_label("Saturation line", show_in_legend)
            params = {"label": label, "color": "black"}
            self._graphic_saturation_line = self._plot_or_update_line(
                axes,
                x,
                y,
                "saturation_line",
                **params,
            )
        else:
            self._set_visibility(axes, "saturation_line", False)

        # Spinodal line
        if plot_spinodal_line:
            if self.spdl_liq is None or self.spdl_vap is None:
                self.spdl_liq, self.spdl_vap = utils.compute_spinodal_line(self, N)
            x = self.spdl_liq[x_variable] + self.spdl_vap[x_variable]
            y = self.spdl_liq[y_variable] + self.spdl_vap[y_variable]
            label = self._get_label("Spinodal line", show_in_legend)
            params = {
                "label": label,
                "color": spinodal_line_color,
                "linewidth": spinodal_line_width,
            }
            self._graphic_spinodal_line = self._plot_or_update_line(
                axes,
                x,
                y,
                "spinodal_line",
                **params,
            )
        else:
            self._set_visibility(axes, "spinodal_line", False)

        # Plot pseudocritical line
        if plot_pseudocritical_line:
            if self.pseudo_critical_line is None:
                self.pseudo_critical_line = utils.compute_pseudocritical_line(self)
            x = self.pseudo_critical_line[x_variable]
            y = self.pseudo_critical_line[y_variable]
            label = self._get_label("Pseudocritical line", show_in_legend)
            params = {
                "label": label,
                "color": "black",
                "linestyle": "--",
                "linewidth": 0.75,
            }
            self._graphic_pseudocritical_line = self._plot_or_update_line(
                axes,
                x,
                y,
                "pseudocritical_line",
                **params,
            )
        else:
            self._set_visibility(axes, "pseudocritical_line", False)

        # Plot quality isolines
        if plot_quality_isolines:
            if self.q_mesh is None:
                self.q_mesh = utils.compute_quality_grid(self, N, quality_levels)
            x = self.q_mesh[x_variable]
            y = self.q_mesh[y_variable]
            _, m = np.shape(x)
            z = np.tile(quality_levels, (m, 1)).T
            params = {"colors": "black", "linestyles": ":", "linewidths": 0.75}
            self._graphics_q_lines = self._plot_or_update_contours(
                axes,
                x,
                y,
                z,
                quality_levels,
                "quality_isolines",
                **params,
            )

            if quality_labels:
                axes.clabel(self._graphics_q_lines, fontsize=9, rightside_up=True)

        else:
            # Remove existing contour lines if they exist
            if "quality_isolines" in self.graphic_elements.get(axes, {}):
                for coll in self.graphic_elements[axes]["quality_isolines"].collections:
                    coll.remove()
                del self.graphic_elements[axes]["quality_isolines"]

        # Plot critical point
        params = {
            "color": "black",
            "marker": "o",
            "markersize": 4.5,
            "markerfacecolor": "w",
        }
        if plot_critical_point:
            x = self.critical_point[x_variable]
            y = self.critical_point[y_variable]
            label = self._get_label("Critical point", show_in_legend)
            self._graphic_critical_point = self._plot_or_update_line(
                axes,
                x,
                y,
                "critical_point",
                label=label,
                **params,
            )
        else:
            self._set_visibility(axes, "critical_point", False)

        # Plot liquid triple point
        if plot_triple_point_liquid:
            x = self.triple_point_liquid[x_variable]
            y = self.triple_point_liquid[y_variable]
            label = self._get_label("Triple point liquid", show_in_legend)
            self._graphic_triple_point_liquid = self._plot_or_update_line(
                axes,
                x,
                y,
                "triple_point_liquid",
                label=label,
                **params,
            )
        else:
            self._set_visibility(axes, "triple_point_liquid", False)

        # Plot vapor triple point
        if plot_triple_point_vapor:
            x = self.triple_point_vapor[x_variable]
            y = self.triple_point_vapor[y_variable]
            label = self._get_label("Triple point vapor", show_in_legend)
            self._graphic_triple_point_vapor = self._plot_or_update_line(
                axes,
                x,
                y,
                "triple_point_vapor",
                label=label,
                **params,
            )
        else:
            self._set_visibility(axes, "triple_point_vapor", False)

        return axes

    def _get_label(self, label, show_in_legend):
        """Returns the appropriate label value based on whether it should be shown in the legend."""
        return label if show_in_legend else "_no_legend_"

    def _plot_or_update_line(self, axes, x_data, y_data, line_name, **plot_params):
        # Ensure there is a dictionary for this axes
        if axes not in self.graphic_elements:
            self.graphic_elements[axes] = {}

        # Check if the line exists for this axes
        if line_name in self.graphic_elements[axes]:
            line = self.graphic_elements[axes][line_name]
            line.set_data(x_data, y_data)
            # Update line properties
            for param, value in plot_params.items():
                setattr(line, param, value)
            line.set_visible(True)
        else:
            # Create a new line with the provided plot parameters
            (line,) = axes.plot(x_data, y_data, **plot_params)
            self.graphic_elements[axes][line_name] = line
        return line

    def _plot_or_update_contours(
        self, axes, x_data, y_data, z_data, contour_levels, line_name, **contour_params
    ):
        # Ensure there is a dictionary for this axes
        if axes not in self.graphic_elements:
            self.graphic_elements[axes] = {}

        # Check if the contour exists for this axes
        if line_name in self.graphic_elements[axes]:
            for coll in self.graphic_elements[axes][line_name].collections:
                coll.remove()  # Remove the old contour collections

        # Create a new contour
        contour = axes.contour(x_data, y_data, z_data, contour_levels, **contour_params)
        self.graphic_elements[axes][line_name] = contour
        return contour

    def _set_visibility(self, axes, line_name, visible):
        if axes in self.graphic_elements and line_name in self.graphic_elements[axes]:
            self.graphic_elements[axes][line_name].set_visible(visible)


# ------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------ #


class FluidState:
    """
    An immutable class representing the thermodynamic state of a fluid.

    This class is designed to provide a read-only representation of a fluid's state,
    with properties accessible through both attribute-style and dictionary-style access.
    The state of a fluid is defined at initialization and cannot be altered thereafter,
    ensuring the integrity of the data.

    Instances of this class store fluid properties and the fluid's name, providing
    methods to access these properties but not to modify them.

    Parameters
    ----------
    properties : dict
        A dictionary where keys are property names (as strings) and values are the
        corresponding properties of the fluid. This dictionary is converted to an
        immutable internal representation.
    fluid_name : str
        The name of the fluid.

    Attributes
    ----------
    properties : dict
        An internal dictionary storing the properties of the fluid state. This attribute
        is immutable and cannot be modified after the object's initialization.
    fluid_name : str
        The name of the fluid. Immutable after initialization.

    Methods
    -------
    to_dict():
        Returns a copy of the fluid properties as a dictionary.
    keys():
        Returns the keys of the fluid properties.
    items():
        Returns the items (key-value pairs) of the fluid properties.
    """

    __slots__ = ("_properties", "fluid_name")  # Define allowed attributes

    def __init__(self, properties, fluid_name):
        # Convert keys to strings and store properties in an internal attribute
        # Ensure the properties dictionary is immutable (e.g., by using a frozendict if modifications are a concern)
        object.__setattr__(
            self, "_properties", {str(k): v for k, v in properties.items()}
        )
        object.__setattr__(self, "fluid_name", fluid_name)

    def __getattr__(self, name):
        # Allows attribute-style access to the fluid properties. If the property does not exist, raises AttributeError.
        try:
            return self._properties[name]
        except KeyError:
            raise AttributeError(f"'{name}' not found in FluidState properties")

    def __getitem__(self, key):
        # Allows dictionary-style access to the fluid properties. If the key does not exist, raises KeyError.
        try:
            return self._properties[str(key)]
        except KeyError:
            raise KeyError(f"'{key}' not found in FluidState properties")

    def __setattr__(self, key, value):
        # Prevents modifications to the instance attributes, ensuring immutability.
        raise AttributeError(
            "Cannot modify properties of an immutable FluidState class"
        )

    def __setitem__(self, key, value):
        # Prevents modifications to the fluid properties through dictionary-style access, ensuring immutability.
        raise AttributeError(
            "Cannot modify properties of an immutable FluidState class"
        )

    def __repr__(self):
        # Returns a string representation of the FluidState instance, including its class name, properties, and fluid name.
        return f"{self.__class__.__name__}({self._properties}, '{self.fluid_name}')"

    def __str__(self):
        prop_str = "\n   ".join([f"{k}: {v:6e}" for k, v in self._properties.items()])
        return f"FluidState:\n   {prop_str}"

    def __iter__(self):
        # Iterate over the object like a dictionary
        return iter(self._properties)

    def to_dict(self):
        # Return a copy of the properties to ensure immutability
        return dict(self._properties)

    def keys(self):
        # Return the keys of the properties
        return self._properties.keys()

    def values(self):
        # Return the values of the properties
        return self._properties.values()

    def items(self):
        # Return the items of the properties
        return self._properties.items()

    def get(self, key, default=None):
        return self._properties.get(key, default)


# ------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------ #


def compute_saturation_line(fluid, N_points=100):
    # Initialize objects to store properties
    prop_names = fluid._properties.keys()
    saturation_liq = {name: [] for name in prop_names}
    saturation_vap = {name: [] for name in prop_names}

    # Define temperature array with refinement close to the critical point
    R = 1 - fluid.triple_point_liquid.T / fluid.critical_point.T
    t1 = np.logspace(np.log10(1 - 0.9999), np.log10(R / 10), int(np.ceil(N_points / 2)))
    t2 = np.logspace(np.log10(R / 10), np.log10(R), int(np.floor(N_points / 2)))
    T_sat = (1 - np.concatenate([t1, t2])) * fluid.critical_point.T

    # Loop over temperatures and property names
    for T in T_sat:
        state_liquid = fluid.get_state(CP.QT_INPUTS, 0.00, T)
        state_vapor = fluid.get_state(CP.QT_INPUTS, 1.00, T)
        for name in prop_names:
            saturation_liq[name].append(state_liquid[name])
            saturation_vap[name].append(state_vapor[name])

    # Add critical point as part of the spinodal line
    for name in prop_names:
        saturation_liq[name] = [fluid.critical_point[name]] + saturation_liq[name]
        saturation_vap[name] = [fluid.critical_point[name]] + saturation_vap[name]

    # Re-format for easy concatenation
    for name in prop_names:
        saturation_liq[name] = list(reversed(saturation_liq[name]))

    return saturation_liq, saturation_vap


def compute_spinodal_line(fluid, N=25):
    # Initialize objects to store properties
    prop_names = fluid._properties.keys()
    spinodal_liq = {name: [] for name in prop_names}
    spinodal_vap = {name: [] for name in prop_names}

    # Temperature array with refinement close to the critical point
    ratio = 1 - 1.1 * fluid.triple_point_liquid.T / fluid.critical_point.T
    t1 = np.logspace(np.log10(1 - 0.999), np.log10(ratio / 10), int(np.ceil(N / 2)))
    t2 = np.logspace(np.log10(ratio / 10), np.log10(ratio), int(np.floor(N / 2)))
    T_spinodal = (1 - np.concatenate([t1, t2])) * fluid.critical_point.T

    # Initial properties for the first temperature point
    props_liq = props.compute_spinodal_point(T_spinodal[0], fluid._AS, "liquid")
    props_vap = props.compute_spinodal_point(T_spinodal[0], fluid._AS, "vapor")

    # Loop over temperatures and property names
    for T in T_spinodal:
        props_liq = props.compute_spinodal_point(T, fluid._AS, "liquid")
        props_vap = props.compute_spinodal_point(T, fluid._AS, "vapor")
        # props_liq = functions.compute_spinodal_point(T, fluid._AS, 'liquid', rho_guess=props_liq["rho"])
        # props_vap = functions.compute_spinodal_point(T, fluid._AS, 'vapor', rho_guess=props_vap["rho"])
        for name in prop_names:
            spinodal_liq[name].append(props_liq[name])
            spinodal_vap[name].append(props_vap[name])

    # Add critical point as part of the spinodal line
    for name in prop_names:
        spinodal_liq[name] = [fluid.critical_point[name]] + spinodal_liq[name]
        spinodal_vap[name] = [fluid.critical_point[name]] + spinodal_vap[name]

    # Re-format for easy concatenation
    for name in prop_names:
        spinodal_liq[name] = list(reversed(spinodal_liq[name]))

    return spinodal_liq, spinodal_vap


def compute_pseudocritical_line(fluid, N_points=100):
    # Initialize objects to store properties
    prop_names = fluid._properties.keys()
    pseudocritical_line = {name: [] for name in prop_names}

    # Define temperature array with refinement close to the critical point
    tau = np.logspace(np.log10(1e-3), np.log10(1), N_points)
    T_range = (1 + tau) * fluid.critical_point.T

    # Loop over temperatures and compute pseudocritical properties
    for T in T_range:
        for name in prop_names:
            state = fluid.get_state(DmassT_INPUTS, fluid.critical_point.d, T)
            pseudocritical_line[name].append(state[name])

    return pseudocritical_line


def compute_quality_grid(fluid, num_points, quality_levels):
    # Define temperature levels
    t1 = np.logspace(np.log10(1 - 0.9999), np.log10(0.1), int(num_points / 2))
    t2 = np.logspace(
        np.log10(0.1),
        np.log10(1 - (fluid.triple_point_liquid.T) / fluid.critical_point.T),
        int(num_points / 2),
    )
    temperature_levels = (1 - np.hstack((t1, t2))) * fluid.critical_point.T

    # Calculate property grid
    quality_grid = []
    for q in quality_levels:
        row = []
        for T in temperature_levels:
            row.append(fluid.get_state(CP.QT_INPUTS, q, T))
        quality_grid.append(row)

    return states_to_dict_2d(quality_grid)



def compute_property_grid(
    fluid,
    input_pair,
    range_1,
    range_2,
):
    """
    Compute fluid properties over a specified range and store them in a dictionary.

    This function creates a meshgrid of property values based on the specified ranges and input pair,
    computes the properties of the fluid at each point on the grid, and stores the results in a
    dictionary where each key corresponds to a fluid property.

    Parameters
    ----------
    fluid : Fluid object
        An instance of the Fluid class.
    input_pair : tuple
        The input pair specifying the property type (e.g., PT_INPUTS for pressure-temperature).
    range1 : tuple
        The range linspace(min, max, n) for the first property of the input pair.
    range2 : tuple
        The range linspace(min, max, n) for the second property of the input pair.

    Returns
    -------
    properties_dict : dict
        A dictionary where keys are property names and values are 2D numpy arrays of computed properties.
    grid1, grid2 : numpy.ndarray
        The meshgrid arrays for the first and second properties.
    """

    # Create the meshgrid
    grid1, grid2 = np.meshgrid(range_1, range_2)

    # Initialize dictionary to store properties and pre-allocate storage
    properties_dict = {key: np.zeros_like(grid1) for key in fluid._properties}

    # Compute properties at each point
    for i in range(len(range_2)):
        for j in range(len(range_1)):
            # Set state of the fluid
            state = fluid.get_state(
                input_pair,
                grid1[i, j],
                grid2[i, j],
            )

            # Store the properties
            for key in state:
                properties_dict[key][i, j] = state[key]

    return properties_dict


def compute_property_grid_rhoT(
    fluid,
    rho_array,
    T_array,
):

    # Calculate property grid
    states_meta = []
    for T in T_array:
        row = []
        for rho in rho_array:
            row.append(props.compute_properties_metastable_rhoT(rho, T, fluid._AS))
        states_meta.append(row)

    # Convert nested list of dictionaries into dictionary of 2D arrays
    states_meta = utils.states_to_dict_2d(states_meta)

    return states_meta


# TODO  Implement class to calculate intersection with saturation line?





def get_state_Qs(fluid, Q, s):
    # Define the residual equation
    def get_residual(p):
        state = fluid.get_state(PSmass_INPUTS, p, s, generalize_quality=True)
        residual = Q - state.Q
        return residual

    # Solve the scalar equation
    p_triple = 1.0 * fluid.triple_point_liquid.p
    p_critical = 1.25 * fluid.critical_point.p
    bounds = [p_triple, p_critical]
    sol = scipy.optimize.root_scalar(get_residual, bracket=bounds, method="brentq")

    # Check if the solver has converged
    if not sol.converged:
        raise ValueError("The root-finding algorithm did not converge!")

    # Compute the outlet state
    state = fluid.get_state(PSmass_INPUTS, sol.root, s, generalize_quality=True)

    return state


def get_isentropic_saturation_state(fluid, s_in):
    # Calculate saturation sate
    if s_in < fluid.critical_point.s:
        state_sat = get_state_Qs(fluid, Q=0.00, s=s_in)
    else:
        state_sat = get_state_Qs(fluid, Q=1.00, s=s_in)

    return state_sat





# # ------------------------------------------------------------------------------------ #
# # Sonic point calculations
# # ------------------------------------------------------------------------------------ #


# class _SonicStateProblem(psv.NonlinearSystemProblem):
#     """ """

#     def __init__(self, fluid, property_pair, prop_1, prop_2):
#         # Calculate the thermodynamic state
#         self.fluid = fluid
#         self.state = fluid.get_state(property_pair, prop_1, prop_2)

#         # Initial guess based in input sstate
#         self.initial_guess = [self.state.d, self.state.T]

#         # # Initial guess based on perfect gass relations
#         # gamma = self.state.gamma
#         # d_star = (2/(gamma + 1)) ** (1/(gamma-1)) * self.state.rho
#         # T_star =  (2/(gamma + 1)) * self.state.T
#         # self.initial_guess = [d_star, T_star]

#     def get_values(self, x):
#         # Ensure x can be indexed and contains exactly two elements
#         if not hasattr(x, "__getitem__") or len(x) != 2:
#             raise ValueError(
#                 "Input x must be a list, tuple or numpy array containing exactly two elements: density and temperature."
#             )

#         # Calculate state for the current density-temperature pair
#         crit_state = self.fluid.get_state(DmassT_INPUTS, x[0], x[1])

#         # Calculate the sonic state residual
#         residual = np.asarray(
#             [
#                 1.0 - (crit_state.h + 0.5 * crit_state.a**2) / self.state.h,
#                 1.0 - crit_state.s / self.state.s,
#             ]
#         )

#         return residual


# class _SonicStateProblem2(psv.OptimizationProblem):
#     """ """

#     def __init__(self, fluid, property_pair, prop_1, prop_2):
#         # Calculate the thermodynamic state
#         self.fluid = fluid
#         self.state = fluid.get_state(property_pair, prop_1, prop_2)

#         # Initial guess based in input sstate
#         self.initial_guess = [self.state.d, self.state.T * 0.9]

#         # # Initial guess based on perfect gass relations
#         # gamma = self.state.gamma
#         # d_star = (2/(gamma + 1)) ** (1/(gamma-1)) * self.state.rho
#         # T_star =  (2/(gamma + 1)) * self.state.T
#         # self.initial_guess = [d_star, T_star]

#     def get_values(self, x):
#         """
#         Compute the residuals for the given density and temperature.

#         Parameters
#         ----------
#         x : list
#             List containing the values for density and temperature.

#         Returns
#         -------
#         np.ndarray
#             Array containing residuals (difference) for the two properties.
#         """

#         # Ensure x can be indexed and contains exactly two elements
#         if not hasattr(x, "__getitem__") or len(x) != 2:
#             raise ValueError(
#                 "Input x must be a list, tuple or numpy array containing exactly two elements: density and temperature."
#             )

#         # Calculate state for the current density-temperature pair
#         crit_state = self.fluid.get_state(DmassT_INPUTS, x[0], x[1])

#         # Calculate the sonic state residual
#         residual = [
#             1.0 - crit_state.s / self.state.s,
#         ]

#         # Objective function
#         self.f = crit_state.d**2 * (self.state.h - crit_state.h)
#         self.f = -self.f / (self.state.d * self.state.a) ** 2

#         # Equality constraints
#         self.c_eq = residual

#         # No inequality constraints given for this problem
#         self.c_ineq = []

#         # Combine objective function and constraints
#         objective_and_constraints = self.merge_objective_and_constraints(
#             self.f, self.c_eq, self.c_ineq
#         )

#         return objective_and_constraints

#     def get_bounds(self):
#         bound_density = (
#             self.fluid.triple_point_vapor.d * 1.5,
#             self.fluid.critical_point.d * 3,
#         )
#         bound_temperature = (
#             self.fluid.triple_point_vapor.T * 1,
#             self.fluid.critical_point.T * 3,
#         )
#         # return [bound_density, bound_temperature]
#         return None

#     def get_n_eq(self):
#         return self.get_number_of_constraints(self.c_eq)

#     def get_n_ineq(self):
#         return self.get_number_of_constraints(self.c_ineq)
    



