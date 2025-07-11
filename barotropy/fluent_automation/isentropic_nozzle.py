import numpy as np
import scipy.optimize
import CoolProp.CoolProp as cp

from .. import properties as props


def compute_supersonic_exit(T_in, p_in, area_ratio, fluid):
    """Compute the static pressure at the exit of a supersonic convergent-divergent nozzle"""

    # Compute inlet state
    fluid.set_state(cp.PT_INPUTS, p_in, T_in)
    state_inlet = fluid.properties
    
    # Compute sonic state for constant entropy
    state_sonic = get_isentropic_sonic_state(state_inlet, fluid.fluid_name)

    # Compute outlet state for constant entropy
    state_exit = get_isentropic_outlet_state(area_ratio, state_inlet, state_sonic, fluid.fluid_name)

    return state_exit, state_sonic, state_inlet


def get_isentropic_sonic_state(state_in, fluid_name):

    # Define residual function
    def get_energy_eq_residual(p, s_in):
        _fluid.get_state(cp.PSmass_INPUTS, p, s_in)
        h = _fluid._properties["h"]
        a = _fluid._properties["speed_sound"]
        res = state_in["h"] -  h - a**2 / 2
        return res

    # Define a new fluid object to avoid mixing states
    _fluid = props.Fluid(fluid_name)

    # Define the initial guess
    _fluid.get_state(cp.PSmass_INPUTS, state_in["p"]*0.528, state_in["s"])
    s_in = state_in["s"]
    p_in = state_in["p"]

    # Solve the nonlinear equation
    sol = scipy.optimize.root_scalar(get_energy_eq_residual, args=(s_in, ), x0=0.528*p_in, bracket=[0.25*p_in, 0.75*p_in], method='brentq')

    # Check if the solver has converged
    if not sol.converged:
        raise ValueError("The root-finding algorithm did not converge!")

    # Compute the outlet state
    _fluid.get_state(cp.PSmass_INPUTS, sol.root, state_in["s"])

    return _fluid._properties


def get_isentropic_outlet_state(area_ratio, state_in, state_sonic, fluid_name):

    # Define residual function
    def get_energy_eq_residual(p, s_in):
        _fluid.get_state(cp.PSmass_INPUTS, p, s_in)
        h_out = _fluid._properties["h"]
        v_out = np.sqrt(2*(state_in["h"] - h_out))
        rho_out = _fluid._properties["rho"]
        res = (rho_out/state_sonic["rho"]) * (v_out/state_sonic["speed_sound"]) * area_ratio - 1
        return res

    # Define a new fluid object to avoid mixing states
    _fluid = props.Fluid(fluid_name)

    # Define the initial guess
    s_in = state_in["s"]
    p_sonic = state_sonic["p"]

    # Solve the nonlinear equation
    sol = scipy.optimize.root_scalar(get_energy_eq_residual, args=(s_in, ), bracket=[0.1*p_sonic, p_sonic], method='brentq')

    # Check if the solver has converged
    if not sol.converged:
        raise ValueError("The root-finding algorithm did not converge!")
    
    # Compute the outlet state
    _fluid.get_state(cp.PSmass_INPUTS, sol.root, state_in["s"])

    return _fluid._properties


if __name__ == "__main__":

    # Define fluid
    fluid_name = "water"
    fluid = props.Fluid(fluid_name)

    # Define inlet state
    T_in = 150+273.15
    p_in = 101325
    area_ratio = 2.00

    # Compute exit pressure
    state_exit, state_sonic, state_in = compute_supersonic_exit(T_in, p_in, area_ratio, fluid)

    # Display results
    print(f"Sonic pressure ratio: {state_sonic['p']/state_in['p']}")
    print(f"Exit pressure ratio: {state_exit['p']/state_in['p']}")
