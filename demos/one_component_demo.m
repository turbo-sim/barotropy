%% Barotropic fluid model calculation minimal working example
% Clear the workspace
clear all
close all
clc

% Define plot settings
set_plot_options()

% Define inlet and outlet conditions
T_in = 22.00 + 273.15;
p_in = 2000.0e3;
p_out = 98.6e3;

% Define polytropic efficiency of the process
eta_poly = 1.0;  % Isentropic process

% Initialize the fluid objects
py.importlib.import_module('CoolProp.CoolProp');
fluid = py.CoolProp.CoolProp.AbstractState('HEOS', 'Water');

% Calculate the initial condition for enthalpy
fluid.update(py.CoolProp.CoolProp.PT_INPUTS, p_in, T_in);
h_in = fluid.hmass;

% Solve the ODE system
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);  % Tight integration tolerance
ode_handle = @(p, grad_hT) polytropic_expansion_ode(p, grad_hT, fluid, eta_poly);
[p, hT] = ode45(ode_handle, [p_in, p_out], [h_in, T_in], options);
props = postprocess_ode(p, hT, ode_handle);

% Plot evolution of temperature
figure(); hold on; box on;
xlabel('Pressure (bar)')
ylabel('Temperature ($^{\circ}$C)')
plot(p/1e5, props.T-273.15, 'k-')

% Plot evolution of entropy
figure(); hold on; box on;
xlabel('Pressure (bar)')
ylabel('Entropy (J/kg/K)')
plot(p/1e5, props.smass-props.smass(1), 'k-')


function [grad_hT, props_out] = polytropic_expansion_ode(p, hT, fluid, eta_poly)
    
    % Rename variables
    h = hT(1);
    T = hT(2);

    % Update thermodynamic state
    fluid.update(py.CoolProp.CoolProp.PT_INPUTS, p, T);
    
    % Compute additional thermodynamic properties
    cp = fluid.cpmass;
    rho = fluid.rhomass;
    dhdp_T = (1 - T*fluid.isobaric_expansion_coefficient)/fluid.rhomass;

    % Compute the slope of the polytropic process
    dhdp = eta_poly/rho;
    dTdp = (dhdp - dhdp_T)/cp;
    grad_hT = [dhdp; dTdp];  % Right hand side of ODE system

    % Prepare additional properties for export
    props_out.p = p;
    props_out.T = T;
    props_out.hmass = fluid.hmass;
    props_out.hmass_error = (h - fluid.hmass);
    props_out.smass = fluid.smass;
    props_out.rhomass = fluid.rhomass;
    props_out.speed_sound = fluid.speed_sound;

end