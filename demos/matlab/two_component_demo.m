%% Barotropic fluid model calculation minimal working example
% Clear the workspace
clear all
close all
clc

% Define inlet and outlet conditions
T_in = 22.00 + 273.15;
p_in = 2000.0e3;
p_out = 98.6e3;

% Define the mixture composition    
R = 67.00;
y_1 = R./(1 + R);
y_2 = 1./(1 + R);

% Define polytropic efficiency of the process
eta_poly = 1.00;  % Isentropic process

% Initialize the fluid objects
py.importlib.import_module('CoolProp.CoolProp');
fluid_1 = py.CoolProp.CoolProp.AbstractState('HEOS', 'Water');
fluid_2 = py.CoolProp.CoolProp.AbstractState('HEOS', 'Nitrogen');

% Calculate the initial condition for enthalpy
fluid_1.update(py.CoolProp.CoolProp.PT_INPUTS, p_in, T_in);
fluid_2.update(py.CoolProp.CoolProp.PT_INPUTS, p_in, T_in);
h_in = y_1*fluid_1.hmass + y_2*fluid_2.hmass;

% Solve the ODE system
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);  % Tight integration tolerance
ode_handle = @(p, grad_hT) polytropic_expansion_ode(p, grad_hT, y_1, y_2, fluid_1, fluid_2, eta_poly);
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


function [grad_hT, props_out] = polytropic_expansion_ode(p, hT, y_1, y_2, fluid_1, fluid_2, eta_poly)
    
    % Rename variables
    h = hT(1);
    T = hT(2);

    % Update thermodynamic state
    fluid_1.update(py.CoolProp.CoolProp.PT_INPUTS, p, T);
    fluid_2.update(py.CoolProp.CoolProp.PT_INPUTS, p, T);
    
    % Compute additional thermodynamic properties
    cp = y_1*fluid_1.cpmass + y_2*fluid_2.cpmass;
    rho = 1/(y_1/fluid_1.rhomass + y_2/fluid_2.rhomass);
    dhdp_T_1 = (1 - T*fluid_1.isobaric_expansion_coefficient)/fluid_1.rhomass;
    dhdp_T_2 = (1 - T*fluid_2.isobaric_expansion_coefficient)/fluid_2.rhomass;
    dhdp_T = y_1*dhdp_T_1 + y_2*dhdp_T_2;

    % Compute the slope of the polytropic process
    dhdp = eta_poly/rho;
    dTdp = (dhdp - dhdp_T)/cp;
    grad_hT = [dhdp; dTdp];  % Right hand side of ODE system

    % Prepare additional properties for export
    props_out.p = p;
    props_out.T = T;
    props_out.mass_frac_1 = y_1;
    props_out.mass_frac_2 = y_2;
    props_out.vol_frac_1 = rho/fluid_1.rhomass*y_1;
    props_out.vol_frac_2 = rho/fluid_2.rhomass*y_2;
    props_out.mixture_ratio = y_1/y_2;
    props_out.hmass_1 = fluid_1.hmass;
    props_out.hmass_2 = fluid_2.hmass;
    props_out.hmass = h;
    props_out.hmass_bis = y_1*fluid_1.hmass + y_2*fluid_2.hmass;
    props_out.hmass_error = (h - y_1*fluid_1.hmass - y_2*fluid_2.hmass);
    props_out.smass_1 = fluid_1.smass;
    props_out.smass_2 = fluid_2.smass;
    props_out.smass = y_1*fluid_1.smass + y_2*fluid_2.smass;
    props_out.rhomass_1 = fluid_1.rhomass;
    props_out.rhomass_2 = fluid_2.rhomass;
    props_out.rhomass = rho;
    props_out.viscosity_1 = fluid_1.viscosity;
    props_out.viscosity_2 = fluid_2.viscosity;
    props_out.viscosity = props_out.vol_frac_1*fluid_1.viscosity + ...
                          props_out.vol_frac_2*fluid_2.viscosity;
    props_out.speed_sound_1 = fluid_1.speed_sound;
    props_out.speed_sound_2 = fluid_2.speed_sound;
    bulk_modulus_1 = fluid_1.rhomass*fluid_1.speed_sound^2;
    bulk_modulus_2 = fluid_2.rhomass*fluid_2.speed_sound^2;
    bulk_modulus = (props_out.vol_frac_1/bulk_modulus_1 + props_out.vol_frac_2/bulk_modulus_2)^(-1);
    props_out.speed_sound = sqrt(bulk_modulus/rho);


end