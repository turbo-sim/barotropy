%% Initialize script
% Clear the workspace
clear all
close all
clc

% Create folder to save results
results_path = 'results';
if not(isfolder(results_path))
    mkdir(results_path)
end

% Define plot settings
set_plot_options()
save_figures = false;
defaultColors = get(groot, 'factoryAxesColorOrder');


%% Plot the phase envelope
% Import CoolProp
py.importlib.import_module('CoolProp.CoolProp');

% Define the fluid mixture as 'Water&Nitrogen'
fluid_mixture = 'Water';

% Define the state using the low-level interface
fluid = py.CoolProp.CoolProp.AbstractState('REFPROP', fluid_mixture);

% Plot the T-s diagram
fig = figure(); ax = gca; hold on; box on; grid off
ax.XLabel.String = "Entropy (J/kg/K)";
ax.YLabel.String = "Temperature (K)";
ax.XScale = "linear";
ax.YScale = "linear";
plot_phase_diagram('smass', 'T', fluid)


%% Plot the isentrope with entropy inputs
% Define the polytropic efficiency of the process
eta_poly = 1.00;

% Define the inlet state
p_in = 10e5;
T_in = 500;
fluid.update(py.CoolProp.CoolProp.PT_INPUTS, p_in, T_in);
s_in = fluid.smass;
h_in = fluid.hmass;
rho_in = fluid.rhomass;

% Define the outlet state
p_out = 1e5;
s_out = s_in;
fluid.update(py.CoolProp.CoolProp.PSmass_INPUTS, p_out, s_out);
T_out = fluid.T;
rho_out = fluid.rhomass;

% Plot the isentrope
plot([s_in, s_out], [T_in, T_out], 'ro-', MarkerFaceColor='white', DisplayName="Entropy inputs")


%% Plot the isentrope as ODE integration
p = linspace(p_in, p_out, 100);

% Compute density/temperature along isentrope
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);
[p, h] = ode45(@(p, h) polytropic_expansion_ode(p, h, fluid, eta_poly), p, h_in, options);

% Compute the entropy to verify calculations
s = p*0;
T = p*0;
for i = 1:numel(p)
    fluid.update(py.CoolProp.CoolProp.HmassP_INPUTS, h(i), p(i))
    s(i) = fluid.smass;
    T(i) = fluid.T;
end

% Plot the isentrope
plot(s, T, 'ko', MarkerFaceColor="w", MarkerSize=2, DisplayName="ODE integration")

% Compute the polytropic process slope
function dh_dp = polytropic_expansion_ode(p, h, fluid, eta_poly)
    fluid.update(py.CoolProp.CoolProp.HmassP_INPUTS, h, p);
    dh_dp = eta_poly/fluid.rhomass;
end
