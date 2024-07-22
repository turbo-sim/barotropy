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


%% Define flow mixture
% Import CoolProp
py.importlib.import_module('CoolProp.CoolProp');

% Define the fluid mixture as 'Water&Nitrogen'
fluid_mixture = 'Water&Nitrogen';

% Define mass fractions of each component in the mixture
% In this example, the mixture contains 70% water and 30% nitrogen by mass
mass_fractions = [0.05, 0.95];

% Define the state using the low-level interface
fluid = py.CoolProp.CoolProp.AbstractState('REFPROP', fluid_mixture);

% Set the mass fractions
fluid.set_mass_fractions(mass_fractions);

% Define the input pressure and temperature
% Pressure in Pascals, Temperature in Kelvin
pressure = 101325;  % 1 atm
temperature = 300;  % Room temperature


%% Compare different input arguments
% Calculate properties based on pressure and temperature
fluid.update(py.CoolProp.CoolProp.PT_INPUTS, pressure, temperature);
entropy = fluid.smass;
density = fluid.rhomass;
temperature = fluid.T;

% Calculate properties based on density and pressure
fluid.update(py.CoolProp.CoolProp.PSmass_INPUTS, pressure, entropy)
entropy_bis = fluid.smass;
density_bis = fluid.rhomass;
temperature_bis = fluid.T;


% Display a header
fprintf('%15s %15s %15s %15s\n', 'Property', 'PT call', 'PS call', 'Deviation');
fprintf('---------------------------------------------------------------\n');

% Calculate and display the error (difference) for each property
fprintf('%15s %15.4f %15.4f %15.4e\n', 'Entropy', entropy, entropy_bis, entropy - entropy_bis);
fprintf('%15s %15.4f %15.4f %15.4e\n', 'Density', density, density_bis, density - density_bis);
fprintf('%15s %15.4f %15.4f %15.4e\n', 'Temperature', temperature, temperature_bis, temperature - temperature_bis);


%% Plot the phase envelope
% Build the phase envelope
fluid.build_phase_envelope('');

% Retrieve the phase envelope data
phase_envelope = fluid.get_phase_envelope_data();

% Plot the phase envelope
figure(); ax = gca; hold on; box on;
xlabel("Density (kg/m$^3$")
ylabel("Pressure (Pa)")
plot(phase_envelope.rhomolar_liq, phase_envelope.p, 'bo-', MarkerfaceColor='w', Displayname='Liquid line')
plot(phase_envelope.rhomolar_vap, phase_envelope.p, 'ro-', MarkerfaceColor='w', Displayname='Vapor line')
legend(Location="northeast")

