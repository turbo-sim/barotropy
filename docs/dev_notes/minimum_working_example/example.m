%% Minimum working example of the barotropic model
% Clear the workspace
clear all %#ok
close all
clc

% Add path to utility functions
addpath(genpath("dependencies"))

% Set plot options
save_figures = false;
set_plot_options()

% Define the fluid name and create the fluid object
fluid_name = 'CO2';
fluid = FluidCoolProp_2Phase('HEOS', fluid_name);

% Define the inlet state
T_in = 300; % Kelvin degrees
p_in = 90e5; % Pascals

% Extrapolate equation of state between the saturation and spinodal lines
include_metastable = true;

% Define the order of the polynomials
polynomial_order = 4;

% Compute properties along isentrope
isentrope_segments = create_barotropic_model(T_in, p_in, fluid, ...
                                             include_metastable=include_metastable, ...
                                             polynomial_order=polynomial_order, ...
                                             properties={'rhomass'});

% Plot barotropic process in p-s diagram
plot_barotropic_process_ps(fluid, isentrope_segments, ...
                           save_figures=save_figures, ...
                           plot_spinodal_line=true, ...
                           spinodal_line_method='robust', ...
                           plot_quality_isolines=true, ...
                           quality_labels=true, ...
                           show_in_legend=false)

% Plot barotropic model fitting
plot_barotropic_polynomials(fluid, isentrope_segments, save_figures=save_figures)

% Export fitted polynomials
export_fluent_expressions(fluid, isentrope_segments)


