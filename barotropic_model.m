%% Generate barotropic fluid property model
% This script evaluate the thermodynamic properties of a fluid along a line
% of constant entropy. The properties between the saturation line and the
% spinodal line can be computed in 2 different ways:
%  
%   1. Assuming phase equilibrium
%   2. Extrapolating the equation of state into the metastable region
%
% Depending on the case, there can be up to 3 isentrope segments:
%   
%   1. Single-phase region
%   2. Metastable region
%   3. Two-phase region (equilibrium)
%
% Each of the isentrope segments is fitted with a polynomial with the form:
% 
%    y = sum(a_i*(p/p_crit)^i) with i = 0, 1,..., n
%
% Where n is the order of the polynomial specified by the user and the
% independent variable is the reduced pressure (i.e., p/p_crit) to improve 
% the conditioning of the system.
%
% The piece-wise polynomials can be exported as Fluent Expressions that
% to be used in conjuction with the pressure-based solver. The expressions
% represent the polynomials in Horner's form to reduce the computational 
% cost and the improve numerical stability of the model.
%
%
% Author: Roberto Agromayor
% Date: 14.08.2023
%

% Clear the workspace
clear all %#ok
close all
clc

% Add path to utility functions
% addpath(genpath("dependencies"))

% Set plot options
set_plot_options()
save_figures = false;

% Define the case parameters
fluid_name = 'CO2';
fluid = FluidCoolProp_2Phase('HEOS', fluid_name);
T_in = 300;
p_in = 90e5;
polynomial_order = 6;
include_metastable = true;
prop_names = ["rhomass", "viscosity"];
% prop_names = {'rhomass', 'speed_sound', 'viscosity', 'conductivity', 'cpmass'};

% Define case name
case_name = strrep(sprintf('results/%s_P%0.2fbar_T%0.2fK', fluid.fluid_name, p_in/1e5, T_in), ".", "_");
if not(isfolder(case_name))
    mkdir(case_name)
end

% Compute properties along isentrope
isentrope_segments = create_barotropic_model(T_in, p_in, fluid, ...
                                             include_metastable=include_metastable, ...
                                             polynomial_order=polynomial_order, ...
                                             properties=prop_names);

% Plot barotropic process in p-s diagram
plot_barotropic_process_ps(fluid, isentrope_segments, case_name, ...
                           save_figures=save_figures, ...
                           plot_spinodal_line=true, ...
                           spinodal_line_method='robust', ...
                           show_in_legend=false)

% Plot barotropic model fitting
plot_barotropic_polynomials(fluid, isentrope_segments, case_name, ...
                            save_figures=save_figures, ...
                            smoothing_factor=0.00)

% Export fitted polynomials
export_fluent_expressions(case_name, isentrope_segments)


