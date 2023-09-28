%% Verify model implementation with enthapy and entropy checks
% Clear the workspace
clear all
close all
clc

% Define plot settings
set_plot_options()
show_figures = true;
save_figures = false;
defaultColors = get(groot, 'factoryAxesColorOrder');

% Create folder to save results
output_dir = sprintf("figures/sensitivity_composition");
if not(isfolder(output_dir))
    mkdir(output_dir)
end

% Import CoolProp
py.importlib.import_module('CoolProp.CoolProp');

% Create separate fluid objects for water and nitrogen
fluid_name = "water_nitrogen";
fluid_label = "Water-Nitrogen";
fluid_1 = py.CoolProp.CoolProp.AbstractState('HEOS', 'Water');
fluid_2 = py.CoolProp.CoolProp.AbstractState('HEOS', 'Nitrogen');

% Load experimental data
exp_data = readtable("../data_experimental/experimental_cases.xlsx");
MixtureRatio = exp_data.mixture_ratio;

% Inlet and outlet conditions
T_in = 22 + 273.15;
p_in = 2000e3;
p_out = 98.6e3;

% Define efficiency for barotropic model
eta_poly = 1.00;

% Initialize figure handles
figs = struct();

% Define which properties to plot
properties_plot = {'temperature', 'density', 'viscosity', 'viscosity_kinematic', 'speed_sound', 'void_fraction', 'entropy', 'enthalpy'};

% Loop over the different mixture ratios
for i = 1:numel(MixtureRatio)

    % Mixture composition    
    R = MixtureRatio(i);
    y_1 = R./(1 + R);
    y_2 = 1./(1 + R);

    % Compute properties along polytropic process
    props = evaluate_barotropic_model_two_components(T_in, ...
                                                     p_in, ...
                                                     p_out, ...
                                                     y_1, ...
                                                     y_2, ...
                                                     fluid_1, ...
                                                     fluid_2, ...
                                                     eta_poly, ...
                                                     N_points=100, ...
                                                     p_min=p_out, ...
                                                     p_max=p_in);

    % Plot properties as a function of pressure
    alpha = 0.25 + (i-1)/(numel(MixtureRatio)-1) * (1.00-0.25);
    label = ['$R=', num2str(props.mixture_ratio(end), '%0.2f'), '$'];
    figs = plot_barotropic_model_two_components(props, ...
                                                alpha=alpha, ...
                                                figure_handles=figs, ...
                                                properties_plot=properties_plot, ...
                                                label=label, ...
                                                show_figures=show_figures);

end

% Save figures once all lines were added
prop_names = fieldnames(figs);
for k = 1:numel(prop_names)
    exportgraphics(figs.(prop_names{k}), fullfile(output_dir, ['sensitivity_mixture_ratio_', prop_names{k} ,'.png']), Resolution=500)
end

