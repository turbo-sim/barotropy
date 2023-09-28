%% Create barotropic models for CFD simulations
% Clear the workspace
clear all
close all
clc

% Define plot settings
set_plot_options()
show_figures = false;
save_figures = true;
defaultColors = get(groot, 'factoryAxesColorOrder');

% Create folder to save results
output_dir = fullfile("output");
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
case_data = readtable("simulation_cases.xlsx");
% case_data = case_data(1:4, :);

% Define properties to create barotropic models for
property_names = {'p'; 'density'; 'viscosity'; 'speed_sound'; 'void_fraction'; 'mass_fraction'};

% Polynomial order
order = 8;

% Loop over experimental cases
for i = 1:numel(case_data.index)

    % Mixture composition    
    R = case_data.mixture_ratio(i);
    y_1 = R./(1 + R);
    y_2 = 1./(1 + R);

    % Inlet and outlet conditions
    T_in = case_data.T0_in(i);
    p_in = case_data.p0_in(i);
    p_out = case_data.p_out(i);

    % Polytropic efficiency of the process
    eta_poly = 1.00;

    % Define case name based on index and identifier
    case_index = case_data.index(i);
    case_tag = case_data.tag{i};
    if isempty(case_tag)
        format_str = 'case_%03d';
        case_name = sprintf(format_str, case_index);
    else
        format_str = 'case_%03d_%s';
        case_name = sprintf(format_str, case_index, case_tag);
    end

    % Create folder to save results
    output_dir = fullfile(output_rootdir, case_name, 'barotropic_model');
    fig_dir = output_dir;
%     fig_dir = fullfile(output_dir, "figures");
    if not(isfolder(fig_dir))
        mkdir(fig_dir)
    end

    % Calculate thermodynamic properties
    props = evaluate_barotropic_model_two_components(T_in, ...
                                                     p_in, ...
                                                     p_out, ...
                                                     y_1, ...
                                                     y_2, ...
                                                     fluid_1, ...
                                                     fluid_2, ...
                                                     eta_poly, ...
                                                     N_points=100, ...
                                                     p_min=0.5*p_out, ...
                                                     p_max=2.0*p_in);

    % Create barotropic model
    barotropicModel = struct();
    barotropicModel.fluid.fluid_name = strrep(strrep(lower(fluid_name), ' ', '_'), '-', '_');
    barotropicModel.p_low = min(props.p);
    barotropicModel.p_high = max(props.p);
    barotropicModel.p_scaling = p_in;
    barotropicModel.property_names = property_names;
    for k = 1:numel(property_names)
        p_norm = props.p/p_in;
        values = props.(property_names{k});
        barotropicModel.property_values.(property_names{k}) = values;
        barotropicModel.property_values.label = fluid_label';
        barotropicModel.polynomial_coefficients.(property_names{k}) = polyfit(p_norm, values, order);
    end

    % Plot barotropic model fitting
    plot_barotropic_model(barotropicModel, ...
                          case_name=case_name, ...
                          output_dir=fig_dir, ...
                          save_figures=save_figures, ...
                          plot_extrapolation=true, ...
                          show_figures=show_figures)
    
    % Export barotropic model
    export_fluent_expressions(barotropicModel, ...
                              case_name=case_name, ...
                              output_dir=output_dir)


end






