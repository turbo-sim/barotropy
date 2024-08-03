%% Sensitivity of barotropic model to polytropic efficiency
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
output_dir = sprintf("figures/sensitivity_efficiency");
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

% Inlet and outlet conditions
T_in = 22 + 273.15;
p_in = 2000e3;
p_out = 98.6e3;

% Define mixture composition
R = 67.8841;
y_1 = R./(1 + R);
y_2 = 1./(1 + R);

% Define efficiency for barotropic model
eta_poly = linspace(0, 1, 6);

% Initialize figure handles
fig_handles_abs = struct();
fig_handles_rel = struct();

% Define which properties to plot
properties_plot = {'temperature', 'density', 'viscosity', 'speed_sound', 'void_fraction', 'entropy', 'enthalpy'};
% properties_plot = {'temperature', 'density', 'viscosity'};

% Define figure dimensions
desired_width = 640;
desired_height = desired_width * (10/16);

% Compute isentropic properties
props_isentropic = evaluate_barotropic_model_two_components(T_in, ...
                                                            p_in, ...
                                                            p_out, ...
                                                            y_1, ...
                                                            y_2, ...
                                                            fluid_1, ...
                                                            fluid_2, ...
                                                            1.00, ...
                                                            N_points=100, ...
                                                            p_min=p_out, ...
                                                            p_max=p_in);

% Loop over the different mixture ratios
for i = 1:numel(eta_poly)

    % Compute properties along polytropic process
    props = evaluate_barotropic_model_two_components(T_in, ...
                                                     p_in, ...
                                                     p_out, ...
                                                     y_1, ...
                                                     y_2, ...
                                                     fluid_1, ...
                                                     fluid_2, ...
                                                     eta_poly(i), ...
                                                     N_points=100, ...
                                                     p_min=p_out, ...
                                                     p_max=p_in);

    % Plot properties as a function of pressure
    alpha = 0.25 + (i-1)/(numel(eta_poly)-1) * (1.00-0.25);
    blended_1 = alpha * defaultColors(1, :) + (1 - alpha) * [1 1 1];
    label = ['$\eta_p=', num2str(eta_poly(i)*100, '%0.0f'), '$\%'];
    fig_handles_abs = plot_barotropic_model_two_components(props, ...
                                                alpha=alpha, ...
                                                figure_handles=fig_handles_abs, ...
                                                properties_plot=properties_plot, ...
                                                label=label, ...
                                                show_figures=show_figures);

    
    % Plot relative density change
    var = 'density';
    if ismember(var, properties_plot)
        if isfield(fig_handles_rel, var)  % Existing figure
            figure(fig_handles_rel.(var));
        else  % New figure
            fig_handles_rel.(var) = figure(); hold on; box on;
            fig_handles_rel.(var).Position(3) = desired_width;
            fig_handles_rel.(var).Position(4) = desired_height;
            xlabel({''; 'Pressure (bar)'})
            ylabel({'Relative density change (\%)', ''})
            ytickformat('%0.2f')
            ylim([-0.5, 0.1])
        end
        plot(props.p/1e5, (props.rhomass./props_isentropic.rhomass-1)*100, Color=blended_1, Linewidth = 1.00, DisplayName=label)
%         legend(location='eastoutside', NumColumns=1)
        legend(location='southeast', NumColumns=1)
    end

    % Plot relative viscosity change
    var = 'viscosity';
    if ismember(var, properties_plot)
        if isfield(fig_handles_rel, var)  % Existing figure
            figure(fig_handles_rel.(var));
        else  % New figure
            fig_handles_rel.(var) = figure(); hold on; box on;
            fig_handles_rel.(var).Position(3) = desired_width;
            fig_handles_rel.(var).Position(4) = desired_height;
            xlabel({''; 'Pressure (bar)'})
            ylabel({'Relative viscosity change (\%)', ''})
            ytickformat('%0.2f')
        end
        plot(props.p/1e5, (props.viscosity./props_isentropic.viscosity-1)*100, Color=blended_1, Linewidth = 1.00, DisplayName=label)
%         legend(location='eastoutside', NumColumns=1)
        legend(location='southeast', NumColumns=1)
    end
 

end


% Save figures once all lines were added
prop_names = fieldnames(fig_handles_abs);
for k = 1:numel(prop_names)
    exportgraphics(fig_handles_abs.(prop_names{k}), fullfile(output_dir, ['sensitivity_efficiency_', prop_names{k} ,'.png']), Resolution=500)
    export_fig(fig_handles_abs.(prop_names{k}), fullfile(output_dir, ['sensitivity_efficiency_', prop_names{k} ,'.svg']), '-svg');
end

prop_names = fieldnames(fig_handles_rel);
for k = 1:numel(prop_names)
    exportgraphics(fig_handles_rel.(prop_names{k}), fullfile(output_dir, ['sensitivity_efficiency_relative_', prop_names{k} ,'.png']), Resolution=500)
    export_fig(fig_handles_rel.(prop_names{k}), fullfile(output_dir, ['sensitivity_efficiency_relative_', prop_names{k} ,'.svg']), '-svg');
end
