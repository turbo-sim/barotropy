%% Plot fluid properties for different mixture ratios
% Clear the workspace
clear all
close all
clc

% Define plot settings
set_plot_options()
show_figures = true;
save_figures = true;
defaultColors = get(groot, 'factoryAxesColorOrder');

% Create folder to save results
output_dir = sprintf("figures/model_verification");
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
eta_poly = 1.00;

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

% Plot the error in enthalpy and entropy
fig = figure(); hold on; box on;
xlabel({''; 'Pressure (bar)'})
ylabel({'Normalized numerical error', ''})
ax = gca;
ax.YScale = 'log';
ax.YLim = 10.^[-16, -10];
plot(props.p/1e5, abs((props.hmass - props.hmass_bis)./props.hmass_bis), Linewidth = 1.0, DisplayName='$\frac{|h_\mathrm{mix}-h_\mathrm{ode}|}{h_\mathrm{ode}}$')
plot(props.p/1e5, abs((props.smass - props.smass(1))./props.smass(1)), Linewidth = 1.0, DisplayName='$\frac{|s_\mathrm{mix}-s_\mathrm{in}|}{s_\mathrm{in}}$')
legend(Location="northeast", FontSize=15, NumColumns=2)

% Save figures once all lines were added
if save_figures
    exportgraphics(fig, fullfile(output_dir, 'barotropic_model_verification.png'), Resolution=500)
end

