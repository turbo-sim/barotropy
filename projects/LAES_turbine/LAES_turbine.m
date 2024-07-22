% Clear the workspace
clear all %#ok
close all
clc

% Set plot options
set_plot_options()
save_figures = true;
defaultColors = get(groot, 'factoryAxesColorOrder');

% Create folder to save results
results_path = './data_simulations';
if not(isfolder(results_path))
    mkdir(results_path)
end

% Load experimental cases
case_table = readtable("./simulation_cases.xlsx", NumHeaderLines=0);
case_indices = case_table.index;
case_table = case_table(ismember(case_table.index, case_indices), :);

% Create the fluid
fluid_nane = "nitrogen";
fluid = FluidCoolProp_2Phase("HEOS", fluid_nane);

% Compute critical point
fluid.abstractstate.update(py.CoolProp.DmassT_INPUTS, fluid.abstractstate.rhomass_critical, fluid.abstractstate.T_critical)
T_critical = fluid.abstractstate.T;
p_critical = fluid.abstractstate.p;
s_critical = fluid.abstractstate.smass;

% Compute triple point
fluid.abstractstate.update(py.CoolProp.QT_INPUTS, 0.0, fluid.abstractstate.Ttriple)
p_triple = fluid.abstractstate.p;

% Barotropic model parameters
polynomial_order = 6;
include_metastable = false;
prop_names = ["rhomass", "viscosity", "speed_sound"];
% prop_names = ["rhomass"];

% Create figure
fig = figure(); ax = gca; hold on; box on; grid off
pbaspect([1.0, 1, 1])
ax.Title.String = {'Barotropic model'; ""};
ax.XLabel.String = {''; "$s/s_{\mathrm{crit}}$ -- Reduced entropy"};
ax.YLabel.String = {"$p/p_{\mathrm{crit}}$ -- Reduced pressure"; ''};
xtickformat('%0.2f')
ytickformat('%0.2f')
ax.XLim = [0.6, 1.1];
ax.YLim = [0, 3];


prop_x = "smass";
prop_y = "p";

% Plot phase diagram
plot_phase_diagram(prop_x, prop_y, fluid.abstractstate, ...
                   plot_critical_point=true, ...
                   plot_saturation_line=true, ...
                   plot_spinodal_line=true, ...
                   spinodal_line_color=0.5*[1,1,1], ...
                   spinodal_line_width=1.25, ...
                   spinodal_line_method='robust', ...
                   plot_quality_isolines=true, ...
                   show_in_legend=true)

% Plot barotropic model
for j = 1:height(case_table)

    % Rename variables
    index = case_table.index(j);
    T_in = case_table.T0_in(j);
    p_in = case_table.p0_in(j);
    p_out = case_table.p_out(j);

    % Define case name
    case_name = sprintf('case_%03d', index);
    fprintf('%s\n', case_name)

    % Create directory to save results
    output_folder = fullfile(results_path, case_name, "barotropic_model");
    if not(isfolder(output_folder))
        mkdir(output_folder)
    end

    % Compute properties along isentrope
    barotropic_model = create_barotropic_model(T_in, p_in, fluid, ...
                                                 include_metastable=include_metastable, ...
                                                 polynomial_order=polynomial_order, ...
                                                 properties=prop_names, ...
                                                 p_high=max(2*p_critical, 1.5*p_in), ...
                                                 p_low=1.5*p_triple, ... 
                                                 N_points=200);
    
    % Plot the isentropic expansion
    for i = 1:numel(barotropic_model.property_values)
        x = barotropic_model.property_values(i).(prop_x);
        y = barotropic_model.property_values(i).(prop_y);
        label = barotropic_model.property_values(i).label;
        if i == 1
            text(ax, x(1), y(1), 0, ['$\;\;$',num2str(case_table.index(j))], FontSize=5, HorizontalAlignment="left")
        end
        if j == 1
            plot(ax, x([1 end]), y([1 end]), Color=defaultColors(i,:), LineWidth=1.00, LineStyle='-', Marker='o', Markersize=3.5, MarkerFaceColor='w', HandleVisibility="off", DisplayName=label)
        else
            plot(ax, x([1 end]), y([1 end]), Color=defaultColors(i,:), LineWidth=1.00, LineStyle='-', Marker='o', Markersize=3.5, MarkerFaceColor='w', HandleVisibility="off", DisplayName=label)
        end
    end

    % Plot barotropic polynomials
    plot_barotropic_model(barotropic_model, case_name=case_name, output_dir=output_folder, save_figures=save_figures, show_figures=true, plot_extrapolation=true, x_scale="log")
    
    % Export fitted polynomials
    export_fluent_expressions(barotropic_model, fluid_name=fluid_nane, case_name=case_name, output_dir=output_folder)
    export_cfx_expressions(barotropic_model, fluid_name=fluid_nane, case_name=case_name, output_dir=output_folder)


end

% Display legend
legend(ax, Location="northeast")

% Rescale units
items = allchild(findall(fig, type='axes'));
for i = 1:numel(items)

    if isequal(class(items(i)), 'matlab.graphics.chart.primitive.Line') || isequal(class(items(i)), 'matlab.graphics.chart.primitive.Contour') 
        items(i).XData =  items(i).XData/s_critical;
        items(i).YData =  items(i).YData/p_critical;

    elseif isequal(class(items(i)), 'matlab.graphics.primitive.Text')
        items(i).Position(1) = items(i).Position(1)/s_critical;
        items(i).Position(2) = items(i).Position(2)/p_critical;
    end

end

% Save p-s diagram
if save_figures
    exportgraphics(fig, fullfile(results_path, 'barotropic_models_ps.png'), Resolution=500)
    savefig(fullfile(results_path, 'barotropic_models_ps.fig'))
end

