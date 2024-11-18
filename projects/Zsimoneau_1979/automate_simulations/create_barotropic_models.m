% Clear the workspace
clear all %#ok
close all
clc

% Set plot options
set_plot_options()
save_figures = true;
defaultColors = get(groot, 'factoryAxesColorOrder');

% Create folder to save results
results_path = '../data_simulations';
if not(isfolder(results_path))
    mkdir(results_path)
end

% Load experimental cases
case_table = readtable("../cases_summary.xlsx", NumHeaderLines=0);
case_indices = [88:1:95 100:1:117];
case_indices = 112;
case_table = case_table(ismember(case_table.Case, case_indices), :);

% Create the fluid
fluid_nane = "nitrogen";
fluid = FluidCoolProp_2Phase("HEOS", fluid_nane);

% Compute critical point
fluid.abstractstate.update(py.CoolProp.DmassT_INPUTS, fluid.abstractstate.rhomass_critical, fluid.abstractstate.T_critical)
T_critical = fluid.abstractstate.T;
p_critical = fluid.abstractstate.p;
s_critical = fluid.abstractstate.smass;

% Compute triple
fluid.abstractstate.update(py.CoolProp.QT_INPUTS, 0.0, fluid.abstractstate.Ttriple)
p_triple = fluid.abstractstate.p;

% Barotropic model parameters
polynomial_order = 6;
include_metastable = false;
% prop_names = ["rhomass", "viscosity", "speed_sound"];
prop_names = ["rhomass"];

% Create figure
fig = figure(); ax = gca; hold on; box on; grid off
pbaspect([2.0, 1, 1])
ax.Title.String = {'Barotropic models for the cases from Simoneau Hendricks (1979)'; ""};
ax.XLabel.String = {''; "$s/s_{\mathrm{crit}}$ -- Reduced entropy"};
ax.YLabel.String = {"$p/p_{\mathrm{crit}}$ -- Reduced pressure"; ''};
xtickformat('%0.2f')
ytickformat('%0.2f')
ax.XLim = [0.925, 1.075];
ax.YLim = [0, 3];

% Plot phase diagram
plot_phase_diagram('smass', "p", fluid.abstractstate, ...
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
    name = case_table.Case(j);
    T_in = case_table.T_0_in(j);
    p_in = case_table.P_0_in(j);
    p_out = case_table.P_out(j);

    % Define case name
    case_name = strrep(sprintf('case_%03d_%s_T%0.2fK_P%0.2fbar', name, fluid.fluid_name, T_in, p_in/1e5), ".", "_");
    case_index = sprintf('case_%03d', name);
    fprintf('%s\n', case_name)

    % Create directory to save results
    output_folder = fullfile(results_path, case_name, "barotropic_model");
    if not(isfolder(output_folder))
        mkdir(output_folder)
    end

    % Compute properties along isentrope
    isentrope_segments = create_barotropic_model(T_in, p_in, fluid, ...
                                                 include_metastable=include_metastable, ...
                                                 polynomial_order=polynomial_order, ...
                                                 properties=prop_names, ...
                                                 p_high=max(2*p_critical, 1.5*p_in), ...
                                                 p_low=1.5*p_triple, ... 
                                                 N_points=100);
    
    % Plot the isentropic expansion
    for i = 1:numel(isentrope_segments)
        x = isentrope_segments(i).(prop_x);
        y = isentrope_segments(i).(prop_y);
        if i == 1
            text(ax, x(1), y(1), 0, ['$\;\;$',num2str(case_table.Case(j))], FontSize=5, HorizontalAlignment="left")
        end
        if j == 1
            plot(ax, x([1 end]), y([1 end]), Color=defaultColors(i,:), LineWidth=1.00, LineStyle='-', Marker='o', Markersize=3.5, MarkerFaceColor='w', HandleVisibility="off", DisplayName=isentrope_segments(i).label)
        else
            plot(ax, x([1 end]), y([1 end]), Color=defaultColors(i,:), LineWidth=1.00, LineStyle='-', Marker='o', Markersize=3.5, MarkerFaceColor='w', HandleVisibility="off", DisplayName=isentrope_segments(i).label)
        end
    end

    % Plot barotropic polynomials
    plot_barotropic_polynomials(fluid, isentrope_segments, case_index, output_dir=output_folder, save_figures=save_figures, show_figures=true)
    
    % Export fitted polynomials
    % Use case index instead of full case name because the file paths
    % become too long that the files could not be read with Python
    export_fluent_expressions(fluid, isentrope_segments, case_index, output_dir=output_folder)


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

