%% Plot the spinodal line in the pressure-temperature plane
% Clear workspace
clear all
close all
clc

% Create folder to save results
results_path = '../';
if not(isfolder(results_path))
    mkdir(results_path)
end

% Set plot options
set_plot_options()
save_figures = true;
defaultColors = get(groot, 'factoryAxesColorOrder');

% Load experimental cases
cases = readtable("../cases_summary.xlsx", NumHeaderLines=1);

% Create the fluid
fluid_nane = "nitrogen";
fluid = FluidCoolProp_2Phase("HEOS", fluid_nane);
prop_x = 'smass';
prop_y = "p";

% Compute critical point
fluid.abstractstate.update(py.CoolProp.DmassT_INPUTS, fluid.abstractstate.rhomass_critical, fluid.abstractstate.T_critical)
T_critical = fluid.abstractstate.T;
p_critical = fluid.abstractstate.p;
s_critical = fluid.abstractstate.smass;

% Create figure
fig = figure(); ax = gca; hold on; box on; grid off;
% pbaspect([1.4, 1, 1])
% ax.Layer="top";
ax.Title.String = {'Nitrogen nozzle cases from Simoneau and Hendricks (1979)'; ""};
ax.XLabel.String = {''; "$s/s_\mathrm{crit}$ -- Reduced entropy"};
ax.YLabel.String = {"$p/p_\mathrm{crit}$ -- Reduced pressure"; ' '};
ax.XScale = "linear";
ax.YScale = "linear";
xtickformat('%.2f');
ytickformat('%.2f')
ax.XLim = [0.7, 1.1];
ax.YLim = [0, 2.2];


% Plot the phase diagram
plot_phase_diagram(prop_x, prop_y, fluid.abstractstate, ...
                   plot_critical_point=true, ...
                   plot_saturation_line=true, ...
                   plot_spinodal_line=true, ...
                   spinodal_line_color='red', ...
                   spinodal_line_width=1.25, ...
                   spinodal_line_method='robust', ...
                   plot_quality_isolines=true)

% Plot boundary conditions of each case
for i = 1:numel(cases.P_0_in)

    % Inlet state
    fluid.set_prop_PT(cases.P_0_in(i), cases.T_0_in(i))
    props_in = fluid.fluid_properties;

    % Outlet isentropic state
    fluid.set_prop_Ps(cases.P_out(i), props_in.smass)
    props_out = fluid.fluid_properties;

    % Plot states
    plot([props_in.(prop_x), props_out.(prop_x)], [props_in.(prop_y), props_out.(prop_y)], 'k-', LineWidth=0.1, Marker='o', MarkerSize=1.5, MarkerFaceColor='w')
    text(props_in.(prop_x), props_in.(prop_y), 0, ['$\quad$', num2str(cases.Case(i))], FontSize=4)

end



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


if save_figures
    exportgraphics(fig, fullfile(results_path, sprintf('cases_summary.png')), Resolution=500)
%     savefig(fig, fullfile(results_path, 'cases_summary.fig'))
end