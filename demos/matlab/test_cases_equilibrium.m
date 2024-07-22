%% Test the barotropic model for different isentropess
% Clear the workspace
clear all
close all
clc

% Add path to utility functions
% addpath(genpath("../dependencies"))

% Set plot options
set_plot_options()
save_figures = true;

% Define the fluid
fluid_name = 'CO2';
fluid = FluidCoolProp_2Phase('HEOS', fluid_name);
fluid.throw_exceptions = false;

% Define case parameters
polynomial_order = 4;
smoothing_factor = 0.0;
include_metastable = false;
prop_names = ["rhomass", "viscosity", "speed_sound"];
% prop_names = {'rhomass', 'speed_sound', 'viscosity', 'conductivity', 'cpmass'};

% Compute critical state
T_critical = fluid.abstractstate.T_critical;
p_critical = fluid.abstractstate.p_critical;
rho_critical = fluid.abstractstate.rhomass_critical;
fluid.abstractstate.update(py.CoolProp.DmassT_INPUTS, rho_critical, T_critical);
s_critical = fluid.abstractstate.smass;

% Compute triple pressure
fluid.abstractstate.update(py.CoolProp.QT_INPUTS, 0, fluid.abstractstate.Ttriple);
p_triple = fluid.abstractstate.p;

% Define the inlet states (entropy levels)
p_in = 2*p_critical;
s_frac_array = [0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6];

% Create figure
fig = figure(); ax = gca; hold on; box on; grid off
ax.Title.String = {sprintf("Barotropic process for %s", fluid.fluid_name); ''};
ax.XLabel.String = {''; "$s/s_{\mathrm{crit}}$ -- Reduced entropy"};
ax.YLabel.String = {"$p/p_{\mathrm{crit}}$ -- Reduced pressure"; ''};
xtickformat('%0.1f')
ytickformat('%0.1f')
% xlim([0.5, 2.0])
% ylim([0.5, 2.0])

% Plot phase diagram
prop_x = 'smass';
prop_y = 'p';
plot_phase_diagram(prop_x, prop_y, fluid.abstractstate, ...
                   plot_saturation_line=true, ...
                   plot_critical_point=true, ...
                   plot_spinodal_line=false, ...
                   spinodal_line_method='standard', ...
                   plot_quality_isolines=true, ...
                   show_in_legend=false)

% Plot barotropic model
for j = 1:numel(s_frac_array)
    
    % Define case name
    fprintf("Isentrope %i out of %i\n", j, numel(s_frac_array))
    case_name = strrep(sprintf('results/%s_equilibrium_%0.0fScrit', fluid.fluid_name, 100*s_frac_array(j)), ".", "_");
    if not(isfolder(case_name))
        mkdir(case_name)
    end

    % Compute inlet temperature
    s_in = s_critical*s_frac_array(j);
    fluid.abstractstate.update(py.CoolProp.PSmass_INPUTS, p_in, s_in);
    T_in = fluid.abstractstate.T;

    % Compute properties along isentrope assuming equilibirum
    isentrope_segments = create_barotropic_model(T_in, p_in, fluid, ...
                                                 include_metastable=include_metastable, ...
                                                 polynomial_order=polynomial_order, ...
                                                 properties=prop_names, ...
                                                 p_high=p_in, ...
                                                 p_low=p_triple, ...
                                                 N_points=100);
    
    % Plot the isentropic expansion
    defaultColors = get(groot, 'factoryAxesColorOrder');
    for i = 1:numel(isentrope_segments)
        x = isentrope_segments.property_values(i).(prop_x);
        y = isentrope_segments.property_values(i).(prop_y);
        label = isentrope_segments.property_values(i).label;
        if j == 1
            plot(ax, x([1 end]), y([1 end]), Color=defaultColors(i,:), LineWidth=1.00, LineStyle='-', Marker='o', Markersize=3.5, MarkerFaceColor='w', HandleVisibility="on", DisplayName=label)
        else
            plot(ax, x([1 end]), y([1 end]), Color=defaultColors(i,:), LineWidth=1.00, LineStyle='-', Marker='o', Markersize=3.5, MarkerFaceColor='w', HandleVisibility="off", DisplayName=label)
        end
    end

    % Plot barotropic polynomials
    plot_barotropic_polynomials(fluid, isentrope_segments, case_name, save_figures=save_figures, show_figures=false)


end

% Add legend
legend(Location="northeast")

% Rescale units
items = allchild(findall(fig, type='axes'));
for i = 1:numel(items)
    items(i).XData =  items(i).XData/s_critical;
    items(i).YData =  items(i).YData/p_critical;
end

% Adjust y-limits
items = allchild(findall(fig, type='axes'));
max_y = 0;
for i = 1:numel(items)
    max_y = max(max_y, max(items(i).YData, [], 'all'));
end
ax.XLim = [0.25 1.75];
ax.YLim = [0 max(1.10, 1.1*max_y)];

% Save p-s diagram
if save_figures
    exportgraphics(fig, fullfile('results', sprintf('%s_ps_diagram_equilibrium.png', fluid_name)), Resolution=500)
end

  
