%% Initialize script
% Clear the workspace
clear all
close all
clc

% Create folder to save results
results_path = 'results';
if not(isfolder(results_path))
    mkdir(results_path)
end

% Define plot settings
set_plot_options()
save_figures = false;
defaultColors = get(groot, 'factoryAxesColorOrder');


%% Compute mixture properties
% Import CoolProp
py.importlib.import_module('CoolProp.CoolProp');

% Create separate fluid objects for water and nitrogen
fluid = py.CoolProp.CoolProp.AbstractState('HEOS', 'Water');

% figure(); hold on;
% plot_phase_diagram("rhomass", "gibbsmass", fluid, show_in_legend=true, plot_quality_isolines=true, plot_spinodal_line=true)


figure(); hold on;
plot_phase_diagram("p", "gibbsmass", fluid, show_in_legend=true, plot_quality_isolines=true, plot_spinodal_line=false)


% figure(); hold on;
% plot_phase_diagram("p", "gibbsmass", fluid, show_in_legend=true, plot_quality_isolines=true, plot_saturation_line=false)

T = 300+273.15;
fluid.update(py.CoolProp.QT_INPUTS, 0, T)
d_liq = fluid.rhomass;
fluid.update(py.CoolProp.QT_INPUTS, 1, T)
d_vap = fluid.rhomass;
d_array = linspace(0.5*d_vap, 1.5*d_vap, 50);
d_array = linspace(0.99*d_liq, 1.01*d_liq, 50);


for i = 1:numel(d_array)
    
    % Computation according to Helmholtz energy equation of state
    props = compute_properties_metastable_Td(T, d_array(i), fluid);
    p_metastable(i) = props.p;
    s_metastable(i) = props.smass;
    g_metastable(i) = props.gibbsmass;

end

d_array = linspace(0.5*d_vap, 1.01*d_liq, 50);

for i = 1:numel(d_array)
    % Computation according to Coolprop including 2-phase equilibrium
    fluid.update(py.CoolProp.CoolProp.DmassT_INPUTS, d_array(i), T)
    p_equilibrium(i) = fluid.p;
    s_equilibrium(i) = fluid.smass;
    g_equilibrium(i) = fluid.gibbsmass;

end

plot(p_metastable, g_metastable, 'b+-')
plot(p_equilibrium, g_equilibrium, 'r+-')


fluid.update(py.CoolProp.QT_INPUTS, 0, T)
p_sat = fluid.p

fluid.update(py.CoolProp.PT_INPUTS, p_sat*1.2, T)
g_liquid = fluid.gibbsmass;
p_liquid = fluid.p
plot(p_liquid, g_liquid, 'bo')


% p = linspace(0.5, 1.5, 100)*1e7;
% 
% for i = 1:numel(p)
% 
%     % Computation according to Coolprop including 2-phase equilibrium
%     fluid.update(py.CoolProp.CoolProp.PT_INPUTS, p(i), T)
%     p_equilibrium(i) = fluid.p;
%     g_equilibrium(i) = fluid.gibbsmass;
% 
% end
% 
% plot(p_equilibrium, g_equilibrium, 'ko')

% T = 300;
% fluid.update(py.CoolProp.QT_INPUTS, 1, T);
% p_sat = fluid.p;
% 
% p = linspace(0.5, 0.99999, 100)*p_sat;
% for i = 1:numel(p)
% 
%     fluid.update(py.CoolProp.PT_INPUTS, p(i), T)
%     gibbs(i) = fluid.gibbsmass;
% end
% plot(T+0*p, gibbs, 'r-+')
% 
% 
% p = linspace(1.00001, 10.5, 100)*p_sat;
% for i = 1:numel(p)
% 
%     fluid.update(py.CoolProp.PT_INPUTS, p(i), T)
%     gibbs(i) = fluid.gibbsmass;
% end
% plot((T)+0*p, gibbs, 'b-+')


xlabel("pressure")
ylabel("gibbs energy")



% %% Plot spinodal line on thermodynamic diagrams
% % Plot the spinodal line on the pressure-density diagram
% fig_1 = figure(); ax = gca; hold on; box on; grid off
% ax.Title.String = sprintf("Maxwell loop and spinodal points for %s", fluid_name);
% ax.XLabel.String = "Density (kg/m$^3$)";
% ax.YLabel.String = "Pressure (bar)";
% ax.XScale = "linear";
% ax.YScale = "linear";
% % ylim([0, 2*p_critical]/1e5)
% plot(sat_vap.rhomass, sat_vap.gibbsmass, 'r', DisplayName='Vapor saturation')
% plot(spinodal_vap.rhomass, spinodal_vap.gibbsmass, 'r:', DisplayName='Vapor spinodal')
% plot(sat_liq.rhomass, sat_liq.gibbsmass, 'b', DisplayName='Liquid saturation')
% plot(spinodal_liq.rhomass, spinodal_liq.gibbsmass, 'b:', DisplayName='Liquid spinodal')
% plot(d_array, g_equilibrium, 'k:', DisplayName='Equilibrium isotherm')
% plot(d_array, g_metastable, 'k-', DisplayName='Helmholtz EoS isotherm')
% plot(spinodal_point_vap.rhomass, spinodal_point_vap.gibbsmass, Color='black', Marker='o', MarkerSize=4, MarkerFaceColor='white', HandleVisibility='off')
% plot(spinodal_point_liq.rhomass, spinodal_point_liq.gibbsmass, Color='black', Marker='o', MarkerSize=4, MarkerFaceColor='white', HandleVisibility='off')
% legend(Location='northwest', FontSize=9)


% T = 100+273.15; 
% Q = 0.0;
% fluid.update(py.CoolProp.QT_INPUTS, Q, T)
% fluid.gibbsmass
% % fluid_1.gibbsmolar
% 
% 
% T = 100+273.15; 
% Q = 1;
% fluid.update(py.CoolProp.QT_INPUTS, Q, T)
% fluid.gibbsmass
% % fluid_1.gibbsmolar
% 
% 
% T = 100+273.15; 
% Q = 0.5;
% fluid.update(py.CoolProp.QT_INPUTS, Q, T)
% fluid.gibbsmass
% % fluid_1.gibbsmolar
% 
% 
% p = 1e5;
% Q = 0.0;
% fluid.update(py.CoolProp.PQ_INPUTS, p, Q)
% s_in = 0.2*fluid.smass;
% 
% p = 1e5;
% Q = 1.0;
% fluid.update(py.CoolProp.PQ_INPUTS, p, Q)
% s_out = 1.2*fluid.smass;
% 
% 
% s = linspace(s_in, s_out, 50)
% 
% for i = 1:numel(s)
% 
%     fluid.update(py.CoolProp.PSmass_INPUTS, p, s(i))
%     gibbs(i) = fluid.gibbsmass;
% 
% end
% 
% 
% figure(); hold on; box on;
% 
% plot(s, gibbs, 'k+')