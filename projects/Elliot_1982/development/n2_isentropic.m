%% Initialize script
% Clear the workspace
clear all
close all
clc

% Define plot settings
set_plot_options()
save_figures = true;
defaultColors = get(groot, 'factoryAxesColorOrder');

% Create folder to save results
output_dir = fullfile("output", "validation_plots");
if not(isfolder(output_dir))
    mkdir(output_dir)
end

%% Validation at different mixture ratios
% Import CoolProp
py.importlib.import_module('CoolProp.CoolProp');

% Create separate fluid objects for water and nitrogen
fluid_name = "water_nitrogen";
fluid = py.CoolProp.CoolProp.AbstractState('HEOS', 'Nitrogen');

% Inlet and outlet conditions
T_in = 22 + 273.15;
p_in = 2000e3;
p_out = 98.6e3;
A_throat = 135*1e-6;
A_exit = 597*1e-6;

cases = table();

cases.T0_in = [94.7; 96.6];
cases.p0_in = [780; 1560]*1e3;
cases.p_out = [210; 960]*1e3;

for i = 1:numel(cases.T0_in)

    fluid.update(py.CoolProp.PT_INPUTS, cases.p0_in(i), cases.T0_in(i))
    cases.s_in(i) = fluid.smass;
    cases.h0_in(i) = fluid.hmass;
    fluid.update(py.CoolProp.PSmass_INPUTS, cases.p_out(i), cases.s_in(i))
    cases.h_out(i) = fluid.hmass;
    cases.v_out(i) = sqrt(2*(cases.h0_in(i) - cases.h_out(i)));


end

cases
% % Load CFD data
% data_cfd_1 = readtable("simulation_cases.xlsx");
% data_cfd_2 = readtable("output/simulation_results.xlsx");
% common_vars = intersect(data_cfd_1.Properties.VariableNames, data_cfd_2.Properties.VariableNames);
% data_cfd_1(:, common_vars) = [];
% data_cfd = [data_cfd_1, data_cfd_2];
% 
% % Load experimental data
% % exp_data = readtable("experimental_data.xlsx");
% data_exp = readtable("../experimental_data/experimental_cases.xlsx");
% 
% % Compute isentropic quantities
% data_isentropic = struct();
% for i = 1:numel(data_exp.index)
% 
%     % Mixture composition    
%     R = data_exp.mixture_ratio(i);
%     y_1 = R./(1 + R);
%     y_2 = 1./(1 + R);
% 
%     % Isentropic calculations
%     eta_poly = 1.00;
%     props = evaluate_barotropic_model_two_components(data_exp.T0_in(i), ...
%                                                           data_exp.p0_in(i), ...
%                                                           data_exp.p_out(i), ...
%                                                           y_1, ...
%                                                           y_2, ...
%                                                           fluid_1, ...
%                                                           fluid, ...
%                                                           eta_poly, ...
%                                                           N_points=100, ...
%                                                           p_min=p_out, ...
%                                                           p_max=p_in);
% 
%     % Compute the isentropic velocity
%     h_in = props.hmass(1);
%     h_out_s = props.hmass(end);
%     data_isentropic.mixture_ratio(i, 1) = R;
%     data_isentropic.density(i, 1) = props.rhomass(end);
%     data_isentropic.density_liquid(i, 1) = props.rhomass_1(end);
%     data_isentropic.density_gas(i, 1) = props.rhomass_2(end);
%     data_isentropic.void_fraction(i, 1) = props.void_fraction(end);
%     data_isentropic.speed_sound(i, 1) = props.speed_sound(end);
%     data_isentropic.exit_velocity(i, 1) = sqrt(2*(h_in-h_out_s));
%     data_isentropic.exit_mach(i, 1) = data_isentropic.exit_velocity(i)/data_isentropic.speed_sound(i);
%     data_isentropic.isentropic_mach(i, 1) = data_isentropic.exit_velocity(i)/data_isentropic.speed_sound(i);
% 
% end
% 
% 
% data_isentropic.exit_velocity
% 
% % %%
% % % Postprocess experimental data
% % data_exp.nozzle_efficiency = (data_exp.exit_velocity./data_isentropic.exit_velocity).^2;
% % data_exp.exit_mach = (data_exp.exit_velocity./data_isentropic.speed_sound);
% % 
% % % Postprocess CFD data to account for exit pressure thrust
% % % (same wrong approach as done in the experiments)
% % data_cfd.thrust_velocity = data_cfd.mass_flow_rate_inlet.*data_cfd.velocity_outlet;
% % data_cfd.thrust_pressure = data_cfd.area_outlet.*(data_cfd.pressure_near_outlet-data_cfd.pressure_outlet);
% % data_cfd.thrust_total = data_cfd.thrust_velocity +data_cfd.thrust_pressure;
% % data_cfd.apparent_velocity = data_cfd.thrust_total./data_cfd.mass_flow_rate_inlet;
% % tags = unique(data_cfd.tag);
% % for i = 1:numel(tags)
% %     idx = strcmp(data_cfd.tag, tags{i});
% %     data_cfd.isentropic_velocity(idx) = data_isentropic.exit_velocity;
% %     data_cfd.mixture_ratio(idx) = data_isentropic.mixture_ratio;
% % end
% % data_cfd.nozzle_efficiency = (data_cfd.apparent_velocity./data_cfd.isentropic_velocity).^2;
% % 
% % 
% % %% Create figures
% % % Exit velocity (apparent)
% % fig_1 = figure(); ax_1 = gca; hold on; box on;
% % xlabel({''; "Mixture ratio ($y_\mathrm{water}/y_\mathrm{nitrogen}$)"})
% % ylabel({"Nozzle exit velocity (m/s)", ''})
% % ax_1.YLim = [50, 200];
% % 
% % % Isentropic efficiency
% % fig_2 = figure(); ax_2 = gca; hold on; box on;
% % xlabel({''; "Mixture ratio ($y_\mathrm{water}/y_\mathrm{nitrogen}$)"})
% % ylabel({"Nozzle isentropic efficiency (\%)", ''})
% % ax_2.YLim = [50, 100];
% % 
% % % Plot experimental data
% % plot(ax_1, data_exp.mixture_ratio, data_exp.exit_velocity, 'o-', markerfacecolor='w', color=defaultColors(2,:), DisplayName="Experimental data")
% % plot(ax_2, data_exp.mixture_ratio, data_exp.nozzle_efficiency*100, 'o-', markerfacecolor='w', color=defaultColors(2,:), DisplayName="Experimental data")
% % 
% % % Plot isentropic data
% % plot(ax_1, data_isentropic.mixture_ratio, data_isentropic.exit_velocity, color=defaultColors(1,:), DisplayName="Isentropic calculation")
% % 
% % % Plot CFD data
% % tags = unique(data_cfd.tag);
% % tags = {'a'; 'drag'}
% % labels = {'Mixture viscosity'; 'Liquid viscosity'}
% % for i = 1:length(tags)
% %     alpha = 0.25 + (i-1)/(numel(tags)-1) * (1.00-0.25);
% %     blended_1 = alpha * [0,0,0] + (1 - alpha) * [1 1 1];
% %     idx = strcmp(data_cfd.tag, tags{i});
% %     roughness = data_cfd.wall_roughness_height(idx);
% % %     label = ['CFD with $k_s=',num2str(1e6*roughness(1), '%0.2f'), '$ $\mu$m'];
% %     label = labels{i};
% %     plot(ax_1, data_cfd.mixture_ratio(idx), data_cfd.apparent_velocity(idx), linestyle="-", marker='o', markersize=2, color=blended_1, DisplayName=label)
% %     plot(ax_2, data_cfd.mixture_ratio(idx), data_cfd.nozzle_efficiency(idx)*100, linestyle="-", marker='o', markersize=2, color=blended_1, DisplayName=label)
% % end
% % 
% % 
% % % Add legends
% % legend(ax_1, location="northeast", NumColumns=1, FontSize=10)
% % legend(ax_2, location="southwest", NumColumns=1, FontSize=10)
% % 
% % % Add a text box in the bottom right corner of the figure
% % str = '$\eta_s = \frac{h_\mathrm{0,in}-h_\mathrm{out}}{h_\mathrm{0,in}-h_{\mathrm{out},s}} = \left(\frac{v_\mathrm{out}}{v_{\mathrm{out},s}}\right)^2$';
% % h = annotation('textbox', 'String', str, ...
% %                'EdgeColor', 'k', 'BackgroundColor', [1, 1, 1], ...
% %                'HorizontalAlignment', 'center', ...
% %                'VerticalAlignment', 'middle', ...
% %                'FitBoxToText', 'on');
% % pos = h.Position;
% % pos(1) = 0.66;
% % pos(2) = 0.15;
% % h.Position = pos;
% % 
% % % Save figures
% % if save_figures
% %     exportgraphics(fig_1, fullfile(output_dir, 'validation_exit_velocity_vs_mixture_ratio.png'), Resolution=500)
% %     exportgraphics(fig_2, fullfile(output_dir, 'validation_nozzle_efficiency_vs_mixture_ratio.png'), Resolution=500)
% % end
% % 
% % 
% % 
% % 
% % % % The isentropic and real properties almost do not change as is the case
% % % % for liquids (which is the majority of the mass)
% % % % 
% % % % - density
% % % % - viscosity
% % % % - speed of sound
% % % % 
% % % % for gases the properties might change a bit due to friction, but the
% % % % effect is not really important in most cases
% % % % 
% % % % Friction and viscous dissipation has the most important effect on
% % % % reducing the speed of the flow (and therefore the mass flow rate?)
% % % % 
% % % 
% % % 
% % % % I have to check if the formula mdot = A*rho*v checks out for the CFD
% % % % simulations. Make several planes at different x locations?
% % % 
% % % % % Slip correction factor estimation
% % % % % velocity_ratio = 4.5;
% % % % % A_out = 0.0005969887023908;
% % % % % mass_flow = data_cfd.mass_flow_rate_inlet;
% % % % % rho_slip = data.DensityLiquid.*data.VolFractionLiquid + velocity_ratio*data.DensityGas.*data.VolFractionGas;
% % % % % v_out = mass_flow/A_out./data.Density;
% % % % % v_out_corrected = (y_1 + velocity_ratio*y_2)*mass_flow/A_out./rho_slip;
% % % % % scaling = v_out./v_out_corrected
% % % % % data_cfd.velocity_out_mas
% % % % % s_avg = data_cfd.velocity_out_mass_avg./scaling;
% % % % 
% % 
