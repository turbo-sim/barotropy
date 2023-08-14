function plot_barotropic_process_ps(fluid, isentrope_segments, case_name, NameValueArgs)

    arguments
        fluid
        isentrope_segments
        case_name (1, 1) string = sprintf('barotropic_model_%s', fluid.fluid_name)
        NameValueArgs.save_figures (1, 1) logical = false
        NameValueArgs.plot_saturation_line (1, 1) logical = true
        NameValueArgs.plot_critical_point (1, 1) logical = true
        NameValueArgs.plot_spinodal_line (1, 1) logical = false
        NameValueArgs.spinodal_line_method (1, 1) string = 'standard'
        NameValueArgs.plot_quality_isolines (1, 1) logical  = false
        NameValueArgs.plot_pseudocritical_line (1, 1) logical = false
        NameValueArgs.quality_levels (:, 1) double = 0.0:0.1:0.9
        NameValueArgs.quality_labels (:, 1) double = false
        NameValueArgs.show_in_legend (1, 1) logical = false
    end

    if not(isfolder(case_name))
        mkdir(case_name)
    end

    % Create figure
    fig = figure(); ax = gca; hold on; box on; grid off
    ax.Title.String = {sprintf("Barotropic process for %s", fluid.fluid_name); ''};
    ax.XLabel.String = {''; "$s/s_{\mathrm{crit}}$ -- Reduced entropy"};
    ax.YLabel.String = {"$p/p_{\mathrm{crit}}$ -- Reduced pressure"; ''};
    xtickformat('%0.1f')
    ytickformat('%0.1f')
    % xlim([0.5, 2.0])
    % ylim([0.5, 2.0])
    
    % Compute critical state
    T_critical = fluid.abstractstate.T_critical;
    p_critical = fluid.abstractstate.p_critical;
    fluid.abstractstate.update(py.CoolProp.PT_INPUTS, p_critical, T_critical);
    s_critical = fluid.abstractstate.smass;
   
    % Plot phase diagram
    prop_x = 'smass';
    prop_y = 'p';
    plot_phase_diagram(prop_x, prop_y, fluid.abstractstate, ...
                       plot_saturation_line=NameValueArgs.plot_saturation_line, ...
                       plot_critical_point=NameValueArgs.plot_critical_point, ...
                       plot_spinodal_line=NameValueArgs.plot_spinodal_line, ...
                       spinodal_line_method=NameValueArgs.spinodal_line_method, ...
                       plot_quality_isolines=NameValueArgs.plot_quality_isolines, ...
                       quality_levels=NameValueArgs.quality_levels, ...
                       quality_labels=NameValueArgs.quality_labels, ...
                       show_in_legend=NameValueArgs.show_in_legend)
    
    % Plot the isentropic expansion
    defaultColors = get(groot, 'factoryAxesColorOrder');
    for i = 1:numel(isentrope_segments)
        x = isentrope_segments(i).(prop_x);
        y = isentrope_segments(i).(prop_y);
%         plot(x, y, Color=defaultColors(i,:), HandleVisibility="off", DisplayName=isentrope_segments(i).label)
        plot(x([1 end]), y([1 end]), Color=defaultColors(i,:), LineWidth=1.00, LineStyle='-', Marker='o', Markersize=3.5, MarkerFaceColor='w', HandleVisibility="on", DisplayName=isentrope_segments(i).label)
    end
    
    % Rescale units
    items = allchild(fig.Children(:));
    for i = 1:numel(items)
        items(i).XData =  items(i).XData/s_critical;
        items(i).YData =  items(i).YData/p_critical;
    end

    % Adjust y-limits
    max_y = 0;
    for i = 1:numel(items)
        max_y = max(max_y, max(items(i).YData, [], 'all'));
    end
    ylim([0 max(1.10, 1.1*max_y)])
    
    % Display legend
    legend(Location="northeast", FontSize=9)
    
    % Save figure
    subname = split(case_name, '/'); subname = subname{end};
    if NameValueArgs.save_figures
        exportgraphics(fig, fullfile(case_name, sprintf('%s_ps_diagram.png', subname)), Resolution=500)
    end

end