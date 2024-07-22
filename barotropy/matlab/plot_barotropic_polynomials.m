function plot_barotropic_polynomials(fluid, isentrope_segments, case_name, NameValueArgs)

    arguments
        fluid
        isentrope_segments
        case_name (1, 1) string = sprintf('barotropic_model_%s', fluid.fluid_name)
        NameValueArgs.smoothing_factor (1, 1) double = 0.00
        NameValueArgs.save_figures (1, 1) logical = true
        NameValueArgs.show_figures (1, 1) logical = true
    end

    if not(isfolder(case_name))
        mkdir(case_name)
    end

    % Use correct syntax for figure visibility
    if NameValueArgs.show_figures
        show_figures = 'on';
    else 
        show_figures = 'off';
    end

    % Redefine property names
    prop_names = isentrope_segments.property_names;
    prop_names = setdiff(prop_names, ["p", "T", "smass"]);

    % Compute critical pressure
    p_critical = fluid.abstractstate.p_critical;

    % Plot thermodynamic properties along expansion
    for j = 1:numel(prop_names)
    
        % Get current property name
        prop_name = prop_names{j};
        
        % Plot property values
        fig = figure('Visible', show_figures);
        ax = subplot(2,1,1); hold on; box on; axis square
        pbaspect([2.5, 1, 1])
        defaultColors = get(groot, 'factoryAxesColorOrder');
        ax.Title.String = {sprintf("Isentropic expansion of %s", fluid.fluid_name); ''};
        ax.XLabel.String = {''; "Reduced pressure"};
        ax.YLabel.String = {[strrep(prop_name, '_', ' '),  ' value']; ''};
        xtickformat('%0.1f')
        ytickformat('%0.2f')
        for i = 1:numel(isentrope_segments)
            pressure = isentrope_segments(i).p;
            prop = isentrope_segments(i).(prop_name);
            plot(pressure/p_critical, prop, Color=defaultColors(i,:), DisplayName=isentrope_segments(i).label)
            plot(pressure([1 end])/p_critical, prop([1 end]), Color=defaultColors(i,:), LineStyle='none', Marker='o', Markersize=3.5, MarkerFaceColor='w', HandleVisibility="off")
        end
        
        % Evaluate polynomial at the same points
        p_poly = [];
        prop_value = [];
        for i = 1:numel(isentrope_segments)
            prop_value = [prop_value, isentrope_segments(i).(prop_name)];
            p_poly = [p_poly, isentrope_segments(i).p];
        end
    %     p_poly = linspace(0, 2*p_in/p_critical, 10000)*p_critical;
        prop_poly = evaluate_isentrope_polynomials(p_poly, prop_name, isentrope_segments, smoothing_factor=NameValueArgs.smoothing_factor);
        plot(ax, p_poly/p_critical, prop_poly, 'k:', DisplayName='Polynomial fit')
        legend(Location="southeast", FontSize=9)
    
        % Plot error distribution
        ax = subplot(2,1,2); hold on; box on; axis square
    %     xtickformat('%0.1f')
    %     ytickformat('%0.2f')
        pbaspect([2.5, 1, 1])
        ax.XLabel.String = {''; "Reduced pressure"};
        ax.YLabel.String = {[strrep(prop_name, '_', ' '), ' rel. error']; ''};
        ax.XScale = "linear";
        ax.YScale = "linear";
        plot(p_poly/p_critical, (prop_value-prop_poly)./prop_value, 'k', LineWidth=0.75)
    
        % Save the figures
        subname = split(case_name, '/'); subname = subname{end};
        if NameValueArgs.save_figures
            exportgraphics(fig, fullfile(case_name, sprintf('%s_%s_polynomial.png', subname, prop_name)), Resolution=500)
        end
         
    end

end


