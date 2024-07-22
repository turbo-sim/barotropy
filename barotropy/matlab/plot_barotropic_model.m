function plot_barotropic_model(barotropic_model, NameValueArgs)

    arguments
        barotropic_model
        NameValueArgs.fluid_name (1, 1) string = barotropic_model.fluid.fluid_name
        NameValueArgs.case_name (1, 1) string = sprintf('barotropic_model_%s', barotropic_model.fluid.fluid_name)
        NameValueArgs.output_dir (1, 1) string = sprintf('barotropic_model_%s', barotropic_model.fluid.fluid_name)
        NameValueArgs.save_figures (1, 1) logical = true
        NameValueArgs.show_figures (1, 1) logical = true
        NameValueArgs.plot_extrapolation (1, 1) logical = false
    end

    if not(isfolder(NameValueArgs.output_dir))
        mkdir(NameValueArgs.output_dir)
    end

    % Use correct syntax for figure visibility
    if NameValueArgs.show_figures
        show_figures = 'on';
    else 
        show_figures = 'off';
    end

    % Rename variables
    p_low = barotropic_model.p_low;
    prop_names = barotropic_model.property_names;
    prop_names = setdiff(prop_names, ["p", "T", "smass"]); % Remove items
    prop_values = barotropic_model.property_values;
    prop_coeffs = barotropic_model.polynomial_coefficients;
    case_name = NameValueArgs.case_name;
    output_dir = NameValueArgs.output_dir;

    % Get scaling pressure
    p_scaling = barotropic_model.p_scaling;

    % Plot thermodynamic properties along expansion
    for k = 1:numel(prop_names)
    
        % Get current property name
        prop_name = prop_names{k};
        
        % Create plot for property values
        fig = figure('Visible', show_figures);
        ax_1 = subplot(2,1,1); hold on; box on; axis square
        pbaspect([2.5, 1, 1])
        defaultColors = get(groot, 'factoryAxesColorOrder');
        ax_1.Title.String = {strrep(sprintf("Isentropic expansion of %s", NameValueArgs.fluid_name), '_', ' '); ''};
        ax_1.XLabel.String = {''; "Reduced pressure"};
        ax_1.YLabel.String = {[strrep(prop_name, '_', ' '),  ' value']; ''};
        xtickformat('%0.1f')
        ytickformat('%0.2f')

        % Create plot for error distribution
        ax_2 = subplot(2,1,2); hold on; box on; axis square
        pbaspect([2.5, 1, 1])
        ax_2.XLabel.String = {''; "Reduced pressure"};
        ax_2.YLabel.String = {[strrep(prop_name, '_', ' '), ' rel. error']; ''};
        ax_2.XScale = "linear";
        ax_2.YScale = "linear";

        % Plot barotropic model
        for i = 1:numel(prop_values)
            p_norm = prop_values(i).p/p_scaling;
            property = prop_values(i).(prop_name);
            property_poly = horner_polyval(prop_coeffs(i).(prop_name), p_norm);
            plot(ax_1, p_norm, property, Color=defaultColors(i,:), DisplayName=prop_values(i).label)
            plot(ax_1, p_norm([1 end]), property([1 end]), Color=defaultColors(i,:), LineStyle='none', Marker='o', Markersize=3.5, MarkerFaceColor='w', HandleVisibility="off")      
            if i == 1
                h = plot(ax_1, p_norm, property_poly, 'k:', DisplayName='Polynomial fit');
            else
                plot(ax_1, p_norm, property_poly, 'k:', HandleVisibility='off')
            end
            plot(ax_2, p_norm, (property-property_poly)./property, 'k', LineWidth=0.75)
        end
        uistack(h,'top');

        if NameValueArgs.plot_extrapolation
            prop_low = polyval(prop_coeffs(end).(prop_names{k}), p_low/p_scaling);
            slope_low = polyval(polyder(prop_coeffs(end).(prop_names{k})), p_low/p_scaling);
            pressure_extrapolation = linspace(-0.0, p_low/p_scaling, 100);
            property_extrapolation = prop_low*exp((pressure_extrapolation-p_low/p_scaling)*slope_low/prop_low);
            plot(ax_1, pressure_extrapolation, property_extrapolation, 'r-', DisplayName='Extrapolation')
        end
        legend(ax_1, Location="southeast", FontSize=9)
    
        % Save the figures
        if NameValueArgs.save_figures
            exportgraphics(fig, fullfile(output_dir, sprintf('%s_%s_polynomial.png', case_name, prop_name)), Resolution=500)
            savefig(fig, fullfile(output_dir, sprintf('%s_%s_polynomial.fig', case_name, prop_name)))
        end
         
    end

end


function y = horner_polyval(coefficients, x)

    % Evaluate polynomial using Horner's rule
    y = coefficients(1);
    for i = 2:numel(coefficients)
        y = coefficients(i) + y.*x;
    end

end

