function export_fluent_expressions(barotropic_model, NameValueArgs)

    arguments
        barotropic_model
        NameValueArgs.fluid_name (1, 1) string = barotropic_model.fluid.fluid_name
        NameValueArgs.case_name (1, 1) string = sprintf('barotropic_model_%s', barotropic_model.fluid.fluid_name)
        NameValueArgs.output_dir (1, 1) string = sprintf('barotropic_model_%s', barotropic_model.fluid.fluid_name)
    end

    if not(isfolder(NameValueArgs.output_dir))
        mkdir(NameValueArgs.output_dir)
    end

    % Rename properties
    prop_names = barotropic_model.property_names;
    prop_names = setdiff(prop_names, ["p", "T", "smass"]); % Remove properties from list
    prop_values = barotropic_model.property_values;
    prop_coeffs = barotropic_model.polynomial_coefficients;
    p_low = barotropic_model.p_low;
    p_scaling = barotropic_model.p_scaling;
    case_name = NameValueArgs.case_name;
    output_dir = NameValueArgs.output_dir;

    % Define the units of the fluid properties
    units = struct();
    units.rhomass = "kg/m^3";
    units.viscosity = "Pa*s";
    units.speed_sound = "m/s";
    units.conductivity = "W/m/K";
    units.cpmass = "J/kg/K";

    % Define additional variables
    polynomial_form = 'horner';

    % Write data to file in a try-catch environment
    fid = fopen(fullfile(output_dir, sprintf('%s_fluent_expressions.txt', case_name)),'wt');
    try

        % Create file header
        fprintf(fid, 'Fluent expressions of %s properties along an isentrope\n', NameValueArgs.fluid_name);
        fprintf(fid, 'Creation datetime: %s\n\n', datetime("now"));
        
        % Create expressions for each property
        for k = 1:numel(prop_names)
    
            % Create piecewise function recursively
            expression = polynomial_expression(prop_coeffs(1).(prop_names{k}), p_scaling, polynomial_form);
            for j = 2:numel(prop_values)
                expression_bis = polynomial_expression(prop_coeffs(j).(prop_names{k}), p_scaling, polynomial_form);
                expression = if_expression(expression, expression_bis, prop_values(j).p(1));
            end

            % Safeguard for negative pressures
            prop_low = polyval(prop_coeffs(end).(prop_names{k}), p_low/p_scaling);
            slope_low = polyval(polyder(prop_coeffs(end).(prop_names{k})), p_low/p_scaling);
            safeguard_expression = sprintf('%0.6e*exp((AbsolutePressure - %0.6e [Pa])/%0.6e [Pa])', prop_low, p_low, (slope_low/prop_low/p_scaling)^-1);
            expression = if_expression(expression, safeguard_expression, p_low);

            % Add fluid property unit
            expression = sprintf('%s * 1 [%s]', expression, units.(prop_names{k}));
            
            % Print expression
            fprintf(fid, '%s_%s\n', case_name, prop_names{k});
            fprintf(fid, '%s\n', expression);
            fprintf(fid, '\n');
    
        end

        % Close file
        fclose(fid);

    catch exception

        % If an error occurs, catch it and display the error message
        fclose(fid);
        error(['Error occurred while writing to file: ' exception.message]);

    end
    
end

function polynomial_string = polynomial_expression(coefficients, p_scaling, polynomial_form)

    arguments
        coefficients (:, 1) double
        p_scaling (1, 1) double
        polynomial_form (1, 1) string = "horner"
    end
    
    if strcmp(polynomial_form, "horner")

        % Create summation terms
        polynomial_string = sprintf('%+0.6e', coefficients(1));
        for i = 2:numel(coefficients)
            polynomial_string = sprintf('%+0.6e + (AbsolutePressure/%0.6e [Pa])*\n(%s)', coefficients(i), p_scaling, polynomial_string);
        end

    elseif strcmp(polynomial_form, "standard")
    
        % Create polynomial terms
        terms = cell(1, numel(coefficients));
        for i = 1:numel(coefficients)
            terms{i} = sprintf('%+0.6e*(AbsolutePressure/%0.6e [Pa])^%0.2f', coefficients(i), p_scaling, i-1);
        end
    
        % Concatenate the polynomial terms
        polynomial_string = strjoin(terms, ' \n');

    else
        error("The polynomial form must be 'horner' or 'standard'")
    end



end

function if_statement = if_expression(expression_1, expression_2, transition_pressure)

    if_statement = sprintf('IF(AbsolutePressure>=%s [Pa], \n%s, \n%s)', num2str(transition_pressure), expression_1, expression_2);

end

