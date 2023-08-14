function export_fluent_expressions(fluid, isentrope_segments, case_name)

    arguments
        fluid
        isentrope_segments
        case_name (1, 1) string = sprintf('barotropic_model_%s', fluid.fluid_name)
    end

    if not(isfolder(case_name))
        mkdir(case_name)
    end

    % Redefine property names
    prop_names = isentrope_segments.property_names;
    prop_names = setdiff(prop_names, ["p", "T", "smass"]);

    % Define the units of the fluid properties
    units = struct();
    units.rhomass = "kg/m^3";
    units.viscosity = "Pa*s";
    units.speed_sound = "m/s";
    units.conductivity = "W/m/K";
    units.cpmass = "J/kg/K";

    % Define additional variables
    polynomial_form = 'horner';
    p_critical = isentrope_segments(1).p_critical;
    p_low = isentrope_segments(end).p(end);

    % Write data to file in a try-catch environment
    subname = split(case_name, '/'); subname = subname{end};
    fid = fopen(fullfile(case_name, sprintf('%s_fluent_expressions.txt', subname)),'wt');
    try

        % Create file header
        fprintf(fid, 'Fluent expressions of %s properties along an isentrope\n', fluid.fluid_name);
        fprintf(fid, 'Creation datetime: %s\n\n', datetime("now"));
        
        % Create expressions for each property
        for k = 1:numel(prop_names)
    
            % Create piecewise function recursively
            expression = polynomial_expression(isentrope_segments(1), prop_names{k}, polynomial_form);
            for j = 2:numel(isentrope_segments)
                expression_bis = polynomial_expression(isentrope_segments(j), prop_names{k}, polynomial_form);
                expression = if_expression(expression, expression_bis, isentrope_segments(j).p(1));
            end
    
            % Safeguard for negative pressures
            prop_lowest = evaluate_isentrope_polynomials(p_low, prop_names{k}, isentrope_segments, smoothing_factor=0);
            safeguard_expression = sprintf('%0.6e*exp((AbsolutePressure - %0.6e [Pa])/%0.6e [Pa])', prop_lowest, p_low, p_critical);
            expression = if_expression(expression, safeguard_expression, p_low);
    
            % Add fluid property unit
            expression = sprintf('%s * 1 [%s]', expression, units.(prop_names{k}));
            
            % Print expression
            fprintf(fid, '%s_%s\n', prop_names{k}, case_name);
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

function polynomial_string = polynomial_expression(isentrope_segment, prop_name, polynomial_form)

    arguments
        isentrope_segment
        prop_name (1, 1) string
        polynomial_form (1, 1) string = "horner"
    end

    % Get critical pressure for scaling
    p_critical = isentrope_segment.p_critical;
    
    if strcmp(polynomial_form, "horner")

        % Get polynomial coefficients as array
        coefficients = isentrope_segment.(strcat('coeff_', prop_name));

        % Create summation terms
        polynomial_string = sprintf('%+0.6e', coefficients(1));
        for i = 2:numel(coefficients)
            polynomial_string = sprintf('%+0.6e + (AbsolutePressure/%0.6e [Pa])*\n(%s)', coefficients(i), p_critical, polynomial_string);
        end

    elseif strcmp(polynomial_form, "standard")

        % Get polynomial coefficients as array
        % Flip because Matlab decreasing exponent polynomial order
        coefficients = flip(isentrope_segment.(strcat('coeff_', prop_name)));
    
        % Get critical pressure for scaling
        p_critical = isentrope_segment.p_critical;
    
        % Create polynomial terms
        terms = cell(1, numel(coefficients));
        for i = 1:numel(coefficients)
            terms{i} = sprintf('%+0.6e*(AbsolutePressure/%0.6e [Pa])^%0.2f', coefficients(i), p_critical, i-1);
        end
    
        % Concatenate the polynomial terms
        polynomial_string = strjoin(terms, ' \n');

    else
        error("The polynomial form must be 'horner' or 'standard'")
    end



end

function if_statement = if_expression(expression_1, expression_2, transition_pressure)

    if_statement = sprintf('IF(AbsolutePressure>=%s [Pa], \n%s, \n%s)', transition_pressure, expression_1, expression_2);

end

