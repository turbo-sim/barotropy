function prop = evaluate_isentrope_polynomials(pressure, property_name, isentrope_segments, NameValueArgs)

    arguments
        pressure (1, :) double
        property_name (1, 1) string
        isentrope_segments
        NameValueArgs.smoothing_factor (1, 1) double = 0.;
    end

    % Rescale the query pressure
    p_crit = isentrope_segments(1).p_critical;
    x = pressure/p_crit;

    % Combine all segments
    prop = 0; 
    for i = 1:numel(isentrope_segments)
        
        % Evaluate polyomials
        prop_segment = horner_polyval(isentrope_segments(i).(strcat('coeff_', property_name)), x);
 
        % Get the limit points of the interval
        x_1 = isentrope_segments(i).p(1)/p_crit;
        x_2 = isentrope_segments(i).p(end)/p_crit;

        % Piecewise polynomial evaluation
        if NameValueArgs.smoothing_factor == 0
    
            % Use in-line 'if' to evaluate the polynomials at the correct range
            prop = prop + prop_segment.*(x <= x_1).*(x >= x_2);

        % Blended polynomial evaluation (smooth and continuous)
        else
        
            % Use sigmoid function to blend the polynomials
            sigma =  (1 + tanh((x-x_1)/NameValueArgs.smoothing_factor))/2;
            if i == 1
                prop = prop_segment;
            else
                prop = (sigma).*prop + (1-sigma).*prop_segment;
            end           
            
    
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