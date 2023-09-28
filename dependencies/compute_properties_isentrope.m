function isentrope = compute_properties_isentrope(state_inlet, p_out, fluid, NameValueArgs)

    arguments
        state_inlet
        p_out
        fluid
        NameValueArgs.N_points (1, 1)
        NameValueArgs.force_helmholtz (1, 1) logical = false
    end

    % Retrieve inlet state
    p_in = state_inlet.p;
    s_in = state_inlet.smass;
    property_values = state_inlet;
    property_names = fieldnames(property_values);

    % Properties from state 1 to state 2
    p_array = linspace(p_in, p_out, NameValueArgs.N_points);
    for i = 1:numel(p_array)

        % Compute properties
        if NameValueArgs.force_helmholtz
            property_values = compute_properties_metastable_ps(p_array(i), s_in, fluid.abstractstate, property_values.T, property_values.rhomass);
        else
            fluid.set_prop_Ps(p_array(i), s_in)
            property_values = fluid.fluid_properties;
        end

        % Store properties
        for j = 1:numel(property_names)
            isentrope.(property_names{j})(i) = property_values.(property_names{j});
        end

    end

end
