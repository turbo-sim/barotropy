function isentrope = create_barotropic_model_segment(state_inlet, p_out, fluid, NameValueArgs)

    arguments
        state_inlet
        p_out
        fluid
        NameValueArgs.N_points (1, 1)
        NameValueArgs.polynomial_order (1, 1) double = 4
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

    % Fit polynomial to data
    p_reduced = isentrope.p/fluid.abstractstate.p_critical;
    for i = 1:numel(property_names)
        isentrope.(['coeff_', property_names{i}]) = polyfit(p_reduced, isentrope.(property_names{i}), NameValueArgs.polynomial_order);
    end

    % Store other variables in the structure
    isentrope.p_critical = fluid.abstractstate.p_critical;
    isentrope.molecular_weight = fluid.abstractstate.molar_mass;
    isentrope.property_names = property_names;

end
