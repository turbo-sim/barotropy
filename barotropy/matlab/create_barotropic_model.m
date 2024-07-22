function barotropic_model = create_barotropic_model(T_in, p_in, fluid, NameValueArgs)

    arguments
        T_in (1, 1) double
        p_in (1, 1) double
        fluid

        NameValueArgs.N_points (1, 1) = 100
        NameValueArgs.polynomial_order (1, 1) double = 4
        NameValueArgs.p_low (1, 1) double = -1
        NameValueArgs.p_high (1, 1) double = -1
        NameValueArgs.p_scale (1, 1) double = -1 
        NameValueArgs.properties (1, :) string = ["rhomass", "cpmass", "speed_sound", "viscosity", "conductivity"]
        NameValueArgs.include_metastable (1, 1) logical = false
        NameValueArgs.spinodal_point_method (1, 1) string = 'robust'
    end

    
    %% Preliminary definitions
    % Define list of property names
    prop_names = [["p", "T", "smass"], NameValueArgs.properties];

    % Compute triple and critical pressures
    fluid.abstractstate.update(py.CoolProp.QT_INPUTS, 0, fluid.abstractstate.Ttriple);
    p_triple = fluid.abstractstate.p;
    p_critical = fluid.abstractstate.p_critical;
        
    % Define the low pressure limit
    if NameValueArgs.p_low ~= -1
        p_low = NameValueArgs.p_low;
    else
        p_low = p_triple;        
    end

    % Define the high pressure limit
    if NameValueArgs.p_high ~= -1
        p_high = NameValueArgs.p_high;
    else
        p_high = p_in;
    end

    % Define the scaling pressure
    if NameValueArgs.p_scale ~= -1
        p_scale = NameValueArgs.p_scale;
    else
        p_scale = p_critical;        
    end


    %% Compute inlet/saturation/spinodal states
    % Compute inlet state (posssibly at a higher pressure than p_in)
    fluid.abstractstate.update(py.CoolProp.PT_INPUTS, p_in, T_in)
    s_in = fluid.abstractstate.smass;
    fluid.set_prop_Ps(p_high, s_in)
    for i = 1:numel(prop_names)
        state_inlet.(prop_names{i}) = fluid.fluid_properties.(prop_names{i});
    end

    % Compute the p-s limits of the saturation line
    [s1_sat, ~, s2_sat, ~] = get_saturation_line_ps_limits(fluid.abstractstate);

    % Compute the p-s limits of the spinodal line
    if NameValueArgs.include_metastable
        [s1_spinodal, ~, s2_spinodal, ~] = get_spinodal_line_ps_limits(fluid.abstractstate);
    end

    % Compute saturation point if entropy is within limits
    if (s_in >= s1_sat) && (s_in <= s2_sat)
        saturation_props = compute_saturation_point_entropy(s_in, fluid.abstractstate);
        fluid.set_prop_rho_T(saturation_props.rhomass, saturation_props.T)
        for i = 1:numel(prop_names)
            state_saturation.(prop_names{i}) = fluid.fluid_properties.(prop_names{i});
        end
    end

    % Compute spinodal point if entropy is within limits
    if NameValueArgs.include_metastable
        if (s_in >= s1_spinodal) && (s_in <= s2_spinodal)
            spinodal_props = compute_spinodal_point_entropy(s_in, fluid.abstractstate, method=NameValueArgs.spinodal_point_method);
            for i = 1:numel(prop_names)
                state_spinodal.(prop_names{i}) = spinodal_props.(prop_names{i});
            end
        end
    end

    
    %% Compute isentrope segments excluding metastable states
    if ~NameValueArgs.include_metastable
        if s_in <= s1_sat
            error(['The inlet entropy must be higher than the liquid entropy at the triple point\n' ...
                   '\t s_in = %0.3f J/kg/K\n' ...
                   '\t s_triple_liq = %0.3f J/kg/K\n'], s_in, s1_sat)
    
        elseif (s_in > s1_sat) && (s_in <= s2_sat)
    
            % Properties from inlet state to saturation
            % Segment ends above the saturation line to prevent common points
            p_sat_plus = (1+1e-6)*state_saturation.p;
            property_segments(1) = compute_properties_isentrope(state_inlet, p_sat_plus, fluid, ...
                N_points=NameValueArgs.N_points, force_helmholtz=true);
        
            % Properties from saturation to the lowest pressure
            property_segments(2) = compute_properties_isentrope(state_saturation, p_low, fluid, ...
                N_points=NameValueArgs.N_points, force_helmholtz=false);
            
            % Give names to each segment
            property_segments(1).label = 'Single-phase region';
            property_segments(2).label = 'Two-phase region';
    
        else % (s_in > s2_sat)
    
            % Properties from inlet state to the lowest pressure
            property_segments(1) = compute_properties_isentrope(state_inlet, p_low, fluid, ...
                N_points=NameValueArgs.N_points, force_helmholtz=true);
        
            % Give names to each segment
            property_segments(1).label = 'Single-phase region';
    
        end
    end


    %% Compute the isentrope segments including metastable states
    if NameValueArgs.include_metastable
        if s_in <= s1_sat
            error(['The inlet entropy must be higher than the liquid entropy at the triple point\n' ...
                   '\t s_in = %0.3f J/kg/K\n' ...
                   '\t s_triple_liq = %0.3f J/kg/K\n'], s_in, s1_sat)
    
        elseif (s_in > s1_sat) && (s_in <= s1_spinodal)
    
            % Properties from inlet state to saturation
            % Segment ends above the saturation line to prevent common points
            p_sat_plus = (1+1e-6)*state_saturation.p;
            property_segments(1) = compute_properties_isentrope(state_inlet, p_sat_plus, fluid, ...
                N_points=NameValueArgs.N_points, force_helmholtz=true);
        
            % Properties from saturation to the lowest pressure
            % Segment ends above the spinodal line to prevent common points
            property_segments(2) = compute_properties_isentrope(state_saturation, p_low, fluid, ...
                N_points=NameValueArgs.N_points, force_helmholtz=true);
        
            % Give names to each segment
            property_segments(1).label = 'Single-phase region';
            property_segments(2).label = 'Metastable region';
    
    
        elseif (s_in > s1_spinodal) && (s_in <= s2_spinodal)
    
            % Properties from inlet state to saturation
            % Segment ends above the saturation line to prevent common points
            p_sat_plus = (1+1e-6)*state_saturation.p;
            property_segments(1) = compute_properties_isentrope(state_inlet, p_sat_plus, fluid, ...
                N_points=NameValueArgs.N_points, force_helmholtz=true);
        
            % Properties from saturation to spinodal point
            % Segment ends above the spinodal line to prevent common points
            p_spinodal_plus = (1+1e-6)*state_spinodal.p;
            property_segments(2) = compute_properties_isentrope(state_saturation, p_spinodal_plus, fluid, ...
                N_points=NameValueArgs.N_points, force_helmholtz=true);
        
            % Properties from spinodal to the lowest pressure
            property_segments(3) = compute_properties_isentrope(state_spinodal, p_low, fluid, ...
                N_points=NameValueArgs.N_points, force_helmholtz=false);
            
            % Give names to each segment
            property_segments(1).label = 'Single-phase region';
            property_segments(2).label = 'Metastable region';
            property_segments(3).label = 'Two-phase region';
    
        elseif (s_in > s2_spinodal) && (s_in <= s2_sat)
    
            % Properties from inlet state to saturation
            % Segment ends above the saturation line to prevent common points
            p_sat_plus = (1+1e-6)*state_saturation.p;
            property_segments(1) = compute_properties_isentrope(state_inlet, p_sat_plus, fluid, ...
                N_points=NameValueArgs.N_points, force_helmholtz=true);
        
            % Properties from saturation to the lowest pressure
            property_segments(2) = compute_properties_isentrope(state_saturation, p_low, fluid, ...
                N_points=NameValueArgs.N_points, force_helmholtz=true);
        
            % Give names to each segment
            property_segments(1).label = 'Single-phase region';
            property_segments(2).label = 'Metastable region';
    
    
        else % (s_in > s2_sat)
    
            % Properties from inlet state to the lowest pressure
            property_segments(1) = compute_properties_isentrope(state_inlet, p_low, fluid, ...
                N_points=NameValueArgs.N_points, force_helmholtz=true);
        
            % Give names to each segment
            property_segments(1).label = 'Single-phase region';
    
        end

    end

    % Fit polynomials to data
    coeffs(numel(property_segments)) = struct();
    order = NameValueArgs.polynomial_order;
    for j = 1:numel(property_segments)
        p_norm = property_segments(j).p/p_scale;
        for i = 1:numel(prop_names)
            values = property_segments(j).(prop_names{i});
            coeffs(j).(prop_names{i}) = polyfit(p_norm, values, order);
        end
    end

    % Export barotropic model
    barotropic_model = struct();
    barotropic_model.fluid = fluid;
    barotropic_model.p_low = p_low;
    barotropic_model.p_high = p_high;
    barotropic_model.p_scaling = p_scale;
    barotropic_model.property_names = prop_names;
    barotropic_model.property_values = property_segments;
    barotropic_model.polynomial_coefficients = coeffs;


end


function [s1, p1, s2, p2] = get_saturation_line_ps_limits(fluid)

    % Triple point liquid entropy
    fluid.update(py.CoolProp.QT_INPUTS, 0.00, fluid.Ttriple);
    s1 = fluid.smass;
    p1 = fluid.p;
    
    % Triple point vapor entropy
    fluid.update(py.CoolProp.QT_INPUTS, 1.00, fluid.Ttriple);
    s2 = fluid.smass;
    p2 = fluid.p;

end


function [s1, p1, s2, p2] = get_spinodal_line_ps_limits(fluid)


    % Compute the spinodal lines
    [spinodal_liq, spinodal_vap] = compute_spinodal_line(fluid, N_points=100, method='robust');

    % Get triple pressure
    fluid.update(py.CoolProp.QT_INPUTS, 0.00, fluid.Ttriple);
    p_triple = fluid.p;

    if spinodal_liq.p(1) <= p_triple
        [~, index] = min(abs(p_triple - spinodal_liq.p));
        s1 = spinodal_liq.smass(index);
        p1 = spinodal_liq.p(index);
    else
        s1 = spinodal_liq.smass(1);
        p1 = spinodal_liq.p(1);
    end

    if spinodal_vap.p(end) <= p_triple
        [~, index] = min(abs(p_triple - spinodal_vap.p));
        s2 = spinodal_vap.smass(index);
        p2 = spinodal_vap.p(index);
    else
        s2 = spinodal_vap.smass(end);
        p2 = spinodal_vap.p(end);
    end

end


function y = horner_polyval(coefficients, x)

    % Evaluate polynomial using Horner's rule
    y = coefficients(1);
    for i = 2:numel(coefficients)
        y = coefficients(i) + y.*x;
    end

end
    
