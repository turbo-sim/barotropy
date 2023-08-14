function properties = compute_spinodal_point_entropy(s_in, fluid, NameValueArgs)

    % This function calculates the spinodal point corresponding to the
    % input entropy specified as argument.
    % 
    % The function verifies if the input entropy is within the limits of 
    % the spinodal line and throws an error if it is outside the range
    %
    % The calculation procedure is based on solving the nonlinear system:
    %
    %   1. s(T,rho) - s_in = 0
    %   2. isothermal_bulk_modulus(T,rho) / p(T,rho) = 0
    %
    % where the properties are evaluated using temperature-entropy function
    % calls to the Helmholtz energy equation of state.
    % 
    % The initial guess to solve the nonlinear system is obtained by 
    % precomputing the spinodal line and using the point with closest 
    % entropy as initial guess.
    %

    arguments
        s_in
        fluid
        NameValueArgs.T_margin (1, 1) double = 0.00
    end

    % Check that the input entropy is within limits
    [spinodal_liq, spinodal_vap] = compute_spinodal_line(fluid, N_points=50, method='standard');
    if s_in < spinodal_liq.smass(1)
        error('Input entropy is lower than the liquid spinodal entropy at the triple point: s_in=%0.2f, s_min=%0.2f', s_in, spinodal_liq.smass(1))
    end    
    if s_in > spinodal_vap.smass(end)
        error('Input entropy is higher than the vapor spinodal entropy at the triple point: s_in=%0.2f, s_max=%0.2f', s_in, spinodal_vap.smass(end))
    end
    
    % Estimate the initial guess
    s_spinodal = [spinodal_liq.smass, spinodal_vap.smass];
    T_spinodal = [spinodal_liq.T, spinodal_vap.T];
    rho_spinodal = [spinodal_liq.rhomass, spinodal_vap.rhomass];
    [~, index] = min(abs(s_in-s_spinodal));
    rhoT_guess = [rho_spinodal(index), T_spinodal(index)];

    % Solve residual equation for isothermal bulk modulus and entropy
    options = optimoptions('fsolve', ...
                           'Algorithm','trust-region-dogleg', ...
                           'MaxIterations', 2000, ...
                           'MaxFunctionEvaluations', 10000, ...
                           'FunctionTolerance', 1e-16, ...
                           'OptimalityTolerance', 1e-16, ...
                           'StepTolerance', 1e-16, ...
                           'FiniteDifferenceType','forward', ...
                           'Display', 'none');
    [x, ~, exitflag, ~]  = fsolve(@(rhoT_pair)get_spinodal_residual(rhoT_pair, s_in, fluid), rhoT_guess, options);
    if exitflag <= 0
        error('Spinodal point calculation did not converge')
    end

    % Compute thermodynamic properties at the spinodal point
    properties = compute_properties_metastable_Td(x(2)+NameValueArgs.T_margin, x(1), fluid);


end

function res = get_spinodal_residual(rhoT_pair, s_in, fluid)

    % Update fluid state
    props = compute_properties_metastable_Td(rhoT_pair(2), rhoT_pair(1), fluid);

    % Compute state residual
    % The bulk modulus is normalized by pressure to scale the problem
    % This is necessary to have a well-conditioned system of equations and
    % achieve tigh
    B_res = props.isothermal_bulk_modulus/props.p;
    s_res = (props.smass - s_in)/s_in;
    res = [B_res; s_res];

end
