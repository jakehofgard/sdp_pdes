% Implementation of energy minimization optimal control problem in
% dimension two.
solvers = ['sedumi', 'mosek', 'scs', 'sdpt3', 'sdpnal'];
mset('yalmip', true)

% Set SCS solver settings
alphas = linspace(1.02, 1.98, 49);
trials = zeros(49, 2);
for i = 1:49
    % Primal scale factor is 1e-6 (better convergence than default 1e-3)
    scaling = 1e-6;
    alpha = alphas(i);
    mset(sdpsettings( ...
        'solver', 'scs', ...
        'verbose', 1, ...
        'showprogress', 1, ...
        'debug', 1, ...
        'savesolveroutput', 1, ...
        'scs.alpha', alpha, ... % DR relaxation. Default value is 1.5, recommended between 1-2.
        'scs.rho_x', scaling, ... % Primal scale factor. Default value is 1e-3.
        'scs.scale', 5.0) ... % Dual scale factor. YALMIP Default factor is 5.0. 
    );
    
    % Get total runtime
    tic;
    
    % Parameters
    
    d = 6; % Relaxation degree
    n = 2; % Dimension of problem
    lmbda = 10; % Regularization parameter
    
    % Surface normal vector
    nu = repmat([1, -1], 1, n);
    
    % Interior occupation measure
    mpol('x', n + 1);
    mpol('z', n);
    m = meas([x; z]);
    
    % Measure over Boundary 1
    mpol('x1', n + 1);
    mpol('z1', n);
    m1 = meas([x1; z1]);
    
    % Measure over Boundary 2
    mpol('x2', n + 1);
    mpol('z2', n);
    m2 = meas([x2; z2]);
    
    % Measure over Boundary 3
    mpol('x3', n + 1);
    mpol('z3', n);
    m3 = meas([x3; z3]);
    
    % Measure over Boundary 4
    mpol('x4', n + 1);
    mpol('z4', n);
    m4 = meas([x4; z4]);
    
    % Interior test function
    phi_int = mmon(x, d);
    
    % Boundary test functions
    phi_1 = mmon(x1, d);
    phi_2 = mmon(x2, d);
    phi_3 = mmon(x3, d);
    phi_4 = mmon(x4, d);
    
    
    % Compile support and moment constraints
    P = msdp(min(mom(x(3)^2 / 2)), 0 <= x(1), x(1) <= 1, 0 <= x(2), x(2) <= 1,...
        0 <= x1(1), x1(1) <= 1, x1(2) == 0, ...
        0 <= x2(1), x2(1) <= 1, x2(2) == 1, ...
        0 <= x3(2), x3(2) <= 1, x3(1) == 0, ...
        0 <= x4(2), x4(2) <= 1, x4(1) == 1, ...
        mom(diff(phi_int, x(1)) + diff(phi_int, x(n+1)) * z(1))...
            - mom(phi_1) + mom(-phi_2) == 0, ...
        mom(diff(phi_int, x(2)) + diff(phi_int, x(n+1)) * z(2))...
            - mom(phi_3) - mom(-phi_4) == 0,...
        mom(phi_int * (x(1) - lmbda * x(1)^2 - lmbda * x(2)^2)) ...
            + mom(phi_1 * z1(1)) + mom(-phi_2 * z2(2)) ...
            + mom(phi_3 * z3(1)) + mom(-phi_4 * z4(2)) == ...
            mom((diff(phi_int, x(1)) + diff(phi_int, x(n+1)) * z(1)) * z(1)) ...
            + mom((diff(phi_int, x(2)) + diff(phi_int, x(n+1)) * z(2)) * z(2)),...
        mom(phi_2 * (1 + x2(1)^2)) == 0);
    
    % Save YALMIP output
    [F, obj, s] = myalmip(P);
    [model, recoverymodel] = export(F, obj, sdpsettings('solver', 'scs'));
    
    % Solve problem using predefined solver
    [status, obj] = msol(P);
    
    % Rescale true objective according to primal scale factor
    true_obj = obj / scaling;
    trials(i, :) = [alpha, true_obj];
    writematrix(trials, "outputs/hjb_energy_low_dim_d=10.csv");

    % Display resclaed objective value
    disp(["Relaxation degree " int2str(d) ": lower bound = " num2str(true_obj)]);
    toc;
end