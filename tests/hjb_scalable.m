% Set solver
solvers = ['sedumi', 'mosek', 'scs', 'sdpt3', 'sdpnal'];
mset('yalmip', true)

mset( ...
    sdpsettings( ...
    'solver', ...
    'scs', ...
    'verbose', 2, ...
    'showprogress', 1, ...
    'debug', 1) ...
    );

tic;
% Parameters

d = 4; % Relaxation degree
n = 3; % Dimension of problem
lmbda = 10; % Regularization parameter

% Interior occupation measure
mpol('x', n + 1);
mpol('z', n);
m = meas([x; z]);

% Boundary measures, defined as a vector of measures
mpol('xb', 2 * n * (n + 1));
mpol('zb', 2 * n^2);
mb = meas([xb; zb]);

% Determine support constraints
K = [];

% Support constraints for occupation measure
for i = 1:n
    K = [K, 0 <= x(i), x(i) <= 1];
end

disp(K);

% Support constraints for boundary measures
for i = 1:2 * n
    for j = 1:n
        ind = (i - 1) * (n + 1) + j;
        if j == ceil(i / 2)
            K = [K, xb(ind) == mod(i, 2)]; 
        else
            K = [K, 0 <= xb(ind), xb(ind) <= 1];
        end
    end
end

disp(K);

% Define interior test functions
phi_int = mmon(x, d);

% Determine moment constraints

M = [];
% Add constraint (C1)
for i = 1:n
    ind_1 = 2 * (i - 1) * (n + 1) + 1;
    ind_2 = ind_1 + n;
    phi_1 = mmon(xb(ind_1:ind_2), d);
    phi_2 = mmon(xb(ind_2 + 1:ind_2 + n + 1), d);
    M = [M, mom(diff(phi_int, x(i)) + diff(phi_int, x(n+1)) * z(i)) == 0];
end

% Add constraint (C2)
norm_term = x(1)^2;
int_term = mom((diff(phi_int, x(1))...
    + diff(phi_int, x(n+1)) * z(1)) * z(1));
bd_term = mom(-mmon(xb(1:n+1), d) * zb(1))...
    + mom(mmon(xb(n+2:2*n+2), d) * zb(n + 2));

for i = 2:n
    norm_term = norm_term + x(i)^2;
    int_term = int_term + mom((diff(phi_int, x(i))...
        + diff(phi_int, x(n+1)) * z(i)) * z(i));
end

first_ord_term = mom(phi_int * (x(1) - lmbda * norm_term));

M = [M, first_ord_term + bd_term == int_term];

% Add constraint (C3)
ft_ind = n + 2;
norm_sq = xb(ft_ind)^2;
phi_f = mmon(xb(ft_ind:ft_ind + n), d);
for i = 2:n
    norm_sq = norm_sq + xb(ft_ind + i)^2;
end
% M = [M, mom(phi_f * norm_sq * (1 - norm_sq)) == 0];

P = msdp(K, M);

[status, obj] = msol(P);

disp(["LMI " int2str(d) " lower bound = " num2str(obj)]);
toc;