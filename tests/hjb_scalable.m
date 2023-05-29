% Set solver
solvers = ['sedumi', 'mosek', 'scs', 'sdpt3', 'sdpnal'];
mset('yalmip', true)
mset(sdpsettings('solver', 'sdpt3', 'verbose', 1, 'showprogress', 1, 'debug', 1));

tic;
% Parameters

d = 2; % Relaxation degree
n = 2; % Dimension of problem
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
    M = [M, mom(diff(phi_int, x(i)) + diff(phi_int, x(n+1)) * z(i)) ...
            - mom(phi_1, d) - mom(-phi_2, d) == 0];
end

% Add constraint (C2)
norm_term = x(1)^2;
int_term = mom((diff(phi_int, x(1))...
    + diff(phi_int, x(n+1)) * z(1)) * z(1));
for i = 2:n
    norm_term = norm_term + x(i)^2;
    int_term = int_term + mom((diff(phi_int, x(i))...
        + diff(phi_int, x(n+1)) * z(i)) * z(i));
end
first_order_term = mom(phi_int * (x(1) - lmbda * norm_term));
bd_term = mom(-mmon(xb(1:n+1), d) * zb(1))...
    + mom(mmon(xb(n+2:2*n+2), d) * zb(n + 2));

for i = 2:2:2 * n
    disp(i);
end

% Add constraint (C3)
ft_ind = n + 2;
norm_sq = xb(ft_ind)^2;
phi_f = mmon(xb(ft_ind:ft_ind + n), d);
for i = 2:n
    norm_sq = norm_sq + xb(ft_ind + i)^2;
end
M = [M, mom(phi_f * norm_sq * (1 - norm_sq)) == 0];

% P = msdp(0 <= x(1), x(1) <= 1, 0 <= x(2), x(2) <= 1,...
%     0 <= x1(1), x1(1) <= 1, x1(2) == 0, ...
%     0 <= x2(1), x2(1) <= 1, x2(2) == 1, ...
%     0 <= x3(2), x3(2) <= 1, x3(1) == 0, ...
%     0 <= x4(2), x4(2) <= 1, x4(1) == 1, ...
%     mom(diff(phi_int, x(1)) + diff(phi_int, x(n+1)) * z(1))...
%         - mom(phi_1) + mom(-phi_2) == 0, ...
%     mom(diff(phi_int, x(2)) + diff(phi_int, x(n+1)) * z(2))...
%         - mom(phi_3) - mom(-phi_4) == 0,...
%     mom(phi_int * (x(1) - lmbda * x(1)^2 - lmbda * x(2)^2)) ...
%         + mom(phi_1 * z1(1)) + mom(-phi_2 * z2(2)) ...
%         + mom(phi_3 * z3(1)) + mom(-phi_4 * z4(2)) == ...
%         mom((diff(phi_int, x(1)) + diff(phi_int, x(n+1)) * z(1)) * z(1)) ...
%         + mom((diff(phi_int, x(2)) + diff(phi_int, x(n+1)) * z(2)) * z(2)),...
%     mom(phi_2 * x2(1)^2 * (1 + x2(1)^2)) == 0);


% % Measure over Boundary 1
% mpol('x1', n + 1);
% mpol('z1', n);
% m1 = meas([x1; z1]);
% 
% % Measure over Boundary 2
% mpol('x2', n + 1);
% mpol('z2', n);
% m2 = meas([x2; z2]);
% 
% % Measure over Boundary 3
% mpol('x3', n + 1);
% mpol('z3', n);
% m3 = meas([x3; z3]);
% 
% % Measure over Boundary 4
% mpol('x4', n + 1);
% mpol('z4', n);
% m4 = meas([x4; z4]);

% Support constraints
% K = [0 <= x(1), x(1) <= 1, 0 <= x(2), x(2) <= 1, ...
% 0 <= x1(1), x1(1) <= 1, x1(2) == 0, ...
% 0 <= x2(1), x2(1) <= 1, x2(2) == 1, ...
% 0 <= x3(2), x3(2) <= 1, x3(1) == 0, ...
% 0 <= x4(2), x4(2) <= 1, x4(1) == 1
% ];
% 
% % Interior test function
% phi_int = mmon(x, d);
% 
% % Boundary test functions
% phi_1 = mmon(x1, d);
% phi_2 = mmon(x2, d);
% phi_3 = mmon(x3, d);
% phi_4 = mmon(x4, d);
% 
% P = msdp(0 <= x(1), x(1) <= 1, 0 <= x(2), x(2) <= 1,...
%     0 <= x1(1), x1(1) <= 1, x1(2) == 0, ...
%     0 <= x2(1), x2(1) <= 1, x2(2) == 1, ...
%     0 <= x3(2), x3(2) <= 1, x3(1) == 0, ...
%     0 <= x4(2), x4(2) <= 1, x4(1) == 1, ...
%     mom(diff(phi_int, x(1)) + diff(phi_int, x(n+1)) * z(1)) == 0, ...
%     mom(diff(phi_int, x(2)) + diff(phi_int, x(n+1)) * z(2)) == 0, ...
%     mom(phi_int * (x(1) - lmbda * x(1)^2 - lmbda * x(2)^2)) == ...
%         mom((diff(phi_int, x(1)) + diff(phi_int, x(n+1)) * z(1)) * z(1)) ...
%         + mom((diff(phi_int, x(2)) + diff(phi_int, x(n+1)) * z(2)) * z(2)), ...
%     mom(phi_1 * (1 + x1(2)^2)) == 0 ...
%     );
% 
% [status, obj] = msol(P);
% 
% disp(["LMI " int2str(d) " lower bound = " num2str(obj)]);
% toc;