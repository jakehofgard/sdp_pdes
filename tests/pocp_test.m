% system variables
mpol x 2; % state vector of size 2
mpol u; % scalar input
% state and input matrices
A = [0 1; 0 0];
B = [0; 1];
% define optimal control problem
prob = pocp( ...
'state', x, ...
'input', u, ...
'dynamics', A*x+B*u);
% set further properties
prob = set(prob, 'idirac', x, [1; 1], ... % initial condition
'fdirac', x, [0; 0], ... % final condition
'tconstraint', [x(2)>=-1; u>=-1; u<=1], ... % constraints
'scost', 1); % setting integral cost h to 1
% call to the solver
[status, cost] = solvepocp(prob, 14); % 14 = degree of moment matrix
disp('The lower bound on the minimum time is');
disp(cost);