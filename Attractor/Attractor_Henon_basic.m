% Global attractor approximation for the discrete-time Henon map.
% Corbinian Schlosser, Milan Korda, April 2020

close all
clear all

addpath(genpath('./functions'))


% Parameters
n = 2;        % Dimension of the state space
d = 8;        % Degree of polynomials v1, v2 and w
alpha = 0.05; % Discount factor

%  Dynamics
a = 1.4; b = 0.3;
Q = diag([1/1.5,1]); % Scaling
f_henon = @(x)( [  1 - a*x(1)^2 + x(2) ; b*x(1)  ] );
f_henon = @(x)(Q*f_henon(inv(Q)*x));

% Yalmip sdpvar representing the state
x = sdpvar(n,1);

% Yalmipize dynamics
f = f_henon(x);

% Constraint set X = [-1,1]^2
gx1 = 1 - x(1)^2;
gx2 = 1 - x(2)^2;

% Polynomials v1, v2, w
[v1, cv1] = polynomial(x,d);
[v2, cv2] = polynomial(x,d);
[w, cw] = polynomial(x,d);

% Sum of squares multipliers
[s1, c1] = polynomial(x,d-2);
[s2, c2] = polynomial(x,d-2);
[s3, c3] = polynomial(x,d-2);
[s4, c4] = polynomial(x,d-2);
[s5, c5] = polynomial(x,d-2);
[s6, c6] = polynomial(x,d-2);
[s7, c7] = polynomial(x,d-2);
[s8, c8] = polynomial(x,d-2);

% alpha * v1 \circ f <= v1 on X
con = [ sos(v1 - alpha*replace(v1,x,f) - gx1*s1 - gx2*s2);
    sos(s1) ; sos(s2) ];

% alpha * v1 \circ f <= v1 on X backward in time
con = [ con ; sos(replace(v2,x,f) - alpha*v2 - gx1*s3 - gx2*s4 );
    sos(s3) ; sos(s4) ];

% w >= V + 1 on X
con = [ con ; sos(w - v1 - v2 - 1 - gx1*s5 - gx2*s6 ) ;
    sos(s5) ; sos(s6) ];

% w >= 0 on X
con = [ con ; sos(w - gx1*s7 - gx2*s8  ) ;
    sos(s7) ; sos(s8)];

% Moments of the Lebesgue measure over X
leb_box = getLebesgueMomentsNew(d,[-1 -1 ; 1 1],1);

% Objective
obj = cw'*leb_box;

% Solver parameters
SDPsolver = lower('mosek');
sosmodel = 2; % 1 - kernel (faster), 2 - image representation (higher accuracy)
options = getSolverParams(SDPsolver,sosmodel);

% Solve
[sol,v,~,res] = solvesos(con,obj,options,[cw;cv1;cv2;c1;c2;c3;c4;c5;c6;c7;c8]);

fprintf('\nMax SOS equality constraint residual = %f \n', norm(res,Inf));

% Retrieve coefficients of the polynomials v1, v2
cv1 = double(cv1);
cv2 = double(cv2);



%% Plot

disp('Starting to plot...')

% Generate a trajectory for comparison
X = Q*[1;0];
for i = 1:40000
    X = [X f_henon(X(:,end))];
end

% Grid
grid = 1500; % 3000 was used in the paper. Takes a bit more time.
gridx = linspace(-1,1,grid); gridy = linspace(-1,1,grid);
[X1,X2] = meshgrid(gridx,gridy);

% Evaluate on the grid
v_vals = polyval_grid2D({cv1,cv2},X1,X2);

% min(v1,v2)
plot_val = min(v_vals{1},v_vals{2});

% Plot
[~,] = contourf(X1,X2,plot_val,[0 0],'linecolor','none'); hold on
colormap([ 0.7*[1 1 1] ]);
h = scatter(X(1,2000:end),X(2,2000:end),'MarkerFaceColor',[0.5 0 0],'MarkerEdgeColor','none'); hold on
h.SizeData = 10;
axis([-1,1,-0.5,0.5])