% Global attractor approximation for the discrete-time Henon map. Takes
% intersection of the outer approximations for several values of the degree and the discount factor

% Corbinian Schlosser, Milan Korda, April 2020

close all
clear all

addpath(genpath('./functions'))


n = 2; % Dimension of the state

%  Dynamics
a = 1.4; b = 0.3;
Q = diag([1/1.5,1]);
f_henon = @(x)( [  1 - a*x(1)^2 + x(2) ; b*x(1)  ] );
f_henon = @(x)(Q*f_henon(inv(Q)*x));

% Yalmip sdpvar representing the state
x = sdpvar(n,1);

% Yalmipize dynamics
f = f_henon(x);

% Constraint set X = [-1,1]^2
gx1 = 1 - x(1)^2;
gx2 = 1 - x(2)^2;


% Discount factors
ALPHA = [0.01 0.02 0.03 0.04 0.05 0.075 0.1 0.2 0.3 0.4 0.5 0.6];

% Degrees
D =  [6,8,10];

% Solver parameters
SDPsolver = lower('mosek');
sosmodel = 2; % 1 - kernel (faster), 2 - image representation (higher accuracy)
options = getSolverParams(SDPsolver,sosmodel);

%% Solve for each degree d in D and each alpha in ALPHA (could be sped up significantly by avoiding reparsing)
CVs = {};
for i = 1:numel(D)
    for j = 1:numel(ALPHA)
        fprintf('degree = %d, alpha = %f \n', D(i), ALPHA(j))        
        
        d = D(i); % Degree of the polynomials
        alpha = ALPHA(j); % Discount factor
        
        % Polynomials V1, V2, w
        [V1, cV1] = polynomial(x,d);
        [V2, cV2] = polynomial(x,d);
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
        
        % alpha * V1 \circ f <= V1 on X
        con = [ sos(V1 - alpha*replace(V1,x,f) - gx1*s1 - gx2*s2);
            sos(s1) ; sos(s2) ];
        
        % alpha * V1 \circ f <= V1 on X backward in time
        con = [ con ; sos(replace(V2,x,f) - alpha*V2 - gx1*s3 - gx2*s4 );
            sos(s3) ; sos(s4) ];
        
        % w >= V + 1 on X
        con = [ con ; sos(w - V1 - V2 - 1 - gx1*s5 - gx2*s6 ) ;
            sos(s5) ; sos(s6) ];
        
        % w >= 0 on X
        con = [ con ; sos(w - gx1*s7 - gx2*s8  ) ;
            sos(s7) ; sos(s8)];
        
        % Moments of the Lebesgue measure over X
        leb_box = getLebesgueMomentsNew(d,[-1 -1 ; 1 1],1);
        
        % Objective
        obj = cw'*leb_box;
        
        % Solve
        [sol,v,~,res] = solvesos(con,obj,options,[cw;cV1;cV2;c1;c2;c3;c4;c5;c6;c7;c8]);
        
        fprintf('\nMax SOS equality constraint residual = %f \n', norm(res,Inf));
        
        % Retrieve coefficients of the polynomials v1, v2
        cV1 = double(cV1);
        cV2 = double(cV2);
        
        % Store
        CVs{end+1} = cV1;
        CVs{end+1} = cV2;
    end
end


%% Plot
disp('Starting to plot...')

% Generate a trajectory for comparison
X = Q*[1;0];
for i = 1:400000
    X = [X f_henon(X(:,end))];
end

% Grid
gridsize = 1500; %3000 was used in the paper. Takes a bit more time.
gridx = linspace(-1,1,gridsize); gridy = linspace(-1,1,gridsize);
[X1,X2] = meshgrid(gridx,gridy);

% Take the running minimum
V_vals = polyval_grid2D(CVs,X1,X2);
plot_val = V_vals{1};
for i = 2:numel(V_vals)
    plot_val = min(plot_val,V_vals{i});
end


% Plot
[~,] = contourf(X1,X2,plot_val,[0 0],'linecolor','none'); hold on
colormap([ 0.7*[1 1 1] ]);
h = scatter(X(1,2000:end),X(2,2000:end),'MarkerFaceColor',[0.5 0 0],'MarkerEdgeColor','none'); hold on
h.SizeData = 10;
axis([-1,1,-0.5,0.5])
