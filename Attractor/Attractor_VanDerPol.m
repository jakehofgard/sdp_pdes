% Global attractor approximation for the continuous-time Van der Pol system.
% Corbinian Schlosser, Milan Korda, April 2020

clear all
close all

addpath(genpath('./functions'))

%  Parameters
n = 2;      % Dimension of the state space
d = 12;     % Degree of polynomials v1, v2 and w
beta = 0.2; % Discount factor

% Dynamics
f = @(t,x)([ 2*x(2,:) ; -0.8*x(1,:) - 10*(x(1,:).^2-0.21).*x(2,:) ]); %Van der Pol
f_ode = f;

% Yalmip sdpvar representing the state
x = sdpvar(n,1);

% Yalmipize the dynamics
f = f(0,x);

% Polynomials v1, v2, w

[w, cw] = polynomial(x,d);
[V1, cV1] = polynomial(x,d);
[V2, cV2] = polynomial(x,d);


% Constraint set X = [-1.1,1.1]^2 \ {||x|| <= 0.4}
xb = [1.1 1.1];
gx1 = xb(1)^2 - x(1)^2;
gx2 = xb(2)^2 - x(2)^2;
r = 0.4;
gx_l = x'*x - r^2;

% Lebesgue moments over X
leb_box = getLebesgueMomentsNew( d,[-xb ; xb], 1 );
leb_small_ball = momsphere(monpowers(n,d),r);
leb_X = leb_box - leb_small_ball;


% Sum of squares multipliers
[s1, c1] = polynomial(x,d-2);
[s2, c2] = polynomial(x,d-2);
[s3, c3] = polynomial(x,d-2);
[s4, c4] = polynomial(x,d-2);
[s5, c5] = polynomial(x,d-2);
[s6, c6] = polynomial(x,d-2);
[s7, c7] = polynomial(x,d-2);
[s8, c8] = polynomial(x,d-2);
[s9, c9] = polynomial(x,d-2);
[s10, c10] = polynomial(x,d-2);
[s11, c11] = polynomial(x,d-2);
[s12, c12] = polynomial(x,d-2);




%% Constraints
% V_x <= beta V on X x U (positime time)
con = [ sos(beta*V1 - jacobian(V1,x)*f - gx1*s1 - gx2*s2 - gx_l*s3 );
    sos(s1) ; sos(s2); sos(s3)  ];
% V_x <= beta V on X x U (negative time)
con = [ con ; sos(beta*V2 - jacobian(V2,x)*(-f) - gx1*s4 - gx2*s5 - gx_l*s6 );
    sos(s4) ; sos(s5) ; sos(s6) ];
% w >= V + 1 on X
con = [ con ; sos(w - V1 - V2 - 1 - gx1*s7 - gx2*s8 - gx_l*s9 ) ;
    sos(s7) ; sos(s8) ; sos(s9) ];
% w >= 0 on X
con = [ con ; sos(w - gx1*s10 - gx2*s11 - gx_l*s12 ) ;
    sos(s10) ; sos(s11) ; sos(s12)  ];

% Objective
obj = cw'*leb_X;

%% Solve
SDPsolver = lower('mosek');

sosmodel = 2; % 1 - kernel (faster), 2 - image representation (higher accuracy)
switch SDPsolver
    case 'sedumi'
        options = sdpsettings('solver','sedumi','sedumi.maxiter',100,'sedumi.eps',1e-16,'sos.model',sosmodel);
    case 'mosek'
        options = sdpsettings('solver','mosek-sdp','sos.model',sosmodel);
    otherwise
        options = sdpsettings('solver',SDPsolver,'sos.model',sosmodel);
end


% Solve
[~,~,~,res] = solvesos(con,obj,options,[cw;cV1;cV2;c1;c2;c3;c4;c5;c6;c7;c8;c9;c10;c11;c12]);

fprintf('\nMax SOS equality constraint residual = %f \n', norm(res,Inf));

% Retrieve coefficients of the polynomials v1, v2 and w
cw = double(cw);
cV1 = double(cV1);
cV2 = double(cV2);


%% Print

fprintf('\nPlotting...')

% Get a trajectory for comparison
[T,X] = ode45(f_ode,[0:0.01:1000],[0.486;-0.2608]); %solving for long time period to create the global attractor
X = X(round(size(X,1)/2):end,:);

% Grid
gridsize = 1500;
gridx = linspace(-xb(1),xb(1),gridsize); gridy = linspace(-xb(2),xb(2),gridsize);
[X1,X2] = meshgrid(gridx,gridy);

% Evaluate on the grid
v_vals = polyval_grid2D({cV1,cV2},X1,X2);

% min(v1,v2)
plot_val = min(v_vals{1},v_vals{2});

% Plot
[~,] = contourf(X1,X2,plot_val,[0 0],'linecolor','none'); hold on
colormap([ 0.7*[1 1 1] ]);
plot(X(:,1),X(:,2),'linewidth',2,'color', 	[0.6350, 0.0780, 0.1840])
axis([-1.1,1.1,-1.1,1.1])
