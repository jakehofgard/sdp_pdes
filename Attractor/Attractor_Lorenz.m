% Global attractor approximation for the continuous-time Lorenz system.
% Corbinian Schlosser, Milan Korda, April 2020

clear all
close all

addpath(genpath('./functions'))

d = 8; % degree of the polynomial
n = 3; % dimension of the state

x = sdpvar(3,1);
[w, cw] = polynomial(x,d);

[V1, cV1] = polynomial(x,d);
[V2, cV2] = polynomial(x,d);

% Dynamics
f = @(t,x)([10*(x(2,:)-x(1,:)) ; x(1,:).*(28-x(3,:)) - x(2,:) ; x(1,:).*x(2,:) - (8/3)*x(3,:) ] );
Q = diag([1/25,1/30,1/50]); % Scale to unit box
f = @(t,x)(Q*f(0,inv(Q)*x));
f_ode = f;


% Yalmipize the dynamics
f = f(0,x);

% Constraint box
xb = 1*[1 1 1];

gx1 = xb(1)^2 - x(1)^2;
gx2 = xb(2)^2 - x(2)^2;
gx3 = xb(3)^2 - x(3)^2;

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


% Moments of Lebesgue measure
leb_box = getLebesgueMomentsNew( d,[-xb ; xb], 1 );


% Constraints
beta = 1; % Discount factor

% nabla V1 . f <= beta V on X x U (positive time)
con = [ sos(beta*V1 - jacobian(V1,x)*f - gx1*s1 - gx2*s2 - gx3*s3);
    sos(s1) ; sos(s2) ; sos(s3)  ];
% nabla V2 . (-f) <= beta V on X x U (negative time)
con = [ con ; sos(beta*V2 - jacobian(V2,x)*(-f) - gx1*s10 - gx2*s11 - gx3*s12);
    sos(s10) ; sos(s11) ; sos(s12)  ];
% w >= V1 + V2 + 1 on X
con = [ con ; sos(w - V1 - V2 - 1 - gx1*s4 - gx2*s5 - gx3*s6) ;
    sos(s4) ; sos(s5) ; sos(s6)];
% w >= 0 on X
 con = [ con ; sos(w - gx1*s7 - gx2*s8 - gx3*s9 ) ;
         sos(s7) ; sos(s8) ; sos(s9) ];


% Objective
obj = cw'*leb_box;

% SDP solver selection and parameters
SDPsolver = lower('mosek');

sosmodel = 1; % 1 - kernel (faster), 2 - image representation (higher accuracy)
switch SDPsolver
    case 'sedumi'
        options = sdpsettings('solver','sedumi','sedumi.maxiter',100,'sedumi.eps',1e-12,'sos.model',sosmodel);
    case 'mosek'
        options = sdpsettings('solver','mosek-sdp','sos.model',sosmodel);
    otherwise
         options = sdpsettings('sos.model',sosmodel);
end

%% Solve
[~,~,~,res] = solvesos(con,obj,options,[cw;cV1;cV2;c1;c2;c3;c4;c5;c6;c7;c8;c9;c10;c11;c12]);

 % Retrieve coefficients of the polynomials v1, v2 and w
cw = double(cw);
cV1 = double(cV1);
cV2 = double(cV2);

% Print info
fprintf('\nDegree = %2d, Solver = %s\n',d,SDPsolver);
fprintf('Max SOS equality constraint residual = %d  \n', norm(res,Inf));


%% Plot

disp('Starting to plot. This may take some time...')

% Simulate a trajectory for comparison
x0 = [0.5;0.5;0.6];
[T,X] = ode45(f_ode,[0,10000],x0); X = X(round(size(X,1) / 5):end,:);

grid_size = 300; % Used 600 for the plot in the paper - takes a longer time.
gridx = linspace(-1,1,grid_size); gridy = linspace(-1,1,grid_size); gridz = linspace(0,1,grid_size);
[X1,X2,X3] = meshgrid(gridx,gridy,gridz);

% Evaluate v1 and v2 on the grid
[V1_val,V2_val] = eval_v1v2_grid3D(cV1,cV2,d,X1,X2,X3);

% min(v1,v2)
plot_val = min(V1_val,V2_val);

% Plot (The artifacts are due to the plotting grid being too coarse)
h = patch(isosurface(X1,X2,X3,plot_val,0),'FaceColor',[1 0 0 ],'EdgeColor','none'); hold on
isonormals(X1,X2,X3,plot_val,h);
view(3); camlight; lighting phong;
axis equal tight vis3d;
alpha(0.3)
camlight left
view(-44,8)
view(-21.5522,   61.5357)
box off; axis off
 
hold on
plot3(X(:,1),X(:,2),X(:,3),'-k')



