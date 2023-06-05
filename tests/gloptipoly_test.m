% Test method of moments on one-dimensional method-of-moments problem,
% using double-integrator example from Lasserre.

x0 = [1; 1]; 
u0 = 0;
d = 20;

if x0(1) >= -(x0(2)^2-2)/2
    tmin = 1+x0(1)+x0(2)+x0(2)^2/2;
elseif x0(1) >= -x0(2)^2/2*sign(x0(2))
    tmin = 2*sqrt(x0(1)+x0(2)^2/2)+x0(2);
else
    tmin = 2*sqrt(-x0(1)+x0(2)^2/2)-x0(2);
end

mpol x1 2
mpol u1
m1 = meas([x1; u1]);

mpol x2 2
m2 = meas(x2);

scaling = tmin;
f = scaling * [x1(2); u1];

g1 = mmon(x1, d);
g2 = mmon(x2, d);

assign([x1; u1],[x0; u0]);
g0 = double(g1);

P = msdp(min(mass(m1)),...
    u1^2 <= 1,...
    x1(2) >= -1,...
    x2'*x2 <= 0,...
    mom(g2) - g0 == mom(diff(g1, x1) * f)); 

[status, obj] = msol(P);
obj = scaling * obj;

disp(["Minimum time = " num2str(tmin)]);
disp(["LMI " int2str(d) " lower bound = " num2str(obj)])