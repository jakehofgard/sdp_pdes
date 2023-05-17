x0 = [1; 1]; 
u0 = 0;
d = 10;

mpol x 2
mpol y 1
mpol z 2
m = meas([x; y; z]);

m1 = meas

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