import numpy as np
import cvxpy as cp
import ncpol2sdpa
import picos
from math import sqrt
from sympy import integrate, N
import matplotlib.pyplot as plt


def J(x):
    return -2*abs(1-2*x)*sqrt(x/(1+x))


def Jk(x, coeffs):
    return sum(ci*x**i for i, ci in enumerate(coeffs))

level = 4
x = ncpol2sdpa.generate_variables('x')[0]
y = ncpol2sdpa.generate_variables('y', 2)
f = (1-2*x)*(y[0] + y[1])

gamma = [integrate(x**i, (x, 0, 1)) for i in range(1, 2*level+1)]
marginals = ncpol2sdpa.flatten([[x**i-N(gamma[i-1]), N(gamma[i-1])-x**i]
                    for i in range(1, 2*level+1)])

inequalities = [x*y[0]**2 + y[1]**2 - x,  - x*y[0]**2 - y[1]**2 + x,
                y[0]**2 + x*y[1]**2 - x,  - y[0]**2 - x*y[1]**2 + x,
                1-x, x]
sdp = ncpol2sdpa.SdpRelaxation(ncpol2sdpa.flatten([x, y]))
sdp.get_relaxation(level, objective=f, momentinequalities=marginals,
                   inequalities=inequalities)
sdp.solve()
coeffs = [sdp.extract_dual_value(0, range(len(inequalities)+1))]
coeffs += [sdp.y_mat[len(inequalities)+1+2*i][0][0] - sdp.y_mat[len(inequalities)+1+2*i+1][0][0]
           for i in range(len(marginals)//2)]

x_domain = [i/100. for i in range(100)]
plt.plot(x_domain, [J(xi) for xi in x_domain], linewidth=2.5)
plt.plot(x_domain, [Jk(xi, coeffs) for xi in x_domain], linewidth=2.5)
plt.show()
