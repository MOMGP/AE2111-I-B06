import math
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy.integrate import quad
import sympy as sp
from sympy import Symbol, integrate

#extra stuff
#x - y-direction along the wingspan
#y - internal shear force

#span, b
b = 53.57

shear_force_distribution = []
internal_shear_force = []


def half_span(b):
    L = b/2
    return L

#from XFLR5
#using x-dir as along the span
def integrand(x):
    return 2*x

#integrating the loading distribution
x = Symbol('x')
def shear_integration(t):
    return integrate(t, x)
#-------------********-------------*********------------*********-------------********-------------*********------------*********#

#appending to lists
y_position = np.arange(0, half_span(b), 1).tolist()

x = sp.Symbol('x')
equation = shear_integration(integrand(x))

for i in range(len(y_position)):
    var_sub = i
    internal_shear_force.append(equation.subs(x, var_sub))


#print(xpos)
plt.plot(y_position, internal_shear_force, label="Shear Force")
plt.xlabel("Y-position")
plt.ylabel("Shear Force")
plt.title("Shear Force along the span")
plt.legend()
plt.show()


