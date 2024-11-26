import math
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy.integrate import quad
import sympy as sp
from sympy import Symbol, integrate
from Aero_loading_XFLR5 import Lift_distribution_for_any_load_case
#extra stuff
#x - y-direction along the wingspan
#y - internal shear force

CL_d = 3
rho = 1.225
V = 62
#span, b

shear_force_distribution = []
internal_shear_force = []


def half_span(b):
    L = b/2
    return L

#from XFLR5
#using x-dir as along the span
#def integrand(x):
 #   return Lift_distribution_for_any_load_case

#integrating the loading distribution
x = Symbol('x')
def shear_integration(t):
    return integrate(t, x)
#-------------********-------------*********------------*********-------------********-------------*********------------*********#

#appending to lists
y_position = np.arange(0, half_span(b), 1).tolist()

x = sp.Symbol('x')
shearequation = shear_integration(Lift_distribution_for_any_load_case(CL_d,rho,V))

for i in range(len(y_position)):
    var_sub = i
    internal_shear_force.append(shearequation.subs(x, var_sub))


#print(xpos)
plt.plot(y_position, internal_shear_force, label="Shear Force")
plt.xlabel("Y-position")
plt.ylabel("Shear Force")
plt.title("Shear Force along the span")
plt.legend()
plt.show()


