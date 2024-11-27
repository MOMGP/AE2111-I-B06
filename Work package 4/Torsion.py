import math
import numpy as np
import sympy as sp
from sympy import Symbol, integrate
from Aero_loading_XFLR5 import Moment_distribution_for_any_load_case
#from geometry_WP4_2 import centroid

Torsional_distribution = []
Internal_torsion = []

CL_d = 3
rho = 1.225
V = 62
Weight_engine = 9630*9.81 #N
Thrust_per_engine = 467060 #N
Lambda_LE = 0.54 #rad
C_r = 7.63 #m
b = 53.57 #m
taper = 0.3
Lambda_c2 = math.atan(math.tan(Lambda_LE) - 0.5*((2*C_r)/b)*(1-taper)) #rad
Thrust_per_engine_perpendicular = Thrust_per_engine * math.cos(Lambda_c2) #N
x_thrust = b/2 * 0.35
x_engine_weight = b/2 * 0.35

x = Symbol('x')
y_position = np.arange(0,b/2,1).tolist()

#Is the moment distribution from the aero_loading_xflr5 file the same as the lift distribution torque around the shear center?

def Torsion_integration(t):
    return integrate(t,x)

Lift_torsion = Torsion_integration(Moment_distribution_for_any_load_case(CL_d,rho,V))

for i in range(len(y_position)):
    var_sub = i
    Internal_torsion.append(Lift_torsion.subs(x, var_sub))

print(Internal_torsion)
