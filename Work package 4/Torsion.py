import math
import numpy as np
import scipy as sp
import sympy as smp
from scipy import integrate
from sympy import Symbol, integrate
from Aero_loading_XFLR5 import Lift_for_integrating
#from geometry_WP4_2 import centroid

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

#Creating an array list of all the spanwise locations
span_loc = []
for i in range(0,26786):
    span_loc.append(i/1000)
span_loc = np.array(span_loc)

#Creating an array list of all the chord lengths at eatch spanwise location
chord_at_span_loc = C_r*(1-((1-taper)*(span_loc/(b/2))))

#Creating a list of the moment arm of the lift with respect to the assumed shear center position
moment_arm_lift = chord_at_span_loc/4 #assuming lift at c/4 of unswept and centroid at c/2 of unswept
#moment_arm_lift = np.array(moment_arm_lift)

#Create te list of the integrated lift at each spanwise location
total_lift_list = []
for i in range(0,26786):
    x = i/1000
    total_lift,L_error=sp.integrate.quad(Lift_for_integrating,x,26.785,args=(0.7,0.31641,241.9574))
    total_lift_list.append(total_lift)

#Create the array list of the torque of the lift around the assumed shear center at each spanwise location
torque_lift_list = []
for i in range(0,26786):
    torque_lift_list.append(total_lift_list[i]*moment_arm_lift[i])
print(torque_lift_list)

#span_wise_lift = []
#for x in span_loc:
#    total_lift,L_error=sp.integrate.quad(torque_lift_distribution,x,26.785)
#    span_wise_lift.append(total_lift)
#print(span_wise_lift)




