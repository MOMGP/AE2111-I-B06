import math
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import sympy as smp
from scipy import integrate
from sympy import Symbol, integrate
from Aero_loading_XFLR5 import Lift_for_integrating, lift_dist_spanwise
#from geometry_WP4_2 import centroid

CL_d = 3
rho = 1.225
V = 62
weight_engine = 9630*9.81 #N
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

# ACTING TORQUES AROUND SHEAR CENTER (NOT INTERNAL TORQUE) ---------------------------------------------------------
#Creating a list of the moment arm of the lift with respect to the assumed shear center position
moment_arm_lift = chord_at_span_loc/4 #assuming lift at c/4 of unswept and centroid at c/2 of unswept

#Creating a list of the lift torque around the shear center along the span
Lift_torque_spanwise = moment_arm_lift * lift_dist_spanwise

#Add contributions of the thrust of the engine
thrust_force = []
moment_arm_thrust = []
n=-1
for i in span_loc:
    n += 1
    if i >= b / 2 * 0.35:
        for s in range(0,26786):
            if s == 9375:
                thrust_force.append(Thrust_per_engine_perpendicular)
                moment_arm_thrust.append(2.085) #based on technical drawing I am assuming that the thrust acts at center of engine which is assumed to be one radius of the engine from the center of the wingbox, which is 2.085 m. Centroid is assumed at c/2 on the camber line.
            else:
                thrust_force.append(0)
                moment_arm_thrust.append(0)
        break
Thrust_torque_spanwise = np.array(thrust_force) * np.array(moment_arm_thrust)

#Add contributions of the weight of the engine
weight_force_engine = []
moment_arm_weight_engine = []
k=-1
for i in span_loc:
    k += 1
    if i >= b / 2 * 0.35:
        for s in range(0,26786):
            if s == 9375:
                weight_force_engine.append(weight_engine)
                moment_arm_weight_engine.append((chord_at_span_loc[9375]/2) + 3.5) #based on technical drawing I am assuming that the weight of the engine acts at 3.5 meter in front of LE, centroid is assumed at c/2
            else:
                weight_force_engine.append(0)
                moment_arm_weight_engine.append(0)
        break
Weight_engine_torque_spanwise = np.array(weight_force_engine) * np.array(moment_arm_weight_engine)

#List of summation of acting torques at each spanwise location
Torque_spanwise = Lift_torque_spanwise + Thrust_torque_spanwise - Weight_engine_torque_spanwise


plt.plot(span_loc, Torque_spanwise)
plt.xlabel("spanwise location")
plt.ylabel("Torque")
plt.xlim(0,b/2)
plt.gca().set_aspect(1/100000, adjustable='box')
plt.show()
# -----------------------------------------------------------------------------------------------------------------------------
'''
# INTERNAL TORQUE OF THE LIFT AROUND SHEAR CENTER CALCULATION METHOD 1---------------------------------------------------------

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
internal_torque_lift_list = []
for i in range(0,26786):
    internal_torque_lift_list.append(total_lift_list[i]*moment_arm_lift[i])

# -----------------------------------------------------------------------------------------------------------------------------

plt.subplot(1,2,1)
plt.plot(span_loc, internal_torque_lift_list)
plt.xlabel("spanwise location")
plt.ylabel("Lift Torque")
plt.xlim(0,b/2)
plt.title("Lift torque dist.")
plt.gca().set_aspect(1/100000, adjustable='box')

plt.subplot(1,2,2)
plt.plot(span_loc, total_lift_list)
plt.xlabel("spanwise location")
plt.ylabel("Lift")
plt.xlim(0,b/2)
plt.title("Lift dist.")

plt.show()
'''

#Is lift_dist_spanwise not the function I need to use to get the torque of the lift? No it is the force acting at each point but not the internal torque
#I Have a list of the forces for each spanwise location, so what do I even need to integrate? To get the internal torque there needs to be integration
#Check previous courses on the internal moment from a distributed shear force.

