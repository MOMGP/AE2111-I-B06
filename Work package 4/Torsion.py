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

#Creating a list of all the lift values for integration for each spanwise location
Lift_for_integrating_list = []
for i in range(0,26786):
    x = i/1000
    Lift_for_integrating_list.append(Lift_for_integrating(x,0.7,0.31641,241.9574))

#Creating a list of the product of the moment arm and the lift values for integration for each spanwise location
#torque_lift_distribution_list = []
#def torque_lift_distribution():
#    for i in range(0,26786):
#        torque_lift_distribution_list.append(Lift_for_integrating_list[i]*moment_arm_lift[i])
#    return torque_lift_distribution_list

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



#---------------------------------------------------------------------------------------------------
#x = smp.Symbol('x')
#y_position = np.arange(0,b/2,1).tolist()

#moment distribution is around c/4

#def Torsion_integration(t):
#    return smp.integrate(t,x)

#Torsional_distribution = []
#Internal_torsion = []
#Lift_dist_spanwise = []
#Span_loc = []

#Lift_dist_spanwise, Span_loc = Lift_distribution_for_any_load_case(CL_d, rho, V)

# Convert lists to numpy arrays for numerical calculations
#Lift_dist_spanwise = np.array(Lift_dist_spanwise)
#Span_loc = np.array(Span_loc)
#Chord_at_x = C_r*(1-((1-taper)*(Span_loc/(b/2))))


#y_centroid = []
#for i in range(26786):
#    y_centroid.append(0.1)
#y_centroid = np.array(y_centroid)

#Calculate Moment About the Centroid
#Moment: ∫ (y - y_centroid) * L(y) dy
#print(Span_loc)
#Lift_moment_arm = Span_loc - y_centroid  # y - y_centroid
#moment_integrand = Lift_moment_arm * Lift_dist_spanwise
#print(moment_integrand)
#Lift_torque_centroid = Torsion_integration(moment_integrand)


#print(Internal_torsion)
