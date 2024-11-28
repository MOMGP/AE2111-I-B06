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

span_loc = []
for i in range(0,26786):
    span_loc.append(i/1000)

span_wise_lift = []
for x in span_loc:
    total_lift,L_error=sp.integrate.quad(Lift_for_integrating,x,26.785,args=(0.7,0.31641,241.9574))
    span_wise_lift.append(total_lift)
#If lift is c/4 on unswept wing, then moment arm is equal to distance between c/4 and the centroid, which I assume to be c/2?
print(span_wise_lift)


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
