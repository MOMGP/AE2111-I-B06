import math
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import sympy as smp
from scipy import integrate
from sympy import Symbol, integrate
from Aero_loading_XFLR5 import Lift_for_integrating, lift_dist_spanwise, normal_force_for_integrating, Lift_distribution_for_any_load_case
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


# INTERNAL TORQUE OF THE LIFT AROUND SHEAR CENTER CALCULATION ---------------------------------------------------------

#TORQUE ONLY FOR NORMAL AERODYNAMIC FORCE
lift, span_loc = Lift_distribution_for_any_load_case(0.7,0.31641,241.9574)

def moment_arm_normal_spanwise(x):
    chord_at_span_loc = C_r*(1-((1-taper)*(x/(b/2))))
    moment_arm_normal_torque = chord_at_span_loc/4 #assuming normal force acting at c/4 and centroid at c/2
    return moment_arm_normal_torque
'''
total_normal_torque,err_normal_torque=sp.integrate.quad(lambda x,CL_d,rho,V,n: normal_force_for_integrating(x,CL_d,rho,V,n) * moment_arm_normal_spanwise(x),0,26.78,args=(0.7,0.31641,241.9574,1),limit=50, epsabs=100)

torque_list = []
torque_error_list = []

fout=open("torque.txt", "w")
for i in np.arange(0,26.78,0.01):
    torque_result,torque_error_result=sp.integrate.quad(lambda x,CL_d,rho,V,n: normal_force_for_integrating(x,CL_d,rho,V,n) * moment_arm_normal_spanwise(x),0,i,args=(0.7,0.31641,241.9574,1),limit=50, epsabs=100)
    torque_result=(-total_normal_torque+torque_result) #Nm
    torque_list.append(torque_result)
    torque_result=str(torque_result)
#    fout.write(torque_result)
#    fout.write('\n')
    torque_error_list.append(torque_error_result)

plt.figure()
plt.plot(span_loc, torque_list, label="Torque",color='purple')
plt.xlabel("Spanwise Location [m]")
plt.ylabel("Torque [Nm]")
plt.title("Torque distribution due to normal force")
#plt.show()
'''
#TORQUE FOR ALL ACTING FORCES
torque_engine_thrust = Thrust_per_engine_perpendicular * 2.085 #based on technical drawing I am assuming that the thrust acts at center of engine which is assumed to be one radius of the engine from the center of the wingbox, which is 2.085 m. Centroid is assumed at c/2 on the camber line.
torque_engine_weight = weight_engine * ((C_r*(1-((1-taper)*((b/2 * 0.35)/(b/2)))))/2 + 3.5) #based on technical drawing I am assuming that the weight of the engine acts at 3.5 meter in front of LE, centroid is assumed at c/2

'''
total_torque = total_normal_torque + torque_engine_thrust - torque_engine_weight

torque_list = []
torque_error_list = []

for i in np.arange(0,26.78,0.001):
    torque_result,torque_error_result=sp.integrate.quad(lambda x,CL_d,rho,V,n: normal_force_for_integrating(x,CL_d,rho,V,n) * moment_arm_normal_spanwise(x),0,i,args=(0.7,0.31641,241.9574,1),limit=50, epsabs=100)
    if i >= 9.37:
        torque_result = torque_result + torque_engine_thrust - torque_engine_weight
    torque_result=(-total_torque+torque_result) #Nm
    torque_list.append(torque_result)
    torque_result=str(torque_result)
    fout.write(torque_result)
    fout.write('\n')
    torque_error_list.append(torque_error_result)

plt.figure()
plt.plot(span_loc, torque_list, label="Torque",color='purple')
plt.xlabel("Spanwise Location [m]")
plt.ylabel("Torque [Nm]")
plt.title("Total Torque distribution")
#plt.show()
'''

#In order for a better approximation of the centroid position, what is the coordinate axis system for the wingbox centroid_x and centroid_y?
#What is the exact position of the centroid of the engine?

def internal_torque_at_x(x, CL_d, rho, V, n):
    torque_list = []
    torque_error_list = []
    total_normal_torque, err_normal_torque = sp.integrate.quad(lambda x, CL_d, rho, V, n: normal_force_for_integrating(x, CL_d, rho, V, n) * moment_arm_normal_spanwise(x), 0,26.78, args=(CL_d, rho, V, n), limit=50, epsabs=100)
    total_torque = n * total_normal_torque + torque_engine_thrust - torque_engine_weight
    for i in np.arange(0,x,0.01): #x between 0 and 26.78 with two decimals
        if x > b/2:
            break
        torque_result, torque_error_result = sp.integrate.quad(lambda x, CL_d, rho, V, n: normal_force_for_integrating(x, CL_d, rho, V, n) * moment_arm_normal_spanwise(x), 0, i, args=(CL_d, rho, V, n), limit=50, epsabs=100)
        if i >= 9.37:
            torque_result = torque_result + torque_engine_thrust - torque_engine_weight
        torque_result = (-total_torque + torque_result) #Nm
        torque_list.append(torque_result)
        torque_error_list.append(torque_error_result)
    if x > b/2:
        print('invalid x entry, half-span is between x = 0 and x = 26.785 m')
        return
    print(torque_list)
    if x == 0:
        return -total_torque
    return torque_list[-1]



# -----------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------



'''
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