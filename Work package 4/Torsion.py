import math
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy import integrate
from Aero_loading_XFLR5 import normal_force_for_integrating, chord_length_interpolation, Moment_for_integrating
#from geometry_WP4_2 import centroid_x, centroid_y

#Centroid coordinate system and implement into the code
#Change internal torque to negative if necessary

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
for i in range(0,2678):
    span_loc.append(i/100)
span_loc = np.array(span_loc)

#Creating an array list of all the chord lengths at eatch spanwise location
chord_at_span_loc = C_r*(1-((1-taper)*(span_loc/(b/2))))

#lift, span_loc = Lift_distribution_for_any_load_case(0.7,0.31641,241.9574)

def moment_arm_normal_spanwise(x):
    chord_at_span_loc = C_r*(1-((1-taper)*(x/(b/2))))
    moment_arm_normal_torque = chord_at_span_loc/4 #assuming normal force acting at c/4 and centroid at c/2
    return moment_arm_normal_torque


torque_engine_thrust = Thrust_per_engine_perpendicular * 2.085 #based on technical drawing I am assuming that the thrust acts at center of engine which is assumed to be one radius of the engine from the center of the wingbox, which is 2.085 m. Centroid is assumed at c/2 on the camber line.
torque_engine_weight = weight_engine * ((C_r*(1-((1-taper)*((b/2 * 0.35)/(b/2)))))/2 + 3.5) #based on technical drawing I am assuming that the weight of the engine acts at 3.5 meter in front of LE, centroid is assumed at c/2

#torque_engine_weight = weight_engine * (3.5 + centroid_x)
#torque_engine_thrust = Thrust_per_engine_perpendicular * (2.085 + centroid_y)

#def moment_arm_normal_spanwise(x):
#    chord_at_span_loc = C_r*(1-((1-taper)*(x/(b/2))))
#    moment_arm_normal_torque = centroid_x - (chord_at_span_loc/4) # assuming normal force acting at c/4 and centroid at centroid_x from LE
#    return moment_arm_normal_torque


# Function returns the internal torque due to the quarter chord pitching moment
def internal_quarter_chord_torque_spanwise(y, CL_d, rho, V):
    total_quarter_pitching_moment_coefficient, err_total_quarter_pitching_moment_coefficient = sp.integrate.quad(pitching_moment_distribution_any_CL,0,y,args=(CL_d),limit=50, epsabs=100)
    chord_at_y = chord_length_interpolation(y)
    q = 0.5 * rho * V**2
    total_quarter_chord_pitching_moment = total_quarter_pitching_moment_coefficient * q * chord_at_y
    return total_quarter_chord_pitching_moment


# Function returns the internal torque at each spanwise location x
'''
def internal_torque_at_x(x, CL_d, rho, V, n):
    torque_list = []
    torque_error_list = []
    total_normal_torque, err_normal_torque = sp.integrate.quad(lambda x, CL_d, rho, V, n: normal_force_for_integrating(x, CL_d, rho, V, n) * moment_arm_normal_spanwise(x), 0,26.78, args=(CL_d, rho, V, n), limit=50, epsabs=100)
    total_torque = total_normal_torque + internal_quarter_chord_torque_spanwise(26.78, CL_d, rho, V) + torque_engine_thrust - torque_engine_weight
    if x == 0:
        return -total_torque
    if x > b/2:
        print('invalid x entry, half-span is between x = 0 and x = 26.785 m')
        return
    for i in np.arange(0,x,0.01): #x between 0 and 26.78 with two decimals
        if x > b/2:
            break
        torque_result, torque_error_result = sp.integrate.quad(lambda x, CL_d, rho, V, n: normal_force_for_integrating(x, CL_d, rho, V, n) * moment_arm_normal_spanwise(x), 0, i, args=(CL_d, rho, V, n), limit=50, epsabs=100)
        torque_result = torque_result + internal_quarter_chord_torque_spanwise(i, CL_d, rho, V)
        if i >= 9.37:
            torque_result = torque_result + torque_engine_thrust - torque_engine_weight
        torque_result = (total_torque - torque_result) #Nm
        torque_list.append(torque_result)
        torque_error_list.append(torque_error_result)
    return torque_list[-1]
'''
def internal_torque_for_plotting(CL_d,rho,V,n):
    torque=[]
    total_torque=0
    ok=0
    print(CL_d)
    for i in np.arange(0, 26.785, 0.01):
        total_torque+=normal_force_for_integrating(i,CL_d,rho,V,n)*moment_arm_normal_spanwise(i)*0.01
        total_torque+=Moment_for_integrating(i,CL_d,rho,V)*0.01
    total_torque=total_torque+n*torque_engine_thrust-n*torque_engine_weight
    for i in np.arange(0,26.785,0.01):
        torque.append(total_torque/1000)
        total_torque-=normal_force_for_integrating(i,CL_d,rho,V,n)*moment_arm_normal_spanwise(i)*0.01
        total_torque-=Moment_for_integrating(i,CL_d,rho,V)*0.01
        if i>=9.37 and ok==0:
            total_torque=total_torque-n*torque_engine_thrust+n*torque_engine_weight
            ok=1
    return torque

'''
# Generates the internal torque diagram for the entered conditions
def internal_torque_diagram (CL_d, rho, V, n):
    total_normal_torque, err_normal_torque = sp.integrate.quad(lambda x, CL_d, rho, V, n: normal_force_for_integrating(x, CL_d, rho, V, n) * moment_arm_normal_spanwise(x), 0, 26.78, args=(CL_d, rho, V, n), limit=50, epsabs=100)
    total_torque = total_normal_torque + internal_quarter_chord_torque_spanwise(26.78, CL_d, rho, V) + torque_engine_thrust - torque_engine_weight
    torque_list = []
    torque_error_list = []
    fout = open("torque.txt", "w")
    for i in np.arange(0, 26.78, 0.01):
        torque_result, torque_error_result = sp.integrate.quad(lambda x, CL_d, rho, V, n: normal_force_for_integrating(x, CL_d, rho, V, n) * moment_arm_normal_spanwise(x), 0, i, args=(CL_d, rho, V, n), limit=50, epsabs=100)
        torque_result += internal_quarter_chord_torque_spanwise(i, CL_d, rho, V)
        if i >= 9.37:
            torque_result = torque_result + torque_engine_thrust - torque_engine_weight
        torque_result = (total_torque - torque_result)  # Nm # If sign convention must be changed to negative internal torque: torque_result = (-total_torque + torque_result)
        torque_list.append(torque_result/1000)
        torque_result = str(torque_result)
        fout.write(torque_result)
        fout.write('\n')
        torque_error_list.append(torque_error_result)
    plt.figure()
    plt.plot(span_loc, torque_list, label="Torque", color='purple')
    plt.xlabel("Spanwise Location [m]")
    plt.ylabel("Torque [kNm]")
    plt.title("Internal Torque distribution function")
    plt.show()
#    return torque_list
'''
#internal_torque_diagram(0.7,0.31641,241.9574,1)
internal_torque_diagram(0.08428454065091404,1.225,258.97,2.5) #critical load case
internal_torque_diagram(0.1896402164645566,1.225,241.96,-1) #critical load case
#internal_torque_diagram(0.08428454065091404,1.225,362.94,2.5) #critical load case
#internal_torque_diagram(0.1896402164645566,1.225,241.96,-1) #critical load case

#plt.figure()
#plt.plot(span_loc, torque_list_n_pos, label="n = 2.5", color='red')
#plt.plot(span_loc, torque_list_n_neg, label="n = -1", color='purple')
#plt.xlabel("Spanwise Location [m]")
#plt.ylabel("Torque [kNm]")
#plt.legend()
#plt.grid(True)
#plt.savefig("Torsional distribution of critical load cases", format="pdf")
#plt.show()




# -------------------------------------------------------------------------------------------------------------------------


'''
#Internal torque diagram for a specific condition

total_normal_torque, err_normal_torque = sp.integrate.quad(lambda x, CL_d, rho, V, n: normal_force_for_integrating(x, CL_d, rho, V, n) * moment_arm_normal_spanwise(x), 0,26.78, args=(0.7,0.31641,241.9574,1), limit=50, epsabs=100)
total_torque = total_normal_torque + internal_quarter_chord_torque_spanwise(26.78, 0.7,0.31641,241.9574) + torque_engine_thrust - torque_engine_weight

torque_list = []
torque_error_list = []

for i in np.arange(0,26.78,0.01):
    torque_result,torque_error_result=sp.integrate.quad(lambda x,CL_d,rho,V,n: normal_force_for_integrating(x,CL_d,rho,V,n) * moment_arm_normal_spanwise(x),0,i,args=(0.7,0.31641,241.9574,1),limit=50, epsabs=100)
    torque_result += internal_quarter_chord_torque_spanwise(i, 0.7,0.31641,241.9574)
    if i >= 9.37:
        torque_result = torque_result + torque_engine_thrust - torque_engine_weight
    torque_result=(-total_torque+torque_result) #Nm
    torque_list.append(torque_result)
    torque_result=str(torque_result)
#    fout.write(torque_result)
#    fout.write('\n')
    torque_error_list.append(torque_error_result)

plt.figure()
plt.plot(span_loc, torque_list, label="Torque",color='purple')
plt.xlabel("Spanwise Location [m]")
plt.ylabel("Torque [Nm]")
plt.title("Total Torque distribution")
plt.show()
'''
# ------------------------------------------------------------------------------------------------------------------------
'''
#TORQUE ONLY FOR NORMAL AERODYNAMIC FORCE

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
# ------------------------------------------------------------------------------------------------------------------------

CL_vals = np.load("Work package 4\\Load_case_arrays\\CL_crit.npy")
rho_vals = np.load("Work package 4\\Load_case_arrays\\Load_factor_crit.npy")
V_vals = np.load("Work package 4\\Load_case_arrays\\Rho_crit.npy")
n_vals = np.load("Work package 4\\Load_case_arrays\\V_crit.npy")

# cases= ["n_crit", "rho_crit", "V_crit", "CL_crit"]
# for i in range(len(cases)):
#     torsion = []
#     for x in np.arange(0, 26.785, 0.5):
#         torsion.append(internal_torque_at_x(x, CL_vals[i], rho_vals[i], V_vals[i], n_vals[i]))
#         print(cases[i] + " is done for length "+str(x))
#     np.save("Work package 4\\Torsions_diff_cases\\" + cases[i]+".npy", np.array(torsion))

