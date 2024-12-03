import numpy as np
import math
import time
import scipy
from matplotlib import pyplot as plt
from scipy.integrate import dblquad, quad, quad_vec
from get_points import get_points, get_geom_from_points
from geometry_WP4_2 import centroid, moments_of_inertia,scaled_chord,scaled_length,centroid_distance
from Aero_loading_XFLR5 import normal_force_for_integrating
# from bending import bending_moment, S

C_r = 7.63 #m
taper = 0.3 
b = 53.57 #m
sigma_y = 450000000 # Pa
sigma_ult = 485000000 #Pa
G = 28000000000 #Pa
E = 72400000000
rho = 2780 # kg/m3
I_xx_r = 1 # dummy value 
stringer_area = None
airfoil_xyy = np.load("Airfoil_geom.npy")
possible_t = np.array([3.665, 3.264, 2.906, 2.588, 2.305, 2.053, 1.628, 1.291, 1.024, .812, .644, .511, .405, .312]) # in mm
possible_x = np.arange(0.2, 0.7, 0.05)
pos_string_num = np.arange(20, 40, 4) #based on Jeroen resource
graphs = False
x_y_y = np.array(get_points(0.2, 0.45, 0.65, 1))
geom = get_geom_from_points(x_y_y, [0.01 for i in range(7)])
x_vals = []
y_vals = []
for i in geom:
    x_vals.append(i[0][0])
    x_vals.append(i[1][0])
    y_vals.append(i[0][1])
    y_vals.append(i[1][1])
# plt.plot(x_vals, y_vals)
# plt.plot(airfoil_xyy[:,0], airfoil_xyy[:,1])
# plt.plot(airfoil_xyy[:,0], airfoil_xyy[:,2])
# plt.show()


def bending_stress(bending_moment, y_max, I_xx):
    sigma = (bending_moment * y_max)/I_xx
    return(sigma)     

def I_xx(I_xx_root,x): #Assuming I_xx scales linearly with length, as thickness presumably stays the same
    result = scaled_length(I_xx_root,scaled_chord(x))
    return result

shear = []
bending = []
n_points = 1 # decrease number for higher amount of points to integrate over (a.k.a. increase accuracy
with open("shear.txt", 'r') as f:
    count = 0
    for i in f:
        line = i.strip('\n')
        count += 1
        if count % n_points == 0 or count == 1:
            shear.append(line)
with open("bending copy.txt", 'r') as f:
    count = 0
    for i in f:
        line = i.strip('\n')
        count += 1
        if count % n_points == 0 or count == 1:
            bending.append(line)          
x_vals = np.linspace(0,b/2,len(shear)) # This doesn't give completely accurate numbers but its slighty off since i didnt find a better method
bending_x_vals = np.arange(0,b/2,1) # taken from bending.py (assuming since only 27 values are given in bending.txt)
shear_interp_func = scipy.interpolate.interp1d(x_vals, shear, kind='quadratic', fill_value="extrapolate") #Unsure if quadratic is appropriate
# bending_interp_func = scipy.interpolate.interp1d(bending_x_vals, bending, kind='quadratic', fill_value="extrapolate")

start = time.time()
def shear_force(x): 
    result = shear_interp_func(x)
    # result, _ = quad(load_distribution, 0, x, limit=200)  # Shear force is the integral of load
    return result
def bending_moments(x): 
    result, _ = quad(shear_force,x,b/2)
    return result # + bending_interp_func #when i get it!!
aof_0, _ = quad(lambda xi: bending_moments(xi) / (E * I_xx(I_xx_r,xi)), 0,b/2, limit=50)
def angle_of_rotation(x): # 
    result, _ = quad(lambda xi: bending_moments(xi) / (E * I_xx(I_xx_r,xi)), x,b/2, limit=50)
    return result - aof_0
def_0, _ = quad(angle_of_rotation,0,b/2,limit=50)
def deflection(x):
    result, _ = quad(angle_of_rotation,x,b/2)
    return result - def_0

print('Maximum deflection is:', abs(deflection(b/2)), "; which is", abs(deflection(b/2)/(b/2) * 100), "% of the wingspan")
end = time.time()
print("integration took", end - start, "seconds")
print('Integrated over', len(x_vals), 'intervals of x')
print(bending)


# #output = [mass, twist_ang, deflection, spar_pts, thicknesses, num_stringer, truncated, end_third]
# #inputs = spar_pos, thicknesses, num_stringer, truncated, end_third
# pos_combs = []
# # for i in np.arange(0.2, 0.5, 0.05):
# #     for j in np.arange(i, 0.65, 0.05):
# #         for k in np.arange(j, 0.75, 0.05):



# Following lines are not necessary for iterating the design, this is purely for analytical purposes
if graphs == True:
    x_vals = np.linspace(0,b/2,10)
    shear_vals = [shear_force(x) for x in x_vals]
    moment_vals = [bending_moments(x) for x in x_vals]
    angle_of_rotation_vals = [angle_of_rotation(x) for x in x_vals]
    deflection_vals = [deflection(x) for x in x_vals] 
    I_xx = [I_xx(I_xx_r,x) for x in x_vals]
    end = time.time()
    print(end - start, "seconds")


    plt.figure(figsize=(12, 6))

    # I_xx vals
    plt.subplot(5, 1, 5)
    plt.plot(x_vals, I_xx, label="I_xx", color="yellow")
    plt.ylabel("I_xx")
    plt.grid(alpha=0.3)
    plt.legend()

    # Shear force
    plt.subplot(5, 1, 1)
    plt.plot(x_vals, shear_vals, label="Shear Force V(x)", color="green")
    plt.ylabel("Shear Force (N)")
    plt.grid(alpha=0.3)
    plt.legend()

    # Bending moment
    plt.subplot(5, 1, 2)
    plt.plot(x_vals, moment_vals, label="Bending Moment M(x)", color="red")
    plt.ylabel("Moment (Nm)")
    plt.grid(alpha=0.3)
    plt.legend()

    plt.subplot(5, 1, 3) # bending moment can be added later aswell if needed :3
    plt.plot(x_vals, angle_of_rotation_vals, label="Slope dv/dx", color="purple")
    plt.xlabel("Position along the beam (m)")
    plt.ylabel("Slope (m)")
    # plt.grid(alpha=0.3)
    plt.legend()



    # Deflection
    plt.subplot(5, 1, 4) #bending moment can be added later aswell if needed :3
    plt.plot(x_vals, deflection_vals, label="Deflection w(x)", color="purple")
    plt.xlabel("Position along the beam (m)")
    plt.ylabel("Deflection (m)")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.show()

    
