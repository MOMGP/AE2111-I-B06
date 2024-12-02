import numpy as np
import math
import time
import scipy
from matplotlib import pyplot as plt
from scipy.integrate import dblquad, quad, quad_vec
from get_points import get_points, get_geom_from_points
from geometry_WP4_2 import centroid, moments_of_inertia
# from bending import bending_moment, S

C_r = 7.63 #m
taper = 0.3 
b = 53.57 #m
sigma_y = 450000000 # Pa
sigma_ult = 485000000 #Pa
G = 28000000000 #Pa
E = 72400000000
rho = 2780 # kg/m3
#I_xx = 1000 # dummy Ixx value
#I_xx, I_yy = moments_of_inertia(wing_box,stringers)
airfoil_xyy = np.load("Airfoil_geom.npy")
possible_t = np.array([3.665, 3.264, 2.906, 2.588, 2.305, 2.053, 1.628, 1.291, 1.024, .812, .644, .511, .405, .312]) # in mm
possible_x = airfoil_xyy[:,0]

x_y_y = np.array(get_points(0.2, 0.45, 0.65, 1))
geom = get_geom_from_points(x_y_y)
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

def check_reqs():
    return 0

def centroid_distance(centroid,airfoil_points):
    distance = []
    for coordinate in airfoil_points:
        y = math.dist(centroid, coordinate)
        distance.append(y)
    max_y = max(distance)    
    return(max_y)
# need to make tip deflection formula
# def v_y(M_x,E,I_xx):
#     # M_x = bending_moment(S)
#     h = M_x / (E*I_xx) 
#     g = scipy.integrate.quad(h,0,b/2)
#     v = scipy.integrate.quad(g,0,b/2) # tip deflection
#     return(v) #idk if this works at all lol?
# print(v_y(bending_moment(S),E,I_xx))



# Integrate the bending moment to find the deflection w(x)
start = time.time()
def I_xx(x): # dummy ixx
    return 0.5 - 0.015 * x 
def load_distribution(x): #dummy load distribution
    if x <= 10:
        return 50000
    else:
        return 40000 - 10 * x ** 2
def shear_force(x): 
    result, _ = quad(load_distribution, 0, x, limit=500)  # Shear force is the integral of load
    return result
def bending_moment(x): 
    result, _ = quad(shear_force, 0, x, limit=500)  # Bending moment is the integral of shear force
    return result 
def angle_of_rotation(x): # 
    result, _ = quad(lambda xi: bending_moment(xi) / (E * I_xx(x)), 0, x, limit=500)
    return result
# def deflection(x):
    result, _ = quad(angle_of_rotation,0,x)
    return result

# Calculate and plot results
x_vals = np.linspace(0, b/2, 10) #create # amount of values evenly spaced
load_vals = [load_distribution(x) for x in x_vals]
shear_vals = [shear_force(x) for x in x_vals]
moment_vals = [bending_moment(x) for x in x_vals]
angle_of_rotation_vals = [angle_of_rotation(x) for x in x_vals]
# deflection_vals = [deflection(x) for x in x_vals] 


print('Angle of rotation is: ', angle_of_rotation(b/2), 'degrees or rad idk')
end = time.time()
print(end - start, "seconds")


plt.figure(figsize=(8, 4))

# Load distribution
plt.subplot(5, 1, 1)
plt.plot(x_vals, load_vals, label="Load Distribution q(x)", color="blue")
plt.ylabel("Load (N/m)")
plt.grid(alpha=0.3)
plt.legend()

# Shear force
plt.subplot(5, 1, 2)
plt.plot(x_vals, shear_vals, label="Shear Force V(x)", color="green")
plt.ylabel("Shear Force (N)")
plt.grid(alpha=0.3)
plt.legend()

# Bending moment
plt.subplot(5, 1, 3)
plt.plot(x_vals, moment_vals, label="Bending Moment M(x)", color="red")
plt.ylabel("Moment (Nm)")
plt.grid(alpha=0.3)
plt.legend()

plt.subplot(5, 1, 4) # bending moment can be added later aswell if needed :3
plt.plot(x_vals, angle_of_rotation_vals, label="Deflection w(x)", color="purple")
plt.xlabel("Position along the beam (m)")
plt.ylabel("Deflection (m)")
plt.grid(alpha=0.3)
plt.legend()
plt.show()



# Deflection
#plt.subplot(1, 1, 1) bending moment can be added later aswell if needed :3
# plt.plot(x_vals, deflection_vals, label="Deflection w(x)", color="purple")
# plt.xlabel("Position along the beam (m)")
# plt.ylabel("Deflection (m)")
# plt.grid(alpha=0.3)
# plt.legend()
# plt.show()


def bending_stress(bending_moment, y_max, I_xx):
    sigma = (bending_moment * y_max)/I_xx
    return(sigma) 

def scaled_length(length,chord):
    scaled_length = length*chord/C_r
    return(scaled_length)
