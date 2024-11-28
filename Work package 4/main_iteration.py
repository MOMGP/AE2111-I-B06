import numpy as np
import math
import scipy
from matplotlib import pyplot as plt
from scipy.integrate import dblquad, quad
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
plt.plot(x_vals, y_vals)
plt.plot(airfoil_xyy[:,0], airfoil_xyy[:,1])
plt.plot(airfoil_xyy[:,0], airfoil_xyy[:,2])
plt.show()

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

def I_xx(x): # dummy ixx :3
    return 0.1 - 0.002 * x 

def load_distribution(x):
    return 1000 * x**2

def shear_force(x):
    result, _ = quad(load_distribution, 0, x)  # Shear force is the integral of load
    return result
def bending_moment(x): #dummy bending moment formula
    return 800000 * x 
def angle_of_rotation(x):
    result, _ = quad(lambda xi: bending_moment(xi) / (E * I_xx(x)), 0, x)
    return result
def deflection(x):
    result, _ = quad(angle_of_rotation,0,x)
    return result

# Calculate and plot results
x_vals = np.linspace(0, b/2, 100) #create # amount of values evenly spaced
deflection_vals = [deflection(x) for x in x_vals] 

plt.figure(figsize=(8, 4))
# Deflection
#plt.subplot(1, 1, 1) bending moment can be added later aswell if needed :3
plt.plot(x_vals, deflection_vals, label="Deflection w(x)", color="purple")
plt.xlabel("Position along the beam (m)")
plt.ylabel("Deflection (m)")
plt.grid(alpha=0.3)
plt.legend()
plt.show()

print('Angle of rotation is: ', angle_of_rotation(b/2 * 180 / math.pi), 'degrees')
def bending_stress(bending_moment, y_max, I_xx):
    sigma = (bending_moment * y_max)/I_xx
    return(sigma) 

def scaled_length(length,chord):
    scaled_length = length*chord/C_r
    return(scaled_length)
