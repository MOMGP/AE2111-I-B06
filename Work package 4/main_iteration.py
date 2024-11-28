import numpy as np
import math
import scipy
from matplotlib import pyplot as plt
import scipy.integrate
from get_points import get_points, get_geom_from_points
from geometry_WP4_2 import centroid, moments_of_inertia
from bending import bending_moment

C_r = 7.63 #m
taper = 0.3 
b = 53.57 #m
sigma_y = 450000000 # Pa
sigma_ult = 485000000 #Pa
G = 28000000000 #Pa
E = 72400000000
rho = 2780 # kg/m3
stringer_area = None
airfoil_xyy = np.load("Airfoil_geom.npy")
possible_t = np.array([3.665, 3.264, 2.906, 2.588, 2.305, 2.053, 1.628, 1.291, 1.024, .812, .644, .511, .405, .312]) # in mm
possible_x = np.arange(0.2, 0.7, 0.05)
pos_string_num = np.arange(20, 40, 4) #based on Jeroen resource

x_y_y = np.array(get_points(0.2, 0.45, 0.65, 1))
geom = get_geom_from_points(x_y_y, [0.01 for i in range(7)])
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

def centroid_distance(centroid,airfoil_points):
    distance = []
    for coordinate in airfoil_points:
        y = math.dist(centroid, coordinate)
        distance.append(y)
    max_y = max(distance)    
    return(max_y)
# need to make tip deflection formula
def v_y(M_x,E,I_xx):
    h = M_x / (E*I_xx) 
    g = scipy.integrate.quad(h,0,b/2)
    v = scipy.integrate.quad(g,0,b/2) # tip deflection
    return(v) #idk if this works at all lol?

def bending_stress(bending_moment, y_max, I_xx):
    sigma = (bending_moment * y_max)/I_xx
    return(sigma) 

def scaled_length(length,chord):
    scaled_length = length*chord/C_r
    return(scaled_length)

#output = [mass, twist_ang, deflection, spar_pts, thicknesses, num_stringer, truncated, end_third]
#inputs = spar_pos, thicknesses, num_stringer, truncated, end_third
pos_combs = []
for i in np.arange(0.2, 0.5, 0.05):
    for j in np.arange(i, 0.65, 0.05):
        for k in np.arange(j, 0.75, 0.05):
            
