import numpy as np
import math
import scipy
from matplotlib import pyplot as plt
import scipy.integrate
from get_points import get_points
from geometry_WP4_2 import centroid, moments_of_inertia

C_r = 7.63 #m
taper = 0.3 
b = 53.57 #m
sigma_y = 450000000 # Pa
sigma_ult = 485000000 #Pa
G = 28000000000 #Pa
E = 72400000000
rho = 2780 # kg/m3
airfoil_xyy = np.load("Airfoil_geom.npy")
possible_t = np.array([3.665, 3.264, 2.906, 2.588, 2.305, 2.053, 1.628, 1.291, 1.024, .812, .644, .511, .405, .312]) # in mm
possible_x = airfoil_xyy[:,0]

x_y_y = np.array(get_points(0.2, 0.6, 0.8, 1))
x_pts = np.array([x_y_y[:, 0], x_y_y[:, 0]]).flatten()
y_pts = np.array([x_y_y[:, 1], x_y_y[:, 2]]).flatten()
plt.scatter(x_pts, y_pts)
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
def v_y(M_x,E,I_xx):
    h = M_x / (E*I_xx)
    g = scipy.integrate.quad(h,0,b/2)
    v = scipy.integrate.quad(g,0,b/2)
    return(v) #idk if this works at all lol?


    


def bending_stress(bending_moment, y_max, I_xx):
    sigma = (bending_moment * y_max)/I_xx
    return(sigma) 

def scaled_chord(spanwise_dist):
    chord = C_r - C_r * (1 - taper) * (spanwise_dist / (b / 2))
    return(chord)

def scaled_length(length,chord):
    scaled_length = length*chord/C_r
    return(scaled_length)


plt.plot(airfoil_xyy[:,0], airfoil_xyy[:,1])
plt.plot(airfoil_xyy[:,0], airfoil_xyy[:,2])
plt.show()


# Calculate sweep at particular position