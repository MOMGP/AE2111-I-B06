import numpy as np
from matplotlib import pyplot as plt
from get_points import get_points
from geometry_WP4_2 import centroid, moments_of_inertia


sigma_y = 450000000 # Pa
sigma_ult = 485000000 #Pa
G = 28000000000 #Pa
E = 72400000000 #Pa
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
    # Bending
    

    # Torsion


    return 0

plt.plot(airfoil_xyy[:,0], airfoil_xyy[:,1])
plt.plot(airfoil_xyy[:,0], airfoil_xyy[:,2])
plt.show()

# ADD WING TIP DEFLECTION AND TWIST 
# Calculate sweep at particular position