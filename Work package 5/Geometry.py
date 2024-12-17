import numpy as np


C_r = 7.63 #m
taper = 0.3 
b = 53.57 #m
dihedral = np.deg2rad(4.75) #deg
sigma_y = 450000000 # Pa
sigma_ult = 485000000 #Pa
G = 28000000000 #Pa
E = 72400000000
rho = 2780 # kg/m3
AR = 10.82
M_CR = 0.82
stringer_area = 2*10**(-4) #m^3



geometry = []
b = 53.57



def scaled_chord(spanwise_dist):
    chord = C_r - C_r * (1 - taper) * (spanwise_dist / (b / 2))
    return(chord)


def get_design(design_number):
    if design_number == 0:          #first design
        spar1 = 0.2
        spar2 = 0.35
        spar3 = 0.65
        thickness_sides = 0.008
        thickness_top_bottom = 0.014
    elif design_number == 1:        #second design
        spar1 = 0.2
        spar2 = 0.4
        spar3 = 0.7
        thickness_sides = 0.01
        thickness_top_bottom = 0.016
    elif design_number == 2:        #third design
        spar1= 0.2
        spar2 = 0.35
        spar3 = 0.7
        thickness_sides = 0.01
        thickness_top_bottom = 0.022
    geometry.append((spar1, spar2, spar3))
    geometry.append((thickness_sides, thickness_top_bottom))
    end_second_cell = 0.175*b
    geometry.append(end_second_cell)
    return geometry

path_airfoil = "Airfoil_geom.npy"

def centroid_z(type, x_over_c1, x_over_c2, y_end, thickness, top=False):
    centroid_z=0
    if type == "spar":
        centroid_z=0
    elif type == "flange" or type =="skin":
        airfoil_geometry = np.load("Airfoil_geom.npy")
        starting_index = -1
        ending_index = -1
        for i in range(len(airfoil_geometry)-1):  #input values for length 1
            j = airfoil_geometry[i][0]
            k = airfoil_geometry[i+1][0]
            if j<x_over_c1<k:
                starting_index = i
            if j<x_over_c2<k:
                ending_index = i
                break
        if top:
            angle = np.arctan((airfoil_geometry[ending_index][1]-airfoil_geometry[starting_index][1])/(airfoil_geometry[ending_index][0]-airfoil_geometry[starting_index][0]))
            for i in range(starting_index, ending_index):
                x_mid = (airfoil_geometry[i][0]+airfoil_geometry[i+1][0])*scaled_chord(y_end)
                y_mid = (airfoil_geometry[i][1]+airfoil_geometry[i+1][1])*scaled_chord(y_end)
                centroid_z+=np.sqrt(((airfoil_geometry[i][0]-airfoil_geometry[i+1][0])*scaled_chord(y_end))**2+((airfoil_geometry[i][1]-airfoil_geometry[i+1][1])*scaled_chord(y_end))**2)*thickness*np.sqrt((x_mid*np.cos(angle)-y_mid*np.sin(angle)-x_mid*np.cos(angle))**2+(y_mid*np.cos(angle)+x_mid*np.sin(angle))**2)
                #todo add stringers
            
        else:
            angle = np.arctan((airfoil_geometry[ending_index][2]-airfoil_geometry[starting_index][2])/(airfoil_geometry[ending_index][0]-airfoil_geometry[starting_index][0]))
            for i in range(starting_index, ending_index):
                x_mid = (airfoil_geometry[i][0]+airfoil_geometry[i+1][0])*scaled_chord(y_end)
                y_mid = (airfoil_geometry[i][2]+airfoil_geometry[i+1][2])*scaled_chord(y_end)
                centroid_z+=np.sqrt(((airfoil_geometry[i][0]-airfoil_geometry[i+1][0])*scaled_chord(y_end))**2+((airfoil_geometry[i][2]-airfoil_geometry[i+1][2])*scaled_chord(y_end))**2)*thickness*np.sqrt((x_mid*np.cos(angle)-y_mid*np.sin(angle)-x_mid*np.cos(angle))**2+(y_mid*np.cos(angle)+x_mid*np.sin(angle))**2)
                #todo add stringers
    else:
        print("centroid_z no valid type")
    




    return centroid_z

def I_xx(type, x_over_c1, x_over_c2, y_begin, y_end, Top=False): #Type can be 'spar', 'flange', or 'skin'
    
    if type=="spar":
        1

    elif type == "flange":
        1

    elif type == "skin":
        1
    else:
        print("centroid_z no valid type")

        


# Main Design Philosophy Minimum weight Buckling Ultimate load
# Thickness top and bottom [m] 0.014 0.016 0.022
# Thickness sides [m] 0.008 0.01 0.01
# Number of stringers [-] 32 40 20
# Mass [kg] 15568 20345 22631
# Critical tip deflection [m] 8.01 6.78 7.98 1
# Critical tip twist [deg] -9.93 -9.74 -9.69 2
# Front spar [x1/c] 0.2 0.2 0.2
# Middle spar [x2/c] 0.35 0.4 0.35
# End spar [x3/c] 0.65 0.7 0.7
# End of secondary cell 0.175b 0.175b 0.175b
# Truncated Yes Yes Yes