import numpy as np
from get_points import get_points
# format for wing_box = [(start), (end), thickness]
#format of stringers = [pos, area]
"""
Geometry:
 ______________2_____________ ______________5_____________ 
|                            |                            |
|                            |                            |
|                            |                            |
1                            3                            7
|                            |                            |
|                            |                            |
|_____________4______________|_____________6______________| 


"""
C_r = 7.63 #m
taper = 0.3 
b = 53.57 #m
dihedral = 4.75 #deg
sigma_y = 450000000 # Pa
sigma_ult = 485000000 #Pa
G = 28000000000 #Pa
E = 72400000000
rho = 2780 # kg/m3
AR = 10.82
M_CR = 0.82

def scaled_chord(spanwise_dist):
    chord = C_r - C_r * (1 - taper) * (spanwise_dist / (b / 2))
    return(chord)

def centroid(wing_box, stringers):
    centroid_sum_x =0
    centroid_sum_y = 0
    area_sum = 0
    for i in range(len(wing_box)):
        start= wing_box[i,0]
        end= wing_box[i,1]
        length = np.sqrt((start[0]-end[0])**2+(start[1]-end[1])**2)
        thickness = wing_box[i,2]
        centroid_sum_x+=(start[0]+end[0])/2*length*thickness
        centroid_sum_y+=(start[1]+end[1])/2*length*thickness
        area_sum += length*thickness
    
    for i in range(stringers):
        centroid_sum_x+= stringers[i,0,0]*stringers[i,1]
        centroid_sum_y+= stringers[i,0,1]*stringers[i,1]
        area_sum+=stringers[i, 1]
    
    centroid_x = centroid_sum_x/area_sum
    centroid_y = centroid_sum_y/area_sum
    return centroid_x, centroid_y
        
def moments_of_inertia(wing_box, stringers):
    centroid_x, centroid_y = centroid(wing_box, stringers)
    I_xx = 0
    I_yy = 0
    I_xy = 0
    for i in range(len(wing_box)):
        start= wing_box[i,0]
        end= wing_box[i,1]
        midpoint = np.array([(start[0]+end[0])/2, (start[1]+end[1])/2])
        length = np.sqrt((start[0]-end[0])**2+(start[1]-end[1])**2)
        angle = np.arctan((start[1]-end[1])/(start[0]-end[0]))
        thickness = wing_box[i,2]
        I_xx+=thickness*length**3*np.sin(angle)**2/12 #reg
        I_xx+=thickness*length*(midpoint[1]-centroid_y)**2
        I_yy+=thickness*length**3*np.cos(angle)**2/12 #reg
        I_yy+=thickness*length*(midpoint[0]-centroid_x)**2
        #PAT for I_xy?
        

    return I_xx, I_yy

def Lambda_n(n_percent):
    Lambda_c4 = np.rad2deg(np.arccos(1.16 / (M_CR + 0.5)))
    taper = 0.2 * (2 - np.deg2rad(Lambda_c4))
    return np.rad2deg(np.arctan(np.tan(np.deg2rad(Lambda_c4))-(4/AR)*((n_percent-25)/100*(1-taper)/(1+taper))))

def get_mass(wing_box_root, wing_box_mid, wing_box_tip, pos_mid, stringers_root, stringers_tip):
    mass = 0
    for i in range(len(wing_box_root)):
        start_root = wing_box_root[i, 0]
        end_root = wing_box_root[i, 1]
        length_root = np.sqrt((start_root[0] - end_root[0]) ** 2 + (start_root[1] - end_root[1]) ** 2)
        thickness_root = wing_box_root[i, 2]
        area_root = length_root*thickness_root

        start_mid = wing_box_mid[i, 0]
        end_mid = wing_box_mid[i, 1]
        length_mid = np.sqrt((start_mid[0] - end_mid[0]) ** 2 + (start_mid[1] - end_mid[1]) ** 2)
        thickness_mid = wing_box_mid[i, 2]
        area_mid = length_mid * thickness_mid

        sweep = np.deg2rad((Lambda_n(start_root / scaled_chord(0))+Lambda_n(end_root / scaled_chord(0))) /2)
        y = b * pos_mid / np.tan(sweep) / np.tan(dihedral)
        mass += (area_root+area_mid)/2 * y * rho



    for i in range(len(wing_box_tip)):
        start_tip = wing_box_tip[i, 0]
        end_tip = wing_box_tip[i, 1]
        length_tip = np.sqrt((start_tip[0] - end_tip[0]) ** 2 + (start_tip[1] - end_tip[1]) ** 2)
        thickness_tip = wing_box_tip[i, 2]
        area_tip = length_tip * thickness_tip

        start_mid = wing_box_mid[i, 0]
        end_mid = wing_box_mid[i, 1]
        length_mid = np.sqrt((start_mid[0] - end_mid[0]) ** 2 + (start_mid[1] - end_mid[1]) ** 2)
        thickness_mid = wing_box_mid[i, 2]
        area_mid = length_mid * thickness_mid

        sweep = np.deg2rad((Lambda_n(start_tip / scaled_chord(b))+Lambda_n(end_tip / scaled_chord(b))) / 2)
        y = b * (1-pos_mid) / np.tan(sweep) / np.tan(dihedral)
        mass += (area_tip + area_mid) / 2 * y * rho

    for i in range(len(stringers_root)):

        start_root = stringers_root[i, 0]
        area_root = stringers_root[1]

        sweep = np.deg2rad((Lambda_n(start_root / scaled_chord(0))) / 2)
        y = b * (1 - pos_mid) / np.tan(sweep) / np.tan(dihedral)
        mass += area_root * y * rho

    for i in range(len(stringers_tip)):

        start_tip = stringers_tip[i, 0]
        area_tip = stringers_tip[1]

        sweep = np.deg2rad((Lambda_n(start_tip / scaled_chord(b))) / 2)
        y = b * (1 - pos_mid) / np.tan(sweep) / np.tan(dihedral)
        mass += area_tip  * y * rho

    return mass

# format for wing_box = [(start), (end), thickness]
#format of stringers = [pos, area]
def get_points_along_spanwise(norm_wing_box_root, norm_stringers, y, end_third_spar, trunctated = False):
    chord = scaled_chord(y)
    geometry = []
    stringers = []
    if trunctated:
        if y<end_third_spar:
            for i in norm_wing_box_root:
                plate = np.array([i[0]*chord,i[1]*chord, i[2]])
                geometry.append(plate)
        else:
            for i in range(4):
                plate = np.array([norm_wing_box_root[i,0]*chord, norm_wing_box_root[i,1]*chord, norm_wing_box_root[i,2]])
                geometry.append(plate)
    else:
        if y<end_third_spar:
            for i in range(4):
                plate = np.array([norm_wing_box_root[i,0]*chord, norm_wing_box_root[i,1]*chord, norm_wing_box_root[i,2]])
                geometry.append(plate)
            x_pos_third_spar = (norm_wing_box_root[6,1,0]*(1-y) + norm_wing_box_root[3, 0, 0]*y)/y
            x_y_y = np.array(get_points(-1, -1, x_pos_third_spar, 1))*chord
            geometry.append(np.array([[norm_wing_box_root[2,0]*chord, [x_y_y[0], x_y_y[1]]], norm_wing_box_root[4,2]]))
            geometry.append(np.array([[norm_wing_box_root[3,0]*chord, [x_y_y[0], x_y_y[2]]], norm_wing_box_root[5,2]]))
            geometry.append(np.array([[norm_wing_box_root[3,0]*chord, [x_y_y[0], x_y_y[2]]], norm_wing_box_root[6,2]]))


        else:
            for i in range(4):
                plate = np.array([norm_wing_box_root[i,0]*chord, norm_wing_box_root[i,1]*chord, norm_wing_box_root[i,2]])
                geometry.append(plate)

    for i in norm_stringers:
        stringers.append(i[0]*chord, i[1])

    stringers = np.array(stringers)
    geometry = np.array(geometry)
    return geometry


profile_1 = np.array([[(0,0),(0,1),0.01],[(0,1),(1,1),0.01],[(1,1),(1,0),0.01],[(1,0),(0,0),0.01],[(1,1),(2,1),0.01],[(2,1),(2,0),0.01],[(2,0),(1,0),0.01]])
profile_2 = np.array([[(0,0),(0,1),0.01],[(0,1),(1,1),0.01],[(1,1),(1,0),0.01],[(1,0),(0,0),0.01],[(1,1),(2,1),0.01],[(2,1),(2,0),0.01],[(2,0),(1,0),0.01]])
profile_3 = np.array([[(0,0),(0,1),0.01],[(0,1),(1,1),0.01],[(1,1),(1,0),0.01],[(1,0),(0,0),0.01]])
pos_profile_3 = 0.5
stringers_1 = np.array([[(0.5,0.05),0.5],[(1.5,0.05),0.5],[(0.5,0.95),0.5],[(1.5,0.95),0.5]])
stringers_2 = np.array([[(0.5,0.05),0.5],[(1.5,0.05),0.5]])

total_mass = get_mass(profile_1, profile_2, profile_3, pos_profile_3, stringers_1, stringers_2)
print(total_mass)