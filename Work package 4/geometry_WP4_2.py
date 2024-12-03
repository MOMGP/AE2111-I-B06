import math
import numpy as np
from get_points import get_points, get_geom_from_points, get_airfoil
from matplotlib import pyplot as plt
# format for wing_box = [(start), (end), thickness]
#format of stringers = [pos, area]
"""
Geometry:
 ______________4_____________ ______________5_____________ 
|                            |                            |
|                            |                            |
|                            |                            |
1                            3                            6
|                            |                            |
|                            |                            |
|_____________2______________|_____________7______________| 


"""
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

def scaled_chord(spanwise_dist):
    chord = C_r - C_r * (1 - taper) * (spanwise_dist / (b / 2))
    return(chord)

def scaled_length(length,chord):
    scaled_length = length*chord/C_r
    return(scaled_length)

def centroid_distance(centroid,airfoil_points):
    distance = []
    for coordinate in airfoil_points:
        y = math.dist(centroid, coordinate)
        distance.append(y)
    max_y = max(distance)    
    return(max_y)

def centroid(wing_box, stringers):
    centroid_sum_x =0
    centroid_sum_y = 0
    area_sum = 0
    for i in range(len(wing_box)):
        start= wing_box[i][0]
        end= wing_box[i][1]
        length = np.sqrt((start[0]-end[0])**2+(start[1]-end[1])**2)
        thickness = wing_box[i][2]
        centroid_sum_x+=(start[0]+end[0])/2*length*thickness
        centroid_sum_y+=(start[1]+end[1])/2*length*thickness
        area_sum += length*thickness
    
    for i in range(len(stringers)):
        centroid_sum_x+= stringers[i][0][0]*stringers[i][1]
        centroid_sum_y+= stringers[i][0][1]*stringers[i][1]
        area_sum+=stringers[i][1]
    
    centroid_x = centroid_sum_x/area_sum
    centroid_y = centroid_sum_y/area_sum
    return centroid_x, centroid_y
        
def moments_of_inertia(wing_box, stringers):
    centroid_x, centroid_y = centroid(wing_box, stringers)
    I_xx = 0
    I_yy = 0
    I_xy = 0
    for i in range(len(wing_box)):
        start= wing_box[i][0]
        end= wing_box[i][1]
        midpoint = np.array([(start[0]+end[0])/2, (start[1]+end[1])/2])
        length = np.sqrt((start[0]-end[0])**2+(start[1]-end[1])**2)
        if start[0]!=end[0]:
            angle = np.arctan((start[1]-end[1])/(start[0]-end[0]))
        else:
            angle = np.pi/2
        thickness = wing_box[i][2]
        I_xx+=thickness*length**3*np.sin(angle)**2/12 #reg
        I_xx+=thickness*length*(midpoint[1]-centroid_y)**2
        I_yy+=thickness*length**3*np.cos(angle)**2/12 #reg
        I_yy+=thickness*length*(midpoint[0]-centroid_x)**2
    
    for i in range(len(stringers)):
        I_xx+=stringers[i][1]*(stringers[i][0][1]-centroid_y)**2
        I_yy+=stringers[i][1]*(stringers[i][0][0]-centroid_x)**2
        
        #PAT for I_xy?
        

    return I_xx, I_yy

def Lambda_n(n_percent):
    Lambda_c4 = np.rad2deg(np.arccos(1.16 / (M_CR + 0.5)))
    taper = 0.2 * (2 - np.deg2rad(Lambda_c4))
    return np.rad2deg(np.arctan(np.tan(np.deg2rad(Lambda_c4))-(4/AR)*((n_percent-25)/100*(1-taper)/(1+taper))))

def get_mass(wing_box_root, wing_box_mid, wing_box_tip, pos_mid, stringers_root):
    mass = 0

    for i in range(4): #calc mass for main box
        start_root = wing_box_root[i][0]
        end_root = wing_box_root[i][1]
        length_root = np.sqrt((start_root[0] - end_root[0]) ** 2 + (start_root[1] - end_root[1]) ** 2)
        thickness_root = wing_box_root[i][2][0]
        area_root = length_root*thickness_root
        start_tip = wing_box_tip[i] [0]
        end_tip = wing_box_tip[i][1]
        length_tip = np.sqrt((start_tip[0] - end_tip[0]) ** 2 + (start_tip[1] - end_tip[1]) ** 2)
        thickness_tip = wing_box_tip[i] [2][0]
        area_tip = length_tip * thickness_tip
        sweep = np.deg2rad((Lambda_n(start_root[0] / scaled_chord(0))+Lambda_n(end_root[0] / scaled_chord(0))) /2)
        y = b / np.cos(sweep) / np.cos(dihedral)
        mass += (area_root+area_tip)/2 * y * rho

    if len(wing_box_mid) == 7: #calc mass if the 2nd wingbox ends at random spar
        for i in range(4,7):
            start_root = wing_box_root[i][0]
            end_root = wing_box_root[i][1]
            length_root = np.sqrt((start_root[0] - end_root[0]) ** 2 + (start_root[1] - end_root[1]) ** 2)
            thickness_root = wing_box_root[i][2][0]
            area_root = length_root * thickness_root
            start_mid = wing_box_mid[i][0]
            end_mid = wing_box_mid[i][1]
            length_mid = np.sqrt((start_mid[0] - end_mid[0]) ** 2 + (start_mid[1] - end_mid[1]) ** 2)
            thickness_mid = wing_box_mid[i][2][0]
            area_mid = length_mid * thickness_mid
            sweep = np.deg2rad((Lambda_n(start_root[0] / scaled_chord(0)) + Lambda_n(end_root[0] / scaled_chord(0))) / 2)
            y = b * pos_mid / np.cos(sweep) / np.cos(dihedral)
            mass += (area_root + area_mid) / 2 * y * rho

    elif len(wing_box_mid) == 4: #calc mass if the 2nd wingbox merges
        for i in range(4,7):
            start_root = wing_box_root[i][0]
            end_root = wing_box_root[i][1]
            length_root = np.sqrt((start_root[0] - end_root[0]) ** 2 + (start_root[1] - end_root[1]) ** 2)
            thickness_root = wing_box_root[i][2][0]
            area_root = length_root * thickness_root
            start_mid = wing_box_mid[i][0]
            end_mid = wing_box_mid[i][1]
            length_mid = np.sqrt((start_mid[0] - end_mid[0]) ** 2 + (start_mid[1] - end_mid[1]) ** 2)
            thickness_mid = wing_box_mid[i][2][0]
            area_mid = length_mid * thickness_mid
            if i == 4 or i == 5:
                area_mid = 0
            sweep = np.deg2rad((Lambda_n(start_root[0] / scaled_chord(0)) + Lambda_n(end_root[0] / scaled_chord(0))) / 2)
            y = b * pos_mid / np.cos(sweep) / np.cos(dihedral)
            mass += (area_root + area_mid) / 2 * y * rho

    for i in range(len(stringers_root)): #calc mass for spars in main wing box

        start_root = stringers_root[i][0]
        area_root = stringers_root[i][1]

        sweep = np.deg2rad((Lambda_n(start_root[0] / scaled_chord(0))) / 2)
        y = b / np.cos(sweep) / np.cos(dihedral)
        mass += area_root[0] * y * rho

    return mass

# format for wing_box = [(start), (end), thickness]
#format of stringers = [pos, area]
#TODO - this currently wrong - make it for lists
def get_points_along_spanwise(norm_wing_box_root, norm_stringers, y, end_third_spar, trunctated=False):
    chord = scaled_chord(y)
    geometry = []
    stringers = []
    
    if trunctated:
        if y <= end_third_spar:
            for i in norm_wing_box_root:
                plate = ([[i[0][0]*chord, i[0][1]*chord], [i[1][0]*chord, i[1][1]*chord], i[2]])
                geometry.append(plate)
        else:
            for i in range(4):
                plate = [[norm_wing_box_root[i][0][0]*chord,norm_wing_box_root[i][0][1]*chord], [norm_wing_box_root[i][1][0]*chord, norm_wing_box_root[i][1][1]*chord], norm_wing_box_root[i][2]]
                geometry.append(plate)
    else:
        if y <= end_third_spar:
            for i in range(4):
                plate = [[norm_wing_box_root[i][0][0]*chord, norm_wing_box_root[i][0][1]*chord], [norm_wing_box_root[i][1][0]*chord, norm_wing_box_root[i][1][1]*chord], norm_wing_box_root[i][2]]
                geometry.append(plate)
            x_pos_third_spar = (norm_wing_box_root[6][1][0]*(end_third_spar - y) + norm_wing_box_root[2][0][0]*y) / end_third_spar
            x_y_y = get_points(-1, -1, x_pos_third_spar, 1) 
            
            geometry.append([[norm_wing_box_root[2][0][0]*chord,norm_wing_box_root[2][0][1]*chord], [x_y_y[0][0]* chord, x_y_y[0][1]* chord], norm_wing_box_root[4][2]])
            geometry.append([[x_y_y[0][0]* chord, x_y_y[0][1]* chord], [x_y_y[1][0]* chord, x_y_y[1][1]* chord], norm_wing_box_root[5][2]])
            geometry.append([[norm_wing_box_root[3][0][0]*chord,norm_wing_box_root[3][0][1]*chord], [x_y_y[0][0]* chord, x_y_y[0][1]* chord], norm_wing_box_root[6][2]])
            
        else:
            for i in range(4):
                plate = [[norm_wing_box_root[i][0][0]*chord,norm_wing_box_root[i][0][1]*chord], [norm_wing_box_root[i][1][0]*chord,norm_wing_box_root[i][1][1]*chord], norm_wing_box_root[i][2]]
                geometry.append(plate)

    for i in norm_stringers:
        stringers.append([[i[0][0]*chord,i[0][1]*chord], i[1]]) 

    return geometry, stringers

def get_stringer_geom_norm(geom, num_of_stringers_tot):
    num_of_stringers_side = int(num_of_stringers_tot/2)
    top_plate_arr = np.array([geom[3][0], geom[3][1]])
    direction_top = top_plate_arr[1]-top_plate_arr[0]
    bot_plate_arr = np.array([geom[1][0], geom[1][1]])
    direction_bot = bot_plate_arr[1]-bot_plate_arr[0]
    stringers = [[[top_plate_arr[0][0]+direction_top[0]*str_num_top/(num_of_stringers_side+1), 
                  top_plate_arr[0][1]+direction_top[1]*str_num_top/(num_of_stringers_side+1)],stringer_area] for str_num_top in range(1, num_of_stringers_side+1)]
    for str_num in range(1, num_of_stringers_side+1):
        stringers.append([[bot_plate_arr[0][0]+direction_bot[0]*str_num/(num_of_stringers_side+1), 
        bot_plate_arr[0][1]+direction_bot[1]*str_num/(num_of_stringers_side+1)], stringer_area])
    return stringers
    

def plot3d_geom(geom_root, stringers, end_third_spar, truncated=False, plot_with_airfoil=False, full_wing=False, plot_stringers = True):
    points_root, stringers_root = get_points_along_spanwise(geom_root, stringers, 0, end_third_spar, trunctated=truncated)
    points_end_third_spar, _ = get_points_along_spanwise(geom_root, stringers, end_third_spar, end_third_spar, trunctated=truncated)
    # print(points_end_third_spar)
    points_tip, stringers_tip = get_points_along_spanwise(geom_root, stringers, b/2, end_third_spar, trunctated=truncated)
    lambda_LE = np.deg2rad(Lambda_n(0))

    # x_vals = []
    # y_vals = []
    # for i in points_end_third_spar:
    #     x_vals.append(i[0][0])
    #     x_vals.append(i[1][0])
    #     y_vals.append(i[0][1])
    #     y_vals.append(i[1][1])
    # plt.plot(x_vals, y_vals)
    # plt.show()
    
    front_spar_x = np.array([[points_root[0][0][0], points_root[0][1][0]], [points_tip[0][0][0] + np.sin(lambda_LE) * b / 2, points_tip[0][1][0] + np.sin(lambda_LE) * b / 2]])
    front_spar_y = np.array([[0, 0], [b/2, b/2]])
    front_spar_z = np.array([[points_root[0][0][1], points_root[0][1][1]],[points_tip[0][0][1] + np.sin(np.deg2rad(dihedral)) * b / 2, points_tip[0][1][1] + np.sin(np.deg2rad(dihedral)) * b / 2]])
    
    snd_spar_x = np.array([[points_root[2][0][0], points_root[2][1][0]], [points_tip[2][0][0] + np.sin(lambda_LE) * b / 2, points_tip[2][1][0] + np.sin(lambda_LE) * b / 2]])
    snd_spar_y = np.array([[0, 0], [b/2, b/2]])
    snd_spar_z = np.array([[points_root[2][0][1], points_root[2][1][1]], [points_tip[2][0][1] + np.sin(np.deg2rad(dihedral)) * b / 2, points_tip[2][1][1] + np.sin(np.deg2rad(dihedral)) * b / 2]])


    thd_spar_x = np.array([[points_root[5][0][0], points_root[5][1][0]], [points_end_third_spar[5][0][0] + np.sin(lambda_LE) * end_third_spar, points_end_third_spar[5][1][0] + np.sin(lambda_LE) * end_third_spar]])
    thd_spar_y = np.array([[0, 0], [end_third_spar, end_third_spar]])
    thd_spar_z = np.array([[points_root[5][0][1], points_root[5][1][1]],[points_end_third_spar[5][0][1] + np.sin(np.deg2rad(dihedral)) * end_third_spar, points_end_third_spar[5][1][1] + np.sin(np.deg2rad(dihedral)) * end_third_spar]])
    # Ensure the axes are properly set for 3D plotting
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    if plot_with_airfoil:
        x_y_y=np.load("Airfoil_geom.npy")
        root_airfoil = x_y_y*scaled_chord(0)
        tip_airfoil = x_y_y*scaled_chord(b/2)

        tip_airfoil[:,0] += np.sin(lambda_LE) * b / 2
        tip_airfoil[:,1] += np.sin(np.deg2rad(dihedral)) * b / 2
        tip_airfoil[:,2] += np.sin(np.deg2rad(dihedral)) * b / 2


        all_x = np.tile(np.vstack((root_airfoil[:,0], tip_airfoil[:,0])).T, [2,1])
        all_y = np.array([0, b/2])
        # print(all_x.shape)
        # print(all_y.shape)
    
        all_Z_upper = np.vstack((root_airfoil[:, 1], tip_airfoil[:, 1])).T
        all_Z_lower = np.vstack((root_airfoil[:, 2], tip_airfoil[:, 2])).T
        # print(all_Z_lower.shape)
        # Tile the Z values to match the meshgrid shape
        Z_upper = np.tile(all_Z_upper, (len(all_y), 1))
        Z_lower = np.tile(all_Z_lower, (len(all_y), 1))

        # Plot root airfoil
        ax.plot_surface(all_x, all_y, Z_upper, color='red', alpha=0.2)
        ax.plot_surface(all_x, all_y, Z_lower, color='red', alpha=0.2)
        if full_wing:
            ax.plot_surface(all_x, -all_y, Z_upper, color='red', alpha=0.2)
            ax.plot_surface(all_x, -all_y, Z_lower, color='red', alpha=0.2)

    
    # Plot surfaces for front and second spar
    ax.plot_surface(front_spar_x, front_spar_y, front_spar_z)
    ax.plot_surface(snd_spar_x, snd_spar_y, snd_spar_z)
    ax.plot_surface(thd_spar_x, thd_spar_y, thd_spar_z)

    #plot stringers
    if plot_stringers:
        for i in range(len(stringers_root)):
            ax.plot([stringers_root[i][0][0], stringers_tip[i][0][0]+np.sin(lambda_LE) * b / 2], [0, b/2], [stringers_root[i][0][1], stringers_tip[i][0][1]+np.sin(np.deg2rad(dihedral)) * b / 2], linewidth=0.3)
    
    if full_wing:
        ax.plot_surface(front_spar_x, -front_spar_y, front_spar_z)
        ax.plot_surface(snd_spar_x, -snd_spar_y, snd_spar_z)
        ax.plot_surface(thd_spar_x, -thd_spar_y, thd_spar_z)
        if plot_stringers:
            for i in range(len(stringers_root)):
                ax.plot([stringers_root[i][0][0], stringers_tip[i][0][0]+np.sin(lambda_LE) * b / 2], [0, -b/2], [stringers_root[i][0][1], stringers_tip[i][0][1]+np.sin(np.deg2rad(dihedral)) * b / 2], linewidth=0.3)
        ax.set_xlim3d(-25, 25)
        ax.set_ylim3d(-25, 25)
        ax.set_zlim3d(-25, 25)

    else:
        ax.set_xlim3d(0, 25)
        ax.set_ylim3d(0, 25)
        ax.set_zlim3d(0, 25)
        
    # plt.show()

profile_root = [
    [[0, 0], [0, 1], [0.01]],
    [[0, 1], [1, 1], [0.01]],
    [[1, 1], [1, 0], [0.01]],
    [[1, 0], [0, 0], [0.01]],
    [[1, 1], [2, 1], [0.01]],
    [[2, 1], [2, 0], [0.01]],
    [[2, 0], [1, 0], [0.01]]
]

profile_mid = [
    [[0, 0], [0, 1], [0.01]],
    [[0, 1], [1, 1], [0.01]],
    [[1, 1], [1, 0], [0.01]],
    [[1, 0], [0, 0], [0.01]],
    [[1, 1], [2, 1], [0.01]],
    [[2, 1], [2, 0], [0.01]],
    [[2, 0], [1, 0], [0.01]]
]

profile_tip = [
    [[0, 0], [0, 1], [0.01]],
    [[0, 1], [1, 1], [0.01]],
    [[1, 1], [1, 0], [0.01]],
    [[1, 0], [0, 0], [0.01]]
]

pos_profile_3 = 0.5

stringers_root = [
    [[0.5, 0.05], [0.01]],
    [[1.5, 0.05], [0.01]],
    [[0.5, 0.95], [0.01]],
    [[1.5, 0.95], [0.01]]
]


total_mass = get_mass(profile_root, profile_mid, profile_tip, pos_profile_3, stringers_root)
# print(total_mass)

spar1_x=0.2
spar2_x=0.5
spar3_x=0.7
x_y_y = get_points(spar1_x, spar2_x, spar3_x, 1)
root_geom = get_geom_from_points(x_y_y, [.644 for i in range(7)])
stringers_norm = get_stringer_geom_norm(root_geom, 32)
# print(root_geom)
x_vals = []
y_vals = []
# for i in root_geom:
#     x_vals.append(i[0][0])
#     x_vals.append(i[1][0])
#     y_vals.append(i[0][1])
#     y_vals.append(i[1][1])
# airfoil_xyy = np.load("Airfoil_geom.npy")
# plt.plot(airfoil_xyy[:,0], airfoil_xyy[:,1])
# plt.plot(airfoil_xyy[:,0], airfoil_xyy[:,2])
# plt.plot(x_vals, y_vals)
# plt.show()

plot3d_geom(root_geom, stringers_norm, b/6, truncated=False, plot_with_airfoil=True, full_wing=True)

root_pts, root_stringers = get_points_along_spanwise(root_geom, stringers_norm, 0, b/6)
print(moments_of_inertia(root_pts, root_stringers)[0])
