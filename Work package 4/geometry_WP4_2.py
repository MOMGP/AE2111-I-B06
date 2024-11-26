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