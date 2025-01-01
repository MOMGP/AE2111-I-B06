from operator import length_hint

import numpy as np
from matplotlib import pyplot as plt
from numpy.ma.core import divide
from pandas.core.interchange.from_dataframe import primitive_column_to_ndarray

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
stringer_base = 0.037
stringer_length = 0.052
stringer_thickness = 0.0023



geometry = []
b = 53.57



def scaled_chord(spanwise_dist):
    chord = C_r - C_r * (1 - taper) * (spanwise_dist / (b / 2))
    return(chord)


def get_design(design_number):
    # define variables for function
    spar1 = 0
    spar2 = 0
    spar3 = 0
    thickness_sides = 0
    thickness_top_bottom = 0
    #minimum weight design
    if design_number == 0:          #first design
        spar1 = 0.2
        spar2 = 0.35
        spar3 = 0.65
        thickness_sides = 0.008
        thickness_top_bottom = 0.014
    #buckling design
    elif design_number == 1:        #second design
        spar1 = 0.2
        spar2 = 0.4
        spar3 = 0.7
        thickness_sides = 0.01
        thickness_top_bottom = 0.016
    #ultimate load design
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

def weight(divided_geometry):
    weight = 0
    airfoil_geometry = np.load('Airfoil_geom.npy')  # get airfoil geometry
    for i in range(len(divided_geometry)):
        print(divided_geometry[i])
        specific = divided_geometry[i][0]
        if specific == "spar":

            #get height
            height = 0
            for n in range(len(airfoil_geometry) - 1):  # input values for length 1
                j = airfoil_geometry[n][0]
                k = airfoil_geometry[n + 1][0]
                if j <= divided_geometry[i][1][1] < k:
                    l = (divided_geometry[i][1][1] - j) / (k - j)
                    top = (1 - l) * airfoil_geometry[n][1] + l * airfoil_geometry[n + 1][1]
                    bottom = (1 - l) * airfoil_geometry[n][2] + l * airfoil_geometry[n + 1][2]

                    #scale height with chord
                    height = (top-bottom)*scaled_chord((divided_geometry[i][1][0]+divided_geometry[i][2][0])/2)

            #add spar weight
            width = divided_geometry[i][2][0]-divided_geometry[i][1][0]
            thickness = divided_geometry[i][3]
            #weight += width * height * thickness * rho

            #add stringer weight
            #weight += stringer_area * width * rho * divided_geometry[i][4]

        elif specific == "flange":

            #calculate chord
            average_chord = (scaled_chord(divided_geometry[i][1][0])+scaled_chord(divided_geometry[i][2][0]))/2
            chord = (divided_geometry[i][2][1]-divided_geometry[i][1][1]) * average_chord

            #add flange weight
            width = divided_geometry[i][2][0] - divided_geometry[i][1][0]
            thickness = divided_geometry[i][3]
            weight += width * chord * thickness * rho

            #add stringer weight
            #weight += stringer_area * width * rho * divided_geometry[i][4]

        elif specific == "skin":
            #todo add skin weight
            weight += 0
    weight = 2*weight
    return weight

path_airfoil = "Airfoil_geom.npy"

def centroid_z(type, x_over_c1, x_over_c2, y_end, n_stringer, thickness, top=False):
    Q_z=0
    area_sum = 0
    if type == "spar":
        return 0
    elif type == "flange" or type =="skin":
        airfoil_geometry = np.load("Airfoil_geom.npy")
        starting_index = -1
        ending_index = -1
        for i in range(len(airfoil_geometry)-1):  #input values for length 1
            j = airfoil_geometry[i][0]
            k = airfoil_geometry[i+1][0]
            if j<x_over_c1<=k:
                starting_index = i
            if j<x_over_c2<=k:
                ending_index = i
                break
        if top:
            angle = np.arctan((airfoil_geometry[ending_index][1]-airfoil_geometry[starting_index][1])/(airfoil_geometry[ending_index][0]-airfoil_geometry[starting_index][0]))
            stringers_placed = 0
            for i in range(starting_index, ending_index): #no plus one for ending index because we add one below, so we go through all panels. 
                x_mid = (airfoil_geometry[i][0]+airfoil_geometry[i+1][0])*scaled_chord(y_end)/2-airfoil_geometry[starting_index][0]*scaled_chord(y_end)
                y_mid = (airfoil_geometry[i][1]+airfoil_geometry[i+1][1])*scaled_chord(y_end)/2-airfoil_geometry[starting_index][1]*scaled_chord(y_end)
                Q_z+= np.sqrt(((airfoil_geometry[i][0]-airfoil_geometry[i+1][0])*scaled_chord(y_end))**2+((airfoil_geometry[i][1]-airfoil_geometry[i+1][1])*scaled_chord(y_end))**2)*thickness*np.sqrt((x_mid-x_mid*np.cos(angle))**2+(y_mid-x_mid*np.sin(angle))**2)    
                area_sum+= np.sqrt(((airfoil_geometry[i][0]-airfoil_geometry[i+1][0])*scaled_chord(y_end))**2+((airfoil_geometry[i][1]-airfoil_geometry[i+1][1])*scaled_chord(y_end))**2)*thickness
                if stringers_placed<(i-starting_index)/(ending_index-starting_index)*(n_stringer+1) and n_stringer!=0:
                    area_sum+=stringer_base*stringer_thickness+(stringer_length-stringer_thickness)*stringer_thickness
                    Q_z+=stringer_base*stringer_thickness*np.sqrt((x_mid-x_mid*np.cos(angle))**2+(y_mid-thickness/2-stringer_thickness/2-x_mid*np.sin(angle))**2)         #stringer base contribution
                    Q_z+=(stringer_length-stringer_thickness)*stringer_thickness*np.sqrt((x_mid-x_mid*np.cos(angle))**2+(y_mid-thickness/2-stringer_thickness-stringer_length/2-x_mid*np.sin(angle))**2)*(y_mid-thickness/2-stringer_thickness-stringer_length/2-x_mid*np.sin(angle))/np.abs((y_mid-thickness/2-stringer_thickness-stringer_length/2-x_mid*np.sin(angle))) #stringer web contribution
                    stringers_placed+=1
                                
        else:
            angle = np.arctan((airfoil_geometry[ending_index][2]-airfoil_geometry[starting_index][2])/(airfoil_geometry[ending_index][0]-airfoil_geometry[starting_index][0]))
            stringers_placed = 0
            for i in range(starting_index, ending_index):
                x_mid = (airfoil_geometry[i][0]+airfoil_geometry[i+1][0])*scaled_chord(y_end)/2-airfoil_geometry[starting_index][0]*scaled_chord(y_end)
                y_mid = (airfoil_geometry[i][2]+airfoil_geometry[i+1][2])*scaled_chord(y_end)/2-airfoil_geometry[starting_index][2]*scaled_chord(y_end)
                Q_z+= np.sqrt(((airfoil_geometry[i][0]-airfoil_geometry[i+1][0])*scaled_chord(y_end))**2+((airfoil_geometry[i][2]-airfoil_geometry[i+1][2])*scaled_chord(y_end))**2)*thickness*np.sqrt((x_mid-x_mid*np.cos(angle))**2+(y_mid-x_mid*np.sin(angle))**2)*(y_mid-x_mid*np.sin(angle))/np.abs(y_mid-x_mid*np.sin(angle))
                area_sum+= np.sqrt(((airfoil_geometry[i][0]-airfoil_geometry[i+1][0])*scaled_chord(y_end))**2+((airfoil_geometry[i][2]-airfoil_geometry[i+1][2])*scaled_chord(y_end))**2)*thickness

                if stringers_placed<(i-starting_index)/(ending_index-starting_index)*(n_stringer+1) and n_stringer!=0:
                    area_sum+=(stringer_base*stringer_thickness+(stringer_length-stringer_thickness)*stringer_thickness)
                    # Q_z+=stringer_base*stringer_thickness*((y_mid+thickness/2+stringer_thickness/2)*np.cos(angle)+x_mid*np.sin(angle))
                    # Q_z+=(stringer_length-stringer_thickness)*stringer_thickness*((y_mid+thickness/2-stringer_thickness+stringer_length/2)*np.cos(angle)+x_mid*np.sin(angle))
                    Q_z+=stringer_base*stringer_thickness*(np.sqrt((x_mid-x_mid*np.cos(angle))**2+(y_mid-x_mid*np.sin(angle))**2)+thickness/2+stringer_thickness/2)         #stringer base contribution
                    Q_z+=(stringer_length-stringer_thickness)*stringer_thickness*(np.sqrt((x_mid-x_mid*np.cos(angle))**2+(y_mid-x_mid*np.sin(angle))**2)+thickness/2+stringer_thickness+stringer_length/2) #stringer web contribution. Latter part of the term is to determine whether contributes positively or negatively
                    stringers_placed+=1
    else:
        print("not a valid type")
        return 0
    return Q_z/area_sum
    

def centroid_z_test(type, x_over_c1, x_over_c2, y_end, n_stringer, thickness, top=False):
    centroid_val = centroid_z(type, x_over_c1, x_over_c2, y_end, n_stringer, thickness, top)
    if type == "flange" or type =="skin":
        airfoil_geometry = np.load("Airfoil_geom.npy")
        starting_index = -1
        ending_index = -1
        for i in range(len(airfoil_geometry)-1):  #input values for length 1
            j = airfoil_geometry[i][0]
            k = airfoil_geometry[i+1][0]
            if j<x_over_c1<=k:
                starting_index = i
            if j<x_over_c2<=k:
                ending_index = i
                break
        if top:
            angle = np.arctan((airfoil_geometry[ending_index][1]-airfoil_geometry[starting_index][1])/(airfoil_geometry[ending_index][0]-airfoil_geometry[starting_index][0]))
            plt.plot([(airfoil_geometry[starting_index][0])*scaled_chord(y_end)+np.sin(angle)*centroid_val, (airfoil_geometry[ending_index][0])*scaled_chord(y_end)+ np.sin(angle)*centroid_val ], [(airfoil_geometry[starting_index][1])*scaled_chord(y_end)+ np.cos(angle)*centroid_val, (airfoil_geometry[ending_index][1])*scaled_chord(y_end) + np.cos(angle)*centroid_val])
            stringers_placed = 0
            for i in range(starting_index, ending_index): #no plus one for ending index because we add one below, so we go through all panels. 
                x_mid = (airfoil_geometry[i][0]+airfoil_geometry[i+1][0])*scaled_chord(y_end)/2-airfoil_geometry[starting_index][0]*scaled_chord(y_end)
                y_mid = (airfoil_geometry[i][1]+airfoil_geometry[i+1][1])*scaled_chord(y_end)/2-airfoil_geometry[starting_index][1]*scaled_chord(y_end)
                plt.plot([ (airfoil_geometry[i][0])*scaled_chord(y_end), airfoil_geometry[i+1][0]*scaled_chord(y_end)], [airfoil_geometry[i][1]*scaled_chord(y_end),airfoil_geometry[i+1][1]*scaled_chord(y_end)])

                if stringers_placed<(i-starting_index)/(ending_index-starting_index)*(n_stringer+1) and n_stringer!=0:
                    # plot stringers (or not, honestly idc, it looks fine without them)
                    stringers_placed+=1

            plt.show()
                    
                                
        else:
            angle = np.arctan((airfoil_geometry[ending_index][2]-airfoil_geometry[starting_index][2])/(airfoil_geometry[ending_index][0]-airfoil_geometry[starting_index][0]))
            plt.plot([(airfoil_geometry[starting_index][0])*scaled_chord(y_end)+np.sin(angle)*centroid_val, (airfoil_geometry[ending_index][0])*scaled_chord(y_end)+ np.sin(angle)*centroid_val], [(airfoil_geometry[starting_index][2])*scaled_chord(y_end)+ np.cos(angle)*centroid_val, (airfoil_geometry[ending_index][2])*scaled_chord(y_end) + np.cos(angle)*centroid_val])
            stringers_placed = 0
            for i in range(starting_index, ending_index): #no plus one for ending index because we add one below, so we go through all panels. 
                x_mid = (airfoil_geometry[i][0]+airfoil_geometry[i+1][0])*scaled_chord(y_end)/2-airfoil_geometry[starting_index][0]*scaled_chord(y_end)
                y_mid = (airfoil_geometry[i][2]+airfoil_geometry[i+1][2])*scaled_chord(y_end)/2-airfoil_geometry[starting_index][2]*scaled_chord(y_end)
                plt.plot([ (airfoil_geometry[i][0])*scaled_chord(y_end), airfoil_geometry[i+1][0]*scaled_chord(y_end)], [airfoil_geometry[i][2]*scaled_chord(y_end),airfoil_geometry[i+1][2]*scaled_chord(y_end)])

                if stringers_placed<(i-starting_index)/(ending_index-starting_index)*(n_stringer+1) and n_stringer!=0:
                    #todo - plot stringers
                    stringers_placed+=1

            plt.show()
    else:
        print("no test for that type")



def I_xx(type, x_over_c1, x_over_c2, y_end, n_stringer, thickness, top=False): #Type can be 'spar', 'flange', or 'skin'
    I_xx = 0
    centroid_z_loc = centroid_z(type, x_over_c1, x_over_c2, y_end, n_stringer, thickness, top)
    if type=="spar":
        airfoil_geometry = np.load("Airfoil_geom.npy")
        index = -1
        for i in range(len(airfoil_geometry)-1):  #input values for length 1
            j = airfoil_geometry[i][0]
            k = airfoil_geometry[i+1][0]
            if j<x_over_c1<=k:
                index = i
                break
        top_z = airfoil_geometry[index][1]*scaled_chord(y_end)
        bot_z = airfoil_geometry[index][2]*scaled_chord(y_end)
        I_xx+=(top_z-bot_z)**3*thickness/12
        for i in range(n_stringer):
            I_xx+= stringer_base**3*stringer_thickness/12+ ((i+1)*(top_z-bot_z)/(n_stringer+1)+bot_z)**2*stringer_thickness*stringer_base #flange term
            I_xx+= ((i+1)*(top_z-bot_z)/(n_stringer+1)+bot_z- stringer_base/2)**2*stringer_thickness*stringer_length#web term - assuming thin so i only compute PAT

    elif type == "flange" or "skin":
        airfoil_geometry = np.load("Airfoil_geom.npy")
        starting_index = -1
        ending_index = -1
        for i in range(len(airfoil_geometry)-1):  #input values for length 1
            j = airfoil_geometry[i][0]
            k = airfoil_geometry[i+1][0]
            if j<x_over_c1<=k:
                starting_index = i
            if j<x_over_c2<=k:
                ending_index = i
                break
        if top:
            angle = np.arctan((airfoil_geometry[ending_index][1]-airfoil_geometry[starting_index][1])/(airfoil_geometry[ending_index][0]-airfoil_geometry[starting_index][0]))
            stringers_placed = 0
            for i in range(starting_index, ending_index): #no plus one for ending index because we add one below, so we go through all panels. 
                x_mid = (airfoil_geometry[i][0]+airfoil_geometry[i+1][0])*scaled_chord(y_end)/2-airfoil_geometry[starting_index][0]*scaled_chord(y_end)
                y_mid = (airfoil_geometry[i][1]+airfoil_geometry[i+1][1])*scaled_chord(y_end)/2-airfoil_geometry[starting_index][1]*scaled_chord(y_end)
                local_length = np.sqrt((airfoil_geometry[i][0]-airfoil_geometry[i+1][0])**2+(airfoil_geometry[i][1]-airfoil_geometry[i+1][1])**2)*scaled_chord(y_end)
                local_angle = np.arctan((airfoil_geometry[i][1]-airfoil_geometry[i+1][1])/(airfoil_geometry[i][0]-airfoil_geometry[i+1][0]))
                displacement = np.sqrt((x_mid-x_mid*np.cos(angle))**2+(y_mid-x_mid*np.sin(angle))**2) * (y_mid-x_mid*np.sin(angle))/np.abs(y_mid-x_mid*np.sin(angle))-centroid_z_loc #todo - this is wrong, missing the centroid?!?!?
                I_xx+= local_length**3*thickness*np.sin(local_angle)**2/12+displacement**2*local_length*thickness

                if stringers_placed<(i-starting_index)/(ending_index-starting_index)*(n_stringer+1) and n_stringer!=0:
                    I_xx+=stringer_base**3*stringer_thickness*np.sin(local_angle)**2/12 + (displacement-thickness/2-stringer_thickness/2)**2*stringer_base*stringer_thickness
                    I_xx+=stringer_length**3*stringer_thickness*np.cos(local_angle)**2/12+ (displacement -thickness/2-stringer_thickness-stringer_length/2)**2*stringer_length*stringer_thickness
                    stringers_placed+=1
                                
        else:
            angle = np.arctan((airfoil_geometry[ending_index][2]-airfoil_geometry[starting_index][2])/(airfoil_geometry[ending_index][0]-airfoil_geometry[starting_index][0]))
            stringers_placed = 0
            for i in range(starting_index, ending_index):
                x_mid = (airfoil_geometry[i][0]+airfoil_geometry[i+1][0])*scaled_chord(y_end)/2-airfoil_geometry[starting_index][0]*scaled_chord(y_end)
                y_mid = (airfoil_geometry[i][2]+airfoil_geometry[i+1][2])*scaled_chord(y_end)/2-airfoil_geometry[starting_index][2]*scaled_chord(y_end)
                local_length = np.sqrt((airfoil_geometry[i][0]-airfoil_geometry[i+1][0])**2+(airfoil_geometry[i][2]-airfoil_geometry[i+1][2])**2)*scaled_chord(y_end)
                local_angle = np.arctan((airfoil_geometry[i][2]-airfoil_geometry[i+1][2])/(airfoil_geometry[i][0]-airfoil_geometry[i+1][0]))
                displacement = np.sqrt((x_mid-x_mid*np.cos(angle))**2+(y_mid-x_mid*np.sin(angle))**2) * (y_mid-x_mid*np.sin(angle))/np.abs(y_mid-x_mid*np.sin(angle))-centroid_z_loc #todo - this is wrong, missing the centroid?!?!?
                I_xx+= local_length**3*thickness*np.sin(local_angle)**2/12+displacement**2*local_length*thickness

                if stringers_placed<(i-starting_index)/(ending_index-starting_index)*(n_stringer+1) and n_stringer!=0:
                    I_xx+=stringer_base**3*stringer_thickness*np.sin(local_angle)**2/12 + (displacement+thickness/2+stringer_thickness/2)**2*stringer_base*stringer_thickness
                    I_xx+=stringer_length**3*stringer_thickness*np.cos(local_angle)**2/12+ (displacement +thickness/2+stringer_thickness+stringer_length/2)**2*stringer_length*stringer_thickness
                    stringers_placed+=1

    else:
        print("centroid_z no valid type")
        return None
    return I_xx 








# MAIN TESTING PLACE
# centroid_z_test("flange", 0.2, 0.35, 3, 4, 0.02, False)
# centroid_z_test("flange", 0.2, 0.35, 3, 4, 0.02, True)
# centroid_z_test("flange", 0.2, 0.35, 20, 4, 0.02, False)
# centroid_z_test("flange", 0.2, 0.35, 20, 4, 0.02, True)
"centroid checks out"

# print(I_xx("spar", 0.3, 0.3, 0, 0, 1, False)) # Only spar - Correct
# print(I_xx("spar", 0.3, 0.3, 0, 1, 0, False)) # Only stringer - Checks out
# print(I_xx("flange", 0.2, 0.35, 0, 0, 0.02, True)) #Looks correct
# print(I_xx("flange", 0.2, 0.35, 0, 0, 0.02, False)) #checks out
# print(I_xx("flange", 0.2, 0.35, 0, 1, 0, False)) #seems right
"I_xx appears correct"


