from operator import length_hint

import numpy as np
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


