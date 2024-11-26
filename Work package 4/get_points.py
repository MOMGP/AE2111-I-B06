import math
import csv
airfoil_chords_list = 'AE2111-I-B06-Cracked-AF\Work package 4\Airfoil Chord List.csv'
import pandas as pd
import numpy as np

airfoil_geometry = []
wingbox_chords = []

def get_airfoil(input_file):
    # Read airfoil points into a numpy array
    airfoil_points = pd.read_csv(input_file).to_numpy()
    
    # Number of rows
    num_rows = airfoil_points.shape[0]
    
    # Divide into top and bottom
    mid_index = num_rows // 2
    top_points = airfoil_points[:mid_index]
    bottom_points = airfoil_points[mid_index:]
    
    # Prepare the output array
    x_y_y = np.zeros((mid_index, 3))
    
    # Fill x, top y, and bottom y values
    x_y_y[:, 0] = top_points[:, 0]  # x-coordinates
    x_y_y[:, 1] = top_points[:, 1]  # top y-coordinates
    x_y_y[:, 2] = bottom_points[:, 1]  # bottom y-coordinates
    
    return x_y_y



#np.save("Airfoil_geom.npy",get_airfoil(airfoil_chords_list))


def get_points(length_1, length_2, length_3, chord):
    airfoil_geometry = np.load('Airfoil_geom.npy') #get airfoil geometry

    length_1 = length_1/chord
    length_2 = length_2/chord
    length_3 = length_3/chord

    for i in range(len(airfoil_geometry)-1):  #input values for length 1
        j = airfoil_geometry[i][0]
        k = airfoil_geometry[i+1][0]
        if j <= length_1 < k:
            l = (length_1 - j) / (k - j)
            top_surface = (1 - l) * airfoil_geometry[i][1] + l * airfoil_geometry[i + 1][1]
            bottom_surface = (1 - l) * airfoil_geometry[i][2] + l * airfoil_geometry[i + 1][2]
            wingbox_chords.append(top_surface)
            wingbox_chords.append(bottom_surface)

        elif j <= length_2 < k:
            l = (length_2 - j) / (k - j)
            top_surface = (1 - l) * airfoil_geometry[i][1] + l * airfoil_geometry[i + 1][1]
            bottom_surface = (1 - l) * airfoil_geometry[i][2] + l * airfoil_geometry[i + 1][2]
            wingbox_chords.append(top_surface)
            wingbox_chords.append(bottom_surface)

        elif j  <= length_3 < k and not length_3 == -1:
            l = (length_3 - j) / (k - j)
            top_surface = (1 - l) * airfoil_geometry[i][1] + l * airfoil_geometry[i + 1][1]
            bottom_surface = (1 - l) * airfoil_geometry[i][2] + l * airfoil_geometry[i + 1][2]
            wingbox_chords.append(top_surface)
            wingbox_chords.append(bottom_surface)
    return wingbox_chords