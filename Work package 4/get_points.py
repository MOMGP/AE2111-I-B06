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


def get_points(length_1, length_2, length_3):
    airfoil_geometry = np.load('C:\\Users\\koppe\\PycharmProjects\\AE2111-I-B06-Cracked-AF\\Airfoil_geom.npy') #get airfoil geometry

    for i in range(len(airfoil_geometry)-1):  #input values for length 1
        j = airfoil_geometry[i][0]
        k = airfoil_geometry[i+1][0]
        if j <= length_1 <= k:
            l = (length_1 - j) / (k - j)
            top_surface = (1 - l) * airfoil_geometry[i][1] + l * airfoil_geometry[i + 1][1]
            bottom_surface = (1 - l) * airfoil_geometry[i][2] + l * airfoil_geometry[i + 1][2]
            wingbox_chords.append(top_surface)
            wingbox_chords.append(bottom_surface)

        elif j <= length_2 <= k:
            l = (length_2 - j) / (k - j)
            top_surface = (1 - l) * airfoil_geometry[i][1] + l * airfoil_geometry[i + 1][1]
            bottom_surface = (1 - l) * airfoil_geometry[i][2] + l * airfoil_geometry[i + 1][2]
            wingbox_chords.append(top_surface)
            wingbox_chords.append(bottom_surface)

        elif j  <= length_3 <= k and not length_3 == -1:
            l = (length_3 - j) / (k - j)
            top_surface = (1 - l) * airfoil_geometry[i][1] + l * airfoil_geometry[i + 1][1]
            bottom_surface = (1 - l) * airfoil_geometry[i][2] + l * airfoil_geometry[i + 1][2]
            wingbox_chords.append(top_surface)
            wingbox_chords.append(bottom_surface)

    intersect = False #set intersect false

    for i in range(0, len(wingbox_chords)-2):
        y1 = wingbox_chords[i]
        y2 = wingbox_chords[i+2]
        if i <= 1:
            x1 = length_1
            x2 = length_2
            a = (x2 - x1) / (y2 - y1)
            b = y1 - a * x1
            x = np.arange(np.round(length_1, 2), np.round(length_2, 2), 0.01)
            for j in range(len(x)):
                y = a * x[j] + b
                yref1 = airfoil_geometry[np.where(airfoil_geometry[:,0] == x[j])][1]
                yref2 = airfoil_geometry[np.where(airfoil_geometry[:,0] == x[j])][2]
                if not yref2 < y <yref1:
                    intersect = True


        elif i >= 2:
            x1 = length_2
            x2 = length_3
            a = (x2 - x1) / (y2 - y1)
            b = y1 - a * x1
            x = np.arange(np.round(length_2, 2), np.round(length_3, 2), 0.01)
            for j in range(x):
                y = a * x[j] + b
                yref1 = airfoil_geometry[np.where(airfoil_geometry[:,0] == x[j])][1]
                yref2 = airfoil_geometry[np.where(airfoil_geometry[:,0] == x[j])][2]
                if not yref2 < y < yref1:
                    intersect = True

    return wingbox_chords(), intersect

winbox_chords, intersect  = get_points(0.2, 0.5, -1)
print(winbox_chords)
print(intersect)
