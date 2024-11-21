import math
import csv
airfoil_chords_list = 'AE2111-I-B06-Cracked-AF\Work package 4\Airfoil Chord List.csv'
import pandas as pd
import numpy as np

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



np.save("Airfoil_geom.npy",get_airfoil(airfoil_chords_list))



