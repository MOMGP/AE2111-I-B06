import numpy as np
from Geometry import get_design
from Divide_Geom import design

geometry = get_design(0)

#longer list if iteration doesnt take all day
skin_thickness = [2.05, 1.8, 1.63, 1.4, 1.29, 1.1, 1.02, 0.91, 0.81, 0.71, 0.64, 0.58, 0.51, 0.46, 0.43]

#shorter list if iteration does take all day
#skin_thickness = [2.05, 1.8, 1.4, 1.1, 0.81, 0.64, 0.51, 0.43]

stringers_spars = np.arange(5,10,1) #is per surface, so total is *3 per half wing
stringers_cell1 = np.arange(9,20,1) #is per surface, so total is *2 per half wing
stringers_cell2 = np.arange(14,30,1) #is per surface, so total is *2 per half wing
number_of_ribs = np.arange(7,20,1) #is per half wing, so total is *2

stringers = [stringers_spars[0], stringers_cell1[0], stringers_cell2[0]]

list = design(geometry, skin_thickness[0], stringers, number_of_ribs[0])

print(list)