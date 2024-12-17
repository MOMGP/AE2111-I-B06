import numpy as np
from Geometry import get_design
from Divide_Geom import design

geometry = get_design(0)

skin_thickness = np.arange(0.001,0.01,0.001)
stringers_spars = np.arange(3,33,3)
stringers_cell1 = np.arange(2,52,2)
stringers_cell2 = np.arange(2,52,2)
number_of_ribs = np.arange(3,50,1)

stringers = [stringers_spars[0], stringers_cell1[0], stringers_cell2[0]]

list = design(geometry, skin_thickness[0], stringers, number_of_ribs[0])

print(list)