import numpy as np
from Geometry import get_design

spar_1, spar_2, spar_3, thickness_top_bottom, thickness_sides, end_second_cell = get_design(1)

skin_thickness = np.arange(0.001,0.01,0.001)
stringers_spars = np.arange(3,33,3)
stringers_cell1 = np.arange(2,52,2)
stringers_cell2 = np.arange(2,52,2)
number_of_ribs = np.arange(2,50,1)