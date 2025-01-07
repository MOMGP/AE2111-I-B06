import numpy as np
from Geometry import get_design, scaled_chord, I_xx_global
from Divide_Geom import design
from get_points import get_points


designs = np.load("Work package 5\\Results\\Design_2.npy")
filtered_designs = designs[designs[:, 1] == 2]
designs_ult = np.load("Work package 5\\Results\\Design_ulti.npy")
index = np.argmin(filtered_designs, axis=0)
# print(index)
# print(filtered_designs[index[0]]*1.07)# gives 1.034, 1.036

# print(f"for ulti load {designs_ult[np.argmin(designs_ult, axis=0)[0]]}")

design_1 = get_design(1)
part_split = design(design_1, 0.0014, [30, 7, 4], 17)
# print()
