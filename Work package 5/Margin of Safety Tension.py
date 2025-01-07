from Geometry import get_design_final, I_xx_global, scaled_chord
from Divide_Geom import design
from get_points import get_points
import numpy as np
b = 53.57

design_1 = get_design_final(1)
part_split_1 = design(design_1[0], design_1[1], design_1[2], design_1[3])
I_xx_1 = I_xx_global(part_split_1)

design_2 = get_design_final(2)
part_split_2 = design(design_2[0], design_2[1], design_2[2], design_2[3])
I_xx_2 = I_xx_global(part_split_2)

design_3 = get_design_final(3)
part_split_3 = design(design_3[0], design_3[1], design_3[2], design_3[3])
I_xx_3 = I_xx_global(part_split_3)

failure_stress = 349.2*(10**6) #Pa

def y_max(y):
    spars_coordinates = np.array(get_points(0.2, 0.4, 0.7, 1)) * scaled_chord(y)
    y_max = spars_coordinates[2,1]
    return y_max


'''
spanwise_stress = []

moment = moment_at_full_position(CL_d,rho,V,n)


for i in range(0,26.78,0.01):
    stress = (moment[i] * y_max[i]) / I_xx[i]
    spanwise_stress.append(stress)

MOS = fail_stress / spanwise_stress
'''

#y_max at each location