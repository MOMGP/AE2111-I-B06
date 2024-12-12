# from main_iteration import deflection, angle_of_rotation
# from geometry_WP4_2 import moments_of_inertia
import numpy as np

bending = []
with open("Work package 4\Bending_critical_loadcases\loadcase21.txt", 'r') as f:
    count = 0
    for i in f:
        line = i.strip('\n')
        count += 1
        if count % 50 == 0 or count == 1:
            bending.append(line)
print(bending, len(bending))