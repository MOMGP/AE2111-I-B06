import math


W_TO = 1 #kg
W_MTO = 1 #kg

D_fus = 1 #m
L_fus = 1 #m

k_s = 0.447 #constant
k_w = 6.67*10**-3 #constant


structure_to_take_off_ratio = k_s*math.sqrt(W_MTO)*((D_fus**2*L_fus)/W_TO)**0.24



