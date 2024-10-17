import math


W_TO = 1 #kg
W_MTO = 1 #kg
W_ZF = 1 #kg

D_fus = 1 #m
L_fus = 1 #m
b_wing = 1 #m
b_ref = 1.905 #m

S_wing = 1 #m2
S_htail = 1 #m2
S_vtail = 1 #m2

V_d = 1 #m/s cruise speed in terms of EAS

quarter_chord_sweep = 1 #degree
half_chord_sweep = 1 #degree
htail_sweep = 1 #degree
vtail_sweep = 1 #degree

thickness_ratio = 1 #constant

k_s = 0.447 #constant
k_w = 6.67*10**-3 #constant
k_h = 1 #constant
k_v = 1 #constant

structure_to_take_off_weight_ratio = k_s*math.sqrt(W_MTO)*((D_fus**2*L_fus)/W_TO)**0.24

b_s = b_wing/math.cos(half_chord_sweep)

W_wing = W_ZF*k_w*b_s**0.75*(1+math.sqrt(b_ref/b_s))*W_MTO**0.55*((b_s/thickness_ratio)/(W_ZF/S_wing))**0.3

W_htail = k_h*((S_htail**0.2*V_d)/(math.cos(htail_sweep)))

W_vtail = k_v*((S_vtail**0.2*V_d)/(math.cos(vtail_sweep)))




