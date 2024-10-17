import math


W_MTO = 259403.4 #kg
W_TO = W_MTO*0.8933 #kg
W_ZF = W_MTO*0.5886 #kg

rho_cruise = 0.316459 #kg/m3

D_fus = 5.727 #m
L_fus = 55.982 #m
b_wing = 62.72 #m
b_ref = 1.905 #m constant
l_wing_to_tail = 1 #m

S_wing = 363.54 #m2
S_htail = 41.5176 #m2
S_vtail = 17.70 #m2

V_cruise_TAS = 241.9574 #m/s
V_cruise_EAS = V_cruise_TAS * math.sqrt(rho_cruise/1.225) #m/s
V_dive = V_cruise_EAS*1.3125 #m/s

quarter_chord_sweep = math.radians(28.5) #degree
half_chord_sweep = math.radians(26.26) #degree
htail_sweep = math.radians(32.18) #degree
vtail_sweep = math.radians(35.2) #degree

thickness_ratio = 0.1 #constant

k_s = 0.447 #constant
k_w = 6.67*10**-3 #constant
k_h = 1 #constant
k_v = 1 #constant
k_f = 0.23 #constant

#structure_to_take_off_weight_ratio = k_s*math.sqrt(W_MTO)*((D_fus**2*L_fus)/W_TO)**0.24

b_s = b_wing/math.cos(half_chord_sweep)

W_wing = W_ZF*k_w*b_s**0.75*(1+math.sqrt(b_ref/b_s))*W_MTO**0.55*((b_s/thickness_ratio)/(W_ZF/S_wing))**0.3

W_htail = k_h*((S_htail**0.2*V_dive)/(math.cos(htail_sweep)))

W_vtail = k_v*((S_vtail**0.2*V_dive)/(math.cos(vtail_sweep)))

W_fus = k_f * math.sqrt(V_dive*(l_wing_to_tail)/(2*D_fus))*(math.pi*(D_fus/2)**2)**1.2

print(b_s, W_wing, W_htail, W_vtail)
