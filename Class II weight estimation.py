import math

W_MTO = 259403.4 * 2.205 #lb
W_ZF = W_MTO*0.5886 * 2.205 #lb

D_fus = 5.727 * 3.2808399 #ft
L_fus = 55.982 * 3.2808399 #ft
b_wing = 62.72 * 3.2808399 #ft
l_wing_to_tail = 31.1233 * 3.2808399 #ft

S_wing = 363.54 * 10.7639104 #ft2
S_control_surface = 109.76658 * 10.7639104 #ft2
S_htail = 86.482 * 10.7639104 #ft2
S_vtail = 69.45789 * 10.7639104 #ft2

V_cruise_TAS = 241.9574 * 3.2808399 #ft/s

quarter_chord_sweep = math.radians(28.5) #degree
half_chord_sweep = math.radians(26.26) #degree
htail_sweep = math.radians(32.18) #degree
vtail_sweep = math.radians(35.2) #degree

aspect_ratio = b_wing**2/S_wing #constant
thickness_ratio = 0.125 #constant
load_factor = 2.5 # constant
ultimate_load_factor = 1.5*load_factor #constant
taper_ratio = 0.3 #constant

W_wing = 0.0051 * (W_MTO * load_factor)**0.557 * S_wing**0.649 * aspect_ratio**0.5 * thickness_ratio**-0.4 (1+)**0.1 * math.cos(quarter_chord_sweep)**-1*




