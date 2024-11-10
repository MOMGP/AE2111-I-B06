import math

#weights
W_MTO = 259403.4 * 2.205 #lb
W_ZF = W_MTO*0.5886 * 2.205 #lb
W_l=0.7 * W_MTO#lb
W_en= 9.630 * 2.205 #lb
W_apu =730.0 #could not find apu_uninstalled so I took this value from reference apu
W_uav=1200.0 #lb
W_c=49442 * 2.205 #lb

#fuselage lengths
D_fus = 5.727 * 3.2808399 #ft
L_fus = 55.982 * 3.2808399 #ft
b_wing = 62.72 * 3.2808399 #ft
l_wing_to_tail = 31.1233 * 3.2808399 #ft

#wing areas and spans
S_wing = 363.54 * 10.7639104 #ft2
S_control_surface = 109.76658 * 10.7639104 #ft2
S_htail = 86.482 * 10.7639104 #ft2
S_vtail = 69.45789 * 10.7639104 #ft2
S_f =819.22 *10.7639104#ft2
S_cs= S_control_surface#ft2

b_htail = 19.21 * 3.2808399#ft
b_vtail = 11.02 * 3.2808399#ft

#speeds
V_cruise_TAS = 241.9574 * 3.2808399 #ft/s
V_stall=(67.6)* 3.2808399 #ft/s

quarter_chord_sweep = math.radians(28.5) #degree
half_chord_sweep = math.radians(26.26) #degree
htail_sweep = math.radians(32.18) #degree
vtail_sweep = math.radians(35.2) #degree

#wing characteristics
aspect_ratio = b_wing**2/S_wing #constant
A_htail=(b_htail**2)/S_htail #unitless
A_vtail=(b_vtail**2)/S_vtail#unitless
thickness_ratio = 0.125 #constant
load_factor = 2.5 # constant
ultimate_load_factor = 1.5*load_factor #constant
taper_ratio = 0.3 #constant
tc_vtail=0.1 #constant

#landing gear
L_m=3.97*39.3701 #in
L_n =L_m #in

#engine
N_lt= 7*3.2808399#ft
N_w=4.17*3.2808399#ft
S_n=105.36*10.7639104 #ft2
L_ec=40*3.2808399#ft assumed wing position 20m in and 10 m tobith sides to each engine

#fuel
V_t=130.1158546*264.172 #gallons
V_i=0.33*V_t#gallons based on reference
V_p=0.67*V_t#gallons based on reference

#controls
N_f=6.0 #based on the raymer statement
N_m=2.0 #based on raymer
I_y= (((b_wing+ L_fus)/2)**2)* (W_MTO/4)*(0.46**2) #lb ft2
L_a= 2*L_ec #ft
V_pr = 100

#constants
K_uht =1.0 #constant
K_y= 0.3 * l_wing_to_tail #constant
H_tH_v =0.0 #for no Ttail
K_z = l_wing_to_tail #constant
K_door =1.0
K_lg = 1.0
K_ws=0.75*((1+2*taper_ratio)/(1+taper_ratio))*(b_wing*math.tan(quarter_chord_sweep)/L_fus)
K_mp=1.0
N_l= 3*1.5 #constant
N_mw=8.0 #constant
N_mss= 2.0 #assumed based on lg geometry
N_nw=2.0 #constant
K_np=1.0 #constant
K_ng=1.017 #constant
K_p= 1.0#constant
K_tr=1.18#constant
N_en = 2.0 #constant
N_t = 5.0 #constant based on reference(left, right, center, and in wings)
K_r=1.0#raymer
K_tp=1.0#raymer
N_c=8.0 #based on reference
R_kva=50.0#raymer
N_gen=N_en#raymer correlation
N_p= 280 + N_c #based on calcs by fuselage team


W_wing = 0.0051 * (W_MTO * load_factor)**0.557 * S_wing**0.649 * aspect_ratio**0.5 * thickness_ratio**-0.4 (1+S_control_surface)**0.1 * math.cos(quarter_chord_sweep)**-1.0




