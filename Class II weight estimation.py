import math
from Matching_diagram_Andrei import *
W_fin.append(259403.4*0.5886-27669)
#weights
W_MTO = 259403.4 * 2.205 #lb
W_ZF = W_MTO*0.5886 #lb
W_l=0.7 * W_MTO#lb
W_en= 9630 * 2.205 #lb
W_apu = 730.0 #could not find apu_uninstalled, so I took this value from reference apu
W_uav=1200.0 #lb
W_c=27669 * 2.205 #lb (changed this number because in our scenario we are not at max payload)
print("W_en",W_en)
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
S_fus =819.22 *10.7639104#ft2
S_cs= S_control_surface#ft2

S_elevator = 0.33*S_htail #ft2

b_htail = 19.21 * 3.2808399#ft
b_vtail = 11.02 * 3.2808399#ft

#speeds
V_cruise_TAS = 241.9574 * 3.2808399 #ft/s
V_stall= 67.6* 3.2808399 #ft/s

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
L_over_D = 19.65 #unitless

#landing gear
L_m=3.97*39.3701 #in
L_n =L_m #in

#engine
N_lt= 7*3.2808399#ft
N_w=4.17*3.2808399#ft
S_n=105.36*10.7639104 #ft2
L_ec=40*3.2808399#ft assumed wing position 20m in and 10 m to both sides to each engine

#fuel
V_t=130.1158546*264.172 #gallons
V_i=0.33*V_t#gallons based on reference
V_p=0.67*V_t#gallons based on reference

#controls
N_f=6.0 #based on the Raymer statement
N_m=2.0 #based on Raymer
I_y= (((b_wing+ L_fus)/2)**2)* (W_MTO/4)*(0.46**2) #lb ft2
L_a= 2*L_ec #ft
V_pr = 920.823*35.315 #ft^2

#constants
K_uht =1.0 #constant
K_y= 0.3 * l_wing_to_tail #constant
H_tH_v =0.0 #for no T-tail
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


W_wing = 0.0051 * (W_MTO * ultimate_load_factor)**0.557 * S_wing**0.649 * aspect_ratio**0.5 * thickness_ratio**-0.4 * (1+taper_ratio)**0.1 * math.cos(quarter_chord_sweep)**-1.0 * S_control_surface**0.1
print("W_wing",W_wing)
W_htail = 0.0379 * K_uht * (1+ D_fus / b_wing)**-0.25 * W_MTO**0.639 * ultimate_load_factor**0.10 * S_htail**0.75 * l_wing_to_tail**-1.0 * K_y**0.704 * math.cos(htail_sweep)**-1.0 * A_htail**0.166 * (1+S_elevator/S_htail)**0.1
print("W_htail",W_htail)
W_vtail = 0.0026 * (1 + H_tH_v)**0.225 * W_MTO**0.556 * ultimate_load_factor**0.536 * l_wing_to_tail**-0.5 * S_vtail**0.5 * K_z**0.875 * math.cos(vtail_sweep)**-1.0 * A_vtail**0.35 * tc_vtail**-0.5
print("W_vtail",W_vtail)
W_fuselage = 0.3280 * K_door * K_lg * (W_MTO * ultimate_load_factor)**0.5 * L_fus**0.25 * S_fus**0.302 * (1+K_ws)**0.04 * L_over_D**0.10
print("W_fuselage",W_fuselage)
W_mlg = 0.0106 * K_mp * W_l**0.888 * N_l**0.25 * L_m**0.4 * N_mw**0.321 * N_mss**-0.5 ** V_stall**0.1
print("W_mlg",W_mlg)
W_nlg = 0.032 * K_np * W_l**0.646 * N_l**0.25 * L_n**0.5 **N_nw**0.45
print("W_nlg",W_nlg)
W_ec = 2.331 * W_en**0.901 * K_p * K_tr #prelim for equation below
print("W_ec",W_ec)
W_nacelle = 0.6724 * K_ng * N_lt**0.1 * N_w**0.294 * ultimate_load_factor**0.119 * W_ec**0.611 * N_en**0.984 * S_n**0.224
print("W_nacelle",W_nacelle)
W_engine_controls = 5*N_en + 0.8*L_ec
print("W_engine_controls",W_engine_controls)
W_starter = 49.19 * (N_en * W_en/1000)**0.541
print("W_starter",W_starter)
W_fuel_system = 2.405 * V_t**0.606 * (1+V_i/V_t)**-1.0 * (1+V_p/V_t)*N_t**0.5
print("W_fuel_system",W_fuel_system)
W_flight_controls = 145.9 * N_f**0.554 * (1+N_m/N_f)**-1.0 * S_cs**0.20 * (I_y * 10**-6)**0.07
print("W_flight_controls",W_flight_controls)
W_apu = 2.2 * W_apu
print("W_apu",W_apu)
W_instruments = 4.509 * K_r * K_tp * N_c**0.541 * N_en * (L_fus + b_wing)**0.5
print("W_instruments",W_instruments)
W_hydraulics = 0.2673 * N_f * (L_fus + b_wing)**0.937
print("W_hydraulics",W_hydraulics)
W_electrical = 7.291 * R_kva**0.782 * L_a **0.346 * N_gen**0.10
print("W_electrical",W_electrical)
W_avionics = 1.73 * W_uav**0.983
print("W_avionics",W_avionics)
W_furnishings = 0.0577 * N_c**0.1 * W_c**0.393 * S_fus**0.75
print("W_furnishings",W_furnishings)
W_air_conditioning = 62.36 * N_p**0.25 * (V_pr/1000)**0.604 * W_uav**0.10
print("W_air_conditioning",W_air_conditioning)
W_anti_ice = 0.002* W_MTO
print("W_anti_ice",W_anti_ice)
W_handling_gear = 3.0 * 10**-4 * W_MTO
print("W_handling_gear",W_handling_gear)
#W_military_cargo_handling_system ?

print("old weight: ", (W_ZF-W_c)/2.205, " in kg")
W_OEW = W_wing + W_htail + W_vtail + W_fuselage + W_mlg + W_nlg + W_nacelle + W_engine_controls + W_starter + W_fuel_system + W_flight_controls + W_apu + W_instruments + W_hydraulics + W_electrical + W_avionics + W_furnishings + W_air_conditioning + W_anti_ice + W_handling_gear
print("new weight: ", W_OEW/2.205, " in kg")
