import numpy as np
import matplotlib.pyplot as plt
#from main import *
from Functions_for_isa_and_alphat import *
#from Drag_polar_estimation_for_diff_configurations import *

#Plotting settings
font_size_full_scr=14
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

#Reference aircraft data points
Wing_loading_ref=[6565,
6565,
5930,
6627,
6471,
6822,
7470,
6943,]
TW_ref_aircraft=[0.261160395,
0.265372659,
0.250480946,
0.256142723,
0.298883713,
0.285765682,
0.223044353,
0.220532531]



#start of the matching diagram

#Wing loading range

Wing_loading=np.arange(100,10100,100)
assumed_landing_mass_fraction=0.7
assumed_max_CL=2.5#when landing
assumed_mass_fraction_cruise=0.93
Ar=10.2
e=0.78#Oswald efficiency factor as from REQ-WNG-02
CD_0=0.0156#during cruise as from REQ-WNG-03
assumed_CL_to=2.1
assumed_maximum_takeoff_mass=259033 #kg


#Approach speed requirement

#Approach speed requirement
def approach_speed_req(landing_mass_fraction, CL_max_landing):
    approach_speed = 70 # approach speed
    rho_0=1.225
    return 1/2*rho_0/landing_mass_fraction*(approach_speed/1.23)**2*CL_max_landing
approach_speed_req_val = approach_speed_req(assumed_landing_mass_fraction, assumed_max_CL)


#Landing field length requirement
def landing_field_req(landing_mass_fraction, CL_max_landing):
    landing_field_len=1865 #m
    rho_0=1.225
    assumed_landing_field_coeff=0.45 #Fromn adsee book, eq 7.9, for CS/FAR-25
    return 1/landing_mass_fraction*(landing_field_len)/assumed_landing_field_coeff*rho_0*CL_max_landing/2
landing_field_req_val =landing_field_req(assumed_landing_mass_fraction, assumed_max_CL)

#Cruise speed requirement
def cruise_req_thrust_over_weight(AR, e, CD_0):
    mass_fraction_cruise = 0.93
    cruise_alt = 11887.2
    M_cr=0.82
    V_cruise = 241.956
    Wing_loading = np.arange(100,10100,100)
    cruise_req=[]
    lapse_rate_req=thrust_lapse_calc(cruise_alt,M_cr)
    p,T,rho=isa_calc(cruise_alt)
    for i in range(0,len(Wing_loading)):
        req=(mass_fraction_cruise/lapse_rate_req*((CD_0*rho*(V_cruise)**2)/(mass_fraction_cruise*2*Wing_loading[i])+mass_fraction_cruise*Wing_loading[i]*2/(math.pi*AR*e*rho*V_cruise**2)))
        cruise_req.append(req)
    return cruise_req
cruise_req_thrust_over_weight_val = cruise_req_thrust_over_weight(Ar, e, CD_0)

#Climb rate requirement
def climb_req_thrust_over_weight(e, AR,  CD_0):
    climb_alt_req=12000
    Wing_loading=np.arange(100,10100,100) 
    mass_fraction_climb=0.93
    p,T,rho=isa_calc(climb_alt_req)
    gamma=1.4
    R_air = 287.05
    ROC=0.5#m/s requirement
    CL_highest_climb=math.sqrt(CD_0*math.pi*AR*e)
    T=T+15 #Hot conditions
    a=math.sqrt(gamma*R_air*T)
    climb_req=[]
    for i in range (0,len(Wing_loading)):
        V_climb=math.sqrt(Wing_loading[i]*2/(rho*CL_highest_climb))
        M_climb=V_climb/a
        lapse_rate_req=thrust_lapse_calc(climb_alt_req,M_climb)
        req=mass_fraction_climb/lapse_rate_req*(math.sqrt(ROC**2*rho/(2*mass_fraction_climb*Wing_loading[i])*CL_highest_climb)+2*math.sqrt(CD_0/(math.pi*AR*e)))
        climb_req.append(req)
    return climb_req
climb_req_thrust_over_weight_val= climb_req_thrust_over_weight(e, Ar, CD_0)

#Climb gradient requirements

#Climb gradient CS25.119
#segment=landing flaps=landing gear=down
def climb_gradient_CS25119(AR, e, CD_0 = 0.0156):
    Wing_loading=np.arange(100,10100,100)
    gamma=1.4
    R_air = 287.05
    flap_deflection_angle_landing=40
    cd0_landing_gear_diff=0.0100
    cd0_diff_landing=flap_deflection_angle_landing*0.0013
    cd_0_landing_landing_gear_on=CD_0+cd0_diff_landing+cd0_landing_gear_diff
    oswald_diff_landing=flap_deflection_angle_landing*0.0026
    e_landing=e+oswald_diff_landing
    CL_119=math.sqrt(cd_0_landing_landing_gear_on*math.pi*AR*e_landing)
    cs_far_119_req_list = []
    p,T,rho=isa_calc(0)
    c_119=3.2/100
    a=math.sqrt(gamma*R_air*T)
    for i in range (0,len(Wing_loading)):
        V_grad= math.sqrt(Wing_loading[i] * 2 / (rho * CL_119))
        M_grad= V_grad / a
        lapse_rate_req=thrust_lapse_calc(0,M_grad)
        #We use one instead of other mass fractions because the requirements specify that for CS FAR 119 and CS FAR 121 (a-c) should be taken at max takeoff mass
        req=1/lapse_rate_req*(c_119+2*math.sqrt(cd_0_landing_landing_gear_on/(math.pi*AR*e_landing)))
        cs_far_119_req_list.append(req)
    return cs_far_119_req_list
cs_far_119_req_val = climb_gradient_CS25119(Ar, e, CD_0 = 0.0156)

#Climb gradient CS 25.121(a)
def climb_gradient_CS25121a(AR, e, CD_0):
    Wing_loading=np.arange(100,10100,100)
    R_air = 287.05
    gamma=1.4
    flap_deflection_angle_take_off=20
    oswald_diff_take_off=flap_deflection_angle_take_off*0.0026
    cd0_diff_take_off=flap_deflection_angle_take_off*0.0013
    cd0_landing_gear_diff=0.0100
    e_take_off=e+oswald_diff_take_off
    cd_0_take_off_landing_gear_on=CD_0+cd0_diff_take_off+cd0_landing_gear_diff
    CL_121_a=math.sqrt(cd_0_take_off_landing_gear_on*math.pi*AR*e_take_off)
    p,T,rho=isa_calc(0)
    a=math.sqrt(gamma*R_air*T)
    cs_far_121_a_req_list=[]
    c_121_a=0
    for i in range (0,len(Wing_loading)):
        V_grad=math.sqrt(Wing_loading[i] * 2 / (rho * CL_121_a))
        M_grad= V_grad / a
        lapse_rate_req=thrust_lapse_calc(0,M_grad)
        req=2*1/lapse_rate_req*(c_121_a+2*math.sqrt((cd_0_take_off_landing_gear_on)/(math.pi*AR*e_take_off)))
        cs_far_121_a_req_list.append(req)
    return cs_far_121_a_req_list
cs_far_121_a_req_val =climb_gradient_CS25121a(Ar,e, CD_0)

#Climb gradient CS 25.121(b)
def climb_gradient_CS25121b(AR, e, CD_0):
    Wing_loading=np.arange(100,10100,100)
    R_air = 287.05
    gamma=1.4
    flap_deflection_angle_take_off=20
    oswald_diff_take_off=flap_deflection_angle_take_off*0.0026
    cd0_diff_take_off=flap_deflection_angle_take_off*0.0013
    e_take_off=e+oswald_diff_take_off
    cd_0_take_off_landing_gear_off=CD_0+cd0_diff_take_off
    CL_121_b=math.sqrt(cd_0_take_off_landing_gear_off*math.pi*AR*e_take_off)
    p,T,rho=isa_calc(0)
    a=math.sqrt(gamma*R_air*T)
    cs_far_121_b_req_list=[]
    c_121_b=2.4/100
    for i in range (0,len(Wing_loading)):
        V_grad=math.sqrt(Wing_loading[i] * 2 / (rho * CL_121_b))
        M_grad= V_grad / a
        lapse_rate_req=thrust_lapse_calc(0,M_grad)
        req=2*1/lapse_rate_req*(c_121_b+2*math.sqrt((cd_0_take_off_landing_gear_off)/(math.pi*AR*e_take_off)))
        cs_far_121_b_req_list.append(req)
    return cs_far_121_b_req_list
cs_far_121_b_req_val =climb_gradient_CS25121b(Ar, e, CD_0)

#Climb gradient CS 25.121(c)
def climb_gradient_CS25121c(AR, e, CD_0):
    Wing_loading=np.arange(100,10100,100)
    R_air = 287.05
    gamma=1.4
    CL_121_c=math.sqrt(CD_0*math.pi*AR*e)
    p,T,rho=isa_calc(0)
    a=math.sqrt(gamma*R_air*T)
    cs_far_121_c_req_list = []
    c_121_c=1.2/100
    for i in range (0,len(Wing_loading)):
        V_grad=math.sqrt(Wing_loading[i] * 2 / (rho * CL_121_c))
        M_grad= V_grad / a
        lapse_rate_req=thrust_lapse_calc(0,M_grad)
        req=2*1/lapse_rate_req*(c_121_c+2*math.sqrt((CD_0)/(math.pi*AR*e)))
        cs_far_121_c_req_list.append(req)
    return cs_far_121_c_req_list
cs_far_121_c_req_val =climb_gradient_CS25121c(Ar, e, CD_0)

#Climb gradient CS 25.121(d)
def climb_gradient_CS25121d(AR, e, CD_0, landing_mass_fraction):
    flap_deflection_angle_landing=40
    oswald_diff_landing=flap_deflection_angle_landing*0.0026
    cd0_diff_landing=flap_deflection_angle_landing*0.0013
    e_landing=e+oswald_diff_landing
    cd_0_landing_landing_gear_off=CD_0+cd0_diff_landing
    Wing_loading=np.arange(100,10100,100)
    R_air = 287.05
    gamma=1.4
    CL_121_d=math.sqrt(cd_0_landing_landing_gear_off*math.pi*AR*e_landing)
    p,T,rho=isa_calc(0)
    a=math.sqrt(gamma*R_air*T)
    c_121_d=2.1/100
    cs_far_121_d_req_list=[]
    for i in range (0,len(Wing_loading)):
        V_grad=math.sqrt(Wing_loading[i] * 2 / (rho * CL_121_d))
        M_grad= V_grad / a
        lapse_rate_req=thrust_lapse_calc(0,M_grad)
        req=2*landing_mass_fraction/lapse_rate_req*(c_121_d+2*math.sqrt((cd_0_landing_landing_gear_off)/(math.pi*AR*e_landing)))
        cs_far_121_d_req_list.append(req)
    return cs_far_121_d_req_list
cs_far_121_d_req_val= climb_gradient_CS25121d(Ar, e, CD_0, assumed_landing_mass_fraction)

#Landing field length req
def landing_field_length_req(AR, e, CL_to):
    flap_deflection_angle_take_off=20
    oswald_diff_take_off=flap_deflection_angle_take_off*0.0026
    e_take_off=e+oswald_diff_take_off
    kt=0.85 #for jet airplanes
    Wing_loading=np.arange(100,10100,100)
    R_air = 287.05
    gamma=1.4
    g_0=9.80665
    p,T,rho=isa_calc(0)
    take_off_len_req=[]
    obstacle_h=11#m, assumed
    take_off_field_len=2790#m
    a=math.sqrt(gamma*R_air*T)
    for i in range (0,len(Wing_loading)):
        V_lf=math.sqrt(Wing_loading[i] * 2 / (rho * CL_to))
        M_lf= V_lf / a
        lapse_rate_req=thrust_lapse_calc(0,M_lf)
        req=1.15*math.sqrt(Wing_loading[i]/(take_off_field_len*kt*rho*g_0*math.pi*AR*e_take_off))+4*obstacle_h/take_off_field_len
        take_off_len_req.append(req)
    return take_off_len_req
take_off_len_req_val =landing_field_length_req(Ar, e, assumed_CL_to)

#Finding the design point
def returns_for_optimization(approach_speed_req_val, landing_field_req_val, cruise_req_thrust_over_weight_val, climb_req_thrust_over_weight_val, cs_far_119_req_val, cs_far_121_b_req_val, cs_far_121_c_req_val, cs_far_121_d_req_val, take_off_len_req_val):
    bound_right=float(min(approach_speed_req_val,landing_field_req_val))
    bound_right=int(bound_right/100)*100
    pos=int(bound_right/100)-1
    bound_low=max(cruise_req_thrust_over_weight_val[pos],climb_req_thrust_over_weight_val[pos],cs_far_119_req_val[pos],cs_far_121_a_req_val[pos],cs_far_121_b_req_val[pos],cs_far_121_c_req_val[pos],cs_far_121_d_req_val[pos],take_off_len_req_val[pos])
    bound_low=round(bound_low,3)
    return (bound_low,bound_right)
bound_low,bound_right= returns_for_optimization(approach_speed_req_val, landing_field_req_val, cruise_req_thrust_over_weight_val, climb_req_thrust_over_weight_val, cs_far_119_req_val, cs_far_121_b_req_val, cs_far_121_c_req_val, cs_far_121_d_req_val, take_off_len_req_val)
#Ploting the req curves
import matplotlib.pyplot as plt

# Create a figure with a specific size (width, height) in inches
fig, ax = plt.subplots(figsize=(12, 6))  # Adjust the size as needed
font_size_full_scr=19
# Plotting data
ax.axvline(x=approach_speed_req_val, color='r', label='Approach Speed')
ax.axvline(x=landing_field_req_val, color='g', label='Landing Field')
ax.plot(Wing_loading, cruise_req_thrust_over_weight_val, color='b', label='Cruise Speed')
ax.plot(Wing_loading, climb_req_thrust_over_weight_val, color='y', label='Climb Rate')
ax.plot(Wing_loading, cs_far_119_req_val, color='purple', label='CS Far 119')
ax.plot(Wing_loading, cs_far_121_a_req_val, color='black', label='CS Far 121(a)')
ax.plot(Wing_loading, cs_far_121_b_req_val, color='orange', label='CS Far 121(b)')
ax.plot(Wing_loading, cs_far_121_c_req_val, color='brown', label='CS Far 121(c)')
ax.plot(Wing_loading, cs_far_121_d_req_val, color='pink', label='CS Far 121(d)')
ax.plot(Wing_loading, take_off_len_req_val, color='grey', label='Take off length')
ax.plot(Wing_loading_ref, TW_ref_aircraft, 'o', label='Reference aircrafts', color='black')
ax.plot(bound_right, bound_low, 'o', label="Design point")

# Set x and y axis limits
ax.set_xlim([0, 10000])
ax.set_ylim([0, 0.5])

# Set x and y ticks font size
ax.tick_params(axis='x', labelsize=font_size_full_scr)
ax.tick_params(axis='y', labelsize=font_size_full_scr)

# Set legend with adjusted position and number of columns
ax.legend(loc="upper left", bbox_to_anchor=(1, 1.05), ncol=1, fancybox=True, shadow=True, fontsize=font_size_full_scr)

# Set axis labels
ax.set_xlabel(r"W/S $\left[{N}/{m^2}\right]$", fontsize=font_size_full_scr)
ax.set_ylabel(r"T/W $\left[{N}/{N}\right]$", fontsize=font_size_full_scr)

# Adjust layout to ensure everything fits well
plt.subplots_adjust(right=0.65, top=0.9, bottom=0.16)  # Adjust this value to reduce right space

# Save the figure
plt.savefig("TW-WS Diagram.pdf", format='pdf')

# Close the plot
plt.close()

# Print additional information
print('Take off thrust is equal to', bound_low * assumed_maximum_takeoff_mass * 9.81)
print('Wing area is equal to', assumed_maximum_takeoff_mass * 9.81 / bound_right)