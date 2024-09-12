import numpy as np
import matplotlib.pyplot as plt
from main import *
from Functions_for_calculating_shit import *
from Drag_polar_estimation_for_diff_configurations import *

#Plotting settings
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
cruise_req_thrust_over_weight=[]
climb_req_thrust_over_weight=[]
cs_far_119_req=[]
cs_far_121_a_req=[]
cs_far_121_b_req=[]
cs_far_121_c_req=[]
cs_far_121_d_req=[]
take_off_len_req=[]

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

Wing_loading=np.arange(100,12100,100)



#Approach speed requirement

approach_speed_req=1/2*rho_0/assumed_landing_mass_fraction*(assumed_aproach_speed/1.23)**2*assumed_max_CL


#Landing field length requirement

landing_field_req=1/assumed_landing_mass_fraction*(landing_field_len)/assumed_landing_field_coeff*rho_0*assumed_max_CL/2


#Cruise speed requirement

lapse_rate_req=thrust_lapse_calc(cruise_alt,M_cr)
p,T,rho=isa_calc(cruise_alt)
for i in range(0,len(Wing_loading)):
    req=(assumed_mass_fraction_cruise/lapse_rate_req*((CD_0*rho*(V_cruise)**2)/(assumed_mass_fraction_cruise*2*Wing_loading[i])+assumed_mass_fraction_cruise*Wing_loading[i]*2/(math.pi*Ar*e*rho*V_cruise**2)))
    cruise_req_thrust_over_weight.append(req)

#Climb rate requirement

p,T,rho=isa_calc(climb_alt_req)
a=math.sqrt(gamma*R_air*T)
for i in range (0,len(Wing_loading)):
    V_climb=math.sqrt(Wing_loading[i]*2/(rho*CL_highest_climb))
    M_climb=V_climb/a
    lapse_rate_req=thrust_lapse_calc(climb_alt_req,M_climb)
    req=assumed_mass_fraction_climb/lapse_rate_req*(math.sqrt(ROC**2*rho/(2*assumed_mass_fraction_climb*Wing_loading[i])*math.sqrt(CD_0*math.pi*Ar*e))+2*math.sqrt(CD_0/(math.pi*Ar*e)))
    climb_req_thrust_over_weight.append(req)


#Climb gradient requirements

#Climb gradient CS25.119
#segment=landing flaps=landing gear=down
CL_119=math.sqrt(cd_0_landing_landing_gear_on*math.pi*Ar*e_landing)
p,T,rho=isa_calc(0)
a=math.sqrt(gamma*R_air*T_0)
for i in range (0,len(Wing_loading)):
    V_grad= math.sqrt(Wing_loading[i] * 2 / (rho * CL_119))
    M_grad= V_grad / a
    lapse_rate_req=thrust_lapse_calc(0,M_grad)
    #We use one instead of other mass fractions because the requirements specify that for CS FAR 119 and CS FAR 121 (a-c) should be taken at max takeoff mass
    req=1/lapse_rate_req*(c_119+2*math.sqrt(cd_0_landing_landing_gear_on/(math.pi*Ar*e_landing)))
    cs_far_119_req.append(req)

#Climb gradient CS 25.121(a)
CL_121_a=math.sqrt(cd_0_take_off_landing_gear_on*math.pi*Ar*e_take_off)
p,T,rho=isa_calc(0)
a=math.sqrt(gamma*R_air*T_0)
for i in range (0,len(Wing_loading)):
    V_grad=math.sqrt(Wing_loading[i] * 2 / (rho * CL_121_a))
    M_grad= V_grad / a
    lapse_rate_req=thrust_lapse_calc(0,M_grad)
    req=2*1/lapse_rate_req*(c_121_a+2*math.sqrt((cd_0_take_off_landing_gear_on)/(math.pi*Ar*e_take_off)))
    cs_far_121_a_req.append(req)

#Climb gradient CS 25.121(b)
CL_121_b=math.sqrt(cd_0_take_off_landing_gear_off*math.pi*Ar*e_take_off)
p,T,rho=isa_calc(0)
a=math.sqrt(gamma*R_air*T_0)
for i in range (0,len(Wing_loading)):
    V_grad=math.sqrt(Wing_loading[i] * 2 / (rho * CL_121_b))
    M_grad= V_grad / a
    lapse_rate_req=thrust_lapse_calc(0,M_grad)
    req=2*1/lapse_rate_req*(c_121_b+2*math.sqrt((cd_0_take_off_landing_gear_off)/(math.pi*Ar*e_take_off)))
    cs_far_121_b_req.append(req)

#Climb gradient CS 25.121(c)
CL_121_c=math.sqrt(CD_0*math.pi*Ar*e)
p,T,rho=isa_calc(0)
a=math.sqrt(gamma*R_air*T_0)
for i in range (0,len(Wing_loading)):
    V_grad=math.sqrt(Wing_loading[i] * 2 / (rho * CL_121_c))
    M_grad= V_grad / a
    lapse_rate_req=thrust_lapse_calc(0,M_grad)
    req=2*1/lapse_rate_req*(c_121_c+2*math.sqrt((CD_0)/(math.pi*Ar*e)))
    cs_far_121_c_req.append(req)

#Climb gradient CS 25.121(d)
CL_121_d=math.sqrt(cd_0_landing_landing_gear_off*math.pi*Ar*e_landing)
p,T,rho=isa_calc(0)
a=math.sqrt(gamma*R_air*T_0)
for i in range (0,len(Wing_loading)):
    V_grad=math.sqrt(Wing_loading[i] * 2 / (rho * CL_121_d))
    M_grad= V_grad / a
    lapse_rate_req=thrust_lapse_calc(0,M_grad)
    req=2*assumed_landing_mass_fraction/lapse_rate_req*(c_121_d+2*math.sqrt((cd_0_landing_landing_gear_off)/(math.pi*Ar*e_landing)))
    cs_far_121_d_req.append(req)

#Landing field length req
for i in range (0,len(Wing_loading)):
    req=1.15*math.sqrt(Wing_loading[i]/(take_off_field_len*kt*rho_0*g_0*math.pi*Ar*e_take_off))+4*obstacle_h/take_off_field_len
    take_off_len_req.append(req)




#Finding the design point
bound_right=float(min(approach_speed_req,landing_field_req))
print(bound_right)
bound_right=int(bound_right/100)*100
pos=int(bound_right/100)-1
bound_low=max(cruise_req_thrust_over_weight[pos],climb_req_thrust_over_weight[pos],cs_far_119_req[pos],cs_far_121_a_req[pos],cs_far_121_b_req[pos],cs_far_121_c_req[pos],cs_far_121_d_req[pos],take_off_len_req[pos])
print(bound_low)
bound_low=round(bound_low,3)
print(bound_low,bound_right)
#Ploting the req curves
plt.figure(figsize=(8,6))
plt.axvline(x=approach_speed_req,color='r',label='Approach Speed')
plt.axvline(x=landing_field_req,color='g',label='Landing Field')
plt.plot(Wing_loading,cruise_req_thrust_over_weight,color='b',label='Cruise Speed')
plt.plot(Wing_loading,climb_req_thrust_over_weight,color='y',label='Climb Rate')
plt.plot(Wing_loading,cs_far_119_req,color='purple',label='CS Far 119')
plt.plot(Wing_loading,cs_far_121_a_req,color='black',label='CS Far 121(a)')
plt.plot(Wing_loading,cs_far_121_b_req,color='orange',label='CS Far 121(b)')
plt.plot(Wing_loading,cs_far_121_c_req,color='brown',label='CS Far 121(c)')
plt.plot(Wing_loading,cs_far_121_d_req,color='pink',label='CS Far 121(d)')
plt.plot(Wing_loading,take_off_len_req,color='grey',label='Take off length')
plt.plot(bound_right,bound_low,'o')
plt.plot(Wing_loading_ref,TW_ref_aircraft,'o',label='Reference aircrafts',color='black')
plt.xlim([0,10000])
plt.ylim([0,0.5])
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),
          ncol=3, fancybox=True, shadow=True)
plt.xlabel('W/S-[N/m2]')
plt.ylabel('T/W=[N/N]')
plt.show()

print('Take off thrust is equal to',bound_low*assumed_maximum_takeoff_mass*9.81)
print('Wing area is equal to',assumed_maximum_takeoff_mass*9.81/bound_right)

