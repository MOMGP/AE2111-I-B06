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
#start of the matching diagram

#Wing loading range

Wing_loading=np.arange(100,10100,100)



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
print(CL_119,cd_0_landing_landing_gear_on*2)
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
print(CL_121_a,cd_0_take_off_landing_gear_on*2)
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
print(CL_121_b,cd_0_take_off_landing_gear_off*2)
a=math.sqrt(gamma*R_air*T_0)
for i in range (0,len(Wing_loading)):
    V_grad=math.sqrt(Wing_loading[i] * 2 / (rho * CL_121_b))
    M_grad= V_grad / a
    lapse_rate_req=thrust_lapse_calc(0,M_grad)
    req=2*1/lapse_rate_req*(c_121_b+2*math.sqrt((cd_0_take_off_landing_gear_off)/(math.pi*Ar*e_take_off)))
    cs_far_121_b_req.append(req)

#Climb gradient CS 25.121(c)
CL_121_c=math.sqrt(CD_0*math.pi*Ar*e)
print(CL_121_c,CD_0*2)
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
print(CL_121_d,cd_0_landing_landing_gear_off*2)
a=math.sqrt(gamma*R_air*T_0)
for i in range (0,len(Wing_loading)):
    V_grad=math.sqrt(Wing_loading[i] * 2 / (rho * CL_121_d))
    M_grad= V_grad / a
    lapse_rate_req=thrust_lapse_calc(0,M_grad)
    req=2*assumed_landing_mass_fraction/lapse_rate_req*(c_121_d+2*math.sqrt((cd_0_landing_landing_gear_off)/(math.pi*Ar*e_landing)))
    cs_far_121_d_req.append(req)
    print(req)

#Landing field length req

#Ploting the req curves


plt.figure()
plt.axvline(x=approach_speed_req,color='r',label='Approach Speed Requirement')
plt.axvline(x=landing_field_req,color='g',label='Landing Field Requirement')
plt.plot(Wing_loading,cruise_req_thrust_over_weight,color='b',label='Cruise speed requirement')
plt.plot(Wing_loading,climb_req_thrust_over_weight,color='y',label='Climb rate requirement')
plt.plot(Wing_loading,cs_far_119_req,color='purple',label='CS Far 119 requirement')
plt.plot(Wing_loading,cs_far_121_a_req,color='black',label='CS Far 121 a_requirement')
plt.plot(Wing_loading,cs_far_121_b_req,color='orange',label='CS Far 121 b_requirement')
plt.plot(Wing_loading,cs_far_121_c_req,color='brown',label='CS Far 121 c_requirement')
plt.plot(Wing_loading,cs_far_121_d_req,color='pink',label='CS Far 121 d_requirement')
plt.xlim([0,10000])
plt.ylim([0,0.7])
plt.show()

