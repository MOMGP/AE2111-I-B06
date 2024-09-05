import numpy as np
import matplotlib.pyplot as plt
from main import *
from Functions_for_calculating_shit import*
#Plotting settings
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
cruise_req_thrust_over_weight=[]

#start of the matching diagram

#Wing loading range

Wing_loading=np.arange(100,10100,100)
#Approach speed requirement
approach_speed_req=1/2*rho_0/assumed_landing_mass_fraction*(assumed_aproach_speed/1.23)**2*assumed_max_CL

#Landing field length requirement
landing_field_req=1/assumed_landing_mass_fraction*(landing_field_len)/assumed_landing_field_coeff*rho_0*assumed_max_CL/2
print(approach_speed_req,landing_field_req)
#Cruise speed requirement
lapse_rate_req=thrust_lapse_calc(cruise_alt,M_cr)
p,T,rho=isa_calc(cruise_alt)
for i in range(0,len(Wing_loading)):
    req=(assumed_mass_fraction_cruise/lapse_rate_req*((CD_0*rho*(V_cruise)**2)/(assumed_mass_fraction_cruise*2*Wing_loading[i])+assumed_mass_fraction_cruise*Wing_loading[i]*2/(math.pi*Ar*e*rho*V_cruise**2)))
    cruise_req_thrust_over_weight.append(req)


plt.figure()
plt.axvline(x=approach_speed_req,color='r',label='Approach Speed Requirement')
plt.axvline(x=landing_field_req,color='g',label='Landing Field Requirement')
plt.plot(Wing_loading,cruise_req_thrust_over_weight,color='b',label='Cruise speed requirement')
plt.xlim([0,10000])
plt.ylim([0,0.7])
plt.show()

