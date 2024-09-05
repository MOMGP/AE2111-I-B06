import numpy as np
import matplotlib.pyplot as plt
from main import *


#start of the matching diagram

#Approach speed requirement
approach_speed_req=1/2*rho_0/assumed_landing_mass_fraction*(assumed_aproach_speed/1.23)**2*assumed_max_CL

#Landing field length requirement
landing_field_req=1/assumed_landing_mass_fraction*(landing_field_len)/assumed_landing_field_coeff*rho_0*assumed_max_CL/2






plt.figure()
plt.axvline(x=approach_speed_req,color='r',label='Approach Speed Requirement')
plt.axvline(x=landing_field_req,color='g',label='Landing Field Requirement')
plt.xlim([0,10000])
plt.show()

