import math
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy.integrate import quad, dblquad
from Aero_loading_XFLR5 import Lift_for_integrating
from Aero_loading_XFLR5 import lift_distribution_any_CL
from Aero_loading_XFLR5 import chord_length_interpolation
from Aero_loading_XFLR5 import normal_force_for_integrating
from Aero_loading_XFLR5 import Lift_distribution_for_any_load_case

#--------------------------------Initial values--------------------------------#

n = 1
Pe = 9.81*9630*n        # Engine weight [N]
ye = 9.37             # Engine spanwise location [m]
Me = Pe*ye            # Engine bending moment around the root [Nm]
b = 53.57             # Wing span [m]
CL_d = 3
rho = 1.225
V = 62
#-------------------------------------------------------------------------------#


lift, span_loc = Lift_distribution_for_any_load_case(0.7,0.31641,241.9574)

moment = []                                                                             # Open list for moment


fout=open("bending.txt", "w")                                                           # Open text file

def moment_inner_int(i):
    return quad(lambda x: normal_force_for_integrating(x, 0.7, 0.31641, 241.9574, 1),
        i, 26.785,limit=50, epsabs=100)[0]


for i in np.arange(0,26.78,0.01):
    moment_result, error = quad(lambda y: moment_inner_int(y),i,26.78,limit=50, epsabs=100) + n*Pe*i + Me
    
    # n*Pe*i comes from a single integral of the engine weight (Pe) over the distance i. Me is the engine bending moment (no integral). n is the load factor 

    if i <= 9.37:                                           # Not including the engine contribution
        moment_result = moment_result - n*Pe*i - Me

    moment.append(moment_result)
    fout.write(str(moment_result))
    fout.write('\n')

    
plt.figure()
plt.plot(span_loc, moment, label="Bending Moment",color='red')
plt.xlabel("spanwise location")
plt.ylabel("Bending Moment")
plt.title("Bending Moment Dist.")
plt.show()
