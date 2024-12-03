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
from Shear_force import shear_force_for_integrating

#--------------------------------Initial values--------------------------------#
n = 1
Pe = 9.81*9630/1000*n        # Engine weight [N]
ye = 9.37             # Engine spanwise location [m]
Me = Pe*ye            # Engine bending moment around the root [Nm]
b = 53.57             # Wing span [m]
CL_d = 3
rho = 1.225
V = 62
#-------------------------------------------------------------------------------#

fout=open('bending.txt','w')
moment = []  # Open list for moment

def moment_inner_int(x,CL_d,rho,V,n):
    return quad(lambda x: normal_force_for_integrating(x, 0.7, 0.31641, 241.9574, 1),
        i, 26.785,limit=500, epsabs=100)[0]

def moment_at_position(x,CL_d,rho,V,n):
    moment_result, error =sp.integrate.quad(shear_force_for_integrating,x,26.785,args=(CL_d,rho,V,n),limit=20,epsabs=100)
    moment_result=moment_result*(-1)
    moment.append(moment_result)
    print(moment_result)
    return moment_result


#START OF SHEAR CALCULATION PAPA MARIO TAKE THIS IN YOUR CODE

def moment_at_full_position(CL_d,rho,V,n):
    moment_total=0
    position=[]
    moment=[]
    shear=[]
    for i in np.arange(0,26.785,0.01):
        res=shear_force_for_integrating(i, 0.7, 0.31641, 241.9574, 1)
        moment_total-=res*0.01
        position.append(i)
        shear.append(res)

    for i in np.arange(0,26.785,0.01):
        moment.append(moment_total)
        pos=int((i*100))
        moment_total+=shear[pos]*0.01
        fout.write(str(moment_total))
        fout.write('\n')
moment_at_full_position(0.7,0.31641,241.9574, 1)
'''
def moment_at_position_updated(x):
    pos=int(x*100-1)
    moment_result=moment[pos]+shear[pos]*0.01
    position.append(x)
    moment.append(moment_result)
    print(pos,shear[pos],moment_result)
for i in np.arange(0.01,26.78,0.01):
    moment_at_position_updated(i)
plt.figure()
plt.plot(position,moment)
plt.show()
plt.figure()
plt.plot(position,shear)
plt.show()

pos=[]

for i in np.arange(0,26.78,1):
    print(i)
    pos.append(i)
    moment_at_position(i,0.7,0.31641,241.9574, 1)
plt.figure()
plt.plot(pos, moment, label="Bending Moment",color='red')
plt.xlabel("spanwise location")
plt.ylabel("Bending Moment")
plt.title("Bending Moment Dist.")
plt.show()
'''