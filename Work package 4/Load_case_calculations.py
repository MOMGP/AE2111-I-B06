import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy import integrate
import Aero_loading_XFLR5
import math

CL_crit=[]
V_crit=[]
Load_factor_crit=[]
Rho_crit=[]
#critical load cases
S=265.2
rho_load_case=[]
Speed=np.genfromtxt('Load_cases',skip_header=1,usecols=(1),dtype=float,max_rows=30)
Weight=np.genfromtxt('Load_cases',skip_header=1,usecols=(2),dtype=float,max_rows=30)
Load_factor=np.genfromtxt('Load_cases',skip_header=1,usecols=(3),dtype=float,max_rows=30)
Alt=np.genfromtxt('Load_cases',skip_header=1,usecols=(4),dtype=str,max_rows=30)
CL_load_case=[]
def CL_calc(W,V,rho):
    CL_d=W*2/(rho*V**2*S)
    return CL_d


for i in range(len(Speed)):
    if(Alt[i]=='FL0'):
        rho_load_case.append(1.225)
    else:
        rho_load_case.append(0.31641)
    Speed[i]=Speed[i]*math.sqrt(1.225/rho_load_case[i])
    CL_load_case.append(CL_calc(Weight[i],Speed[i],rho_load_case[i]))
    #print(CL_load_case[i],Speed[i],Load_factor[i],rho_load_case[i])


#crit load cases
for i in range(len(Speed)):
    if i==20 or i==21 or i==25 or i==26:
        CL_crit.append(CL_load_case[i])
        V_crit.append(Speed[i])
        Load_factor_crit.append(Load_factor[i])
        Rho_crit.append(rho_load_case[i])
CL_crit=np.array(CL_crit)
V_crit=np.array(V_crit)
Load_factor_crit=np.array(Load_factor_crit)
Rho_crit=np.array(Rho_crit)

np.save('Load_case_arrays\\CL_crit.npy',CL_crit)
np.save('Load_case_arrays\\V_crit.npy',V_crit)
np.save('Load_case_arrays\\Load_factor_crit.npy',Load_factor_crit)
np.save('Load_case_arrays\\Rho_crit.npy',Rho_crit)