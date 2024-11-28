import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy import integrate
import Aero_loading_XFLR5
import math
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
print(Speed)
print(Weight)
print(Load_factor)
print(Alt)
for i in range(len(Speed)):
    if(Alt[i]=='FL0'):
        rho_load_case.append(1.225)
    else:
        rho_load_case.append(0.31641)
    Speed[i]=Speed[i]*math.sqrt(1.225/rho_load_case[i])
    CL_load_case.append(CL_calc(Weight[i],Speed[i],rho_load_case[i]))
for i in range (len(CL_load_case)):
    CL_d=CL_load_case[i]
    L=[]
    span_loc=[]
    L,span_loc=
