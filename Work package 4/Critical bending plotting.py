import math
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from Aero_loading_XFLR5 import normal_force_for_integrating
from bending import  moment_at_full_position
rho_load_case2=[]
#critical load cases
S=265.2
Speed2=np.genfromtxt('Load_cases2',skip_header=1,usecols=(1),dtype=float,max_rows=30)
Weight2=np.genfromtxt('Load_cases2',skip_header=1,usecols=(2),dtype=float,max_rows=30)
Load_factor2=np.genfromtxt('Load_cases2',skip_header=1,usecols=(3),dtype=float,max_rows=30)
Alt2=np.genfromtxt('Load_cases2',skip_header=1,usecols=(4),dtype=str,max_rows=30)
CL_load_case2=[]
def CL_calc(W,V,rho):
    CL_d=W*2/(rho*V**2*S)
    return CL_d

for i in range(len(Speed2)):
    if(Alt2[i]=='FL0'):
        rho_load_case2.append(1.225)
    else:
        rho_load_case2.append(0.31641)
    Speed2[i]=Speed2[i]*math.sqrt(1.225/rho_load_case2[i])
    CL_load_case2.append(CL_calc(Weight2[i],Speed2[i],rho_load_case2[i]))
fout=open("bending_for_plotting.txt", "w")
for k in range (len(CL_load_case2)):
    CL_d=CL_load_case2[k]
    V=Speed2[k]
    n=Load_factor2[k]
    rho=rho_load_case2[k]
    shear_force=[]
    loc=[]
    moment=moment_at_full_position(CL_d,rho,V,n)
    fout.write("Load case ")
    fout.write(str(k+1))
    fout.write("\n")
    for i in range(0,len(moment)):
        fout.write(str(moment[i]))
        fout.write("\n")

