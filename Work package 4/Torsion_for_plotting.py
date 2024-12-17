from Torsion import internal_torque_for_plotting
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import math
from Aero_loading_XFLR5 import AOA
rho_load_case2=[]
def CL_calc(W,V,rho):
    CL_d=W*2/(rho*V**2*S)
    return CL_d
#critical load cases
S=265.2
Speed2=np.genfromtxt('Load_cases2',skip_header=1,usecols=(1),dtype=float,max_rows=30)
Weight2=np.genfromtxt('Load_cases2',skip_header=1,usecols=(2),dtype=float,max_rows=30)
Load_factor2=np.genfromtxt('Load_cases2',skip_header=1,usecols=(3),dtype=float,max_rows=30)
Alt2=np.genfromtxt('Load_cases2',skip_header=1,usecols=(4),dtype=str,max_rows=30)
CL_load_case2=[]
for i in range(len(Speed2)):
    if(Alt2[i]=='FL0'):
        rho_load_case2.append(1.225)
    else:
        rho_load_case2.append(0.31641)
    Speed2[i]=Speed2[i]*math.sqrt(1.225/rho_load_case2[i])
    CL_load_case2.append(CL_calc(Weight2[i],Speed2[i],rho_load_case2[i]))
fout=open("torsion_for_plotting.txt", "w")
for k in range (len(CL_load_case2)):
    CL_d=CL_load_case2[k]
    V=Speed2[k]
    n=Load_factor2[k]
    rho=rho_load_case2[k]
    shear_force=[]
    fout.write("Load case ")
    fout.write(str(k+1))
    fout.write("\n")
    print("LOAD CASE",k+1,AOA(CL_d))
    torque=internal_torque_for_plotting(CL_d,rho,V,n)
    for i in range (len(torque)):
        fout.write(str(torque[i]))
        fout.write("\n")

