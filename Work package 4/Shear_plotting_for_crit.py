import math
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from Load_case_calculations import *
from Shear_force import shear_force_for_integrating
from Aero_loading_XFLR5 import normal_force_for_integrating
import math
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
    if i==25:
        print(CL_load_case2[i],rho_load_case2[i],Speed2[i],Load_factor2[i])
fout=open("Shear_for_plotting.txt", "w")
for k in range (len(CL_load_case2)):
    fout.write("Load case ")
    fout.write(str(k+1))
    fout.write("\n")
    CL_d=CL_load_case2[k]
    V=Speed2[k]
    n=Load_factor2[k]
    rho=rho_load_case2[k]
    shear_force=[]
    if k==21:
        print(k+1,CL_d,rho,V,n)
    if k==25:
        print(k+1,CL_d,rho,V,n)
        '''
    loc=[]
    for i in np.arange(0, 26.785, 0.01):
        normal_int, err_normal = sp.integrate.quad(normal_force_for_integrating, i, 26.785, args=(CL_d, rho, V, n),limit=50, epsabs=100)
        shear_result = (-normal_int) / 1000
        if i <= 9.37:
            shear_result = (shear_result + 94.4703 * n)
        fout.write(str(shear_result))
        fout.write('\n')
'''
        '''
    if(k==20):
        plt.plot(loc, shear_force,label='LC-21',color='red')
    elif(k==21):
        plt.plot(loc, shear_force,label='LC-22',color='green')
    elif(k==25):
        plt.plot(loc, shear_force,label='LC-26',color='black')
    elif(k==26):
        plt.plot(loc, shear_force,label='LC-27',color='blue')
    else:
        plt.plot(loc, shear_force,color='purple')

plt.legend(loc='lower right')
plt.savefig("Shear load case.pdf", format='pdf')
plt.show()
'''
