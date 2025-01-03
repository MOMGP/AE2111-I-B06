import math
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy.integrate import quad
from Aero_loading_XFLR5 import Lift_for_integrating
from Aero_loading_XFLR5 import Lift_distribution_for_any_load_case
from Aero_loading_XFLR5 import normal_force_for_integrating
#extra stuff
#x - y-direction along the wingspan
#y - internal shear force

CL_d = 3
rho = 1.225
V = 62
#span, b
b = 53.57
y_engine = 9.37 #[m]
W_engine = 9630*9.81/1000

shear_force_distribution = []
internal_shear_force = []


def half_span(b):
    L = b/2
    return 
    
#-------------********-------------*********------------*********-------------********-------------*********------------*********#

lift, span_loc = Lift_distribution_for_any_load_case(0.7,0.31641,241.9574)
shear_force_1 = []
shear_force_2 = []
L_error_list = []
normal_list =[]
total_lift,err=sp.integrate.quad(Lift_for_integrating,0,26.785,args=(0.7,0.31641,241.9574),limit=50, epsabs=100)
total_normal,err_normal=sp.integrate.quad(normal_force_for_integrating,0,26.785,args=(0.7,0.31641,241.9574,1),limit=50, epsabs=100)


fout=open("shear_n=-1.txt", "w")
for i in np.arange(0,26.78,0.01):
    normal_int_1,err_normal=sp.integrate.quad(normal_force_for_integrating,i,26.785,args=(0.7,0.31641,241.9574,-1),limit=50, epsabs=100)

    shear_result_1= (-normal_int_1)/1000

    if i <= 9.37:
        shear_result_1 = (shear_result_1 + 94.4703)

    shear_force_1.append(shear_result_1)
    shear_result_1=str(shear_result_1)
    fout.write(shear_result_1)
    fout.write('\n')
plt.plot(span_loc, shear_force_1, label="n = -1",color='purple') 

fout=open("shear_n=2.5.txt", "w")
for i in np.arange(0,26.78,0.01):
    normal_int_2,err_normal=sp.integrate.quad(normal_force_for_integrating,i,26.785,args=(0.7,0.31641,241.9574,2.5),limit=50, epsabs=100)

    shear_result_2= (-normal_int_2)/1000

    if i <= 9.37:
        shear_result_2 = (shear_result_2 + 94.4703)

    shear_force_2.append(shear_result_2)
    shear_result_2=str(shear_result_2)
    fout.write(shear_result_2)
    fout.write('\n')
plt.plot(span_loc, shear_force_2, label="n = 2.5",color='red') 


def shear_force_for_integrating(x,CL_d,rho,V,n):
    span_locs = x
    normal_int,err_normal=sp.integrate.quad(normal_force_for_integrating,x,26.785,args=(CL_d,rho,V,n),limit=50, epsabs=100)
    shear_result= (-normal_int)/1000

    if x <= 9.37:
        shear_result = (shear_result + 94.4703)
    return shear_result


def shear_force_full_values(CL_d,rho,V,n):
    fout = open("shear.txt", "w")
    for i in np.arange(0,26.785,0.01):
        normal_int, err_normal = sp.integrate.quad(normal_force_for_integrating,i, 26.785, args=(CL_d, rho, V, n),limit=50, epsabs=100)
        shear_result= (-normal_int)/1000
        if i <= 9.37:
            shear_result = (shear_result + 94.4703*n)
        fout.write(str(shear_result))
        fout.write('\n')


plt.xlabel("Spanwise Location [m]")
plt.ylabel("Shear Force [kN]")
plt.grid(True)
plt.legend()
plt.savefig('Shear_Force.pdf')
plt.show()
