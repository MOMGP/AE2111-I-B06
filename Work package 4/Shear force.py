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
W_engine = 9630*9.81 

shear_force_distribution = []
internal_shear_force = []


def half_span(b):
    L = b/2
    return L

#from XFLR5
#using x-dir as along the span
#def integrand(x):
 #   return Lift_distribution_for_any_load_case

'''
#integrating the loading distribution
x = Symbol('x')
def shear_integration(t):
    return integrate(t, x)
    '''
#-------------********-------------*********------------*********-------------********-------------*********------------*********#
'''
def step(y):          # Shear force step function            
    if y_position<0 or y_position>half_span(b):
        S_engine = 0
    elif y_position<y_engine:
        S_engine = 0
    else:
        S_engine = W_engine
    return S_engine

#appending to lists
y_position = np.arange(0, half_span(b), 1).tolist()

x = sp.Symbol('x')
#Lift_distribution_for_any_load_case(CL_d,rho,V)
shearequation = shear_integration(2*x)
'''
'''for i in range(len(y_position)):
    if y_position[i]<0 or y_position[i]>half_span(b):
       
    elif y_position[i]<y_engine:
        
    else:
        internal_shear_force.append(shearequation.subs(x, var_sub)-step(y_position))
'''
"""for i in range(len(y_position)):
    var_sub = i
    if y_position[i]<0 or y_position[i]>half_span(b):
       internal_shear_force.append(shearequation.subs(x, var_sub))
    elif y_position[i]<y_engine:
        internal_shear_force.append(shearequation.subs(x, var_sub))
    else:
        internal_shear_force.append(shearequation.subs(x, var_sub)-step(y_position))
        """
    


# print(Lift_for_integrating(26.785,0.7,0.31641,241.9574))

lift, span_loc = Lift_distribution_for_any_load_case(0.7,0.31641,241.9574)
shear_force = []
L_error_list = []
total_lift,err=sp.integrate.quad(Lift_for_integrating,0,26.785,args=(0.7,0.31641,241.9574),limit=50, epsabs=100)
total_normal,err_normal=sp.integrate.quad(normal_force_for_integrating,0,26.785,args=(0.7,0.31641,241.9574),limit=50, epsabs=100)
#print(total_lift)
fout=open("shear.txt", "w")
for i in np.arange(0,26.78,0.01):
    shear_result,L_error_result=sp.integrate.quad(Lift_for_integrating,0,i,args=(0.7,0.31641,241.9574),limit=50, epsabs=100)
    shear_result= 94.4703+(-total_normal+shear_result)/1000

    if i >= 9.37:
        shear_result = (shear_result - 94.4703)

    shear_force.append(shear_result)
    shear_result=str(shear_result)
    fout.write(shear_result)
    fout.write('\n')
    L_error_list.append(L_error_result)


def shear_force_for_integrating(x,CL_d,rho,V,n):
    span_locs = x
    shear_result,L_error_result=sp.integrate.quad(Lift_for_integrating,0,x,args=(CL_d,rho,V,n),limit=50, epsabs=100)
    shear_result= 94.4703+(-total_normal+shear_result)/1000
    if i >= 9.37:
            shear_result = (shear_result - 94.4703)
    return shear_result

print(shear_force_for_integrating(10,0.7,0.31641,241.9574))


plt.figure()
plt.plot(span_loc, shear_force, label="Shear Force",color='purple')
plt.xlabel("spanwise location")
plt.ylabel("shear force")
plt.title("shear force dist.")
plt.show()
