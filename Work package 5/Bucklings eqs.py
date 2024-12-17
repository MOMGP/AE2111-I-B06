import numpy as np
import math

#Parameters
k_v = -1 # ASSUMED factor of average shear compared to max 
k_s = -1 # ASSUMED FROM GRAPH (figure 16 in reader)
# k_s might become a function if we 
E = 100 # youngs modulus
v = 0.33 # poisson ratio for our material

#For shear, only the webs are considered, not flanges :P 
def shear(h1,t1,h2,t2,h3,t3,V,q1,q2):
    tau_average = V / (h1 * t1 + h2 * t2 + h3 * t3)
    tau_max = tau_average * k_v # and xxx * q1,q2 --> need to super impose shear flow into the shear from normal force. 
    return tau_max
#Critical buckling from shear for the web
def crit_web_shear(k_s,t,b):
    k_s = placeholder_eq(a,b,y)
    crit_web_shear = math.pi() ** 2 * k_s * E / (12 * (1 - v ** 2)) * (t / b) **2
    return crit_web_shear