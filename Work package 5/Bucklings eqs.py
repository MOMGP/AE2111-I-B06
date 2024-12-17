import numpy as np
import math
from #Kc_function.py import
#Parameters
k_v = -1 # ASSUMED factor of average shear compared to max 
K = 1 # End condition factor of the column, assumed to be pinned; 1. This is because we assume the stringers to be continuous.
A = 
E = 72400000000 # youngs modulus
v = 0.33 # poisson ratio for our material

#For shear, only the webs are considered, not flanges :P 
def shear(h1,t1,h2,t2,h3,t3,V,q1,q2):
    tau_average = V / (h1 * t1 + h2 * t2 + h3 * t3)
    tau_max = tau_average * k_v # and xxx * q1,q2 --> need to super impose shear flow into the shear from normal force. 
    return tau_max
#Critical buckling from shear for the web
def crit_web_shear(t,a,b):
    #k_s = ks(a,b)
    crit_web_shear = math.pi() ** 2 * k_s * E / (12 * (1 - v ** 2)) * (t / b) **2
    return crit_web_shear

def crit_skin_stress(t,a,b):
    #k_c = kc(a,b)
    crit_stress = math.pi() ** 2 * k_c * E / (12 * (1 - v ** 2)) * (t / b) **2
    return crit_stress

def crit_stringer_stress(K,L,I,A): #L= length, I = mmoi, A = area
    #I could be function dunno :3
    crit_column_stress = K * math.pi() * E * I / (L ** 2 * A)
    return crit_column_stress

