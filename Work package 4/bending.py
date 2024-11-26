import math
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.integrate import quad

# Initial values
Pe = 9.81*9630        # Engine weight [N]
ye = 9.37             # Engine spanwise location [m]
Me = Pe*ye            # Engine bending moment around the root [Nm]
b = 53.57             # Wing span [m]
S = -1                # PLACEHOLDER DUMMY for shear
# Functions

def step(y):                   # Moment step function               
    if y<0 or y>b:
        print('span error')
        M = -1
    elif y<ye:
        M = 0
    else:
        M = Me
    return M

def bending_moment(S):          # Bending moment at a location y
    
    M = []
    span_loc = []

    for i in range (0,26786):
        y = i/1000                     # spanwise location
        moment =-step(y)        # moment as a function of shear + step
        M.append(0.001*moment)
        span_loc.append(y)
    
    return M,span_loc

M, span_loc = bending_moment(S)

plt.plot(span_loc,M)
plt.xlabel("Spanwise Position [m]")
plt.ylabel("Bending Moment [kNm]")
plt.title("Spanwise Bending Moment")
plt.legend()
plt.show()