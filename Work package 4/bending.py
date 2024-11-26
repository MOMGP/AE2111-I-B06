import math
import numpy as np
import matplotlib.pyplot as plt

# Initial values
Pe = 9.81*9630        #Engine weight [N]
ye = 9.37             #Engine spanwise location [m]
Me = Pe*ye            #Engine bending moment around the root [Nm]
b = 26.785            #Wing half span [m]
S = -1                #PLACEHOLDER DUMMY for shear
# Functions

def step(y):                    #Moment step function               
    if y<0 or y>b:
        print('span error')
        M=-1
    elif y<ye:
        M=0
    else:
        M=Me
    return M

def bending_moment(S):         #Bending moment at a location y
    
    M=[]
    span_loc=[]

    for i in range (0,26786):
        y=i/1000                     #spanwise location
        moment=1-step(y)             #moment as a function of shear + step
        M.append(moment)
        span_loc.append(y)
    
    return M,span_loc

M, span_loc = bending_moment(S)

plt.plot(span_loc,M)
plt.xlim([0,27])
plt.show()