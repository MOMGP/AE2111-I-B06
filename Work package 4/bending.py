import math
import numpy as np

# Initial values
Pe = 9.81*9630        #Engine weight [N]
ye = 9.37             #Engine spanwise location [m]
Me = Pe*ye            #Engine bending moment around the root [Nm]
b = 26.785            #Wing half span [m]

# Functions

def step(y):           #Moment step function               
    if y<0 or y>b:
        print('span error')
        M=-1
    elif y<ye:
        M=0
    else:
        M=Me
    return M

