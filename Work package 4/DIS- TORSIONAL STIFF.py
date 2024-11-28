from math import cos,sin
import scipy.integrate
from sympy import symbols, Eq, solve # type: ignore
import matplotlib as mp 
import numpy as np
import scipy as scipy 
import matplotlib as plt
#torsional stiffenes function
#So the torsional stiffenes is defined as ts
#                                           c
#                    c1                                      c2
#       ------------------------------------------------------------------------------
#       |                                   |                                        |
#       |                                   |                                        |
#       |                h1                 |                 h2                     |
#       | <-------------------------------> | <------------------------------------> | b 
#     a |                                   |                                        |
#       |                                   |                                        |
#       |                                   |                                        |
#       |                                   |                                        |
#       -----------------------------------------------------------------------------
#                     d1                                       d2
#
#       <---------------------------------------------------------------------------->
#                                           h
#
#reference area
def get_A(a,b,h):
    A = 0.5 * (a+b)*h 
    return A  

#single cell 
def get_J_SS(h,a,b,c,d,t1,t2,t3,t4):
    sum = a/t1 + b/t2 + c/t3 + d/t4
    A = get_A(a,b,h) 
    J = (4* A**2)/sum
    return J 

#rot = rate of twist 
#multi cell
def get_J_MS(a, b1, b, h1, h2, G, t1, t2, t3, t4, t5, t6, t7, c1, c2, d1, d2):
    # Define symbolic variables
    rot, q1, q2 = symbols('rot q1 q2')

    # Left side
    A1 = get_A(a, b1, h1)  # Assuming get_A is a provided function
    eq1 = Eq(rot, (1 / (2 * A1)) * (((q1 - q2) * b1 / (G * t4)) + (q1 * c1 / (G * t5)) + (q1 * a / (G * t6)) + (q1 * d1 / (G * t7))))

    # Right side
    A2 = get_A(b1, b, h2)
    eq2 = Eq(rot, (1 / (2 * A2)) * ((q2 * b / (G * t2)) + (q2 * c2 / (G * t1)) + ((q2 - q1) * b1 / (G * t4)) + (q2 * d2 / (G * t3))))

    # Third formula
    eq3 = Eq(1, 2 * A1 * q1 + 2 * A2 * q2)

    # Solve the system of equations
    solutions = solve([eq1, eq2, eq3], (rot, q1, q2))

    # Extract rot for the final J
    rot_solution = solutions[rot]
    J = (1 / G) * rot_solution
    return J


#now that we got the J's we need to find a way to let the program determain when to switch to a single cell J calculation 
def Torsional_Stiffness(h,h1,h2,G,a,b,c,d,b1,c1,c2,d1,d2,t1, t2, t3, t4, t5, t6, t7):
    if h1 > 0:
        tor_stiff = G*get_J_MS(a, b1, b, h1, h2, G, t1, t2, t3, t4, t5, t6, t7, c1, c2, d1, d2)
    else: 
        tor_stiff = G*get_J_SS(h,a,b,c,d,t1,t2,t3,t4)
    return tor_stiff 



# T(x) = tor_stiff * theta'
# theta' = T(x)/tor stiff 

# def get_twist_angle(tor_stiff,x):
#    hws = 26.785
#    T = lambda x: #...    
#    theta = 1/tor_stiff * scipy.integrate.quad(T,0,x)
#    return theta


print(Torsional_Stiffness(1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))

#graphing 
hws = 26.785
x = np.arange(0,hws,0.01)
 