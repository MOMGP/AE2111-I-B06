from math import cos,sin
from sympy import symbols, Eq, solve # type: ignore
import numpy as np
import scipy as scipy 
from matplotlib import pyplot as plt
from geometry_WP4_2 import get_points_along_spanwise
from get_points import get_points,get_geom_from_points

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
"""
Geometry:
 ______________4_____________ ______________5_____________ 
|                            |                            |
|                            |                            |
|                            |                            |
1                            3                            6
|                            |                            |
|                            |                            |
|_____________2______________|_____________7______________| 

List_output = [[(x_s, y_s), (x_e,y_e), t]]
"""
def translate(val_list):
    if len(val_list)==4:

        h = np.abs(val_list[0][0][0] - val_list[2][0][0])
        a = np.sqrt((val_list[0][0][0]-val_list[0][1][0])**2+ (val_list[0][0][1]-val_list[0][1][1])**2)
        b = np.sqrt((val_list[2][0][0]-val_list[2][1][0])**2+ (val_list[2][0][1]-val_list[2][1][1])**2)
        c = np.sqrt((val_list[3][0][0]-val_list[3][1][0])**2+ (val_list[3][0][1]-val_list[3][1][1])**2)
        d = np.sqrt((val_list[1][0][0]-val_list[1][1][0])**2+ (val_list[1][0][1]-val_list[1][1][1])**2)
        t1 = val_list[3][2]
        t2 = val_list[2][2]
        t3 = val_list[1][2]
        t4 = val_list[0][2]
        h1 = 0
    

        return [h1,h,a,b,c,d,t1,t2,t3,t4]
    else:
        h1 = np.abs(val_list[0][0][0] - val_list[2][0][0])
        h2 = np.abs(val_list[2][0][0] - val_list[5][0][0])
        a = np.sqrt((val_list[0][0][0]-val_list[0][1][0])**2+ (val_list[0][0][1]-val_list[0][1][1])**2)
        b = np.sqrt((val_list[5][0][0]-val_list[5][1][0])**2+ (val_list[5][0][1]-val_list[5][1][1])**2)
        c1 = np.sqrt((val_list[3][0][0]-val_list[3][1][0])**2+ (val_list[3][0][1]-val_list[3][1][1])**2)
        c2 = np.sqrt((val_list[4][0][0]-val_list[4][1][0])**2+ (val_list[4][0][1]-val_list[4][1][1])**2) 
        d1 = np.sqrt((val_list[1][0][0]-val_list[1][1][0])**2+ (val_list[1][0][1]-val_list[1][1][1])**2) 
        d2 = np.sqrt((val_list[6][0][0]-val_list[6][1][0])**2+ (val_list[6][0][1]-val_list[6][1][1])**2)
        b1 = np.sqrt((val_list[2][0][0]-val_list[2][1][0])**2+ (val_list[2][0][1]-val_list[2][1][1])**2)

        t1 = val_list[4][2]
        t2 = val_list[5][2]
        t3 = val_list[6][2]
        t4 = val_list[2][2]
        t5 = val_list[3][2]
        t6 = val_list[0][2]
        t7 = val_list[1][2]


        return [h1, a, b1, b, h2, t1, t2, t3, t4, t5, t6, t7, c1, c2, d1, d2]
               
#a = translate([
#     [[0, 0], [0, 1], [0.01]],
#     [[0, 1], [1, 1], [0.01]],
#     [[1, 1], [1, 0], [0.01]],
#     [[1, 0], [0, 0], [0.01]],
#     [[1, 1], [2, 1], [0.01]],
#     [[2, 1], [2, 0], [0.01]],
#     [[2, 0], [1, 0], [0.01]]
# ])
# print(a[1])

# val_list = [
#     [[0, 0], [0, 1], 0.01],
#     [[0, 1], [1, 1], 0.01],
#     [[1, 1], [1, 0], 0.01],
#     [[1, 0], [0, 0], 0.01],
#     [[1, 1], [2, 1], 0.01],
#     [[2, 1], [2, 0], 0.01],
#     [[2, 0], [1, 0], 0.01]
# ]


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
def get_J_MS(h1, a, b1, b, h2, G, t1, t2, t3, t4, t5, t6, t7, c1, c2, d1, d2):
    # Define symbolic variables
    rot, q1, q2 = symbols('rot q1 q2')

    # Left side
    A1 = get_A(a, b1, h1)
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
    J = 1 / (G * rot_solution)
    return J


#now that we got the J's we need to find a way to let the program determain when to switch to a single cell J calculation 


def get_torr_stiff_list(norm_wing_box_root, norm_stringers, end_third_spar, cond, x):

    y = []
    # for i in range(0,26785,500):
    G = 28*10**9
    #     step = i*0.001
    val_list = get_points_along_spanwise(norm_wing_box_root, norm_stringers, x, end_third_spar, trunctated=cond)[0]
    values = translate(val_list)
    if values[0] > 0:
        tor_stiff = G * get_J_MS(values[0], values[1],values[2],values[3],values[4],G,values[5],values[6],values[7],values[8],values[9],values[10],values[11],values[12],values[13],values[14],values[15])
    else: 
        tor_stiff = G * get_J_SS(values[1], values[2],values[3],values[4],values[5],values[6],values[7],values[8],values[9])#[h1,h,a,b,c,d,t1,t2,t3,t4]
    return tor_stiff
    #     y.append(tor_stiff)#h,a,b,c,d,t1,t2,t3,t4
    # return y
    return 
span = 2*26.785
spar1_x=0.2
spar2_x=0.5
spar3_x=0.7
x_y_y = get_points(spar1_x, spar2_x, spar3_x, 1)
root_geom = get_geom_from_points(x_y_y, [0.0005 for i in range(7)])
#print(get_points_along_spanwise(root_geom, [[[0,0],0]], span/6, True))
# print(get_torr_stiff_list(root_geom, [[[0,0],0]], span/6, True))









# T(x) = tor_stiff * theta'
# theta' = T(x)/tor stiff 

# def get_twist_angle(tor_stiff,x):
#    hws = 26.785
#    T = lambda x: #...    
#    theta = 1/tor_stiff * scipy.integrate.quad(T,0,x)
#    return theta
#graphing 
hws = 26.785
x = np.arange(0,hws,0.5)
y = []
for i in x:
    y.append(get_torr_stiff_list(root_geom, [[[0,0],0]], span/6, True, i))

y = np.array(y)
plt.plot(x,y)
plt.plot([0,1], [1,2])
plt.show()
print("I'm going to kms if you dont run this, BITCH")
