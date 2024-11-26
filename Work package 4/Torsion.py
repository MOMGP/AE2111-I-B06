import math

Weight_engine = 9630*9.81 #N
Thrust_per_engine = 467060 #N
Lambda_LE = 0.54 #rad
C_r = 7.63 #m
b = 53.57 #m
taper = 0.3

Lambda_c2 = math.atan(math.tan(Lambda_LE) - 0.5*((2*C_r)/b)*(1-taper))
