import numpy
from matplotlib import pyplot as plt
#CONSTANTS
lapse_rate = -0.0065 #K m^-1
rho_0 = 1.225 #kg m^-3
R_air = 287.05 #J kg^-1 K^-1
T_0 = 288.15 #K

#Airplane constraints
M_p_max = 49442 #kg
M_cr = 0.82
h_cr = 39000 #ft
d_to = 2790 #m
d_ld = 1856 #m 
range_dm = 13797 #km
payload_dm = 27669 #kg
range_MTOW_full_fuel = 13983 #km
M_p_max_full_fuel = 26308 #km
r_ferry = 15811 #km



#ASSUMPTIONS
assumed_max_CL=2
assumed_maximum_takeoff_mass=150000 #kg
assumed_aproach_speed=74.163 #aproach speed
assumed_landing_mass_fraction=0.8 #just assumed, no calculations
