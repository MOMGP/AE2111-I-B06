
import numpy
import math
from matplotlib import pyplot as plt
#Plotting things
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

#CONSTANTS
lapse_rate = -0.0065 #K m^-1
rho_0 = 1.225 #kg m^-3
R_air = 287.05 #J kg^-1 K^-1
T_0 = 288.15 #K
p_0=101325 #Pa
gamma=1.4
theta_t_break=1.07


#Airplane constraints
M_p_max = 49442 #kg
M_cr = 0.77
speed_sound_cr=299.46 #m/s
h_cr = 39000 #ft
cruise_alt=10000#m
V_cruise=M_cr*speed_sound_cr
d_to = 2790 #m
d_ld = 1856 #m
range_dm = 13797 #km
payload_dm = 27669 #kg
range_MTOW_full_fuel = 13983 #km
M_p_max_full_fuel = 26308 #km
r_ferry = 15811 #km
landing_field_len=1480 #m
Bypass_ratio=11
Ar=10.6#Aspect ratio
e=0.78#Oswald efficiency factor
CD_0=0.018#during cruise
assumed_mass_fraction_cruise=0.95
ROC=4.7#m/s requirement
climb_alt_req=7400#m
CL_highest_climb=math.sqrt(CD_0*math.pi*Ar*e)
CD_highest_climb=2*CD_0
#ASSUMPTIONS
assumed_max_CL=2.3#when landing
assumed_maximum_takeoff_mass=150000 #kg
assumed_aproach_speed=68 #aproach speed
assumed_landing_mass_fraction=0.87 #just assumed, no calculations
assumed_landing_field_coeff=0.45 #Fromn adsee book, eq 7.9, for CS/FAR-25
assumed_mass_fraction_climb=0.95
