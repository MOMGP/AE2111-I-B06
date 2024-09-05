import numpy
from matplotlib import pyplot as plt
from Prelim_w_est import *
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
M_cr = 0.82
speed_sound_cr=295.07 #m/s
h_cr = 39000 #ft
cruise_alt=11887.2#m
d_to = 2790 #m
d_ld = 1856 #m
range_dm = 13797 #km
payload_dm = 27669 #kg
range_MTOW_full_fuel = 13983 #km
M_p_max_full_fuel = 26308 #km
r_ferry = 15811 #km
landing_field_len=1856 #m
V_cr=speed_sound_cr*M_cr

#ASSUMPTIONS
assumed_max_CL=2
C_f=0.0026
Swet_to_S=6
CD_0=CD_0_est
Bypass_ratio=7.4
Ar=9.1#Aspect ratio
oswald_e=0.7829#Oswald efficiency factor
TSFC=TSFC_est
C_L=C_L_est
C_D=C_D_est
assumed_mass_fraction_cruise=0.95
assumed_maximum_takeoff_mass=150000 #kg
assumed_aproach_speed=74.163 #aproach speed
assumed_landing_mass_fraction=0.8 #just assumed, no calculations
assumed_landing_field_coeff=0.45 #Fromn adsee book, eq 7.9, for CS/FAR-25



