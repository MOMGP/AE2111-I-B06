
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
g_0=9.81


#Airplane constraints
M_p_max = 49442 #kg
M_cr = 0.82
speed_sound_cr=295.07#m/s
h_cr = 39000 #ft
cruise_alt=11887#m
V_cruise=M_cr*speed_sound_cr
d_to = 2790 #m
d_ld = 1856 #m
range_dm = 13797 #km
payload_dm = 27669 #kg
range_MTOW_full_fuel = 13983 #km
M_p_max_full_fuel = 26308 #km
r_ferry = 15811 #km
landing_field_len=1865 #m
Bypass_ratio=7.4 # As from REQ-PRP-01
Ar=9.1#Aspect ratio as from REQ-WNG-01
e=0.78#Oswald efficiency factor as from REQ-WNG-02
CD_0=0.0156#during cruise as from REQ-WNG-03
assumed_mass_fraction_cruise=0.93
ROC=0.5#m/s requirement
climb_alt_req=12000#m
CL_highest_climb=math.sqrt(CD_0*math.pi*Ar*e)
CD_highest_climb=2*CD_0
take_off_field_len=2790#m
CL_to=2.1
#ASSUMPTIONS
assumed_max_CL=2.5#when landing
assumed_maximum_takeoff_mass=359339 #kg
assumed_aproach_speed=70 #aproach speed
assumed_landing_mass_fraction=0.66
assumed_landing_field_coeff=0.45 #Fromn adsee book, eq 7.9, for CS/FAR-25
assumed_mass_fraction_climb=0.93

#CS FAR 119 and 121(a-d) climb gradient req
c_119=3.2/100
c_121_a=0
c_121_b=2.4/100
c_121_c=1.2/100
c_121_d=2.1/100

#take off assumptions
kt=0.85 #for jet airplanes
obstacle_h=11#m, assumed
