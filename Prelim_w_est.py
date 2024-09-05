import sklearn.linear_model
from main import *
import numpy as np
import scipy as sp
import sklearn
font_size_full_scr=12
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

C_f=0.0026
Swet_to_S=6
Bypass_ratio=7.4
Ar=9.1#Aspect ratio
oswald_e=0.7829#Oswald efficiency factor
CD_0_est=C_f*Swet_to_S #Eq 6.15
C_L_est=(np.pi*Ar*oswald_e*CD_0_est)**0.5 #Eq 6.11
C_D_est=2*CD_0_est #Eq 6.11
TSFC_est=22*Bypass_ratio**(-0.19) #Eq 6.23
specific_energy_eff_est=V_cruise/TSFC_est #Eq 6.21 - making e_f*\eta_j=e_f,eff


#Whole main (proper)

"""
import numpy
from matplotlib import pyplot as plt
from Prelim_w_est import *
#Plotting things


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

"""
MTOM_vals=np.array([242000, 242000, 227900, 254700, 186900, 297500, 276500, 257000])
OEM_vals=np.array([120600, 129400, 120000, 129000, 90000, 138100, 131000, 120228])
Payload_vals=np.array([49400, 45600, 41100, 53000, 43800, 59000, 48150, 49400])

MTOM_vals = MTOM_vals.reshape(-1, 1)
OEM_vals = OEM_vals.reshape(-1, 1)
Payload_vals = Payload_vals.reshape(-1, 1)

lin_MTOM_OEM=sklearn.linear_model.LinearRegression()
reg_MTOM_OEM=lin_MTOM_OEM.fit(MTOM_vals, OEM_vals)
lin_MTOM_pay=sklearn.linear_model.LinearRegression()
reg_MTOM_pay=lin_MTOM_pay.fit(MTOM_vals, Payload_vals)
lin_pay_OEM=sklearn.linear_model.LinearRegression()
reg_pay_OEM=lin_pay_OEM.fit(Payload_vals, OEM_vals)
"""
#MTOM - OEM
ax=plt.subplot()
ax.scatter(MTOM_vals, OEM_vals, rasterized=True)
ax.plot(
    [MTOM_vals.min(), MTOM_vals.max()],
    [reg_MTOM_OEM.coef_[0] * MTOM_vals.min() + reg_MTOM_OEM.intercept_,
     reg_MTOM_OEM.coef_[0] * MTOM_vals.max() + reg_MTOM_OEM.intercept_],
    '--', color='black')
formula_text = f"Regression: OEW = {reg_MTOM_OEM.coef_[0][0]:.3f} * MTOW + {reg_MTOM_OEM.intercept_[0]:.0f}\n$r^2$ = {np.round(sp.stats.pearsonr(MTOM_vals.flatten(), OEM_vals.flatten())[0]**2,2)}"
ax.text(0.025, 0.975, formula_text, transform=ax.transAxes, fontsize=font_size_full_scr, verticalalignment='top')
ax.set_xlabel("MTOW [kg]", fontsize=font_size_full_scr)
ax.set_ylabel("OEW [kg]", fontsize=font_size_full_scr)
ax.grid(linestyle = '-', linewidth=0.1)
plt.savefig("MTOW vs. OEW Prelim.pdf", format="pdf")
plt.close()

#MTOM - payload
ax=plt.subplot()
ax.scatter(MTOM_vals, Payload_vals, rasterized=True)
ax.plot(
    [MTOM_vals.min(), MTOM_vals.max()],
    [reg_MTOM_pay.coef_[0] * MTOM_vals.min() + reg_MTOM_pay.intercept_,
     reg_MTOM_pay.coef_[0] * MTOM_vals.max() + reg_MTOM_pay.intercept_],
    '--', color='black')
formula_text = "Regression: $M_{pl}$" +f"= {reg_MTOM_pay.coef_[0][0]:.3f} * MTOW + {reg_MTOM_pay.intercept_[0]:.0f}\n$r^2$ = {np.round(sp.stats.pearsonr(MTOM_vals.flatten(), Payload_vals.flatten())[0]**2,2)}"
ax.text(0.025, 0.975, formula_text, transform=ax.transAxes, fontsize=font_size_full_scr, verticalalignment='top')
ax.set_xlabel("MTOW [kg]", fontsize=font_size_full_scr)
ax.set_ylabel("$M_{pl}$ [kg]", fontsize=font_size_full_scr)
ax.grid(linestyle = '-', linewidth=0.1)
plt.savefig("MTOW vs. $M_{pl}$ Prelim.pdf", format="pdf")
plt.close()

#Payload - OEM
ax=plt.subplot()
ax.scatter(Payload_vals, OEM_vals, rasterized=True)
ax.plot(
    [Payload_vals.min(), Payload_vals.max()],
    [reg_pay_OEM.coef_[0] * Payload_vals.min() + reg_pay_OEM.intercept_,
     reg_pay_OEM.coef_[0] * Payload_vals.max() + reg_pay_OEM.intercept_],
    '--', color='black')
formula_text = "Regression: OEW" +f"= {reg_pay_OEM.coef_[0][0]:.3f} *"+ "$M_{pl}$"+f" + {reg_pay_OEM.intercept_[0]:.0f}\n$r^2$ = {np.round(sp.stats.pearsonr(Payload_vals.flatten(), OEM_vals.flatten())[0]**2,2)}"
ax.text(0.025, 0.975, formula_text, transform=ax.transAxes, fontsize=font_size_full_scr, verticalalignment='top')
ax.set_xlabel("$M_{pl}$ [kg]", fontsize=font_size_full_scr)
ax.set_ylabel("OEW [kg]", fontsize=font_size_full_scr)
ax.grid(linestyle = '-', linewidth=0.1)
plt.savefig("$M_{pl}$ vs. OEW Prelim.pdf", format="pdf")
plt.close()
"""