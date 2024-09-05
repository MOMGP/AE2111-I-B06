from main import *
import numpy as np
CD_0_est=C_f*Swet_to_S #Eq 6.15
C_L_est=(np.pi*Ar*oswald_e*CD_0_est)**0.5 #Eq 6.11
C_D_est=2*CD_0_est #Eq 6.11
TSFC_est=22*Bypass_ratio**(-0.19) #Eq 6.23
specific_energy_eff_est=V_cr/TSFC_est #Eq 6.21 - making e_f*\eta_j=e_f,eff




#ASSUMPTIONS - put in main later
"""
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
"""