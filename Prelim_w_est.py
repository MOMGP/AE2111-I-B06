from main import *
import numpy as np
CD_0_est=C_f*Swet_to_S #Eq 6.15
C_L_est=(np.pi*Ar*oswald_e*CD_0_est)**0.5 #Eq 6.11
C_D_est=2*CD_0_est #Eq 6.11
TSFC_est=22*Bypass_ratio**(-0.19) #Eq 6.23
specific_energy_eff=V_cr/TSFC_est #Eq 6.21 - making e_f*\eta_j=e_f,eff



