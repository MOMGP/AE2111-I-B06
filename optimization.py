#Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#from Matching_diagram_Andrei #FINISH 

#Data
file_path_CL_alpha = "C:\\Users\\Equipo\\Downloads\\Polar_Graph_cl_alpha.csv"
C_l_v_alpha = pd.read_csv(file_path_CL_alpha).to_numpy()
C_l_v_alpha = np.array([C_l_v_alpha[:,0], C_l_v_alpha[:,1]])

file_path_CL_v_CD = "C:\\Users\\Equipo\\Downloads\\Polar_Graph_CL_CD.csv"
CL_v_CD = pd.read_csv(file_path_CL_v_CD).to_numpy()
CL_v_CD = np.array([CL_v_CD[:,0], CL_v_CD[:,1]])

file_path_CM_v_alpha = "C:\\Users\\Equipo\\Downloads\\Polar_Graph_CM.csv"
CM_v_alpha = pd.read_csv(file_path_CM_v_alpha).to_numpy()
CM_v_alpha = np.array([CM_v_alpha[:,0], CM_v_alpha[:,1]])

#Parameters - constants
M_CR = 0.82
B = 10
C_f = 0.0026
gamma=1.4
R=287.05
S_wet_to_S_w = 6
C_D0 = C_f*S_wet_to_S_w
g=9.80665
t_E = 1800
h_cr = 39000*0.3048
T_to_W = 0.329
wing_loading = 7000
T_CR = 288.15-6.5*10**(-3)*min(h_cr, 11000)
V_CR = M_CR*np.sqrt(T_CR*gamma*R) #m s-1
p_CR =101325*(T_CR/288.15)**(g/(6.5*10**(-3)*R))*np.exp(-g/(R*T_CR)*(h_cr-11000))
rho_CR = p_CR/(T_CR*R)
rho_kerosene = 820
C_T = 22*B**(-0.19)
e_f = 43
eta_j = V_CR/(C_T*e_f)
OEM_to_MTOM = 0.48203
L_to_D_CR = 17.16 #  CHANGE ONCE GET VALUES
R_nom = np.array([13797000, 13983000, 15811000])
R_lost = 1/0.7*L_to_D_CR*(h_cr+V_CR**2/(2*g))
f_con = 0.05
R_div = 400000
R_endu = V_CR*t_E
R_eq = (R_nom+R_lost)*(1+f_con)+1.2*R_div+R_endu
m_pl = np.array([27669, 26308, 0])
m_f_to_MTOM = 1-np.exp(-R_eq*g/(eta_j*e_f*10**6*L_to_D_CR))
MTOM = np.array([m_pl[0]/(1-m_f_to_MTOM[0]-OEM_to_MTOM), m_pl[1]/(1-m_f_to_MTOM[1]-OEM_to_MTOM), m_pl[1]/(1-m_f_to_MTOM[1]-OEM_to_MTOM)-m_pl[1]])
m_f = m_f_to_MTOM*MTOM
V_f_req = np.max(m_f)/rho_kerosene
Lambda_c4 = np.rad2deg(np.arccos(1.16/(M_CR+0.5)))
taper = 0.2*(2-np.deg2rad(Lambda_c4))
eta_tank = 0.55
C_l_max = 2.12
C_L_max_to_C_l_max = 0.8
C_L_max = C_L_max_to_C_l_max*C_l_max
C_d_cruise = None
alpha_CR = None
print(V_CR)
#Formulas
def SAR(S, C_D):
    return V_CR/(0.5*rho_CR*V_CR**2*S*C_D*C_T)

def Lambda_n(AR, n_percent):
    return np.rad2deg(np.arctan(np.tan(np.deg2rad(Lambda_c4))-(4/AR)*((n_percent-25)/100*(1-taper)/(1+taper))))

def wingspan_eff(AR, Lambda_LE):
    return 4.61*(1-0.045*AR**0.68)*(np.cos(np.deg2rad(Lambda_LE))**0.15)-3.1

def wing_fuel_volume(eta_tank, t_to_c, S_w, AR):
    return 0.9*eta_tank*t_to_c*S_w**1.5*AR**(-0.5)

def C_L_alpha(AR, Lambda_c2):
    beta = np.sqrt(1-M_CR**2)
    return 2*np.pi*AR/(2+np.sqrt((AR*beta/0.95)**2*(1+np.tan(Lambda_c2)**2/beta**2)+4))

def Oswald_eff_factor(AR, Lambda_LE):
    return 4.61*(1-0.045*AR**0.68)*(np.cos(np.deg2rad(Lambda_LE)))**0.15-3.1

#Variables
AR = 10.2
S = 363.53
e = 0.764
Lambda_c2 = Lambda_n(AR, 50)
Lambda_LE = Lambda_n(AR, 0)
Lambda_TE = Lambda_n(AR, 100)

#Requirements
reqs=[]
current=[]

file_path_CL_alpha = "C:\\Users\\Equipo\\Downloads\\Polar_Graph_cl_alpha.csv"
C_l_v_alpha = pd.read_csv(file_path_CL_alpha).to_numpy()
C_l_v_alpha = np.array([C_l_v_alpha[:,0], C_l_v_alpha[:,1]])


# plt.plot(C_l_v_alpha[0], C_l_v_alpha[1])
# plt.show()

# plt.plot(CL_v_CD[0], CL_v_CD[1])
# plt.show()

# plt.plot(CM_v_alpha[0], CM_v_alpha[1])
# plt.show()


for i in range(len(CL_v_CD[0])-1):
    if (CL_v_CD[1,i]-0.857)*(CL_v_CD[1,i+1]-0.857)<0:
        C_d_cruise=CL_v_CD[0,i] #gives cruise airfoil drag
        break

for i in range(len(C_l_v_alpha[0])-1):
    if (C_l_v_alpha[1,i]-0.857)*(C_l_v_alpha[1,i+1]-0.857)<0:
        alpha_CR= (C_l_v_alpha[0,i]) #gives cruise angle
        break