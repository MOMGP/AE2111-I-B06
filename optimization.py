#Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Matching_diagram_Andrei import approach_speed_req, climb_gradient_CS25119, climb_gradient_CS25121a, climb_gradient_CS25121b, climb_gradient_CS25121c, climb_gradient_CS25121d, landing_field_req, climb_req_thrust_over_weight, cruise_req_thrust_over_weight, landing_field_length_req, returns_for_optimization #FINISH 
font_size_full_scr=14
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

#Data
'''
file_path_CL_alpha = "C:\\Users\\Equipo\\Downloads\\Polar_Graph_cl_alpha.csv"
C_l_v_alpha = pd.read_csv(file_path_CL_alpha).to_numpy()
C_l_v_alpha = np.array([C_l_v_alpha[:,0], C_l_v_alpha[:,1]])

file_path_CL_v_CD = "C:\\Users\\Equipo\\Downloads\\Polar_Graph_CL_CD.csv"
CL_v_CD = pd.read_csv(file_path_CL_v_CD).to_numpy()
CL_v_CD = np.array([CL_v_CD[:,0], CL_v_CD[:,1]])

file_path_CM_v_alpha = "C:\\Users\\Equipo\\Downloads\\Polar_Graph_CM.csv"
CM_v_alpha = pd.read_csv(file_path_CM_v_alpha).to_numpy()
CM_v_alpha = np.array([CM_v_alpha[:,0], CM_v_alpha[:,1]])
'''
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
C_T = 22*B**(-0.19)*10**(-6)
e_f = 43
eta_j = V_CR/(C_T*10**(6)*e_f)
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
MTOM_crit_case = np.max(MTOM)
V_f_req = np.max(m_f)/rho_kerosene
Lambda_c4 = np.rad2deg(np.arccos(1.16/(M_CR+0.5)))
taper = 0.2*(2-np.deg2rad(Lambda_c4))
eta_tank = 0.55
C_l_max = 2.12
C_L_max_to_C_l_max = 0.8
C_L_max = C_L_max_to_C_l_max*C_l_max
C_d_cruise = None #determined below
alpha_CR = None #determined below  
landing_mass_fraction = np.max(1-m_f_to_MTOM+(1-np.exp(-(1.2*R_div+R_endu)*g/(eta_j*e_f*10**6*L_to_D_CR)))*(1-m_f_to_MTOM+1-np.exp(-(1.2*R_div+R_endu)*g/(eta_j*e_f*10**6*L_to_D_CR))))
Delta_CL_max = 0.804
CL_max_landing =  C_L_max+Delta_CL_max
CL_max_TO = C_L_max+Delta_CL_max*0.6
t_to_c = 0.3*((1-((5+M_CR**2)/(5+1.15**2))**3.5)*np.sqrt(1-M_CR**2)/M_CR**2)**(2/3)
N_eng = 2


#Formulas
def SAR(S, C_D):
    return V_CR/(0.5*rho_CR*V_CR**2*S*C_D*C_T)

def Lambda_n(AR, n_percent):
    return np.rad2deg(np.arctan(np.tan(np.deg2rad(Lambda_c4))-(4/AR)*((n_percent-25)/100*(1-taper)/(1+taper))))

def wing_fuel_volume(eta_tank, t_to_c, S_w, AR):
    return 0.9*eta_tank*t_to_c*S_w**1.5*AR**(-0.5)

def C_L_alpha(AR, Lambda_c2):
    beta = np.sqrt(1-M_CR**2)
    return 2*np.pi*AR/(2+np.sqrt((AR*beta/0.95)**2*(1+np.tan(Lambda_c2)**2/beta**2)+4))

def Oswald_eff_factor(AR):
    return ((1 + 0.12 * M_CR**6) * (1 + (0.142 + (0.005 * (1 + 1.5 * (taper - 0.6)**2)) * AR * (10 * t_to_c)**0.33) / (np.cos(np.radians(Lambda_c4))**2) + 0.1 * (3 * N_eng + 1) / (4 + AR)**0.8))**-1

def get_wing_area_and_TtoW(AR, e, CD_0, CL_max_landing, CL_to, landing_mass_fraction,m_crit):
    #Andrei formulas
    approach_speed_result = approach_speed_req(landing_mass_fraction, CL_max_landing)
    landing_field_result = landing_field_req(landing_mass_fraction, CL_max_landing)
    cruise_req_result = cruise_req_thrust_over_weight(AR, e, CD_0)
    climb_req_result = climb_req_thrust_over_weight(e, AR, CD_0)
    cs_far_119_result = climb_gradient_CS25119(AR, e, CD_0)
    cs_far_121_b_result = climb_gradient_CS25121b(AR, e, CD_0)
    cs_far_121_c_result = climb_gradient_CS25121c(AR, e, CD_0)
    cs_far_121_d_result = climb_gradient_CS25121d(AR, e, CD_0, landing_mass_fraction)
    take_off_len_result = landing_field_length_req(AR, e, CL_to)
    
    bound_low, bound_right = returns_for_optimization(
        approach_speed_result, landing_field_result, cruise_req_result, 
        climb_req_result, cs_far_119_result, cs_far_121_b_result, 
        cs_far_121_c_result, cs_far_121_d_result, take_off_len_result
    )
    
    S = m_crit * 9.81 / bound_right
    T_to_W = bound_low
    return S, T_to_W

def cruise_lift(S,m_crit,m_f):
    return 1.1/(0.5*rho_CR*V_CR**2)*0.5*((m_crit-m_f[0]*0.05)*g/S+(m_crit-(1-0.118)*m_f[0])*g/S)

def drag_coeff(AR, e, CD_0, C_L):
    return CD_0+C_L**2/(np.pi*AR*e)

def lift_to_drag(C_L, C_D):
    return C_L/C_D

def wingspan(S, AR):
    return np.sqrt(S*AR)
def V2_val(C_L_maxTO, S, TtoW):
    #return np.sqrt(TtoW*MTOM_crit_case*g/(1.225*drag_coeff(AR, e, C_D0, C_L_maxTO)*S))
    return  82.19
def rollrate(b, V2):
    return 0.0454/0.4347*np.pi/9*V2/b*2


#Variables
AR = 10.2
S = 363.53
e = Oswald_eff_factor(AR)
b = wingspan(S, AR)
Lambda_c2 = Lambda_n(AR, 50)
Lambda_LE = Lambda_n(AR, 0)
Lambda_TE = Lambda_n(AR, 100)
CL_des = cruise_lift(S,259403.4,m_f)
C_D = drag_coeff(AR, e, C_D0, CL_des)
V2 = V2_val(C_D, S, T_to_W)
P = rollrate(b, V2)
#Requirements
reqs=np.array([130.1, 0.71, 338.7, 16.4, 60/11, 0]) #V, e, S, L/D, P, SAR
current=[]

#file_path_CL_alpha = "C:\\Users\\Equipo\\Downloads\\Polar_Graph_cl_alpha.csv"
#C_l_v_alpha = pd.read_csv(file_path_CL_alpha).to_numpy()
#C_l_v_alpha = np.array([C_l_v_alpha[:,0], C_l_v_alpha[:,1]])


# plt.plot(C_l_v_alpha[0], C_l_v_alpha[1])
# plt.show()

# plt.plot(CL_v_CD[0], CL_v_CD[1])
# plt.show()

# plt.plot(CM_v_alpha[0], CM_v_alpha[1])
# plt.show()

'''
for i in range(len(CL_v_CD[0])-1):
    if (CL_v_CD[1,i]-0.857)*(CL_v_CD[1,i+1]-0.857)<0:
        C_d_cruise=CL_v_CD[0,i] #gives cruise airfoil drag
        break

for i in range(len(C_l_v_alpha[0])-1):
    if (C_l_v_alpha[1,i]-0.857)*(C_l_v_alpha[1,i+1]-0.857)<0:
        alpha_CR= (C_l_v_alpha[0,i]) #gives cruise angle
        break
'''
print("SAR original is "+str(SAR(S, C_D)))
"""
AR_range=np.arange(5, 15, 0.01)
vals =[]
for i in AR_range:
    e= Oswald_eff_factor(i)
    S, T_to_W = get_wing_area_and_TtoW(i, e, C_D0, CL_max_landing, CL_max_TO, 0.7)
    C_L_cruise = cruise_lift(S)
    C_D = drag_coeff(i, e, C_D0, C_L_cruise)
    L_over_D = lift_to_drag(C_L_cruise, C_D)
    b = wingspan(S, i)
    V2 = V2_val(CL_max_TO, S, T_to_W)
    P= rollrate(b, V2)*180/np.pi
    V_fuel = wing_fuel_volume(eta_tank, t_to_c, S, i)
    SAR_val = SAR(S, C_D)
    current = np.array([V_fuel, e, S, L_over_D, P, SAR_val])
    vals.append(current)
meet_reqs = []
for i in range(len(vals)):
    for j in range(len(vals[i])):
        if j==len(vals[i])-1 and vals[i][j]>reqs[j]:
            meet_reqs.append(True)
        elif vals[i][j]>reqs[j]:
            continue
        else:
            meet_reqs.append(False)
            break
vals_array = np.array(vals)
valid_AR_vals =[]
valid_SAR_vals = []
invalid_AR_vals1 =[]
invalid_SAR_vals1 = []
invalid_AR_vals2 =[]
invalid_SAR_vals2 = []
last_valid = 0
encountered_valid = False
for i in range(len(meet_reqs)):
    if meet_reqs[i]:
        valid_AR_vals.append(AR_range[i])
        valid_SAR_vals.append(vals[i][-1])
        encountered_valid = True
    else:
        if encountered_valid:
            invalid_AR_vals2.append(AR_range[i])
            invalid_SAR_vals2.append(vals[i][-1])
        else:
            invalid_AR_vals1.append(AR_range[i])
            invalid_SAR_vals1.append(vals[i][-1])
print("The maximum SAR is "+ str(max(valid_SAR_vals))+" for a corresponding AR of "+str(valid_AR_vals[valid_SAR_vals.index(max(valid_SAR_vals))]))


plt.plot(invalid_AR_vals1, invalid_SAR_vals1, "r", label = "Do not meet\nrequirements")
plt.plot(invalid_AR_vals2, invalid_SAR_vals2, "r")
plt.plot(valid_AR_vals, valid_SAR_vals, "k", label = "Meet requirements")
plt.plot([AR_range[288],AR_range[288]], [0, 200], "r--")
plt.plot([AR_range[583],AR_range[583]], [0, 200], "r--")
plt.xlabel(r"Aspect Ratio", fontsize = font_size_full_scr)
plt.ylabel(r"SAR $\left[m\ kg^{-1}\right]$", fontsize= font_size_full_scr)
plt.xticks(fontsize = font_size_full_scr)
plt.yticks(fontsize = font_size_full_scr)
plt.legend(fontsize = font_size_full_scr)
plt.xlim(5,15)
plt.ylim(0, 200)
plt.grid()
plt.savefig("Optimization_plot.pdf", format = "pdf")
plt.close()
"""
AR=10.82
e= Oswald_eff_factor(AR)
S, T_to_W = get_wing_area_and_TtoW(AR, e, C_D0, CL_max_landing, CL_max_TO, 0.7,259403.4)
C_L_cruise = cruise_lift(S,MTOM_crit_case,m_f)
C_D = drag_coeff(AR, e, C_D0, C_L_cruise)
L_over_D = lift_to_drag(C_L_cruise, C_D)
b = wingspan(S, AR)
V2 = V2_val(CL_max_TO, S, T_to_W)
P= rollrate(b, V2)*180/np.pi
V_fuel = wing_fuel_volume(eta_tank, t_to_c, S, AR)
SAR_val = SAR(S, C_D)
print("Oswald Efficiency Factor (e):", e)
print("Wing Area (S) [m²]:", S)
print("Thrust-to-Weight Ratio (T/W):", T_to_W)
print("Cruise Lift Coefficient (C_L_cruise):", C_L_cruise)
print("Drag Coefficient (C_D):", C_D)
print("Lift-to-Drag Ratio (L/D):", L_over_D)
print("Wingspan (b) [m]:", b)
print("V2 Speed [m/s]:", V2)
print("Roll Rate (P) [deg/s]:", P)
print("Wing Fuel Volume (V_fuel) [m³]:", V_fuel)
print("Specific Air Range (SAR) [m/kg]:", SAR_val)
