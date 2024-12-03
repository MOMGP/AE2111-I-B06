import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy import integrate
import math
#MISC functions
def q_calc(rho,V):
    q=1/2*rho*V**2
    return q
C_L0_original = 0.576972
C_L10_original= 1.394834
C_D0=0.022 #FROM WP3
Ar=10.82 #FROM WP3
e=0.71 #FROM WP3

#DATA FROM 0 deg AOA

#spanwise location
span_0=np.genfromtxt('alpha=0',skip_header=40,usecols=(0),dtype=float,max_rows=19)
#chord length
chord_l0=np.genfromtxt('alpha=0',skip_header=40,usecols=(1),dtype=float,max_rows=19)
#induced angle
ai_0=np.genfromtxt('alpha=0',skip_header=40,usecols=(2),dtype=float,max_rows=19)
#lift coefficient
CL_0=np.genfromtxt('alpha=0',skip_header=40,usecols=(3),dtype=float,max_rows=19)
#induced drag
ICd_0=np.genfromtxt('alpha=0',skip_header=40,usecols=(5),dtype=float,max_rows=19)
#pitching moment about the quarter chord
CmAir_quarter_0=np.genfromtxt('alpha=0',skip_header=40,usecols=(7),dtype=float,max_rows=19)

#DATA FROM 10 deg AOA

span_10=np.genfromtxt('alpha=10',skip_header=40,usecols=(0),dtype=float,max_rows=19)
#chord length
chord_l10=np.genfromtxt('alpha=10',skip_header=40,usecols=(1),dtype=float,max_rows=19)
#induced angle
ai_10=np.genfromtxt('alpha=10',skip_header=40,usecols=(2),dtype=float,max_rows=19)
#lift coefficient
CL_10=np.genfromtxt('alpha=10',skip_header=40,usecols=(3),dtype=float,max_rows=19)
#induced drag
ICd_10=np.genfromtxt('alpha=10',skip_header=40,usecols=(5),dtype=float,max_rows=19)
#pitching moment about the quarter chord
CmAir_quarter_10=np.genfromtxt('alpha=10',skip_header=40,usecols=(7),dtype=float,max_rows=19)

#Functions for interpolating the values above
#FOR 0 AOA
lift_coeff_interpolation_0=sp.interpolate.interp1d(span_0,CL_0,kind='linear',fill_value="extrapolate")
induced_drag_interpolation_0=sp.interpolate.interp1d(span_0,ICd_0,kind='linear',fill_value="extrapolate")
picth_mon_quarter_chord_interpolation_0=sp.interpolate.interp1d(span_0,CmAir_quarter_0,kind='linear',fill_value="extrapolate")
chord_length_interpolation=sp.interpolate.interp1d(span_0,chord_l0,kind='linear',fill_value="extrapolate")

#FOR 10 AOA
lift_coeff_interpolation_10=sp.interpolate.interp1d(span_10,CL_10,kind='linear',fill_value="extrapolate")
induced_drag_interpolation_10=sp.interpolate.interp1d(span_10,ICd_10,kind='linear',fill_value="extrapolate")
picth_mon_quarter_chord_interpolation_10=sp.interpolate.interp1d(span_10,CmAir_quarter_10,kind='linear',fill_value="extrapolate")

#DISTRIBUTION FOR 0 AOA
def lift_spanwise_0(y,q): #calculate lift per unit span at y position
    lift=lift_coeff_interpolation_0(y)*chord_length_interpolation(y)*q
    return lift
def induced_drag_spanwise_0(y,q):
    Idrag=induced_drag_interpolation_0(y)*q*chord_length_interpolation(y)
    return Idrag
def moment_spanwise_0(y,q):
    moment=picth_mon_quarter_chord_interpolation_0(y)*q*chord_length_interpolation(y)**2
    return moment

#DISTRIBUTION FOR 10 AOA

def lift_spanwise_10(y,q): #calculate lift per unit span at y position
    lift=lift_coeff_interpolation_10(y)*chord_length_interpolation(y)*q
    return lift
def induced_drag_spanwise_10(y,q):
    Idrag=induced_drag_interpolation_10(y)*q*chord_length_interpolation(y)
    return Idrag
def moment_spanwise_10(y,q):
    moment=picth_mon_quarter_chord_interpolation_10(y)*q*chord_length_interpolation(y)**2
    return moment

#INTERPOLATION OF THE DISTRIBUTIONS AT DIFFERENT CL VALUES (BASED ON LOADING CASES)

def lift_distribution_any_CL(CL_d,y):
    CL0=lift_coeff_interpolation_0(y)
    CL10=lift_coeff_interpolation_10(y)
    CL_real=CL0+(CL_d-C_L0_original)/(C_L10_original-C_L0_original)*(CL10-CL0)
    return CL_real

def drag_distribution_any_CL(CL_d,y):
    CD0=induced_drag_interpolation_0(y)
    CD10=induced_drag_interpolation_10(y)
    CD_real=CD0+(CL_d-C_L0_original)/(C_L10_original-C_L0_original)*(CD10-CD0)
    return CD_real
def pitching_moment_distribution_any_CL(CL_d,y):
    CM0=picth_mon_quarter_chord_interpolation_0(y)
    CM10=picth_mon_quarter_chord_interpolation_10(y)
    Cm_real=CM0+(CL_d-C_L0_original)/(C_L10_original-C_L0_original)*(CM10-CM0)
    return Cm_real

def AOA(CL_d):#degrees
    angle=(CL_d-C_L0_original)/(C_L10_original-C_L0_original)*10
    return angle
def Lift_distribution_for_any_load_case(CL_d,rho,V):
    L = []
    span_loc=[]
    q=q_calc(rho,V)
    for i in range (0,2678):
        span_location=i/100
        CL_value=lift_distribution_any_CL(CL_d,span_location)
        lift=CL_value*q*chord_length_interpolation(span_location)
        L.append(lift)
        span_loc.append(span_location)
    return L,span_loc
def Drag_distribution_for_any_load_case(CL_d,rho,V):
    D=[]
    q=q_calc(rho,V)
    for i in range (0,2678):
        span_location=i/100
        CL_value=lift_distribution_any_CL(CL_d,span_location)
        CD_value=C_D0+CL_value**2/np.pi/e/Ar
        drag=CD_value*q*chord_length_interpolation(span_location)
        D.append(drag)
    return D
def Moment_distribution_for_any_load_case(CL_d,rho,V):
    M=[]
    q=q_calc(rho,V)
    for i in range (0,2678):
        span_location=i/100
        Cm_value=pitching_moment_distribution_any_CL(CL_d,span_location)
        moment=Cm_value*q*chord_length_interpolation(span_location)**2
        M.append(moment)
    return M


lift_dist_spanwise,span_loc=Lift_distribution_for_any_load_case(0.7,0.31641,241.9574)
drag_dist=Drag_distribution_for_any_load_case(0.7,0.31641,241.9574)
moment_dist=Moment_distribution_for_any_load_case(0.7,0.31641,241.9574)

#ploting the spanwise lift distribution
plt.plot(span_loc,lift_dist_spanwise)
plt.xlim([-27,27])
plt.gca().set_aspect(1/4000, adjustable='box')
#plt.show()
''''
#ploting the spanwise drag distribution
plt.plot(span_loc,drag_dist)
plt.xlim([-27,27])
plt.gca().set_aspect(1/400, adjustable='box')
#plt.show()

#ploting the spanwise moment distribution
plt.plot(span_loc,moment_dist)
plt.xlim([-27,27])
#plt.show()
'''
#NEW FUNCTIONS FOR INTEGRATING
def Lift_for_integrating(x,CL_d,rho,V):
    q=q_calc(rho,V)
    span_location=x
    CL_value=lift_distribution_any_CL(CL_d,span_location)
    lift=CL_value*q*chord_length_interpolation(span_location)
    return lift
def Drag_for_integrating(x,CL_d,rho,V):
    q = q_calc(rho, V)
    span_location = x
    CL_value = lift_distribution_any_CL(CL_d, span_location)
    CD_value = C_D0 + CL_value ** 2 / np.pi / e / Ar
    drag = CD_value * q * chord_length_interpolation(span_location)
    return drag
def Moment_for_integrating(x,CL_d,rho,V):
    q=q_calc(rho,V)
    span_location=x
    Cm_value=pitching_moment_distribution_any_CL(CL_d,span_location)
    moment=Cm_value*q*chord_length_interpolation(span_location)**2
    return moment

#NORMAL FORCE CALCULATION (I'M DUMB)
normal_force=[]
def normal_force_for_integrating(x,CL_d,rho,V,n):
    q=q_calc(rho,V)
    span_location=x
    CL_value = lift_distribution_any_CL(CL_d, span_location)
    CD_value=C_D0 + CL_value ** 2 / np.pi / e / Ar
    drag =n* CD_value * q * chord_length_interpolation(span_location)
    lift = n*CL_value * q * chord_length_interpolation(span_location)
    aoa=AOA(CL_d)
    aoa=math.radians(aoa)
    normal=math.cos(aoa)*lift+math.sin(aoa)*drag
    return normal
'''
total_lift,L_error=sp.integrate.quad(Lift_for_integrating,0,26.786,args=(0.7,0.31641,241.9574),limit=50, epsabs=100)
total_drag,D_error=sp.integrate.quad(Drag_for_integrating,0,26.786,args=(0.7,0.31641,241.9574),limit=50,epsabs=100)
total_moment,M_error=sp.integrate.quad(Moment_for_integrating,0,26.786,args=(0.7,0.31641,241.9574),limit=50,epsabs=100)
'''

'''''
for i in np.arange(0,26.78,0.01):
    normal_force_for_integrating(i, 0.7,0.31641,241.9574, 1)
    

#Testing shit
#print(total_lift,L_error)
plt.plot(span_loc,normal_force)
plt.xlim([-27,27])
plt.gca().set_aspect(1/4000, adjustable='box')
plt.show()
#print(total_drag,D_error)
#print(total_moment,M_error)
'''