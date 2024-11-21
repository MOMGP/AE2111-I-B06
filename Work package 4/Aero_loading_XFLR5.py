import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy as sp

#MISC functions
def q_calc(rho,V):
    q=1/2*rho*V**2
    return q
C_L0_original = 0.576972
C_L10_original= 1.394834
C_D0=0.022 #FROM WP3
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

#DISTRIBUTION FOR 0 AOA

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
    CM0=CmAir_quarter_0(y)
    CM10=CmAir_quarter_10(y)
    Cm_real=CM0+(CL_d-C_L0_original)/(C_L10_original-C_L0_original)*(CM10-CM0)
    return Cm_real

def AOA(CL_d):#degrees
    angle=(CL_d-C_L0_original)/(C_L10_original-C_L0_original)*10
    return angle
def Lift_distribution_for_any_load_case(CL_d,rho,V):
    L = []
    D = []
    M = []
    span_loc=[]
    q=q_calc(rho,V)
    for i in range (0,26786):
        span_location=i/1000
        CL_value=lift_distribution_any_CL(CL_d,span_location)
        lift=CL_value*q*chord_length_interpolation(span_location)
        L.append(lift)
        L.append(lift)
        span_loc.append(span_location)
        span_loc.append(-span_location)
    return L,span_loc

lift_dist_spanwise,span_loc=Lift_distribution_for_any_load_case(0.7,0.31641,241.9574)
plt.plot(span_loc,lift_dist_spanwise)
plt.xlim([0,27])
plt.show()
print(AOA(0.6))
print(lift_dist_spanwise)