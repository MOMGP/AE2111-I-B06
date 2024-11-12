from optimization import *
from Matching_diagram_Andrei import *
C_D0_fuselage=0.003857389
C_D0_wave_drag=0.001904762
C_D0_fuselage_upsweep=0.000538676
D_fuselage=7.35077
def C_D0_wing_calc(S,Cf,FF,IF,S_ref):
    C_D0_wing=1.07*2*S*Cf*FF*IF/S_ref
    return C_D0_wing
def C_D0_fuselage_base_area_calc(D_fuselage_base,S_ref):
    C_D0_fuselage_base_area=D_fuselage/S_ref
    return C_D0_fuselage_base_area
def C_D0_HTail_calc(S,Cf,FF,IF,S_ref):
    C_D0_HTail=Cf*S*FF*IF/S_ref
    return C_D0_HTail

def C_DO_VTail_calc(S,Cf,FF,IF,S_ref):
    C_DO_VTail=Cf*S*FF*IF/S_ref
    return C_DO_VTail
