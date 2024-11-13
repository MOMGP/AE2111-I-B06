from optimization import *
from Matching_diagram_Andrei import *

#Constant drags(not affected by class 2)
C_D0_fuselage=0.003857389
C_D0_wave_drag=0.001904762
C_D0_fuselage_upsweep=0.000538676
#Friction coefficients
Cf_wing=0.002128259
Cf_HTail=0.002219171
Cf_VTail=0.002360944
Cf_nacelle=0.001992739

#Influence factors
IF_wing=1
IF_HTail=1.045
IF_VTail=1.045
IF_nacelle=1

#Form factors
FF_wing=1.337009317
FF_HTail=1.28841584
FF_VTail=1.272896358
FF_nacelle=1.587529976
def C_D0_wing_calc(S,Cf,FF,IF,S_ref):
    C_D0_wing=1.07*2*S*Cf*FF*IF/S_ref
    return C_D0_wing

def C_D0_HTail_calc(S,Cf,FF,IF,S_ref):
    C_D0_HTail=Cf*S*FF*IF/S_ref
    return C_D0_HTail

def C_DO_VTail_calc(S,Cf,FF,IF,S_ref):
    C_DO_VTail=Cf*S*FF*IF/S_ref
    return C_DO_VTail

def C_D0_Nacelle_calc(S,Cf,FF,IF,S_ref):
    C_D0_Nacelle=Cf*S*FF*IF/S_ref
    return C_D0_Nacelle