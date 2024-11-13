import math
#import mass fractions
#import mac length + S,B wing

OEW_M = 0.482026757 #need import
Fuel_M = 0.411309244 #need import
Payload_M = 0.106663999 #need import

B_wing = 62.7172366 #need import
S_wing = 363 #need import
MAC = 6.3553836 #need import
YLEMAC = 7 #need import
Chord_root = 8.914109483 #need import


xc_OEW = 0.25
Mi = 0.250339286
Xi = 25.04135279
MiXi = 25.04135279
Mj = 0.205179167
def X_Lemac_calc(MAC,MiXi,Mj,Mi,xc_OEW):
    X_LEMAC = MiXi + MAC * (0.4 * (Mj / Mi) - xc_OEW * (1 + Mj / Mi))
    return X_LEMAC
def cg(MAC,MiXi,Mj,Mi,xc_OEW,OEW_M,Fuel_M,Payload_M):
    Xj = 0.4 * MAC
    MjXj = 0.109 * Xj

    X_LEMAC = X_Lemac_calc(MAC,MiXi,Mj,Mi,xc_OEW)
    X_TEMAC = X_LEMAC + MAC
    X_cg_OEW = X_LEMAC + MAC*xc_OEW

    OEW_X = X_cg_OEW
    Payload_X = 24.275
    Fuel_X = X_LEMAC + Xj

    cg_aft = max((OEW_M*OEW_X + Payload_M*Payload_X)/(OEW_M+Payload_M), (OEW_M*OEW_X + Payload_M*Payload_X + Fuel_M*Fuel_X)/(OEW_M+Payload_M+Fuel_M), (OEW_M*OEW_X + Fuel_M*Fuel_X)/(OEW_M+Fuel_M))
    cg_for = min((OEW_M*OEW_X + Payload_M*Payload_X)/(OEW_M+Payload_M), (OEW_M*OEW_X + Payload_M*Payload_X + Fuel_M*Fuel_X)/(OEW_M+Payload_M+Fuel_M), (OEW_M*OEW_X + Fuel_M*Fuel_X)/(OEW_M+Fuel_M))
    return cg_aft, cg_for

def S_h_calc(MAC, S_wing, cg_aft):
    S_h = (0.982*MAC*S_wing)/(52.5-cg_aft)
    return S_h

def B_h_calc(S_h):
    B_h = math.sqrt(4.27*S_h)
    return B_h

def Cr_h_calc(S_h, B_h):
    Cr_h = (2*S_h)/(B_h*(1.35))
    return Cr_h

def LE_sweep_h_calc(Cr_h, B_h):
    LE_sweep_h = math.degrees(math.atan(math.tan(math.radians(32.18))-Cr_h/(2*B_h)*(0.35-1)))
    return LE_sweep_h

def X_LEMAC_h_calc(LE_sweep_h, Cr_h):
    X_LEMAC_h = (1+2*0.35)/12*4.27*math.tan(math.radians(LE_sweep_h))*Cr_h
    return X_LEMAC_h

def L_wing_to_tail_calc_h(X_LEMAC, MAC):
    L_wing_to_tail_h = 52.5-MAC*0.25-X_LEMAC
    return L_wing_to_tail_h

def S_v_calc(S_wing, B_wing, cg_aft):
    S_v = (0.079*S_wing*B_wing)/(52.1-cg_aft)
    return S_v

def B_v_calc(S_v):
    B_v = math.sqrt(1.75*S_v)
    return B_v
def Cr_v_calc(S_v,B_v):
    Cr_v = 2*S_v/(B_v*1.5)
    return Cr_v
def LE_Sweep_V_calc(Cr_v,B_v):
    LE_Sweep_v = math.degrees(math.atan(math.tan(math.radians(40.41)) - Cr_v / (2 * B_v) * (0.5 - 1)))
    return LE_Sweep_v
def X_LEMAC_v_calc(LE_sweep_v,Cr_v):
    X_LEMAC_v = (1 + 2 * 0.5) / 12 * 1.75 * math.tan(math.radians(LE_sweep_v)) * Cr_v
    return X_LEMAC_v
def L_wing_to_tail_calc_v(X_LEMAC, MAC):
    L_wing_to_tail_v = 50-MAC*0.25-X_LEMAC
    return L_wing_to_tail_v
cg_aft, cg_for = cg(MAC,MiXi,Mj,Mi,xc_OEW,OEW_M,Fuel_M,Payload_M)
S_h = S_h_calc(MAC, S_wing, cg_aft)
B_h = B_h_calc(S_h)
Cr_h = Cr_h_calc(S_h, B_h)
LE_sweep_h = LE_sweep_h_calc(Cr_h, B_h)
X_LEMAC=X_Lemac_calc(MAC,MiXi,Mj,Mi,xc_OEW)
X_LEMAC_h = X_LEMAC_h_calc(LE_sweep_h, Cr_h)
L_wing_to_tail_h = L_wing_to_tail_calc_h(X_LEMAC, MAC)
S_v = S_v_calc(S_wing, B_wing, cg_aft)
B_v = B_v_calc(S_v)
Cr_v = Cr_v_calc(S_v, B_v)
LE_Sweep_V = LE_Sweep_V_calc(Cr_v, B_v)
X_LEMAC_v = X_LEMAC_v_calc(LE_Sweep_V, Cr_v)
L_wing_to_tail_v=L_wing_to_tail_calc_v(X_LEMAC,MAC)

print(L_wing_to_tail_v,L_wing_to_tail_h)