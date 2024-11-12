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

def cg(MAC):
    Xj = 0.4 * MAC
    MjXj = 0.109 * Xj

    X_LEMAC = MiXi+MAC*(0.4*(Mj/Mi)-xc_OEW * (1+Mj/Mi))
    X_TEMAC = X_LEMAC + MAC
    X_cg_OEW = X_LEMAC + MAC*xc_OEW

    OEW_X = X_cg_OEW
    Payload_X = 24.275
    Fuel_X = X_LEMAC + Xj

    cg_aft = max((OEW_M*OEW_X + Payload_M*Payload_X)/(OEW_M+Payload_M), (OEW_M*OEW_X + Payload_M*Payload_X + Fuel_M*Fuel_X)/(OEW_M+Payload_M+Fuel_M), (OEW_M*OEW_X + Fuel_M*Fuel_X)/(OEW_M+Fuel_M))
    cg_for = min((OEW_M*OEW_X + Payload_M*Payload_X)/(OEW_M+Payload_M), (OEW_M*OEW_X + Payload_M*Payload_X + Fuel_M*Fuel_X)/(OEW_M+Payload_M+Fuel_M), (OEW_M*OEW_X + Fuel_M*Fuel_X)/(OEW_M+Fuel_M))
    return cg_aft, cg_for

def S_h(MAC, S_wing, cg_aft):
    S_h = (0.982*MAC*S_wing)/(52.5-cg_aft)
    return S_h

def B_h(S_h):
    B_h = math.sqrt(4.27*S_h)
    return B_h

def Cr_h(S_h, B_h):
    Cr_h = (2*S_h)/(B_h*(1.35))
    return Cr_h

def LE_sweep_h(Cr_h, B_h):
    LE_sweep_h = math.degrees(math.atan(math.tan(math.radians(32.18))-Cr_h/(2*B_h)*(0.35-1)))
    return LE_sweep_h

def X_LEMAC_h(LE_sweep_h, Cr_h):
    X_LEMAC_h = (1+2*0.35)/12*4.27*math.tan(math.radians(LE_sweep_h))*Cr_h
    return X_LEMAC_h

def L_wing_to_tail(X_LEMAC_h, Cr_h, X_LEMAC, Chord_root, YLEMAC):
    L_wing_to_tail = 52.5-X_LEMAC_h*0.25-X_LEMAC_h+0.25*Cr_h-(X_LEMAC-YLEMAC+Chord_root*0.25)
    return L_wing_to_tail

def S_v(S_wing, B_wing, cg_aft):
    S_v = (0.079*S_wing*B_wing)/(52.1-cg_aft)
    return S_v

def B_v(S_v):
    B_v = math.sqrt(1.75*S_v)
    return B_v

cg_aft, cg_for = cg(MAC)
S_h = S_h(MAC, S_wing, cg_aft)
B_h = B_h(S_h)
Cr_h = Cr_h(S_h, B_h)
LE_sweep_h = LE_sweep_h(Cr_h, B_h)
X_LEMAC_h = X_LEMAC_h(LE_sweep_h, Cr_h)
L_wing_to_tail = L_wing_to_tail(X_LEMAC_h, Cr_h, X_LEMAC, Chord_root, YLEMAC)
S_v = S_v(S_wing, B_wing, cg_aft)
B_v = B_v(S_v)
