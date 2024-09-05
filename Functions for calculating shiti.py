import math
from main import *
def isa_calc(altitude):

    #constants
    e=2.71828
    ft=0.3048
    FL=100
    b=[0,11000,20000,32000,47000,51000,71000,86000]
    a=[0,-0.0065,0,0.0010,0.0028,0,-0.0028,-0.0020]
    #sea level conditions
    g=9.80665
    R=287
    T=288.15
    rho=1.225
    p=101325

    for i in range (1,8):
        c=b[i];
        b[i]=min(altitude,b[i])
        if c!=b[i]: sf=i; break
        T_0=T
        p_0=p
        T=T+a[i]*(b[i]-b[i-1])
        if(a[i]!=0):p=p_0*(T/T_0)**(-g/(R*a[i]))
        else:p=p_0*(e**(-g/(R*T)*(b[i]-b[i-1])))


    if(a[sf]==0):
        p_0=p
        p = p_0 * (e ** (-g / (R * T) * (b[sf] - b[sf-1])))
    else:
        p_0=p
        T_0=T
        T = T + a[sf] * (b[sf]-b[sf-1])
        p = p_0 * (T / T_0) ** (-g / (R * a[sf]))

    rho = p / (R * T)
    return p,T,rho

def thrust_lapse_calc(T,p,M):
        Temp=T*(1+(gamma-1)/2*(M**2))
        Pressure=p(1+(gamma-1)/2*(M**2))**(gamma/(gamma-1))
        sigma_t=Pressure/p_0
        theta_t=Temp/T_0
        if Bypass_ratio<5:
            if theta_t<=theta_t_break:
                lapse_rate=sigma_t
            else:
                lapse_rate=sigma_t*(1-2.1*(theta_t-theta_t_break)/theta_t)
        if Bypass_ratio>=5 and Bypass_ratio<15:
            if theta_t<=theta_t_break:
                lapse_rate=sigma_t*(1-(0.43+0.014*Bypass_ratio)*math.sqrt(M))
            else:
                lapse_rate=sigma_t*(1-(0.43+0.014*Bypass_ratio)*math.sqrt(M)-3*(theta_t-theta_t_break)/(1.5+M))
        return lapse_rate
