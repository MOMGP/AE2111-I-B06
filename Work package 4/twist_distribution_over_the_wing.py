#twist distribution over the wing 
import scipy as sp
import numpy as np
from scipy import integrate
from Torsion import internal_torque_for_plotting
from get_points import get_points,get_geom_from_points
from DIS_TORSIONAL_STIFF import get_torr_stiff_list
from matplotlib import pyplot as plt
import time
cases = ["CL", "n", "rho", "V"]
#𝑑θ/𝑑𝑦 = 𝑇(𝑦)/𝐺𝐽(𝑦)
def rot(x, CL_d,rho,V,n,norm_wing_box_root, norm_stringers, end_third_spar, cond):
    int_torque_x = internal_torque_for_plotting(CL_d,rho,V,n)
    torque=int_torque_x[int(x*100)]*1000
    rot = torque/get_torr_stiff_list(norm_wing_box_root, norm_stringers, end_third_spar, cond, x)
    return rot
LC_22=[0.1896402164645566,1.225,241.96,-1.0] #CL_d,rho,V,n
LC_26=[0.640875724061477,0.31641,258.97916748380965,2.5]

def twist_angle(CL_d,rho,V,n, norm_wing_box_root, norm_stringers, end_third_spar, cond):
    G = 28*10**9
    twist = []
    
    sum = 0
    previous_val =0
    for i in range(0,26785,500):
        step = i*0.001
        temp = rot(step,CL_d,rho,V,n,norm_wing_box_root, norm_stringers, end_third_spar, cond)        
        sum+=(temp +previous_val)*0.5*0.5/1.9
        previous_val = temp
        # result =sp.integrate.quad(lambda x, CL_d, rho, V, n,norm_wing_box_root, norm_stringers, end_third_spar, cond: rot(x, CL_d, rho, V, n,norm_wing_box_root, norm_stringers, end_third_spar, cond), 0,step,args=(CL_d, rho, V, n,norm_wing_box_root, norm_stringers, end_third_spar, cond))
        # twist.append(result)
        twist.append(sum)
    return twist
span = 2*26.785
spar1_x=0.2
spar2_x=0.35
spar3_x=0.7
t_sides = 0.01
t_tb = 0.022
x_y_y = get_points(spar1_x, spar2_x, spar3_x, 1)
norm_wing_box_root = get_geom_from_points(x_y_y, [t_sides, t_tb, t_sides, t_tb, t_tb, t_sides, t_tb])
norm_stringers = [[[0,0],0]]
end_third_spar = span*0.175 
cond = True 
hws = 26.785
x = np.arange(0,hws,0.5)
#y = get_torr_stiff_list(root_geom, [[[0,0],0]], span/6, True)
#--------------------------------------------------------------------------------------------------------------
CL_d = 0.7
rho = 0.31641
V = 241.9574
n = 1


#ploting 
hws = 26.785
x = np.arange(0,hws,0.5)
y = twist_angle(LC_22[0],LC_22[1],LC_22[2],LC_22[3],norm_wing_box_root, norm_stringers, end_third_spar, cond)
y1 = twist_angle(LC_26[0],LC_26[1],LC_26[2],LC_26[3],norm_wing_box_root, norm_stringers, end_third_spar, cond)
#y2 = twist_angle('rho' ,norm_wing_box_root, norm_stringers, end_third_spar, cond)
#y3 = twist_angle('V' ,norm_wing_box_root, norm_stringers, end_third_spar, cond)
plt.plot(x,y,label='n = -1',color = 'blue')
plt.plot(x, y1, label="n= 2.5", color="red", linestyle="-")
#plt.plot(x, y2, label="n = -1", color="Black", linestyle="-")
#plt.plot(x, y3, label="Line 1", color="pink", linestyle="-")
plt.grid(True)
plt.xlabel('Spanwise location, y [m]')
plt.ylabel('Angle of twist [rad]')
plt.legend(loc='upper right')
plt.savefig("Twist_dist_D3.pdf", format="pdf")
plt.show()





    
