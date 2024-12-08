#twist distribution over the wing 
import scipy as sp
import numpy as np
from scipy import integrate
#from Torsion import internal_torque_at_x
from get_points import get_points,get_geom_from_points
from DIS_TORSIONAL_STIFF import get_torr_stiff_list
from matplotlib import pyplot as plt
import time
cases = ["CL", "n", "rho", "V"]
#𝑑θ/𝑑𝑦 = 𝑇(𝑦)/𝐺𝐽(𝑦)
def rot(x, case,norm_wing_box_root, norm_stringers, end_third_spar, cond):
    int_torque_x = np.load("Work package 4\\Torsions_diff_cases\\"+case+"_crit.npy")[int(x*2)]*100*0.35
    rot = int_torque_x/get_torr_stiff_list(norm_wing_box_root, norm_stringers, end_third_spar, cond, x)
    return rot


def twist_angle(case, norm_wing_box_root, norm_stringers, end_third_spar, cond):
    G = 28*10**9
    twist = []
    
    sum = 0
    previous_val =0
    for i in range(0,26785,500):
        step = i*0.001
        temp = rot(step, case,norm_wing_box_root, norm_stringers, end_third_spar, cond)        
        sum+=(previous_val+temp)*0.5*0.5
        previous_val = temp
        # result =sp.integrate.quad(lambda x, CL_d, rho, V, n,norm_wing_box_root, norm_stringers, end_third_spar, cond: rot(x, CL_d, rho, V, n,norm_wing_box_root, norm_stringers, end_third_spar, cond), 0,step,args=(CL_d, rho, V, n,norm_wing_box_root, norm_stringers, end_third_spar, cond))
        # twist.append(result)
        twist.append(sum)
    return twist
span = 2*26.785
spar1_x=0.2
spar2_x=0.5
spar3_x=0.7
x_y_y = get_points(spar1_x, spar2_x, spar3_x, 1)
norm_wing_box_root = get_geom_from_points(x_y_y, [0.0005 for i in range(7)])
norm_stringers = [[[0,0],0]]
end_third_spar = span/6 
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
# x = np.arange(0,hws,0.5)
# y = twist_angle("n",norm_wing_box_root, norm_stringers, end_third_spar, cond)
# plt.plot(x,y)
# plt.show()





    
