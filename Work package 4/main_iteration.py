import numpy as np
import math
import time
import scipy
from numba import cuda, float32
from matplotlib import pyplot as plt
from scipy.integrate import dblquad, quad, quad_vec
from get_points import get_points, get_geom_from_points
from geometry_WP4_2 import centroid, moments_of_inertia, get_stringer_geom_norm,scaled_chord,scaled_length,centroid_distance, get_points_along_spanwise, get_mass
from Aero_loading_XFLR5 import normal_force_for_integrating
#from bending import moment_at_full_position
#from Shear_force import shear_force_full_values
from twist_distribution_over_the_wing import twist_angle
import os

C_r = 7.63 #m
taper = 0.3 
b = 53.57 #m
sigma_y = 450000000 # Pa
sigma_ult = 485000000 #Pa
G = 28000000000 #Pa
E = 72400000000
rho = 2780 # kg/m3
stringer_area = 0.0002
possible_t = np.arange(2, 24, 2)*0.001
possible_t_side = np.arange(0.5, 20, 0.5)*0.001
possible_x = np.arange(0.2, 0.7, 0.05)
pos_string_num = np.arange(2, 6, 2) #based on Jeroen resource
graphs = False

#Cross section plotting below
airfoil_xyy = np.load("Airfoil_geom.npy")
x_y_y = np.array(get_points(0.35, 0.55, 0.65, 1)) #EDIT THIS
geom = get_geom_from_points(x_y_y, [0.01 for i in range(7)])
stringers_for_plotting = get_stringer_geom_norm(geom, 6) #EDIT THIS (for caelan <3)
stringers_x=[]
stringers_y =[]
for i in range(len(stringers_for_plotting)):
    stringers_x.append(stringers_for_plotting[i][0][0])
    stringers_y.append(stringers_for_plotting[i][0][1])
x_vals = []
y_vals = []
for i in geom:
    x_vals.append(i[0][0])
    x_vals.append(i[1][0])
    y_vals.append(i[0][1])
    y_vals.append(i[1][1])
plt.figure(figsize=(20,4))
plt.plot(x_vals, y_vals)
plt.scatter(stringers_x, stringers_y, s=100, cmap="o")
plt.xlabel("x/c")
plt.ylabel("t/c")
plt.xlim(-0.05, 1.05)
plt.plot(airfoil_xyy[:,0], airfoil_xyy[:,1])
plt.plot(airfoil_xyy[:,0], airfoil_xyy[:,2])
plt.savefig("Cross_section_min_W.pdf", format="pdf")


def bending_stress(bending_moment, y_max, I_xx):
    sigma = (bending_moment * y_max)/I_xx
    return(sigma)     
# Only works for 2 sparred wingbox
# shear = []
# bending = []
# n_points = 50 # decrease number for higher amount of points to integrate over (a.k.a. increase accuracy
# with open("shear.txt", 'r') as f:
#     count = 0
#     for i in f:
#         line = i.strip('\n')
#         count += 1
#         if count % n_points == 0 or count == 1:
#             shear.append(line)
# with open("Work package 4\\bending.txt", 'r') as f:
#     count = 0
#     for i in f:
#         line = i.strip('\n')
#         count += 1
#         if count % n_points == 0 or count == 1:
#             bending.append(line)         
# shear = np.array(shear) #todo - load array made in bending
skip = 15
bending_x_vals = np.arange(0,26.785,0.01*skip)
#I used interpolation functions to simplify the integrations + dont know any other method to form a function out of the data
start = time.time()
# def shear_force(x): #Just defines the interpolation function for the shear force. I
#     return -1000*shear[int(x/b*2*len(shear)+0.02)]
#     # result, _ = quad(load_distribution, 0, x, limit=200)  # Shear force is the integral of load
# def bending_moments(x): #Integrates the interpolated function of shear, and adds the interpolation functions of bending provided by bending.py
#     # result, _ = quad(shear_force,x,b/2)
#     return 1000*bending[int(x/b*2*len(bending)+0.02)]
def angle_of_rotation(root_geom, root_str, end_third_spar, truncated, case): # does the integral of M(x) / E * I(x) which is the slope 
    bending = np.load("Work package 4\\Bending_vals\\"+case+"_crit.npy")[::skip]
    angle_of_rotation = []
    sum = 0
    diff = bending_x_vals[2]-bending_x_vals[1]
    for i in range(bending.size):
        if (i==0):
            angle_of_rotation.append(sum)
            geom, stringers = get_points_along_spanwise(root_geom, root_str, 0, end_third_spar, truncated)
            prev_val = 1000*bending[i]/(E * moments_of_inertia(geom, stringers)[0])
        else:
            geom, stringers  =  get_points_along_spanwise(root_geom, root_str, bending_x_vals[i], end_third_spar, truncated)
            current_val = 1000*bending[i]/(E * moments_of_inertia(geom, stringers)[0])
            sum+=(current_val+prev_val)*0.5*diff
            angle_of_rotation.append(sum)
            prev_val = current_val
    return angle_of_rotation

def deflection(case, root_geom, root_str, end_third_spar, truncated): # integrates slope to obtain deflection
    bending = np.load("Work package 4\\Bending_vals\\"+case+"_crit.npy")[::skip]
    AOR = angle_of_rotation(root_geom, root_str, end_third_spar, truncated, case)
    deflection =[]
    sum = 0
    diff = bending_x_vals[2]-bending_x_vals[1]
    for i in range(bending.size):
        if (i==0):
            deflection.append(sum)
        else:
            sum+=(AOR[i]+AOR[i-1])*0.5*diff
            deflection.append(sum)
    return deflection

# pwease dont remove mario >///<
# print('Maximum deflection is:', abs(deflection(b/2)), "m; which is", abs(deflection(b/2)/(b/2) * 100), "% of the wingspan")
# end = time.time()
# print("integration took", end - start, "seconds")
# print('Interpolated over', len(x_vals), 'intervals of x')






# Following lines are not necessary for iterating the design, this is purely for analytical purposes
# if graphs == True:
#     x_vals = np.linspace(0,b/2,100)
#     shear_vals = [shear_force(x) for x in x_vals]
#     moment_vals = [bending_moments(x) for x in x_vals]
#     angle_of_rotation_vals = [angle_of_rotation(x) for x in x_vals]
#     deflection_vals = [deflection(x) for x in x_vals] 
#     I_xx = [I_xx(I_xx_r,x) for x in x_vals]
#     end = time.time()
#     print(end - start, "seconds")


#     plt.figure(figsize=(12, 6))

#     # I_xx vals
#     plt.subplot(5, 1, 5)
#     plt.plot(x_vals, I_xx, label="I_xx", color="yellow")
#     plt.ylabel("I_xx")
#     plt.grid(alpha=0.3)
#     plt.legend()

#     # Shear force
#     plt.subplot(5, 1, 1)
#     plt.plot(x_vals, shear_vals, label="Shear Force V(x)", color="green")
#     plt.ylabel("Shear Force (N)")
#     plt.grid(alpha=0.3)
#     plt.legend()

#     # Bending moment
#     plt.subplot(5, 1, 2)
#     plt.plot(x_vals, moment_vals, label="Bending Moment M(x)", color="red")
#     plt.ylabel("Moment (Nm)")
#     plt.grid(alpha=0.3)
#     plt.legend()

#     plt.subplot(5, 1, 3) # bending moment can be added later aswell if needed :3
#     plt.plot(x_vals, angle_of_rotation_vals, label="Slope dv/dx", color="purple")
#     plt.xlabel("Position along the beam (m)")
#     plt.ylabel("Slope (m)")
#     # plt.grid(alpha=0.3)
#     plt.legend()



#     # Deflection
#     plt.subplot(5, 1, 4) #bending moment can be added later aswell if needed :3
#     plt.plot(x_vals, deflection_vals, label="Deflection w(x)", color="purple")
#     plt.xlabel("Position along the beam (m)")
#     plt.ylabel("Deflection (m)")
#     plt.grid(alpha=0.3)
#     plt.legend()
#     plt.show()

def bending_stress(bending_moment, y_max, I_xx):
    sigma = (bending_moment * y_max)/I_xx
    return(sigma) 

def scaled_length(length,chord):
    scaled_length = length*chord/C_r
    return(scaled_length)


#output = [mass, twist_ang, deflection, spar_pts, thicknesses, num_stringer, truncated, end_third]
#inputs = spar_pos, thicknesses, num_stringer, end_third, truncated
pos_combs = []
cases = ["n", "rho", "V", "CL"]
#TEST 
# x_y_y = get_points(0.2, 0.4, 0.6, 1)
# root_geom = get_geom_from_points(x_y_y, [0.001 for i in range(7)])
# root_stringer = get_stringer_geom_norm(root_geom, 4)
# tip_twists = []
# tip_deflections = []
# plt.plot(bending_x_vals, deflection("CL",root_geom, root_stringer, b/6, False))
# plt.title("x v. def")
# plt.show()
# plt.plot(np.arange(0,b/2,0.5), twist_angle("CL", root_geom, root_stringer, b/6, False))
# plt.title("x v. twist")
# plt.show()
pos_steps = 0.05
count = 0
skip_t_sides = False
skip_str = False
skip_thd_end = False
orig_time = time.time()
for i in np.arange(0.2, 0.50, pos_steps): # iterating front spar pos (skipping every 0.05 space)
    for j in np.arange(i+pos_steps, 0.65, pos_steps): #iterating second spar pos (skipping every 0.05 space)
        for k in np.arange(j+pos_steps, 0.75, pos_steps): # iterating third spar pos (skipping every 0.05 space)
            print("currently in k", count, f"time = {time.time()- orig_time}")
            combinations = []
            for truncated in [False, True]:
                for t_tb in possible_t:
                    for t_sides in possible_t_side:
                        for thd_end in np.arange(0.05*b, b/4, 0.025*b): #iterating over end of third spar pos
                            for n in pos_string_num:
                                x_y_y = get_points(i, j, k, 1)
                                root_geom = get_geom_from_points(x_y_y, [t_sides, t_tb, t_sides, t_tb, t_tb, t_sides, t_tb]) #todo - change this to account for top, bot, and sides
                                root_stringer = get_stringer_geom_norm(root_geom, n)
                                tip_twists = []
                                tip_deflections = []
                                for c in cases:
                                    tip_deflections.append(np.abs(deflection(c,root_geom, root_stringer, thd_end, truncated)[-1]))
                                for c in cases:
                                    tip_twists.append(np.abs(twist_angle(c, root_geom, root_stringer, thd_end, truncated)[-1]))
                                crit_tip_twist = max(tip_twists)
                                tip_deflection = max(tip_deflections)
                                if (tip_deflection<b*0.15 and crit_tip_twist<np.deg2rad(10)):
                                    mass= get_mass(get_points_along_spanwise(root_geom, root_stringer, 0, thd_end, truncated)[0], get_points_along_spanwise(root_geom, root_stringer, thd_end, thd_end, truncated)[0], get_points_along_spanwise(root_geom, root_stringer, b/2, thd_end, truncated)[0], thd_end, get_points_along_spanwise(root_geom, root_stringer, 0, thd_end, truncated)[1])
                                    combination = np.array([mass,
                                            crit_tip_twist,
                                            tip_deflection,
                                            i, #first spar
                                            j, #snd spar
                                            k, #3rd spar
                                            t_tb,
                                            t_sides,
                                            n,
                                            truncated,
                                            thd_end,
                                            count
                                            ]) 
                                    combinations.append(combination)
                                    #If already met, more thickness, stringers, or a longer yehudi are just going to add weight
                                    skip_t_sides=True
                                    skip_str=True
                                    skip_thd_end=True
                                    
                                count+=1
                                if skip_str:
                                        skip_str=False
                                        break
                            if skip_thd_end:
                                    skip_thd_end=False
                                    break
                        if skip_t_sides:
                            skip_t_sides=False
                            break
            combinations = np.array(combinations)
            # print(combinations)
            # print(combinations.size)
            if combinations.size!=0:
                lowest_mass = np.min(combinations[:,0])
                directory = f"Work package 4\\Combinations\\{np.round(i, 2)}\\{np.round(j, 2)}\\{np.round(k, 2)}"
                file_path = os.path.join(directory, f"min_W_{np.round(lowest_mass, 0)}.npy")

                # Create directories if they don't exist
                os.makedirs(directory, exist_ok=True)

                # Save the file
                np.save(file_path, combinations)

                print(f"File saved to: {file_path}")
                                
                            
#considerations - if smth works with a particular thickness, dont increase
# if smth works with a particular num of stringers, dont increase
# if smth dont work due to twist, skip over all stringers 
#if smth work with small yehudi (thd_end), skip next
