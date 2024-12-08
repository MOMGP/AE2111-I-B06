import numpy as np
from matplotlib import pyplot as plt
from get_points import get_points
from geometry_WP4_2 import get_geom_from_points, get_stringer_geom_norm, get_points_along_spanwise, get_mass, moments_of_inertia, plot3d_geom, Lambda_n, scaled_chord
from twist_distribution_over_the_wing import twist_angle 

best_comb = np.load("Work package 4\\Combinations\\0.2\\0.35\\0.65\\min_W_15568.0.npy")
snd_best = np.load("Work package 4\\Combinations\\0.2\\0.4\\0.7\\min_W_17797.0.npy")
# print(snd_best)

skip = 15
bending_x_vals = np.arange(0,26.785,0.01*skip)
E = 72400000000
b = 53.57 #m
dihedral = np.deg2rad(4.75) #deg

def angle_of_rotation(root_geom, root_str, end_third_spar, truncated, case): # does the integral of M(x) / E * I(x) which is the slope 
    bending = np.load("Work package 4\\Bending_vals\\"+case+"_crit.npy")[::skip]*100*0.5
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
    bending = np.load("Work package 4\\Bending_vals\\"+case+"_crit.npy")[::skip]*1000
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


t_tb = 0.022
t_sides = 0.01
truncated = True
# thd_end = b*0.175
n=20
cases = ["n", "rho", "V", "CL"]
i=0.2 #CAELAN - set to .3
j=0.4 #CAELAN - set to .5
k=0.7
thd_end = 0.175*b
combinations = []
x_y_y = get_points(i, j, k, 1)
root_geom = get_geom_from_points(x_y_y, [t_sides, t_tb, t_sides, t_tb, t_tb, t_sides, t_tb])
root_stringer = get_stringer_geom_norm(root_geom, n)
tip_twists = []
tip_deflections = []
for c in cases:
    tip_deflections.append(np.abs(deflection(c,root_geom, root_stringer, thd_end, truncated)[-1]))
for c in cases:
    tip_twists.append(np.abs(twist_angle(c, root_geom, root_stringer, thd_end, truncated)[-1]))
crit_tip_twist = max(tip_twists)
tip_deflection = max(tip_deflections)
print(tip_deflection, crit_tip_twist)
# for n in np.arange(2, 8,2):

#CAELAN - comment this and un-tab
# for t_tb in possible_t_tb_vals*10**(-3):
#     for t_sides in np.arange(0.4, 2, 0.1)*10**(-3):
#CALEAN - end here

# # for thd_end in np.arange(0.05*b, b/4, 0.025*b):
# # for i in np.arange(0.25, 0.5, 0.05):
# #     for j in np.arange(0.25, 0.65, 0.05):
# #         if (j<=i):
# #             combinations.append([i,j, 0, 0, 0])
# #             continue

#CAELAN - un-tab
# x_y_y = get_points(i, j, k, 1)
# possible_t_tb_vals = np.arange(0.008, 0.028, 0.004)
# # for t_tb in possible_t_tb_vals:
# #     for t_sides in np.arange(0.006, 0.02, 0.002):
# for n in np.arange(20, 44, 4):
#     for thd_end in np.arange(0.1*b, b/4, 0.025*b):
#         for t_tb in np.arange(0.01, 0.026,  0.002):
#             for t_sides in np.arange(0.008, 0.024, 0.002):
#                 x_y_y = get_points(i, j, k, 1)
#                 root_geom = get_geom_from_points(x_y_y, [t_sides, t_tb, t_sides, t_tb, t_tb, t_sides, t_tb])
#                 root_stringer = get_stringer_geom_norm(root_geom, n)
#                 tip_twists = []
#                 tip_deflections = []
#                 for c in cases:
#                     tip_deflections.append(np.abs(deflection(c,root_geom, root_stringer, thd_end, truncated)[-1]))
#                 for c in cases:
#                     tip_twists.append(np.abs(twist_angle(c, root_geom, root_stringer, thd_end, truncated)[-1]))
#                 crit_tip_twist = max(tip_twists)
#                 tip_deflection = max(tip_deflections)
#                 print(tip_deflection, crit_tip_twist)
#                 if (tip_deflection<b*0.15 and crit_tip_twist<np.deg2rad(10)): #CAELAN - remove 1.5 from BOTH
#                     mass= get_mass(get_points_along_spanwise(root_geom, root_stringer, 0, thd_end, truncated)[0], get_points_along_spanwise(root_geom, root_stringer, thd_end, thd_end, truncated)[0], get_points_along_spanwise(root_geom, root_stringer, b/2, thd_end, truncated)[0], thd_end, get_points_along_spanwise(root_geom, root_stringer, 0, thd_end, truncated)[1])
#                     combinations.append([t_tb, t_sides, thd_end, n,  mass, tip_deflection, crit_tip_twist])
#                     # print([i, j, mass, tip_deflection, crit_tip_twist])
#                 else:
#                     combinations.append([t_tb,t_sides, thd_end, n, np.inf, tip_deflection, crit_tip_twist])
#     print("finished with "+str(n)+"number of combinations is "+str(len(combinations)))

# combinations = np.array(combinations)

# np.save("Work Package 4\\Testing\\testing_t.npy", combinations)
combinations = np.load("Work Package 4\\Testing\\testing_t.npy")
mass_column = combinations[:, 4]
index_of_lowest_mass = np.argmin(mass_column)

# print("Critical twist (rad) is "+str(np.radians(10)))
# print("Critical bending is "+str(0.15*b))
print(combinations[index_of_lowest_mass])
# x_y_y = get_points(i, j, j+0.2, 1)
# geom_root = get_geom_from_points (x_y_y, [t_sides, t_tb, t_sides, t_tb, t_tb, t_sides, t_tb])
# stringers = get_stringer_geom_norm(geom_root, 2)
#plot3d_geom(geom_root, stringers, b*0.1, truncated=False, plot_with_airfoil=True, full_wing=True, plot_stringers = True)

# TESTING for i,j

# fig = plt.figure(figsize = (12,10))
# ax = plt.axes(projection='3d')

# X, Y = combinations[:,0].reshape(-1, np.arange(0.25, 0.55, 0.05).size), combinations[:,1].reshape(-1, np.arange(0.25, 0.55, 0.05).size)
# ax.scatter(X, Y, combinations[:,2].reshape(-1, np.arange(0.25, 0.55, 0.05).size))
# plt.xlabel("i")
# plt.ylabel("j")
# plt.title("pos v. mass")
# plt.show()
# ax = plt.axes(projection='3d')
# ax.scatter(X,Y, combinations[:,3].reshape(-1, np.arange(0.25, 0.55, 0.05).size))
# plt.xlabel("i")
# plt.ylabel("j")
# plt.title("pos v. def")
# plt.show()

# ax = plt.axes(projection='3d')
# ax.scatter(X,Y, combinations[:,4].reshape(-1, np.arange(0.25, 0.55, 0.05).size))
# plt.xlabel("i")
# plt.ylabel("j")
# plt.title("pos v. twist")
# plt.show()



# #TESTING - for THD end
# plt.plot(combinations[:,0], combinations[:,1])
# plt.title("Mass for varying thd_end")
# plt.show()

# plt.plot(combinations[:,0], combinations[:,2])
# plt.title("Def for varying thd_end")
# plt.show()

# plt.plot(combinations[:,0], combinations[:,3])
# plt.title("Twist for varying thd_end")
# plt.show()


#TESTING - for n
# plt.plot(combinations[:,0], combinations[:,1])
# plt.title("Mass for varying n")
# plt.show()

# plt.plot(combinations[:,0], combinations[:,2])
# plt.title("Def for varying n")
# plt.show()

# plt.plot(combinations[:,0], combinations[:,3])
# plt.title("Twist for varying n")
# plt.show()


# TESTING for t
# fig = plt.figure(figsize = (12,10))
# ax = plt.axes(projection='3d')

# X, Y = combinations[:,0].reshape(-1, possible_t_tb_vals.size), combinations[:,1].reshape(-1, possible_t_tb_vals.size)
# ax.scatter(X, Y, combinations[:,2].reshape(-1, possible_t_tb_vals.size))
# plt.xlabel("t_tb")
# plt.ylabel("t_sides")
# plt.title("thicknesses v. mass")
# plt.show()
# ax = plt.axes(projection='3d')
# ax.scatter(X,Y, combinations[:,3].reshape(-1, possible_t_tb_vals.size))
# plt.xlabel("t_tb")
# plt.ylabel("t_sides")
# plt.title("thicknesses v. def")
# plt.show()

# ax = plt.axes(projection='3d')
# ax.scatter(X,Y, combinations[:,4].reshape(-1, possible_t_tb_vals.size))
# plt.xlabel("t_tb")
# plt.ylabel("t_sides")
# plt.title("thicknesses v. twist")
# plt.show()

"""combination = np.array([mass,
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
                                            ]) """
"""[ 1.86851201e+03  1.71066055e-01 -1.11641910e+00  4.00000000e-01
   6.00000000e-01  6.50000000e-01  2.00000000e-03  5.00000000e-04
   2.00000000e+01  0.00000000e+00  2.67850000e+00  6.99950000e+04]"""




# x_y_y = get_points(i, j, k, 1)
# root_geom = get_geom_from_points(x_y_y, [t_sides, t_tb, t_sides, t_tb, t_tb, t_sides, t_tb]) #todo - change this to account for top, bot, and sides
# root_stringer = get_stringer_geom_norm(root_geom, n)
# tip_twists = []
# tip_deflections = []
# for c in cases:
#     tip_deflections.append(np.abs(deflection(c,root_geom, root_stringer, thd_end, truncated)[-1]))
# for c in cases:
#     tip_twists.append(np.abs(twist_angle(c, root_geom, root_stringer, thd_end, truncated)[-1]))
# crit_tip_twist = max(tip_twists)
# tip_deflection = max(tip_deflections)
# if (tip_deflection<b*0.15 and crit_tip_twist<np.deg2rad(10)):
#     mass= get_mass(get_points_along_spanwise(root_geom, root_stringer, 0, thd_end, truncated)[0], get_points_along_spanwise(root_geom, root_stringer, thd_end, thd_end, truncated)[0], get_points_along_spanwise(root_geom, root_stringer, b/2, thd_end, truncated)[0], thd_end, get_points_along_spanwise(root_geom, root_stringer, 0, thd_end, truncated)[1])
#     combinations.append([i, j, mass, tip_deflection, crit_tip_twist])
# else:
#     combinations.append([i,j, 0, 0, 0])