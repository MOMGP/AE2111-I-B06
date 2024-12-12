import numpy as np
from matplotlib import pyplot as plt
from get_points import get_points
from geometry_WP4_2 import get_geom_from_points, get_stringer_geom_norm, get_points_along_spanwise, get_mass, moments_of_inertia, plot3d_geom, Lambda_n, scaled_chord, moments_of_inertia
from twist_distribution_over_the_wing import twist_angle 

best_comb = np.load("Work package 4\\Combinations\\0.4\\0.6\\0.65\\min_W_1869.0.npy")
snd_best = np.load("Work package 4\\Combinations\\0.35\\0.55\\0.6\\min_W_1884.0.npy")
# print(best_comb)
# print(np.degrees(0.171))

skip = 15
bending_x_vals = np.arange(0,26.785,0.01*skip)
E = 72400000000
b = 53.57 #m
dihedral = np.deg2rad(4.75) #deg

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


t_tb = 0.001
t_sides = 0.0007
truncated = False
thd_end = 0.1*b
n=2
cases = ["n", "rho", "V", "CL"]
i=0.35
j=0.55
k=0.7
combinations = []

#print(np.load("Work package 4\\Combinations\\0.2\\0.3\\0.4\\min_W_4203.0.npy"))
# for n in np.arange(2, 8,2):


# for t_tb in np.arange(0.8, 2, 0.1)*10**(-3):
#     for t_sides in np.arange(0.4, 2, 0.1)*10**(-3):


# # for thd_end in np.arange(0.05*b, b/4, 0.025*b):
# # for i in np.arange(0.25, 0.5, 0.05):
# #     for j in np.arange(0.25, 0.65, 0.05):
# #         if (j<=i):
# #             combinations.append([i,j, 0, 0, 0])
# #             continue

#         x_y_y = get_points(i, j, j+0.15, 1)
#         root_geom = get_geom_from_points(x_y_y, [t_sides, t_tb, t_sides, t_tb, t_tb, t_sides, t_tb]) #todo - change this to account for top, bot, and sides
#         root_stringer = get_stringer_geom_norm(root_geom, n)
#         tip_twists = []
#         tip_deflections = []
#         for c in cases:
#             tip_deflections.append(np.abs(deflection(c,root_geom, root_stringer, thd_end, truncated)[-1]))
#         for c in cases:
#             tip_twists.append(np.abs(twist_angle(c, root_geom, root_stringer, thd_end, truncated)[-1]))
#         crit_tip_twist = max(tip_twists)
#         tip_deflection = max(tip_deflections)
#         if (tip_deflection<b*0.15/1.5 and crit_tip_twist<np.deg2rad(10)/1.5):
#             mass= get_mass(get_points_along_spanwise(root_geom, root_stringer, 0, thd_end, truncated)[0], get_points_along_spanwise(root_geom, root_stringer, thd_end, thd_end, truncated)[0], get_points_along_spanwise(root_geom, root_stringer, b/2, thd_end, truncated)[0], thd_end, get_points_along_spanwise(root_geom, root_stringer, 0, thd_end, truncated)[1])
#             combinations.append([t_tb, t_sides, mass, tip_deflection, crit_tip_twist])
#             # print([i, j, mass, tip_deflection, crit_tip_twist])
#         else:
#             combinations.append([t_tb,t_sides, np.inf, np.inf, np.inf])
# combinations = np.array(combinations)

# np.save("Work Package 4\\Testing\\testing_t.npy", combinations)
combinations = np.load("Work Package 4\\Testing\\testing_t.npy")
mass_column = combinations[:, 2]
index_of_lowest_mass = np.argmin(mass_column)

print("Critical twist (rad) is "+str(np.radians(10)))
print("Critical bending is "+str(0.15*b))
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
fig = plt.figure(figsize = (12,10))
ax = plt.axes(projection='3d')

X, Y = combinations[:,0].reshape(-1, np.arange(0.8, 2, 0.1).size), combinations[:,1].reshape(-1, np.arange(0.8, 2, 0.1).size)
ax.scatter(X, Y, combinations[:,2].reshape(-1, np.arange(0.8, 2, 0.1).size))
plt.xlabel("t_tb")
plt.ylabel("t_sides")
plt.title("thicknesses v. mass")
plt.show()
ax = plt.axes(projection='3d')
ax.scatter(X,Y, combinations[:,3].reshape(-1, np.arange(0.8, 2, 0.1).size))
plt.xlabel("t_tb")
plt.ylabel("t_sides")
plt.title("thicknesses v. def")
plt.show()

ax = plt.axes(projection='3d')
ax.scatter(X,Y, combinations[:,4].reshape(-1, np.arange(0.8, 2, 0.1).size))
plt.xlabel("t_tb")
plt.ylabel("t_sides")
plt.title("thicknesses v. twist")
plt.show()

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


def plot_def_in_aircraft(geom_root, stringers, end_third_spar, truncated=False, plot_with_airfoil=False, full_wing=False, plot_stringers = True):
    points_root, stringers_root = get_points_along_spanwise(geom_root, stringers, 0, end_third_spar, trunctated=truncated)
    points_end_third_spar, _ = get_points_along_spanwise(geom_root, stringers, end_third_spar, end_third_spar, trunctated=truncated)
    # print(points_end_third_spar)
    points_tip, stringers_tip = get_points_along_spanwise(geom_root, stringers, b/2, end_third_spar, trunctated=truncated)
    lambda_LE = np.deg2rad(Lambda_n(0))

    # x_vals = []
    # y_vals = []
    # for i in points_end_third_spar:
    #     x_vals.append(i[0][0])
    #     x_vals.append(i[1][0])
    #     y_vals.append(i[0][1])
    #     y_vals.append(i[1][1])
    # plt.plot(x_vals, y_vals)
    # plt.show()
    
    front_spar_x = np.array([[points_root[0][0][0], points_root[0][1][0]], [points_tip[0][0][0] + np.sin(lambda_LE) * b / 2, points_tip[0][1][0] + np.sin(lambda_LE) * b / 2]])
    front_spar_y = np.array([[0, 0], [b/2, b/2]])
    front_spar_z = np.array([[points_root[0][0][1], points_root[0][1][1]],[points_tip[0][0][1] + np.sin(np.deg2rad(dihedral)) * b / 2, points_tip[0][1][1] + np.sin(np.deg2rad(dihedral)) * b / 2]])
    
    snd_spar_x = np.array([[points_root[2][0][0], points_root[2][1][0]], [points_tip[2][0][0] + np.sin(lambda_LE) * b / 2, points_tip[2][1][0] + np.sin(lambda_LE) * b / 2]])
    snd_spar_y = np.array([[0, 0], [b/2, b/2]])
    snd_spar_z = np.array([[points_root[2][0][1], points_root[2][1][1]], [points_tip[2][0][1] + np.sin(np.deg2rad(dihedral)) * b / 2, points_tip[2][1][1] + np.sin(np.deg2rad(dihedral)) * b / 2]])


    thd_spar_x = np.array([[points_root[5][0][0], points_root[5][1][0]], [points_end_third_spar[5][0][0] + np.sin(lambda_LE) * end_third_spar, points_end_third_spar[5][1][0] + np.sin(lambda_LE) * end_third_spar]])
    thd_spar_y = np.array([[0, 0], [end_third_spar, end_third_spar]])
    thd_spar_z = np.array([[points_root[5][0][1], points_root[5][1][1]],[points_end_third_spar[5][0][1] + np.sin(np.deg2rad(dihedral)) * end_third_spar, points_end_third_spar[5][1][1] + np.sin(np.deg2rad(dihedral)) * end_third_spar]])
    # Ensure the axes are properly set for 3D plotting
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    if plot_with_airfoil:
        x_y_y=np.load("Airfoil_geom.npy")
        root_airfoil = x_y_y*scaled_chord(0)
        tip_airfoil = x_y_y*scaled_chord(b/2)

        tip_airfoil[:,0] += np.sin(lambda_LE) * b / 2
        tip_airfoil[:,1] += np.sin(np.deg2rad(dihedral)) * b / 2
        tip_airfoil[:,2] += np.sin(np.deg2rad(dihedral)) * b / 2


        all_x = np.tile(np.vstack((root_airfoil[:,0], tip_airfoil[:,0])).T, [2,1])
        all_y = np.array([0, b/2])
        # print(all_x.shape)
        # print(all_y.shape)
    
        all_Z_upper = np.vstack((root_airfoil[:, 1], tip_airfoil[:, 1])).T
        all_Z_lower = np.vstack((root_airfoil[:, 2], tip_airfoil[:, 2])).T
        # print(all_Z_lower.shape)
        # Tile the Z values to match the meshgrid shape
        Z_upper = np.tile(all_Z_upper, (len(all_y), 1))
        Z_lower = np.tile(all_Z_lower, (len(all_y), 1))

        # Plot root airfoil
        ax.plot_surface(all_x, all_y, Z_upper, color='red', alpha=0.2)
        ax.plot_surface(all_x, all_y, Z_lower, color='red', alpha=0.2)
        if full_wing:
            ax.plot_surface(all_x, -all_y, Z_upper, color='red', alpha=0.2)
            ax.plot_surface(all_x, -all_y, Z_lower, color='red', alpha=0.2)

    
    # Plot surfaces for front and second spar
    ax.plot_surface(front_spar_x, front_spar_y, front_spar_z)
    ax.plot_surface(snd_spar_x, snd_spar_y, snd_spar_z)
    ax.plot_surface(thd_spar_x, thd_spar_y, thd_spar_z)

    #plot stringers
    if plot_stringers:
        for i in range(len(stringers_root)):
            ax.plot([stringers_root[i][0][0], stringers_tip[i][0][0]+np.sin(lambda_LE) * b / 2], [0, b/2], [stringers_root[i][0][1], stringers_tip[i][0][1]+np.sin(np.deg2rad(dihedral)) * b / 2], linewidth=0.3)
    
    if full_wing:
        ax.plot_surface(front_spar_x, -front_spar_y, front_spar_z)
        ax.plot_surface(snd_spar_x, -snd_spar_y, snd_spar_z)
        ax.plot_surface(thd_spar_x, -thd_spar_y, thd_spar_z)
        if plot_stringers:
            for i in range(len(stringers_root)):
                ax.plot([stringers_root[i][0][0], stringers_tip[i][0][0]+np.sin(lambda_LE) * b / 2], [0, -b/2], [stringers_root[i][0][1], stringers_tip[i][0][1]+np.sin(np.deg2rad(dihedral)) * b / 2], linewidth=0.3)
        ax.set_xlim3d(-25, 25)
        ax.set_ylim3d(-25, 25)
        ax.set_zlim3d(-25, 25)

    else:
        ax.set_xlim3d(0, 25)
        ax.set_ylim3d(0, 25)
        ax.set_zlim3d(0, 25)
        
    plt.show()

    x_y_y = get_points(i, j, j+0.15, 1)
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
        combinations.append([i, j, mass, tip_deflection, crit_tip_twist])
    else:
        combinations.append([i,j, 0, 0, 0])
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


# t_tb = 0.014
# t_sides = 0.008
# truncated = True
# # thd_end = b*0.175
# n=32
# cases = ["n", "rho", "V", "CL"]
# i=0.2 #CAELAN - set to .3
# j=0.35 #CAELAN - set to .5
# k=0.65
# thd_end = 0.175*b
# combinations = []
# x_y_y = get_points(i, j, k, 1)
# root_geom = get_geom_from_points(x_y_y, [t_sides, t_tb, t_sides, t_tb, t_tb, t_sides, t_tb])
# root_stringer = get_stringer_geom_norm(root_geom, n)
# tip_twists = []
# tip_deflections = []
# for c in cases:
#     tip_deflections.append(np.abs(deflection(c,root_geom, root_stringer, thd_end, truncated)[-1]))
# for c in cases:
#     tip_twists.append(np.abs(twist_angle(c, root_geom, root_stringer, thd_end, truncated)[-1]))
# crit_tip_twist = max(tip_twists)
# tip_deflection = max(tip_deflections)
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
# for i in np.arange(0.2, 0.65, 0.05):
#     for j in np.arange(0.25, 0.7, 0.05):
# for t_tb in np.arange(10, 32, 2)*0.001:
#     for t_sides in np.arange(2, 22, 2)*0.001:
#         # if (i>=j):
#         #     combinations.append([i, j, i, j,  np.inf, np.inf, np.inf])
#         #     continue
#         x_y_y = get_points(i, j, k, 1)
#         root_geom = get_geom_from_points(x_y_y, [t_sides, t_tb, t_sides, t_tb, t_tb, t_sides, t_tb])
#         root_stringer = get_stringer_geom_norm(root_geom, n)
#         tip_twists = []
#         tip_deflections = []
#         for c in cases:
#             tip_deflections.append(np.abs(deflection(c,root_geom, root_stringer, thd_end, truncated)[-1]))
#         for c in cases:
#             tip_twists.append(np.abs(twist_angle(c, root_geom, root_stringer, thd_end, truncated)[-1]))
#         crit_tip_twist_abs_ind = tip_twists.index(max(tip_twists))
#         tip_deflection_abs_ind = tip_deflections.index(max(tip_deflections))
#         crit_tip_twist = twist_angle(cases[crit_tip_twist_abs_ind], root_geom, root_stringer, thd_end, truncated)[-1]
#         tip_deflection = deflection(cases[tip_deflection_abs_ind],root_geom, root_stringer, thd_end, truncated)[-1]
#         print(tip_deflection, crit_tip_twist)
#         # if (tip_deflection<b*0.15 and crit_tip_twist<np.deg2rad(10)): #CAELAN - remove 1.5 from BOTH
#         mass= get_mass(get_points_along_spanwise(root_geom, root_stringer, 0, thd_end, truncated)[0], get_points_along_spanwise(root_geom, root_stringer, thd_end, thd_end, truncated)[0], get_points_along_spanwise(root_geom, root_stringer, b/2, thd_end, truncated)[0], thd_end, get_points_along_spanwise(root_geom, root_stringer, 0, thd_end, truncated)[1])
#         combinations.append([t_tb, t_sides, i, j,  mass, -tip_deflection, np.degrees(crit_tip_twist)])
        # print([i, j, mass, tip_deflection, crit_tip_twist])
    # else:
    #     combinations.append([i,j, thd_end, n, np.inf, tip_deflection, crit_tip_twist])
# print("finished with "+str(n)+"number of combinations is "+str(len(combinations)))

# combinations = np.array(combinations)

# np.save("Work Package 4\\Testing\\testing_t.npy", combinations)
# combinations = np.load("Work Package 4\\Testing\\testing_t.npy")
# print("Critical twist (rad) is "+str(np.radians(10)))
# print("Critical bending is "+str(0.15*b))
# x_y_y = get_points(i, j, j+0.2, 1)
# geom_root = get_geom_from_points (x_y_y, [t_sides, t_tb, t_sides, t_tb, t_tb, t_sides, t_tb])
# stringers = get_stringer_geom_norm(geom_root, 2)
#plot3d_geom(geom_root, stringers, b*0.1, truncated=False, plot_with_airfoil=True, full_wing=True, plot_stringers = True)






# TESTING for i,j
# Assuming 'combinations' is defined with appropriate data

# Reshape the data
# X = combinations[:, 0].reshape(-1, np.arange(0.2, 0.65, 0.05).size)
# Y = combinations[:, 1].reshape(-1, np.arange(0.2, 0.65, 0.05).size)
# Z_mass = combinations[:, 4].reshape(-1, np.arange(0.2, 0.65, 0.05).size)

# fig, ax = plt.subplots(figsize=(5, 5))

# # Create 2D colormap
# c = ax.contourf(X, Y, Z_mass, cmap='viridis')
# ax.set_xlabel("$x_1$ [x/c]")
# ax.set_ylabel("$x_2$ [x/c]")

# # Add colorbar
# fig.colorbar(c, label="Mass [kg]")

# plt.savefig("i-j mass sensitivity.pdf", format="pdf")
# plt.close()


# Z_deflection = combinations[:, 5].reshape(-1, np.arange(0.2, 0.65, 0.05).size)

# fig, ax = plt.subplots(figsize=(5, 5))

# # Create 2D colormap
# c = ax.contourf(X, Y, Z_deflection, cmap='viridis')
# ax.set_xlabel("$x_1$ [x/c]")
# ax.set_ylabel("$x_2$ [x/c]")

# # Add colorbar
# fig.colorbar(c, label="Tip Deflection, v(b/2) [m]")

# plt.savefig("i-j def sensitivity.pdf", format="pdf")
# plt.close()



# Z_twist = combinations[:, 6].reshape(-1, np.arange(0.2, 0.65, 0.05).size)

# fig, ax = plt.subplots(figsize=(5, 5))

# # Create 2D colormap
# c = ax.contourf(X, Y, Z_twist, cmap='viridis')
# ax.set_xlabel("$x_1$ [x/c]")
# ax.set_ylabel("$x_2$ [x/c]")

# # Add colorbar
# fig.colorbar(c, [label="Tip Twist, $\\theta(b/2)$ [deg]"])

# plt.savefig("i-j twist sensitivity.pdf", format="pdf")
# plt.close()





#Testing k
# plt.figure(figsize = (6,5))
# plt.plot(combinations[:,0], combinations[:,4]) #weight
# plt.xlabel("$x_3$ [x/c]")
# plt.ylabel("Mass [kg]")
# plt.savefig("k mass sensitivity.pdf", format="pdf")
# plt.close()

# plt.figure(figsize = (5,5))
# plt.plot(combinations[:,0], combinations[:,5]) #weight
# plt.xlabel("$x_3$ [x/c]")
# plt.ylabel("Tip Deflection, v(b/2) [m]")
# plt.savefig("k def sensitivity.pdf", format="pdf")
# plt.close()
# plt.figure(figsize = (5,5))
# plt.plot(combinations[:,0], combinations[:,6]) #weight
# plt.xlabel("$x_3$ [x/c]")
# plt.ylabel("Tip Twist, $\\theta(b/2)$ [deg]")
# plt.savefig("k twist sensitivity.pdf", format="pdf")
# plt.close()





# TESTING for n, thd_end
# Reshape the data
# grid_shape = (np.arange(20, 44, 4).size, np.arange(0.05*b, b/4+0.025*b, 0.025*b).size)
# X = combinations[:, 0].reshape(grid_shape)
# Y = combinations[:, 1].reshape(grid_shape)

# # Mass plot
# Z_mass = combinations[:, 4].reshape(grid_shape).T
# fig, ax = plt.subplots(figsize=(5, 5))
# c = ax.imshow(Z_mass, cmap='viridis', origin='lower', 
#               extent=[X.min(), X.max(), Y.min(), Y.max()], aspect='auto')
# ax.set_xlabel("n"), ax.set_ylabel("$y_{3,end}$"), fig.colorbar(c, label="Mass [kg]")
# plt.savefig("n-3e mass sensitivity.pdf", format="pdf"), plt.close()

# # Deflection plot
# Z_deflection = combinations[:, 5].reshape(grid_shape).T
# fig, ax = plt.subplots(figsize=(5, 5))
# c = ax.imshow(Z_deflection, cmap='viridis', origin='lower', 
#               extent=[X.min(), X.max(), Y.min(), Y.max()], aspect='auto')
# ax.set_xlabel("n"), ax.set_ylabel("$y_{3,end}$"), fig.colorbar(c, label="Tip Deflection, v(b/2) [m]")
# plt.savefig("n-3e def sensitivity.pdf", format="pdf"), plt.close()

# # Twist plot
# Z_twist = combinations[:, 6].reshape(grid_shape).T
# fig, ax = plt.subplots(figsize=(5, 5))
# c = ax.imshow(Z_twist, cmap='viridis', origin='lower', 
#               extent=[X.min(), X.max(), Y.min(), Y.max()], aspect='auto')
# ax.set_xlabel("n"), ax.set_ylabel("$y_{3,end}$"), fig.colorbar(c, label="Tip Twist, $\\theta(b/2)$ [deg]")
# plt.savefig("n-3e twist sensitivity.pdf", format="pdf"), plt.close()


#Testing thicknesses
# Reshape the data
# X = combinations[:, 0].reshape(-1, np.arange(2, 22,2).size)*1000
# Y = combinations[:, 1].reshape(-1, np.arange(2, 22,2).size)*1000
# Z_mass = combinations[:, 4].reshape(-1, np.arange(2, 22,2).size)

# fig, ax = plt.subplots(figsize=(5, 5))

# # Create 2D colormap
# c = ax.contourf(X, Y, Z_mass, cmap='viridis')
# ax.set_xlabel("$t_{fl}$ [mm]")
# ax.set_ylabel("$t_{spar}$ [mm]")

# # Add colorbar
# fig.colorbar(c, label="Mass [kg]")

# plt.savefig("t mass sensitivity.pdf", format="pdf")
# plt.close()


# Z_deflection = combinations[:, 5].reshape(-1, np.arange(2, 22,2).size)

# fig, ax = plt.subplots(figsize=(5, 5))

# # Create 2D colormap
# c = ax.contourf(X, Y, Z_deflection, cmap='viridis')
# ax.set_xlabel("$t_{fl}$ [mm]")
# ax.set_ylabel("$t_{spar}$ [mm]")

# # Add colorbar
# fig.colorbar(c, label="Tip Deflection, v(b/2) [m]")

# plt.savefig("t def sensitivity.pdf", format="pdf")
# plt.close()



# Z_twist = combinations[:, 6].reshape(-1, np.arange(2, 22,2).size)

# fig, ax = plt.subplots(figsize=(5, 5))

# # Create 2D colormap
# c = ax.contourf(X, Y, Z_twist, cmap='viridis')
# ax.set_xlabel("$t_{fl}$ [mm]")
# ax.set_ylabel("$t_{spar}$ [mm]")

# # Add colorbar
# fig.colorbar(c, label="Tip Twist, $\\theta(b/2)$ [deg]")

# plt.savefig("t twist sensitivity.pdf", format="pdf")
# plt.close()


#Testing truncation


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
# n=32
# truncated = truncated
# root_geom = get_geom_from_points(x_y_y, [t_sides, t_tb, t_sides, t_tb, t_tb, t_sides, t_tb])
# root_stringer = get_stringer_geom_norm(root_geom, n)
# tip_twists = []
# tip_deflections = []
# for c in cases:
#     tip_deflections.append(np.abs(deflection(c,root_geom, root_stringer, thd_end, truncated)[-1]))
# for c in cases:
#     tip_twists.append(np.abs(twist_angle(c, root_geom, root_stringer, thd_end, truncated)[-1]))
# crit_tip_twist_abs_ind = tip_twists.index(max(tip_twists))
# tip_deflection_abs_ind = tip_deflections.index(max(tip_deflections))
# crit_tip_twist_1 = twist_angle(cases[crit_tip_twist_abs_ind], root_geom, root_stringer, thd_end, truncated)[-1]
# tip_deflection_1 = deflection(cases[tip_deflection_abs_ind],root_geom, root_stringer, thd_end, truncated)[-1]
# mass1= get_mass(get_points_along_spanwise(root_geom, root_stringer, 0, thd_end, truncated)[0], get_points_along_spanwise(root_geom, root_stringer, thd_end, thd_end, truncated)[0], get_points_along_spanwise(root_geom, root_stringer, b/2, thd_end, truncated)[0], thd_end, get_points_along_spanwise(root_geom, root_stringer, 0, thd_end, truncated)[1])

# x_y_y = get_points(i, j, k, 1)
# n=32
# root_geom = get_geom_from_points(x_y_y, [t_sides, t_tb, t_sides, t_tb, t_tb, t_sides, t_tb])
# root_stringer = get_stringer_geom_norm(root_geom, n)
# tip_twists = []
# tip_deflections = []
# for c in cases:
#     tip_deflections.append(np.abs(deflection(c,root_geom, root_stringer, thd_end, truncated)[-1]))
# for c in cases:
#     tip_twists.append(np.abs(twist_angle(c, root_geom, root_stringer, thd_end, truncated)[-1]))
# crit_tip_twist_abs_ind = tip_twists.index(max(tip_twists))
# tip_deflection_abs_ind = tip_deflections.index(max(tip_deflections))
# crit_tip_twist_2 = twist_angle(cases[crit_tip_twist_abs_ind], root_geom, root_stringer, thd_end, truncated)[-1]
# tip_deflection_2 = deflection(cases[tip_deflection_abs_ind],root_geom, root_stringer, thd_end, truncated)[-1]
# mass2= get_mass(get_points_along_spanwise(root_geom, root_stringer, 0, thd_end, truncated)[0], get_points_along_spanwise(root_geom, root_stringer, thd_end, thd_end, truncated)[0], get_points_along_spanwise(root_geom, root_stringer, b/2, thd_end, truncated)[0], thd_end, get_points_along_spanwise(root_geom, root_stringer, 0, thd_end, truncated)[1])
# print("derivative in twist is "+str(np.abs(np.degrees(crit_tip_twist_2-crit_tip_twist_1)/2)))
# print("derivative of def is "+str(np.abs((tip_deflection_2-tip_deflection_1)/2)))
# print("derivative of mass is "+str(np.abs((mass2-mass1)/2)))
