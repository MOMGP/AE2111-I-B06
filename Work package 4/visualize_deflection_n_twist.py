import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
from get_points import get_points
from geometry_WP4_2 import get_geom_from_points, get_stringer_geom_norm, get_points_along_spanwise, get_mass, moments_of_inertia, plot3d_geom, Lambda_n, scaled_chord
from twist_distribution_over_the_wing import twist_angle 
from data_analysis import deflection, angle_of_rotation
b=53.57
dihedral = np.deg2rad(4.75) #rad
skip = 15
bending_x_vals = np.arange(0,26.785,0.01*skip)
twist_y_vals = np.arange(0,26785,500)*0.001
E = 72400000000
b = 53.57 #m
dihedral = np.deg2rad(4.75) #deg
cases = ["n", "rho", "V", "CL"]

#geometry definition
i=0.2
j=0.35
k=0.7
n=32
t_tb = 0.014
t_sides = 0.008
thd_end = b*0.175
truncated = True


def plot_def_in_aircraft(geom_root, stringers, end_third_spar, truncated=True, plot_with_airfoil=False, full_wing=False, plot_stringers = True, show_deflection = False, show_twist = False, view_engine = False, zoom=1):
    points_root, stringers_root = get_points_along_spanwise(geom_root, stringers, 0, end_third_spar, trunctated=truncated)
    points_end_third_spar, _ = get_points_along_spanwise(geom_root, stringers, end_third_spar, end_third_spar, trunctated=truncated)
    # print(points_end_third_spar)
    points_tip, stringers_tip = get_points_along_spanwise(geom_root, stringers, b/2, end_third_spar, trunctated=truncated)
    lambda_LE = np.deg2rad(Lambda_n(0))
    root_geom = geom_root
    root_stringer = stringers
    tip_deflections = []
    tip_twists = []
    for c in cases:
        tip_deflections.append(np.abs(deflection(c,root_geom, root_stringer, thd_end, truncated)[-1]))
        tip_twists.append(np.abs(twist_angle(c, root_geom, root_stringer, thd_end, truncated)[-1]))
    index_deflections = tip_deflections.index(max(tip_deflections))
    index_twists = tip_twists.index(max(tip_twists))
    deflections = deflection(cases[index_deflections],root_geom, root_stringer, thd_end, truncated)
    for i in range(len(deflections)):
        deflections[i] = -deflections[i]
    twists = twist_angle(cases[index_twists], root_geom, root_stringer, thd_end, truncated)
    for i in range(len(twists)):
        twists[i]=-twists[i]
    fig = plt.figure(figsize = (5, 5))
    ax = fig.add_subplot(111, projection='3d')

    front_spar_x = []
    front_spar_y = []
    front_spar_z = []

    snd_spar_x = []
    snd_spar_y = []
    snd_spar_z = []

    thd_spar_x = []
    thd_spar_y = []
    thd_spar_z = []
    
    color_vals12 = []
    color_vals3 = []
    if show_twist:
        for i in range(twist_y_vals.size):
            prev_geom, prev_str = get_points_along_spanwise(root_geom, root_stringer, twist_y_vals[i], thd_end, truncated)
            current_geom, current_str = get_points_along_spanwise(root_geom, root_stringer, twist_y_vals[i], thd_end, truncated)
            top_x_prev_1 = prev_geom[0][0][0]
            bot_x_prev_1 = prev_geom[0][1][0]
            top_x_now_1 = current_geom[0][0][0]
            bot_x_now_1 = current_geom[0][1][0]
            top_x_prev_2 = prev_geom[2][0][0]
            bot_x_prev_2 = prev_geom[2][1][0]
            top_x_now_2 = current_geom[2][0][0]
            bot_x_now_2 = current_geom[2][1][0]

            top_z_prev_1 = prev_geom[0][0][1]
            bot_z_prev_1 = prev_geom[0][1][1]
            top_z_now_1 = current_geom[0][0][1]
            bot_z_now_1 = current_geom[0][1][1]
            top_z_prev_2 = prev_geom[2][0][1]
            bot_z_prev_2 = prev_geom[2][1][1]
            top_z_now_2 = current_geom[2][0][1]
            bot_z_now_2 = current_geom[2][1][1]

            # Front Spar X
            front_spar_x.append(np.array([
                np.cos(twists[i]) * top_x_prev_1 - np.sin(twists[i]) * top_z_prev_1 + np.sin(lambda_LE) * twist_y_vals[i],
                np.cos(twists[i]) * bot_x_prev_1 - np.sin(twists[i]) * bot_z_prev_1 + np.sin(lambda_LE) * twist_y_vals[i],
            ]))

            # Front Spar Y
            front_spar_y.append(np.array([
                twist_y_vals[i],
                twist_y_vals[i],
            ]))
            color_vals12.append(twists[i])

            # Front Spar Z
            front_spar_z.append(np.array([
                np.sin(twists[i]) * top_x_prev_1 + np.cos(twists[i]) * top_z_prev_1 + np.sin(dihedral) * twist_y_vals[i],
                np.sin(twists[i]) * bot_x_prev_1 + np.cos(twists[i]) * bot_z_prev_1 + np.sin(dihedral) * twist_y_vals[i],
            ]))

            # Second Spar X
            snd_spar_x.append(np.array([
                np.cos(twists[i]) * top_x_prev_2 - np.sin(twists[i]) * top_z_prev_2 + np.sin(lambda_LE) * twist_y_vals[i],
                np.cos(twists[i]) * bot_x_prev_2 - np.sin(twists[i]) * bot_z_prev_2 + np.sin(lambda_LE) * twist_y_vals[i],
            ]))

            # Second Spar Y
            snd_spar_y.append(np.array([
                twist_y_vals[i],
                twist_y_vals[i],
            ]))

            # Second Spar Z
            snd_spar_z.append(np.array([
                np.sin(twists[i]) * top_x_prev_2 + np.cos(twists[i]) * top_z_prev_2 + np.sin(dihedral) * twist_y_vals[i],
                np.sin(twists[i]) * bot_x_prev_2 + np.cos(twists[i]) * bot_z_prev_2 + np.sin(dihedral) * twist_y_vals[i],
            ]))

            if (twist_y_vals[i]<thd_end):
                top_x_prev_3 = prev_geom[5][0][0]
                bot_x_prev_3 = prev_geom[5][1][0]
                top_x_now_3 = current_geom[5][0][0]
                bot_x_now_3 = current_geom[5][1][0]
                top_z_prev_3 = prev_geom[5][0][1]
                bot_z_prev_3 = prev_geom[5][1][1]
                top_z_now_3 = current_geom[5][0][1]
                bot_z_now_3 = current_geom[5][1][1]

                # Third Spar X
                thd_spar_x.append(np.array([
                    np.cos(twists[i]) * top_x_prev_3 - np.sin(twists[i]) * top_z_prev_3 + np.sin(lambda_LE) * twist_y_vals[i],
                    np.cos(twists[i]) * bot_x_prev_3 - np.sin(twists[i]) * bot_z_prev_3 + np.sin(lambda_LE) * twist_y_vals[i],  # Switched
                ]))
                color_vals3.append(twists[i])


                # Third Spar Y
                thd_spar_y.append(np.array([
                    twist_y_vals[i],
                    twist_y_vals[i],
                ]))

                # Third Spar Z
                thd_spar_z.append(np.array([
                    np.sin(twists[i]) * top_x_prev_3 + np.cos(twists[i]) * top_z_prev_3 + np.sin(dihedral) * twist_y_vals[i],
                    np.sin(twists[i]) * bot_x_prev_3 + np.cos(twists[i]) * bot_z_prev_3 + np.sin(dihedral) * twist_y_vals[i],  # Switched
                ]))

    elif show_deflection:
        # Front Spar X
        for i in range(bending_x_vals.size): 
            prev_geom, prev_str = get_points_along_spanwise(root_geom, root_stringer, bending_x_vals[i], thd_end, truncated)
            current_geom, current_str = get_points_along_spanwise(root_geom, root_stringer, bending_x_vals[i], thd_end, truncated)

            # === Front Spar ===
            top_x_prev_1 = prev_geom[0][0][0]
            bot_x_prev_1 = prev_geom[0][1][0]
            top_z_prev_1 = prev_geom[0][0][1]
            bot_z_prev_1 = prev_geom[0][1][1]

            # Front Spar X
            front_spar_x.append(np.array([
                top_x_prev_1 + np.sin(lambda_LE) * bending_x_vals[i],
                bot_x_prev_1 + np.sin(lambda_LE) * bending_x_vals[i],
            ]))

            # Front Spar Y
            front_spar_y.append(np.array([
                bending_x_vals[i],
                bending_x_vals[i],
            ]))
            color_vals12.append(deflections[i])
            # color_vals12.append(-np.load("Work package 4\\Bending_vals\\n_crit.npy")[::skip][i]*1000)

            # Front Spar Z (with deflections)
            front_spar_z.append(np.array([
                top_z_prev_1 + np.sin(dihedral) * bending_x_vals[i] + deflections[i],
                bot_z_prev_1 + np.sin(dihedral) * bending_x_vals[i] + deflections[i],
            ]))

            # === Second Spar ===
            top_x_prev_2 = prev_geom[2][0][0]
            bot_x_prev_2 = prev_geom[2][1][0]
            top_z_prev_2 = prev_geom[2][0][1]
            bot_z_prev_2 = prev_geom[2][1][1]

            # Second Spar X
            snd_spar_x.append(np.array([
                top_x_prev_2 + np.sin(lambda_LE) * bending_x_vals[i],
                bot_x_prev_2 + np.sin(lambda_LE) * bending_x_vals[i],
            ]))

            # Second Spar Y
            snd_spar_y.append(np.array([
                bending_x_vals[i],
                bending_x_vals[i],
            ]))

            # Second Spar Z (with deflections)
            snd_spar_z.append(np.array([
                top_z_prev_2 + np.sin(dihedral) * bending_x_vals[i] + deflections[i],
                bot_z_prev_2 + np.sin(dihedral) * bending_x_vals[i] + deflections[i],
            ]))

            # === Third Spar (if within range) ===
            if bending_x_vals[i] < thd_end:
                top_x_prev_3 = prev_geom[5][0][0]
                bot_x_prev_3 = prev_geom[5][1][0]
                top_z_prev_3 = prev_geom[5][0][1]
                bot_z_prev_3 = prev_geom[5][1][1]

                # Third Spar X
                thd_spar_x.append(np.array([
                    top_x_prev_3 + np.sin(lambda_LE) * bending_x_vals[i],
                    bot_x_prev_3 + np.sin(lambda_LE) * bending_x_vals[i],
                ]))
                color_vals3.append(deflections[i])
                # color_vals3.append(-np.load("Work package 4\\Bending_vals\\n_crit.npy")[::skip][i]*1000)

                # Third Spar Y
                thd_spar_y.append(np.array([
                    bending_x_vals[i],
                    bending_x_vals[i],
                ]))

                # Third Spar Z (with deflections)
                thd_spar_z.append(np.array([
                    top_z_prev_3 + np.sin(dihedral) * bending_x_vals[i] + deflections[i],
                    bot_z_prev_3 + np.sin(dihedral) * bending_x_vals[i] + deflections[i],
                ]))

    
    else:
        front_spar_x = np.array([[points_root[0][0][0], points_root[0][1][0]], [points_tip[0][0][0] + np.sin(lambda_LE) * b / 2, points_tip[0][1][0] + np.sin(lambda_LE) * b / 2], [points_end_third_spar[0][0][0] + np.sin(lambda_LE) * thd_end, points_end_third_spar[0][1][0] + np.sin(lambda_LE) * thd_end]])
        front_spar_y = np.array([[0, 0], [b/2, b/2]])
        front_spar_z = np.array([[points_root[0][0][1], points_root[0][1][1]],[points_tip[0][0][1] + np.sin(np.deg2rad(dihedral)) * b / 2, points_tip[0][1][1] + np.sin(np.deg2rad(dihedral)) * b / 2], [points_end_third_spar[0][0][1] + np.sin(np.deg2rad(dihedral)) * thd_end, points_end_third_spar[0][1][1] + np.sin(np.deg2rad(dihedral)) * thd_end]])
        
        snd_spar_x = np.array([[points_root[2][0][0], points_root[2][1][0]], [points_tip[2][0][0] + np.sin(lambda_LE) * b / 2, points_tip[2][1][0] + np.sin(lambda_LE) * b / 2],  [points_end_third_spar[2][0][0] + np.sin(lambda_LE) * thd_end, points_end_third_spar[2][1][0] + np.sin(lambda_LE) * thd_end]])
        snd_spar_y = np.array([[0, 0], [b/2, b/2]])
        snd_spar_z = np.array([[points_root[2][0][1], points_root[2][1][1]], [points_tip[2][0][1] + np.sin(np.deg2rad(dihedral)) * b / 2, points_tip[2][1][1] + np.sin(np.deg2rad(dihedral)) * b / 2], [points_end_third_spar[2][0][1] + np.sin(np.deg2rad(dihedral)) * thd_end, points_end_third_spar[2][1][1] + np.sin(np.deg2rad(dihedral)) *  thd_end]])
        color_vals12 = [1 , 1]

        thd_spar_x = np.array([[points_root[5][0][0], points_root[5][1][0]], [points_end_third_spar[5][0][0] + np.sin(lambda_LE) * end_third_spar, points_end_third_spar[5][1][0] + np.sin(lambda_LE) * end_third_spar]])
        thd_spar_y = np.array([[0, 0], [end_third_spar, end_third_spar]])
        thd_spar_z = np.array([[points_root[5][0][1], points_root[5][1][1]],[points_end_third_spar[5][0][1] + np.sin(np.deg2rad(dihedral)) * end_third_spar, points_end_third_spar[5][1][1] + np.sin(np.deg2rad(dihedral)) * end_third_spar]])
        color_vals3 = [1 , 1]
    # Ensure the axes are properly set for 3D plotting
        if plot_stringers:
            for i in range(len(stringers_root)):
                ax.plot([stringers_root[i][0][0], stringers_tip[i][0][0]+np.sin(lambda_LE) * b / 2], [0, b/2], [stringers_root[i][0][1], stringers_tip[i][0][1]+np.sin(np.deg2rad(dihedral)) * b / 2], linewidth=0.3, color='r')
            if full_wing:
                if plot_stringers:
                    for i in range(len(stringers_root)):
                        ax.plot([stringers_root[i][0][0], stringers_tip[i][0][0]+np.sin(lambda_LE) * b / 2], [0, -b/2], [stringers_root[i][0][1], stringers_tip[i][0][1]+np.sin(np.deg2rad(dihedral)) * b / 2], linewidth=0.3)
                
    if plot_with_airfoil:
        x_y_y=np.load("Airfoil_geom.npy")
        if show_twist:
            # Initialize 2D arrays for the top and bottom surfaces
            n_twists = twist_y_vals.size
            n_points = x_y_y.shape[0]

            # Create empty 2D arrays to store coordinates for the top and bottom surfaces
            airfoil_x_top = np.zeros((n_twists, n_points))
            airfoil_y_top = np.zeros((n_twists, n_points))
            airfoil_z_top = np.zeros((n_twists, n_points))

            airfoil_x_bottom = np.zeros((n_twists, n_points))
            airfoil_y_bottom = np.zeros((n_twists, n_points))
            airfoil_z_bottom = np.zeros((n_twists, n_points))

            for i in range(n_twists):
                # Scale the geometry for the current twist
                prev_geom = x_y_y * scaled_chord(twist_y_vals[i])

                for j in range(n_points):
                    # Top surface coordinates
                    airfoil_x_top[i, j] = prev_geom[j, 0] * np.cos(twists[i]) - prev_geom[j, 1] * np.sin(twists[i]) + np.sin(lambda_LE) * twist_y_vals[i]
                    airfoil_y_top[i, j] = twist_y_vals[i]
                    airfoil_z_top[i, j] = prev_geom[j, 0] * np.sin(twists[i]) + prev_geom[j, 1] * np.cos(twists[i]) + np.sin(dihedral) * twist_y_vals[i]

                    # Bottom surface coordinates
                    airfoil_x_bottom[i, j] = prev_geom[j, 0] * np.cos(twists[i]) - prev_geom[j, 2] * np.sin(twists[i]) + np.sin(lambda_LE) * twist_y_vals[i]
                    airfoil_y_bottom[i, j] = twist_y_vals[i]
                    airfoil_z_bottom[i, j] = prev_geom[j, 0] * np.sin(twists[i]) + prev_geom[j, 2] * np.cos(twists[i]) + np.sin(dihedral) * twist_y_vals[i]

            # Plot the top surface
            ax.plot_surface(airfoil_x_top, airfoil_y_top, airfoil_z_top, color='#0B0B0B', alpha=0.05)

            # Plot the bottom surface
            ax.plot_surface(airfoil_x_bottom, airfoil_y_bottom, airfoil_z_bottom, color='#0B0B0B', alpha=0.05)
            if full_wing:
                ax.plot_surface(airfoil_x_top, -airfoil_y_top, airfoil_z_top, color='#0B0B0B', alpha=0.05)

                # Plot the bottom surface
                ax.plot_surface(airfoil_x_bottom, -airfoil_y_bottom, airfoil_z_bottom, color='#0B0B0B', alpha=0.05)
        elif show_deflection:
            n_twists = bending_x_vals.size
            n_points = x_y_y.shape[0]

            # Create empty 2D arrays to store coordinates for the top and bottom surfaces
            airfoil_x_top = np.zeros((n_twists, n_points))
            airfoil_y_top = np.zeros((n_twists, n_points))
            airfoil_z_top = np.zeros((n_twists, n_points))

            airfoil_x_bottom = np.zeros((n_twists, n_points))
            airfoil_y_bottom = np.zeros((n_twists, n_points))
            airfoil_z_bottom = np.zeros((n_twists, n_points))

            for i in range(n_twists):
                # Scale the geometry for the current twist
                prev_geom = x_y_y * scaled_chord(bending_x_vals[i])

                for j in range(n_points):
                    # Top surface coordinates
                    airfoil_x_top[i, j] = prev_geom[j, 0] + np.sin(lambda_LE) * bending_x_vals[i]
                    airfoil_y_top[i, j] = bending_x_vals[i]
                    airfoil_z_top[i, j] = prev_geom[j, 1]  + np.sin(dihedral) * bending_x_vals[i] +  deflections[i]

                    # Bottom surface coordinates
                    airfoil_x_bottom[i, j] = prev_geom[j, 0] + np.sin(lambda_LE) * bending_x_vals[i]
                    airfoil_y_bottom[i, j] = bending_x_vals[i]
                    airfoil_z_bottom[i, j] = prev_geom[j, 2]  + np.sin(dihedral) * bending_x_vals[i] +  deflections[i]

            # Plot the top surface
            ax.plot_surface(airfoil_x_top, airfoil_y_top, airfoil_z_top, color='#0B0B0B', alpha=0.15)

            # Plot the bottom surface
            ax.plot_surface(airfoil_x_bottom, airfoil_y_bottom, airfoil_z_bottom, color='#0B0B0B', alpha=0.15)
            if full_wing:
                ax.plot_surface(airfoil_x_top, -airfoil_y_top, airfoil_z_top, color='#0B0B0B', alpha=0.15)

                # Plot the bottom surface
                ax.plot_surface(airfoil_x_bottom, -airfoil_y_bottom, airfoil_z_bottom, color='#0B0B0B', alpha=0.15)

        else:
            tip_airfoil = x_y_y*scaled_chord(b/2)
            root_airfoil = x_y_y*scaled_chord(0)
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
            ax.plot_surface(all_x, all_y, Z_upper, color='#0B0B0B', alpha=0.2)
            ax.plot_surface(all_x, all_y, Z_lower, color='#0B0B0B', alpha=0.2)
            if full_wing:
                ax.plot_surface(all_x, -all_y, Z_upper, color='#0B0B0B', alpha=0.2)
                ax.plot_surface(all_x, -all_y, Z_lower, color='#0B0B0B', alpha=0.2)

    
    # Plot surfaces for front and second spar
    front_spar_x = np.array(front_spar_x)
    front_spar_y = np.array(front_spar_y)
    front_spar_z = np.array(front_spar_z)
    snd_spar_x = np.array(snd_spar_x)
    snd_spar_y = np.array(snd_spar_y)
    snd_spar_z = np.array(snd_spar_z)
    thd_spar_x = np.array(thd_spar_x)
    thd_spar_y = np.array(thd_spar_y)
    thd_spar_z = np.array(thd_spar_z)
    

    color_vals12=np.array(color_vals12)
    color_vals3 = np.array(color_vals3)

    if show_twist:
        color_vals12=np.degrees(np.array(color_vals12))
        color_vals3 = np.degrees(np.array(color_vals3))
        color_vals12=-color_vals12
        color_vals3=-color_vals3

    norm = colors.Normalize(np.min(color_vals12), np.max(color_vals12))
    
        
    color_vals12 = np.tile(color_vals12, [2,1]).T.reshape(-1, 2)
    color_vals3 = np.tile(color_vals3, [2,1]).T.reshape(-1, 2)

    # ax.plot_surface(front_spar_x[:2], front_spar_y[:2], front_spar_z[:2], rstride=1, cstride=1, facecolors =cm.rainbow(norm(color_vals12)))
    # ax.plot_surface(snd_spar_x[:2], snd_spar_y[:2], snd_spar_z[:2], rstride=1, cstride=1, facecolors = cm.rainbow(norm(color_vals12)))
    # ax.plot_surface(np.array([[snd_spar_x[0][0], snd_spar_x[-1][0]], [thd_spar_x[0][0],  thd_spar_x[-1][0]]]), 
    #                 np.array([[0, thd_end], [0, thd_end]]),
    #                 np.array([[snd_spar_z[0][0], snd_spar_z[len(thd_spar_z)-1][0] ], [thd_spar_z[0][0], thd_spar_z[-1][0]]]), color = 'b', alpha = 0.6)
    if show_twist:
        fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cm.rainbow), ax=ax, label="Twist, $\\theta$ [deg]", shrink = 0.7)
    elif show_deflection:
        # fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cm.rainbow), ax=ax, label="Bending moment distribution, [kN/m]", shrink=0.7)
        fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cm.rainbow), ax=ax, label="Wing deflection, [m]", shrink=0.7)
    #plot stringers
    if show_twist:
        if plot_stringers:
            for i in range(len(stringers)):
                stringer_x=[]
                stringer_y=[]
                stringer_z=[]
                for j in range(twist_y_vals.size):
                    prev_str = get_points_along_spanwise(geom_root, stringers, twist_y_vals[j], end_third_spar, trunctated=truncated)[1]
                    stringer_x.append(prev_str[i][0][0]+np.sin(lambda_LE)* twist_y_vals[j])
                    stringer_y.append(twist_y_vals[j])
                    stringer_z.append(prev_str[i][0][1]+np.sin(dihedral)*twist_y_vals[j])
                ax.plot(stringer_x, stringer_y, stringer_z, cmap = cm.rainbow(norm(color_vals12[:,1])), linewidth=0.3, alpha = 0.3)
            if full_wing:
                stringer_x=[]
                stringer_y=[]
                stringer_z=[]
                for j in range(twist_y_vals.size):
                    prev_str = get_points_along_spanwise(geom_root, stringers, twist_y_vals[j], end_third_spar, trunctated=truncated)[1]
                    stringer_x.append(prev_str[i][0][0]+np.sin(lambda_LE)* twist_y_vals[j])
                    stringer_y.append(-twist_y_vals[j])
                    stringer_z.append(prev_str[i][0][1]+np.sin(dihedral)*twist_y_vals[j])
                ax.plot(stringer_x, stringer_y, stringer_z, "r", linewidth=0.3, alpha = 0.3)
    if view_engine:
        radius = 4.17/2
        length = 10
        center = (7.5-length*0.7, 9.37, 9.37*np.sin(dihedral)-radius*1.5)
        theta = np.linspace(0, 2 * np.pi, 25)  # Angle around the cylinder
        x = np.linspace(0, length, 25)          # Length of the cylinder
        theta, x = np.meshgrid(theta, x)        # Create a meshgrid for the angles and length
        curve = -(x-length/2)**2/(length**2)-(x-length/2)/length/3+1
        # Parametric equations for the cylinder pointing in the x-direction
        y = center[1] + radius * np.cos(theta)*curve
        z = center[2] + radius * np.sin(theta)*curve
        x = center[0] + x

        ax.plot_surface(x, y, z, color='k', alpha=0.4)

        theta_cap = np.linspace(0, 2 * np.pi, 25)
        theta_cap, r = np.meshgrid(theta_cap, np.linspace(0, radius, 25))  # Circular coordinates

        y_cap = center[1] + r * np.cos(theta_cap)*0.5
        z_cap = center[2] + r * np.sin(theta_cap)*0.5

        # Cap at the start of the cylinder
        x_cap_start = np.full_like(y_cap, center[0])
        ax.plot_surface(x_cap_start, y_cap, z_cap, color='k', alpha=0.4)
        ax.plot_surface(x_cap_start+length, y_cap, z_cap, color='k', alpha=0.4)
        if full_wing:
            ax.plot_surface(x, -y, z, color='k', alpha=0.4)
            ax.plot_surface(x_cap_start, -y_cap, z_cap, color='k', alpha=0.4)
            ax.plot_surface(x_cap_start+length, -y_cap, z_cap, color='k', alpha=0.4)

    if full_wing:
        ax.plot_surface(front_spar_x, -front_spar_y, front_spar_z, rstride=1, cstride=1, facecolors =cm.rainbow(norm(color_vals12)))
        ax.plot_surface(snd_spar_x, -snd_spar_y, snd_spar_z, rstride=1, cstride=1, facecolors = cm.rainbow(norm(color_vals12)))
        ax.plot_surface(thd_spar_x, -thd_spar_y, thd_spar_z, rstride=1, cstride=1, facecolors = cm.rainbow(norm(color_vals3)))
        lims = 25/zoom
        x_offset = 4
        ax.set_xlim3d(-lims+x_offset, lims+x_offset)
        ax.set_ylim3d(-lims, lims)
        ax.set_zlim3d(-lims, lims)

    else:
        x_offset = 4
        lims = 25/zoom
        ax.set_xlim3d(x_offset, lims+x_offset)
        ax.set_ylim3d(0, lims)
        ax.set_zlim3d(0, lims)
    ax.view_init(elev=90, azim=0)
    # Hide grid lines
    # ax.set_xticks([])
    # ax.set_yticks([])
    # ax.set_zticks([])
    # plt.axis('off')
    plt.tight_layout()
    plt.show()

i=0.2
j=0.35
k=0.65
x_y_y = get_points(i, j, k, 1)
root_geom = get_geom_from_points(x_y_y, [t_sides, t_tb, t_sides, t_tb, t_tb, t_sides, t_tb])
root_str = get_stringer_geom_norm(root_geom, n)
plot_def_in_aircraft(root_geom, root_str, thd_end, truncated=truncated, plot_with_airfoil=True, full_wing=True, plot_stringers=True, show_deflection=True, show_twist=False, view_engine=True, zoom=1)