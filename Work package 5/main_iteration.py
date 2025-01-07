import numpy as np
from Geometry import get_design, I_xx, centroid_z, scaled_chord, get_extra_weight
from Ks_function import linear_interpolation_Ks
from Divide_Geom import design
from Bucklings_eqs import max_shear, crit_web_shear, crit_skin_stress, crit_stringer_stress
from DIS_TORSIONAL_STIFF import get_J_MS, translate, get_J_SS
from geometry_WP4_2 import get_points_along_spanwise, get_geom_from_points, get_mass
from get_points import get_points

b = 53.57 #m
k_v = 1.5 # ASSUMED factor of average shear compared to max.
#-> for k_v, we are using 1.5, which is the factor for rectangular cross sections (webs are rectangular)
K = 1 # End condition factor of the column, assumed to be pinned; 1. This is because we assume the stringers to be continuous.
A = 0.000200 # m^2 area of the stringer thats set by some constraint :P
E = 72400000000 # youngs modulus
kt = 0.5
v = 0.33 # poisson ratio for our material
G = 28*10**9 #Pa

def skin_area(i, j, k, y_end):
    airfoil_geometry = np.load('Airfoil_geom.npy')
    area = 0
    for x in range(len(airfoil_geometry)-1):
        if airfoil_geometry[x,0]<i:
            area+=0.5*(scaled_chord(0)+scaled_chord(b/2))*np.sqrt((airfoil_geometry[x,0]-(airfoil_geometry[x+1,0]))**2+(airfoil_geometry[x,1]-(airfoil_geometry[x+1,1]))**2)*b/2
            area+=0.5*(scaled_chord(0)+scaled_chord(b/2))*np.sqrt((airfoil_geometry[x,0]-(airfoil_geometry[x+1,0]))**2+(airfoil_geometry[x,2]-(airfoil_geometry[x+1,2]))**2)*b/2
        if  j<=airfoil_geometry[x,0]<=k:
            area+=np.sqrt((airfoil_geometry[x,0]-(airfoil_geometry[x+1,0]))**2+(airfoil_geometry[x,1]-(airfoil_geometry[x+1,1]))**2)*(b/2-y_end)*0.5*(scaled_chord(b/2)+scaled_chord(y_end))
            area+=np.sqrt((airfoil_geometry[x,0]-(airfoil_geometry[x+1,0]))**2+(airfoil_geometry[x,2]-(airfoil_geometry[x+1,2]))**2)*(b/2-y_end)*0.5*(scaled_chord(b/2)+scaled_chord(y_end))
        else:
            area+=0.5*(scaled_chord(0)+scaled_chord(b/2))*np.sqrt((airfoil_geometry[x,0]-(airfoil_geometry[x+1,0]))**2+(airfoil_geometry[x,1]-(airfoil_geometry[x+1,1]))**2)*b/2
            area+=0.5*(scaled_chord(0)+scaled_chord(b/2))*np.sqrt((airfoil_geometry[x,0]-(airfoil_geometry[x+1,0]))**2+(airfoil_geometry[x,2]-(airfoil_geometry[x+1,2]))**2)*b/2
    return area

shear_vals  = np.genfromtxt("Work package 5\\shear_n=2.5.txt")*1000*kt
int_torque_x = np.load("Work package 4\\Torsions_diff_cases\\n_crit.npy")*1000#remember to use [int(x*2)]
bending_moment = np.load("Work package 4\\Bending_vals\\n_crit.npy")*1000
def get_margin_factors(geometry, t_skin, stringers, n_ribs):
    design_breakdown = design(geometry, t_skin, stringers, n_ribs)
    x_y_y = get_points(geometry[0][0], geometry[0][1], geometry[0][2], 1)
    root_geom = get_geom_from_points(x_y_y, [geometry[1][0], geometry[1][1], geometry[1][0], geometry[1][1], geometry[1][1], geometry[1][0], geometry[1][1]])
    bay_length = b/2/(n_ribs-1)
    geom_i_j_k = geometry[0]
    geom_tsides_tb = geometry[1]
    geom_y_end = geometry[2]
    design_margins = np.zeros(2) #1: shear, 2: column buckling 
    design_margins[:]=np.inf
    for i in range(n_ribs-1):
        wing_pos = (i+1)*bay_length
        current_geom = get_points_along_spanwise(root_geom, [], wing_pos, geometry[2], True)
        values = translate(current_geom[0])
        if (wing_pos<25.1):
            current_shear = shear_vals[int(wing_pos*100)]
        else:
            current_shear = 0
        current_bending_mom  = bending_moment[int(wing_pos*100)]
        if wing_pos <geom_y_end:
            height_sp1 = values[1]
            height_sp2 = values[2]
            height_sp3 = values[3]
            _, q1, q2 = get_J_MS(values[0], values[1],values[2],values[3],values[4],G,values[5],values[6],values[7],values[8],values[9],values[10],values[11],values[12],values[13],values[14],values[15])
            q1 *=int_torque_x[int(wing_pos*2)]
            q2 *=int_torque_x[int(wing_pos*2)]
            tau_avg = current_shear/(height_sp1*geom_tsides_tb[0]+height_sp2*geom_tsides_tb[0]+height_sp3*geom_tsides_tb[0])
            sp1_margin = np.abs(crit_web_shear(geom_tsides_tb[0], b/(n_ribs-1), height_sp1)/(tau_avg*k_v+height_sp1*q1))
            sp2_margin = np.abs(crit_web_shear(geom_tsides_tb[0], b/(n_ribs-1), height_sp2)/(tau_avg*k_v+height_sp2*q2))
            sp3_margin = np.abs(crit_web_shear(geom_tsides_tb[0], b/(n_ribs-1), height_sp3)/(tau_avg*k_v+height_sp3*(q1-q2)))
            min_margin_curr = min(sp1_margin, sp2_margin, sp3_margin)
            if min_margin_curr<design_margins[0]:
                design_margins[0]=min_margin_curr
        
        else: 
            height_sp1 = values[2]
            height_sp2 = values[3]
            avg_shear = current_shear/(height_sp1*geom_tsides_tb[0]+height_sp2*geom_tsides_tb[0])
            q_torque = int_torque_x[int(wing_pos*2)]/(0.5*(height_sp1+height_sp2)*values[1])
            sp1_margin = np.abs(crit_web_shear(geom_tsides_tb[0], b/(n_ribs-1), height_sp1)/(avg_shear*k_v+height_sp1*q_torque))
            sp2_margin = np.abs(crit_web_shear(geom_tsides_tb[0], b/(n_ribs-1), height_sp2)/(avg_shear*k_v+height_sp2*q_torque))
            min_margin_curr = min(sp1_margin, sp2_margin)
            if min_margin_curr<design_margins[0]:
                design_margins[0]=min_margin_curr
    from_type = None
    for curr_part in range(len(design_breakdown)):
        if design_breakdown[curr_part][0]=="spar":
            continue
        elif design_breakdown[curr_part][0]=="flange":
            z_cent = centroid_z(design_breakdown[curr_part][0], design_breakdown[curr_part][1][1], design_breakdown[curr_part][2][1], design_breakdown[curr_part][2][0], design_breakdown[curr_part][4], geom_tsides_tb[1], design_breakdown[curr_part][5])
            I_xx_val = I_xx(design_breakdown[curr_part][0], design_breakdown[curr_part][1][1], design_breakdown[curr_part][2][1], design_breakdown[curr_part][2][0], design_breakdown[curr_part][4], geom_tsides_tb[1], design_breakdown[curr_part][5])
            area = np.abs(design_breakdown[curr_part][2][1]-design_breakdown[curr_part][1][1])*scaled_chord(wing_pos)*geom_tsides_tb[1]+design_breakdown[curr_part][4]*A
            crit_stress = crit_stringer_stress(K, bay_length, I_xx_val, area)
            sigma_max = z_cent*current_bending_mom/I_xx_val*1000
            if (sigma_max<0):
                margin_column_stress = crit_stress/np.abs(sigma_max)
                if margin_column_stress<design_margins[1]:
                    design_margins[1]=margin_column_stress
                    from_type = design_breakdown[curr_part][0]       
        else:
            z_cent = centroid_z(design_breakdown[curr_part][0], design_breakdown[curr_part][1][1], design_breakdown[curr_part][2][1], design_breakdown[curr_part][2][0], design_breakdown[curr_part][4], t_skin, design_breakdown[curr_part][5])
            I_xx_val = I_xx(design_breakdown[curr_part][0], design_breakdown[curr_part][1][1], design_breakdown[curr_part][2][1], design_breakdown[curr_part][2][0], design_breakdown[curr_part][4], t_skin, design_breakdown[curr_part][5])
            area = np.abs(design_breakdown[curr_part][2][1]-design_breakdown[curr_part][1][1])*scaled_chord(wing_pos)*geom_tsides_tb[1]+design_breakdown[curr_part][4]*A
            crit_stress = crit_stringer_stress(K, bay_length, I_xx_val, area)
            sigma_max = z_cent*current_bending_mom/I_xx_val*1000
            if (sigma_max<0):
                margin_column_stress = crit_stress/np.abs(sigma_max)
                if margin_column_stress<design_margins[1]:
                    design_margins[1]=margin_column_stress
                    from_type = design_breakdown[curr_part][0]
    return design_margins
#longer list if iteration doesnt take all day
skin_thickness = np.array([4, 3, 2.05, 1.8, 1.63, 1.4, 1.29, 1.1])*10**(-3) #these were never going to happen: 1.1, 1.02, 0.91, 0.81, 0.71

#shorter list if iteration does take all day
#skin_thickness = [2.05, 1.8, 1.4, 1.1, 0.81, 0.64, 0.51, 0.43]

stringers_skin = np.arange(22,31,1) #is per surface, so total is *3 per half wing
stringers_cell1 = np.arange(4,12,2) #is per surface, so total is *2 per half wing
stringers_cell2 = np.arange(7,14,2) #is per surface, so total is *2 per half wing
number_of_ribs = np.arange(10,18,1) #is per half wing, so total is *2

# stringers = [stringers_skin[0], stringers_cell1[0], stringers_cell2[0]]

# geometry_0 = get_design(0)
# list = design(geometry_0, skin_thickness[0], stringers, number_of_ribs[0])
# print(list)
skip_skin_stringer = False
skip_c1_stringer = False
skip_c2_stringer = False
skip_rib = False


for design_num in range(3):
    if design_num==0:
        continue
    print("the design number is "+str(design_num))
    all_designs = []
    geometry = get_design(design_num)
    tot_skin_area = skin_area(geometry[0][0], geometry[0][1], geometry[0][2], geometry[2])
    x_y_y = get_points(geometry[0][0], geometry[0][1], geometry[0][2], 1)
    thd_end = geometry[2]
    truncated = True
    root_geom = get_geom_from_points(x_y_y, [geometry[1][0], geometry[1][1], geometry[1][0], geometry[1][1], geometry[1][1], geometry[1][0], geometry[1][1]])
    for t_skin in skin_thickness:
        for n_ribs in number_of_ribs:
            if skip_rib:
                skip_rib=False
                break
            for str_skin in stringers_skin:
                if skip_skin_stringer or skip_rib:
                    skip_skin_stringer=False
                    break
                for str_c1 in stringers_cell1:
                    if skip_c1_stringer or skip_rib:
                        skip_c1_stringer=False
                        break
                    for str_c2 in stringers_cell2:
                        if skip_c2_stringer or skip_rib:
                            skip_c2_stringer=False
                            break
                        stringers_curr = [str_skin, str_c1, str_c2]
                        margin_factors = get_margin_factors(geometry, t_skin, stringers_curr, n_ribs)
                        mass_orig= get_mass(get_points_along_spanwise(root_geom, [], 0, thd_end, truncated)[0], get_points_along_spanwise(root_geom, [], thd_end, thd_end, truncated)[0], get_points_along_spanwise(root_geom, [], b/2, thd_end, truncated)[0], thd_end, [])
                        extra_mass = get_extra_weight(tot_skin_area, t_skin, (str_skin+str_c1+str_c2), n_ribs)*2
                        tot_mass = mass_orig+extra_mass
                        if (margin_factors[0]>1 and margin_factors[1]>1):
                            skip_skin_stringer=True
                            skip_c1_stringer=True
                            skip_c2_stringer=True
                            all_designs.append([tot_mass, design_num, t_skin, n_ribs, str_skin, str_c1, str_c2]) #todo - finish
                        elif margin_factors[0]<1 and margin_factors[1]>1:
                            skip_rib=True
                        elif margin_factors[0]<0.85 or margin_factors[1]<0.85:
                            skip_c1_stringer=True
                            skip_c2_stringer=True
                        print("for "+str([tot_mass, design_num, t_skin, n_ribs, str_skin, str_c1, str_c2, margin_factors[0], margin_factors[1]])+", margins are "+str(margin_factors)+" and mass is "+str(tot_mass))
        np.save(f"Work package 5//Results//Design_{design_num}.npy", all_designs)       
  

