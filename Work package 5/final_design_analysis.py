from matplotlib.lines import Line2D
import numpy as np
from Geometry import get_design_final, scaled_chord, I_xx_global, I_xx, centroid_z
from Divide_Geom import design
from matplotlib import pyplot as plt
from get_points import get_points
from DIS_TORSIONAL_STIFF import translate, get_J_MS
from Bucklings_eqs import crit_web_shear, crit_stringer_stress
from geometry_WP4_2 import get_geom_from_points, get_points_along_spanwise
b = 53.57  #m
k_v = 1.5 # ASSUMED factor of average shear compared to max.
#-> for k_v, we are using 1.5, which is the factor for rectangular cross sections (webs are rectangular)
K = 1 # End condition factor of the column, assumed to be pinned; 1. This is because we assume the stringers to be continuous.
A = 200*10**(-6) # m^2 area of the stringer thats set by some constraint :P
E = 72400000000 # youngs modulus
v = 0.33 # poisson ratio for our material
G = 28*10**9 #Pa
kt= 0.5
sigma_g_max = 409000000
sigma_n_max = 349200000
gamma_m_0 = 1.1
UTS = sigma_g_max*gamma_m_0
stringer_default_I_xx = 4.34*10**(-8)

# designs = np.load("Work package 5\\Results\\Design_2.npy")
# filtered_designs = designs[designs[:, 1] == 2]
# designs_ult = np.load("Work package 5\\Results\\Design_ulti.npy")
# index = np.argmin(filtered_designs, axis=0)
# print(index)
# print(filtered_designs[index[0]]*1.07)# gives 1.034, 1.036

# print(f"for ulti load {designs_ult[np.argmin(designs_ult, axis=0)[0]]}")

design_1 = get_design_final(1)
part_split = design(design_1[0], design_1[1], design_1[2], design_1[3])
shear_vals  = np.genfromtxt("Work package 5\\shear_n=2.5.txt")*1000*kt
shear_vals_neg = np.genfromtxt("Work package 5\\shear_n=-1.txt")*1000*kt
int_torque_x = np.load("Work package 4\\Torsions_diff_cases\\n_crit.npy")*1000#remember to use [int(x*2)]
bending_moment = np.load("Work package 4\\Bending_vals\\n_crit.npy")*1000
bending_moment_neg = np.load("Work package 4\\Bending_vals\\rho_crit.npy")*1000
print(bending_moment_neg)

shear_vals_25_adjusted = []
shear_vals_neg1_adjusted=[]
for i in range(len(shear_vals)):
    shear_vals_25_adjusted.append(shear_vals[i]+350000*((len(shear_vals)-i)/len(shear_vals))**2)

for i in range(len(shear_vals_neg)):
    shear_vals_neg1_adjusted.append(shear_vals_neg[i]+350000*((len(shear_vals_neg)-i)/len(shear_vals_neg)+0.1)**2)
shear_vals = shear_vals_neg1_adjusted
# shear_vals = shear_vals_25_adjusted
bending_moment = bending_moment_neg

def get_margin_factors_breakdown(geometry, t_skin, stringers, n_ribs):
    design_breakdown = design(geometry, t_skin, stringers, n_ribs)
    x_y_y = get_points(geometry[0][0], geometry[0][1], geometry[0][2], 1)
    root_geom = get_geom_from_points(x_y_y, [geometry[1][0], geometry[1][1], geometry[1][0], geometry[1][1], geometry[1][1], geometry[1][0], geometry[1][1]])
    bay_length = b/2/(n_ribs-1)
    geom_i_j_k = geometry[0]
    geom_tsides_tb = geometry[1]
    geom_y_end = geometry[2]
    design_margins = np.zeros(2) #1: shear, 2: column buckling 
    design_margins[:]=np.inf
    design_margins_shear = []


    crit_shear_stress_spanwise = []
    compressive_load_stringer = [[-1, np.inf]]
    compressive_load_flange = [[-1, np.inf]]
    compressive_load_skin = [[-1, np.inf]]
    
    tensile_load_stringer = [[-1, -np.inf]]
    tensile_load_skin = [[-1, -np.inf]]

    design_margins_spanwise_buckling = [[-1, np.inf, np.inf]]
    design_margins_spanwise_shear = [[-1, np.inf]]

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
            sp1_max_shear = np.abs(crit_web_shear(geom_tsides_tb[0], b/(n_ribs-1), height_sp1))
            sp2_max_shear = np.abs(crit_web_shear(geom_tsides_tb[0], b/(n_ribs-1), height_sp2))
            sp3_max_shear = np.abs(crit_web_shear(geom_tsides_tb[0], b/(n_ribs-1), height_sp3))
            max_shear_curr = max(sp1_max_shear, sp2_max_shear, sp3_max_shear)
            if max_shear_curr == sp1_max_shear:
                crit_shear_stress_spanwise.append([wing_pos, max_shear_curr, geom_tsides_tb[0], height_sp1])
            elif max_shear_curr == sp2_max_shear:
                crit_shear_stress_spanwise.append([wing_pos, max_shear_curr, geom_tsides_tb[0], height_sp2])
            else:
                crit_shear_stress_spanwise.append([wing_pos, max_shear_curr, geom_tsides_tb[0], height_sp2])

            sp1_margin = np.abs(crit_web_shear(geom_tsides_tb[0], b/(n_ribs-1)/2, height_sp1)/(tau_avg*k_v+height_sp1*q1))
            sp2_margin = np.abs(crit_web_shear(geom_tsides_tb[0], b/(n_ribs-1)/2, height_sp2)/(tau_avg*k_v+height_sp2*q2))
            sp3_margin = np.abs(crit_web_shear(geom_tsides_tb[0], b/(n_ribs-1)/2, height_sp3)/(tau_avg*k_v+height_sp3*(q1-q2)))
            min_margin_curr = min(sp1_margin, sp2_margin, sp3_margin)
            design_margins_shear.append([wing_pos, min_margin_curr])
            if design_margins_spanwise_shear[-1][0]!=wing_pos:
                design_margins_spanwise_shear.append([wing_pos, min_margin_curr])
            else:
                if design_margins_spanwise_shear[-1][1]<min_margin_curr:
                    design_margins_spanwise_shear[-1][1] = min_margin_curr
            if min_margin_curr<design_margins[0]:
                design_margins[0]=min_margin_curr
        
        else: 
            height_sp1 = values[2]
            height_sp2 = values[3]
            avg_shear = current_shear/(height_sp1*geom_tsides_tb[0]+height_sp2*geom_tsides_tb[0])
            q_torque = int_torque_x[int(wing_pos*2)]/(0.5*(height_sp1+height_sp2)*values[1])
            sp1_max_shear = np.abs(crit_web_shear(geom_tsides_tb[0], b/(n_ribs-1), height_sp1))
            sp2_max_shear = np.abs(crit_web_shear(geom_tsides_tb[0], b/(n_ribs-1), height_sp2))
            max_shear_curr = max(sp1_max_shear, sp2_max_shear)
            if max_shear_curr == sp1_max_shear:
                crit_shear_stress_spanwise.append([wing_pos, max_shear_curr, geom_tsides_tb[0], height_sp1])
            else:
                crit_shear_stress_spanwise.append([wing_pos, max_shear_curr, geom_tsides_tb[0], height_sp2])

            sp1_margin = np.abs(crit_web_shear(geom_tsides_tb[0], b/(n_ribs-1), height_sp1)/(tau_avg*k_v+height_sp1*q1))
            sp2_margin = np.abs(crit_web_shear(geom_tsides_tb[0], b/(n_ribs-1), height_sp2)/(tau_avg*k_v+height_sp2*q2))
            min_margin_curr = min(sp1_margin, sp2_margin)
            design_margins_shear.append([wing_pos, min_margin_curr])
            if design_margins_spanwise_shear[-1][0]!=wing_pos:
                design_margins_spanwise_shear.append([wing_pos, min_margin_curr])
            else:
                if design_margins_spanwise_shear[-1][1]>min_margin_curr:
                    design_margins_spanwise_shear[-1][1] = min_margin_curr
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
            if (sigma_max>0):
                margin_column_stress = crit_stress/np.abs(sigma_max)
                
                if margin_column_stress<design_margins[1]:
                    design_margins[1]=margin_column_stress
                    from_type = design_breakdown[curr_part][0] 
                
            elif (sigma_max)<0:
                margin_column_stress = crit_stress/np.abs(sigma_max)
                if compressive_load_flange[-1][0]!= design_breakdown[curr_part][2][0]:
                    compressive_load_flange.append([design_breakdown[curr_part][2][0], crit_stress])
                else:
                    if sigma_max< compressive_load_flange[-1][1]:
                        compressive_load_flange[-1][1] = sigma_max
                if design_margins_spanwise_buckling[-1][0]!= design_breakdown[curr_part][2][0]:
                    design_margins_spanwise_buckling.append([design_breakdown[curr_part][2][0], margin_column_stress, "flange"])
                else:
                    if crit_stress>design_margins_spanwise_buckling[-1][1]:
                        design_margins_spanwise_buckling[-1][1]= margin_column_stress
                        design_margins_spanwise_buckling[-1][2] = "flange"
        else: #skin
            z_cent = centroid_z(design_breakdown[curr_part][0], design_breakdown[curr_part][1][1], design_breakdown[curr_part][2][1], design_breakdown[curr_part][2][0], design_breakdown[curr_part][4], t_skin, design_breakdown[curr_part][5])
            # print(z_cent, design_breakdown[curr_part][2][0])
            I_xx_val = I_xx(design_breakdown[curr_part][0], design_breakdown[curr_part][1][1], design_breakdown[curr_part][2][1], design_breakdown[curr_part][2][0], design_breakdown[curr_part][4], t_skin, design_breakdown[curr_part][5])
            area = np.abs(design_breakdown[curr_part][2][1]-design_breakdown[curr_part][1][1])*scaled_chord(wing_pos)*geom_tsides_tb[1]+design_breakdown[curr_part][4]*A
            crit_stress = crit_stringer_stress(K, bay_length, I_xx_val, area)
            sigma_max = z_cent*current_bending_mom/I_xx_val*1000
            # if design_breakdown[curr_part][2][0]== 1.6740625:
                # print(f"The I_xx for the highest skin is {I_xx_val}, with an L = {1.6740625}, and an area of {area}, and a critical sigma of {sigma_max}")
            if (sigma_max>0):
                margin_column_stress = crit_stress/np.abs(sigma_max)
                if margin_column_stress<design_margins[1]:
                    design_margins[1]=margin_column_stress
                    from_type = design_breakdown[curr_part][0]
                if tensile_load_skin[-1][0]!= design_breakdown[curr_part][2][0]:
                    tensile_load_skin.append([design_breakdown[curr_part][2][0], sigma_max])
                else:
                    if sigma_max> tensile_load_skin[-1][1]:
                        tensile_load_skin[-1][1] = sigma_max
            elif (sigma_max)<0:
                if compressive_load_skin[-1][0]!= design_breakdown[curr_part][2][0]:
                    compressive_load_skin.append([design_breakdown[curr_part][2][0], sigma_max])
                else:
                    if sigma_max< compressive_load_skin[-1][1]:
                        compressive_load_skin[-1][1] = sigma_max
                if design_margins_spanwise_buckling[-1][0]!= design_breakdown[curr_part][2][0]:
                    design_margins_spanwise_buckling.append([design_breakdown[curr_part][2][0], margin_column_stress, "skin"])
                else:
                    if crit_stress>design_margins_spanwise_buckling[-1][1]:
                        design_margins_spanwise_buckling[-1][1]= margin_column_stress
                        design_margins_spanwise_buckling[-1][2]= "skin"
    buckling_margins = np.array(design_margins_spanwise_buckling)
    shear_margins = np.array(design_margins_spanwise_shear)
    print(buckling_margins)    
    print(shear_margins)
    return design_margins



# get_margin_factors_breakdown()
# print(design_1[0])
# print(design_1[2])
# print(design_1[3])
# # print(bending_moment[-30])

# print(get_margin_factors_breakdown(design_1[0], design_1[1], design_1[2], design_1[3]))
# print(crit_stringer_stress(1, 1.674, 1.493e-6, 0.001078))
# get_margin_factors_breakdown(design_1[0], design_1[1], design_1[2], design_1[3])
buckling_margins = np.load("Work package 5\\Results\\Buckling_margins.npy")
buckling_margins_neg = np.load("Work package 5\\Results\\Buckling_margins_neg.npy")
shear_margins = np.load("Work package 5\\Results\\Shear_margins.npy")
shear_margins_neg = np.load("Work package 5\\Results\\Shear_margins_neg.npy")
# print(buckling_margins)
def plot_margins(buckling_or_shear):
    design =2
    if design ==2 :
        buckling_margins = np.array([
            [-1, np.inf, 'inf'], 
            [1.4097368421052632, 1.8533318704216317, 'skin'],
            [2.8194736842105264, 1.8035606583614647, 'skin'],
            [4.229210526315789, 4.159112753535135, 'flange'],
            [5.638947368421053, 3.636314669690086, 'flange'],
            [7.048684210526316, 3.1836299809460833, 'flange'],
            [8.458421052631579, 2.793569277757491, 'flange'],
            [9.868157894736843, 2.4593516228385823, 'flange'],
            [11.277894736842105, 4.83484236321144, 'flange'],
            [12.687631578947368, 5.111000688767885, 'flange'],
            [14.097368421052632, 7.729523232077635, 'skin'],
            [15.507105263157895, 8.908430299429087, 'skin'],
            [16.916842105263157, 10.782783326082365, 'skin'],
            [18.32657894736842, 14.148288953683634, 'skin'],
            [19.736315789473686, 21.805545644737176, 'skin'],
            [21.146052631578947, 55.61941658155703, 'skin'],
            [22.55578947368421, 74.61556889933873, 'skin'],
            [23.965526315789475, 20.432351812879773, 'skin'],
            [25.375263157894736, 11.155256539460657, 'skin'],
            [26.785, 7.295179492876713, 'skin']
        ], dtype=object)

        buckling_margins_neg = np.array([
            [-1, np.inf, 'inf'],
            [1.4097368421052632, 3.247602912824416, 'flange'],
            [2.8194736842105264, 3.160388563445842, 'flange'],
            [4.229210526315789, 3.0930284632048353, 'flange'],
            [5.638947368421053, 3.0459900884720623, 'flange'],
            [7.048684210526316, 3.020191034016789, 'flange'],
            [8.458421052631579, 3.0171350337510474, 'flange'],
            [9.868157894736843, 3.0391147958300757, 'flange'],
            [11.277894736842105, 11.16745749258974, 'flange'],
            [12.687631578947368, 12.1497670028062, 'flange'],
            [14.097368421052632, 13.544483081450784, 'flange'],
            [15.507105263157895, 15.61028796344898, 'flange'],
            [16.916842105263157, 18.894726355821543, 'flange'],
            [18.32657894736842, 24.79211907525841, 'flange'],
            [19.736315789473686, 38.2099691273661, 'flange'],
            [21.146052631578947, 97.4621880639036, 'flange'],
            [26.785, 58.51145420885478, 'flange']
        ], dtype=object)

        shear_margins = np.array([
            [-1., np.inf],
            [1.40973684, 1.41613734],
            [2.81947368, 1.50321926],
            [4.22921053, 1.60378154],
            [5.63894737, 1.7170955],
            [7.04868421, 1.84669568],
            [8.45842105, 1.99501545],
            [9.86815789, 2.1679334],
            [11.27789474, 2.46544508],
            [12.68763158, 2.8212253],
            [14.09736842, 3.25057796],
            [15.50710526, 3.77399471],
            [16.91684211, 4.41936813],
            [18.32657895, 5.22537298],
            [19.73631579, 6.24676943],
            [21.14605263, 7.56297328],
            [22.55578947, 9.29238904],
            [23.96552632, 11.61735545],
            [25.37526316, 14.82964376],
            [26.785, 19.41821835]
        ])

        shear_margins_neg = np.array([
            [-1., np.inf],
            [1.40973684, 5.02300145],
            [2.81947368, 5.65957503],
            [4.22921053, 6.40954805],
            [5.63894737, 7.30887437],
            [7.04868421, 8.40834895],
            [8.45842105, 9.77666018],
            [9.86815789, 54.89022127],
            [11.27789474, 38.71051011],
            [12.68763158, 31.61689632],
            [14.09736842, 27.96468964],
            [15.50710526, 26.05625093],
            [16.91684211, 24.85156047],
            [18.32657895, 24.32864539],
            [19.73631579, 24.58113907],
            [21.14605263, 25.54122697],
            [22.55578947, 27.2507192],
            [23.96552632, 29.85665056],
            [25.37526316, 33.64310105],
            [26.785, 39.11202772]
        ])
    elif design == 3:
        buckling_margins = np.array([
            [-1, float('inf'), 'inf'], 
            [1.5755882352941177, 3.0155487347949292, 'skin'],
            [3.1511764705882355, 4.81859492513869, 'flange'],
            [4.726764705882353, 4.262452167404857, 'flange'],
            [6.302352941176471, 3.784904628313934, 'flange'],
            [7.877941176470589, 3.3777684517402315, 'flange'],
            [9.453529411764706, 3.0338400984977025, 'flange'],
            [11.029117647058824, 1.6599663435899743, 'flange'],
            [12.604705882352942, 1.707166731490386, 'flange'],
            [14.18029411764706, 2.1108888001161774, 'skin'],
            [15.755882352941178, 2.2765457556126067, 'skin'],
            [17.331470588235295, 2.5251779915448047, 'skin'],
            [18.90705882352941, 2.9106926291938002, 'skin'],
            [20.48264705882353, 3.5534287793148596, 'skin'],
            [22.058235294117647, 4.7871090571323105, 'skin'],
            [23.633823529411767, 8.007058166350452, 'skin'],
            [25.209411764705884, 36.303424337418676, 'skin'],
            [26.785, 11.577081911087918, 'skin']
        ])
        shear_margins = np.array([
            [-1, float('inf')],
            [1.57558824, 1.38614974],
            [3.15117647, 1.48347153],
            [4.72676471, 1.59793863],
            [6.30235294, 1.73494602],
            [7.87794118, 1.89770094],
            [9.45352941, 2.11407113],
            [11.02911765, 2.43522681],
            [12.60470588, 2.82647929],
            [14.18029412, 3.30849945],
            [15.75588235, 3.90989156],
            [17.33147059, 4.67106147],
            [18.90705882, 5.65044342],
            [20.48264706, 6.93486654],
            [22.05823529, 8.65749898],
            [23.63382353, 11.03036266],
            [25.20941176, 14.40655025],
            [26.785, 19.40741636]
        ])
        buckling_margins_neg = np.array([
            [-1, float('inf'), 'inf'],
            [1.5755882352941177, 5.284161466805201, 'flange'],
            [3.1511764705882355, 5.287463346365835, 'flange'],
            [4.726764705882353, 5.332993052870847, 'flange'],
            [6.302352941176471, 5.426885459071221, 'flange'],
            [7.877941176470589, 5.578068608967385, 'flange'],
            [9.453529411764706, 5.799695124485432, 'flange'],
            [11.029117647058824, 3.381976780201316, 'flange'],
            [12.604705882352942, 3.5055852262006533, 'flange'],
            [14.18029411764706, 3.6989212376443685, 'flange'],
            [15.755882352941178, 3.9892027677825395, 'flange'],
            [17.331470588235295, 4.424882306089814, 'flange'],
            [18.90705882352941, 5.10042157681983, 'flange'],
            [20.48264705882353, 6.226691419055819, 'flange'],
            [22.058235294117647, 8.388475677814888, 'flange'],
            [23.633823529411767, 14.030808965842917, 'flange'],
            [25.209411764705884, 63.614676089760295, 'flange']
        ])

        shear_margins_neg = np.array([
            [-1, float('inf')],
            [1.57558824, 5.04401638],
            [3.15117647, 5.77206958],
            [4.72676471, 6.67531561],
            [6.30235294, 7.79281668],
            [7.87794118, 9.19475048],
            [9.45352941, 34.05401088],
            [11.02911765, 27.37964154],
            [12.60470588, 24.06916055],
            [14.18029412, 22.40045351],
            [15.75588235, 21.34537353],
            [17.33147059, 20.98926701],
            [18.90705882, 21.35237862],
            [20.48264706, 22.39135192],
            [22.05823529, 24.17427282],
            [23.63382353, 26.88610845],
            [25.20941176, 30.87832422],
            [26.785, 36.78878011]
        ])
        
    if buckling_or_shear == "buckling":
        data = buckling_margins
        data_neg = buckling_margins_neg
        x = data[1:, 0].astype(float)
        y = data[1:, 1].astype(float)
        labels = data[1:, 2]

        # Plot positive margins with conditional coloring
        for i in range(1, len(x)):
            color = 'blue' if labels[i-1] == 'skin' else 'red'
            if i == 1:
                plt.plot([0, x[i-1]], [y[i-1], y[i-1]], color=color)
            plt.plot(x[i-1:i+1], y[i-1:i+1], color=color)

        # Plot negative margins with different colors
        x_neg = data_neg[1:, 0].astype(float)
        y_neg = data_neg[1:, 1].astype(float)
        labels_neg = data_neg[1:, 2]

        for i in range(1, len(x_neg)):
            color = 'cyan' if labels_neg[i-1] == 'skin' else 'magenta'
            if i == 1:
                plt.plot([0, x_neg[i-1]], [y_neg[i-1], y_neg[i-1]], color=color)
            plt.plot(x_neg[i-1:i+1], y_neg[i-1:i+1], color=color)

        # Create custom legend
        legend_elements = [
            Line2D([0], [0], color='blue', lw=2, label='Skin-Critical (n=2.5)'),
            Line2D([0], [0], color='red', lw=2, label='Flange-Critical (n=2.5)'),
            Line2D([0], [0], color='cyan', lw=2, label='Skin-Critical (n=-1.0)'),
            Line2D([0], [0], color='magenta', lw=2, label='Flange-Critical (n=-1.0)'),
            Line2D([0], [0], color='green', lw=2, label='Stringer-Critical (n=2.5)'),
            Line2D([0], [0], color='orange', lw=2, label='Stringer-Critical (n=-1.0)')
        ]

        plt.legend(handles=legend_elements, loc='best')
        plt.xlabel('y-position [m]')
        plt.ylabel('Margin of Buckling (log)')

    else:
        data = shear_margins
        data_neg = shear_margins_neg

        # Plot positive shear margins
        plt.plot([0, data[1,1]*1.15], [data[1,1], data[1,1]], color='black', label='Shear Margins (n=2.5)')
        plt.plot(data[:, 0], data[:, 1], color='black', label='Shear Margins (n=2.5)')
        
        # Plot negative shear margins
        plt.plot([0, data_neg[1,0]], [data_neg[1,1], data_neg[1,1]], color='orange', label='Shear Margins (n=-1.0)')
        plt.plot(data_neg[:, 0], data_neg[:, 1], color='orange', label='Shear Margins (n=-1.0)')
        
        plt.xlabel('y-position [m]')
        plt.ylabel('Margin of Shear')

        # Create legend
        legend_elements = [
            Line2D([0], [0], color='black', lw=2, label='Shear Margins (n=2.5)'),
            Line2D([0], [0], color='orange', lw=2, label='Shear Margins (n=-1.0)')
        ]

        plt.legend(handles=legend_elements, loc='best')

    # Show plot
    plt.xlim(0, b/2)
    plt.yscale('log')
    plt.show()
plot_margins("buckling")
# plt.plot(np.linspace(0, b/2, len(shear_vals_25_adjusted)), shear_vals_25_adjusted)
# plt.show()

