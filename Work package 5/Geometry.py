geometry = []

def get_design(design_number):
    if design_number == 1:          #first design
        spar1 = 0.2
        spar2 = 0.35
        spar3 = 0.65
        thickness_top_bottom = 0.014
        thickness_sides = 0.008
    elif design_number == 2:        #second design
        spar1 = 0.2
        spar2 = 0.4
        spar3 = 0.7
        thickness_top_bottom = 0.016
        thickness_sides = 0.01
    elif design_number == 3:        #third design
        spar1= 0.2
        spar2 = 0.35
        spar3 = 0.7
        thickness_top_bottom = 0.022
        thickness_sides = 0.01
    geometry.append(spar1)
    geometry.append(spar2)
    geometry.append(spar3)
    geometry.append(thickness_top_bottom)
    geometry.append(thickness_sides)
    b = 53.57
    end_second_cell = 0.175*b
    geometry.append(end_second_cell)
    return geometry


# Main Design Philosophy Minimum weight Buckling Ultimate load
# Thickness top and bottom [m] 0.014 0.016 0.022
# Thickness sides [m] 0.008 0.01 0.01
# Number of stringers [-] 32 40 20
# Mass [kg] 15568 20345 22631
# Critical tip deflection [m] 8.01 6.78 7.98 1
# Critical tip twist [deg] -9.93 -9.74 -9.69 2
# Front spar [x1/c] 0.2 0.2 0.2
# Middle spar [x2/c] 0.35 0.4 0.35
# End spar [x3/c] 0.65 0.7 0.7
# End of secondary cell 0.175b 0.175b 0.175b
# Truncated Yes Yes Yes