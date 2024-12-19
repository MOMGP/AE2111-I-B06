import numpy as np
from Geometry import get_design

b = 53.57

def design(geometry, skin_thickness, stringers, ribs):
    full_list = []
    HLD_pos = [3,8.6,9,12.5]

    #add spars and flanges
    span_pos = np.arange(0,b/2+b/(2*(ribs-1)),b/(2*(ribs-1)))
    for i in range(len(span_pos)-1):
        if span_pos[i] >= geometry[2]: #check if double cell
            #add spars
            for j in range(len(geometry[0])-1):
                full_list.append(("spar", (float(span_pos[i]), geometry[0][j]),(float(span_pos[i+1]), geometry[0][j]), geometry[1][0], int(stringers[0])))

            #add flanges
            for j in range(len(geometry[0])-2):
                full_list.append(("flange", (float(span_pos[i]), geometry[0][j]), (float(span_pos[i+1]), geometry[0][j+1]), geometry[1][1], int(stringers[j+1]), True))
                full_list.append(("flange", (float(span_pos[i]), geometry[0][j]), (float(span_pos[i+1]), geometry[0][j+1]), geometry[1][1], int(stringers[j+1]), False))

            #add skin
            if span_pos[i + 1] <= HLD_pos[0] or span_pos[i] >= HLD_pos[1] and span_pos[i + 1] <= HLD_pos[2] or span_pos[i] >= HLD_pos[3]:
                full_list.append(("skin", (float(span_pos[i]), geometry[0][1]), (float(span_pos[i + 1]), 1), float(skin_thickness), True))
                full_list.append(("skin", (float(span_pos[i]), geometry[0][1]), (float(span_pos[i + 1]), 1), float(skin_thickness), False))

        else:
            #add spars
            for j in range(len(geometry[0])):
                full_list.append(("spar", (float(span_pos[i]), geometry[0][j]),(float(span_pos[i+1]/(ribs-1)), geometry[0][j]), geometry[1][0], int(stringers[0])))

            #add flanges
            for j in range(len(geometry[0])-1):
                full_list.append(("flange", (float(span_pos[i]), geometry[0][j]),(float(span_pos[i + 1]), geometry[0][j + 1]), geometry[1][1], int(stringers[j+1]), True))
                full_list.append(("flange", (float(span_pos[i]), geometry[0][j]),(float(span_pos[i + 1]), geometry[0][j + 1]), geometry[1][1], int(stringers[j+1]), False))

            #add skin
            if span_pos[i+1] <= HLD_pos[0] or span_pos[i] >= HLD_pos[1] and span_pos[i+1] <= HLD_pos[2] or span_pos[i] >= HLD_pos[3]:
                full_list.append(("skin", (float(span_pos[i]), geometry[0][2]), (float(span_pos[i+1]), 1),float(skin_thickness) ,True))
                full_list.append(("skin", (float(span_pos[i]), geometry[0][2]), (float(span_pos[i+1]), 1),float(skin_thickness) ,False))

    return full_list

def weight():
    weight = 0
    return weight