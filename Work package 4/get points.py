airfoil_geometry = []
wingbox_chords = []


def get_points(length_1, length_2, length_3):
    airfoil_geometry = [[0.5, 1, 2],[0.51, 1.1, 1.9],[0.52, 1.2, 1.8]] #get airfoil geometry

    for i in range(len(airfoil_geometry)):  #input values for length 1
        j = airfoil_geometry[i][0]
        if length_1 == j:
            wingbox_chords.append(airfoil_geometry[j][1], airfoil_geometry[j][2])

    for i in range(len(airfoil_geometry)):  #input values for length 2
        j = airfoil_geometry[i][0]
        if length_2 == j:
            wingbox_chords.append(airfoil_geometry[j][1], airfoil_geometry[j][2])

    for i in range(len(airfoil_geometry)):  # input values for length 3
        j = airfoil_geometry[i][0]
        if length_3 == j:
            wingbox_chords.append(airfoil_geometry[j][1], airfoil_geometry[j][2])

    return wingbox_chords()