#from bending import moment

spanwise_stress = []
I_xx = []
y_max = []

for i in range(0,26.78,0.01):
    stress = (moment[i] * y_max[i]) / I_xx[i]
    spanwise_stress.append(stress)