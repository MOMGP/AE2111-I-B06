from bending import moment_at_full_position

spanwise_stress = []
I_xx = []
y_max = []
moment = moment_at_full_position(CL_d,rho,V,n)

for i in range(0,26.78,0.01):
    stress = (moment[i] * y_max[i]) / I_xx[i]
    spanwise_stress.append(stress)