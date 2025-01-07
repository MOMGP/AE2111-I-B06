from Geometry import get_design, I_xx_global
from Divide_Geom import design


design_1 = get_design(1)
print(I_xx_global(design_1))

spanwise_stress = []
I_xx = []
y_max = []
moment = moment_at_full_position(CL_d,rho,V,n)
fail_stress = []

for i in range(0,26.78,0.01):
    stress = (moment[i] * y_max[i]) / I_xx[i]
    spanwise_stress.append(stress)

MOS = fail_stress / spanwise_stress
