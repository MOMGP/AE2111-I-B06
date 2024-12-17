#Kc funtion
from matplotlib import pyplot as plt
import numpy as np
# from math import e
e = 2.71828182846
x = np.arange(0,5,0.1)
y = []

for i in x: 
    if i <= 1.0: 
        value = e**(-14.2*i+8) + 4
    elif  1.0 < i and i <= 1.4: 
        value = (5/4)*i+2.75
    elif 1.4 < i and i <= 2.0:
        value = -(5/6)*i + 17/3
    elif 2.0 < i and i < 2.4:
        value = 0.5*i + 3
    elif 2.4 < i and i <= 3.0: 
        value = -(1/3)*i + 5
    elif i > 3.0: 
        value = 4 
    y.append(value)
plt.xlim(0, 5)  # Set the range for the x-axis
plt.ylim(0, 7)
plt.plot(x,y, label = 'C')
plt.legend(loc = 'lower right')
plt.grid(True, which='major', linestyle='-', linewidth=0.5, color='gray')  # Major grid
plt.minorticks_on()  # Turn on minor ticks
plt.grid(True, which='minor', linestyle=':', linewidth=0.5, color='gray')
plt.ylabel("$K_c$")
plt.xlabel('$a/b$')
plt.show()
#graphing
