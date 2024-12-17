import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
plt.figure(figsize=(9,5))
loc=np.arange(0, 26.785, 0.01)
skip=1
'''
for i in range(0,30):
    shear=np.genfromtxt('Shear_for_plotting.txt',skip_header=skip,usecols=(0),dtype=float,max_rows=2679)

    if i==29:
        plt.plot(loc,shear,color='green',label='Other load cases')
    else:
       plt.plot(loc, shear, color='green')
    skip+=2680
skip=1
'''
for i in range(0,30):
    shear=np.genfromtxt('Shear_for_plotting.txt',skip_header=skip,usecols=(0),dtype=float,max_rows=2679)
    if i==20:
        plt.plot(loc,shear,color='red',label='LC-21')
    elif i==21:
        plt.plot(loc,shear,color='blue',label='LC-22')
    '''
    elif i==25:
        plt.plot(loc,shear,color='purple',label='LC-26')
    elif i==26:
        plt.plot(loc,shear,color='black',label='LC-27')
    '''
    skip+=2680
plt.legend(loc='lower right')
plt.grid(alpha=0.9)
plt.xlim([0,27])
plt.xlabel('Spanwise position [m]',fontsize=14)
plt.ylabel('Shear [kN]',fontsize=14)
plt.xticks(np.arange(0, 28, 2))
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig("Shear_distribution.pdf", format='pdf')
plt.show()
