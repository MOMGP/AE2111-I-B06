import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
plt.figure(figsize=(9,4.5))
loc=np.arange(0, 26.785, 0.01)
skip=1
'''
for i in range(0,30):
    torsion=np.genfromtxt('torsion_for_plotting.txt',skip_header=skip,usecols=(0),dtype=float,max_rows=2679)
    if i==29:
        plt.plot(loc,torsion,color='green',label='Other load cases')
    else:
       plt.plot(loc, torsion, color='green')
    skip+=2680
skip=1
'''
for i in range(0,30):
    torsion=np.genfromtxt('torsion_for_plotting.txt',skip_header=skip,usecols=(0),dtype=float,max_rows=2679)
    if i==20:
        print(1)
        #plt.plot(loc,torsion,color='red',label='LC-21')
    elif i==21:
        plt.plot(loc,torsion,color='blue',label='LC-22')
    elif i==25:
        plt.plot(loc,torsion,color='purple',label='LC-26')
    elif i==26:
        print(1)
        #plt.plot(loc,torsion,color='black',label='LC-27')
    skip+=2680

plt.legend(loc='upper right')
plt.grid(alpha=0.9)
plt.xlim([0,27])
plt.xlabel('Spanwise position [m]',fontsize=13)
plt.ylabel('Torque [kNm]',fontsize=13)
plt.xticks(np.arange(0, 28, 2))
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig("Torque_distribution.pdf", format='pdf')
plt.show()


