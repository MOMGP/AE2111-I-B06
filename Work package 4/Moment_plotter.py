import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
plt.figure(figsize=(9,5))
loc=np.arange(0, 26.785, 0.01)
skip=1
'''
for i in range(0,30):
    moment=np.genfromtxt('bending_for_plotting.txt',skip_header=skip,usecols=(0),dtype=float,max_rows=2679)
    if i==29:
        plt.plot(loc,moment,color='green',label='Other load cases')
    else:
       plt.plot(loc, moment, color='green')
    skip+=2680
'''
skip=1
for i in range(0,30):
    moment=np.genfromtxt('bending_for_plotting.txt',skip_header=skip,usecols=(0),dtype=float,max_rows=2679)
    if i==20:
        plt.plot(loc,moment,color='red',label='LC-21')
    elif i==21:
        plt.plot(loc,moment,color='blue',label='LC-22')
    '''
    elif i==25:
        plt.plot(loc,moment,color='purple',label='LC-26')
    elif i==26:
        plt.plot(loc,moment,color='black',label='LC-27')
    '''
    skip+=2680
plt.legend(loc='upper right')
plt.grid(alpha=0.9)
plt.xlim([9,27])
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel('Spanwise position [m]',fontsize=14)
plt.ylabel('Moment [kNm]',fontsize=14)
plt.xticks(np.arange(0, 28, 2))
plt.savefig("Moment_distribution.pdf", format='pdf')
plt.show()
