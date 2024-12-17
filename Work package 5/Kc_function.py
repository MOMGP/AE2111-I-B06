#Kc function 
from matplotlib import pyplot as plt
import numpy as np
e = 2.71828182846
def Kc(a,b):
    x = a/b
    if x <= 1.0: 
        value = e**(-14.2*x+8) + 4
    elif  1.0 < x and x <= 1.4: 
        value = (5/4)*x+2.75
    elif 1.4 < x and x <= 2.0:
        value = -(5/6)*x + 17/3
    elif 2.0 < x and x < 2.4:
        value = 0.5*x + 3
    elif 2.4 < x and x <= 3.0: 
        value = -(1/3)*x + 5
    elif x > 3.0: 
        value = 4
    return value

# print(Kc(3.4,5.6))


    
