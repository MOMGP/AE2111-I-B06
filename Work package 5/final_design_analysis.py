import numpy as np
designs = np.load("Work package 5\\Results\\Design_1.npy")
index = np.argmin(designs, axis=0)
print(index)
print(designs[12])