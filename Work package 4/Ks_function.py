import matplotlib.pyplot as plt
#data points
points = [
    [1.0361445783132535, 14.858895705521473],
    [1.0451807228915664, 14.509202453987731],
    [1.0632530120481931, 14.15950920245399],
    [1.0993975903614457, 13.717791411042946],
    [1.1355421686746991, 13.276073619631903],
    [1.1897590361445785, 12.88957055214724],
    [1.2801204819277112, 12.466257668711657],
    [1.3343373493975905, 12.171779141104295],
    [1.4246987951807233, 11.822085889570552],
    [1.569277108433735, 11.435582822085891],
    [1.6867469879518073, 11.049079754601227],
    [1.8132530120481931, 10.754601226993866],
    [1.9578313253012047, 10.404907975460123],
    [2.2560240963855422, 10.036809815950921],
    [2.4457831325301207, 9.85276073619632],
    [2.7530120481927716, 9.70552147239264],
    [3.024096385542169, 9.70552147239264],
    [3.1957831325301207, 9.613496932515337],
    [3.385542168674699, 9.576687116564418],
    [3.566265060240964, 9.521472392638039],
    [3.765060240963856, 9.484662576687118],
    [3.9728915662650603, 9.484662576687118],
    [4.198795180722892, 9.484662576687118],
    [4.388554216867471, 9.503067484662576],
    [4.578313253012049, 9.503067484662576],
    [4.795180722891567, 9.521472392638039],
    [5.003012048192772, 9.503067484662576],
]

def linear_interpolation_Ks(a,b):
    x_input = a/b 
    # Check if x_input is outside the data range
    if x_input <= points[0][0]:
        return points[0][1]  # Return the first y value
    elif x_input >= points[-1][0]:
        return points[-1][1]  # Return the last y value
    
    # Find the interval [x1, x2] where x_input lies
    for i in range(len(points) - 1):
        x1, y1 = points[i]
        x2, y2 = points[i + 1]
        if x1 <= x_input <= x2:
            # Perform linear interpolation
            y_output = y1 + (y2 - y1) * (x_input - x1) / (x2 - x1)
            return y_output
# print(linear_interpolation())
# Separate the data into x and y components
x_values = [point[0] for point in points]
y_values = [point[1] for point in points]

# Plot the points and linear segments
plt.figure(figsize=(8, 6))
plt.plot(x_values, y_values, marker='o', linestyle='-', color='b', label='Linearized Data')
plt.title('Linearized Curve from Data Points')
plt.xlabel('a/b')
plt.ylabel('$k_s$')
plt.grid(True)
plt.legend()
plt.show()
