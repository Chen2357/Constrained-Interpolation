import numpy as np
import whitney as wit
import matplotlib.pyplot as plt
import timeit

np.random.seed(223)

start_time = timeit.default_timer()
plt.rcParams["figure.figsize"] = [5, 5]      

# Tests time of nested func
r = np.arange(3,10,1)
dt = []
for i in r: 
    print(i)
    for j in range(2):
        coordinates = np.concatenate(
            (
                np.random.rand(i,2) * 0.1 + 0.1,
                np.random.rand(i,2) * 0.1 + 0.9,
                np.random.rand(i,2) * 0.1 + np.array([0.1, 0.9])
            ), axis=0
        )
        ## Functions which we are testing
        start_time  = timeit.default_timer()
        root = wit.Hypercube([0, 0], 1, coordinates)
        root.quadDecompose()
        root.compress()
        pointy = root.find_nearest_neighbor(coordinates[0])
    dt.append((timeit.default_timer()-start_time)/(1+i*np.log(i))) # log plot
    #dt.append((timeit.default_timer()-start_time)) # linear plot
    
plt.plot(r[1:], dt[1:]) # log plot
# plt.plot(r*np.log(r), dt) # linear plot
plt.show()