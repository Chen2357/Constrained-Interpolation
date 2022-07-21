import numpy as np
import whitney as wit
import matplotlib.pyplot as plt
import timeit

#   np.random.seed(223)

plt.rcParams["figure.figsize"] = [5, 5]      

# Tests time of nested func
r = np.arange(3,500,50)
dt = []
for i in r: 
    print(i)
    start_time = timeit.default_timer()
    for j in range(5):
        coordinates = np.concatenate(
            (
                np.random.rand(i,2) * 0.1 + 0.1,
                np.random.rand(i,2) * 0.1 + 0.9,
                np.random.rand(i,2) * 0.1 + np.array([0.1, 0.9])
            ), axis=0
        )
        ## Functions which we are testing
        root = wit.Hypercube([0, 0], 1, coordinates)
        root.quadDecompose()
        root.compress()
        seperation_factor = 0.5
        wspairs = root.well_separated_pairs_decomposition(seperation_factor)
        
        
        filtered_pairs = wit.filter_pairs(wspairs, 2)
        nearest_neighbors, distances = wit.find_nearest_neighbor(filtered_pairs, root.points[0], 2)
        
    dt.append((timeit.default_timer()-start_time)/(1+i*np.log(i))) # log plot
    #dt.append((timeit.default_timer()-start_time)) # linear plot
    
plt.plot(r[1:], dt[1:]) # log plot
# plt.plot(r*np.log(r), dt) # linear plot
plt.show()
