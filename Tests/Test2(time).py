import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

import numpy as np
import whitney as wit
import matplotlib.pyplot as plt
import timeit

# np.random.seed(223)

plt.rcParams["figure.figsize"] = [5, 5]      

r = np.arange(100,1000,100)
dt = []
for i in r: 
    print(i)
    start_time = timeit.default_timer()
    for j in range(3):
        # coordinates = np.concatenate(
        #     (
        #         np.random.rand(i,2) * 0.1 + 0.1,
        #         np.random.rand(i,2) * 0.1 + 0.9,
        #         np.random.rand(i,2) * 0.1 + np.array([0.1, 0.9])
        #     ), axis=0
        # )
        coordinates = np.random.rand(i, 2)
        
        root = wit.Hypercube([0, 0], 1, coordinates)
        root.quadDecompose()
        root.compress()
        seperation_factor = 0.5
        ws_pairs = root.well_separated_pairs_decomposition(seperation_factor)
        wit.all_nearest_neighbors(ws_pairs, 1)
        
    dt.append((timeit.default_timer()-start_time)/(1+i*np.log(i)))
    
plt.plot(r[1:], dt[1:])
plt.show()
