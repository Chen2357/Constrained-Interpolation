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

plt.rcParams["figure.figsize"] = [7, 5]      

trials = 3
seperation_factor = 0.5
k = 2

r = np.arange(3, 10000, 250)
time = []
num_pairs = []


# Tests time of nested funcs
for i in r: 
    print(i)
    start_time = timeit.default_timer()
    total_pairs = 0
    for j in range(trials):
        
        # coordinates = np.concatenate(
        #     (
        #         np.random.rand(i//3,2) * 0.1 + 0.1,
        #         np.random.rand(i//3,2) * 0.1 + 0.9,
        #         np.random.rand(i//3,2) * 0.1 + np.array([0.1, 0.9])
        #     ), axis=0
        # )
        
        coordinates = np.random.rand(i,2)
        
        ## Functions which we are testing
        root = wit.Hypercube([0, 0], 1, coordinates)
        root.quadDecompose()
        root.compress()
        ws_pairs = root.well_separated_pairs_decomposition(seperation_factor)
        nearest_neighbors = wit.all_nearest_neighbors(ws_pairs, k)
        
        total_pairs += len(ws_pairs)
        
    time.append((timeit.default_timer() - start_time) / trials / (1 + i * np.log(i)))
    num_pairs.append((total_pairs / trials)/(1 + i))


title = "Computing {} nearest neighbors for each point \n Trials per n: {}, Seperation factor: {}".format(k, trials, seperation_factor)
fig, axs = plt.subplots(2, 1)

axs[0].set_title(title)
axs[0].plot(r[1:], time[1:]) # log plot
axs[0].set_ylabel('Time / nlog(n)')
axs[0].grid(True)

axs[1].plot(r[1:], num_pairs[1:])
axs[1].set_ylabel('#total pairs in WSPD / n')
axs[1].set_xlabel('n')
axs[1].grid(True)

fig.tight_layout()
plt.show()
