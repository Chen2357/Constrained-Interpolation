import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

import numpy as np
import whitney as wit
import matplotlib.pyplot as plt
import timeit
from test_module import *

# np.random.seed(223)

trials = 3
seperation_factor = 0.5
k = 2

r = np.arange(3, 10000, 250)
time = []
num_pairs = []

for i in r:
    print(i)
    start_time = timeit.default_timer()
    total_pairs = 0
    for j in range(trials):
        coordinates = sample_points(i, "random")

        root = wit.Hypercube([0, 0], 1, coordinates)
        root.quad_decompose()
        root.compress()
        ws_pairs = root.well_separated_pairs_decomposition(seperation_factor)

        # total_pairs += len(ws_pairs)

        all_nearest_neighbors = wit.all_nearest_neighbors(ws_pairs, k)

        for test_point_index in range(i):
            test_point = root.points[test_point_index]
            nearest_neighbors = all_nearest_neighbors[test_point.tobytes()]

            distance = wit.metric_distance(test_point, nearest_neighbors[0])
            # wit_square = root.whitney_square(test_point, distance, nearest_neighbors)
            wit_square = root.whitney_square(test_point, distance)

    time.append((timeit.default_timer() - start_time) / trials / (1 + i * np.log(i)))
    num_pairs.append((total_pairs / trials)/(1 + i))


title = "Computing {} nearest neighbors for each point \n Trials per n: {}, Seperation factor: {}".format(k, trials, seperation_factor)
plt.figure(figsize=(7, 5))
fig, axs = plt.subplots(2, 1)

axs[0].set_title(title)
axs[0].plot(r[1:], time[1:]) # log plot
axs[0].set_ylabel('Time / nlog(n)')
axs[0].grid(True)

# axs[1].plot(r[1:], num_pairs[1:])
# axs[1].set_ylabel('#total pairs in WSPD / n')
# axs[1].set_xlabel('n')
# axs[1].grid(True)

fig.tight_layout()
plt.show()
