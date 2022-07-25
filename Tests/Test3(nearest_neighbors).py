import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

import whitney as wit
import numpy as np
from Debug import *
import matplotlib.pyplot as plt

coordinates = sample_points(20, "random")

root = wit.Hypercube([0, 0], 1, coordinates)
root.quad_decompose()
# root.compress()
seperation_factor = 0.5
ws_pairs = root.well_separated_pairs_decomposition(seperation_factor)

test_point_index = 0
test_point = root.points[test_point_index]

nearest_neighbors = wit.all_nearest_neighbors(ws_pairs, 3)[test_point.tobytes()]
print(test_point)
print()
print(nearest_neighbors)
print()

distances = wit.metric_distance(test_point, root.points)
indices = distances.argsort()[1:4]
print(root.points[indices])

fig, ax = plt.subplots(1)
root.plot(ax)
plt.xlim(0,1)
plt.ylim(0,1)
