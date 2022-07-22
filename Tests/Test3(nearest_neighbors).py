import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

import whitney as wit
import numpy as np
from Debug import *

coordinates = sample_points(30, "clusters")

root = wit.Hypercube([0, 0], 1, coordinates)
root.quadDecompose()
root.compress()
seperation_factor = 0.5
ws_pairs = root.well_separated_pairs_decomposition(seperation_factor)

test_point_index = 2
test_point = root.points[test_point_index]

nearest_neighbors = wit.all_nearest_neighbors(ws_pairs, 2)[test_point.tobytes()]
print(test_point)
print()
print(nearest_neighbors)
print()

distances = wit.metric_distance(test_point, root.points)
indices = distances.argsort()[1:3]
print(root.points[indices])
