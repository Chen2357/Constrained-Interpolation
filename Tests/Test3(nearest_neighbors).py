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

test_point_index = 0
test_point = root.points[test_point_index]

nearest_neighbors = wit.all_nearest_neighbors(ws_pairs, 2)[test_point.tobytes()]
print(root.points[0])
print(nearest_neighbors)

nearest_index = np.argmin(wit.metric_distance(np.delete(root.points, test_point_index, 0), test_point))
print(root.points[nearest_index + (1 if nearest_index >= test_point_index else 0)])
