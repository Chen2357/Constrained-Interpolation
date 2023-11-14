import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

import whitney as wit
import numpy as np
from test_module import *
import matplotlib.pyplot as plt
coordinates = sample_points(50, "random")

root = wit.Hypercube([0, 0], 1, coordinates)
root.quad_decompose()
root.compress()

test_point = np.array([0.5, 0.5])
print(root.query_nearest_point(test_point))

distances = wit.metric_distance(test_point, root.points)
index = np.argmin(distances)
print(root.points[index])
