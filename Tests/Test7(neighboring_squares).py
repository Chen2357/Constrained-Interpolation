import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

import whitney as wit
import numpy as np
from Debug import *
import matplotlib.pyplot as plt
coordinates = sample_points(25, "random")

root = wit.Hypercube([0, 0], 1, coordinates)
root.quad_decompose()
# root.compress()
seperation_factor = 0.5
ws_pairs = root.well_separated_pairs_decomposition(seperation_factor)


root2 = wit.Hypercube([0, 0], 1, coordinates)
root2.whitney_decompose()


fig, ax = plt.subplots(1)
# plt.figure(figsize=(5, 5))
root.plot(ax, edgecolor = 'b')


all_nearest_neighbors = wit.all_nearest_neighbors(ws_pairs, 2)
for test_point_index in range(1):
    test_point = root.points[test_point_index]
    print(test_point)
    nearest_neighbors = all_nearest_neighbors[test_point.tobytes()]
    
    distance = wit.metric_distance(test_point, nearest_neighbors[0])
    wit_squares = root.whitney_square2(test_point, distance)
    unique_wit_squares = []
    for square in wit_squares:
        if square in unique_wit_squares:
            continue
        else:
            unique_wit_squares.append(square)

    for square in unique_wit_squares:
        # print(square.pos, square.width)
        square.plot(ax, facecolor = 'r')
        
