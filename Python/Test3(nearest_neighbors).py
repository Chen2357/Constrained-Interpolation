import numpy as np
import whitney as wit

import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from ipywidgets import interact

coordinates = np.concatenate(
    (
        np.random.rand(10,2) * 0.1 + 0.1,
        np.random.rand(10,2) * 0.1 + 0.9,
        np.random.rand(10,2) * 0.1 + np.array([0.1, 0.9])
    ), axis=0
)

trivial_coordinates = np.array([[0.2,0.2], [0.4,0.4], [0.9,0.9]])


root = wit.Hypercube([0, 0], 1, coordinates)
root.quadDecompose()
root.compress()
seperation_factor = 0.5
ws_pairs = root.well_separated_pairs_decomposition(seperation_factor)


filtered_pairs = wit.filter_pairs(ws_pairs, 2)
nearest_neighbors, distances = wit.find_nearest_neighbor(filtered_pairs, root.points[0], 2)
print(root.points[0])
print(nearest_neighbors)
print(distances)
print()

nearest_index = np.argmin(np.max(np.abs(root.points[1:] - root.points[0]),1))
print(root.points[nearest_index + 1])

