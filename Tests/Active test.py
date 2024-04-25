import numpy as np
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.patches import Rectangle

from test_module import *

points = np.array([
    [0.01, 0.9],
    [0.01, 0.8],
    [0.01, 0.7],
    [0.01, 0.6],
    [0.01, 0.5],
    [0.01, 0.4],
    [0.01, 0.01],
    [0.4, 0.01],
    [0.5, 0.01],
    [0.6, 0.01],
    [0.7, 0.01],
    [0.8, 0.01],
    [0.9, 0.01],
])

points = sample_points(25, 'clusters')
root = wit.Hypercube([0, 0], 1, points)
cz = wit.CZ_Decomposition(root)

fig, ax = plt.subplots()
ax.set_aspect('equal')
# ax.set_xlim(0, 1)
# ax.set_ylim(0, 1)
ax.set_xlim(0, 0.3)
ax.set_ylim(0, 0.3)
root.plot(ax)