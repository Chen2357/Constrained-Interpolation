import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

import whitney as wit
from test_module import *

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from ipywidgets import interact

# Coords
def f(x):
   return (x-0.5)**2+0.25

x_vals1 = np.linspace(0.1, 0.4, 10)
x_vals2 = np.linspace(0.6, 0.9, 10)
x_vals = np.hstack([x_vals1, x_vals2, 0.5])

coordinates = None
for i in range(len(x_vals)):
    if i == 0:
        coordinates = [x_vals[i], f(x_vals[i])]
    else:
        coordinates = np.vstack([coordinates, [x_vals[i], f(x_vals[i])]])
coordinates = np.vstack([coordinates, [0.5, 0.75]])

# Creating WSPD
test_point = coordinates[0]
root = wit.Hypercube([0, 0], 1, coordinates)
root.quad_decompose()
root.compress()
seperation_factor = 0.5
ws_pairs = disambiguate_paris(root.well_separated_pairs_decomposition(seperation_factor))
nearest_neighbors = wit.all_nearest_neighbors(ws_pairs, 3)[test_point.tobytes()]


# For graphing
nonTrivialPairs: list[list[wit.Hypercube]] = []
for pair in ws_pairs:
    if len(pair[0].points) > 1 or len(pair[1].points) > 1:
        nonTrivialPairs.append(pair)

# Test point
pairs_including_test_point = []
for pair in nonTrivialPairs:
    if pair[0].contains(test_point) or pair[1].contains(test_point):
        pairs_including_test_point.append(pair)

# Need Pseudocode
@interact(i=(0, max(len(nonTrivialPairs)-1, 0)), j=(0, max(len(ws_pairs)-1,0)), k=(0, max(len(pairs_including_test_point)-1,0)))
def checkPairs(i, j, k):
    fig, ax = plt.subplots(1)
    root.plot(ax)

    if len(nonTrivialPairs) != 0:
        pair = nonTrivialPairs[i]
        rectangles = [
            Rectangle(pair[0].pos, pair[0].width, pair[0].width),
            Rectangle(pair[1].pos, pair[1].width, pair[1].width)
            ]
        collection = PatchCollection(rectangles, alpha=0.5, edgecolor='k',facecolor='r')
        ax.add_collection(collection)

    if len(ws_pairs) != 0:
        pair = ws_pairs[j]
        rectangles = [
            Rectangle(pair[0].pos, pair[0].width, pair[0].width),
            Rectangle(pair[1].pos, pair[1].width, pair[1].width)
            ]
        collection = PatchCollection(rectangles, alpha=0.5, edgecolor='k',facecolor='m')
        ax.add_collection(collection)

    if len(pairs_including_test_point) != 0:
        pair = pairs_including_test_point[k]
        rectangles = [
            Rectangle(pair[0].pos, pair[0].width, pair[0].width),
            Rectangle(pair[1].pos, pair[1].width, pair[1].width)
            ]
        collection = PatchCollection(rectangles, alpha=0.5, edgecolor='k',facecolor='g')
        ax.add_collection(collection)

    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.show()

