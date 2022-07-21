import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

import whitney as wit
import whitney.Debug as wtd

import numpy as np
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

# coordinates = np.array([[0.2,0.2], [0.4,0.4], [0.9,0.9]])

root = wit.Hypercube([0, 0], 1, coordinates)
root.quadDecompose()
root.compress()
seperation_factor = 0.5
ws_pairs = wtd.disambiguate_paris(root.well_separated_pairs_decomposition(seperation_factor))

# For graphing
nonTrivialPairs: list[list[wit.Hypercube]] = []    
for pair in ws_pairs: 
    if len(pair[0].points) > 1 or len(pair[1].points) > 1:
        nonTrivialPairs.append(pair)

# Need Pseudocode
@interact(i=(0, max(len(nonTrivialPairs)-1, 0)), j=(0, max(len(ws_pairs)-1,0)))
def checkPairs(i, j):
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

    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.show()
    
