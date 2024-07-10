# %% Imports
import sys, os
sys.path.insert(1, os.path.join(sys.path[0], '..'))

import whitney as wit
import numpy as np
import queue

# %% Numerical CZ Decomposition helper functions
def ref_axis_angle(points, width, a1, a2):
    if len(points) < 3:
        if len(points) == 2:
            return np.arctan2(points[1,1] - points[0,1], points[1,0] - points[0,0])
        else:
            return 0

    best_rad = None
    minCurvature = a2 / width

    for rad in np.arange(0, np.pi, np.pi/20):
        rotation = np.array([
            [np.cos(rad), -np.sin(rad)],
            [np.sin(rad), np.cos(rad)]
        ])

        rotated_points = points @ rotation
        rotated_points = rotated_points[rotated_points[:,0].argsort()]

        if np.any(np.diff(rotated_points[:,0]) == 0):
            continue
        maxSlope = np.max(np.abs(np.diff(rotated_points[:,1])/np.diff(rotated_points[:,0])))
        maxCurvature = np.max(np.abs(np.diff(np.diff(rotated_points[:,1])/np.diff(rotated_points[:,0]))))

        if maxSlope <= a1 and maxCurvature <= minCurvature:
            minCurvature = maxCurvature
            best_rad = rad

    return best_rad

def get_numerical_CZ(E, a1, a2):
    numerical_CZ = wit.Hypercube(np.array([0, 0]), 1, E)

    def is_ok(cube: wit.Hypercube):
        points = numerical_CZ.search_in(cube.dialated(3))
        angle = ref_axis_angle(points, cube.width, a1, a2)

        return angle is not None

    numerical_CZ._decompose(is_ok)

    return numerical_CZ

def inject_comparison_data(root: wit.Hypercube, a1, a2):
    q = queue.Queue()
    q.put(root)
    while q.qsize() != 0:
        square = q.get()
        square: wit.Hypercube

        points = root.search_in(square.dialated(3))
        angle = ref_axis_angle(points, square.width, a1, a2)

        square.is_CZ = angle is not None # type: ignore
        square.CZ_generation = 0 if angle is None else (1 if square.parent is None else square.parent.CZ_generation + 1) # type: ignore

        for child in square.children:
            if child is not None:
                q.put(child)

def CZ_comparison_fillcolor_map(cube: wit.Hypercube):
    return f"rgba(100, 100, 100, {cube.CZ_generation / 5})" if cube.is_CZ else "rgba(255, 0, 0, 0.75)" # type: ignore
# %% Compute numerical CZ decomposition
E = np.random.rand(20, 2)

a1 = 1
a2 = 0.01

numerical_CZ = get_numerical_CZ(E, a1, a2)
wit.Plotting.plot_hypercube(numerical_CZ).show()

# %% CZ decomposition comparison
root = wit.Hypercube((0, 0), 1, E)
wit.CZ_decompose(root, post_shrinking=0.25)

inject_comparison_data(root, a1, a2)

all_CZ = np.all([leaf.is_CZ for leaf in root.leaves]) # type: ignore
if all_CZ:
    print("All leaves are CZ")
else:
    print("Not all leaves are CZ")

wit.Plotting.plot_hypercube(root, fillcolor_map=CZ_comparison_fillcolor_map, opacity=1).show()