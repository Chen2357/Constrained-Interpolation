import numpy as np

import queue
from test_module import *

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

def assign_is_CZ(root: wit.Hypercube, a1, a2):
    q = queue.Queue()
    q.put(root)
    while q.qsize() != 0:
        square = q.get()
        square: wit.Hypercube

        points = root.search_in(square.dialated(3))
        angle = ref_axis_angle(points, square.width, a1, a2)

        square.ref_axis_angle = angle
        square.is_CZ = angle is not None
        square.CZ_generation = 0 if angle is None else (1 if square.parent is None else square.parent.CZ_generation + 1)

        for child in square.children:
            if child is not None:
                q.put(child)

def CZ_comparison_fillcolor_map(cube: wit.Hypercube):
    return f"rgba(100, 100, 100, {cube.CZ_generation / 5})" if cube.is_CZ else "rgba(255, 0, 0, 0.75)"

# wit.Plotting.plot_hypercube(root, fillcolor_map=fillcolor_map, opacity=1)
# %%

def get_actual_CZ(E, a1, a2):
    actual_CZ = wit.Hypercube(np.array([0, 0]), 1, E)
    q = queue.Queue()

    q.put(actual_CZ)
    while q.qsize() != 0:
        cube: wit.Hypercube = q.get()
        points = actual_CZ.search_in(cube.dialated(3))
        angle = ref_axis_angle(points, cube.width, a1, a2)

        cube.ref_axis_angle = angle

        if angle is not None: continue

        cube.subdivide()
        for child in cube.children:
            q.put(child)

    return actual_CZ