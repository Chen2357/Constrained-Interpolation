# %%
import numpy as np
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
from ipywidgets.widgets import interact

from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.patches import Rectangle

import queue
from test_module import *

# %% Ellipsoid operations

def forward_transformation(x):
    return np.array([
        [1, x[0], x[1]],
        [0, 1, 0],
        [0, 0, 1]
    ])

def pullback(ellipsoid, forward):
    return forward.T @ ellipsoid @ forward

def scale(ellipsoid, c):
    return ellipsoid / c**2

def sum(ellipsoid1, ellipsoid2):
    return np.linalg.inv(np.linalg.inv(ellipsoid1) + np.linalg.inv(ellipsoid2))

def _inv_john_ellipsoid(points: np.ndarray, T = 10):
    m, n = points.shape
    w = np.repeat(n / m, m)

    for k in range(0, T-1):
        g = np.linalg.inv(points.T @ np.diag(w) @ points)
        w = w * np.einsum('ij,jk,ik->i', points, g, points)

    return points.T @ np.diag(w) @ points

def intersection(ellipsoids, T = 10):
    L = np.linalg.cholesky(ellipsoids)
    points = np.concatenate(np.swapaxes(L, 1, 2))
    return _inv_john_ellipsoid(points, T)

def norm(x):
    return np.linalg.norm(x, ord=np.inf)

def diameter(ellipsoids: np.ndarray):
    return 2 / np.sqrt(np.linalg.eig(ellipsoids)[0].min(axis=-1))

# %% Plotting
def _approximate_polygon(boundary_func, n: int):
    angles = np.linspace(0, 2*np.pi, n)
    directions = np.array([np.cos(angles), np.sin(angles)]).T
    points = boundary_func(directions) * directions.T

    return points.T

def _to_boundary_func(elipse: np.ndarray):
    def boundary_func(directions):
        result = []
        for direction in directions:
            result.append(np.sqrt((direction @ direction) / (direction @ elipse @ direction)))
        return np.array(result)

    return boundary_func

def _plot_convex_set(ax: Axes, boundary_func, color=None):
    n = 1000
    points = _approximate_polygon(boundary_func, n)
    if color is not None:
        ax.fill(points[:,0], points[:,1], alpha=0.5, color=color)
    else:
        ax.fill(points[:,0], points[:,1], alpha=0.5)

def plot_ellipse(ax: Axes, ellipse: np.ndarray, color=None):
    _plot_convex_set(ax, _to_boundary_func(ellipse), color)

# %% Sigma finding

thickness = 0.001

def sigma_0(x):
    return pullback(np.array([
        [1/thickness**2, 0, 0],
        [0, 1, 0],
        [0, 0, 1]
    ]), forward_transformation(x))

def ball(delta, x):
    return pullback(np.array([
        [1/delta**4, 0, 0],
        [0, 1/delta**2, 0],
        [0, 0, 1/delta**2]
    ]), forward_transformation(x))

C = 1
dim = 2

def find_sigma(points: np.ndarray, N = 6):
    sigma = np.array([sigma_0(x) for x in points])

    def recrusion(sigma):
        new_sigma = sigma
        for i in range(len(points)):
            x = points[i]

            intersectands = [new_sigma[i]]

            for j in range(len(points)):
                if i == j:
                    continue
                y = points[j]
                distance = norm(x-y)

                intersectands.append(sum(new_sigma[j], scale(ball(distance, x), C)))

            new_sigma[i] = intersection(intersectands)

        return new_sigma

    for _ in range(6):
        sigma = recrusion(sigma)

    return np.array([pullback(sigma[i], forward_transformation(-points[i])) for i in range(len(points))])
# %% Approximate CZ decomposition

def asvoid(arr):
    """View the array as dtype np.void (bytes)
    This collapses ND-arrays to 1D-arrays, so you can perform 1D operations on them.
    https://stackoverflow.com/a/16216866/190597 (Jaime)"""
    arr = np.ascontiguousarray(arr)
    return arr.view(np.dtype((np.void, arr.dtype.itemsize * arr.shape[-1])))

def find_index(arr, x):
    arr_as1d = asvoid(arr)
    x = asvoid(x)
    return np.nonzero(arr_as1d == x)[0]

def approximate_CZ(points: np.ndarray, a: float = 10) -> wit.Hypercube:
    root = wit.Hypercube(np.array([0, 0]), 1, points)
    sigma = find_sigma(points)

    def indices_of_points(p):
        return [find_index(points, q)[0] for q in p]

    def is_good(square: wit.Hypercube):
        p = root.search_in(square.dialated(3))
        i = indices_of_points(p)
        d = diameter(sigma[i])
        return np.all(d >= a * square.width)

    q = queue.Queue()
    q.put(root)

    while q.qsize() != 0:
        current = q.get()
        current: wit.Hypercube
        if is_good(current):
            continue
        else:
            current.subdivide()
            for child in current.children:
                q.put(child)

    return root
# %% Verify CZ decomposition

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

# %% Testing

# E = np.array([
#     [0.01, 0.9],
#     [0.01, 0.8],
#     [0.01, 0.7],
#     [0.01, 0.6],
#     [0.01, 0.5],
#     [0.01, 0.4],
#     [0.01, 0.01],
#     [0.4, 0.01],
#     [0.5, 0.01],
#     [0.6, 0.01],
#     [0.7, 0.01],
#     [0.8, 0.01],
#     [0.9, 0.01],
# ])
np.random.seed(0)
E = sample_points(10, 'clusters')

a1 = 1
a2 = 1
root = approximate_CZ(E, 10)
assign_is_CZ(root, a1, a2)

# %% Plotting
def fillcolor_map(cube: wit.Hypercube):
    return f"rgba(0, 206, 250, {cube.CZ_generation / 5})" if cube.is_CZ else "rgba(255, 0, 0, 0.75)"

wit.Plotting.plot_hypercube(root, fillcolor_map=fillcolor_map, opacity=1)
# %%

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

wit.Plotting.plot_hypercube(actual_CZ)