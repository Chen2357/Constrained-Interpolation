# %%
import numpy as np
import matplotlib.pyplot as plt

from ipywidgets import interact

import queue
from test_module import *

# %%
def second_derivative(three_points, rad):
    rotation = np.array([
        [np.cos(rad), -np.sin(rad)],
        [np.sin(rad), np.cos(rad)]
    ])

    rotated_points = three_points @ rotation
    return 2*np.polyfit(rotated_points[:,0], rotated_points[:,1], 2)[0]

def argsort_points_along_angle(points, rad):
    rotation = np.array([
        [np.cos(rad), -np.sin(rad)],
        [np.sin(rad), np.cos(rad)]
    ])
    rotated_points = points @ rotation
    return rotated_points[:,0].argsort()


def _best(points, width, a1, a2):
    if len(points) < 3:
        if len(points) == 2:
            return (0, np.arctan2(points[1,1] - points[0,1], points[1,0] - points[0,0]))
        else:
            return (0, 0)

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

    return (minCurvature, best_rad)

def ref_axis_angle(points, width, a1, a2):
    return _best(points, width, a1, a2)[1]

def best_curavture(points, width, a1, a2):
    return _best(points, width, a1, a2)[0]

class CZDecomposition:
    def __init__(self, points: np.ndarray, values: np.ndarray, a1: float, a2: float):
        self.points = points
        self.values = values
        self.a1 = a1
        self.a2 = a2

        square = wit.Hypercube([0,0], 1, points)

        q = queue.Queue()
        q.put(square)
        while q.qsize() != 0:
            cube: wit.Hypercube = q.get()
            points = square.search_in(cube.dialated(3))
            angle = ref_axis_angle(points, cube.width, a1, a2)

            cube.ref_axis_angle = angle

            if angle is not None: continue

            cube.subdivide()
            for child in cube.children:
                q.put(child)

        self.root = square

    def worst_point(self, cube: wit.Hypercube):
        worst_point = None
        worst_curvature = -np.inf
        for point in self.root.search_in(cube.parent.dialated(3)):
            if cube.contains(point):
                continue

            curvature = best_curavture(np.concatenate([cube.points, [point]]), cube.parent.width, self.a1, self.a2)
            if curvature > worst_curvature:
                worst_point = point
                worst_curvature = curvature

        return worst_point


# %%

set1 = np.array([[x, x] for x in np.linspace(0.01, 0.49, 100)])
set2 = np.array([[x, 1.5-x] for x in np.linspace(0.51, 0.99, 100)])

points: np.ndarray = np.concatenate([set1, set2])

def f(points):
    return -np.array([(point[0]-0.5)**2 + (point[1]-0.5)**2 for point in points])
values: np.ndarray = points[:,0]

def value_of_point(point: np.ndarray):
    return values[np.where((points == point).all(axis=1))][0]

# %%
czd = CZDecomposition(points, values, 0.999, 1)

@interact(x=(0, 0.9999, 0.01), y=(0, 0.9999, 0.01))
def plot(x, y):
    fig, ax = plt.subplots(1)
    czd.root.plot(ax, edgecolor = 'b', facecolor = "None")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    ax.plot(x, y, 'rx')
    cube = czd.root.search([x, y])

    cube.dialated(1).plot(ax, edgecolor = None, facecolor = "r", alpha = 0.2)

    x = np.linspace(0, 1, 100)
    y = np.tan(cube.ref_axis_angle) * (x-0.5) + 0.5
    ax.plot(x, y, 'r')

    if len(cube.points) == 0:
        return

    worst_point = czd.worst_point(cube)

    if worst_point is not None:
        ax.plot(worst_point[0], worst_point[1], 'ks')

    rad: float = cube.ref_axis_angle
    points_in_cube = cube.points[argsort_points_along_angle(cube.points, rad)]

    rep_index = len(points_in_cube) // 2

    rotation = np.array([
        [np.cos(rad), -np.sin(rad)],
        [np.sin(rad), np.cos(rad)]
    ])

    rotated_points = points_in_cube @ rotation

    perp_diff = value_of_point(worst_point) - value_of_point(points_in_cube[rep_index])
    tangent_diff_quotient = (value_of_point(points_in_cube[rep_index+1]) - value_of_point(points_in_cube[rep_index-1])) / (rotated_points[rep_index+1][0] - rotated_points[rep_index-1][0])

    matrix = np.array([
        worst_point - points_in_cube[rep_index],
        [np.cos(rad), np.sin(rad)]
    ])

    A = np.linalg.inv(matrix) @ np.array([perp_diff, tangent_diff_quotient])

    # plot vector A at rep point
    ax.arrow(*points_in_cube[rep_index], *(A * 0.2), head_width=0.05, head_length=0.1)

    # def P(x):
    #     return value_of_point(points_in_cube[rep_index]) + (x - points_in_cube[rep_index]) @ A

    # x_min = 0
    # x_max = 1
    # x_res = 80
    # x_step = (x_max - x_min) / (x_res - 1)

    # y_min = 0
    # y_max = 1
    # y_res = 80
    # y_step = (y_max - y_min) / (y_res - 1)

    # def x(i):
    #     return x_min + i*x_step

    # def y(j):
    #     return y_min + j*y_step

    # X = x(np.arange(x_res))
    # Y = y(np.arange(y_res))
    # Z_ = np.zeros((x_res, y_res), dtype=float)
    # for i in range(x_res):
    #     for j in range(y_res):
    #         Z_[i, j] = P([x(i), y(j)])
    # Z_ = Z_.T

    # plt.pcolormesh(Y, X, Z_)
    # plt.colorbar()
    # plt.show()

    # [worst_point - Q.rep
    # ref_axis_angle]

    # P_Q(x) = value_{Q.rep} + (x - Q_rep) @ A
    # P_Q(worst_point) = value_{worst_point}
    # d(P_Q)/d(ref_axis_angle) = discrete derivative of the points
# %%
