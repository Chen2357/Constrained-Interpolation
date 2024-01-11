# %%
import numpy as np
import matplotlib.pyplot as plt

import queue

from ipywidgets import interact
import ipywidgets as widgets

from test_module import *

# %%
def plot_quad(ax, three_points, rad, format):
    rad2 = wit.find_orientation(three_points)

    rotation = np.array([
        [np.cos(rad), -np.sin(rad)],
        [np.sin(rad), np.cos(rad)]
    ])

    rotated_points = three_points @ rotation
    p = np.poly1d(np.polyfit(rotated_points[:,0], rotated_points[:,1], 2))

    x = np.linspace(-2, 2, 100)
    y = p(x)

    interp = np.array([x, y]).T @ rotation.T

    ax.plot(interp[:,0], interp[:,1], format)

def second_derivative(three_points, rad):
    rotation = np.array([
        [np.cos(rad), -np.sin(rad)],
        [np.sin(rad), np.cos(rad)]
    ])

    rotated_points = three_points @ rotation
    return 2*np.polyfit(rotated_points[:,0], rotated_points[:,1], 2)[0]

# %%
np.random.seed(123)
points = sample_points(10, "random")

square = wit.Hypercube([0,0], 1, points)
square.quad_decompose()
# root.compress()
seperation_factor = 0.5
ws_pairs = square.well_separated_pairs_decomposition(seperation_factor)

nearest_neighbors = wit.all_nearest_neighbors(ws_pairs, 3)

@interact(x=(0, 1, 0.01), y=(0, 1, 0.01))
def plot(x, y):
    fig, ax = plt.subplots(1)

    square.plot(ax, edgecolor = 'b', facecolor = "None")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal', adjustable='box')

    ax.plot(x, y, 'kx')

    cube = square.search([x,y])
    if len(cube.points) == 0:
        return

    point = cube.points[0]
    # square.query_nearest_point([x,y])
    three_points = np.concatenate([nearest_neighbors[point.tobytes()][0:2], [point]])

    rad = wit.find_orientation(three_points)

    ax.plot(point[0], point[1], 'rx')
    ax.annotate("B", (point[0], point[1]))
    ax.annotate("Q", (three_points[1,0], three_points[1,1]))

    x = np.linspace(0, 1, 100)
    y = np.tan(rad) * (x-0.5) + 0.5
    ax.plot(x, y, 'r--')

    plot_quad(ax, three_points, rad, 'r')

    next_point = three_points[0]
    next_next_point: npt.NDArray

    ax.plot(next_point[0], next_point[1], 'gx')
    ax.annotate("P", (next_point[0], next_point[1]))

    for p in nearest_neighbors[next_point.tobytes()]:
        if not np.all(p == three_points[1]) and not np.all(p == point):
            next_next_point = p
            break

    ax.annotate("R", (next_next_point[0], next_next_point[1]))

    next_three_point = np.array([point, next_point, next_next_point])
    next_rad = wit.find_orientation(next_three_point)

    x = np.linspace(0, 1, 100)
    y = np.tan(next_rad) * (x-0.5) + 0.5
    ax.plot(x, y, 'g--')

    plot_quad(ax, next_three_point, next_rad, 'g')

    print("DIfference in orientation (rad):", np.min([(rad - next_rad) % (np.pi), (next_rad - rad) % (np.pi)]))
    print("Second derivative change of BPQ from red to green orientation:", np.abs(second_derivative(three_points, rad)) - np.abs(second_derivative(three_points, next_rad)))
    print("Second derivative change of BPR from green to red orientation:", np.abs(second_derivative(next_three_point, next_rad)) - np.abs(second_derivative(next_three_point, rad)))

# %%
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

# %%
np.random.seed(123)
# points = sample_points(100, "random")

set1 = np.array([[x,x] for x in np.linspace(0.01, 0.49, 100)])
set2 = np.array([[x, 1.5-x] for x in np.linspace(0.51, 0.99, 100)])

points = np.concatenate([set1, set2])

square = wit.Hypercube([0,0], 1, points)

# biggest = wit.Hypercube([0,0], 0)

a1 = 1
a2 = 1

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

@interact(x=(0, 1, 0.01), y=(0, 1, 0.01))
def plot(x, y):
    fig, ax = plt.subplots(1)
    square.plot(ax, edgecolor = 'b', facecolor = "None")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    ax.plot(x, y, 'rx')
    cube = square.search([x,y])

    cube.dialated(1).plot(ax, edgecolor = None, facecolor = "r", alpha = 0.2)

    x = np.linspace(0, 1, 100)
    y = np.tan(cube.ref_axis_angle) * (x-0.5) + 0.5
    ax.plot(x, y, 'r')

    for point in square.search_in(cube.parent.dialated(3)):
        if point in cube.points:
            continue

        angle = ref_axis_angle(np.concatenate([cube.points, [point]]), cube.parent.width, a1, a2)
        if angle is None:
            ax.plot(point[0], point[1], 'gx')