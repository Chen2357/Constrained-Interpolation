# %%
from multiprocessing import Value
import numpy as np
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
from ipywidgets.widgets import interact
from soupsieve import closest

def forward_transformation(x):
    return np.array([
        [1, x[0], x[1]],
        [0, 1, 0],
        [0, 0, 1]
    ])

def pullback(ellipsoid, forward):
    return forward.T @ ellipsoid @ forward

def ball(delta, x):
    return pullback(np.array([
        [1/delta**2, 0, 0],
        [0, 1/delta**2, 0],
        [0, 0, 1/delta**2]
    ]), forward_transformation(x))

thickness = 0.001

def sigma_0(x):
    return pullback(np.array([
        [1/thickness**2, 0, 0],
        [0, 1, 0],
        [0, 0, 1]
    ]), forward_transformation(x))

def scale(ellipsoid, c):
    return ellipsoid / c**2

def sum(ellipsoid1, ellipsoid2):
    return np.linalg.inv(np.linalg.inv(ellipsoid1) + np.linalg.inv(ellipsoid2))

def intersection(ellipsoid1, ellipsoid2):
    return (ellipsoid1 + ellipsoid2)

def norm(x):
    return np.linalg.norm(x, ord=np.inf)

def approximate_polygon(boundary_func, n: int):
    angles = np.linspace(0, 2*np.pi, n)
    directions = np.array([np.cos(angles), np.sin(angles)]).T
    points = boundary_func(directions) * directions.T

    return points.T

def approximate_polygon_3D(boundary_func, n: int, m: int):
    phi = np.linspace(0, 2*np.pi, n)
    theta = np.linspace(0, np.pi, m)
    directions = []
    for p in phi:
        for t in theta:
            directions.append([np.sin(t)*np.cos(p), np.sin(t)*np.sin(p), np.cos(t)])
    directions = np.array(directions)
    points = boundary_func(directions) * directions.T

    return points.T

def plot_convex_set(ax: Axes, boundary_func, color=None):
    n = 1000
    points = approximate_polygon(boundary_func, n)
    if color is not None:
        ax.fill(points[:,0], points[:,1], alpha=0.5, color=color)
    else:
        ax.fill(points[:,0], points[:,1], alpha=0.5)

def to_boundary_func(elipse: np.ndarray):
    def boundary_func(directions):
        result = []
        for direction in directions:
            result.append(np.sqrt((direction @ direction) / (direction @ elipse @ direction)))
        return np.array(result)

    return boundary_func
# %%
def nearest_k_points(points, k):
    """
    Returns an n by k array of the k nearest points to each point in the n by d array points
    """
    distances = np.sum(points**2, axis=1).reshape(-1, 1) - 2 * points @ points.T + np.sum(points**2, axis=1).reshape(1, -1)
    return np.argsort(distances, axis=1)[:, 1:k+1]

def is_pos_def(x):
    return np.all(np.linalg.eigvals(x) > 0)
# %%
# def is_sin_5pix_over_270_in_sigma(sigma, x):
#     the_function = np.array([0, 5*np.pi/270, 0])
#     return the_function @ pullback(sigma, forward_transformation(-x)) @ the_function < 1

C = 2

dim = 2
E = np.array([
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
# E = np.random.rand(10, 2)

inspection_index = 1

closest_points_info = nearest_k_points(E, 2)

sigma = np.array([sigma_0(x) for x in E])
sigma_history = [sigma.copy()]

artifical_expand = 1.5
artifical_shrink = 0.5

def recrusion(sigma):
    new_sigma = sigma
    for i in range(len(E)):
        x = E[i]

        for j in range(len(E)):
            if i == j:
                continue
            y = E[j]
            distance = norm(x-y)

            new_sigma[i] = intersection(new_sigma[i], sum(sigma[j], scale(ball(distance, x), C)))

        for j in closest_points_info[i]:
            closest_point = E[j]

            closest_point_unit_direction = (closest_point - x) / np.linalg.norm(closest_point - x, ord=2)

            projection = np.outer(closest_point_unit_direction, closest_point_unit_direction)

            shrink = 1/np.sqrt(artifical_expand) * (np.eye(dim) - projection) + 1/np.sqrt(artifical_shrink) * projection

            shrink = np.block([
                [1, np.zeros((1, dim))],
                [np.zeros((dim, 1)), shrink]
            ])

            new_sigma[i] = pullback(shrink @ pullback(new_sigma[i], forward_transformation(-x)) @ shrink, forward_transformation(x))

    return new_sigma

for i in range(10):
    sigma = recrusion(sigma)
    sigma_history.append(sigma.copy())
# %%
# sigma_history = np.array(sigma_history)
# scale = 1/np.linalg.det(sigma_history[:,0])

# j = 0
# size = np.zeros(len(sigma_history))

# for i in range(len(sigma_history)):
#     size[i] = 1/np.linalg.det(pullback(sigma_history[i][j], forward_transformation(-E[j]))[1:,1:])

# fig, ax = plt.subplots()
# ax.plot(np.arange(1, 100), size[:99])
# ax.set_yscale('log')
# ax.set_xscale('log')
# %%
def plot_sigma_at(ax, sigma, x):
    pullback_sigma = pullback(sigma, forward_transformation(-x))
    n = 1000
    points = approximate_polygon(to_boundary_func(pullback_sigma[1:,1:]), n)
    points = points + x
    ax.fill(points[:,0], points[:,1], alpha=0.5)

# for i in range(len(E)):
#     plot_sigma(ax, sigma_history[7][i], E[i])

# i = 0

# for j in range(len(sigma_history)):
#     plot_sigma(ax, sigma_history[j][i], E[i])

@interact(j=(0, len(sigma_history)-1))
def plot_at(j):
    fig, ax = plt.subplots()
    print(pullback(sigma_history[j][0], forward_transformation(-E[0])))
    for i in range(len(E)):
        plot_sigma_at(ax, sigma_history[j][i], E[i])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
# %%
