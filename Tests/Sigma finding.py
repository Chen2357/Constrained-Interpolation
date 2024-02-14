# %%
from multiprocessing import Value
import numpy as np
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
from ipywidgets.widgets import interact

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

def sum_combination(points_1, points_2):
    result = []
    for p in points_1:
        for q in points_2:
            result.append(p + q)

    return np.array(result)

def sum(ellipsoid1, ellipsoid2):
    guess = np.linalg.inv(np.linalg.inv(ellipsoid1) + np.linalg.inv(ellipsoid2))

    bound1 = to_boundary_func(ellipsoid1)
    bound2 = to_boundary_func(ellipsoid2)
    p1 = approximate_polygon_3D(bound1, 20, 10)
    p2 = approximate_polygon_3D(bound2, 20, 10)

    result = []
    for p in p1:
        for q in p2:
            result.append([p, q, p + q])

    for i in range(len(result)):
        p = result[i][2]
        x = result[i][0]
        y = result[i][1]
        if p @ guess/2 @ p > 1:
            print(ellipsoid1, ellipsoid2, p, x, y)
            raise ValueError("Counterexample!!!")

    return guess

def intersection(ellipsoid1, ellipsoid2):
    return ellipsoid1 + ellipsoid2

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
def is_sin_5pix_over_270_in_sigma(sigma, x):
    the_function = np.array([0, 5*np.pi/270, 0])
    return the_function @ pullback(sigma, forward_transformation(-x)) @ the_function < 1

C = 2

E = np.array([
    [0.2, 0.2],
    [0.4, 0.4],
    [0.6, 0.6],
    [0.8, 0.8],
])

sigma = np.array([sigma_0(x) for x in E])
sigma_history = [sigma.copy()]

def recrusion(sigma):
    new_sigma = sigma
    for i in range(len(E)):
        x = E[i]
        for j in range(len(E)):
            if i == j:
                continue
            y = E[j]
            new_sigma[i] = intersection(new_sigma[i], sum(sigma[j], scale(ball(norm(x-y), x), C)))
            if not is_sin_5pix_over_270_in_sigma(new_sigma[i], x):
                print(f"i: {i}, j: {j}, n: {len(sigma_history)}")
                print(new_sigma[i])

    return new_sigma

for i in range(7):
    sigma = recrusion(sigma)
    sigma_history.append(sigma.copy())

# %%
# fig, ax = plt.subplots()

def plot_sigma(ax, sigma, x):
    pullback_sigma = pullback(sigma, forward_transformation(-x))
    plot_convex_set(ax, to_boundary_func(pullback_sigma[1:,1:]))

# for i in range(len(E)):
#     plot_sigma(ax, sigma_history[7][i], E[i])

i = 0

# for j in range(len(sigma_history)):
#     plot_sigma(ax, sigma_history[j][i], E[i])

@interact(j=(0, len(sigma_history)-1))
def plot(j):
    fig, ax = plt.subplots()
    plot_sigma(ax, sigma_history[j][i], E[i])
    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(-0.5, 0.5)
# plot_sigma(ax, sigma_history[60][i], E[i])
# plot_sigma(ax, sigma_history[100][i], E[i])
# %%
pullback(sigma_history[7][2], forward_transformation(-E[2]))
# %%
