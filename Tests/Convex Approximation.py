# %%
import numpy as np
import matplotlib.pyplot as plt
from ipywidgets.widgets import interact
from matplotlib.axes import Axes

# Only for 2D convex set
def approximate_polygon(boundary_func, n: int):
    angles = np.linspace(0, 2*np.pi, n)
    directions = np.array([np.cos(angles), np.sin(angles)]).T
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

def sum_combination(points_1, points_2):
    result = []
    for p in points_1:
        for q in points_2:
            result.append(p + q)

    return np.array(result)


A = np.array([[1, 0], [0, 2]])
B = np.array([[2, 1], [1, 1]])
AB = A @ B
A_bound = to_boundary_func(A)
B_bound = to_boundary_func(B)
AB_bound = to_boundary_func(AB)
fig, ax = plt.subplots()
plot_convex_set(ax, A_bound)
plot_convex_set(ax, B_bound)
plot_convex_set(ax, AB_bound)
# %%

def add_boundary_func(boundary_func1, boundary_func2):
    def boundary_func(directions):
        return boundary_func1(directions) + boundary_func2(directions)

    return boundary_func

A = np.array([[1, 9], [9, 1000]])
B = np.array([[100, 0], [0, 1]])
AB = np.linalg.inv(np.linalg.inv(A) + np.linalg.inv(B))/2
A_bound = to_boundary_func(A)
B_bound = to_boundary_func(B)
AB_bound = to_boundary_func(AB)
fig, ax = plt.subplots()
plot_convex_set(ax, A_bound)
plot_convex_set(ax, B_bound)
plot_convex_set(ax, AB_bound)

sum_points = sum_combination(approximate_polygon(A_bound, 100), approximate_polygon(B_bound, 100))
ax.plot(sum_points[:,0], sum_points[0:,1], "k.", alpha=0.2)

for p in sum_points:
    if p @ AB @ p > 1:
        print("Counterexample!!!")
# %%

@interact(a=(0, 5, 0.1), b=(0, 5, 0.1), c=(0, 5, 0.1), d=(0, 5, 0.1), e=(0, 5, 0.1), f=(0, 5, 0.1))
def plot(a, b, c, d, e, f):
    A = np.array([[a, b], [b, c]])
    B = np.array([[d, e], [e, f]])
    AB = np.linalg.inv(np.linalg.inv(A) + np.linalg.inv(B))
    A_bound = to_boundary_func(A)
    B_bound = to_boundary_func(B)
    AB_bound = to_boundary_func(AB)
    fig, ax = plt.subplots()
    plot_convex_set(ax, A_bound)
    plot_convex_set(ax, B_bound)
    plot_convex_set(ax, AB_bound)
    plot_convex_set(ax, add_boundary_func(A_bound, B_bound))
# %%
# %%

points_1 = np.array([[1, 0], [0, 1], [1, 1]])
points_2 = np.array([[10, 0], [0, 10], [10, 10]])
sum_combination(points_1, points_2)
# %%
np.argmax(sum_points[0])
sum_points[:,0]
# %%
E1 = np.array([
    [1_000_000, 400_000, 400_000],
    [400_000, 160_001, 160_000],
    [400_000, 160_000, 160_001]
])

E2 = np.array([
    [6.25, 1.25, 1.25],
    [1.25, 6.5, 0.25],
    [1.25, 0.25, 6.5]
])

p1 = np.array([0.01009804, 0, -0.02774411])
p2 = np.array([-0.3647, -0.06086, -0.0652])
p12 = p1 + p2

print(p1 @ E1 @ p1)
print(p2 @ E2 @ p2)

guess = np.linalg.inv(np.linalg.inv(E1) + np.linalg.inv(E2))
print(p12 @ guess @ p12)
print(p12 @ guess @ p12)
# %%
