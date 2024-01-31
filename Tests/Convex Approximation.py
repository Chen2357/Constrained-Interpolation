# %%
import numpy as np
import matplotlib.pyplot as plt
from sympy import plot

class ConvexSet:
    # chararteric_func: R^n -> {1, 0}
    # boundary_func: R^n -> R
    def __init__(self, chararteric_func, boundary_func):
        self.chararteric_func = chararteric_func
        self.boundary_func = boundary_func

def intersection(a: ConvexSet, b: ConvexSet):
    return ConvexSet(
        lambda x: a.chararteric_func(x) * b.chararteric_func(x),
        lambda x: np.minimum(a.boundary_func(x), b.boundary_func(x))
    )

def addition(a: ConvexSet, b: ConvexSet):
    return ConvexSet(
        lambda x: (a.chararteric_func(x) + b.chararteric_func(x)) >= np.linalg.norm(x, ord=2),
        lambda x: a.boundary_func(x) + b.boundary_func(x)
    )

# Only for 2D convex set
def approximate_polygon(convex_set: ConvexSet, n: int):
    angles = np.linspace(0, 2*np.pi, n)
    directions = np.array([np.cos(angles), np.sin(angles)]).T
    points = convex_set.boundary_func(directions) * directions.T

    return points

def plot_convex_set(ax: plt.Axes, convex_set: ConvexSet):
    n = 1000
    points = approximate_polygon(convex_set, n)
    ax.fill(points[0], points[1], alpha=0.5)

# %%

square = ConvexSet(
    lambda x: np.all(np.abs(x) <= 1, axis=1),
    lambda x: np.linalg.norm(x, ord=2, axis=1) / np.max(np.abs(x), axis=1)
)

circle = ConvexSet(
    lambda x: np.linalg.norm(x, ord=2, axis=1) <= 1.2,
    lambda x: 1.2
)

fig, ax = plt.subplots(1)
plot_convex_set(ax, square)
plot_convex_set(ax, circle)
plot_convex_set(ax, intersection(square, circle))
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)
# %%