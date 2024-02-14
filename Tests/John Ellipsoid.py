# %%
import numpy as np
import matplotlib.pyplot as plt

# Only for 2D convex set
def approximate_polygon(boundary_func, n: int):
    angles = np.linspace(0, 2*np.pi, n)
    directions = np.array([np.cos(angles), np.sin(angles)]).T
    points = boundary_func(directions) * directions.T

    return points

def plot_convex_set(ax: plt.Axes, boundary_func):
    n = 1000
    points = approximate_polygon(boundary_func, n)
    ax.fill(points[0], points[1], alpha=0.5)

def elipsoid_to_boundary_func(elipsoid: np.ndarray):
    def boundary_func(directions):
        result = []
        for direction in directions:
            result.append(np.sqrt((direction @ direction) / (direction @ elipsoid @ direction)))
        return np.array(result)

    return boundary_func

# def polytope_to_boundary_func(polytope: np.ndarray):
#     def boundary_func(directions):
#         result = []
#         for direction in directions:
#             result.append(np.sqrt((direction @ direction) / (direction @ elipsoid @ direction)))
#         return np.array(result)

#     return boundary_func

# %%