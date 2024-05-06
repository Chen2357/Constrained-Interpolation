# %%
import numpy as np
import matplotlib.pyplot as plt
from ipywidgets.widgets import interact
from matplotlib.axes import Axes

import plotly.graph_objects as go

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

def ellipse_graph_object(boundary_func, color=None):
    n = 1000
    points = approximate_polygon(boundary_func, n)
    return go.Scatter(x=points[:,0], y=points[:,1], fill='toself', fillcolor=color)

def _inv_john_ellipsoid(points: np.ndarray, T = 10):
    m, n = points.shape
    w = np.repeat(n / m, m)

    for k in range(0, T-1):
        g = np.linalg.inv(points.T @ np.diag(w) @ points)
        w = w * np.einsum('ij,jk,ik->i', points, g, points)

    return points.T @ np.diag(w) @ points

def sum(ellipsoid1, ellipsoid2):
    return np.linalg.inv(np.linalg.inv(ellipsoid1) + np.linalg.inv(ellipsoid2))

def intersection(ellipsoids, T = 20, ax: Axes | None = None, fig: go.Figure | None = None):
    eigenvalues, eigenvectors = np.linalg.eig(ellipsoids)
    L = np.einsum("ijk,ik->ijk", eigenvectors, np.sqrt(eigenvalues))
    T_approx = np.einsum("ijk,ilk->ijl", L, L)

    if not np.allclose(T_approx, ellipsoids):
        raise ValueError("!!!")

    points = np.concatenate(np.swapaxes(L, 1, 2))
    if ax is not None:
        ax.plot(points[:,0], points[:,1], 'o')
    if fig is not None:
        fig.add_trace(go.Scatter(x=points[:,0], y=points[:,1], mode='markers', showlegend=False))

    return _inv_john_ellipsoid(points, T)

# fig, ax = plt.subplots()
fig = go.Figure()

# A = np.array([[1, 0], [0, 2]])
# B = np.array([[2, 1], [1, 1]])
A = np.array([[ 2.13, 0], [0, 100*1.473]])
B = np.array([[ 2.1282549, 1.72], [ 1.72, 100*1.350]])
C = np.array([[ 2.10, -0.657 ], [-0.657 ,  100*1.429]])

intersect = intersection([A, B, C])

A_bound = to_boundary_func(A)
B_bound = to_boundary_func(B)
C_bound = to_boundary_func(C)
intersect_bound = to_boundary_func(intersect)
# inv_bound = to_boundary_func(np.linalg.inv(intersect))
# AB_bound = to_boundary_func(AB)
# AB_inv_bound = to_boundary_func(np.linalg.inv(AB))

fig.add_trace(ellipse_graph_object(A_bound))
fig.add_trace(ellipse_graph_object(B_bound))
fig.add_trace(ellipse_graph_object(C_bound))
fig.add_trace(ellipse_graph_object(intersect_bound))
# fig.add_trace(ellipse_graph_object(inv_bound))
# plot_convex_set(ax, A_bound)
# plot_convex_set(ax, B_bound)
# plot_convex_set(ax, C_bound)
# plot_convex_set(ax, intersect_bound)
# plot_convex_set(ax, inv_bound)