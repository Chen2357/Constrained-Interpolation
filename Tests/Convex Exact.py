# %%
import numpy as np
import matplotlib.pyplot as plt
from ipywidgets.widgets import interact
from matplotlib.axes import Axes

# %%
def to_covector_to_vector_func(ellipse: np.ndarray):
    return lambda x: (x @ ellipse)

def sum(covector_to_vector_func1, covector_to_vector_func2):
    def result_func(x):
        vector1 = covector_to_vector_func1(x)
        vector2 = covector_to_vector_func2(x)
        length1 = np.sqrt(np.sum(x * vector1, axis=1))[:,None]
        length2 = np.sqrt(np.sum(x * vector2, axis=1))[:,None]

        vectors = vector1 / length1 + vector2 / length2

        return vectors * np.sum(x * vectors, axis=1)[:,None]

    return result_func

def approximate_polygon(covector_to_vector_func, n: int):
    angles = np.linspace(0, 2*np.pi, n)
    covectors = np.array([np.cos(angles), np.sin(angles)]).T
    vectors = covector_to_vector_func(covectors)
    return vectors / np.sqrt(np.sum(covectors * vectors, axis=1))[:,None]

def plot_convex_set(ax: Axes, covector_to_vector_func, color=None):
    n = 1000
    points = approximate_polygon(covector_to_vector_func, n)
    if color is not None:
        ax.fill(points[:,0], points[:,1], alpha=0.5, color=color)
    else:
        ax.fill(points[:,0], points[:,1], alpha=0.5)

# %%
# A = np.array([[2, 0], [0, 1]])
B = np.array([[1, 0], [0, 1]])
# A_func = to_covector_to_vector_func(A)
A_func = lambda x: np.sign(x) * np.sum(np.abs(x), axis=1)[:,None]
B_func = to_covector_to_vector_func(B)
sum_AB_func = sum(A_func, B_func)
fig, ax = plt.subplots()
plot_convex_set(ax, A_func)
plot_convex_set(ax, B_func)
plot_convex_set(ax, sum_AB_func, 'g')

# %%
def to_vector_to_covector_func(dual_ellipse: np.ndarray):
    return lambda x: (x @ dual_ellipse)

def intersection(vector_to_covector_func1, vector_to_covector_func2):
    def result_func(x):
        covector1 = vector_to_covector_func1(x)
        covector2 = vector_to_covector_func2(x)
        length1 = np.sqrt(np.sum(covector1 * x, axis=1))[:,None]
        length2 = np.sqrt(np.sum(covector2 * x, axis=1))[:,None]

        return np.where(length1 > length2, covector1, covector2)

    return result_func

def approximate_dual_polygon(vector_to_covector_func, n: int):
    angles = np.linspace(0, 2*np.pi, n)
    vectors = np.array([np.cos(angles), np.sin(angles)]).T
    return vectors / np.sqrt(np.sum(vector_to_covector_func(vectors) * vectors, axis=1))[:,None]

def plot_dual_convex_set(ax: Axes, vector_to_covector_func, color=None):
    n = 1000
    points = approximate_dual_polygon(vector_to_covector_func, n)
    if color is not None:
        ax.fill(points[:,0], points[:,1], alpha=0.5, color=color)
    else:
        ax.fill(points[:,0], points[:,1], alpha=0.5)
# %%
# A = np.array([[2, 0], [0, 1]])
B = np.array([[1.5, 0], [0, 2]])
# A_func = to_vector_to_covector_func(A)
A_func = lambda x: np.sign(x) * np.sum(np.abs(x), axis=1)[:,None]
B_func = to_vector_to_covector_func(B)
intersection_AB_func = intersection(A_func, B_func)
fig, ax = plt.subplots()
plot_dual_convex_set(ax, A_func)
plot_dual_convex_set(ax, B_func)
plot_dual_convex_set(ax, intersection_AB_func)
ax.set_aspect('equal')
# %%
