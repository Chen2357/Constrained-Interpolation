import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

import whitney as wit
import numpy as np
import numpy.typing as npt

import plotly.graph_objects as go
from plotly.graph_objs import Figure

def disambiguate_paris(pairs: list[list[wit.Hypercube]]):
    i = 0
    while i < len(pairs):
        for j in range(i+1, len(pairs)):
            if (pairs[i][1] == pairs[j][0] and pairs[i][0] == pairs[j][1]) or pairs[i] == pairs[j]:
                del pairs[j]
                break
        i += 1
    return pairs

def sample_points(count: int, type: str = "random"):
    """Types: random, clusters, worst, circle, L"""

    if type == "random":
        return np.random.rand(count, 2)

    if type == "uniform":
        width = np.floor(np.sqrt(count))
        if width**2 != count:
            raise ValueError("uniform distribution must have square number of points")
        gap = 1 / (width + 1)
        return np.array([(i*gap, j*gap) for i in np.arange(1, (width+1)) for j in np.arange(1, (width+1))])

    if type == "clusters":
        return np.concatenate(
            (
                np.random.rand(count//3,2) * 0.1 + 0.1,
                np.random.rand(count//3,2) * 0.1 + 0.9,
                np.random.rand(count//3 + count % 3,2) * 0.1 + np.array([0.1, 0.9])
            ), axis=0
        )

    if type == "worst":
        array = 2.0 ** np.arange(-1, -count-1, -1)
        return np.vstack((array, array)).T

    if type == "circle":
        return np.array([0.5, 0.5]) + 0.49 * np.array([np.cos(np.linspace(0, 2*np.pi, count, endpoint=False)), np.sin(np.linspace(0, 2*np.pi, count, endpoint=False))]).T

    if type == "L":
        vertical_count = (count - 1) // 2
        horizontal_count = count - 1 - vertical_count

        vertical = np.array([np.repeat(0.01, vertical_count), np.linspace(0.4, 0.9, vertical_count)]).T
        horizontal = np.array([np.linspace(0.4, 0.9, horizontal_count), np.repeat(0.01, horizontal_count)]).T
        return np.concatenate((vertical, horizontal, [np.array([0.01, 0.01])]), axis=0)

    raise ValueError("Unknown type")

def sample_polynomials(points: npt.NDArray, degree):
    """Types: random"""
    if len(points) == 0:
        return {}
    polynomials: wit.Polynomials = {}
    dimension = len(points[0])
    for point in points:
        polynomials[point.tobytes()] = zero_polynomial(dimension, degree)

    return polynomials

def zero_polynomial(dimension, degree):
    array = np.zeros(np.repeat(degree + 1, dimension))
    return array

def approximate_polygon(boundary_func, n: int):
    angles = np.linspace(0, 2*np.pi, n)
    directions = np.array([np.cos(angles), np.sin(angles)]).T
    points = boundary_func(directions) * directions.T

    return points.T

def to_boundary_func(elipse: np.ndarray):
    def boundary_func(directions):
        result = []
        for direction in directions:
            result.append(np.sqrt((direction @ direction) / (direction @ elipse @ direction)))
        return np.array(result)

    return boundary_func

def plot_sigma_at(fig: Figure, sigma, x, n=1000):
    # pullback_sigma = wit._pullback(sigma, wit._forward_transformation(-x))
    points = approximate_polygon(to_boundary_func(sigma[1:,1:]), n)
    points = points + x
    fig.add_trace(go.Scatter(x=points[:,0], y=points[:,1], fill="toself"))

def _plot_sigma_at(ax, sigma, x):
    # pullback_sigma = wit._pullback(sigma, wit._forward_transformation(-x))
    n = 1000
    points = approximate_polygon(to_boundary_func(sigma[1:,1:]), n)
    points = points + x
    ax.fill(points[:,0], points[:,1], alpha=0.5)
