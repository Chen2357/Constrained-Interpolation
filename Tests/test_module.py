import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

import whitney as wit
import numpy as np
import numpy.typing as npt

def disambiguate_paris(pairs: list[list[wit.Hypercube]]):
    i = 0
    while i < len(pairs):
        for j in range(i+1, len(pairs)):
            if (pairs[i][1] == pairs[j][0] and pairs[i][0] == pairs[j][1]) or pairs[i] == pairs[j]:
                del pairs[j]
                break
        i += 1
    return pairs

def sample_points(count: int, type: str):
    """Types: random, clusters, worst"""

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
