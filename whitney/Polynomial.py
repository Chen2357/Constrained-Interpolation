import numpy as np
import numpy.typing as npt
from typing import Union, Dict
from .Hypercube import Hypercube

Polynomials = Dict["bytes", npt.NDArray]

def evaluate_polynomial(polynomial: npt.NDArray, point):
    degree_plus_1 = len(polynomial)

    result = polynomial
    for x in point:
        jet = [x ** n for n in range(degree_plus_1)]
        result = np.tensordot(jet, result, axes = (0,0))

    return result

def _assign_jet(root: Hypercube, cube: Hypercube, jets: npt.NDArray):
    # points = root.search_in(Hypercube(cube.pos - cube.width, 3*cube.width, 1))
    # if len(points) == 0:
    #     points = root.search_in(Hypercube(cube.parent.pos - cube.parent.width, 3*cube.parent.width, 1))

    center = cube.pos + cube.width / 2
    point = root.query_nearest_point(center)
    # i = np.argmin(np.sum((points - center)**2, axis=1))
    # index = np.where(root.points == point)[0][0]
    index = np.where(np.all(root.points == point, axis=1))[0][0]
    cube.pointed_jet = np.concatenate([point, jets[index]])

def assign_jets(root: Hypercube, jets: npt.ArrayLike):
    jets = np.array(jets)
    for leaf in root.leaves:
        _assign_jet(root, leaf, jets)

def evaluate_cube(cube: Hypercube, point: npt.ArrayLike):
    anchor = cube.pointed_jet[0:cube.dimension]
    constant = cube.pointed_jet[cube.dimension]
    gradient = cube.pointed_jet[(cube.dimension+1):]

    gradient_squared = np.sum(gradient**2)

    if gradient_squared < 1e-7:
        return constant
    else:
        a = gradient_squared / (4 * constant)
        x0 = anchor - 2*constant/(gradient_squared) * gradient

        return a*np.sum((point - x0)**2)