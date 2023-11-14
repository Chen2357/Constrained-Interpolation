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


# NOTE This is a temporary solution to the problem of evaluating polynomials.
import sympy as sp

def find_orientation(points: npt.ArrayLike) -> float:
    points = np.array(points)
    if np.shape(points) != (3, 2):
        raise ValueError("Only 3 2D-points are supported")

    A = points[1] - points[0]
    B = points[2] - points[1]
    C = points[0] - points[2]

    x1, y1, x2, y2, x3, y3 = sp.symbols("x_1, y_1, x_2, y_2, x_3, y_3")

    theta = sp.symbols("\\theta")
    z = sp.symbols("z")

    expr: sp.Symbol = (y1*x1**2 + y2*x2**2 + y3*x3**2).subs({
        x1: A[0]*sp.cos(theta) - A[1]*sp.sin(theta),
        y1: A[0]*sp.sin(theta) + A[1]*sp.cos(theta),
        x2: B[0]*sp.cos(theta) - B[1]*sp.sin(theta),
        y2: B[0]*sp.sin(theta) + B[1]*sp.cos(theta),
        x3: C[0]*sp.cos(theta) - C[1]*sp.sin(theta),
        y3: C[0]*sp.sin(theta) + C[1]*sp.cos(theta),
    }).subs({
        sp.cos(theta): (1-z**2)/(1+z**2),
        sp.sin(theta): 2*z/(1+z**2),
    }).expand() * ((1+z**2)**3)

    z_roots = sp.nroots(expr.simplify())
    theta_roots = [2*sp.atan(z) for z in z_roots]

    best_rad = 0
    min_curvature = np.inf

    for rad in theta_roots:
        rad = float(rad)
        rotation = np.array([
            [np.cos(rad), -np.sin(rad)],
            [np.sin(rad), np.cos(rad)]
        ])

        rotated_points = points @ rotation

        num = np.linalg.det([
            [rotated_points[0,1], rotated_points[0,0], 1],
            [rotated_points[1,1], rotated_points[1,0], 1],
            [rotated_points[2,1], rotated_points[2,0], 1]
        ])
        denom = np.linalg.det([
            [rotated_points[0,0]**2, rotated_points[0,0], 1],
            [rotated_points[1,0]**2, rotated_points[1,0], 1],
            [rotated_points[2,0]**2, rotated_points[2,0], 1]
        ])

        curvature = np.abs(num / denom)
        if curvature < min_curvature:
            min_curvature = curvature
            best_rad = rad

    return best_rad

