import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

import whitney as wit
from test_module import *

import numpy as np
import sympy as sp

def find_orientation_symbolic(points: npt.ArrayLike) -> float:
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

n = 100
points = sample_points(n+2)

for i in range(n):
    symbolic_result = find_orientation_symbolic(points[i:i+3]) % np.pi
    result = wit.find_orientation(points[i:i+3]) % np.pi
    if not np.isclose(symbolic_result, result):
        print(f"Symbolic result: {symbolic_result}")
        print(f"Result: {result}")
        print(f"Points: {points[i:i+3]}")
        raise ValueError("Symbolic result and result are not close")

    print(f"Test {i+1}/{n} passed")

print("All tests passed")