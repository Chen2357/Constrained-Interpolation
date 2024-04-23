import numpy as np
import numpy.typing as npt

def _forward_transformation(x: npt.NDArray):
    return np.block([[1, x], [np.zeros((len(x), 1)), np.eye(len(x))]])

def _pullback(ellipsoid: npt.NDArray, forward: npt.NDArray):
    return forward.T @ ellipsoid @ forward

def scale(ellipsoid, c):
    return ellipsoid / c**2

def sum(ellipsoid1: npt.NDArray, ellipsoid2: npt.NDArray):
    return np.linalg.inv(np.linalg.inv(ellipsoid1) + np.linalg.inv(ellipsoid2))

def _inv_john_ellipsoid(points: npt.NDArray, T = 10):
    m, n = points.shape
    w = np.repeat(n / m, m)

    for k in range(0, T-1):
        g = np.linalg.inv(points.T @ np.diag(w) @ points)
        w = w * np.einsum('ij,jk,ik->i', points, g, points)

    return points.T @ np.diag(w) @ points

def intersection(ellipsoids: npt.NDArray, T = 10):
    L = np.linalg.cholesky(ellipsoids)
    points = np.concatenate(np.swapaxes(L, 1, 2))
    return _inv_john_ellipsoid(points, T)

