# %%
import numpy as np
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
from ipywidgets.widgets import interact

from typing import Callable, Tuple
import numpy.typing as npt

from copy import deepcopy

# %%
def _gradient_descent_inverse(func, step=0.001, beta=0.9, tol=1e-3):
    def result(y):
        x = y
        previous_change = None

        # REVIEW: should we set a maximum number of iterations?
        for _ in range(100000):
            value, gradient = func(x)
            diff = value - y

            if np.max(np.abs(diff)) < tol:
                return (x, np.linalg.inv(gradient))

            if previous_change is not None:
                change = np.einsum("ij,ijk->ik", diff, gradient) + beta * previous_change
            else:
                change = np.einsum("ij,ijk->ik", diff, gradient)

            x = x - step * change
            previous_change = change

        raise ValueError(f"Gradient descent did not converge")

    return result

class ConvexSet:
    def __init__(self,
                 boundary: Callable[[npt.NDArray], Tuple[npt.NDArray, npt.NDArray]],
                 inverse_boundary: Callable[[npt.NDArray], Tuple[npt.NDArray, npt.NDArray]] | None = None):
        """
        boundary: (np.ndarray) -> (np.ndarray, np.ndarray)

        boundary maps vectors to covectors and gradients
        """

        self.boundary = boundary

        self.inverse_boundary: Callable[[npt.NDArray], Tuple[npt.NDArray, npt.NDArray]]
        if inverse_boundary is None:
            self.inverse_boundary = _gradient_descent_inverse(self.boundary)
        else:
            self.inverse_boundary = inverse_boundary

    def intersection(self, other: 'ConvexSet'):
        def result(x):
            covector1, gradient1 = self.boundary(x)
            covector2, gradient2 = other.boundary(x)
            length1 = np.sqrt(np.sum(covector1 * x, axis=1))
            length2 = np.sqrt(np.sum(covector2 * x, axis=1))

            return (
                np.where((length1 > length2)[:, np.newaxis], covector1, covector2),
                np.where((length1 > length2)[:, np.newaxis, np.newaxis], gradient1, gradient2)
            )

        return ConvexSet(result)

    def sum(self, other: 'ConvexSet'):
        def result_inverse(x):
            vector1, inverse_gradient1 = self.inverse_boundary(x)
            vector2, inverse_gradient2 = other.inverse_boundary(x)
            length1 = np.sqrt(np.sum(x * vector1, axis=1))
            length2 = np.sqrt(np.sum(x * vector2, axis=1))

            vectors = vector1 / length1[:, np.newaxis] + vector2 / length2[:, np.newaxis]

            lamb = np.sum(x * vectors, axis=1)
            result = vectors * lamb[:, np.newaxis]

            length1_grad = vector1 + np.einsum("ij,ijk->ik", x, inverse_gradient1)
            length2_grad = vector2 + np.einsum("ij,ijk->ik", x, inverse_gradient2)

            vector_grad = \
                inverse_gradient1 / length1[:, np.newaxis, np.newaxis] \
                + inverse_gradient2 / length2[:, np.newaxis, np.newaxis] \
                - vector1[:, :, np.newaxis] * length1_grad[:, np.newaxis, :] \
                / (2 * length1**3)[:, np.newaxis, np.newaxis] \
                - vector2[:, :, np.newaxis] * length2_grad[:, np.newaxis, :] \
                / (2 * length2**3)[:, np.newaxis, np.newaxis]

            lamb_grad = vectors + np.einsum("ij,ijk->ik", x, vector_grad)

            inverse_gradient = \
                vector_grad * lamb[:, np.newaxis, np.newaxis] \
                + vectors[:, :, np.newaxis] * lamb_grad[:, np.newaxis, :]

            return (result, inverse_gradient)

        return ConvexSet(_gradient_descent_inverse(result_inverse), result_inverse)

    def plot(self, ax: Axes, n=1000, **kwargs):
        angles = np.linspace(0, 2*np.pi, n)
        directions = np.array([np.cos(angles), np.sin(angles)]).T
        length = np.sqrt(
            np.sum(directions * self.boundary(directions)[0], axis=1))
        points = directions / length[:, np.newaxis]

        ax.fill(points[:, 0], points[:, 1], **kwargs)

    def _check_integrity(self, n, rtol=1e-3):
        id = np.eye(n)
        N = 100
        x = np.random.rand(100, n)
        ids = np.repeat(id[np.newaxis], 100, axis=0)

        covector, gradient = self.boundary(x)
        vector, inverse_gradient = self.inverse_boundary(covector)
        assert np.allclose(vector, x, rtol = rtol)
        assert np.allclose(np.einsum("ijk,ikl->ijl", gradient, inverse_gradient), ids, rtol = rtol)
# %%

def forward_transformation(x: npt.NDArray):
    return np.array([
        [1, x[0], x[1]],
        [0, 1, 0],
        [0, 0, 1]
    ])

def pushforward(convex_set: ConvexSet, transformation: np.ndarray) -> ConvexSet:
    def result_boundary(x):
        inv_transformation = np.linalg.inv(transformation)
        covector, gradient = convex_set.boundary(x @ inv_transformation.T)

        return (
            covector @ inv_transformation,
            np.einsum("jm,ijk,kn->imn", inv_transformation, gradient, inv_transformation)
        )

    def result_inverse_boundary(x):
        vector, inverse_gradient = convex_set.inverse_boundary(x @ transformation)

        return (
            vector @ transformation.T,
            np.einsum("mj,ijk,nk->imn", transformation, inverse_gradient, transformation)
        )

    return ConvexSet(
        result_boundary,
        result_inverse_boundary
    )

def pullback(convex_set: ConvexSet, transformation: npt.NDArray) -> ConvexSet:
    def result_boundary(x):
        covector, gradient = convex_set.boundary(x @ transformation.T)

        return (
            covector @ transformation,
            np.einsum("jm,ijk,kn->imn", transformation, gradient, transformation)
        )

    def result_inverse_boundary(x):
        inv_transformation = np.linalg.inv(transformation)
        vector, inverse_gradient = convex_set.inverse_boundary(x @ inv_transformation)

        return (
            vector @ inv_transformation.T,
            np.einsum("mj,ijk,nk->imn", inv_transformation, inverse_gradient, inv_transformation)
        )

    return ConvexSet(
        result_boundary,
        result_inverse_boundary
    )

def convex_set_from_ellipsoid(ellipsoid: npt.NDArray) -> ConvexSet:
    def result_boundary(x):
        return (x @ ellipsoid, np.repeat(ellipsoid[np.newaxis], len(x), axis=0))

    inv = np.linalg.inv(ellipsoid)
    def result_inverse_boundary(x):
        return (x @ inv, np.repeat(inv[np.newaxis], len(x), axis=0))

    return ConvexSet(result_boundary, result_inverse_boundary)

def ball(delta: float, x: npt.NDArray) -> ConvexSet:
    ellipsoid = np.diag([1/delta**2, 1/delta**2, 1/delta**2])
    convest_set = convex_set_from_ellipsoid(ellipsoid)
    return pullback(convest_set, forward_transformation(x))

thickness = 0.1

def sigma_0(x: npt.NDArray) -> ConvexSet:
    ellipsoid = np.diag([1/thickness**2, 1, 1])
    convest_set = convex_set_from_ellipsoid(ellipsoid)
    return pullback(convest_set, forward_transformation(x))

def scale(convex_set: ConvexSet, c: float) -> ConvexSet:
    def result_boundary(x):
        covector, gradient = convex_set.boundary(x / c)

        return (
            covector / c,
            gradient / c**2
        )

    def result_inverse_boundary(x):
        vector, inverse_gradient = convex_set.inverse_boundary(x * c)

        return (
            vector * c,
            inverse_gradient * c**2
        )

    return ConvexSet(
        result_boundary,
        result_inverse_boundary
    )

def norm(x: npt.NDArray):
    return np.linalg.norm(x, ord=np.inf)

def nearest_k_points(points: npt.NDArray, k: int):
    """
    Returns an n by k array of the k nearest points to each point in the n by d array points
    """
    distances = np.sum(points**2, axis=1).reshape(-1, 1) - 2 * points @ points.T + np.sum(points**2, axis=1).reshape(1, -1)
    return np.argsort(distances, axis=1)[:, 1:k+1]

def to_boundary_func(convex_set: ConvexSet) -> Callable[[npt.NDArray], npt.NDArray]:
    def boundary_func(directions: npt.NDArray):
        n = len(directions)
        directions = np.concatenate([np.zeros([n, 1]), directions], axis=1)

        covector, _ = convex_set.boundary(directions)
        return directions[:,1:] / np.sqrt(np.sum(covector * directions, axis=1))[:, np.newaxis]

    return boundary_func

def approximate_polygon(boundary_func: Callable[[npt.NDArray], npt.NDArray], n: int):
    angles = np.linspace(0, 2*np.pi, n)
    directions = np.array([np.cos(angles), np.sin(angles)]).T
    return boundary_func(directions)

def plot_boundary_func(ax: Axes, boundary_func: Callable[[npt.NDArray], npt.NDArray], color=None):
    n = 100
    points = approximate_polygon(boundary_func, n)
    if color is not None:
        ax.fill(points[:,0], points[:,1], alpha=0.5, color=color)
    else:
        ax.fill(points[:,0], points[:,1], alpha=0.5)
# %%
C = 2

dim = 2

parabola = np.array([[x, x**2] for x in np.linspace(-1, 1, 40)])

E = np.array([
    [0.01, 0.8],
    [0.01, 0.6],
    [0.01, 0.4],
    [0.01, 0.01],
    [0.4, 0.01],
    [0.6, 0.01],
    [0.8, 0.01],
])

closest_points_info = nearest_k_points(E, 2)

sigma = [sigma_0(x) for x in E]
sigma_history = [deepcopy(sigma)]

debug_sigma_element: ConvexSet

def recrusion(sigma: list[ConvexSet]):
    new_sigma = sigma
    for i in range(len(E)):
        x = E[i]

        for j in range(len(E)):
            if i == j:
                continue
            y = E[j]
            distance = float(norm(x-y))

            # try:
            #     new_sigma[j].sum(scale(ball(distance, x), C))._check_integrity(3)
            # except:
            #     debug_sigma_element = new_sigma[j].sum(scale(ball(distance, x), C))
            #     raise ValueError("Integerity check failed")

            new_sigma[i] = new_sigma[i].intersection(new_sigma[j].sum(scale(ball(distance, x), C)))
            # try:
            #     new_sigma[i]._check_integerity(3)
            # except:
            #     print(i, E[i])
            #     print(j, E[j])
            #     raise ValueError("Integerity check failed")

    return new_sigma

for i in range(1):
    try:
        sigma = recrusion(sigma)
        sigma_history.append(deepcopy(sigma))
    except:
        print("Failed at", i)

# %%

def plot_convex_set_at(ax: Axes, convex_set: ConvexSet, x: npt.NDArray):
    pullback_convex_set = pushforward(convex_set, forward_transformation(x))
    n = 1000
    points = approximate_polygon(to_boundary_func(pullback_convex_set), n)
    points = points + x
    ax.fill(points[:,0], points[:,1], alpha=0.5)

# @interact(j=(0, len(sigma_history)-1))
def plot_at(j = 1):
    fig, ax = plt.subplots()
    ax: Axes
    print(pullback(sigma_history[j][0], forward_transformation(-E[0])))
    for i in range(len(E)):
        plot_convex_set_at(ax, sigma_history[j][i], E[i])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
# plot_at()
# %%
