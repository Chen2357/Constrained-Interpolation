# %%
import numpy as np
import matplotlib.pyplot as plt
from ipywidgets.widgets import interact
from matplotlib.axes import Axes

# %%
def _gradient_descent_inverse(func, step=0.01, tol=1e-3):
    def result(y):
        x = y

        # REVIEW: should we set a maximum number of iterations?
        for _ in range(100):
            value, gradient = func(x)
            diff = value - y
            loss = np.sum(diff**2)
            if loss < tol:
                break
            x = x - step * np.einsum("ij,ijk->ik", diff, gradient)

        return (x, np.linalg.inv(gradient))

    return result


class ConvexSet:
    def __init__(self, boundary, inverse_boundary=None):
        """
        boundary: (np.ndarray) -> (np.ndarray, np.ndarray)

        boundary maps vectors to covectors and gradients
        """
        self.boundary = boundary

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

            result = vectors * np.sum(x * vectors, axis=1)[:, np.newaxis]

            vector_grad = \
                inverse_gradient1 / length1[:, np.newaxis, np.newaxis] \
                + inverse_gradient2 / length2[:, np.newaxis, np.newaxis] \
                - vector1[:, :, np.newaxis] * vector1[:, np.newaxis, :] \
                / (2 * length1**3)[:, np.newaxis, np.newaxis] \
                - vector2[:, :, np.newaxis] * vector2[:, np.newaxis, :] \
                / (2 * length2**3)[:, np.newaxis, np.newaxis]

            inverse_gradient = \
                vector_grad * np.sum(x * vectors, axis=1)[:, np.newaxis, np.newaxis] \
                + vectors[:, :, np.newaxis] * vectors[:, np.newaxis, :] \
                + vectors[:, :, np.newaxis] * x[:, np.newaxis, :] @ vector_grad

            return (result, inverse_gradient)

        return ConvexSet(_gradient_descent_inverse(result_inverse), result_inverse)

    def plot(self, ax: Axes, n=1000, **kwargs):
        angles = np.linspace(0, 2*np.pi, n)
        directions = np.array([np.cos(angles), np.sin(angles)]).T
        length = np.sqrt(
            np.sum(directions * self.boundary(directions)[0], axis=1))
        points = directions / length[:, np.newaxis]

        ax.fill(points[:, 0], points[:, 1], **kwargs)


# %%
disk1 = ConvexSet(
    lambda x: (x, np.repeat(np.eye(2)[np.newaxis], len(x), axis=0)),
    lambda x: (x, np.repeat(np.eye(2)[np.newaxis], len(x), axis=0))
)

disk2 = ConvexSet(
    lambda x: (
        x @ np.diag([2, 0.5]), np.repeat(np.diag([2, 0.5])[np.newaxis], len(x), axis=0)),
    lambda x: (
        x @ np.diag([1/2, 2]), np.repeat(np.diag([1/2, 2])[np.newaxis], len(x), axis=0))
)

fig, ax = plt.subplots()
disk1.plot(ax, alpha=0.5)
disk2.plot(ax, alpha=0.5)
disk1.intersection(disk2).plot(ax, alpha=0.5)
ax.set_aspect('equal')
# %%
