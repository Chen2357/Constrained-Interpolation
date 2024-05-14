import numpy as np
import numpy.typing as npt
from .Convex import _pullback, _forward_transformation, scale, sum, _inv_john_ellipsoid, intersection
from .Hypercube import Hypercube
from .Utility import asvoid, find_index
import queue

class CZ_Decomposition:
    def __init__(self, root: Hypercube, a = 30, thickness = 0.001, N = 6, T = 10):
        """
        root: Hypercube
            The root hypercube of the CZ decomposition
        thickness: float
            The thickness of Ellipsoid
        a: int
            decomposition constant
        N: int
            The number of iterations to approximate sigma
        T: int
            The number of iterations to approximate John Ellipsoid
        """
        self.root = root
        self.thickness = thickness
        self.a = a
        self.N = N
        self.T = T
        self._CZ_decompose()

    @property
    def points(self):
        return self.root.points

    def _CZ_decompose(self, a = 30):
        # root = wit.Hypercube(np.array([0, 0]), 1, points)
        sigma = self._approximate_sigma()

        def indices_of_points(p):
            return [find_index(self.points, q)[0] for q in p]

        def is_good(square: Hypercube):
            p = self.root.search_in(square.dialated(3))
            i = indices_of_points(p)
            diameters = 2 / np.sqrt(np.linalg.eig(sigma[i])[0].min(axis=-1))
            return np.all(diameters >= self.a * square.width)

        q = queue.Queue()
        q.put(self.root)

        while q.qsize() != 0:
            current = q.get()
            current: Hypercube
            if is_good(current):
                continue
            else:
                current.subdivide()
                for child in current.children:
                    q.put(child)

        # return root

    def _sigma_0(self, x):
        return _pullback(np.array([
            [1/self.thickness**2, 0, 0],
            [0, 1, 0],
            [0, 0, 1]
        ]), _forward_transformation(x))

    def _ball(self, delta, x):
        return _pullback(np.array([
            [1/delta**4, 0, 0],
            [0, 1/delta**2, 0],
            [0, 0, 1/delta**2]
        ]), _forward_transformation(x))

    def _approximate_sigma(self, C = 1):
        sigma = np.array([self._sigma_0(x) for x in self.points])

        def recursion(sigma):
            new_sigma = sigma
            for i in range(len(self.points)):
                x = self.points[i]

                intersectands = [sigma[i]]

                for j in range(len(self.points)):
                    if i == j:
                        continue
                    y = self.points[j]
                    distance = np.linalg.norm(x - y, ord=2)

                    intersectands.append(sum(sigma[j], scale(self._ball(distance, x), C)))

                new_sigma[i] = intersection(intersectands)

            return new_sigma

        for _ in range(6):
            sigma = recursion(sigma)

        return np.array([_pullback(sigma[i], _forward_transformation(-self.points[i])) for i in range(len(self.points))])