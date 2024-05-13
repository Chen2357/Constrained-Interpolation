from test_module import *
import numpy as np

class Data:
    def __init__(self, root: wit.Hypercube, thickness = 0.001):
        self.root = root
        self.thickness = thickness

    @property
    def points(self):
        return self.root.points

    def _sigma_0(self, x):
        return wit._pullback(np.array([
            [1/self.thickness**2, 0, 0],
            [0, 1, 0],
            [0, 0, 1]
        ]), wit._forward_transformation(x))

    def _ball(self, delta, x):
        return wit._pullback(np.array([
            [1/delta**4, 0, 0],
            [0, 1/delta**2, 0],
            [0, 0, 1/delta**2]
        ]), wit._forward_transformation(x))

    def _approximate_sigma(self, C = 1):
        sigma = np.array([self._sigma_0(x) for x in self.points])

        def recursion(sigma):
            new_sigma = sigma
            for i in range(len(self.points)):
                x = self.points[i]

                intersectands = [new_sigma[i]]

                for j in range(len(self.points)):
                    if i == j:
                        continue
                    y = self.points[j]
                    distance = np.linalg.norm(x - y, ord=2)

                    intersectands.append(sum(new_sigma[j], wit.scale(self._ball(distance, x), C)))

                new_sigma[i] = wit.intersection(intersectands)

            return new_sigma

        for _ in range(6):
            sigma = recursion(sigma)

        return np.array([wit._pullback(sigma[i], wit._forward_transformation(-self.points[i])) for i in range(len(self.points))])

E = np.array([
    [0.01, 0.9],
    [0.01, 0.85],
    [0.01, 0.8],
    [0.01, 0.75],
    [0.01, 0.7],
    [0.01, 0.65],
    [0.01, 0.6],
    [0.01, 0.55],
    [0.01, 0.5],
    [0.01, 0.45],
    [0.01, 0.4],
    [0.01, 0.01],
    [0.4, 0.01],
    [0.45, 0.01],
    [0.5, 0.01],
    [0.55, 0.01],
    [0.6, 0.01],
    [0.65, 0.01],
    [0.7, 0.01],
    [0.75, 0.01],
    [0.8, 0.01],
    [0.85, 0.01],
    [0.9, 0.01],
])

root = wit.Hypercube((0, 0), 1, E)
Data(root)._approximate_sigma()