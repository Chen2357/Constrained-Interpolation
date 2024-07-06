import numpy as np
import numpy.typing as npt
from .Convex import _pullback, _forward_transformation, scale, sum, _inv_john_ellipsoid, intersection, _sum_with_inverse
from .Hypercube import Hypercube
import queue
from .WSPD import build_wspd

import timeit


from scipy import spatial

def _diam_inf(points):
    try:
        candidates = points[spatial.ConvexHull(points).vertices]
    except:
        return np.max(points[:, 0]) - np.min(points[:, 0]) + np.max(points[:, 1]) - np.min(points[:, 1])

    dist_mat = spatial.distance_matrix(candidates, candidates, p=float('inf')) # type: ignore
    return np.sqrt(2) * np.max(dist_mat)

class CZ_Decomposition:
    def __init__(self, root: Hypercube, a: float = 30, s: float = 2, thickness = 0.001, N = 6, T = 10, post_shrinking = 0.25):
        """
        root: Hypercube
            The root hypercube of the CZ decomposition
        thickness: float
            The thickness of Ellipsoid
        a: float
            decomposition constant
        s: float
            The separation constant for WSPD
        N: int
            The number of iterations to approximate sigma
        T: int
            The number of iterations to approximate John Ellipsoid
        """
        self.root = root
        self.thickness = thickness
        self.a = a
        self.s = s
        self.N = N
        self.T = T
        self.post_shrinking = post_shrinking

        self.group_sigma_counter = float(0)
        self.sigma_temp_counter = float(0)
        self.sigma_bar_counter = float(0)
        self.sigma_prime_counter = float(0)
        self.new_sigma_counter = float(0)

        self._CZ_decompose()

        

    @property
    def points(self):
        return self.root.points

    def _CZ_decompose(self):
        # Produce finer decomposition

        sigma = scale(self._approximate_sigma(), self.post_shrinking)
        def is_good(square: Hypercube):
            i = self.root.indices_search_in(square.dialated(3))
            diameters = 2 / np.sqrt(np.linalg.eig(sigma[i])[0].min(axis=-1))
            return np.all(diameters >= self.a * square.width)

        q = queue.Queue()
        q.put(self.root)

        while q.qsize() != 0:
            current: Hypercube = q.get()
            if is_good(current):
                continue
            else:
                current.subdivide()
                for child in current.children:
                    q.put(child)

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

    def _ball_inverse(self, delta, x):
        return _pullback(np.array([
            [delta**4, 0, 0],
            [0, delta**2, 0],
            [0, 0, delta**2]
        ]), _forward_transformation(-x).T)

    def _approximate_sigma(self):
        sigma = np.array([self._sigma_0(x) for x in self.points])
        rotated_points = np.array([[1-point[1], point[0]] for point in self.points])
        groups, well_separated_pairs_indices = build_wspd(self.points, self.s)

        self._groups = groups
        self._well_separated_pairs_indices = well_separated_pairs_indices
        self._group_sigma: list[np.ndarray]
        self._sigma_bar: list[list[np.ndarray]]
        self._sigma_prime: list[np.ndarray]

        def recursion(self, sigma):
            # For each A in T, define sigma(A) = intersection(sigma(x) + B(x, diam(A)) for x in A)

            """"""
            group_sigma_time = timeit.default_timer()
            """"""
            group_sigma = [
                intersection([
                    _sum_with_inverse(
                        sigma[j],
                        self._ball_inverse(_diam_inf(self.points[groups[i]]), self.points[j])
                    )
                    for j in groups[i]
                ])
                for i in range(len(groups))
            ]
            self._group_sigma = group_sigma
            """"""
            group_sigma_time_diff = timeit.default_timer() - group_sigma_time
            self.group_sigma_counter = self.group_sigma_counter + group_sigma_time_diff
            """"""

            # For each (A, B) in L,
            #  Define sigma_1(A) = sigma(A) + B(x_A, diam(A))
            #  Define sigma_2(B) = sigma(B) + B(x_B, diam(B))
            # For each (A, B) in L, define sigma(A, B) = intersection(sigma_1(A), sum(sigma_2(B) + B(x_A, |x_A - x_B|)))
            """"""
            sigma_temp_time = timeit.default_timer()
            """"""
            sigma_temp = [
                _sum_with_inverse(
                    group_sigma[j],
                    self._ball_inverse(_diam_inf(self.points[groups[j]]), self.points[np.random.choice(groups[j])])
                )
                for j in range(len(groups))
            ]
            """"""
            sigma_temp_time_diff = timeit.default_timer() - sigma_temp_time
            self.sigma_temp_counter = self.sigma_temp_counter + sigma_temp_time_diff
            """"""

            """"""
            sigma_bar_time = timeit.default_timer()
            """"""
            sigma_bar = [[np.empty(0, dtype=float) for _ in range(len(groups))] for _ in range(len(groups))]

            for j, k in well_separated_pairs_indices:
                sigma_bar[j][k] = intersection([
                    sigma_temp[j],
                    _sum_with_inverse(
                        sigma_temp[k],
                        self._ball_inverse(
                            np.linalg.norm(self.points[groups[j][0]] - self.points[groups[k][0]], ord=np.inf),
                            self.points[groups[j][0]]
                        )
                    )
                ])

            self._sigma_bar = sigma_bar
            """"""
            sigma_bar_time_diff = timeit.default_timer() - sigma_bar_time
            self.sigma_bar_counter = self.sigma_bar_counter + sigma_bar_time_diff
            """"""

            # For each A in T, define sigma'(A) = intersection(sigma_bar(A, B) where (A, B) in L)
            """"""
            sigma_prime_time = timeit.default_timer()
            """"""
            sigma_prime = [np.empty(0, dtype=float) for _ in range(len(groups))]
            for i in range(len(groups)):
                intersectands = []
                for j in range(len(groups)):
                    if [i, j] in well_separated_pairs_indices:
                        intersectands.append(sigma_bar[i][j])

                sigma_prime[i] = intersection(intersectands)

            self._sigma_prime = sigma_prime
            """"""
            sigma_prime_time_diff = timeit.default_timer() - sigma_prime_time
            self.sigma_prime_counter = self.sigma_prime_counter + sigma_prime_time_diff
            """"""

            # For each x in E, redefine sigma(x) = simga(x) intersect intersection(sigma(A) for A in T where x in A)
            """"""
            new_sigma_time = timeit.default_timer()
            """"""
            new_sigma = [np.empty(0, dtype=float) for _ in range(len(self.points))]
            for i in range(len(self.points)):
                intersectands = [sigma[i]]
                for j in range(len(groups)):
                    if i in groups[j]:
                        intersectands.append(sigma_prime[j])
                new_sigma[i] = intersection(intersectands)
                
            """"""
            new_sigma_time_diff = timeit.default_timer() - new_sigma_time
            self.new_sigma_counter = self.new_sigma_counter + new_sigma_time_diff
            """"""
            return new_sigma

        for _ in range(self.N):
            sigma = recursion(self, sigma)

        return np.array([_pullback(sigma[i], _forward_transformation(-self.points[i])) for i in range(len(self.points))])
        # sigma = np.array([self._sigma_0(x) for x in self.points])

        # def recursion(sigma):
        #     new_sigma = sigma
        #     for i in range(len(self.points)):
        #         x = self.points[i]

        #         intersectands = [sigma[i]]

        #         for j in range(len(self.points)):
        #             if i == j:
        #                 continue
        #             y = self.points[j]
        #             distance = np.linalg.norm(x - y, ord=2)

        #             intersectands.append(sum(sigma[j], scale(self._ball(distance, x), C)))

        #         new_sigma[i] = intersection(intersectands)

        #     return new_sigma

        # for _ in range(6):
        #     sigma = recursion(sigma)

        # return np.array([_pullback(sigma[i], _forward_transformation(-self.points[i])) for i in range(len(self.points))])