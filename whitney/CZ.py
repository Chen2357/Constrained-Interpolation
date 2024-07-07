import numpy as np
import numpy.typing as npt
from .Convex import _pullback, _forward_transformation, scale, sum, _inv_john_ellipsoid, intersection, _sum_with_inverse
from .Hypercube import Hypercube
import queue
from .WSPD import build_wspd


from scipy import spatial

def _diam_inf(points):
    try:
        candidates = points[spatial.ConvexHull(points).vertices]
    except:
        return np.max(points[:, 0]) - np.min(points[:, 0]) + np.max(points[:, 1]) - np.min(points[:, 1])

    dist_mat = spatial.distance_matrix(candidates, candidates, p=float('inf')) # type: ignore
    return np.sqrt(2) * np.max(dist_mat)

def CZ_decompose(root: Hypercube, a: float = 30, s: float = 2, thickness = 0.001, N = 6, T = 10, post_shrinking = 0.25):
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
    decomposition = CZ_Decomposition(root, a, s, thickness, N, T, post_shrinking)
    decomposition._CZ_decompose()
    return decomposition

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

    def _compute_group_sigma(self, sigma):
        # For each A in T, define sigma(A) = intersection(sigma(x) + B(x, diam(A)) for x in A)
        return [
            intersection([
                _sum_with_inverse(
                    sigma[j],
                    self._ball_inverse(_diam_inf(self.points[self._groups[i]]), self.points[j])
                )
                for j in self._groups[i]
            ])
            for i in range(len(self._groups))
        ]

    def _compute_sigma_temp(self):
        return [
            _sum_with_inverse(
                self._group_sigma[j],
                self._ball_inverse(_diam_inf(self.points[self._groups[j]]), self.points[np.random.choice(self._groups[j])])
            )
            for j in range(len(self._groups))
        ]

    def _compute_sigma_bar(self, sigma_temp):
        # For each (A, B) in L,
        #  Define sigma_1(A) = sigma(A) + B(x_A, diam(A))
        #  Define sigma_2(B) = sigma(B) + B(x_B, diam(B))
        # For each (A, B) in L, define sigma(A, B) = intersection(sigma_1(A), sum(sigma_2(B) + B(x_A, |x_A - x_B|)))
        sigma_bar = [[np.empty(0, dtype=float) for _ in range(len(self._groups))] for _ in range(len(self._groups))]
        # print(np.shape(sigma_bar))

        for j, k in self._well_separated_pairs_indices:
            sigma_bar[j][k] = intersection([
                sigma_temp[j],
                _sum_with_inverse(
                    sigma_temp[k],
                    self._ball_inverse(
                        np.linalg.norm(self.points[self._groups[j][0]] - self.points[self._groups[k][0]], ord=np.inf),
                        self.points[self._groups[j][0]]
                    )
                )
            ])

        return sigma_bar

    def _compute_sigma_prime(self, sigma_bar):
        # For each A in T, define sigma'(A) = intersection(sigma_bar(A, B) where (A, B) in L)
        intersectands = [[] for _ in range(len(self._groups))]

        for j, k in self._well_separated_pairs_indices:
            intersectands[j].append(sigma_bar[j][k])

        sigma_prime = [intersection(intersectands[i]) for i in range(len(self._groups))]

        return sigma_prime

    def _compute_new_sigma(self, sigma, sigma_prime):
        # For each x in E, redefine sigma(x) = simga(x) intersect intersection(sigma(A) for A in T where x in A)
        new_sigma = [np.empty(0, dtype=float) for _ in range(len(self.points))]
        for i in range(len(self.points)):
            intersectands = [sigma[i]]
            for j in range(len(self._groups)):
                if i in self._groups[j]:
                    intersectands.append(sigma_prime[j])
            new_sigma[i] = intersection(intersectands)

        return new_sigma

    def _approximate_sigma(self):
        # Algorithm based on page 92 of Fefferman's paper "Fiiting a Cm-smooth function to data II"
        sigma = np.array([self._sigma_0(x) for x in self.points])
        groups, well_separated_pairs_indices = build_wspd(self.points, self.s)

        self._groups = groups
        self._well_separated_pairs_indices = well_separated_pairs_indices
        self._group_sigma: list[np.ndarray]
        self._sigma_bar: list[list[np.ndarray]]
        self._sigma_prime: list[np.ndarray]

        def recursion(sigma):
            group_sigma = self._compute_group_sigma(sigma)
            self._group_sigma = group_sigma

            sigma_temp = self._compute_sigma_temp()

            sigma_bar = self._compute_sigma_bar(sigma_temp)
            self._sigma_bar = sigma_bar

            sigma_prime = self._compute_sigma_prime(sigma_bar)
            self._sigma_prime = sigma_prime

            new_sigma = self._compute_new_sigma(sigma, sigma_prime)
            return new_sigma

        for _ in range(self.N):
            sigma = recursion(sigma)

        return np.array([_pullback(sigma[i], _forward_transformation(-self.points[i])) for i in range(len(self.points))])