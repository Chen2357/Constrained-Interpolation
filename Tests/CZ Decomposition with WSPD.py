from sympy import group
from test_module import *
from actual_CZ import *
import numpy as np
import plotly.graph_objects as go
from plotly.graph_objs import Figure

from scipy import spatial
import queue


def _diam_inf(points):
    try:
        candidates = points[spatial.ConvexHull(points).vertices]
    except:
        return np.max(points[:, 0]) - np.min(points[:, 0]) + np.max(points[:, 1]) - np.min(points[:, 1])

    dist_mat = spatial.distance_matrix(candidates, candidates, p=np.inf)  # type: ignore
    return np.max(dist_mat)


def diam_inf(points):
    return np.sqrt(2) * _diam_inf(points)

class Data:
    def __init__(self, root: wit.Hypercube, thickness=0.01):
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

    def _ball_inverse(self, delta, x):
        return wit._pullback(np.array([
            [delta**4, 0, 0],
            [0, delta**2, 0],
            [0, 0, delta**2]
        ]), wit._forward_transformation(x))

    def _ball(self, delta, x):
        return wit._pullback(np.array([
            [1/delta**4, 0, 0],
            [0, 1/delta**2, 0],
            [0, 0, 1/delta**2]
        ]), wit._forward_transformation(x))

    def _CZ_decompose(self, a = 30):
        # root = wit.Hypercube(np.array([0, 0]), 1, points)
        sigma = wit.scale(self._approximate_sigma(), 0.25)

        def indices_of_points(p):
            return [wit.find_index(self.points, q)[0] for q in p]

        def is_good(square: wit.Hypercube):
            p = self.root.search_in(square.dialated(3))
            i = indices_of_points(p)
            diameters = 2 / np.sqrt(np.linalg.eig(sigma[i])[0].min(axis=-1))
            return np.all(diameters >= a * square.width)

        q = queue.Queue()
        q.put(self.root)

        while q.qsize() != 0:
            current = q.get()
            current: wit.Hypercube
            if is_good(current):
                continue
            else:
                current.subdivide()
                for child in current.children:
                    q.put(child)

    def _test_group_sigma(self, s=2):
        sigma = np.array([self._sigma_0(x) for x in self.points])
        groups, well_separated_pairs_indices = wit.build_wspd(self.points, s)

        group_sigma = [
            wit.intersection([
                wit._sum_with_inverse(
                    sigma[j],
                    self._ball_inverse(diam_inf(self.points[groups[i]]), self.points[j])
                )
                for j in groups[i]
            ])
            for i in range(len(groups))
        ]

        return (group_sigma, groups, well_separated_pairs_indices)

    def _test_sigma_bar(self, s=2):
        group_sigma, groups, well_separated_pairs_indices = self._test_group_sigma(s)
        # np.random.shuffle(well_separated_pairs_indices)

        sigma_bar = [[np.empty(0, dtype=float) for _ in range(len(groups))] for _ in range(len(groups))]

        sigma_temp = [
            wit._sum_with_inverse(
                group_sigma[j],
                self._ball_inverse(diam_inf(self.points[groups[j]]), self.points[np.random.choice(groups[j])])
            )
            for j in range(len(groups))
        ]

        for j, k in well_separated_pairs_indices:
            print("j, k: ", j, k)
            # sigma1 = wit._sum_with_inverse(
            #     group_sigma[j],
            #     self._ball_inverse(diam_inf(self.points[groups[j]]), self.points[groups[j][0]])
            # )
            # sigma2 = wit._sum_with_inverse(
            #     group_sigma[k],
            #     self._ball_inverse(diam_inf(self.points[groups[k]]), self.points[groups[k][0]])
            # )
            random_point_in_group_j = self.points[groups[j][0]]
            random_point_in_group_k = self.points[groups[k][0]]

            sigma_bar[j][k] = wit.intersection([
                sigma_temp[j],
                wit._sum_with_inverse(
                    sigma_temp[k],
                    self._ball_inverse(
                        np.linalg.norm(random_point_in_group_j -random_point_in_group_k, ord=np.inf),
                        random_point_in_group_j
                    )
                )
            ])

        return (sigma_bar, groups, well_separated_pairs_indices)


    def _approximate_sigma(self, s=2):
        # Assume all the lambda's are singletons (lambda = [A])

        sigma = np.array([self._sigma_0(x) for x in self.points])
        groups, well_separated_pairs_indices = wit.build_wspd(self.points, s)

        def recursion(sigma):
            # For each A in T, define sigma(A) = intersection(sigma(x) + B(x, diam(A)) for x in A)
            group_sigma = [
                wit.intersection([
                    wit._sum_with_inverse(
                        sigma[j],
                        self._ball_inverse(diam_inf(self.points[groups[i]]), self.points[j])
                    )
                    for j in groups[i]
                ])
                for i in range(len(groups))
            ]

            # For each (A, B) in L,
            #  Define sigma_1(A) = sigma(A) + B(x_A, diam(A))
            #  Define sigma_2(B) = sigma(B) + B(x_B, diam(B))
            # For each (A, B) in L, define sigma(A, B) = intersection(sigma_1(A), sum(sigma_2(B) + B(x_A, |x_A - x_B|)))
            sigma_temp = [
                wit._sum_with_inverse(
                    group_sigma[j],
                    self._ball_inverse(diam_inf(self.points[groups[j]]), self.points[np.random.choice(groups[j])])
                )
                for j in range(len(groups))
            ]

            sigma_bar = [[np.empty(0, dtype=float) for _ in range(len(groups))] for _ in range(len(groups))]

            for j, k in well_separated_pairs_indices:
                sigma_bar[j][k] = wit.intersection([
                    sigma_temp[j],
                    wit._sum_with_inverse(
                        sigma_temp[k],
                        self._ball_inverse(
                            np.linalg.norm(self.points[groups[j][0]] - self.points[groups[k][0]], ord=np.inf),
                            self.points[groups[j][0]]
                        )
                    )
                ])

            # For each A in T, define sigma'(A) = intersection(sigma_bar(A, B) where (A, B) in L)
            sigma_prime = [np.empty(0, dtype=float)] * len(groups)
            for i in range(len(groups)):
                intersectands = []
                for j in range(len(groups)):
                    if [i, j] in well_separated_pairs_indices:
                        intersectands.append(sigma_bar[i][j])

                sigma_prime[i] = wit.intersection(intersectands)

            # For each x in E, redefine sigma(x) = simga(x) intersect intersection(sigma(A) for A in T where x in A)
            new_sigma = [np.empty(0, dtype=float)] * len(self.points)
            for i in range(len(self.points)):
                intersectands = [sigma[i]]
                for j in range(len(groups)):
                    if i in groups[j]:
                        intersectands.append(sigma_prime[j])
                new_sigma[i] = wit.intersection(intersectands)

            return new_sigma

        for _ in range(6):
            sigma = recursion(sigma)

        return np.array([wit._pullback(sigma[i], wit._forward_transformation(-self.points[i])) for i in range(len(self.points))])

parabola = np.array([[x, x**2] for x in np.linspace(-1, 1, 40)])

set1 = parabola[20:] * 0.2
set2 = parabola[20:] * 0.5 + np.array([0.2, 0.2])
set3 = (parabola @ np.array([[0, 1], [-1, 0]])) * 0.035 + np.array([0.5, 0.8])
set4 = (parabola @ np.array([[np.cos(np.pi/6), -np.sin(np.pi/6)], [np.sin(np.pi/6), np.cos(np.pi/6)]])) * 0.25 + np.array([0.75, 0.2])
E = np.concatenate([set1, set2, set3, set4])

# E = np.array([
#     [0.01, 0.9],
#     [0.01, 0.85],
#     [0.01, 0.8],
#     [0.01, 0.75],
#     [0.01, 0.7],
#     [0.01, 0.65],
#     [0.01, 0.6],
#     [0.01, 0.55],
#     [0.01, 0.5],
#     [0.01, 0.45],
#     [0.01, 0.4],
#     [0.01, 0.01],
#     [0.4, 0.01],
#     [0.45, 0.01],
#     [0.5, 0.01],
#     [0.55, 0.01],
#     [0.6, 0.01],
#     [0.65, 0.01],
#     [0.7, 0.01],
#     [0.75, 0.01],
#     [0.8, 0.01],
#     [0.85, 0.01],
#     [0.9, 0.01],
# ])

# E = np.array([[a, 1-a] for a in np.arange(0.1, 1, 0.1)])
# E = np.array([[a, b] for b in np.arange(0.1, 1, 0.1) for a in np.arange(0.1, 1, 0.1)])

# E = np.random.rand(20, 2)

E = np.array([point for point in E if 0 <= point[0] <= 1 and 0 <= point[1] <= 1])
np.random.shuffle(E)

root = wit.Hypercube((0, 0), 1, E)
sigma = Data(root)._approximate_sigma()

# %%
def approximate_polygon(boundary_func, n: int):
    angles = np.linspace(0, 2*np.pi, n)
    directions = np.array([np.cos(angles), np.sin(angles)]).T
    points = boundary_func(directions) * directions.T

    return points.T

def to_boundary_func(elipse: np.ndarray):
    def boundary_func(directions):
        result = []
        for direction in directions:
            result.append(np.sqrt((direction @ direction) / (direction @ elipse @ direction)))
        return np.array(result)

    return boundary_func

def plot_sigma_at(fig: Figure, sigma, x, n=1000):
    # pullback_sigma = wit._pullback(sigma, wit._forward_transformation(-x))
    points = approximate_polygon(to_boundary_func(sigma[1:,1:]), n)
    points = points + x
    fig.add_trace(go.Scatter(x=points[:,0], y=points[:,1], fill="toself"))

def _plot_sigma_at(ax, sigma, x):
    # pullback_sigma = wit._pullback(sigma, wit._forward_transformation(-x))
    n = 1000
    points = approximate_polygon(to_boundary_func(sigma[1:,1:]), n)
    points = points + x
    ax.fill(points[:,0], points[:,1], alpha=0.5)

# import matplotlib.pyplot as plt

# fig, ax = plt.subplots()
# ax.plot(E[:,0], E[:,1], 'x')
# for i in range(len(E)):
#     _plot_sigma_at(ax, sigma[i] * 20, E[i])
# ax.set_xlim(0, 1)
# ax.set_ylim(0, 1)
# ax.set_aspect('equal')

fig = go.Figure()
fig.add_trace(go.Scatter(x=E[:,0], y=E[:,1], mode='markers'))
for i in range(len(E)):
    plot_sigma_at(fig, sigma[i] * 20, E[i])

fig.update_layout(
    xaxis=dict(range=[0, 1]),
    yaxis=dict(range=[0, 1])
)
fig.show()

# %%
a1 = 1
a2 = 1
Data(root)._CZ_decompose()
assign_is_CZ(root, a1, a2)
wit.Plotting.plot_hypercube(root, opacity=1, fillcolor_map=CZ_comparison_fillcolor_map)
# %%
wit.Plotting.plot_hypercube(get_actual_CZ(E, a1, a2), opacity=1)