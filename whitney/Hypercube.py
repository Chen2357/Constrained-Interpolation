import numpy as np
import numpy.typing as npt
from typing import Union, Dict, Callable
import queue

class Hypercube:
    """This is a hypercube object that serves as a node in our quadtree.

    Attributes:
        pos:        Position of the corner with smallest value in every direction.
        width:      Side length of hypercube.
        points:     Numpy array containing points in this Hypercube. Axis 0 used to index points, axis 1 contains their coordinates.
        children:   List containing all children hypercubes. See [pdf] for details.
        parent:     Parent of this cube in quadtree.
        level:      How many subdivides to generate this hypercube from quadtree root (root hypercube has level 0).
    """
    def __init__(self,
            pos: npt.ArrayLike,
            width: float,
            points: Union[npt.ArrayLike, None] = None,
            point_indices: Union[npt.ArrayLike, None] = None,
            level: int = 0,
            pointed_jet: Union[npt.ArrayLike, None] = None
        ):
        """Initialize hypercube.

        Inputs:
            pos:    Array, length determines dimension of hypercube.
            width:  Float.
            points: 2D array, axis 1 must match dimesion. Default to None.
            point_indices: 1D array, indices of points in the quadtree. Default to 0, 1, 2, ..., len(points)-1.
            level:  Int Default to 0.
        """
        self.pos = np.array(pos)
        self.width = width
        self.points = np.reshape(points, (-1, len(self.pos))) if points is not None else np.empty([0, len(self.pos)])
        self.point_indices = np.array(point_indices) if point_indices is not None else np.arange(len(self.points))

        self.children: list[Union[Hypercube, None]] = [None for _ in range(2**self.dimension)]
        self.parent: Union[Hypercube, None] = None
        self.level = level
        self.pointed_jet = np.array(pointed_jet if pointed_jet is not None else [])

        if len(self.pos) != np.size(self.points, 1):
           raise ValueError("Dimension mismatch")

    # def __eq__(self, other: "Hypercube") -> bool:
    #     return np.all(self.pos == other.pos) and self.width == other.width

    @property
    def dimension(self):
        """Dimension of hypercube, determined by length of pos attribute."""
        return len(self.pos)

    @property
    def is_leaf(self):
        return all(child is None for child in self.children)

    @property
    def rep(self) -> npt.NDArray | None:
        """"Representative of hypercube. It is the first element in points attribute and None when points are empty."""
        return self.points[0] if len(self.points) else None

    @property
    def root(self):
        """Furthest ancestor of hypercube in quadtree."""
        cube = self
        while cube.parent is not None:
            cube = cube.parent
        return cube

    def subdivide(self):
        """Populates children attribute with 2^dimension equal sized hypercubes. Partitions points into appropriate children."""

        center = self.pos + self.width/2
        quad_info = self.points > center

        for i in range(2**self.dimension):
            kernel = (i // 2**np.arange(self.dimension)) % 2
            # kernel[j] is j-th binary digit of i, 0th digit is unit digit
            is_point_in_quad = np.all(quad_info == kernel, axis=tuple(range(1, self.dimension)))

            cube = Hypercube(
                self.pos + self.width/2 * kernel,
                self.width/2,
                self.points[is_point_in_quad],
                self.point_indices[is_point_in_quad],
                self.level + 1
            )
            cube.parent = self
            self.children[i] = cube

    def decompose(self, method: str):
        if method == "quad":
            self._decompose(lambda cube: len(cube.points) <= 1)
        elif method == "whitney":
            def whitney_criterion(cube: Hypercube) -> np.bool_:
                center = cube.pos + cube.width / 2
                distances = metric_distance(center, self.points)
                num_points_in_dialation = np.sum(distances <= cube.width * 1.5)
                return num_points_in_dialation <= 1

            self._decompose(whitney_criterion)

    def _decompose(self, terminate: Callable[["Hypercube"], Union[bool, np.bool_]]):
        """Construct the quadtree recursively using subdivide until terminate function returns True."""
        q = queue.Queue()
        q.put(self)
        while q.qsize() != 0:
            cube: Hypercube = q.get()
            if terminate(cube): continue

            cube.subdivide()
            for child in cube.children:
                q.put(child)

    def compress(self):
        """Eliminates all internal nodes which contain only one nonempty child. Modifies parent and children attributes accordingly."""
        q = queue.Queue()
        q.put(self)

        while q.qsize() != 0:
            cube: Hypercube
            cube = q.get()

            for i in range(2**self.dimension):
                potential_child = cube.children[i]

                while potential_child is not None:
                    num_points = np.array([len(child.points) if child is not None else 0 for child in potential_child.children])

                    if np.sum(num_points > 0) == 1:
                        j = np.argwhere(num_points > 0)[0][0]
                        potential_child = potential_child.children[j]
                    else:
                        q.put(potential_child)
                        break

                cube.children[i] = potential_child

    @property
    def leaves(self):
        q = queue.Queue()
        q.put(self)

        leaves: list[Hypercube] = []
        while q.qsize() != 0:
            cube: Hypercube
            cube = q.get()

            if cube.is_leaf:
                leaves.append(cube)
            else:
                for child in cube.children:
                    if child is None: continue
                    q.put(child)
        return leaves

    # def well_separated_pairs_decomposition(self, s: float):
    #     """Uses quadtree wih .self as root to return list of tuples of well seperated hypercubes. See [pdf] for details."""
    #     return well_separated_pairs(self, self, s)

    def contains(self, points: npt.ArrayLike) -> npt.NDArray[np.bool_]:
        """Returns array of Booleans to indicate which input points lie within .self, returns False otherwise."""
        points = np.array(points)
        return np.all(self.pos <= points.reshape(-1,self.dimension), axis = 1) & np.all(points.reshape(-1,self.dimension) < self.pos + self.width, axis = 1)

    def search(self, point: npt.ArrayLike):
        """Returns the leaf node that contains query point."""
        if self.is_leaf and self.contains(point):
            return self
        for child in self.children:
            if child is not None and child.contains(point):
                return child.search(point)
        raise ValueError(point, self.pos, self.width, "The point is not in this hypercube.")

    def is_well_separated(self, other: "Hypercube", s: float):
        """Returns True if .self and input other are s-well-seperated, returns False otherwise."""
        if len(self.points) <= 1 and len(other.points) <= 1:
            return True

        assert self.rep is not None
        assert other.rep is not None
        distance = metric_distance(self.rep, other.rep)
        radius = self.width + other.width

        return radius < s * distance

    # def whitney_square(self, point: npt.ArrayLike, target_width: float):
    #     """Returns all whitney squares whose 1.1 dilation contains a point"""
    #     """For test purpose only"""
    #     point = np.array(point)
    #     leaf = self._search_leaf(point)
    #     result = []
    #     level = 0
    #     width = leaf.width
    #     while width >= target_width / 4:
    #         if width > 4 * target_width:
    #             level += 1
    #             width /= 2
    #             continue
    #         for i in [-1, 0, 1]:
    #             for j in [-1, 0, 1]:
    #                 pos = np.array([i * width, j * width]) + leaf.pos
    #                 center = pos + width / 2
    #                 distance = metric_distance(center, point)
    #                 if distance <= (width / 2) * 1.1:
    #                     dilated_cube = Hypercube(pos - width, 3 * width)
    #                     if len(self.search_in(dilated_cube)) <= 1:
    #                         result.append(Hypercube(pos, width))
    #         # if result != []:
    #         #     return result
    #         level += 1
    #         width /= 2
    #       return result

    # def whitney_decompose(self):
    #     """Need to be improved, track points"""
    #     q = queue.Queue()
    #     q.put(self)
    #     while q.qsize() != 0:
    #         cube: Hypercube = q.get()
    #         center = cube.pos + cube.width/2
    #         distances = metric_distance(center, self.points)
    #         sum_squares = np.sum(distances <= cube.width * 1.5)
    #         if sum_squares <= 1:
    #             continue

    #         cube.subdivide()
    #         for child in cube.children:
    #             q.put(child)

    # def cubes_dilation_contains_point(self, point: npt.ArrayLike, dilation_factor: float) -> list["Hypercube"]:
    #     base_cube = self.search(point)

    #     q = queue.Queue()
    #     q.put(self)
    #     candidates = []

    #     while q.qsize() != 0:
    #         cube: Hypercube = q.get()
    #         if base_cube.intersects(Hypercube(cube.pos - cube.width*(dilation_factor-1)/2, cube.width*dilation_factor)):

    #             if cube.is_leaf:
    #                 candidates.append(cube)
    #             else:
    #                 for child in cube.children:
    #                     q.put(child)

    #     result = []
    #     for cube in candidates:
    #         if Hypercube(cube.pos - cube.width*(dilation_factor-1)/2, cube.width*dilation_factor).contains(point):
    #             result.append(cube)

    #     return result

    def intersects(self, other: "Hypercube"):
        "Returns True if two hypercubes intersect, returns False otherwise."
        return np.all(self.pos <= other.pos + other.width) and np.all(other.pos <= self.pos + self.width)

    def is_subset(self, other: "Hypercube"):
        "Returns True if .self is a subset of input other, returns False otherwise."
        return np.all(self.pos >= other.pos) and np.all(other.pos + other.width >= self.pos + self.width)

    def indices_search_in(self, cube: "Hypercube", num_point: int | None = None) -> np.ndarray:
        """Returns all the indices of points (up to a maximum of num_point if not None) in .self that are contained in input cube."""
        q = queue.Queue()
        q.put(self)
        result = []
        while q.qsize() != 0:
            c: Hypercube = q.get()
            if c.is_subset(cube):
                result.append(c.point_indices)
            elif c.is_leaf:
                result.append(c.point_indices[cube.contains(c.points)])
            else:
                for child in c.children:
                    if child is not None and child.intersects(cube):
                        q.put(child)
            if num_point != None and len(result) >= num_point:
                break
        if result == []:
            return np.empty([0, self.dimension])
        return np.concatenate(result)

    def search_in(self, cube: "Hypercube", num_point: int | None = None):
        """Returns all points (up to a maximum of num_point if not None) in .self that are contained in input cube."""
        # q = queue.Queue()
        # q.put(self)
        # result = []
        # while q.qsize() != 0:
        #     c: Hypercube = q.get()
        #     if c.is_subset(cube):
        #         result.append(c.points)
        #     elif c.is_leaf:
        #         result.append(c.points[cube.contains(c.points)])
        #     else:
        #         for child in c.children:
        #             if child is not None and child.intersects(cube):
        #                 q.put(child)
        #     if num_point != None and len(result) >= num_point:
        #         break
        # if result == []:
        #     return np.empty([0, self.dimension])
        # return np.concatenate(result)
        return self.points[self.indices_search_in(cube, num_point)]

    def dialated(self, dialation_factor: float):
        """Returns a new hypercube that is a dilation of .self by the dialation_factor."""
        return Hypercube(self.pos - self.width*(dialation_factor-1)/2, self.width*dialation_factor)

    def query_nearest_point(self, query_point: npt.ArrayLike):
        """Returns the nearest point to query_point in .self."""
        query = np.array(query_point)
        smallest_cube = self
        while not smallest_cube.is_leaf:
            center = smallest_cube.pos + smallest_cube.width/2
            kernel = query > center
            child_index = np.packbits(kernel, bitorder='little')[0]

            if len(smallest_cube.children[child_index].points) >= 1:
                smallest_cube = smallest_cube.children[child_index]
            elif len(smallest_cube.children[child_index].points) == 0:
                break

        candidate_radii = metric_distance(query, smallest_cube.points)
        candidate_index = np.argmin(candidate_radii)
        candidate = smallest_cube.points[candidate_index]
        radius = candidate_radii[candidate_index]

        points = self.search_in(Hypercube(query - radius/2, radius))
        if len(points) == 0:
            return candidate
        nearest_index = np.argmin(metric_distance(query, points))

        return points[nearest_index]

# def well_separated_pairs(u: Hypercube, v: Hypercube, s: float) -> list[tuple[Hypercube, Hypercube]]:
#     """Returns well-seperated pairs for the Cartesian Product of points in u and v as a list of tuples of hypercubes."""
#     if u.isLeaf and v.isLeaf and u == v: return []
#     if u.rep is None or v.rep is None:
#         return []
#     elif u.is_well_separated(v, s):
#         return [(u, v)]
#     else:
#         if not u.isLeaf and not v.isLeaf:
#             swap = u.width < v.width # u.level > v.level
#         else:
#             swap = u.isLeaf

#         if swap:
#             children = v.children
#             v = u
#         else:
#             children = u.children

#         return [pair for child in children if child is not None for pair in well_separated_pairs(child, v, s)]

# def all_nearest_neighbors(well_separated_pairs: list[tuple[Hypercube, Hypercube]], k: int):
#     """Returns dictionary containing nearest k neighbors for all points."""
#     neighbors: Dict["bytes", list[npt.NDArray]] = defaultdict(list[npt.NDArray])
#     for pair in well_separated_pairs:
#         if len(pair[0].points) <= k:
#             for point in pair[0].points:
#                 neighbors[point.tobytes()].append(pair[1].points)
#         if len(pair[1].points) <= k:
#             for point in pair[1].points:
#                 neighbors[point.tobytes()].append(pair[0].points)

#     nearest_neighbors: Dict["bytes", npt.NDArray] = {}
#     for point_bytes in neighbors:
#         point = np.frombuffer(point_bytes)

#         nearest_points = np.concatenate(neighbors[point_bytes])
#         nearest_points = np.unique(nearest_points, axis = 0)
#         distances = metric_distance(point, nearest_points)

#         indices = distances.argpartition(k-1)[:k]
#         order = distances[indices].argsort()

#         nearest_neighbors[point_bytes] = nearest_points[indices[order]]

#     return nearest_neighbors

def metric_distance(point: npt.NDArray, points: npt.NDArray) -> npt.NDArray:
    """Computed distance using L_infinity metric."""
    return np.max(np.abs(points - point).reshape(-1,len(point)), 1)

