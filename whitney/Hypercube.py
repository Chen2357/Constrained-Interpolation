import numpy as np
import numpy.typing as npt
from typing import Union, Dict
import queue
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from collections import defaultdict

class Hypercube:
    def __init__(self, pos: npt.ArrayLike, width: float, points: Union[npt.ArrayLike, None] = None, level: int = 0) -> None:
        self.pos = np.array(pos)
        self.width = width
        self.points = np.array(points) if points is not None else np.empty([0, len(self.pos)])

        self.children: list[Union[Hypercube, None]] = [None for _ in range(2**self.dimension)]
        self.parent: Union[Hypercube, None] = None
        self.level = level
        #if len(self.pos) != np.size(self.points, 1):
        #    raise ValueError("Dimension mismatch")

    def __eq__(self, other: "Hypercube"):
        return np.all(self.pos == other.pos) and self.width == other.width

    @property
    def dimension(self):
        return len(self.pos)

    @property
    def isLeaf(self):
        return len(self.points) <= 1

    @property
    def rep(self):
        return self.points[0] if len(self.points) else None
    
    @property
    def root(self):
        cube = self
        while cube.parent is not None:
            cube = cube.parent
        return cube

    def subdivide(self):
        if self.dimension > 8:
            raise RuntimeError("Dimension of Hypercube cannot be greater than 8")

        center = self.pos + self.width/2
        quad_info = self.points > center

        for i in range(2**self.dimension):
            kernel = np.unpackbits(np.uint8(i), count=self.dimension, bitorder='little')
            # kernel[j] is j-th binary digit of i, 0th digit is unit digit
            is_point_in_quad = np.all(quad_info == kernel, axis=1)

            self.children[i] = Hypercube(
                self.pos + self.width/2 * kernel,
                self.width/2,
                self.points[is_point_in_quad],
                self.level + 1
            )
            self.children[i].parent = self
            
    def quadDecompose(self):
        q = queue.Queue()
        q.put(self)
        while q.qsize() != 0:
            cube: Hypercube = q.get()
            if cube.isLeaf: continue
            
            cube.subdivide()
            for child in cube.children:
                q.put(child)

    def compress(self):
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

            if cube.isLeaf:
                leaves.append(cube)
            else:
                for child in cube.children:
                    if child is None: continue
                    q.put(child)
        return leaves

    def plot(self, ax, edgecolor='k', facecolor=None, alpha=0.5):
        if self.dimension != 2:
            raise ValueError("Dimension must be 2 to use plot")

        plt.scatter(self.points[:,0], self.points[:,1], marker="x")
        pc = PatchCollection(
            [Rectangle(cube.pos, cube.width, cube.width) for cube in self.leaves],
            edgecolor=edgecolor,
            facecolor=facecolor,
            alpha=alpha
        )
        ax.add_collection(pc)

    def well_separated_pairs_decomposition(self, s: float):
        return well_separated_paris(self, self, s)
    
    def contains(self, point: npt.ArrayLike):
        return np.all(self.pos < point) and np.all(point < self.pos + self.width)
    
    def search(self, point: npt.ArrayLike):
        for child in self.children:
            if child is None: continue
            if child.contains(point):
                return child.search(point)
        return self

    def is_well_separated(self, other: "Hypercube", s: float):
        if self.isLeaf and other.isLeaf: return True
        distance = np.max(np.abs(self.rep - other.rep))
        radius = self.width + other.width

        return radius < s * distance

def well_separated_paris(u: Hypercube, v: Hypercube, s: float):
    if u.isLeaf and v.isLeaf and u == v: return []
    if u.rep is None or v.rep is None:
        return []
    elif u.is_well_separated(v, s):
        return [[u, v]]
    else:
        if not u.isLeaf and not v.isLeaf:
            swap = u.width < v.width # u.level > v.level
        else:
            swap = u.isLeaf

        if swap:
            children = v.children
            v = u
        else:
            children = u.children

        return [pair for child in children for pair in well_separated_paris(child, v, s)]


def all_nearest_neighbors(well_separated_pairs: list[list[Hypercube]], k: int):
    neighbors: Dict["bytes", list[npt.NDArray]] = defaultdict(list[npt.NDArray])
    for pair in well_separated_pairs:
        if len(pair[0].points) <= k:
            for point in pair[0].points:
                neighbors[point.tobytes()].append(pair[1].points)
        if len(pair[1].points) <= k:
            for point in pair[1].points:
                neighbors[point.tobytes()].append(pair[0].points)

    nearest_neighbors: Dict["bytes", npt.NDArray] = {}
    for point_bytes in neighbors:
        point = np.frombuffer(point_bytes)
        distances = []
        nearest_points = []

        for p in np.concatenate(neighbors[point_bytes]):
            if any(np.all(p == p2) for p2 in nearest_points):
                continue            
            
            distance = np.max(np.abs(p - point))
            index = np.searchsorted(distances, distance)
            if len(nearest_points) < k:
                distances.insert(index, distance)
                nearest_points.insert(index, p)
            elif index < k:
                distances.insert(index, distance)
                nearest_points.insert(index, p)
                distances.pop()
                nearest_points.pop()

        nearest_neighbors[point_bytes] = np.array(nearest_points)

    return nearest_neighbors

def filter_pairs(pairs: list[list[Hypercube]], k):
    result = []
    for pair in pairs:
        if len(pair[0].points) <= k or len(pair[1].points) <= k:
            result.append(pair)
    return result

def find_nearest_neighbor(filtered_pairs: list[list[Hypercube]], point, k):
    neighbors = []
    for pair in filtered_pairs:
        if pair[0].contains(point) and len(pair[1].points) <= k:
            for p in pair[1].points:
                neighbors.append(p)
        elif pair[1].contains(point) and len(pair[0].points) <= k:
            for p in pair[0].points:
                neighbors.append(p)
                
    if len(neighbors) < k:
        raise ValueError("not enough points")

    nearest_points = []
    distances = []
    
    for p in neighbors:
        if any(np.all(p == p2) for p2 in nearest_points):
            continue            
        
        distance = np.max(np.abs(p - point))
        index = np.searchsorted(distances, distance)
        if len(nearest_points) < k:
            distances.insert(index, distance)
            nearest_points.insert(index, p)
        elif index < k:
            distances.insert(index, distance)
            nearest_points.insert(index, p)
            distances.pop()
            nearest_points.pop()
    
    return np.array(nearest_points), np.array(distances)

    
            



