import numpy as np
import numpy.typing as npt
from typing import Union
import queue
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

class Hypercube:
    def __init__(self, pos: npt.ArrayLike, width: float, points: Union[npt.ArrayLike, None] = None) -> None:
        self.pos = np.array(pos)
        self.width = width
        self.points = np.array(points) if points is not None else np.empty([0, len(self.pos)])

        self.children: list[Union[Hypercube, None]] = [None for _ in range(2**self.dimension)]
        self.parent: Union[Hypercube, None] = None

        if len(self.pos) != np.size(self.points, 1):
            raise ValueError("Dimension mismatch")

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

    def subdivide(self):
        if self.dimension > 8:
            raise RuntimeError("Dimension of Hypercube cannot be greater than 8")

        center = self.pos + self.width/2
        quad_info = self.points > center

        for i in range(2**self.dimension):
            kernel = np.unpackbits(np.uint8(i), count=self.dimension, bitorder='little')

            is_point_in_quad = np.all(quad_info == kernel, axis=1)

            self.children[i] = Hypercube(
                self.pos + self.width/2 * kernel,
                self.width/2,
                self.points[is_point_in_quad]
            )


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
    
def is_well_separated(u: Hypercube, v: Hypercube, s: float):
    if u.isLeaf and v.isLeaf: return True
    distance = np.max(np.abs(u.rep - v.rep))
    radius = u.width + v.width

    return radius < s * distance

def well_separated_paris(u: Hypercube, v: Hypercube, s: float):
    if u.isLeaf and v.isLeaf and u == v: return []
    if u.rep is None or v.rep is None:
        return []
    elif is_well_separated(u, v, s):
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

def disambiguate_paris(pairs: list[list[Hypercube]]):
    i = 0
    while i < len(pairs):
        for j in range(i+1, len(pairs)):
            if (pairs[i][1] == pairs[j][0] and pairs[i][0] == pairs[j][1]) or pairs[i] == pairs[j]:
                del pairs[j]
                break
        i += 1
    return pairs