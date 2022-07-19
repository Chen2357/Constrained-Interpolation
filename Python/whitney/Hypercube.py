import numpy as np
import numpy.typing as npt
from typing import Union
import queue
import matplotlib.pyplot as plt
from . import Utilities as ut
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

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
    
    def contain(self, point: npt.ArrayLike):
        return np.all(self.pos < point) and np.all(point < self.pos + self.width)
    
    def search(self, point: npt.ArrayLike):
        for child in self.children:
            if child is None: continue
            if child.contain(point):
                return child.search(point)
        return self

    def pos_in_level(self):
        pos = (self.pos - self.root.pos) / self.root.width
        width = self.width / self.root.width

        pos /= width
        pos = np.round(pos).astype(int)
        return pos
    
    def follow_path(self, path, test_point):
        cube = self
        level = cube.level
        i = 0
        while i < len(path):
            if cube.isLeaf:
                return cube
            if len(cube.children[path[i]].points) == 0:
                return cube
            if not cube.children[path[i]].contain(test_point):
                return cube.children[path[i]]
            cube = cube.children[path[i]]
            i += cube.level - level
            level = cube.level
        return cube                       
    
    def find_nearest_neighbor(self, point):
        cutoff = 10
        points_to_check = np.empty([0,self.dimension])
        q = queue.Queue()
        cube = self.search(point)
        pos = cube.pos_in_level()
        center = cube.pos + cube.width/2
        for i in range(0,3**self.dimension):
            relative_pos = ut.change_base(i,3,self.dimension)-1
            if np.any(0 > pos + relative_pos) or np.any(pos + relative_pos >= 2 ** cube.level): 
                continue
            if np.all(relative_pos == 0):
                continue
            path = path_from_pos_in_level(pos + relative_pos, cube.level)
            candiate_cube = self.follow_path(path, center + cube.width * relative_pos)
            q.put(candiate_cube)
        
        while q.qsize() != 0:
            cube = q.get()
            if cube.contain(point):
                continue
            if cube.isLeaf:
                points_to_check = np.append(points_to_check, cube.points, 0)
                continue
            if len(cube.points) > cutoff: # FLESH OUT
                points_to_check = np.append(points_to_check, cube.points, 0)
                continue
            else:
                points_to_check = np.append(points_to_check, cube.points, 0)
        
        nearest_index = np.argmin(np.max(np.abs(points_to_check - point),1))
        return points_to_check[nearest_index]        
            
def path_from_pos_in_level(pos, level):
    pos_binary = np.zeros((len(pos), level), int)
    for i in range(len(pos)):
        pos_binary[i]= ut.change_base(pos[i], 2, level)
    path = np.zeros(level, int)
    for j in range(level):
        path[j] = ut.binary_to_int(pos_binary[:,j],True)
    return path                 
    
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




