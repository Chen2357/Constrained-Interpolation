import numpy as np
import numpy.typing as npt
from typing import Union, Dict
import queue
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from collections import defaultdict

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
    def __init__(self, pos: npt.ArrayLike, width: float, points: Union[npt.ArrayLike, None] = None, level: int = 0) -> None:
        """Initialize hypercube.
        
        Inputs:
            pos:    Array, length determines dimension of hypercube.
            width:  Float.
            points: 2D array, axis 1 must match dimesion. Default to None.
            level:  Int Default to 0.
        """
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
        """Dimension of hypercube, determined by length of pos attribute."""
        return len(self.pos)

    @property
    def isLeaf(self):
        """True when number points in this hypercube is <= 1, False otherwise."""
        return len(self.points) <= 1
    
    @property
    def has_no_children(self):
        return all(child is None for child in self.children)

    @property
    def rep(self):
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
        """Populates children attribute with 2^dimension equal sized hypercubes. Partitions points into appropriate children. See [pdf] for details."""
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
            
    def quad_decompose(self):
        """Construct the quadtree recursively using subdivide such that every point is contained within a leaf node."""
        q = queue.Queue()
        q.put(self)
        while q.qsize() != 0:
            cube: Hypercube = q.get()
            if cube.isLeaf: continue
            
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

            if cube.has_no_children: # Make this .isLeaf if looking at quadtree
                leaves.append(cube)
            else:
                for child in cube.children:
                    if child is None: continue
                    q.put(child)
        return leaves

    def plot(self, ax, edgecolor='k', facecolor=None, alpha=0.5):
        if self.dimension != 2:
            raise ValueError("Dimension must be 2 to use plot")

        plt.scatter(self.points[:,0], self.points[:,1], marker = "x")
        pc = PatchCollection(
            [Rectangle(cube.pos, cube.width, cube.width) for cube in self.leaves],
            edgecolor=edgecolor,
            facecolor=facecolor,
            alpha=alpha
        )
        ax.add_collection(pc)

    def well_separated_pairs_decomposition(self, s: float):
        """Uses quadtree wih .self as root to return list of tuples of well seperated hypercubes. See [pdf] for details."""
        return well_separated_pairs(self, self, s)
    
    def contains(self, points: npt.ArrayLike):
        """Returns array of Booleans to indicate which input points lie within .self, returns False otherwise."""
        return np.all(self.pos <= points.reshape(-1,self.dimension), axis = 1) & np.all(points.reshape(-1,self.dimension) < self.pos + self.width,  axis = 1)
    
    def search(self, point: npt.ArrayLike):
        """Returns the leaf node that contains input point."""
        if self.isLeaf and self.contains(point):
            return self
        for child in self.children:
            if child.contains(point):
                return child.search(point)
        raise ValueError(point, self.pos, self.width, "The point is not in this hypercube.")
        

    def is_well_separated(self, other: "Hypercube", s: float):
        """Returns True if .self and input other are s-well-seperated, returns False otherwise."""
        if self.isLeaf and other.isLeaf: return True
        distance = metric_distance(self.rep, other.rep)
        radius = self.width + other.width

        return radius < s * distance
    
    def whitney_square(self, point: npt.ArrayLike, target_width: float):
        """Returns all whitney squares whose 1.1 dilation contains a point"""
        leaf = self.search(point)
        result = []
        level = 0
        width = leaf.width
        while width >= target_width / 4:
            if width > 4 * target_width:
                level += 1 
                width /= 2
                continue
            for i in range(2**level):
                for j in range(2**level):
                    pos = np.array([i * width, j * width]) + leaf.pos
                    center = pos + width / 2 
                    distance = metric_distance(center, point)
                    if distance <= (width / 2) * 1.1:
                        dilated_cube = Hypercube(pos - width, 3 * width)
                        if len(self.search_in(dilated_cube)) <= 1:
                            result.append(Hypercube(pos, width))
            if result != []:
                return result
            level += 1 
            width /= 2
        return result
    
    def whitney_decompose(self):
        """DEPRECATED, DO NOT USE"""
        q = queue.Queue()
        q.put(self)
        while q.qsize() != 0:
            cube: Hypercube = q.get()
            center = cube.pos + cube.width/2
            distances = metric_distance(center, self.points)
            sum_squares = np.sum(distances <= cube.width * 1.5)
            if sum_squares <= 1:
                # print(cube.pos, cube.width)
                continue
        
            cube.subdivide()
            for child in cube.children:
                q.put(child)
                
    def intersects(self, other: "Hypercube"):
        return np.all(self.pos <= other.pos + other.width) and np.all(other.pos <= self.pos + self.width)
    
    def is_subset(self, other: "Hypercube"):
        return np.all(self.pos >= other.pos) and np.all(other.pos + other.width >= self.pos + self.width)
    
    def search_in(self, cube: "Hypercube"):
        q = queue.Queue()
        q.put(self)
        result = []
        while q.qsize() != 0:
            c: Hypercube = q.get() 
            if c.is_subset(cube):
                result.append(c.points)
            elif c.has_no_children:
                result.append(c.points[cube.contains(c.points)])
            else: 
                for child in c.children:
                    if child is not None and child.intersects(cube):
                        q.put(child)
        if result == []:
            return np.empty([0, self.dimension])
        return np.concatenate(result)

        
    def query_nearest_point(self, query: npt.ArrayLike):
        smallest_cube = self
        while not smallest_cube.has_no_children:
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

def well_separated_pairs(u: Hypercube, v: Hypercube, s: float):
    """Returns well-seperated pairs for the Cartesian Product of points in u and v as a list of tuples of hypercubes."""
    if u.isLeaf and v.isLeaf and u == v: return []
    if u.rep is None or v.rep is None:
        return []
    elif u.is_well_separated(v, s):
        return [(u, v)]
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

        return [pair for child in children for pair in well_separated_pairs(child, v, s)]

def all_nearest_neighbors(well_separated_pairs: list[list[Hypercube]], k: int):
    """Returns dictionary containing nearest k neighbors for all points."""
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

        nearest_points = np.concatenate(neighbors[point_bytes])
        nearest_points = np.unique(nearest_points, axis = 0)
        distances = metric_distance(point, nearest_points)

        indices = distances.argpartition(k-1)[:k]
        order = distances[indices].argsort()

        nearest_neighbors[point_bytes] = nearest_points[indices[order]]

    return nearest_neighbors

def metric_distance(point: npt.NDArray, points: npt.NDArray) -> npt.NDArray:
    """Computed distance using L_infinity metric."""
    return np.max(np.abs(points - point).reshape(-1,len(point)), 1)

