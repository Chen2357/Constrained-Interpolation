import wspd
import numpy as np

from scipy import spatial

def diam_inf(points):
    if len(points) <= 2:
        return 0
    candidates = points[spatial.ConvexHull(points).vertices]
    dist_mat = spatial.distance_matrix(candidates, candidates, p=np.inf)  # type: ignore
    return np.sqrt(2) * np.max(dist_mat)

points = np.array(((1, 0), (2, 1), (1, 10), (2, 9), (10,4), (11,5), (5,4), (5.5, 3.5)))
s = 2
n = 8

num_points = len(points)
dim = len(points[0])
wspd_points = [wspd.point(points[i]) for i in range(num_points)]  # type: ignore

L = wspd.build_wspd(num_points, dim, s, wspd_points) # type: ignore

print(L)

big_union = []
for A, B in L:
    for x in A:
        for y in B:
            big_union.append((x, y))
            big_union.append((y, x))

print("Starting (1)")

for i in range(n):
    for j in range(n):
        if i != j and (i, j) not in big_union:
            print("Missing point: ", i, j)

print("Done with (1)")

print("Starting (2)")

for A, B in L:
    for C, D in L:
        if A != C or B != D:
            for x in A:
                for y in B:
                    for z in C:
                        for w in D:
                            if x == z and y == w:
                                print("Duplicate point: ", x, y)

print("Done with (2)")

print("Starting (3)")

def diam(A):
    return diam_inf(points[A])

def dist(A, B):
    dist_mat = spatial.distance_matrix(points[A], points[B], p=np.inf) # type: ignore
    return np.min(dist_mat)

for A, B in L:
    if max(diam(A), diam(B)) >= s * dist(A, B):
        print("Bad separation: ", A, B)

print("Done with (3)")