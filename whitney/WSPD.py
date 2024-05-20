import wspd
import numpy as np

def build_wspd(points: np.ndarray, s: float):
    num_points = len(points)
    dim = len(points[0])
    wspd_points = [wspd.point(points[i]) for i in range(num_points)]  # type: ignore

    L_half = wspd.build_wspd(num_points, dim, s, wspd_points) # type: ignore
    L_half: list[tuple[list[int], list[int]]]

    L = L_half + [(B, A) for A, B in L_half]

    T_with_duplicates = [A for tup in L for A in tup]

    T = []
    for A in T_with_duplicates:
        if A not in T:
            T.append(A)

    new_L = [[T.index(A) for A in tup] for tup in L]

    return (T, new_L)