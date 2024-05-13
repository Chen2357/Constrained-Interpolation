import wspd
import numpy as np

def build_wspd(points: np.ndarray, s: float):
    num_points = len(points)
    dim = len(points[0])
    wspd_points = [wspd.point(points[i]) for i in range(num_points)]  # type: ignore

    L = wspd.build_wspd(num_points, dim, s, wspd_points) # type: ignore
    L: list[tuple[list[int], list[int]]]

    T_with_duplicates = [set(A) for tup in L for A in tup]

    T = []
    for A in T_with_duplicates:
        if A not in T:
            T.append(A)

    return (T, L)