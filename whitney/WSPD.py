import wspd
import numpy as np

def build_wspd(points: np.ndarray, s: float):
    num_points = len(points)
    dim = len(points[0])
    wspd_points = [wspd.point(points[i]) for i in range(num_points)]  # type: ignore

    L_half = wspd.build_wspd(num_points, dim, s, wspd_points) # type: ignore
    L_half = [(np.array(A), np.array(B)) for A, B in L_half]

    groups = []
    group_to_index = {}
    well_separated_pairs_indices = []

    for A, B in L_half:
        A_bytes = A.tobytes()
        B_bytes = B.tobytes()

        if A_bytes not in group_to_index:
            group_to_index[A_bytes] = len(groups)
            groups.append(A)
        if B_bytes not in group_to_index:
            group_to_index[B_bytes] = len(groups)
            groups.append(B)

        A_index = group_to_index[A_bytes]
        B_index = group_to_index[B_bytes]
        well_separated_pairs_indices.append((A_index, B_index))
        well_separated_pairs_indices.append((B_index, A_index))

    return (groups, well_separated_pairs_indices)