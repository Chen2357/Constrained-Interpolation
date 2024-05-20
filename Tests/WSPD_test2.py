import wspd
import numpy as np
import plotly.express as px

def estimate_space_complexity(N: int):
    points = np.random.rand(N, 2)

    s = 2
    num_points = len(points)
    dim = len(points[0])
    wspd_points = [wspd.point(points[i]) for i in range(num_points)]  # type: ignore

    L = wspd.build_wspd(num_points, dim, s, wspd_points) # type: ignore

    return sum(len(A) + len(B) for A, B in L)

mean_complexity = []
max_complexity = []

for N in range(500, 10000, 500):
    print(N)
    mean = 0
    maximum = 0
    for _ in range(10):
        mean += estimate_space_complexity(N)
        maximum = max(maximum, estimate_space_complexity(N))
    mean /= 10

    mean_complexity.append((N, mean))
    max_complexity.append((N, maximum))

mean_complexity = np.array(mean_complexity)
max_complexity = np.array(max_complexity)

fig = px.line(mean_complexity, x=0, y=1)
fig.add_scatter(x=max_complexity[:, 0], y=max_complexity[:, 1], mode='lines')
fig.show()