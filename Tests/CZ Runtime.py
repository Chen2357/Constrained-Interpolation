import numpy as np
from test_module import *
from actual_CZ import *
import timeit

import plotly.graph_objects as go

r = np.arange(300, 400, 20)
time = []
num_pairs = []

for i in r:
    E = sample_points(i, "circle")
    root = wit.Hypercube((0, 0), 1, E)
    
    start_time = timeit.default_timer()

    wit.CZ_Decomposition(root, N=2, post_shrinking=1)

    time_diff = timeit.default_timer() - start_time
    time.append(time_diff)
    print(i, time_diff)

fig = go.Figure(data=go.Line(x=r, y=time / (r**2)))
fig.show()
