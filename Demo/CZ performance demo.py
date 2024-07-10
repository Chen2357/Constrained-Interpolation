# %% Imports
import sys, os
sys.path.insert(1, os.path.join(sys.path[0], '..'))

import whitney as wit
import numpy as np
import timeit

import plotly.express as px

# %% CZ performance demo
# A quadratic growth of sample
n = np.array([120, 160, 220, 300, 400, 520, 660, 820, 1000, 1200])
time = []

for num in n:
    E = np.random.rand(num, 2)

    start_time = timeit.default_timer()

    root = wit.Hypercube((0, 0), 1, E)
    wit.CZ_decompose(root)

    time_diff = timeit.default_timer() - start_time
    time.append(time_diff)
    print("Number of points:", num, "Time:", time_diff)

# %% Plot
px.line(x=n, y=time, markers=True).show()
# %%
