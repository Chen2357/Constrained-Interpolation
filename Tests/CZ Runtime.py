import numpy as np
from test_module import *
from actual_CZ import *
import timeit

import plotly.graph_objects as go
from plotly.subplots import make_subplots


np.random.seed(1)
r = np.arange(10, 40, 50)
time = []
num_pairs = []

# for i in r:
#     E = sample_points(i, "circle")
#     root = wit.Hypercube((0, 0), 1, E)
    
#     start_time = timeit.default_timer()

#     wit.CZ_Decomposition(root, N=2, post_shrinking=1)

#     time_diff = timeit.default_timer() - start_time
#     time.append(time_diff)
#     print(i, time_diff)

group_sigma_counter = []
sigma_temp_counter = []
sigma_bar_counter = []
sigma_prime_counter = []
new_sigma_counter = []

time = []

for i in r:
    E = sample_points(i, "circle")
    root = wit.Hypercube((0, 0), 1, E)
    
    start_time = timeit.default_timer()

    decomposition = wit.CZ_Decomposition(root, N=2, post_shrinking=1)

    time_diff = timeit.default_timer() - start_time
    time.append(time_diff)
    print(i, time_diff)

    group_sigma_counter.append(decomposition.group_sigma_counter)
    sigma_temp_counter.append(decomposition.sigma_temp_counter)
    sigma_bar_counter.append(decomposition.sigma_bar_counter)
    sigma_prime_counter.append(decomposition.sigma_prime_counter)
    new_sigma_counter.append(decomposition.new_sigma_counter)

fig = make_subplots(rows=6, cols=1)

fig.append_trace(go.Line(x=r, y=time), row=1, col=1)
fig.append_trace(go.Line(x=r, y=group_sigma_counter), row=2, col=1)
fig.append_trace(go.Line(x=r, y=sigma_temp_counter), row=3, col=1)
fig.append_trace(go.Line(x=r, y=sigma_bar_counter), row=4, col=1)
fig.append_trace(go.Line(x=r, y=sigma_prime_counter), row=5, col=1)
fig.append_trace(go.Line(x=r, y=new_sigma_counter), row=6, col=1)

fig.update_layout(height=1200, width=600, title_text="Runtime Comparison")
fig.show()