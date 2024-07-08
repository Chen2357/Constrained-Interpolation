import numpy as np
from test_module import *
from actual_CZ import *
import timeit

import wspd

from plotly.subplots import make_subplots
import plotly.graph_objects as go

r = np.arange(60, 1000, 40)
time = []
num_pairs = []

prep_time = []
group_time = []
temp_time = []
bar_time = []
prime_time = []
new_time = []

for i in r:
    E = sample_points(i, "random")
    root = wit.Hypercube((0, 0), 1, E)

    start_time = timeit.default_timer()

    # wit.CZ_decompose(root, N=1, post_shrinking=1)

    decomposition = wit.CZ_Decomposition(root, N=2, post_shrinking=1)

    # Preparation time
    mark = timeit.default_timer()
    decomposition._pre_arpproximate_sigma()
    prep_time.append(timeit.default_timer() - mark)

    # Group sigma
    mark = timeit.default_timer()
    decomposition._compute_group_sigma()
    group_time.append(timeit.default_timer() - mark)

    # Sigma temp
    mark = timeit.default_timer()
    decomposition._compute_sigma_temp()
    temp_time.append(timeit.default_timer() - mark)

    # Sigma bar
    mark = timeit.default_timer()
    decomposition._compute_sigma_bar()
    bar_time.append(timeit.default_timer() - mark)

    # Sigma prime
    mark = timeit.default_timer()
    decomposition._compute_sigma_prime()
    prime_time.append(timeit.default_timer() - mark)

    # New sigma
    mark = timeit.default_timer()
    decomposition._compute_new_sigma()
    new_time.append(timeit.default_timer() - mark)

    time_diff = timeit.default_timer() - start_time
    time.append(time_diff)
    print(i, time_diff)

# %% Plot
def plot(titles, times):
    rows = len(titles)
    fig = make_subplots(rows=rows, cols=1, subplot_titles=titles)

    for i in range(rows):
        fig.add_trace(go.Scatter(x=r, y=times[i]/(r*np.log(r)), mode='lines+markers', name=titles[i]), row=i+1, col=1)

    fig.update_layout(height=150*rows)
    fig.show()

plot([
    "Preparation",
    "Group sigma",
    "Sigma temp",
    "Sigma bar",
    "Sigma prime",
    "New sigma",
    "Total"
], [
    prep_time,
    group_time,
    temp_time,
    bar_time,
    prime_time,
    new_time,
    time
])
# %%
