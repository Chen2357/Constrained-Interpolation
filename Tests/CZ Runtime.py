import numpy as np
from test_module import *
from actual_CZ import *
import timeit

from plotly.subplots import make_subplots
import plotly.graph_objects as go

r = np.arange(60, 300, 40)
time = []
num_pairs = []

group_time = []
temp_time = []
bar_time = []
prime_time = []
new_time = []

for i in r:
    E = sample_points(i, "circle")
    root = wit.Hypercube((0, 0), 1, E)

    start_time = timeit.default_timer()

    # wit.CZ_decompose(root, N=1, post_shrinking=1)

    decomposition = wit.CZ_Decomposition(root, N=2, post_shrinking=1)

    sigma = np.array([decomposition._sigma_0(x) for x in decomposition.points])
    groups, well_separated_pairs_indices = wit.build_wspd(decomposition.points, decomposition.s)

    decomposition._groups = groups
    decomposition._well_separated_pairs_indices = well_separated_pairs_indices

    # Group sigma
    mark = timeit.default_timer()
    decomposition._group_sigma = decomposition._compute_group_sigma(sigma)
    group_time.append(timeit.default_timer() - mark)

    # Sigma temp
    mark = timeit.default_timer()
    _sigma_temp = decomposition._compute_sigma_temp()
    temp_time.append(timeit.default_timer() - mark)

    # Sigma bar
    mark = timeit.default_timer()
    decomposition._sigma_bar = decomposition._compute_sigma_bar(_sigma_temp)
    bar_time.append(timeit.default_timer() - mark)

    # Sigma prime
    mark = timeit.default_timer()
    decomposition._sigma_prime = decomposition._compute_sigma_prime(decomposition._sigma_bar)
    prime_time.append(timeit.default_timer() - mark)

    # New sigma
    mark = timeit.default_timer()
    decomposition._compute_new_sigma(sigma, decomposition._sigma_prime)
    new_time.append(timeit.default_timer() - mark)

    time_diff = timeit.default_timer() - start_time
    time.append(time_diff)
    print(i, time_diff)

# %% Plot
fig = make_subplots(rows=6, cols=1, subplot_titles=("Group sigma", "Sigma temp", "Sigma bar", "Sigma prime", "New sigma", "Total time"))

fig.add_trace(go.Scatter(x=r, y=group_time/(r*np.log(r)), mode='lines+markers', name="Group sigma"), row=1, col=1)
fig.add_trace(go.Scatter(x=r, y=temp_time/(r*np.log(r)), mode='lines+markers', name="Sigma temp"), row=2, col=1)
fig.add_trace(go.Scatter(x=r, y=bar_time/(r*np.log(r)), mode='lines+markers', name="Sigma bar"), row=3, col=1)
fig.add_trace(go.Scatter(x=r, y=prime_time/(r*np.log(r)), mode='lines+markers', name="Sigma prime"), row=4, col=1)
fig.add_trace(go.Scatter(x=r, y=new_time/(r*np.log(r)), mode='lines+markers', name="New sigma"), row=5, col=1)
fig.add_trace(go.Scatter(x=r, y=time/(r*np.log(r)), mode='lines+markers', name="Total time"), row=6, col=1)

fig.update_layout(height=800)
fig.show()
# %%
