import numpy as np
from test_module import *
from actual_CZ import *
import timeit

from plotly.subplots import make_subplots
import plotly.graph_objects as go

r = np.arange(60, 1000, 40)
time = []
num_pairs = []

preperation_time = []
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

    # Preperation time
    mark = timeit.default_timer()
    decomposition._pre_arpproximate_sigma()
    preperation_time.append(timeit.default_timer() - mark)

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
fig = make_subplots(rows=7, cols=1, subplot_titles=("Preparation time", "Group sigma", "Sigma temp", "Sigma bar", "Sigma prime", "New sigma", "Total time"))

fig.add_trace(go.Scatter(x=r, y=preperation_time/(r*np.log(r)), mode='lines+markers', name="Preparation time"), row=1, col=1)
fig.add_trace(go.Scatter(x=r, y=group_time/(r*np.log(r)), mode='lines+markers', name="Group sigma"), row=2, col=1)
fig.add_trace(go.Scatter(x=r, y=temp_time/(r*np.log(r)), mode='lines+markers', name="Sigma temp"), row=3, col=1)
fig.add_trace(go.Scatter(x=r, y=bar_time/(r*np.log(r)), mode='lines+markers', name="Sigma bar"), row=4, col=1)
fig.add_trace(go.Scatter(x=r, y=prime_time/(r*np.log(r)), mode='lines+markers', name="Sigma prime"), row=5, col=1)
fig.add_trace(go.Scatter(x=r, y=new_time/(r*np.log(r)), mode='lines+markers', name="New sigma"), row=6, col=1)
fig.add_trace(go.Scatter(x=r, y=time/(r*np.log(r)), mode='lines+markers', name="Total time"), row=7, col=1)

fig.update_layout(height=1000)
fig.show()
# %%
