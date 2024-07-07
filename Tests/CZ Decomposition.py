import numpy as np
from test_module import *
from actual_CZ import *

import plotly.graph_objects as go

# %%
# parabola = np.array([[x, x**2] for x in np.linspace(-1, 1, 40)])

# set1 = parabola[20:] * 0.2
# set2 = parabola[20:] * 0.5 + np.array([0.2, 0.2])
# set3 = (parabola @ np.array([[0, 1], [-1, 0]])) * 0.035 + np.array([0.5, 0.8])
# set4 = (parabola @ np.array([[np.cos(np.pi/6), -np.sin(np.pi/6)], [np.sin(np.pi/6), np.cos(np.pi/6)]])) * 0.25 + np.array([0.75, 0.2])
# E = np.concatenate([set1, set2, set3, set4])

E = sample_points(20, "L")

# E = np.array([[a, 1-a] for a in np.arange(0.1, 1, 0.1)])
# E = np.array([[a, b] for b in np.arange(0.1, 1, 0.1) for a in np.arange(0.1, 1, 0.1)])

# E = np.random.rand(20, 2)

E = np.array([point for point in E if 0 <= point[0] <= 1 and 0 <= point[1] <= 1])
np.random.shuffle(E)

root = wit.Hypercube((0, 0), 1, E)
decomposition = wit.CZ_decompose(root, post_shrinking=1)
sigma = decomposition._approximate_sigma()
# %% Plotting

fig = go.Figure()
fig.add_trace(go.Scatter(x=E[:,0], y=E[:,1], mode='markers'))
for i in range(len(E)):
    plot_sigma_at(fig, sigma[i] * 20, E[i])

fig.update_layout(
    xaxis=dict(range=[0, 1]),
    yaxis=dict(range=[0, 1])
)
fig.show()

# %%
a1 = 1
a2 = 0.01
assign_is_CZ(root, a1, a2)
wit.Plotting.plot_hypercube(root, opacity=1, fillcolor_map=CZ_comparison_fillcolor_map)
# %%
actual_CZ = get_actual_CZ(E, a1, a2)
print(np.log2(np.min([leaf.width for leaf in root.leaves])))
print(np.all([leaf.is_CZ for leaf in root.leaves])) # type: ignore
print(np.log2(np.min([leaf.width for leaf in actual_CZ.leaves])))
wit.Plotting.plot_hypercube(actual_CZ, opacity=1)
# %%

for i in np.arange(20, 300, 20):
    np.random.seed(i)
    E = np.random.rand(20, 2)
    root = wit.Hypercube((0, 0), 1, E)
    sigma = wit.CZ_Decomposition(root, N=2)._approximate_sigma()
    assign_is_CZ(root, a1, a2)
    actual_CZ = get_actual_CZ(E, a1, a2)
    print(np.all([leaf.is_CZ for leaf in root.leaves])) # type: ignore
    print(np.log2(np.min([leaf.width for leaf in root.leaves])))
    print(np.log2(np.min([leaf.width for leaf in actual_CZ.leaves])))
# %%
