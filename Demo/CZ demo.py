# %% Imports
import sys, os
sys.path.insert(1, os.path.join(sys.path[0], '..'))

import whitney as wit
import numpy as np

# %% Compute CZ Decomposition

# Data points
E = np.random.rand(20, 2)

# Post shrinking factor
# Lower value corresponds to finer decomposition
post_shrinking = 0.25

root = wit.Hypercube((0, 0), 1, E)
wit.CZ_decompose(root, post_shrinking=post_shrinking)

wit.Plotting.plot_hypercube(root).show()