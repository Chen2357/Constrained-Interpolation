import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

import whitney as wit
import numpy as np

point = np.Array([0.2,0.4])
polynomials : wit.Polynomials = {}
polynomials[point.tobytes()] = np.array([[3,2,4],[4,1,0],[1, 0, 0]])

    