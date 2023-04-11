import numpy as np
import numpy.typing as npt
from typing import Union, Dict

Polynomials = Dict["bytes", npt.NDArray]

def evaluate_polynomial(polynomial: npt.NDArray, point):
    degree_plus_1 = len(polynomial)
    
    result = polynomial
    for x in point:
        jet = [x ** n for n in range(degree_plus_1)]
        result = np.tensordot(jet, result, axes = (0,0))
    
    return result