import numpy as np

def change_base(number, base, memory):
    result = np.zeros(memory, int)
    converted = np.array(list(np.base_repr(number, base)), int)
    result[-converted.shape[0]:] = converted
    return result

def binary_to_int(binary_list, reverse: bool = False):
    result = 0
    for digit in reversed(binary_list) if reverse else binary_list:
        result = (result << 1) | digit
    return result

