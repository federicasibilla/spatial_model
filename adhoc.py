import numpy as np

def R0_adhoc():
    R0 =  np.zeros((101, 101, 2))
    R0[:101//2, :, 0] = 100
    R0[101//2:, :, 1] = 100
    return R0

def N0_adhoc():
    N0 =  np.zeros((101, 101, 2))
    N0[:101//2, :, 0] = 100
    N0[101//2:, :, 1] = 100
    return N0

