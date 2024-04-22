import numpy as np


class IsoO_Class:
    
    def __init__(self,data):
        num_g = 1
        m, r, p, T, k, P, s, dP_dr_num, dP_dr_ana = [var[num_g:-(num_g+1)] for var in data[:-1]]
        boundary_i = data[-1]
        num_p = len(m)
        
        self._lo = num_g
        self._hi = num_g + num_p - 1
        self._boundary_i = boundary_i
        self._var = [m, r, p, k, T, P, s]


class FreO_Class:
    
    def __init__(self,data):
        num_g = 1
        m, r, p, T, k, P, s = [var[num_g:-(num_g+1)] for var in data[:-1]]
        boundary_i = data[-1]
        num_p = len(m)
        
        self._lo = num_g
        self._hi = num_g + num_p - 1
        self._boundary_i = boundary_i
        self._var = [m, r, p, k, T, P, s]


class ComO_Class:
    
    def __init__(self,data):
        num_g = 1
        m, r, p, T, k, P, s = [var[num_g:-(num_g+1)] for var in data[:-1]]
        boundary_i = data[-1]
        num_p = len(m)
        
        self._lo = num_g
        self._hi = num_g + num_p - 1
        self._boundary_i = boundary_i
        self._var = [m, r, p, k, T, P, s]

# self._break