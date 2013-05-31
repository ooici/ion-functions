import numpy as np
cimport numpy as np

cimport cython

np.import_array()

ctypedef np.uint8_t uint8_t

def cython_confirm():
    print "Confirmed"


def inner_spike(a,b,acc, N, L):
    cdef np.ndarray[double] x
    cdef np.ndarray[double] y
    cdef np.ndarray[signed char] f
    it = np.nditer([a, b, None],
                   flags=['reduce_ok', 'external_loop', 'buffered', 'delay_bufalloc'],
                   op_flags=[['readonly'],['readonly'],['readwrite','allocate']],
                   op_dtypes=['float64','float64','int8'],
                   op_axes=[None, [0, -1], [0, -1]])
    it.operands[-1][...] = 0
    it.reset()
    for ai, bi, oi in it:
        x = ai
        y = bi
        f = oi
        size = x.shape[0]
        for i in range(size):
            value = (N * np.max([x.max() - x.min(), acc])) > np.abs(y[0] - x.mean())
            f[i] = value

    return it.operands[-1]
