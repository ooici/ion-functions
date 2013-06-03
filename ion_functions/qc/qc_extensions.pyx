import numpy as np
cimport numpy as np

cimport cython

np.import_array()

cdef extern from "stuck.h":
    int stuck(signed char *out, double *dat, size_t len, double reso, int num)

cdef extern from "spike.h":
    int spike(signed char *out, double *dat, size_t len, int L, double N, double acc)
    


@cython.boundscheck(False)
@cython.wraparound(False)
def stuckvalues(dat, reso, num):
    cdef int dat_shape = dat.shape[0]
    cdef np.ndarray[double] x = dat
    cdef np.ndarray[signed char] out = np.zeros([dat_shape], dtype=np.int8)
    out.fill(1)
    stuck(&out[0], &x[0], dat_shape, reso, num)

    return out
            
@cython.boundscheck(False)
@cython.wraparound(False)
def spikevalues(dat, L, N, acc):
    cdef int dat_shape = dat.shape[0]
    cdef np.ndarray[double] x = dat
    cdef np.ndarray[signed char] out = np.zeros([dat_shape], dtype=np.int8)
    out.fill(1)
    spike(&out[0], &x[0], dat_shape, L, N, acc)

    return out

