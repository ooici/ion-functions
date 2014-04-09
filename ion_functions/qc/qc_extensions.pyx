import numpy as np
cimport numpy as np

cimport cython

np.import_array()

cdef extern from "stuck.h":
    int stuck(signed char *out, double *dat, size_t len, double reso, int num)

cdef extern from "spike.h":
    int spike(signed char *out, double *dat, size_t len, int L, double N, double acc)

cdef extern from "gradient.h":
    int gradient(signed char *out, double *dat, double *x, size_t len, double grad_min, double grad_max, double mindx, double startdat, double toldat, double skipped_value)
    
cdef extern from "time_utils.h":
    int ntp_month_vector(short int *out, double *input, size_t len)





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

@cython.boundscheck(False)
@cython.wraparound(False)
def gradientvalues(dat, x, grad_min, grad_max, mindx, startdat, toldat, skipped_value=-99.):
    cdef int dat_shape = dat.shape[0]
    cdef np.ndarray[double] idat = dat
    cdef np.ndarray[double] ix = x
    cdef np.ndarray[signed char] out = np.zeros([dat_shape], dtype=np.int8)
    cdef float _skip = skipped_value
    out.fill(1)
    gradient(&out[0], &idat[0], &ix[0], dat_shape, grad_min, grad_max, mindx, startdat, toldat, _skip)
    return out

@cython.boundscheck(False)
@cython.wraparound(False)
def ntp_to_month(dat):
    cdef int dat_shape = dat.shape[0]
    cdef np.ndarray[double] idat = dat
    cdef np.ndarray[short int] out = np.zeros([dat_shape], dtype=np.int16)
    
    ntp_month_vector(&out[0], &idat[0], dat_shape)
    return out

