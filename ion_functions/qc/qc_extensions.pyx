import numpy as np
cimport numpy as np

cimport cython

np.import_array()



cdef inline double d_abs(double a) : return a if a > 0 else -a
cdef inline int int_min(int a, int b) : return a if a <= b else b

cdef inline int cmp_res(double x, double y, double reso):
    return (d_abs(x-y) < reso)

cdef extern from "stuck.h":
    int stuck(signed char *out, int out_len, double *dat, int dat_len, double reso, int num)


def stuckvalues(dat, reso, num):
    cdef int dat_shape = dat.shape[0]
    cdef np.ndarray[double] x = dat
    cdef np.ndarray[signed char] out = np.zeros([dat_shape], dtype=np.int8)
    out.fill(1)
    stuck(&out[0], dat_shape, &x[0], dat_shape, reso, num)


    return out
            
