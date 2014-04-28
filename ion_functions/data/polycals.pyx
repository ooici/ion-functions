import numpy as np
cimport numpy as np

np.import_array() # Pretty standard stuff

from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free

cdef extern from "polycals.h":
    ctypedef struct coeff_vector:
        size_t N
        double *coeff
    size_t c_polycal "polycal" (double *out, coeff_vector *cals, double *cal_t, size_t cal_len, double *x, double *t, size_t x_len)


cdef class CoeffVector:
    cdef coeff_vector *data
    cdef int coeff_len
    def __cinit__(self, coeffs):
        # TODO: raise exception on null
        # There's a case where coeffs could be none, an empty array or a whole plethora of other not-nice things

        self.data = <coeff_vector*> PyMem_Malloc(sizeof(coeff_vector) * len(coeffs))
        self.coeff_len = len(coeffs)
        if not self.data:
            raise MemoryError()
        print "memory alloc'd"
        for i, coeff in enumerate(coeffs):
            self.data[i].coeff = <double *> PyMem_Malloc(sizeof(double) * len(coeff))
            if not self.data[i].coeff:
                raise MemoryError()
            print "    internal memory alloc'd"
            self.data[i].N = len(coeff)
            for j,c in enumerate(coeff):
                self.data[i].coeff[j] = <double> c
                print self.data[i].coeff[j]


    def __init__(self, coeffs):
        print "__init__"

    def __dealloc__(self):
        for i in range(self.coeff_len):
            PyMem_Free(self.data[i].coeff)
            print "    internal memory freed"
        PyMem_Free(self.data)
        print "Memory was deallocated"






def polycal(coefficients, calibration_times, x, times):
    raise NotImplementedError("polycal is not ready for use")
    cdef CoeffVector v = CoeffVector(coefficients)
    cdef np.ndarray[double] cal_t = calibration_times
    cdef np.ndarray[double] ix = x
    cdef np.ndarray[double] itimes = times
    cdef np.ndarray[double] out = np.zeros(x.shape, dtype=np.float)

    c_polycal(&out[0], v.data, &cal_t[0], cal_t.shape[0], &ix[0], &itimes[0], x.shape[0])

    return out


