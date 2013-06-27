import numpy as np
cimport numpy as np

cimport cython
import os
import datetime

np.import_array()

cdef extern from "wmm.h":
    ctypedef struct WMM_Model:
        int initialized

    int wmm_initialize(char *filename, WMM_Model *model)
    int wmm_free(WMM_Model *model)
    double wmm_declination(WMM_Model *model, double lat, double lon, double z, int year, int month, int day)


cdef class WMM:
    cdef WMM_Model model

    def __init__(self, filename):
        cdef int retval
        cdef char *fname = filename
        if not os.path.exists(filename):
            raise OSError("File does not exist")

        retval = wmm_initialize(fname, &self.model)
        if retval:
            raise RuntimeError("Unable to initialized WMM Model")

    property initialized:
        def __get__(self):
            if self.model.initialized == 1:
                return True
            return False

    def __del__(self):
        cdef int retval
        retval = wmm_free(&self.model)

        if retval:
            raise RuntimeError("Unable to free WMM Model")

    def declination(self, double lat, double lon, double z, date):
        '''
        lat: degrees north
        lon: degrees east
        z: km (MSL)
        date: datetime.date object
        '''
        if not isinstance(date, datetime.date):
            raise TypeError("date is not a datetime.date object")

        cdef retval = wmm_declination(&self.model, lat, lon, z, <int> date.year, <int> date.month, <int> date.day)
        return retval




