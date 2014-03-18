import numpy as np
cimport numpy as np

cimport cython
import os
import datetime
import pkg_resources

np.import_array()

cdef extern from "wmm.h":
    ctypedef struct WMM_Model:
        int initialized
    ctypedef struct velocity_profile:
        size_t len
        double *uu
        double *vv
        double *lat
        double *lon
        double *z
        np.int64_t *timestamp

    int wmm_initialize(char *filename, WMM_Model **model)
    int wmm_free(WMM_Model *model)
    double wmm_declination(WMM_Model *model, double lat, double lon, double z, int year, int month, int day)
    size_t wmm_velocity_correction(velocity_profile *in_vp, WMM_Model *model, velocity_profile *out_vp)

cdef class WMM:
    # CSF change this to us a pointer to a WMM_Model.  wmm_initialize will
    # either perform the same init as before, or, it we are reconfiging to 
    # a previously seen config file, quietly return a pointer to an existing
    # WMM_model

    cdef WMM_Model *model

    def __cinit__(self, filename):
        cdef int retval
        cdef char *fname = filename
        
        # CSF: wmm_initialize already checks for file existance.  Don't check twice,
        # and prefer checking in C
        #if not os.path.exists(filename):
        #    raise OSError("File does not exist")

        retval = wmm_initialize(fname, &self.model)
        if retval:
            raise RuntimeError("Unable to initialized WMM Model")

    property initialized:
        def __get__(self):
            if self.model[0].initialized == 1:
                return True
            return False

    def __dealloc__(self):
        cdef int retval
        retval = wmm_free(self.model)

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


        cdef retval = wmm_declination(<WMM_Model *>self.model, <double>lat, <double>lon, <double>z, <int> date.year, <int> date.month, <int> date.day)
        return retval

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def velocity_correction(self, uu, vv, lat, lon, z, timestamp, zflag=-1):
        '''
        Velocity Correction for magnetic declination
        Returns uu and vv corrected
        Parameters
            uu: Eastern Velocity Vector
            vv: Northern Velocity Vector
            lat: Latitude (degrees north)
            lon: Longitude (degrees east)
            z: Height MSL (km)
            timestamp: UNIX Timestamp
            zflag: 1 = above sea-level, -1 = below sea-level
        '''
        if uu.shape != vv.shape:
            raise TypeError("Vectors are not aligned")
        if np.atleast_1d(lat).shape[0] == 1:
            lat = np.ones(uu.shape, np.float) * lat
        if np.atleast_1d(lon).shape[0] == 1:
            lon = np.ones(uu.shape, np.float) * lon
        z *= zflag
        if np.atleast_1d(z).shape[0] == 1:
            z = np.ones(uu.shape, np.float) * z
        if np.atleast_1d(timestamp).shape[0] == 1:
            timestamp = np.ones(uu.shape, np.int64) * timestamp
        return self._velocity_correction(uu, vv, lat, lon, z, timestamp)


    @cython.boundscheck(False)
    @cython.wraparound(False)
    def _velocity_correction(self, uu, vv, lat, lon, z, timestamp):
        cdef np.ndarray[double] uu_in = uu
        cdef np.ndarray[double] vv_in = vv
        cdef np.ndarray[double] uu_cor = np.empty(uu.shape, np.float)
        cdef np.ndarray[double] vv_cor = np.empty(uu.shape, np.float)
        cdef np.ndarray[double] lat_in = lat
        cdef np.ndarray[double] lon_in = lon
        cdef np.ndarray[double] z_in = z
        cdef np.ndarray[np.int64_t] timestamp_in = timestamp
        cdef size_t retval

        cdef velocity_profile in_vp
        cdef velocity_profile out_vp
        in_vp.len = uu.shape[0]
        in_vp.uu = &uu_in[0]
        in_vp.vv = &vv_in[0]
        in_vp.lat = &lat_in[0]
        in_vp.lon = &lon_in[0]
        in_vp.z = &z_in[0]
        in_vp.timestamp = &timestamp_in[0]

        out_vp.len = uu.shape[0]
        out_vp.uu = &uu_cor[0]
        out_vp.vv = &vv_cor[0]

        retval = wmm_velocity_correction(&in_vp, self.model, &out_vp)
        if retval != uu.shape[0]:
            raise RuntimeError("Failed to Process All Vector Elements")
        return uu_cor, vv_cor

