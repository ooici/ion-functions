#!/usr/bin/env python

"""
@package ion_functions.data.met_functions
@file ion_functions/data/met_functions.py
@author Stuart Pearce
@brief Module containing functions for the met family of instruments
"""
# TODO: Get the height look up values working so that z can be used in
# TODO: the wrapper functions.

import numpy as np
import numexpr as ne

from ion_functions.data.adcp_functions import adcp_magvar
from ion_functions.data.generic_functions import magnetic_declination, wmm_model
from ion_functions.data.wmm import WMM


def windavg_mag_corr_north(uu, vv, lat, lon, timestamp, z=0):
    """
    Corrects the northward wind velocity from a METBK
    instrument for magnetic declination.

    This function is a wrapper around the method "velocity_correction"
    of the "ion_functions.data.wmm.WMM" class.
    """
    wmm = WMM(wmm_model)
    uu = np.asanyarray(uu, dtype=np.float)
    vv = np.asanyarray(vv, dtype=np.float)
    lat = np.asanyarray(lat, dtype=np.float)
    lon = np.asanyarray(lon, dtype=np.float)
    z = np.asanyarray(z, dtype=np.float)/1000.
    timestamp = np.asanyarray(timestamp, dtype=np.int64) - 2208988800
    vv_cor = wmm.velocity_correction(uu, vv, lat, lon, z, timestamp)[1]
    return vv_cor


def windavg_mag_corr_east(uu, vv, lat, lon, timestamp, z=0):
    """
    Corrects the eastward wind velocity from a METBK
    instrument for magnetic declination.

    This function is a wrapper around the method "velocity_correction"
    of the "ion_functions.data.wmm.WMM" class.
    """
    wmm = WMM(wmm_model)
    uu = np.asanyarray(uu, dtype=np.float)
    vv = np.asanyarray(vv, dtype=np.float)
    lat = np.asanyarray(lat, dtype=np.float)
    lon = np.asanyarray(lon, dtype=np.float)
    z = np.asanyarray(z, dtype=np.float)/1000.
    timestamp = np.asanyarray(timestamp, dtype=np.int64) - 2208988800
    uu_cor = wmm.velocity_correction(uu, vv, lat, lon, z, timestamp)[0]
    return uu_cor
