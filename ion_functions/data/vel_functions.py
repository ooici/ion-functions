#!/usr/bin/env python

"""
@package ion_functions.data.vel_functions
@file ion_functions/data/vel_functions.py
@author Stuart Pearce
@brief Module containing velocity family instrument related functions
"""
# TODO: Get the depth look up values working so that z can be used in
# TODO: wrapper functions.

import numpy as np
import numexpr as ne

from ion_functions.data.adcp_functions import adcp_magvar
from ion_functions.data.generic_functions import magnetic_declination, wmm_model
from ion_functions.data.wmm import WMM


def valid_lat(lat):
    if isinstance(lat, np.ndarray):
        if np.any(lat > 90) or np.any(lat < -90):
            return False
        return True
    else:
        return -90 <= lat and lat <= 90


def valid_lon(lon):
    if isinstance(lon, np.ndarray):
        if np.any(lon > 180) or np.any(lon < -180):
            return False
        return True
    else:
        return -180 <= lon and lon <= 180


# wrapper functions for use in ION
def nobska_mag_corr_east(uu, vv, lat, lon, timestamp, z=0):
    """
    Corrects the eastward velocity from a VEL3D
    Nobska MAVS 4 instrument for magnetic declination.

    This function is a wrapper around the method "velocity_correction"
    of the "ion_functions.data.wmm.WMM" class.
    """
    if not valid_lat(lat) or not valid_lon(lon):
        return np.ones(uu.shape, dtype=np.float) * -9999
    wmm = WMM(wmm_model)
    uu = np.asanyarray(uu, dtype=np.float)
    vv = np.asanyarray(vv, dtype=np.float)
    lat = np.asanyarray(lat, dtype=np.float)
    lon = np.asanyarray(lon, dtype=np.float)
    z = np.asanyarray(z, dtype=np.float)/1000.
    timestamp = np.asanyarray(timestamp, dtype=np.int64) - 2208988800
    uu_cor = wmm.velocity_correction(uu, vv, lat, lon, z, timestamp)[0]
    return ne.evaluate('uu_cor/100.')  # convert from cm/s to m/s


def nobska_mag_corr_north(uu, vv, lat, lon, timestamp, z=0):
    """
    Corrects the northward velocity from a VEL3D
    Nobska MAVS 4 instrument for magnetic declination.

    This function is a wrapper around the method "velocity_correction"
    of the "ion_functions.data.wmm.WMM" class.
    """
    if not valid_lat(lat) or not valid_lon(lon):
        return np.ones(uu.shape, dtype=np.float) * -9999
    wmm = WMM(wmm_model)
    uu = np.asanyarray(uu, dtype=np.float)
    vv = np.asanyarray(vv, dtype=np.float)
    lat = np.asanyarray(lat, dtype=np.float)
    lon = np.asanyarray(lon, dtype=np.float)
    z = np.asanyarray(z, dtype=np.float)/1000.
    timestamp = np.asanyarray(timestamp, dtype=np.int64) - 2208988800
    vv_cor = wmm.velocity_correction(uu, vv, lat, lon, z, timestamp)[1]
    return ne.evaluate('vv_cor/100.')  # convert from cm/s to m/s


def nortek_mag_corr_east(uu, vv, lat, lon, timestamp, z=0):
    """
    Corrects the eastward velocity from VEL3D-CD Nortek Vector
    and VELPT Nortek Aquadopp instruments for magnetic declination.

    This function is a wrapper around the method "velocity_correction"
    of the "ion_functions.data.wmm.WMM" class.
    """
    if not valid_lat(lat) or not valid_lon(lon):
        return np.ones(uu.shape, dtype=np.float) * -9999
    wmm = WMM(wmm_model)
    uu = np.asanyarray(uu, dtype=np.float)
    vv = np.asanyarray(vv, dtype=np.float)
    lat = np.asanyarray(lat, dtype=np.float)
    lon = np.asanyarray(lon, dtype=np.float)
    z = np.asanyarray(z, dtype=np.float)/1000.
    timestamp = np.asanyarray(timestamp, dtype=np.int64) - 2208988800
    uu_cor = wmm.velocity_correction(uu, vv, lat, lon, z, timestamp)[0]
    return ne.evaluate('uu_cor/1000.')  # convert from mms/ to m/s


def nortek_mag_corr_north(uu, vv, lat, lon, timestamp, z=0):
    """
    Corrects the northward velocity from VEL3D-CD Nortek Vector
    and VELPT Nortek Aquadopp instruments for magnetic declination.

    This function is a wrapper around the method "velocity_correction"
    of the "ion_functions.data.wmm.WMM" class.
    """
    if not valid_lat(lat) or not valid_lon(lon):
        return np.ones(uu.shape, dtype=np.float) * -9999
    wmm = WMM(wmm_model)
    uu = np.asanyarray(uu, dtype=np.float)
    vv = np.asanyarray(vv, dtype=np.float)
    lat = np.asanyarray(lat, dtype=np.float)
    lon = np.asanyarray(lon, dtype=np.float)
    z = np.asanyarray(z, dtype=np.float)/1000.
    timestamp = np.asanyarray(timestamp, dtype=np.int64) - 2208988800
    vv_cor = wmm.velocity_correction(uu, vv, lat, lon, z, timestamp)[1]
    return ne.evaluate('vv_cor/1000.')  # convert from mms/ to m/s


### the commented out function below is no longer needed due to a ###
### WMM C code implementation and is now deprecated ####

## proper functions
#def vel_mag_correction(uu, vv, lat, lon, timestamp, z, zflag=-1):
#    """
#    Description:
#
#        Magnetic declination correction for velocities referenced to
#        magnetic North. Given the directional velocity components U and
#        V, it finds the magnetic declination for the location, depth,
#        and date, rotates the vectors from relative to magnetic North to
#        true North. It should be called with the directional string to
#        get an individual horizontal vector component.
#
#    Implemented by:
#
#        2013-04-17: Stuart Pearce. Initial code.
#        2013-04-24: S.P. Changed to be general for all velocities
#
#    Usage:
#
#        vel_corr = vel_mag_correction(uu, vv, lat, lon, timestamp,
#                                      z, dirstr='all')
#
#            where
#
#        vel_corr = The corrected velocity components [D/T].  A tuple of
#            (uu_corr,vv_corr) is returned if dirstr == 'all'.
#        uu = Uncorrected input eastward velocity [D/T]
#        uu = Uncorrected input northward velocity [D/T]
#        lat = Instrument latitude (North positive) [decimal degrees]
#        lon = Instrument longitude (East positive) [decimal degrees]
#        timestamp = The data timestamp (NTP timestamp format)
#            [seconds since 1900-01-01]
#        z = Instrument depth relative to sea level (positive
#            values only) [meters].
#        dirstr = String of the direction for the requested component
#            ('east','north','all').  The default is 'all'.
#
#    References:
#
#        OOI (2012). Data Product Specification for Turbulent Point Water
#            Velocity. Document Control Number 1341-00781.
#            https://alfresco.oceanobservatories.org/ (See: Company Home
#            >> OOI >> Controlled >> 1000 System Level >>
#            1341-00781_Data_Product_SPEC_VELPTTU_Nobska_OOI.pdf)
#    """
#    # retrieve magnetic declination
#    theta = magnetic_declination(lat, lon, timestamp, z, zflag=-1)
#
#    # correct the velocities for magnetic declination
#    #   the algorithm for Nobska & Nortek VELPTTU's are the same as
#    #   adcp_magvar
#    magvar = np.vectorize(adcp_magvar)
#    uu_cor, vv_cor = magvar(theta, uu, vv)
#
#    return (uu_cor, vv_cor)
