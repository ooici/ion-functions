#!/usr/bin/env python
"""
@package ion_functions.data.vel_functions
@file ion_functions/data/vel_functions.py
@author Stuart Pearce
@brief Module containing velocity family instrument related functions
"""

import numpy as np
import numexpr as ne

from ion_functions.data.generic_functions import magnetic_declination, magnetic_correction


def valid_lat(lat):
    """
    Checks if inputs are valid latitude values.
    """
    if isinstance(lat, np.ndarray):
        if np.any(lat > 90) or np.any(lat < -90):
            return False
        return True
    else:
        return -90 <= lat and lat <= 90


def valid_lon(lon):
    """
    Checks if inputs are valid longitude values.
    """
    if isinstance(lon, np.ndarray):
        if np.any(lon > 180) or np.any(lon < -180):
            return False
        return True
    else:
        return -180 <= lon and lon <= 180


# wrapper functions for use in ION
def nobska_mag_corr_east(u, v, lat, lon, timestamp, z=0):
    """
    Corrects the eastward velocity from a VEL3D
    Nobska MAVS 4 instrument for magnetic declination.

    This function is a wrapper around the method "vel_mag_correction".
    """
    if not valid_lat(lat) or not valid_lon(lon):
        return np.ones(u.shape, dtype=np.float) * -9999

    u_cor = vel_mag_correction(u, v, lat, lon, timestamp, z)[0]
    u_cor = ne.evaluate('u_cor / 100.')  # convert from cm/s to m/s

    return u_cor


def nobska_mag_corr_north(u, v, lat, lon, timestamp, z=0):
    """
    Corrects the northward velocity from a VEL3D
    Nobska MAVS 4 instrument for magnetic declination.

    This function is a wrapper around the method "vel_mag_correction".
    """
    if not valid_lat(lat) or not valid_lon(lon):
        return np.ones(u.shape, dtype=np.float) * -9999

    v_cor = vel_mag_correction(u, v, lat, lon, timestamp, z)[1]
    v_cor = ne.evaluate('v_cor / 100.')  # convert from cm/s to m/s

    return v_cor


def nortek_mag_corr_east(u, v, lat, lon, timestamp, z=0.0):
    """
    Corrects the eastward velocity from VEL3D-CD Nortek Vector
    and VELPT Nortek Aquadopp instruments for magnetic declination.

    This function is a wrapper around the method "vel_mag_correction".
    """
    if not valid_lat(lat) or not valid_lon(lon):
        return np.ones(u.shape, dtype=np.float) * -9999

    u_cor = vel_mag_correction(u, v, lat, lon, timestamp, z)[0]
    u_cor = ne.evaluate('u_cor / 1000.')  # convert from mms/ to m/s

    return u_cor


def nortek_mag_corr_north(u, v, lat, lon, timestamp, z=0.0):
    """
    Corrects the northward velocity from VEL3D-CD Nortek Vector
    and VELPT Nortek Aquadopp instruments for magnetic declination.

    This function is a wrapper around the method "vel_mag_correction".
    """
    if not valid_lat(lat) or not valid_lon(lon):
        return np.ones(u.shape, dtype=np.float) * -9999

    v_cor = vel_mag_correction(u, v, lat, lon, timestamp, z)[1]
    v_cor = ne.evaluate('v_cor/1000.')  # convert from mms/ to m/s

    return v_cor


def vel_mag_correction(u, v, lat, lon, ntp_timestamp, z=0.0, zflag=-1):
    """
    Description:

        Magnetic variation correction for velocities referenced to
        magnetic North. Given the directional velocity components U and
        V, it finds the magnetic declination for the location, depth,
        and date, rotates the vectors from relative to magnetic North to
        true North.

    Implemented by:

        2013-04-17: Stuart Pearce. Initial code.
        2013-04-24: Stuart Pearce. Changed to be general for all velocity
                    instruments.
        2014-02-05: Christopher Wingard. Edited to use magnetic corrections in
                    the generic_functions module.

    Usage:

        u_cor, v_cor = vel_mag_correction(u, v, lat, lon, ntp_timestamp, z, zflag)

            where

        u_cor = eastward velocity profiles, in earth coordinates, with
            the correction for magnetic variation applied.
        v_cor = northward velocity profiles, in earth coordinates,
            with the correction for magnetic variation applied.

        u = uncorrected eastward velocity profiles in Earth coordinates
        v = uncorrected northward velocity profiles in Earth coordinates
        lat = latitude of the instrument [decimal degrees].  East is
            positive, West negative.
        lon = longitude of the instrument [decimal degrees]. North
            is positive, South negative.
        ntp_timestamp = NTP time stamp from a data particle
            [secs since 1900-01-01].
        z = depth or height of instrument relative to sealevel [meters].
            Positive values only. Default value is 0.
        zflag = indicates whether to use z as a depth or height relative
            to sealevel. -1=depth (i.e. -z) and 1=height (i.e. +z). -1
            is the default

    References:

        OOI (2012). Data Product Specification for Turbulent Point Water
            Velocity. Document Control Number 1341-00781.
            https://alfresco.oceanobservatories.org/ (See: Company Home
            >> OOI >> Controlled >> 1000 System Level >>
            1341-00781_Data_Product_SPEC_VELPTTU_Nobska_OOI.pdf)
    """
    # retrieve the magnetic declination
    theta = magnetic_declination(lat, lon, ntp_timestamp, z, zflag)

    # apply the magnetic variation correction
    magvar = np.vectorize(magnetic_correction)
    u_cor, v_cor = magvar(theta, u, v)

    return u_cor, v_cor
