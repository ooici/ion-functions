#!/usr/bin/env python

"""
@package ion_functions.data.vel_functions
@file ion_functions/data/vel_functions.py
@author Stuart Pearce
@brief Module containing velocity family instrument related functions
"""
# TODO: unittests in test_vel_functions.py, NOTE: no test data for
# Nobska yet

# TODO: Get the depth look up values working so that z can be used in
# TODO: wrapper functions.

import numpy as np

from ion_functions.data.adcp_functions import adcp_magvar
from ion_functions.data.generic_functions import magnetic_declination

# wrapper functions for use in ION
def nobska_mag_corr_east(uu,vv,lat,lon,timestamp,z=0):
    """
    Wrapper function to correct the eastward velocity from a VEL3D
    Nobska instrument for magnetic declination.
    """
    uu_cor = vel_mag_correction(uu,vv,lat,lon,timestamp,z)[0]
    return uu_cor/100.  # convert from cm/s to m/s


def nobska_mag_corr_north(uu,vv,lat,lon,timestamp,z=0):
    """
    Wrapper function to correct the northward velocity from a VEL3D
    Nobska instrument for magnetic declination.  See vel_mag_correction
    function
    """
    vv_cor = vel_mag_correction(uu,vv,lat,lon,timestamp,z)[1]
    return vv_cor/100.  # convert from cms/ to m/s


def nortek_mag_corr_east(uu,vv,lat,lon,timestamp,z=0):
    """
    Wrapper function to correct the eastward velocity from a VEL3D
    Nortek instrument for magnetic declination.  See vel_mag_correction
    function
    """
    uu_cor = vel_mag_correction(uu,vv,lat,lon,timestamp,z)[0]
    return uu_cor/1000.  # convert from mms/ to m/s


def nortek_mag_corr_north(uu,vv,lat,lon,timestamp,z=0):
    """
    Wrapper function to correct the northward velocity from a VEL3D
    Nortek instrument for magnetic declination.  See vel_mag_correction
    function
    """
    vv_cor = vel_mag_correction(uu,vv,lat,lon,timestamp,z)[1]
    return vv_cor/1000.  # convert from mms/ to m/s





# proper functions
def vel_mag_correction(uu, vv, lat, lon, timestamp, z, zflag=-1):
    """
    Description:

        Magnetic declination correction for velocities referenced to
        magnetic North. Given the directional velocity components U and
        V, it finds the magnetic declination for the location, depth,
        and date, rotates the vectors from relative to magnetic North to
        true North. It should be called with the directional string to
        get an individual horizontal vector component.

    Implemented by:

        2013-04-17: Stuart Pearce. Initial code.
        2013-04-24: S.P. Changed to be general for all velocities

    Usage:

        vel_corr = vel_mag_correction(uu, vv, lat, lon, timestamp,
                                      z, dirstr='all')

            where

        vel_corr = The corrected velocity components [D/T].  A tuple of
            (uu_corr,vv_corr) is returned if dirstr == 'all'.
        uu = Uncorrected input eastward velocity [D/T]
        uu = Uncorrected input northward velocity [D/T]
        lat = Instrument latitude (North positive) [decimal degrees]
        lon = Instrument longitude (East positive) [decimal degrees]
        timestamp = The data timestamp (NTP timestamp format)
            [seconds since 1900-01-01]
        z = Instrument depth relative to sea level (positive
            values only) [meters].
        dirstr = String of the direction for the requested component
            ('east','north','all').  The default is 'all'.   
    
    References: 
    
        OOI (2012). Data Product Specification for Turbulent Point Water
            Velocity. Document Control Number 1341-00781.
            https://alfresco.oceanobservatories.org/ (See: Company Home
            >> OOI >> Controlled >> 1000 System Level >>
            1341-00781_Data_Product_SPEC_VELPTTU_Nobska_OOI.pdf)
    """
    # retrieve magnetic declination
    theta = magnetic_declination(lat, lon, timestamp, z, zflag=-1)
    
    # correct the velocities for magnetic declination
    #   the algorithm for Nobska & Nortek VELPTTU's are the same as
    #   adcp_magvar
    magvar = np.vectorize(adcp_magvar)
    uu_cor, vv_cor = magvar(theta, uu, vv)
    
    return (uu_cor, vv_cor)

