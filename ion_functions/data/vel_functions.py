#!/usr/bin/env python

"""
@package ion_functions.data.vel_functions
@file ion_functions/data/vel_functions.py
@author Stuart Pearce
@brief Module containing velocity family instrument related functions
"""
# TODO: unittests in test_vel_functions.py, NOTE: no test data for
# Nobska yet
import numpy as np

from ion_functions.data.adcp_functions import adcp_magvar
from ion_functions.data.generic_functions import magnetic_declination

# wrapper functions for use in ION
def nobska_mag_corr_east(uu,vv,lat,lon,z,timestamp):
    """
    Wrapper function to correct the eastward velocity from a VEL3D
    Nobska instrument for magnetic declination.
    """
    uu_cor = nobska_mag_correction(uu,vv,lat,lon,z,timestamp,dir_='east')
    return uu_cor
    
def nobska_mag_corr_north(uu,vv,lat,lon,z,timestamp):
    """
    Wrapper function to correct the northward velocity from a VEL3D
    Nobska instrument for magnetic declination.
    """
    vv_cor = nobska_mag_correction(uu,vv,lat,lon,z,timestamp,dir_='north')
    return vv_cor

def nobska_mag_corr_up(ww):
    """
    Converts upward velocity from a VEL3D Nobska instrument from cm/s to
    m/s.
    
    References: 
    
        OOI (2012). Data Product Specification for Turbulent Point Water
            Velocity. Document Control Number 1341-00781.
            https://alfresco.oceanobservatories.org/ (See: Company Home
            >> OOI >> Controlled >> 1000 System Level >>
            1341-00781_Data_Product_SPEC_VELPTTU_Nobska_OOI.pdf)
    """
    if np.size(ww) is 1 and isinstance(ww,float):
        ww = np.array([ww])
    ww_cor = ww/100.
    return ww_cor


# proper functions
def nobska_mag_correction(uu, vv, lat, lon, z, timestamp,
                           zflag=-1, dir_='all'):
    """
    Description:

        Magnetic corrections for the VEL3D-B/Nobska velocities. Given
        the directional velocity components U and V, it finds the
        magnetic declination for the location, depth, and date, rotates
        the vectors from relative to magnetic North to true North, and
        then scales the output from cm/s to m/s. It should be called
        with the directional to get an individual horizontal vector
        component.

    Implemented by:

        2013-04-17: Stuart Pearce. Initial code.

    Usage:

        vel_comp = vel3d_b_mag_correction(uu, vv, lat, lon, z,
                                          timestamp,dir_='all')

            where

        vel_comp = the velocity component in the dir_ direction [m/s]
        uu = input eastward velocity [cm/s]
        uu = input northward velocity [cm/s]
        lat = Instrument latitude (North positive) [decimal degrees]
        lon = Instrument longitude (East positive) [decimal degrees]
        z = Instrument depth relative to sea level (positive
            values only) [meters]
        timestamp = The data timestamp (NTP timestamp format)
            [seconds since 1900-01-01]
        dir_ = String of the direction for the requested component
            ('east','north','all').  The default is 'all'.   
    
    References: 
    
        OOI (2012). Data Product Specification for Turbulent Point Water
            Velocity. Document Control Number 1341-00781.
            https://alfresco.oceanobservatories.org/ (See: Company Home
            >> OOI >> Controlled >> 1000 System Level >>
            1341-00781_Data_Product_SPEC_VELPTTU_Nobska_OOI.pdf)
    """
    # retrieve magnetic declination
    theta = magnetic_declination(lat, lon, z, timestamp, zflag=-1)
    
    # correct the velocities for magnetic declination
    #   the algorithm is the same as adcp_magvar
    uu_cor, vv_cor = adcp_magvar(theta, uu, vv)
    
    # convert the velocities from cm/s to m/s
    uu_cor /= 100.  # m/s = 100 cm/s
    vv_cor /= 100.
    
    # return according to dir_ direction
    if dir_ == 'all':
        return (uu_cor, vv_cor)
    elif dir_ == 'east':
        return uu_cor
    elif dir_ == 'north':
        return vv_cor
