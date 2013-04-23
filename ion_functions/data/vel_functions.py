#!/usr/bin/env python

"""
@package ion_functions.data.vel_functions
@file ion_functions/data/vel_functions.py
@author Stuart Pearce
@brief Module containing velocity family instrument related functions
"""

from ion_functions.data.adcp_functions import adcp_magvar
from ion_functions.data.generic_functions import magnetic_declination

def vel3d_b_mag_correction(uu, vv, ww, lat, lon, z, timestamp,
                           zflag=-1, dir_='all'):
    """
    Description:

        A simple wrapper function for the magnetic corrections for the
        VEL3D-B instrument.  Given three directional velocity components
        U, V, and W, it finds the magnetic declination for the
        location and date, provides it to the adcp_magvar function and
        then scales the output from cm/s to m/s.  It should be called
        with the directional to get an individual 3D vector component.


    Implemented by:

        2013-04-17: Stuart Pearce. Initial code.

    Usage:

        vel_comp = vel3d_b_mag_correction(dir_, lat, lon, z, timestamp,
                                  zflag=-1)

            where

        vel_comp = the velocity component in the dir_ direction [m/s]
        uu = input eastward velocity [cm/s]
        uu = input northward velocity [cm/s]
        ww = input upward velocity [cm/s]
        lat = Instrument latitude (North positive) [decimal degrees]
        lon = Instrument longitude (East positive) [decimal degrees]
        z = Instrument depth or height relative to sea level (positive
            values only, use zflag -1 for a depth and zflag +1 for
            height) [meters]
        timestamp = The data timestamp (NTP timestamp format)
            [seconds since 1900-01-01]
        zflag = Flag that sets z to a depth below, or a heigh above sea
            level.  The default is depth (-1 for depth, +1 for height).
        dir_ = String of the direction for the requested component
            ('east','north','up','all').  The default is 'all'.
    
    
    References: 
    
        OOI (2012). Data Product Specification for Turbulent Point Water
            Velocity. Document Control Number 1341-00781.
            https://alfresco.oceanobservatories.org/ (See: Company Home
            >> OOI >> Controlled >> 1000 System Level >>
            1341-00781_Data_Product_SPEC_VELPTTU_Nobska_OOI.pdf)
    """
    # retrieve magnetic declination
    theta = magnetic_declination(lat, lon, z, timestamp, zflag)
    
    # correct the velocities for magnetic declination
    uu_cor, vv_cor = adcp_magvar(theta, uu, vv)
    
    # convert the velocities from cm/s to m/s
    uu_cor /= 100.  # m/s = 100 cm/s
    vv_cor /= 100.
    ww /= 100. 
    
    # return according to dir_ direction
    if dir_ == 'all':
        return (uu_cor, vv_cor, ww)
    elif dir_ == 'east':
        return uu_cor
    elif dir_ == 'north':
        return vv_cor
    elif dir_ == 'up':
        return ww
