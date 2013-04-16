#!/usr/bin/env python

"""
@package ion_functions.data.adcp_functions
@file ion_functions/data/adcp_functions.py
@author Christopher Wingard
@brief Module containing CTD related data-calculations.
"""

from ion_functions.data.generic_functions import magnetic_declination
import numpy as np

# Wrapper functions to create 1:1 outputs for ParameterFunctions in Preload
def adcp_beam_eastward(b1, b2, b3, b4, h, p, r, v, lat, lon, z, dt):
    """
    Wrapper function to compute the Eastward Velocity Profile (VELPROF-VLE)
    from the beam coordinate transformed data.
    """
    # compute the beam to instrument transform
    u, v, w, e = adcp_beam2ins(b1, b2, b3, b4)
    
    # compute the instrument to earth beam transform
    uu, vv, ww = adcp_ins2earth(u, v, w, h, p, r, v)
    
    # calculate the magnetic variation and correct the velocity profiles
    zflag = -1      # sets depth to negative for below sea surface
    theta = magnetic_declination(lat, lon, z, dt, zflag)
    uu_cor, vv_cor = adcp_magvar(theta, uu, vv)
    
    # return the Eastward Velocity Profile
    return uu_cor    


def adcp_beam_northward(b1, b2, b3, b4, h, p, r, v, lat, lon, z, dt):
    """
    Wrapper function to compute the Northward Velocity Profile (VELPROF-VLN)
    from the beam coordinate transformed data.
    """
    # compute the beam to instrument transform
    u, v, w, e = adcp_beam2ins(b1, b2, b3, b4)
    
    # compute the instrument to earth beam transform
    uu, vv, ww = adcp_ins2earth(u, v, w, h, p, r, v)
    
    # calculate the magnetic variation and correct the velocity profiles
    zflag = -1      # sets depth to negative for below sea surface
    theta = magnetic_declination(lat, lon, z, dt, zflag)
    uu_cor, vv_cor = adcp_magvar(theta, uu, vv)
    
    # return the Northward Velocity Profile
    return vv_cor    


def adcp_beam_vertical(b1, b2, b3, b4, h, p, r, v):
    """
    Wrapper function to compute the Upward Velocity Profile (VELPROF-VLU)
    from the beam coordinate transformed data.
    """
    # compute the beam to instrument transform
    u, v, w, e = adcp_beam2ins(b1, b2, b3, b4)
    
    # compute the instrument to earth beam transform
    uu, vv, ww = adcp_ins2earth(u, v, w, h, p, r, v)
    
    # return the Upward Velocity Profile
    return ww    


def adcp_beam_error(b1, b2, b3, b4):
    """
    Wrapper function to compute the Error Velocity (VELPROF-ERR) from the beam
    coordinate transformed data.
    """
    # compute the beam to instrument transform
    u, v, w, e = adcp_beam2ins(b1, b2, b3, b4)    

    # return the Error Velocity Profile
    return e


##### ADCP Beam to Earth Transform and Magnetic Variation Correction Functions 
def adcp_beam2ins(b1, b2, b3, b4):
    """
    Description:

        This function converts the Beam Coordinate transformed velocity
        profiles to the instrument coordinate system. The calculation is
        defined in the Data Product Specification for Velocity Profile and Echo
        Intensity - DCN 1341-00750.

    Implemented by:

        2013-04-10: Christopher Wingard. Initial code.

    Usage:

        u, v, w, e = adcp_beam2ins(b1, b2, b3, b4)

            where

        u = "east" velocity profiles in instrument coordinates [mm s-1]
        v = "north" velocity profiles in instrument coordinates [mm s-1]
        w = "vertical" velocity profiles in instrument coordinates [mm s-1]
        e = "error" velocity profiles [mm s-1]
        
        b1 = "beam 1" velocity profiles in beam coordinates [mm s-1]
        b2 = "beam 2" velocity profiles in beam coordinates [mm s-1]
        b3 = "beam 3" velocity profiles in beam coordinates [mm s-1]
        b4 = "beam 4" velocity profiles in beam coordinates [mm s-1]

    Example:


    References:
    
        OOI (2012). Data Product Specification for Velocity Profile and Echo
            Intensity. Document Control Number 1341-00750.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00050_Data_Product_SPEC_VELPROF_OOI.pdf)
    """

    theta = 20.0 / 180.0 * np.pi
    a = 1.0 / (2.0 * np.sin(theta))
    b = 1.0 / (4.0 * np.cos(theta))
    c = 1.0
    d = a / np.sqrt(2.0)
    
    u = c * a * (b1 - b2)
    v = c * a * (b4 - b3)
    w = b * (b1 + b2 + b3 + b4)
    e = d * (b1 + b2 - b3 - b4)
    
    return (u, v, w, e)

def adcp_ins2earth(u, v, w, heading, pitch, roll, vertical):
    """
    Description:

        This function converts the Instrument Coordinate transformed velocity
        profiles to the Earth coordinate system. The calculation is defined in
        the Data Product Specification for Velocity Profile and Echo Intensity
        - DCN 1341-00750.

    Implemented by:

        2013-04-10: Christopher Wingard. Initial code.

    Usage:

        uu, vu, ww = adcp_ins2earth(u, v, w, heading, pitch, roll, vertical)

            where

        uu = "east" velocity profiles in earth coordinates [mm s-1]
        vv = "north" velocity profiles in earth coordinates [mm s-1]
        ww = "vertical" velocity profiles in earth coordinates [mm s-1]
        
        u = east velocity profiles in instrument coordinates [mm s-1]
        v = north velocity profiles in instrument coordinates [mm s-1]
        w = vertical velocity profiles in instrument coordinates [mm s-1]
        heading = instrument's uncorrected magnetic heading [degrees]
        pitch = instrument pitch [degrees]
        roll = instrument roll [degrees]
        vertical = instrument's vertical orientation (0 = downward looking and
            1 = upward looking)
        
    Example:


    References:
    
        OOI (2012). Data Product Specification for Velocity Profile and Echo
            Intensity. Document Control Number 1341-00750.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00050_Data_Product_SPEC_VELPROF_OOI.pdf)
    """
        
    if vertical == 1:
        R = roll + 180.0
    else:
        R = roll
        
    Rrad = np.radians(R)
    Hrad = np.radians(heading)
    t1rad = np.radians(pitch)
    t2rad = np.radians(roll)
    P = np.arctan(np.tan(t1rad) * np.cos(t2rad))
        
    # create the multiplication matrix
    M1 = np.array([
        [np.cos(Hrad), np.sin(Hrad), 0.0],
        [-np.sin(Hrad), np.cos(Hrad), 0.0],
        [0.0, 0.0, 1.0]
        ])
    M2 = np.array([
        [1.0, 0.0, 0.0],
        [0.0, np.cos(P), -np.sin(P)],
        [0.0, np.sin(P), np.cos(P)]    
        ])
    M3 = np.array([
        [np.cos(Rrad), 0.0, np.sin(Rrad)],
        [0.0, 1.0, 0.0],
        [-np.sin(Rrad), 0.0, np.cos(Rrad)]
        ])
    M = M1 * M2 * M3
    
    # apply the Earth transform
    nbins = len(u)
    uu = vv = ww = np.zeros(nbins)
    for i in range(nbins):
        inst = np.array([u[i], v[i], w[i]]).reshape(3,1)
        vel = np.dot(M,inst).reshape(3)
        uu[i] = vel[0]
        vv[i] = vel[1]
        ww[i] = vel[2]
        
    return (uu, vv, ww)


def adcp_magvar(theta, uu, vv):
    """
    Description:

        This function corrects the velocity profiles for the magnetic
        declination at the measurement location. The calculation is defined in
        the Data Product Specification for Velocity Profile and Echo Intensity
        - DCN 1341-00750. Magnetic declination is obtained from the geomag
        toolbox.

    Implemented by:

        2013-04-10: Christopher Wingard. Initial code.

    Usage:

        uu_cor, vv_cor = adcp_magvar(theta, uu, vv)

            where

        uu_cor = eastward velocity profiles in earth coordinates [mm s-1], with
            the correction for magnetic variation applied.
        vv_cor = northward velocity profiles in earth coordinates [mm s-1],
            with the correction for magnetic variation applied.
        
        theta = magnetic variation based on location (latitude, longitude and
            altitude) and date [degrees]
        uu = uncorrected eastward velocity profiles in Earth coordinates
            [mm s-1]
        vv = uncorrected northward velocity profiles in Earth coordinates
            [mm s-1]
    
    Example:


    References:
    
        OOI (2012). Data Product Specification for Velocity Profile and Echo
            Intensity. Document Control Number 1341-00750.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00050_Data_Product_SPEC_VELPROF_OOI.pdf)
    """
    
    theta_rad = np.radians(theta)
    M = np.array([
        [np.cos(theta_rad), np.sin(theta_rad)],
        [-np.sin(theta_rad), np.cos(theta_rad)]
    ])

    nbins = len(uu)
    uu_cor = vv_cor = np.zeros(nbins)
    for i in range(nbins):
        vel = np.array([uu[i], vv[i]]).reshape(2,1)
        cor = np.dot(M, vel).reshape(2)
        uu_cor[i] = cor[0]
        vv_cor[i] = cor[1]
        
    return (uu_cor, vv_cor)

