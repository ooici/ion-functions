#!/usr/bin/env python

"""
@package ion_functions.data.adcp_functions
@file ion_functions/data/adcp_functions.py
@author Christopher Wingard
@brief Module containing CTD related data-calculations.
"""

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
        e = "error" velocity profiles in instrument coordinates [mm s-1]
        
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
    # Import Numpy 
    import numpy as np

    theta = 20 / 180 * np.pi
    a = 1 / (2 * np.sin(theta))
    b = 1 / (4 * np.cos(theta))
    c = 1
    d = a / np.sqrt(2)
    
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
    # Import Numpy 
    import numpy as np
    
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
    M1 = np.array([[np.cos(Hrad), np.sin(Hrad), 0],
                  [-np.sin(Hrad), np.cos(Hrad), 0],
                  [0, 0, 1]])
    M2 = np.array([[1, 0, 0], [0, np.cos(P), -np.sin(P)],
                  [0, np.sin(P), np.cos(P)]]);
    M3 = np.array([[np.cos(Rrad), 0, np.sin(Rrad)], [0, 1, 0],
                  [-np.sin(Rrad), 0, np.cos(Rrad)]]);
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

def adcp_magvar(decln, uu, vv):
    """
    Description:

        This function applies the magnetic declination correction to the Earth
        coordinate velocity profiles. The calculation is defined in the Data
        Product Specification for Velocity Profile and Echo Intensity
        - DCN 1341-00750.

    Implemented by:

        2013-04-10: Christopher Wingard. Initial code.

    Usage:

    Example:


    References:
    
        OOI (2012). Data Product Specification for Velocity Profile and Echo
            Intensity. Document Control Number 1341-00750.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00050_Data_Product_SPEC_VELPROF_OOI.pdf)
    """
    # Import Numpy 
    import numpy as np
    
    decln_rad = np.radians(dcln)
    M = np.array([
        [np.cos(decln.rad), np.sin(decln_rad)],
        [-np.sin(decln.rad), np.cos(decln_rad)]
    ])

    nbins = len(uu)
    uu_cor = vv_cor = np.zeros(nbins)
    for i in range(nbins):
        vel = np.array([uu[i], vv[i]]).reshape(2,1)
        cor = np.dot(M, vel)
        uu_cor[i] = cor[0]
        vv_cor[i] = cor[1]
        
    return (uu_cor, vv_cor)
