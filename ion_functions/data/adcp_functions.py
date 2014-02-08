#!/usr/bin/env python
"""
@package ion_functions.data.adcp_functions
@file ion_functions/data/adcp_functions.py
@author Christopher Wingard
@brief Module containing CTD related data-calculations.
"""

import numpy as np
import numexpr as ne

from ion_functions.data.generic_functions import magnetic_declination, magnetic_correction
from ion_functions.utils import isscalar


# Wrapper functions to create 1:1 outputs for ParameterFunctions in Preload
def adcp_beam_eastward(b1, b2, b3, b4, h, p, r, vf, lat, lon, z, dt):
    """
    Description:

        Wrapper function to compute the Eastward Velocity Profile (VELPROF-VLE)
        from beam coordinate transformed velocity profiles as defined in the
        Data Product Specification for Velocity Profile and Echo Intensity -
        DCN 1341-00750.

    Implemented by:

        2013-04-10: Christopher Wingard. Initial code.
        2014-02-03: Christopher Wingard. Formatting and adjusting to use
                    magnetic declination values calulated use the WMM 2010.

    Usage:

        uu_cor = adcp_beam_eastward(b1, b2, b3, b4, h, p, r, vf, lat, lon, z, dt)

            where

        uu_corr = east velocity profiles in Earth coordinates corrected for the
                  magnetic declination (VELPROF-VLE_L1) [mm s-1]

        b1 = "beam 1" velocity profiles in beam coordinates (VELPROF-B1_L0) [mm s-1]
        b2 = "beam 2" velocity profiles in beam coordinates (VELPROF-B2_L0) [mm s-1]
        b3 = "beam 3" velocity profiles in beam coordinates (VELPROF-B3_L0) [mm s-1]
        b4 = "beam 4" velocity profiles in beam coordinates (VELPROF-B4_L0) [mm s-1]
        h = instrument's uncorrected magnetic heading [degrees]
        p = instrument pitch [degrees]
        r = instrument roll [degrees]
        vf = instrument's vertical orientation (0 = downward looking and
            1 = upward looking)
        lat = instrument's deployment latitude [decimal degrees]
        lon = instrument's deployment longitude [decimal degrees]
        z = instrument's pressure sensor reading (depth) [dm]
        dt = sample date and time value [seconds since 1900-01-01]
    """
    # Test shapes of inputs to determine if we are dealing with one record or
    # several records. Reshape scalars input as size m arrays to size m x 1
    # arrays.
    if isscalar(h) is False:
        h = h.reshape(h.shape[0], 1)
        p = p.reshape(p.shape[0], 1)
        r = r.reshape(r.shape[0], 1)
        vf = vf.reshape(vf.shape[0], 1)

    # compute the beam to instrument transform
    beam2ins = np.vectorize(adcp_beam2ins)
    u, v, w, e = beam2ins(b1, b2, b3, b4)

    # compute the instrument to earth beam transform
    ins2earth = np.vectorize(adcp_ins2earth)
    uu, vv, ww = ins2earth(u, v, w, h, p, r, vf)

    # scale decimeter depth input to meters
    zm = ne.evaluate('z / 10.')

    # calculate the magnetic declination using the WWM2010 model
    theta = magnetic_declination(lat, lon, dt, zm)
    if isscalar(theta) is False:
        theta = theta.reshape(theta.shape[0], 1)

    magvar = np.vectorize(magnetic_correction)
    uu_cor, vv_cor = magvar(theta, uu, vv)

    # scale eastward velocity to m/s
    uu_cor = ne.evaluate('uu_cor / 1000.')  # mm/s -> m/s

    # return the Eastward Velocity Profile
    return uu_cor


def adcp_beam_northward(b1, b2, b3, b4, h, p, r, vf, lat, lon, z, dt):
    """
    Description:

        Wrapper function to compute the Eastward Velocity Profile (VELPROF-VLN)
        from beam coordinate transformed velocity profiles as defined in the
        Data Product Specification for Velocity Profile and Echo Intensity -
        DCN 1341-00750.

    Implemented by:

        2013-04-10: Christopher Wingard. Initial code.
        2014-02-03: Christopher Wingard. Formatting and adjusting to use
                    magnetic declination values calulated use the WMM 2010.

    Usage:

        uu_cor = adcp_beam_eastward(b1, b2, b3, b4, h, p, r, vf, lat, lon, z, dt)

            where

        uu_corr = east velocity profiles in Earth coordinates corrected for the
                  magnetic declination (VELPROF-VLE_L1) [mm s-1]

        b1 = "beam 1" velocity profiles in beam coordinates (VELPROF-B1_L0) [mm s-1]
        b2 = "beam 2" velocity profiles in beam coordinates (VELPROF-B2_L0) [mm s-1]
        b3 = "beam 3" velocity profiles in beam coordinates (VELPROF-B3_L0) [mm s-1]
        b4 = "beam 4" velocity profiles in beam coordinates (VELPROF-B4_L0) [mm s-1]
        h = instrument's uncorrected magnetic heading [degrees]
        p = instrument pitch [degrees]
        r = instrument roll [degrees]
        vf = instrument's vertical orientation (0 = downward looking and
            1 = upward looking)
        lat = instrument's deployment latitude [decimal degrees]
        lon = instrument's deployment longitude [decimal degrees]
        z = instrument's pressure sensor reading (depth) [dm]
        dt = sample date and time value [seconds since 1900-01-01]
    """
    # Test shapes of inputs to determine if we are dealing with one record or
    # several records. Reshape scalars input as size m arrays to size m x 1
    # arrays.
    if isscalar(h) is False:
        h = h.reshape(h.shape[0], 1)
        p = p.reshape(p.shape[0], 1)
        r = r.reshape(r.shape[0], 1)
        vf = vf.reshape(vf.shape[0], 1)

    # compute the beam to instrument transform
    beam2ins = np.vectorize(adcp_beam2ins)
    u, v, w, e = beam2ins(b1, b2, b3, b4)

    # compute the instrument to earth beam transform
    ins2earth = np.vectorize(adcp_ins2earth)
    uu, vv, ww = ins2earth(u, v, w, h, p, r, vf)

    # scale decimeter depth input to meters
    zm = ne.evaluate('z / 10.')

    # calculate the magnetic declination using the WWM2010 model
    theta = magnetic_declination(lat, lon, dt, zm)
    if isscalar(theta) is False:
        theta = theta.reshape(theta.shape[0], 1)

    magvar = np.vectorize(magnetic_correction)
    uu_cor, vv_cor = magvar(theta, uu, vv)

    # scale northward velocity to m/s
    vv_cor = ne.evaluate('vv_cor / 1000.')  # mm/s -> m/s

    # return the Northward Velocity Profile
    return vv_cor


def adcp_beam_vertical(b1, b2, b3, b4, h, p, r, vf):
    """
    Wrapper function to compute the Upward Velocity Profile (VELPROF-VLU)
    from the beam coordinate transformed data.
    """
    # Test shapes of inputs to determine if we are dealing with one record or
    # several records. Reshape scalars input as size m arrays to size m x 1
    # arrays.
    if isscalar(h) is False:
        h = h.reshape(h.shape[0], 1)
        p = p.reshape(p.shape[0], 1)
        r = r.reshape(r.shape[0], 1)
        vf = vf.reshape(vf.shape[0], 1)

    # compute the beam to instrument transform
    beam2ins = np.vectorize(adcp_beam2ins)
    u, v, w, e = beam2ins(b1, b2, b3, b4)

    # compute the instrument to earth beam transform
    ins2earth = np.vectorize(adcp_ins2earth)
    uu, vv, ww = ins2earth(u, v, w, h, p, r, vf)

    # scale vertical velocity to m/s
    ww = ne.evaluate('ww / 1000.')  # mm/s -> m/s

    # return the Upward Velocity Profile
    return ww


def adcp_beam_error(b1, b2, b3, b4):
    """
    Wrapper function to compute the Error Velocity (VELPROF-ERR) from the beam
    coordinate transformed data.
    """
    # compute the beam to instrument transform
    beam2ins = np.vectorize(adcp_beam2ins)
    u, v, w, e = beam2ins(b1, b2, b3, b4)

    # scale error velocity to m/s
    e = ne.evaluate('e/1000.')  # mm/s

    # return the Error Velocity Profile
    return e


def adcp_earth_eastward(u, v, z, lat, lon, dt):
    """
    Wrapper function to compute the Eastward Velocity Profile (VELPROF-VLE)
    from the Earth coordinate transformed data.
    """
    # scale decimeter depth input to meters
    zm = ne.evaluate('z / 10.')

    # calculate the magnetic declination using the WWM2010 model
    theta = magnetic_declination(lat, lon, dt, zm)
    if isscalar(theta) is False:
        theta = theta.reshape(theta.shape[0], 1)

    magvar = np.vectorize(magnetic_correction)
    uu_cor, vv_cor = magvar(theta, u, v)

    # scale eastward velocity from [mm s-1] to [m s-1]
    uu_cor = ne.evaluate('uu_cor / 1000.')

    # return the Eastward Velocity Profile (VELPROF-VLE_L1)
    return uu_cor


def adcp_earth_northward(u, v, z, lat, lon, dt):
    """
    Wrapper function to compute the Northward Velocity Profile (VELPROF-VLN)
    from the Earth coordinate transformed data.
    """
    # scale decimeter depth input to meters
    zm = ne.evaluate('z / 10.')

    # calculate the magnetic declination using the WWM2010 model
    theta = magnetic_declination(lat, lon, dt, zm)
    if isscalar(theta) is False:
        theta = theta.reshape(theta.shape[0], 1)

    magvar = np.vectorize(magnetic_correction)
    uu_cor, vv_cor = magvar(theta, u, v)

    # scale northward velocity from [mm s-1] to [m s-1]
    vv_cor = ne.evaluate('vv_cor / 1000.')

    # return the Northward Velocity Profile (VELPROF-VLN_L1)
    return vv_cor


##### ADCP Beam to Earth Transforms and Magnetic Variation Corrections
def adcp_beam2ins(b1, b2, b3, b4):
    """
    Description:

        This function converts the Beam Coordinate transformed velocity
        profiles to the instrument coordinate system. The calculations are
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

    u = ne.evaluate('c * a * (b1 - b2)')
    v = ne.evaluate('c * a * (b4 - b3)')
    w = ne.evaluate('b * (b1 + b2 + b3 + b4)')
    e = ne.evaluate('d * (b1 + b2 - b3 - b4)')

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

    References:

        OOI (2012). Data Product Specification for Velocity Profile and Echo
            Intensity. Document Control Number 1341-00750.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00050_Data_Product_SPEC_VELPROF_OOI.pdf)
    """
    if vertical == 1:
        R = ne.evaluate('roll + 180.0')
    else:
        R = roll

    Rrad = np.radians(R)
    Hrad = np.radians(heading)
    t1rad = np.radians(pitch)
    t2rad = np.radians(roll)
    P = np.arctan(np.tan(t1rad) * np.cos(t2rad))

    # create the multiplication matrix
    M1 = np.matrix([
        [np.cos(Hrad), np.sin(Hrad), 0.0],
        [-np.sin(Hrad), np.cos(Hrad), 0.0],
        [0.0, 0.0, 1.0]])
    M2 = np.matrix([
        [1.0, 0.0, 0.0],
        [0.0, np.cos(P), -np.sin(P)],
        [0.0, np.sin(P), np.cos(P)]])
    M3 = np.matrix([
        [np.cos(Rrad), 0.0, np.sin(Rrad)],
        [0.0, 1.0, 0.0],
        [-np.sin(Rrad), 0.0, np.cos(Rrad)]])
    M = M1 * M2 * M3

    # apply the Earth transform
    u = np.atleast_1d(u)
    v = np.atleast_1d(v)
    w = np.atleast_1d(w)
    nbins = len(u)
    uu = np.zeros(nbins)
    vv = np.zeros(nbins)
    ww = np.zeros(nbins)
    for i in range(nbins):
        inst = np.array([u[i], v[i], w[i]]).reshape(3, 1)
        vel = np.dot(M, inst).reshape(3)
        uu[i] = np.array(vel)[0][0]
        vv[i] = np.array(vel)[0][1]
        ww[i] = np.array(vel)[0][2]

    return (uu, vv, ww)
