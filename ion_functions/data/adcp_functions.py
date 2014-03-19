#!/usr/bin/env python
"""
@package ion_functions.data.adcp_functions
@file ion_functions/data/adcp_functions.py
@author Christopher Wingard
@brief Module containing CTD related data-calculations.
"""
import numpy as np
from ion_functions.data.generic_functions import magnetic_declination, magnetic_correction


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
    # force shapes of inputs to arrays of the correct dimensions
    b1 = np.atleast_2d(b1)
    b2 = np.atleast_2d(b2)
    b3 = np.atleast_2d(b3)
    b4 = np.atleast_2d(b4)
    h = np.atleast_1d(h)
    p = np.atleast_1d(p)
    r = np.atleast_1d(r)
    vf = np.atleast_1d(vf)
    z = np.atleast_1d(z) / 10.  # scale decimeter depth input to meters
    lat = np.atleast_1d(lat)
    lon = np.atleast_1d(lon)
    dt = np.atleast_1d(dt)

    # compute the beam to instrument transform
    u, v, w, e = adcp_beam2ins(b1, b2, b3, b4)

    # calculate the magnetic declination using the WWM2010 model
    theta = magnetic_declination(lat, lon, dt, z)

    # iterate through arrays for processing multiple records
    uu_cor = np.empty((b1.shape))
    for indx in range(b1.shape[0]):
        # compute the instrument to earth beam transform
        uu, vv, ww = adcp_ins2earth(u[indx, :], v[indx, :], w[indx, :],
                                    h[indx], p[indx], r[indx], vf[indx])

        # apply the magnetic variation correction
        uu_cor[indx, :], vcor = magnetic_correction(theta[indx], uu, vv)

    # scale eastward velocity to m/s
    uu_cor = uu_cor / 1000.  # mm/s -> m/s

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
    # force shapes of inputs to arrays of the correct dimensions
    b1 = np.atleast_2d(b1)
    b2 = np.atleast_2d(b2)
    b3 = np.atleast_2d(b3)
    b4 = np.atleast_2d(b4)
    h = np.atleast_1d(h)
    p = np.atleast_1d(p)
    r = np.atleast_1d(r)
    vf = np.atleast_1d(vf)
    z = np.atleast_1d(z) / 10.  # scale decimeter depth input to meters
    lat = np.atleast_1d(lat)
    lon = np.atleast_1d(lon)
    dt = np.atleast_1d(dt)

    # compute the beam to instrument transform
    u, v, w, e = adcp_beam2ins(b1, b2, b3, b4)

    # calculate the magnetic declination using the WWM2010 model
    theta = magnetic_declination(lat, lon, dt, z)

    # iterate through arrays for processing multiple records
    vv_cor = np.empty((b1.shape))
    for indx in range(b1.shape[0]):
        # compute the instrument to earth beam transform
        uu, vv, ww = adcp_ins2earth(u[indx, :], v[indx, :], w[indx, :],
                                    h[indx], p[indx], r[indx], vf[indx])

        # apply the magnetic variation correction
        ucor, vv_cor[indx, :] = magnetic_correction(theta[indx], uu, vv)

    # scale northward velocity to m/s
    vv_cor = vv_cor / 1000.  # mm/s -> m/s

    # return the Northward Velocity Profile
    return vv_cor


def adcp_beam_vertical(b1, b2, b3, b4, h, p, r, vf):
    """
    Wrapper function to compute the Upward Velocity Profile (VELPROF-VLU)
    from the beam coordinate transformed data.
    """
    # force shapes of inputs to arrays of the correct dimensions
    b1 = np.atleast_2d(b1)
    b2 = np.atleast_2d(b2)
    b3 = np.atleast_2d(b3)
    b4 = np.atleast_2d(b4)
    h = np.atleast_1d(h)
    p = np.atleast_1d(p)
    r = np.atleast_1d(r)
    vf = np.atleast_1d(vf)

    # compute the beam to instrument transform
    u, v, w, e = adcp_beam2ins(b1, b2, b3, b4)

    # iterate through arrays for processing multiple records
    ww = np.empty((b1.shape))
    for indx in range(b1.shape[0]):
        # compute the instrument to earth beam transform
        uu, vv, ww[indx, :] = adcp_ins2earth(u[indx, :], v[indx, :], w[indx, :],
                                             h[indx], p[indx], r[indx], vf[indx])

    # scale eastward velocity to m/s
    ww = ww / 1000.  # mm/s -> m/s

    # return the Eastward Velocity Profile
    return ww


def adcp_beam_error(b1, b2, b3, b4):
    """
    Wrapper function to compute the Error Velocity (VELPROF-ERR) from the beam
    coordinate transformed data.
    """
    # force input arrays to 2d shape
    b1 = np.atleast_2d(b1)
    b2 = np.atleast_2d(b2)
    b3 = np.atleast_2d(b3)
    b4 = np.atleast_2d(b4)

    # compute the beam to instrument transform
    u, v, w, e = adcp_beam2ins(b1, b2, b3, b4)

    # scale error velocity to m/s
    e = e / 1000.   # mm/s

    # return the Error Velocity Profile
    return e


def adcp_earth_eastward(u, v, z, lat, lon, dt):
    """
    Wrapper function to compute the Eastward Velocity Profile (VELPROF-VLE)
    from the Earth coordinate transformed data.
    """
    # force shapes of inputs to arrays
    u = np.atleast_2d(u)
    v = np.atleast_2d(v)
    z = np.atleast_1d(z) / 10.  # scale decimeter depth input to meters
    lat = np.atleast_1d(lat)
    lon = np.atleast_1d(lon)
    dt = np.atleast_1d(dt)

    # calculate the magnetic declination using the WWM2010 model
    theta = magnetic_declination(lat, lon, dt, z)

    # iterate through arrays for processing multiple records
    uu_cor = np.empty(u.shape)
    for indx in range(u.shape[0]):
        # apply the magnetic variation correction
        uu_cor[indx, :], vcor = magnetic_correction(theta[indx], u[indx, :], v[indx, :])

    # scale eastward velocity to m/s
    uu_cor = uu_cor / 1000.  # mm/s -> m/s

    # return the Eastward Velocity Profile
    return uu_cor


def adcp_earth_northward(u, v, z, lat, lon, dt):
    """
    Wrapper function to compute the Northward Velocity Profile (VELPROF-VLN)
    from the Earth coordinate transformed data.
    """
    # force shapes of inputs to arrays of the correct dimensions
    u = np.atleast_2d(u)
    v = np.atleast_2d(v)
    z = np.atleast_1d(z) / 10.  # scale decimeter depth input to meters
    lat = np.atleast_1d(lat)
    lon = np.atleast_1d(lon)
    dt = np.atleast_1d(dt)

    # calculate the magnetic declination using the WWM2010 model
    theta = magnetic_declination(lat, lon, dt, z)

    # iterate through arrays for processing multiple records
    vv_cor = np.empty(u.shape)
    for indx in range(u.shape[0]):
        # apply the magnetic variation correction
        ucor, vv_cor[indx, :] = magnetic_correction(theta[indx], u[indx, :], v[indx, :])

    # scale northward velocity to m/s
    vv_cor = vv_cor / 1000.  # mm/s -> m/s

    # return the Northward Velocity Profile
    return vv_cor


##### ADCP Beam to Earth Transforms and Magnetic Variation Corrections
def adcp_backscatter(raw, sfactor):
    """
    Description:

        Converts the echo intensity data from counts to dB using a facotry
        specified scale factor (nominally 0.45 dB/count for the Workhorse
        family of ADCPs and 0.61 dB/count for the ExplorerDVL family). As
        defined in the Data Product Specification for Velocity Profile and Echo
        Intensity - DCN 1341-00750.

    Implemented by:

        2014-04-21: Christopher Wingard. Initial code.

    Usage:

        dB = adcp_backscatter(raw, sfactor)

            where

        dB = Relative Echo Intensity (ECHOINT_L1) [dB]

        raw = raw echo intensity (ECHOINT_L0) [count]
        sfactor = factory supplied scale factor, instrument and beam specific [dB/count]

    References:

        OOI (2012). Data Product Specification for Velocity Profile and Echo
            Intensity. Document Control Number 1341-00750.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00050_Data_Product_SPEC_VELPROF_OOI.pdf)
    """
    if np.isscalar(sfactor) is False:
        sfactor = sfactor.reshape(sfactor.shape[0], 1)

    dB = raw * sfactor
    return dB


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
    # c = 1.0   # not used in subsequent calculations (1 * N = N)
    d = a / np.sqrt(2.0)

    u = a * (b1 - b2)
    v = a * (b4 - b3)
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

    References:

        OOI (2012). Data Product Specification for Velocity Profile and Echo
            Intensity. Document Control Number 1341-00750.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00050_Data_Product_SPEC_VELPROF_OOI.pdf)
    """
    # insure we are dealing with array inputs and determine the number of bins
    # for output arrays
    u = np.atleast_1d(u)
    v = np.atleast_1d(v)
    w = np.atleast_1d(w)

    # if the unit is oriented looking up, add 180 degrees
    mask = (vertical == 1)
    R = roll + (180.0 * mask)

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

    # apply the Earth transform for each instrument velocity bin
    uu, vv, ww = np.array(np.dot(M, np.array([u, v, w])))

    return (uu, vv, ww)
