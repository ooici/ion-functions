#!/usr/bin/env python
"""
@package ion_functions.data.adcp_functions
@file ion_functions/data/adcp_functions.py
@author Christopher Wingard
@brief Module containing ADCP related data-calculations.
"""
import numpy as np
from ion_functions.data.generic_functions import magnetic_declination


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
                    magnetic declination values calculated use the WMM 2010.
        2014-04-04: Russell Desiderio. Optimized code performance by replacing the
                    for loops previously used to calculate 2D and 3D vectorized
                    coordinate transformations with calls to np.einsum (numpy
                    Einstein summation function).

    Usage:

        uu_cor = adcp_beam_eastward(b1, b2, b3, b4, h, p, r, vf, lat, lon, z, dt)

            where

        uu_corr = east velocity profiles in Earth coordinates corrected for the
                  magnetic declination (VELPROF-VLE_L1) [m s-1]

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

    # compute the instrument to earth beam transform
    uu, vv, _ = adcp_ins2earth(u, v, w, h, p, r, vf)

    # compute the magnetic variation, and
    theta = magnetic_declination(lat, lon, dt, z)

    # correct for it
    uu_cor, _ = magnetic_correction(theta, uu, vv)

    # scale velocity to m/s
    uu_cor = uu_cor / 1000.  # mm/s -> m/s

    # return the Eastward Velocity Profile
    return uu_cor


def adcp_beam_northward(b1, b2, b3, b4, h, p, r, vf, lat, lon, z, dt):
    """
    Description:

        Wrapper function to compute the Northward Velocity Profile (VELPROF-VLN)
        from beam coordinate transformed velocity profiles as defined in the
        Data Product Specification for Velocity Profile and Echo Intensity -
        DCN 1341-00750.

    Implemented by:

        2013-04-10: Christopher Wingard. Initial code.
        2014-02-03: Christopher Wingard. Formatting and adjusting to use
                    magnetic declination values calculated use the WMM 2010.
        2014-03-28: Russell Desiderio. Corrected documentation only (east->north).
        2014-04-04: Russell Desiderio. Optimized code performance by replacing the
                    for loops previously used to calculate 2D and 3D vectorized
                    coordinate transformations with calls to np.einsum (numpy
                    Einstein summation function).

    Usage:

        vv_cor = adcp_beam_northward(b1, b2, b3, b4, h, p, r, vf, lat, lon, z, dt)

            where

        vv_corr = north velocity profiles in Earth coordinates corrected for the
                  magnetic declination (VELPROF-VLN_L1) [m s-1]

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

    # compute the instrument to earth beam transform
    uu, vv, _ = adcp_ins2earth(u, v, w, h, p, r, vf)

    # compute the magnetic variation, and
    theta = magnetic_declination(lat, lon, dt, z)

    # corect for it
    _, vv_cor = magnetic_correction(theta, uu, vv)

    # scale velocity to m/s
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

    # compute the instrument to earth beam transform
    _, _, ww = adcp_ins2earth(u, v, w, h, p, r, vf)

    # scale upward velocity to m/s
    ww = ww / 1000.  # mm/s -> m/s

    # return the Upward Velocity Profile
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

    This function was written to test intermediate data product values, and
    is not used to calculate the final data products. However, it is used
    in the unit tests, and so should not be deleted.
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

    # correct u and v for magnetic declination
    uu_cor, _ = magnetic_correction(theta, u, v)

    # scale velocity to m/s
    uu_cor = uu_cor / 1000.  # mm/s -> m/s

    # return the Eastward Velocity Profile
    return uu_cor


def adcp_earth_northward(u, v, z, lat, lon, dt):
    """
    Wrapper function to compute the Northward Velocity Profile (VELPROF-VLN)
    from the Earth coordinate transformed data.

    This function was written to test intermediate data product values, and
    is not used to calculate the final data products. However, it is used
    in the unit tests, and so should not be deleted.
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

    # correct u and v for magnetic declination
    _, vv_cor = magnetic_correction(theta, u, v)

    # scale velocity to m/s
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
    c = 1.0   # +1.0 for convex transducer head, -1 for concave
    d = a / np.sqrt(2.0)

    u = c * a * (b1 - b2)
    v = c * a * (b4 - b3)
    w = b * (b1 + b2 + b3 + b4)
    e = d * (b1 + b2 - b3 - b4)

    return (u, v, w, e)


def vadcp_beam2ins(b1, b2, b3, b4, b5):
    """
    Description:

        This function converts the 5 Beam Coordinate transformed velocity
        profiles for the VADCP to the instrument coordinate system. The
        calculations are defined in the Data Product Specification for
        Turbulent Velocity Profile and Echo Intensity - DCN 1341-00760.

    Implemented by:

        2014-06-25: Christopher Wingard. Initial code.

    Usage:

        u, v, w, e = vadcp_beam2ins(b1, b2, b3, b4, b5)

            where

        u = "east" velocity profiles in instrument coordinates [mm s-1]
        v = "north" velocity profiles in instrument coordinates [mm s-1]
        w = "vertical" velocity profiles in instrument coordinates [mm s-1]
        e = "error" velocity profiles [mm s-1]

        b1 = "beam 1" velocity profiles in beam coordinates [mm s-1]
        b2 = "beam 2" velocity profiles in beam coordinates [mm s-1]
        b3 = "beam 3" velocity profiles in beam coordinates [mm s-1]
        b4 = "beam 4" velocity profiles in beam coordinates [mm s-1]
        b5 = "beam 5" velocity profiles in beam coordinates [mm s-1]

    References:

        OOI (2012). Data Product Specification for Velocity Profile and Echo
            Intensity. Document Control Number 1341-00750.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00050_Data_Product_SPEC_VELPROF_OOI.pdf)
    """
    # Setup the transformation matrix based on the following constants
    theta = 20.0 / 180.0 * np.pi        # fixed at 20 degrees for 4 beam unit
    a = 1.0 / (2.0 * np.sin(theta))
    b = 1.0 / (4.0 * np.cos(theta))
    c = 1.0   # +1.0 for convex transducer head, -1 for concave
    d = a / np.sqrt(2.0)
    # for the 5th beam, theta equals 0. Thus, a = 0, d = 0 and b = 0.25

    # The transformation matrix, is applied as follows
    # u = | c * a -c * a     0      0     0 | * beam 1-5
    # v = |   0      0    -c * a  c * a   0 | * beam 1-5
    # w = |   b      b      -b     -b     b | * beam 1-5
    # e = |   d      d      -d     -d     d | * beam 1-5

    # or, in other words
    u = c * a * (b1 - b2)
    v = c * a * (-b3 + b4)
    w = b * (b1 + b2 + b3 + b4) + 0.25 * b5
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
        2014-04-04: Russell Desiderio. Optimized code performance by replacing the for
                    loops previously used to calculate vectorized matrix multiplication
                    products with calls to np.einsum (numpy Einstein summation function).

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
    # insure we are dealing with array inputs
    u = np.atleast_2d(u)
    v = np.atleast_2d(v)
    w = np.atleast_2d(w)

    heading = np.atleast_1d(heading)
    pitch = np.atleast_1d(pitch)
    roll = np.atleast_1d(roll)
    vertical = np.atleast_1d(vertical)

    # if the unit is oriented looking up, add 180 degrees
    mask = (vertical == 1)
    R = roll + (180.0 * mask)

    # roll
    Rrad = np.radians(R)
    cos_R = np.cos(Rrad)
    sin_R = np.sin(Rrad)
    # heading
    Hrad = np.radians(heading)
    cos_H = np.cos(Hrad)
    sin_H = np.sin(Hrad)
    # pitch
    t1rad = np.radians(pitch)
    t2rad = np.radians(roll)
    Prad = np.arctan(np.tan(t1rad) * np.cos(t2rad))
    cos_P = np.cos(Prad)
    sin_P = np.sin(Prad)

    # determine array size
    n_packets = u.shape[0]
    n_uvw = u.shape[1]

    # initialize vectors to be used as matrix elements
    ones = np.ones(n_packets)
    zeros = ones * 0.0

    # the rollaxis calls reorient the matrices so that their lead index is
    # the data packet index
    M1 = np.array([[cos_H, sin_H, zeros],
                   [-sin_H, cos_H, zeros],
                   [zeros, zeros, ones]])
    M1 = np.rollaxis(M1, 2)
    M2 = np.array([[ones, zeros, zeros],
                   [zeros, cos_P, -sin_P],
                   [zeros, sin_P, cos_P]])
    M2 = np.rollaxis(M2, 2)
    M3 = np.array([[cos_R, zeros, sin_R],
                   [zeros, ones, zeros],
                   [-sin_R, zeros, cos_R]])
    M3 = np.rollaxis(M3, 2)

    # construct input array of coordinates (velocities) to be transformed.
    # the basis set is 3D (E,N,U) so that the middle dimension is sized at 3.
    uvw = np.zeros((n_packets, 3, n_uvw))

    # pack the coordinates (velocities) to be transformed into the appropriate slices.
    uvw[:, 0, :] = u
    uvw[:, 1, :] = v
    uvw[:, 2, :] = w

    # the Einstein summation is here configured to do the matrix
    # multiplication MM(i,l) = M1(i,j) * M2(j,k) * M3(k,l) on each slice h.
    MM = np.einsum('hij,hjk,hkl->hil', M1, M2, M3)

    # the Einstein summation is here configured to do the matrix
    # multiplication uvw_earth(i,m) = MM(i,l) * uvw(l,m) on each slice h.
    uvw_earth = np.einsum('hil,hlm->him', MM, uvw)

    # NOTE:
    # these last two executable statements run about a factor of 2
    # faster in the 10000 data packet performance tests versus combining
    # these operations into the one statement:
    #     uvw_earth = np.einsum('hij,hjk,hkl,hlm->him', M1, M2, M3, uvw)

    # break out the coordinate slices and return them
    uu = uvw_earth[:, 0, :]
    vv = uvw_earth[:, 1, :]
    ww = uvw_earth[:, 2, :]

    return (uu, vv, ww)


def magnetic_correction(theta, u, v):
    """
    Description:

        This function corrects velocity profiles for the magnetic variation
        (declination) at the measurement location.  The magnetic declination
        is obtained from the 2010 World Magnetic Model (WMM2010) provided by
        NOAA (see wmm_declination).

        This version handles 'vectorized' input variables without using for
        loops. It was specifically written to handle the case of a 1D array of
        theta values, theta=f(i), with corresponding sets of 'u' and 'v' values
        such that u=f(i,j) and v=f(i,j), where there are j 'u' and 'v' values
        for each theta(i).

    Implemented by:

        2014-04-04: Russell Desiderio. Initial code. This function is used to
                    calculate magnetic corrections by the functions contained
                    in this module instead of the function magnetic_correction
                    found in ion_functions.data.generic_functions.

    Usage:

        u_cor, v_cor = magnetic_correction(theta, u, v)

            where

        u_cor = eastward velocity profiles, in earth coordinates, with
            the correction for magnetic variation applied.
        v_cor = northward velocity profiles, in earth coordinates,
            with the correction for magnetic variation applied.

        theta = magnetic variation based on location (latitude, longitude and
            altitude) and date; units of theta are [degrees]
        u = uncorrected eastward velocity profiles in earth coordinates
        v = uncorrected northward velocity profiles in earth coordinates

    References:

        OOI (2012). Data Product Specification for Velocity Profile and Echo
            Intensity. Document Control Number 1341-00750.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00750_Data_Product_SPEC_VELPROF_OOI.pdf)

        OOI (2013). Data Product Specification for Turbulent Velocity Profile
            and Echo Intensity. Document Control Number 1341-00760.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00760_Data_Product_SPEC_VELPROF_OOI.pdf)
    """
    # force shapes of inputs to arrays
    theta = np.atleast_1d(theta)
    uv = np.atleast_2d(u)
    v = np.atleast_2d(v)

    theta_rad = np.radians(theta)
    cosT = np.cos(theta_rad)
    sinT = np.sin(theta_rad)

    M = np.array([[cosT, sinT],
                  [-sinT, cosT]])

    # roll axes so that the lead index represents data packet #.
    M = np.rollaxis(M, 2)

    # the coordinate system is 2D, so the middle dimension is sized at 2.
    uv = np.zeros((u.shape[0], 2, u.shape[1]))

    # pack the coordinates to be rotated into the appropriate slices
    uv[:, 0, :] = u
    uv[:, 1, :] = v

    # the Einstein summation is here configured to do the matrix
    # multiplication uv_cor(i,k) = M(i,j) * uv(j,k) on each slice h.
    uv_cor = np.einsum('hij,hjk->hik', M, uv)

    # the magnetically corrected u values are:
    u_cor = uv_cor[:, 0, :]

    # the magnetically corrected v values are:
    v_cor = uv_cor[:, 1, :]

    # return corrected u and v values
    return (u_cor, v_cor)
