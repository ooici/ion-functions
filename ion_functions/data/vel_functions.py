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

# NOTE:
#    The previous version of this module had each function return an
#    array of all fill values (-9999) if the lat or lon is invalid. This
#    should not occur. Really lat and lon should be checked and handled
#    in the QAQC functions, and only apply a fill value for the single
#    instance of invalid lat & lon rather than the whole array. I
#    believe they were coded here for testing purposes. The part that
#    returned the fill values now instead just raises the ValueError
#    exception.
#               -Stuart Pearce
#               (this message should be removed if/when the lats & lons
#               are checked in the QCQA functions)
from exceptions import ValueError


# wrapper functions for use in ION
def nobska_mag_corr_east(u, v, lat, lon, timestamp, z=0):
    """
    Corrects the eastward velocity from a VEL3D-B Nobska MAVS 4
    instrument for magnetic declination to produce an L1 VELPTTU-VLE
    OOI data product.

    Given a velocity vector with components u & v in the magnetic East
    and magnetic North directions respectively, this function calculates
    the magnetic declination for the location, depth, and time of the
    vector from the World Magnetic Model (WMM) and transforms the vector
    to a true Earth reference frame.

    This function is a wrapper around the function "vel_mag_correction".

    Usage:

        u_cor = nobska_mag_corr_east(u, v, lat, lon, ntp_timestamp, z=0)

            where

        u_cor = eastward velocity, in true Earth frame, with the
            correction for magnetic declination applied. [m/s]

        u = uncorrected eastward velocity in magnetic Earth frame. [cm/s]
        v = uncorrected northward velocity in magnetic Earth frame. [cm/s]
        lat = latitude of the instrument [decimal degrees].  East is
            positive, West negative.
        lon = longitude of the instrument [decimal degrees]. North
            is positive, South negative.
        ntp_timestamp = NTP time stamp from a data particle
            [secs since 1900-01-01].
        z = depth of instrument relative to sea level [meters].
            Positive values only. Default value is 0.

    References:

        OOI (2012). Data Product Specification for Turbulent Point Water
            Velocity. Document Control Number 1341-00781.
            https://alfresco.oceanobservatories.org/ (See: Company Home
            >> OOI >> Controlled >> 1000 System Level >>
            1341-00781_Data_Product_SPEC_VELPTTU_Nobska_OOI.pdf)
    """
   # Check for valid latitudes & longitudes
    if not valid_lat(lat) or not valid_lon(lon):
        # commented out next line according to NOTE above. -SP
        #return np.ones(u.shape, dtype=np.float) * -9999
        raise ValueError('Latitudes or Longitudes are not within the valid range!')

    # correct for magnetic declination
    u_cor = vel_mag_correction(u, v, lat, lon, timestamp, z)[0]
    u_cor = ne.evaluate('u_cor / 100.')  # convert from cm/s to m/s

    # return true compass referenced East velocity in m/s
    return u_cor


def nobska_mag_corr_north(u, v, lat, lon, timestamp, z=0):
    """
    Corrects the northward velocity from a VEL3D-B Nobska MAVS 4
    instrument for magnetic declination to produce an L1 VELPTTU-VLN
    OOI data product.

    Given a velocity vector with components u & v in the magnetic East
    and magnetic North directions respectively, this function calculates
    the magnetic declination for the location, depth, and time of the
    vector from the World Magnetic Model (WMM) and transforms the vector
    to a true Earth reference frame.

    This function is a wrapper around the function "vel_mag_correction".

    Usage:

        v_cor = nobska_mag_corr_north(u, v, lat, lon, ntp_timestamp, z)

            where

        v_cor = northward velocity, in true Earth frame, with the
            correction for magnetic declination applied. [m/s]

        u = uncorrected eastward velocity in magnetic Earth frame. [cm/s]
        v = uncorrected northward velocity in magnetic Earth frame. [cm/s]
        lat = latitude of the instrument [decimal degrees].  East is
            positive, West negative.
        lon = longitude of the instrument [decimal degrees]. North
            is positive, South negative.
        ntp_timestamp = NTP time stamp from a data particle
            [secs since 1900-01-01].
        z = depth of instrument relative to sealevel [meters].
            Positive values only. Default value is 0.

    References:

        OOI (2012). Data Product Specification for Turbulent Point Water
            Velocity. Document Control Number 1341-00781.
            https://alfresco.oceanobservatories.org/ (See: Company Home
            >> OOI >> Controlled >> 1000 System Level >>
            1341-00781_Data_Product_SPEC_VELPTTU_Nobska_OOI.pdf)
    """
   # check for valid latitudes & longitudes
    if not valid_lat(lat) or not valid_lon(lon):
        # commented out next line according to NOTE above. -SP
        #return np.ones(u.shape, dtype=np.float) * -9999
        raise ValueError('Latitudes or Longitudes are not within the valid range!')

    # correct for magnetic declination
    v_cor = vel_mag_correction(u, v, lat, lon, timestamp, z)[1]
    v_cor = ne.evaluate('v_cor / 100.')  # convert from cm/s to m/s

    # return true compass referenced North velocity in m/s
    return v_cor


def nobska_scale_up_vel(w):
    """
    Converts a Nobska MAVS-4 (VEL3D-B) vertical velocity measurement
    from cm/s to m/s

    Usage:

        w_mps = nobska_scale_up_vel(w_cmps)

            where

        w_mps = Output vertical velocity. [m/s]
        w_cmps = Input vertical velocity. [cm/s]

    References:

        OOI (2012). Data Product Specification for Turbulent Point Water
            Velocity. Document Control Number 1341-00781.
            https://alfresco.oceanobservatories.org/ (See: Company Home
            >> OOI >> Controlled >> 1000 System Level >>
            1341-00781_Data_Product_SPEC_VELPTTU_Nobska_OOI.pdf)
    """
    return w / 100.0


def nortek_mag_corr_east(u, v, lat, lon, timestamp, z=0.0):
    """
    Corrects the eastward velocity from VEL3D-CD Nortek Vector, VEL3D-K
    Nortek Aquadopp II, or VELPT Nortek Aquadopp instruments for
    magnetic declination to produce an L1 VELPTTU-VLE or an L1
    VELPTMN-VLE OOI data product.

    Given a velocity vector with components u & v in the magnetic East
    and magnetic North directions respectively, this function calculates
    the magnetic declination for the location, depth, and time of the
    vector from the World Magnetic Model (WMM) and transforms the vector
    to a true Earth reference frame.

    This function is a wrapper around the function "vel_mag_correction".

    Usage:

        u_cor = nortek_mag_corr_east(u, v, lat, lon, ntp_timestamp, z)

            where

        u_cor = eastward velocity , in true Earth frame, with the
            correction for magnetic declination applied. [m/s]

        u = uncorrected eastward velocity in magnetic Earth frame. [m/s]
        v = uncorrected northward velocity in magnetic Earth frame. [m/s]
        lat = latitude of the instrument [decimal degrees].  East is
            positive, West negative.
        lon = longitude of the instrument [decimal degrees]. North
            is positive, South negative.
        ntp_timestamp = NTP time stamp from a data particle
            [secs since 1900-01-01].
        z = depth of instrument relative to sealevel [meters].
            Positive values only. Default value is 0.

    References:

        OOI (2012). Data Product Specification for Turbulent Point Water
            Velocity. Document Control Number 1341-00780.
            https://alfresco.oceanobservatories.org/ (See: Company Home
            >> OOI >> Controlled >> 1000 System Level >>
            1341-00780_Data_Product_SPEC_VELPTTU_Nortek_OOI.pdf)
    """
   # check for valid latitudes & longitudes
    if not valid_lat(lat) or not valid_lon(lon):
        # commented out next line according to NOTE above. -SP
        #return np.ones(u.shape, dtype=np.float) * -9999
        raise ValueError('Latitudes or Longitudes are not within the valid range!')

    # correct for magnetic declination
    u_cor = vel_mag_correction(u, v, lat, lon, timestamp, z)[0]

    # return true compass referenced East velocity in m/s
    return u_cor


def nortek_mag_corr_north(u, v, lat, lon, timestamp, z=0.0):
    """
    Corrects the northward velocity from VEL3D-CD Nortek Vector, VEL3D-K
    Nortek Aquadopp II, and VELPT Nortek Aquadopp instruments for
    magnetic declination to produce an L1 VELPTTU-VLN or an L1
    VELPTMN-VLN OOI data product.

    Given a velocity vector with components u & v in the magnetic East
    and magnetic North directions respectively, this function calculates
    the magnetic declination for the location, depth, and time of thel
    vector from the World Magnetic Model (WMM) and transforms the vector
    to a true Earth reference frame.

    This function is a wrapper around the function "vel_mag_correction".

    Usage:

        v_cor = nortek_mag_corr_north(u, v, lat, lon, ntp_timestamp, z)

            where

        v_cor = northward velocity, in true Earth frame, with the
            correction for magnetic declination applied. [m/s]

        u = uncorrected eastward velocity in magnetic Earth frame. [m/s]
        v = uncorrected northward velocity in magnetic Earth frame. [m/s]
        lat = latitude of the instrument [decimal degrees].  East is
            positive, West negative.
        lon = longitude of the instrument [decimal degrees]. North
            is positive, South negative.
        ntp_timestamp = NTP time stamp from a data particle
            [secs since 1900-01-01].
        z = depth of instrument relative to sealevel [meters].
            Positive values only. Default value is 0.

    References:

        OOI (2012). Data Product Specification for Turbulent Point Water
            Velocity. Document Control Number 1341-00780.
            https://alfresco.oceanobservatories.org/ (See: Company Home
            >> OOI >> Controlled >> 1000 System Level >>
            1341-00780_Data_Product_SPEC_VELPTTU_Nortek_OOI.pdf)
    """
   # check for valid latitudes & longitudes
    if not valid_lat(lat) or not valid_lon(lon):
        # commented out next line according to NOTE above. -SP
        #return np.ones(u.shape, dtype=np.float) * -9999
        raise ValueError('Latitudes or Longitudes are not within the valid range!')

    # correct for magnetic declination
    v_cor = vel_mag_correction(u, v, lat, lon, timestamp, z)[1]

    # return true compass referenced North velocity in m/s
    return v_cor


def nortek_up_vel(w):
    """
    Returns a Nortek instrument (VEL3D-C,D,K or VELPT) vertical velocity
    measurement. This function is an identity function to return the
    same value, but is required because of how OOINet and OOI have
    designated the system to work.

    Usage:

        w = nortek_up_vel(w)

            where

        w = Output and Input vertical velocity.  [m/s]

    References:

        OOI (2012). Data Product Specification for Turbulent Point Water
            Velocity. Document Control Number 1341-00781.
            https://alfresco.oceanobservatories.org/ (See: Company Home
            >> OOI >> Controlled >> 1000 System Level >>
            1341-00781_Data_Product_SPEC_VELPTTU_Nobska_OOI.pdf)
    """
    return w


def vel3dk_mag_corr_east(u, v, lat, lon, timestamp, z, Vscale):
    """
    Corrects the eastward velocity from VEL3D-K (Nortek Aquadopp II)
    instruments for magnetic declination to produce the L1 VELPTTU-VLE
    OOI data product.

    Takes a velocity in integer counts from a VEL3D-K (Aquadopp II on a McLane
    Profiler(MMP)) with the provided Vscale parameter from an MMP
    A#####.DEC binary data file.

    Given a velocity vector with components u & v in the magnetic East
    and magnetic North directions respectively, this function calculates
    the magnetic declination for the location, depth, and time of the
    vector from the World Magnetic Model (WMM) and transforms the vector
    to a true Earth reference frame.

    This function is a wrapper around the function "vel_mag_correction".

    Usage:
        u_cor = nortek_mag_corr_east(u, v, lat, lon, ntp_timestamp, z, Vscale)

            where

        u_cor = floating point eastward velocity, in true Earth frame,
            with the correction for magnetic declination applied. [m/s]

        u = uncorrected eastward velocity in magnetic Earth
            frame.  Integer scaled by 10^Vscale. [scaled integer distance/s]
        v = uncorrected northward velocity in magnetic Earth
            frame.  Integer scaled by 10^Vscale. [scaled integer distance/s]
        lat = latitude of the instrument [decimal degrees].  East is
            positive, West negative.
        lon = longitude of the instrument [decimal degrees]. North
            is positive, South negative.
        ntp_timestamp = NTP time stamp from a data particle
            [secs since 1900-01-01].
        z = depth of instrument relative to sealevel [meters].
            Positive values only. Default value is 0.
        Vscale = velocity scaling exponent factor.

    References:

        VEL3D-K IDD (2014) (No DPS as of 2014-03-03)
        https://confluence.oceanobservatories.org/display/instruments/
        VEL3D-K__stc_imodem+-+Telemetered

        OOI (2012). Data Product Specification for Turbulent Point Water
            Velocity. Document Control Number 1341-00780.
            https://alfresco.oceanobservatories.org/ (See: Company Home
            >> OOI >> Controlled >> 1000 System Level >>
            1341-00780_Data_Product_SPEC_VELPTTU_Nortek_OOI.pdf)
    """
    # convert from scaled, integer distance/s (as received from the
    # binary data file) to floating point m/s using the Vscale parameter
    # from the MMP binary data file A#####.DEC
    u = u * 10**Vscale
    v = v * 10**Vscale

   # check for valid latitudes & longitudes
    if not valid_lat(lat) or not valid_lon(lon):
        raise ValueError('Latitudes or Longitudes are not within the valid range!')

    # correct for magnetic declination
    u_cor = vel_mag_correction(u, v, lat, lon, timestamp, z)[0]

    # return true compass referenced East velocity in m/s
    return u_cor


def vel3dk_mag_corr_north(u, v, lat, lon, timestamp, z, Vscale):
    """
    Corrects the eastward velocity from VEL3D-K (Nortek Aquadopp II)
    instruments for magnetic declination to produce the L1 VELPTTU-VLE
    OOI data product.

    Takes a velocity in integer counts from a VEL3D-K (Aquadopp II on a McLane
    Profiler(MMP)) with the provided Vscale parameter from an MMP
    A#####.DEC binary data file.

    Given a velocity vector with components u & v in the magnetic East
    and magnetic North directions respectively, this function calculates
    the magnetic declination for the location, depth, and time of the
    vector from the World Magnetic Model (WMM) and transforms the vector
    to a true Earth reference frame.

    This function is a wrapper around the function "vel_mag_correction".

    Usage:
        u_cor = nortek_mag_corr_east(u, v, lat, lon, ntp_timestamp, z, Vscale)

            where

        u_cor = floating point eastward velocity, in true Earth frame,
            with the correction for magnetic declination applied. [m/s]

        u = uncorrected eastward velocity in magnetic Earth
            frame.  Integer scaled by 10^Vscale. [scaled integer distance/s]
        v = uncorrected northward velocity in magnetic Earth
            frame.  Integer scaled by 10^Vscale. [scaled integer distance/s]
        lat = latitude of the instrument [decimal degrees].  East is
            positive, West negative.
        lon = longitude of the instrument [decimal degrees]. North
            is positive, South negative.
        ntp_timestamp = NTP time stamp from a data particle
            [secs since 1900-01-01].
        z = depth of instrument relative to sealevel [meters].
            Positive values only. Default value is 0.
        Vscale = velocity scaling exponent factor.

    References:

        VEL3D-K IDD (2014) (No DPS as of 2014-03-03)
        https://confluence.oceanobservatories.org/display/instruments/
        VEL3D-K__stc_imodem+-+Telemetered

        OOI (2012). Data Product Specification for Turbulent Point Water
            Velocity. Document Control Number 1341-00780.
            https://alfresco.oceanobservatories.org/ (See: Company Home
            >> OOI >> Controlled >> 1000 System Level >>
            1341-00780_Data_Product_SPEC_VELPTTU_Nortek_OOI.pdf)
    """
    # convert from scaled, integer distance/s (as received from the
    # binary data file) to floating point m/s using the Vscale parameter
    # from the MMP binary data file A#####.DEC
    u = u * 10**Vscale
    v = v * 10**Vscale

   # check for valid latitudes & longitudes
    if not valid_lat(lat) or not valid_lon(lon):
        raise ValueError('Latitudes or Longitudes are not within the valid range!')

    # correct for magnetic declination
    v_cor = vel_mag_correction(u, v, lat, lon, timestamp, z)[1]

    # return true compass referenced East velocity in m/s
    return v_cor


def vel3dk_scale_up_vel(w, Vscale):
    """
    Takes an integer vertical velocity in generic distance per second
    units from a VEL3D-K (Aquadopp II on a McLane Profiler(MMP)) with
    the provided Vscale parameter from an MMP A#####.DEC binary data
    file to scale the velocity to a floating point in m/s.

    Usage:
        w_mps = vel3dk_scale_up_vel(w, Vscale)

            where

        w_mps = floating point vertical velocity. [m/s]

        w = integer vertical velocity. Integer scaled by 10^Vscale.
            [scaled integer distance/s]
        Vscale = velocity scaling exponent factor.

    References:

        VEL3D-K IDD (2014) (No DPS as of 2014-03-03)
        https://confluence.oceanobservatories.org/display/instruments/
        VEL3D-K__stc_imodem+-+Telemetered
    """
    # convert from scaled, integer distance/s (as received from the
    # binary data file) to floating point m/s using the Vscale parameter
    # from the MMP binary data file A#####.DEC
    u = u * 10**Vscale
    v = v * 10**Vscale

   # check for valid latitudes & longitudes
    if not valid_lat(lat) or not valid_lon(lon):
        raise ValueError('Latitudes or Longitudes are not within the valid range!')

    # correct for magnetic declination
    v_cor = vel_mag_correction(u, v, lat, lon, timestamp, z)[1]

    # return true compass referenced East velocity in m/s
    return v_cor


##### Sub functions #####
# the main sub-function
def vel_mag_correction(u, v, lat, lon, ntp_timestamp, z=0.0, zflag=-1):
    """
    Description:

        Given a velocity vector U, measured in a sensor frame that is
        referenced to Earth's magnetic field, with components u & v in
        the magnetic East and magnetic North directions respectively;
        vel_mag_correction transforms U to true Earth referenced
        directions by a rotation that removes the magnetic declination.
        Magnetic Declination, theta(x,y,z,t), is the azimuthal angular
        offset between true North and magnetic North as a function of
        Earth referenced location (latitude, longitude, & height/depth)
        and time. Magnetic declination is estimated from the World
        Magnetic Model (WMM) using the location and time of the vector.

    Usage:

        u_cor, v_cor = vel_mag_correction(u, v, lat, lon, ntp_timestamp, z, zflag)

            where

        u_cor = eastward velocity, in true Earth frame, with the
            correction for magnetic declination applied.
        v_cor = northward velocity, in true Earth frame, with the
            correction for magnetic declination applied.

        u = uncorrected eastward velocity in magnetic Earth frame.
        v = uncorrected northward velocity in magnetic Earth frame.
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
            is the default, because Oceanography!

    Implemented by:

        2013-04-17: Stuart Pearce. Initial code.
        2013-04-24: Stuart Pearce. Changed to be general for all velocity
                    instruments.
        2014-02-05: Christopher Wingard. Edited to use magnetic corrections in
                    the generic_functions module.

    References:

        OOI (2012). Data Product Specification for Turbulent Point Water
            Velocity. Document Control Number 1341-00781.
            https://alfresco.oceanobservatories.org/ (See: Company Home
            >> OOI >> Controlled >> 1000 System Level >>
            1341-00781_Data_Product_SPEC_VELPTTU_Nobska_OOI.pdf)

        OOI (2012). Data Product Specification for Turbulent Point Water
            Velocity. Document Control Number 1341-00780.
            https://alfresco.oceanobservatories.org/ (See: Company Home
            >> OOI >> Controlled >> 1000 System Level >>
            1341-00780_Data_Product_SPEC_VELPTTU_Nortek_OOI.pdf)
    """
    # retrieve the magnetic declination
    theta = magnetic_declination(lat, lon, ntp_timestamp, z, zflag)

    # apply the magnetic declination correction
    magvar = np.vectorize(magnetic_correction)
    u_cor, v_cor = magvar(theta, u, v)

    return u_cor, v_cor


# helper sub-functions
def valid_lat(lat):
    """valid_lat(lat) -> boolean

    Checks if inputs are valid latitude values.
    Returns True if value is between -90 and 90,
    False otherwise.
    """
    if isinstance(lat, np.ndarray):
        if np.any(lat > 90) or np.any(lat < -90):
            return False
        return True
    else:
        return -90 <= lat and lat <= 90


def valid_lon(lon):
    """valid_lon(lon) -> boolean

    Checks if inputs are valid longitude values.
    Returns True if value is between -180 and 180,
    False otherwise.
    """
    if isinstance(lon, np.ndarray):
        if np.any(lon > 180) or np.any(lon < -180):
            return False
        return True
    else:
        return -180 <= lon and lon <= 180
