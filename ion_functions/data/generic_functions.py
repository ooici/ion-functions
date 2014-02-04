#!/usr/bin/env python

"""
@package ion_functions.data.generic_functions
@file ion_functions/data/generic_functions.py
@author Christopher Mueller
@brief Module containing generic data-calculation functions.  Primarily
    used for calculating values in Parameter Functions
"""

# common imports
import datetime
import numpy as np
import numexpr as ne
import os
import pkg_resources
import time

from ion_functions.data.wmm import WMM


def magnetic_declination(lat, lon, ntp_timestamp, z=0, zflag=-1):
    """
    Description:

        Wrapper function, vectorizing inputs to wmm_declination

    Implemented by:

        2014-02-02: Christopher Wingard. Initial Code.
    """
    decln = np.vectorize(wmm_declination)
    mag_dec = decln(lat, lon, ntp_timestamp, z, zflag)
    return mag_dec


def set_wmm_model(year):
    """
    Based on year of sample, determine which WMM model coefficients file to
    use. Raises an exception if the file does not exist.
    """
    # set the WMM Coefficients file name based on year input.
    cof_file = 'WMM%4d.COF' % year

    # see if the file exists, if not raise an exception error.
    try:
        wmm_model = pkg_resources.resource_string(__name__, cof_file)
    except pkg_resources.ResolutionError as e:
        print("Error Type %s: Unable to find the WMM%4d.COF Coefficients file" % e, year)
    else:
        return wmm_model


def wmm_declination(lat, lon, ntp_timestamp, z=0, zflag=-1):
    """
    Description:

        Magnetic declination (a.k.a. magnetic variation) as a function
        of location and date from the World Magnetic Model (WMM).

        The magnetic declination correction is used to correct velocity vectors
        in several OOI data product transformations.

    Implemented by:

        2013-03-20: Stuart Pearce. Initial code.
        2013-06:    Luke Campbell. Implemented the WMM C code for speed
                    over the Python geomag library.
        2014-02-02: Christopher Wingard. Adjusted to allow for updates to the
                    WMM coefficients table (WMM.COF), that are updated every 5
                    years. Renamed to wmm_declination in order to keep
                    magnetic_declination name preserved in other modules.

    Usage:

        mag_dec = wmm_declination(lat,lon,z,ntp_timestamp,zflag=-1)

            where

        mag_dec = magnetic declination/variation value [degrees from N].
            Positive values are eastward, negative westward of North.
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

    Example:

        >>> lat = 45.0 # Location is deep water off of Oregon coast
        >>> lon = -128
        >>> z = -1000
        >>> ntp_timestamp = 3574792037.958   # 2013-04-12 14:47:17

        >>> mag_dec = wmm_declination(lat, lon, z,
        >>>                               ntp_timestamp, -1)
        >>> print mag_dec
        16.465045980896086

    References:

        Maus, S., S. Macmillan, S. McLean, B. Hamilton, A. Thomson, M. Nair,
            and C. Rollins, 2010, The US/UK World Magnetic Model for 2010-2015,
            NOAA Technical Report NESDIS/NGDC.
            http://www.ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml

    """
    # convert ntp timestamp to unix timestamp and then a datetime object
    unix_timestamp = ntp_timestamp - 2208988800.
    dates = datetime.datetime.utcfromtimestamp(unix_timestamp).date()

    # determine which WMM model to use.
    if dates.year >= 2010 and dates.year < 2015:
        wmm_model = set_wmm_model(2010)

    if dates.year >= 2015 and dates.year < 2020:
        wmm_model = set_wmm_model(2015)

    if dates.year >= 2020 and dates.year < 2025:
        wmm_model = set_wmm_model(2020)

    if dates.year >= 2025 and dates.year < 2030:
        wmm_model = set_wmm_model(2025)

    if dates.year >= 2030 and dates.year < 2035:
        wmm_model = set_wmm_model(2030)

    wmm = WMM(wmm_model)  # World Magnetic Model calculates magnetic declination

    # set the depth to negative for below sealevel (if needed) and convert from
    # meters to kilometers.
    z = ne.evaluate('z / 1000.')  # m -> km
    if z > 0 & zflag == -1:   # check that depth is a positive number first
        z = ne.evaluate('zflag * z')    # convert z to indicate depth

    # calculate the magnetic declination
    print 'here'
    mag_dec = wmm.declination(lat, lon, z, dates)

    return mag_dec


def ntp_to_unix_time(ntp_timestamp):
    """DEPRECATED. Use 2208988800 seconds as the offset between an NTP
    timestamp and a Unix timestamp instead.

    Description:

        Convert an NTP time stamp (epoch=1900-01-01) to a Unix timestamp
        (epoch=1970-01-01).  NOTE: Will become deprecated when coi-
        services time utilities module is brought into ion-functions.

    Implemented by:

        2013-04-12: Stuart Pearce. Initial code.

    Usage:

        unix_ts = ntp_to_unix_time(ntp_ts)

            where

        ntp_ts = NTP timestamp [seconds since 1900-01-01]
        unix_ts = Unix timestamp [seconds since 1970-01-01]

    Example:

        >>> ntp_ts = 3574792037.958

        >>> unix_ts = ntp_to_unix_time(ntp_ts)
        >>> unix_ts
        1365803237.9580002

    References:

        The NTP FAQ and HOWTO (2006). http://www.ntp.org/ntpfaq/
    """
    SYSTEM_EPOCH = datetime.date(*time.gmtime(0)[0:3])
    NTP_EPOCH = datetime.date(1900, 1, 1)
    NTP_DELTA = (SYSTEM_EPOCH - NTP_EPOCH).total_seconds()

    unix_timestamp = ne.evaluate('ntp_timestamp - NTP_DELTA')
    return unix_timestamp


def extract_parameter(in_array, index):
    """
    Description:

        Extracts and returns a single value from an array. Used, for example,
        to extract the L0 PH434SI from the array holding the 24 sets of 4 light
        measurements made by the Sunburst SAMI-II pH instrument (PHSEN).

    Implemented by:

        2013-04-19: Christopher Wingard. Initial code.

    Usage:

        out_value = extract_parameter(in_array, index)

            where

        out_value = the desired parameter.
        in_array = the input array holding the value.
        index = 0-based index in array to value.

    References:

        None.
    """
    out_value = in_array[index]
    return out_value


def bilinear_interpolation(x, y, points):
    '''
    Interpolate (x,y) from values associated with four points.

    The four points are a list of four triplets:  (x, y, value).
    The four points can be in any order.  They should form a rectangle.

        >>> bilinear_interpolation(12, 5.5,
        ...                        [(10, 4, 100),
        ...                         (20, 4, 200),
        ...                         (10, 6, 150),
        ...                         (20, 6, 300)])
        165.0
    '''
    # See formula at:  http://en.wikipedia.org/wiki/Bilinear_interpolation

    # order points by x, then by y
    pts = np.sort(points.view('f8,f8,f8'), order=['f0', 'f1'])
    (x1, y1, q11), (_x1, y2, q12), (x2, _y1, q21), (_x2, _y2, q22) = points

    if x1 != _x1 or x2 != _x2 or y1 != _y1 or y2 != _y2:
        raise ValueError('points do not form a rectangle')
    if not x1 <= x <= x2 or not y1 <= y <= y2:
        raise ValueError('(x, y) not within the rectangle')

    return (q11 * (x2 - x) * (y2 - y) +
            q21 * (x - x1) * (y2 - y) +
            q12 * (x2 - x) * (y - y1) +
            q22 * (x - x1) * (y - y1)) / ((x2 - x1) * (y2 - y1) + 0.0)


def error(x, y):
    return np.abs(x - y) / np.abs(y)
