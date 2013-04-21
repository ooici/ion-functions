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
import time

import geomag

# Example function from ctd_functions.py
def magnetic_declination(lat, lon, z, ntp_timestamp, zflag=-1):
    """
    Description:

        Magnetic declination (a.k.a. magnetic variation) as a function
        of location and date from the World Magnetic Model (WMM). Uses
        the geomag Python library implementation of the WMM. Declination
        is used in several OOI data product transformations.

    Implemented by:

        2013-03-20: Stuart Pearce. Initial code.

    Usage:

        mag_dec = magnetic_declination(lat,lon,z,ntp_timestamp,zflag=-1)

            where

        mag_dec = magnetic declination/variation value [degrees from N].
            Positive values are eastward, negative westward of North.
        lat = latitude of the instrument [decimal degrees].  East is
            positive, West negative.
        lon = longitude of the instrument [decimal degrees]. North
            is positive, South negative.
        z = depth or height of instrument relative to sealevel [meters].
            Positive values only.
        ntp_timestamp = NTP time stamp from a data particle
            [secs since 1900-01-01].
        zflag = indicates whether to use z as a depth or height relative
            to sealevel. -1=depth (i.e. -z) and 1=height (i.e. +z). -1
            is the default

    Example:

        >>> lat = 45.0 # Location is deep water off of Oregon coast
        >>> lon = -128
        >>> z = -1000
        >>> ntp_timestamp = 3574792037.958   # 2013-04-12 14:47:17

        >>> mag_dec = magnetic_declination(lat, lon, z,
        >>>                               ntp_timestamp, -1)
        >>> print mag_dec
        16.465045980896086

    References:
    
        World Magnetic Model (2010). http://www.ngdc.noaa.gov/geomag/WMM
        /DoDWMM.shtml
    """

    # convert ntp timestamp to unix timestamp and then a datetime object
    unix_timestamp = ntp_to_unix_time(ntp_timestamp)
    
    # the data timestamp is in UTC, and while the input to the WMM does
    # not specify, reasonably it should handle time as a UTC time rather
    # than a local time. If the WMM does expect a local time, since a
    # date is required by the geomag library, the difference of discrete
    # day timesteps results in an average error that is much smaller
    # than the uncertainty of almost all compasses and so local time
    # versus UTC can be ignored.
    datestamp = datetime.datetime.utcfromtimestamp(unix_timestamp).date()
    
    # give the z value the proper vector direction (i.e negative down)
    z = z*zflag
    
    # geomag python library requires depth in feet (tisk tisk) and then
    # in the library converts them right back to meters!  Ridiculous!
    # ALL science shall be in GOD's units ... SI units!!!
    z *= 3.28084  # m*3.28084 ft/m = ft
    
    # TODO: Handling depths; should I use Udunits? Zflag or other?
    
    # the magnetic declination at a given location & time
    mag_dec = geomag.declination(lat, lon, z, datestamp)
    return mag_dec


def ntp_to_unix_time(ntp_timestamp):
    """
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
    
    unix_timestamp = ntp_timestamp - NTP_DELTA
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