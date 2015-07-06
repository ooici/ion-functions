#!/usr/bin/env python
"""
@package ion_functions.data.adcp_functions
@file ion_functions/data/adcp_functions.py
@author Christopher Wingard, Russell Desiderio, Craig Risien
@brief Module containing ADCP related data-calculations.
"""
import numpy as np
from ion_functions.data.generic_functions import magnetic_declination
from ion_functions.data.generic_functions import replace_fill_with_nan

# instrument fill value unprocessed by CI
# (bad beam velocity sentinel output by tRDI ADCP instruments)
ADCP_FILLVALUE = -32768

"""
      **** For instruments programmed in beam coordinates:
           (ADCPS-I,K;  ADCPT-B,D,E)
      adcp_beam_eastward -- calculates VELPROF-VLE_L1
      adcp_beam_northward -- calculates VELPROF-VLN_L1
      adcp_beam_vertical -- calculates VELPROF-VLU_L1
      adcp_beam_error -- calculates VELPROF-ERR_L1

      **** For instruments programmed in earth coordinates:
           (ADCPA;  ADCPS-J,L,N; ADCPT-C,F,G,M)
      adcp_earth_eastward -- calculates VELPROF-VLE_L1
      adcp_earth_northward -- calculates VELPROF-VLN_L1
      adcp_earth_vertical -- calculates VELPROF-VLU_L1
      adcp_earth_error -- calculates VELPROF-ERR_L1

      **** For the VADCP programmed in beam coordinates:
      vadcp_beam_eastward -- calculates VELTURB-VLE_L1
      vadcp_beam_northward -- calculates VELTURB-VLN_L1
      vadcp_beam_vertical_true -- calculates VELTURB-VLU-5BM_L1
      vadcp_beam_vertical_est -- calculates VELTURB-VLU-4BM_L1
      vadcp_beam_error -- calculates VELTURB-ERR_L1

      **** For all tRDI ADCP instruments:
      adcp_backscatter -- calculates ECHOINT-B1_L1,
                          calculates ECHOINT-B2_L1,
                          calculates ECHOINT-B3_L1,
                          calculates ECHOINT-B4_L1.

      **** Base functions used by above functions
      adcp_beam2ins -- applies the beam to instrument transform using a 4
            beam solution for instruments programmed in beam coordinates
      adcp_ins2earth -- applies the instrument to Earth transform for all
            instruments originally programmed in beam coordinates.
      magnetic_correction -- corrects horizontal velocities for the magnetic
            variation (declination) at the measurement location.

      **** Supplementary functions to calculate velocity bin depths:
      adcp_bin_depths -- calculates bin depths for the pd0 output format
                         (virtually all tRDI ADCPs deployed by OOI); uses
                         TEOS-10 functions p_from_z and enthalpy_SSO_0_p.
      adcp_bin_depths_pd8 -- calculates bin depths for the pd8 output format,
                             assuming that (1) the ADCP operator recorded the
                             necessary input variables and (2) these are somehow
                             entered into the CI system.

"""


# Wrapper functions to create the VELPROF L1 data products for instruments
# programmed in beam coordinates by RSN (ADCPS-I,K and ADCPT-B,D,E)
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
        2014-04-04: Russell Desiderio. Optimized code performance by replacing
                    the for loops previously used to calculate 2D and 3D
                    vectorized coordinate transformations with calls to
                    np.einsum (numpy Einstein summation function).
        2014-06-25: Christopher Wingard. Edited to account for units of
                    heading, pitch, roll and depth
        2015-06-10: Russell Desiderio.
                    (a) moved the conditioning of input beam velocities to adcp_beam2inst.
                    (b) moved the conditioning of compass readings to adcp_inst2earth.
                    (c) removed the depth dependence from the magnetic declination.

    Usage:

        uu_cor = adcp_beam_eastward(b1, b2, b3, b4, h, p, r, vf, lat, lon, z, dt)

            where

        uu_corr = east velocity profiles in Earth coordinates corrected for the
                  magnetic declination (VELPROF-VLE_L1) [m s-1]

        b1 = "beam 1" velocity profiles in beam coordinates (VELPROF-B1_L0) [mm s-1]
        b2 = "beam 2" velocity profiles in beam coordinates (VELPROF-B2_L0) [mm s-1]
        b3 = "beam 3" velocity profiles in beam coordinates (VELPROF-B3_L0) [mm s-1]
        b4 = "beam 4" velocity profiles in beam coordinates (VELPROF-B4_L0) [mm s-1]
        h = instrument's uncorrected magnetic heading [cdegrees]
        p = instrument pitch [cdegrees]
        r = instrument roll [cdegrees]
        vf = instrument's vertical orientation (0 = downward looking and
            1 = upward looking)
        lat = instrument's deployment latitude [decimal degrees]
        lon = instrument's deployment longitude [decimal degrees]
        z = instrument's pressure sensor reading (depth) [daPa]
        dt = sample date and time value [seconds since 1900-01-01]
    """
    # force shapes of inputs to arrays of the correct dimensions
    #z = np.atleast_1d(z) / 1000.  # scale daPa depth input to dbar
    #z = z * 1.019716  # use a simple approximation to calculate depth in m
    lat = np.atleast_1d(lat)
    lon = np.atleast_1d(lon)
    dt = np.atleast_1d(dt)

    # compute the beam to instrument transform
    u, v, w, eee = adcp_beam2ins(b1, b2, b3, b4)
    #print eee

    # compute the instrument to earth beam transform
    uu, vv, _ = adcp_ins2earth(u, v, w, h, p, r, vf)

    # compute the magnetic variation, and ...
    theta = magnetic_declination(lat, lon, dt)

    # ... correct for it
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
        2014-03-28: Russell Desiderio. Corrected documentation only.
        2014-04-04: Russell Desiderio. Optimized code performance by replacing
                    the for loops previously used to calculate 2D and 3D
                    vectorized coordinate transformations with calls to
                    np.einsum (numpy Einstein summation function).
        2014-06-25: Christopher Wingard. Edited to account for units of
                    heading, pitch, roll and depth
        2015-06-10: Russell Desiderio.
                    (a) moved the conditioning of input beam velocities to adcp_beam2inst.
                    (b) moved the conditioning of compass readings to adcp_inst2earth.
                    (c) removed the depth dependence from the magnetic declination.

    Usage:

        vv_cor = adcp_beam_northward(b1, b2, b3, b4, h, p, r, vf, lat, lon, z, dt)

            where

        vv_corr = north velocity profiles in Earth coordinates corrected for the
                  magnetic declination (VELPROF-VLN_L1) [m s-1]

        b1 = "beam 1" velocity profiles in beam coordinates (VELPROF-B1_L0) [mm s-1]
        b2 = "beam 2" velocity profiles in beam coordinates (VELPROF-B2_L0) [mm s-1]
        b3 = "beam 3" velocity profiles in beam coordinates (VELPROF-B3_L0) [mm s-1]
        b4 = "beam 4" velocity profiles in beam coordinates (VELPROF-B4_L0) [mm s-1]
        h = instrument's uncorrected magnetic heading [cdegrees]
        p = instrument pitch [cdegrees]
        r = instrument roll [cdegrees]
        vf = instrument's vertical orientation (0 = downward looking and
            1 = upward looking)
        lat = instrument's deployment latitude [decimal degrees]
        lon = instrument's deployment longitude [decimal degrees]
        z = instrument's pressure sensor reading (depth) [daPa]
        dt = sample date and time value [seconds since 1900-01-01]
    """
    # force shapes of inputs to arrays of the correct dimensions
    #z = np.atleast_1d(z) / 1000.  # scale daPa depth input to dbar
    #z = z * 1.019716  # use a simple approximation to calculate depth in m
    lat = np.atleast_1d(lat)
    lon = np.atleast_1d(lon)
    dt = np.atleast_1d(dt)

    # compute the beam to instrument transform
    u, v, w, _ = adcp_beam2ins(b1, b2, b3, b4)

    # compute the instrument to earth beam transform
    uu, vv, _ = adcp_ins2earth(u, v, w, h, p, r, vf)

    # compute the magnetic variation, and ...
    theta = magnetic_declination(lat, lon, dt)

    # ... correct for it
    _, vv_cor = magnetic_correction(theta, uu, vv)

    # scale velocity to m/s
    vv_cor = vv_cor / 1000.  # mm/s -> m/s

    # return the Northward Velocity Profile
    return vv_cor


def adcp_beam_vertical(b1, b2, b3, b4, h, p, r, vf):
    """
    Description:

        Wrapper function to compute the Upward Velocity Profile (VELPROF-VLU)
        from beam coordinate transformed velocity profiles as defined in the
        Data Product Specification for Velocity Profile and Echo Intensity -
        DCN 1341-00750.

    Implemented by:

        2013-04-10: Christopher Wingard. Initial code.
        2014-02-03: Christopher Wingard. Formatting and adjusting to use
                    magnetic declination values calculated using the WMM 2010.
        2014-04-04: Russell Desiderio. Optimized code performance by replacing
                    the for loops previously used to calculate 2D and 3D
                    vectorized coordinate transformations with calls to
                    np.einsum (numpy Einstein summation function).
        2014-06-25: Christopher Wingard. Edited to account for units of
                    heading, pitch, roll and depth
        2015-06-10: Russell Desiderio.
                    (a) moved the conditioning of input beam velocities to adcp_beam2inst.
                    (b) moved the conditioning of compass readings to adcp_inst2earth.

    Usage:

        ww_cor = adcp_beam_vertical(b1, b2, b3, b4, h, p, r, vf)

            where

        ww_cor = vertical velocity profiles (VELPROF-VLU_L1) [m s-1]

        b1 = "beam 1" velocity profiles in beam coordinates (VELPROF-B1_L0) [mm s-1]
        b2 = "beam 2" velocity profiles in beam coordinates (VELPROF-B2_L0) [mm s-1]
        b3 = "beam 3" velocity profiles in beam coordinates (VELPROF-B3_L0) [mm s-1]
        b4 = "beam 4" velocity profiles in beam coordinates (VELPROF-B4_L0) [mm s-1]
        h = instrument's uncorrected magnetic heading [cdegrees]
        p = instrument pitch [cdegrees]
        r = instrument roll [cdegrees]
        vf = instrument's vertical orientation (0 = downward looking and
            1 = upward looking)
    """
    # compute the beam to instrument transform
    u, v, w, _ = adcp_beam2ins(b1, b2, b3, b4)

    # compute the instrument to earth beam transform
    _, _, ww = adcp_ins2earth(u, v, w, h, p, r, vf)

    # scale upward velocity to m/s
    ww = ww / 1000.  # mm/s -> m/s

    # return the Upward Velocity Profile
    return ww


def adcp_beam_error(b1, b2, b3, b4):
    """
    Description:

        Wrapper function to compute the Error Velocity Profile (VELPROF-ERR)
        from beam coordinate transformed velocity profiles as defined in the
        Data Product Specification for Velocity Profile and Echo Intensity -
        DCN 1341-00750.

    Implemented by:

        2013-04-10: Christopher Wingard. Initial code.
        2015-06-10: Russell Desiderio.
                    Moved the conditioning of input beam velocities to adcp_beam2inst.

    Usage:

        ww_cor = adcp_beam_error(b1, b2, b3, b4)

            where

        e = Error velocity profiles (VELPROF-ERR_L1) [m s-1]

        b1 = "beam 1" velocity profiles in beam coordinates (VELPROF-B1_L0) [mm s-1]
        b2 = "beam 2" velocity profiles in beam coordinates (VELPROF-B2_L0) [mm s-1]
        b3 = "beam 3" velocity profiles in beam coordinates (VELPROF-B3_L0) [mm s-1]
        b4 = "beam 4" velocity profiles in beam coordinates (VELPROF-B4_L0) [mm s-1]
    """
    # compute the beam to instrument transform
    _, _, _, e = adcp_beam2ins(b1, b2, b3, b4)

    # scale error velocity to m/s
    e = e / 1000.   # mm/s

    # return the Error Velocity Profile
    return e


# Wrapper functions to create the VELPROF L1 data products for instruments
# programmed in Earth coordinates by CGSN (Pioneer and Endurance) (ADCPA,
# ADCPS-J,L,N and ADCPT-C,F,G,M)
def adcp_earth_eastward(u, v, z, lat, lon, dt):
    """
    Description:

        Wrapper function to compute the Eastward Velocity Profile (VELPROF-VLE)
        from Earth coordinate transformed velocity profiles as defined in the
        Data Product Specification for Velocity Profile and Echo Intensity -
        DCN 1341-00750.

    Implemented by:

        2013-04-10: Christopher Wingard. Initial code.
        2014-02-03: Christopher Wingard. Formatting and adjusting to use
                    magnetic declination values calculated use the WMM 2010.
        2014-04-04: Russell Desiderio. Optimized code performance by replacing
                    the for loops previously used to calculate 2D and 3D
                    vectorized coordinate transformations with calls to
                    np.einsum (numpy Einstein summation function).
        2014-06-25: Christopher Wingard. Edited to account for units of
                    heading, pitch, roll and depth
        2015-06-10: Russell Desiderio.
                    Removed the depth dependence from the magnetic declination.
        2015-06-25: Russell Desiderio. Incorporated int fillvalue -> Nan.

    Usage:

        uu_cor = adcp_earth_eastward(u, v, z, lat, lon, dt)

            where

        uu_cor = eastward velocity profiles in Earth coordinates corrected for
                 the magnetic declination (VELPROF-VLE_L1) [m s-1]

        u = Eastward velocity profiles (VELPROF-VLE_L0) [mm s-1]
        v = Northward velocity profiles (VELPROF-VLN_L0) [mm s-1]
        z = instrument's pressure sensor reading (depth) [daPa]
        lat = instrument's deployment latitude [decimal degrees]
        lon = instrument's deployment longitude [decimal degrees]
        dt = sample date and time value [seconds since 1900-01-01]
    """
    # force shapes of inputs to arrays
    u = np.atleast_2d(u)
    v = np.atleast_2d(v)

    # on input, the elements of u and v are of type int.
    u, v = replace_fill_with_nan(ADCP_FILLVALUE, u, v)

    #z = np.atleast_1d(z) / 1000.  # scale daPa depth input to dbar
    #z = z * 1.019716  # use a simple approximation to calculate depth in m
    lat = np.atleast_1d(lat)
    lon = np.atleast_1d(lon)
    dt = np.atleast_1d(dt)

    # compute the magnetic variation, and ...
    theta = magnetic_declination(lat, lon, dt)

    # ... correct for it
    uu_cor, _ = magnetic_correction(theta, u, v)

    # scale velocity to m/s
    uu_cor = uu_cor / 1000.  # mm/s -> m/s

    # return the Eastward Velocity Profile
    return uu_cor


def adcp_earth_northward(u, v, z, lat, lon, dt):
    """
    Description:

        Wrapper function to compute the Northward Velocity Profile (VELPROF-VLN)
        from Earth coordinate transformed velocity profiles as defined in the
        Data Product Specification for Velocity Profile and Echo Intensity -
        DCN 1341-00750.

    Implemented by:

        2013-04-10: Christopher Wingard. Initial code.
        2014-02-03: Christopher Wingard. Formatting and adjusting to use
                    magnetic declination values calculated use the WMM 2010.
        2014-04-04: Russell Desiderio. Optimized code performance by replacing
                    the for loops previously used to calculate 2D and 3D
                    vectorized coordinate transformations with calls to
                    np.einsum (numpy Einstein summation function).
        2014-06-25: Christopher Wingard. Edited to account for units of
                    heading, pitch, roll and depth
        2015-06-10: Russell Desiderio.
                    Removed the depth dependence from the magnetic declination.
        2015-06-25: Russell Desiderio. Incorporated int fillvalue -> Nan.

    Usage:

        vv_cor = adcp_earth_northward(u, v, z, lat, lon, dt)

            where

        vv_cor = northward velocity profiles in Earth coordinates corrected for
                 the magnetic declination (VELPROF-VLN_L1) [m s-1]

        u = Eastward velocity profiles (VELPROF-VLE_L0) [mm s-1]
        v = Northward velocity profiles (VELPROF-VLN_L0) [mm s-1]
        z = instrument's pressure sensor reading (depth) [daPa]
        lat = instrument's deployment latitude [decimal degrees]
        lon = instrument's deployment longitude [decimal degrees]
        dt = sample date and time value [seconds since 1900-01-01]
    """
    # force shapes of inputs to arrays
    u = np.atleast_2d(u)
    v = np.atleast_2d(v)

    # on input, the elements of u and v are of type int.
    u, v = replace_fill_with_nan(ADCP_FILLVALUE, u, v)

    #z = np.atleast_1d(z) / 1000.  # scale daPa depth input to dbar
    #z = z * 1.019716  # use a simple approximation to calculate depth in m
    lat = np.atleast_1d(lat)
    lon = np.atleast_1d(lon)
    dt = np.atleast_1d(dt)

    # compute the magnetic variation, and ...
    theta = magnetic_declination(lat, lon, dt)

    # ... correct for it
    _, vv_cor = magnetic_correction(theta, u, v)

    # scale velocity to m/s
    vv_cor = vv_cor / 1000.  # mm/s -> m/s

    # return the Northward Velocity Profile
    return vv_cor


def adcp_earth_vertical(w):
    """
    Description:

        Wrapper function to compute the Upward Velocity Profile (VELPROF-VLU)
        from Earth coordinate transformed velocity profiles as defined in the
        Data Product Specification for Velocity Profile and Echo Intensity -
        DCN 1341-00750.

    Implemented by:

        2014-06-25: Christopher Wingard. Initial code.
        2015-06-25: Russell Desiderio. Incorporated int fillvalue -> Nan.

    Usage:

        w_scl = adcp_earth_vertical(w)

            where

        w_scl = scaled upward velocity profiles in Earth coordinates
                (VELPROF-VLN_L1) [m s-1]

        w = upward velocity profiles (VELPROF-VLU_L0) [mm s-1]
    """
    w = replace_fill_with_nan(ADCP_FILLVALUE, w)

    # scale velocity to m/s
    w_scl = w / 1000.  # mm/s -> m/s

    # return the Upward Velocity Profile
    return w_scl


def adcp_earth_error(e):
    """
    Description:

        Wrapper function to compute the Error Velocity Profile (VELPROF-ERR)
        from Earth coordinate transformed velocity profiles as defined in the
        Data Product Specification for Velocity Profile and Echo Intensity -
        DCN 1341-00750.

    Implemented by:

        2014-06-25: Christopher Wingard. Initial code.
        2015-06-25: Russell Desiderio. Incorporated int fillvalue -> Nan.

    Usage:

        e_scl = adcp_earth_vertical(w)

            where

        e_scl = scaled error velocity profiles in Earth coordinates
                (VELPROF-ERR_L1) [m s-1]

        e = error velocity profiles (VELPROF-ERR_L0) [mm s-1]
    """
    e = replace_fill_with_nan(ADCP_FILLVALUE, e)

    # scale velocity to m/s
    e_scl = e / 1000.  # mm/s -> m/s

    # return the scaled Error Velocity Profile
    return e_scl


# Compute the VELTURB_L1 data products for the VADCP instrument deployed by RSN.
def vadcp_beam_eastward(b1, b2, b3, b4, h, p, r, vf, lat, lon, z, dt):
    """
    Description:

        Wrapper function to compute the Eastward Velocity Profile (VELTURB-VLE)
        from beam coordinate transformed velocity profiles as defined in the
        Data Product Specification for Turbulent Velocity Profile and Echo Intensity -
        DCN 1341-00760.

    Implemented by:

        2014-06-25: Christopher Wingard. Initial code, based on existing ADCP
        2015-06-10: Russell Desiderio.
                    (a) moved the conditioning of input beam velocities to adcp_beam2inst.
                    (b) moved the conditioning of compass readings to adcp_inst2earth.
                    (c) removed the depth dependence from the magnetic declination.

    Usage:

        uu_cor = vadcp_beam_eastward(b1, b2, b3, b4, h, p, r, vf, lat, lon, z, dt)

            where

        uu_cor = east velocity profiles in Earth coordinates corrected for the
                  magnetic declination (VELTURB-VLE_L1) [m s-1]

        b1 = "beam 1" velocity profiles in beam coordinates (VELTURB-B1_L0) [mm s-1]
        b2 = "beam 2" velocity profiles in beam coordinates (VELTURB-B2_L0) [mm s-1]
        b3 = "beam 3" velocity profiles in beam coordinates (VELTURB-B3_L0) [mm s-1]
        b4 = "beam 4" velocity profiles in beam coordinates (VELTURB-B4_L0) [mm s-1]
        h = instrument's uncorrected magnetic heading [cdegrees]
        p = instrument pitch [cdegrees]
        r = instrument roll [cdegrees]
        vf = instrument's vertical orientation (0 = downward looking and
            1 = upward looking)
        lat = instrument's deployment latitude [decimal degrees]
        lon = instrument's deployment longitude [decimal degrees]
        z = instrument's pressure sensor reading (depth) [daPa]
        dt = sample date and time value [seconds since 1900-01-01]
    """
    # force shapes of inputs to arrays of the correct dimensions
    #z = np.atleast_1d(z) / 1000.  # scale daPa depth input to dbar
    #z = z * 1.019716  # use a simple approximation to calculate depth in m
    lat = np.atleast_1d(lat)
    lon = np.atleast_1d(lon)
    dt = np.atleast_1d(dt)

    # compute the beam to instrument transform
    u, v, w, _ = adcp_beam2ins(b1, b2, b3, b4)

    # compute the instrument to earth beam transform
    uu, vv, _ = adcp_ins2earth(u, v, w, h, p, r, vf)

    # compute the magnetic variation, and ...
    theta = magnetic_declination(lat, lon, dt)

    # ... correct for it
    uu_cor, _ = magnetic_correction(theta, uu, vv)

    # scale velocity to m/s
    uu_cor = uu_cor / 1000.  # mm/s -> m/s

    # return the Eastward Velocity Profile
    return uu_cor


def vadcp_beam_northward(b1, b2, b3, b4, h, p, r, vf, lat, lon, z, dt):
    """
    Description:

        Wrapper function to compute the Northward Velocity Profile
        (VELTURB-VLN) from beam coordinate transformed velocity profiles as
        defined in the Data Product Specification for Turbulent Velocity
        Profile and Echo Intensity - DCN 1341-00760.

    Implemented by:

        2014-06-25: Christopher Wingard. Initial code, based on existing ADCP
        2015-06-10: Russell Desiderio.
                    (a) moved the conditioning of input beam velocities to adcp_beam2inst.
                    (b) moved the conditioning of compass readings to adcp_inst2earth.
                    (c) removed the depth dependence from the magnetic declination.

    Usage:

        vv_cor = vadcp_beam_northward(b1, b2, b3, b4, h, p, r, vf, lat, lon, z, dt)

            where

        vv_cor = north velocity profiles in Earth coordinates corrected for the
                  magnetic declination (VELTURB-VLN_L1) [m s-1]

        b1 = "beam 1" velocity profiles in beam coordinates (VELTURB-B1_L0) [mm s-1]
        b2 = "beam 2" velocity profiles in beam coordinates (VELTURB-B2_L0) [mm s-1]
        b3 = "beam 3" velocity profiles in beam coordinates (VELTURB-B3_L0) [mm s-1]
        b4 = "beam 4" velocity profiles in beam coordinates (VELTURB-B4_L0) [mm s-1]
        h = instrument's uncorrected magnetic heading [cdegrees]
        p = instrument pitch [cdegrees]
        r = instrument roll [cdegrees]
        vf = instrument's vertical orientation (0 = downward looking and
            1 = upward looking)
        lat = instrument's deployment latitude [decimal degrees]
        lon = instrument's deployment longitude [decimal degrees]
        z = instrument's pressure sensor reading (depth) [dm]
        dt = sample date and time value [seconds since 1900-01-01]
    """
    # force shapes of inputs to arrays of the correct dimensions
    #z = np.atleast_1d(z) / 1000.  # scale daPa depth input to dbar
    #z = z * 1.019716  # use a simple approximation to calculate depth in m
    lat = np.atleast_1d(lat)
    lon = np.atleast_1d(lon)
    dt = np.atleast_1d(dt)

    # compute the beam to instrument transform
    u, v, w, _ = adcp_beam2ins(b1, b2, b3, b4)

    # compute the instrument to earth beam transform
    uu, vv, _ = adcp_ins2earth(u, v, w, h, p, r, vf)

    # compute the magnetic variation, and ...
    theta = magnetic_declination(lat, lon, dt)

    # ... corect for it
    _, vv_cor = magnetic_correction(theta, uu, vv)

    # scale velocity to m/s
    vv_cor = vv_cor / 1000.  # mm/s -> m/s

    # return the Northward Velocity Profile
    return vv_cor


def vadcp_beam_vertical_est(b1, b2, b3, b4, h, p, r, vf):
    """
    Description:

        Wrapper function to compute the "estimated" Upward Velocity Profile
        (VELTURB-VLU-4BM) from the beam coordinate transformed velocity profiles as
        defined in the Data Product Specification for Turbulent Velocity
        Profile and Echo Intensity - DCN 1341-00760. This provides the
        traditional estimate of the vertical velocity component from a 4 beam
        solution, where each beam is facing outward at an angle (20 degrees)
        relative to the vertical.

    Implemented by:

        2014-06-25: Christopher Wingard. Initial code, based on existing ADCP
        2015-06-10: Russell Desiderio.
                    (a) moved the conditioning of input beam velocities to adcp_beam2inst.
                    (b) moved the conditioning of compass readings to adcp_inst2earth.
        2015-06-22: Russell Desiderio. Renamed this data product.

    Usage:

        ww_est = vadcp_beam_vertical_est(b1, b2, b3, b4, h, p, r, vf)

            where

        ww_est = estimated vertical velocity profiles in Earth coordinates
                 (VELTURB-VLU-4BM_L1) [m s-1]

        b1 = "beam 1" velocity profiles in beam coordinates (VELTURB-B1_L0) [mm s-1]
        b2 = "beam 2" velocity profiles in beam coordinates (VELTURB-B2_L0) [mm s-1]
        b3 = "beam 3" velocity profiles in beam coordinates (VELTURB-B3_L0) [mm s-1]
        b4 = "beam 4" velocity profiles in beam coordinates (VELTURB-B4_L0) [mm s-1]
        h = instrument's uncorrected magnetic heading [cdegrees]
        p = instrument pitch [cdegrees]
        r = instrument roll [cdegrees]
        vf = instrument's vertical orientation (0 = downward looking and
            1 = upward looking)
    """
    # compute the beam to instrument transform
    u, v, w, _ = adcp_beam2ins(b1, b2, b3, b4)

    # compute the instrument to earth beam transform
    _, _, ww = adcp_ins2earth(u, v, w, h, p, r, vf)

    # scale upward velocity to m/s
    ww = ww / 1000.  # mm/s -> m/s

    # return the estimated Upward Velocity Profile
    return ww


def vadcp_beam_vertical_true(b1, b2, b3, b4, b5, h, p, r, vf):
    """
    Description:

        Wrapper function to compute the "true" Upward Velocity Profile
        (VELTURB-VLU-5BM) from the beam coordinate transformed velocity profiles as
        defined in the Data Product Specification for Turbulent Velocity
        Profile and Echo Intensity - DCN 1341-00760. This is assumed to provide
        a better estimate of the true vertical velocity component, since beam 5
        is pointing directly up.

    Implemented by:

        2014-06-25: Christopher Wingard. Initial code, based on existing ADCP
        2015-06-10: Russell Desiderio.
                    (a) moved the conditioning of input beam velocities to adcp_beam2inst.
                    (b) moved the conditioning of compass readings to adcp_inst2earth.
        2015-06-22: Russell Desiderio. Renamed this data product.
        2015-06-25: Russell Desiderio. Incorporated b5 int fillvalue -> Nan.

    Usage:

        ww_true = vadcp_beam_vertical_true(b1, b2, b3, b4, b5, h, p, r, vf)

            where

        ww_true = true vertical velocity profiles in Earth coordinates
                  (VELTURB-VLU-5BM_L1) [m s-1]

        b1 = "beam 1" velocity profiles in beam coordinates (VELTURB-B1_L0) [mm s-1]
        b2 = "beam 2" velocity profiles in beam coordinates (VELTURB-B2_L0) [mm s-1]
        b3 = "beam 3" velocity profiles in beam coordinates (VELTURB-B3_L0) [mm s-1]
        b4 = "beam 4" velocity profiles in beam coordinates (VELTURB-B4_L0) [mm s-1]
        b5 = "beam 5" velocity profiles in beam coordinates (VELTURB-B5_L0) [mm s-1]
        h = instrument's uncorrected magnetic heading [cdegrees]
        p = instrument pitch [cdegrees]
        r = instrument roll [cdegrees]
        vf = instrument's vertical orientation (0 = downward looking and
            1 = upward looking)
    """
    # compute the beam to instrument transform
    # fill values in the 4 beams are checked for inside adcp_beam2ins
    u, v, _, _ = adcp_beam2ins(b1, b2, b3, b4)

    # check b5 for the presence of fill values
    b5 = replace_fill_with_nan(ADCP_FILLVALUE, b5)

    # compute the instrument to earth beam transform
    # fill values in the adcp orientation parameters are checked for inside adcp_ins2earth
    _, _, ww = adcp_ins2earth(u, v, b5, h, p, r, vf)

    # scale upward velocity to m/s
    ww = ww / 1000.  # mm/s -> m/s

    # return the true Upward Velocity Profile
    return ww


def vadcp_beam_error(b1, b2, b3, b4):
    """
    Description:

        Wrapper function to compute the Error Velocity Profile (VELTURB-ERR)
        from the beam coordinate transformed velocity profiles as defined in
        the Data Product Specification for Turbulent Velocity Profile and Echo
        Intensity - DCN 1341-00760.

    Implemented by:

        2014-06-25: Christopher Wingard. Initial code, based on existing ADCP
        2015-06-10: Russell Desiderio.
                    Moved the conditioning of input beam velocities to adcp_beam2inst.

    Usage:

        e = vadcp_beam_northward(b1, b2, b3, b4)

            where

        e = error velocity profiles (VELTURB-ERR_L1) [m s-1]

        b1 = "beam 1" velocity profiles in beam coordinates (VELTURB-B1_L0) [mm s-1]
        b2 = "beam 2" velocity profiles in beam coordinates (VELTURB-B2_L0) [mm s-1]
        b3 = "beam 3" velocity profiles in beam coordinates (VELTURB-B3_L0) [mm s-1]
        b4 = "beam 4" velocity profiles in beam coordinates (VELTURB-B4_L0) [mm s-1]
    """
    # compute the beam to instrument transform
    _, _, _, e = adcp_beam2ins(b1, b2, b3, b4)

    # scale error velocity to m/s
    e = e / 1000.   # mm/s

    # return the Error Velocity Profile
    return e


# Calculates ECHOINT_L1 for all tRDI ADCPs
def adcp_backscatter(raw, sfactor):
    """
    Description:

        Converts the echo intensity data from counts to dB using a factory
        specified scale factor (nominally 0.45 dB/count for the Workhorse
        family of ADCPs and 0.61 dB/count for the ExplorerDVL family). As
        defined in the Data Product Specification for Velocity Profile and Echo
        Intensity - DCN 1341-00750.

    Implemented by:

        2014-04-21: Christopher Wingard. Initial code.
        2015-06-25: Russell Desiderio. Incorporated int fillvalue -> Nan.

    Usage:

        dB = adcp_backscatter(raw, sfactor)

            where

        dB = Relative Echo Intensity (ECHOINT_L1) [dB]

        raw = raw echo intensity (ECHOINT_L0) [count]
        sfactor = factory supplied scale factor, instrument and beam specific [dB/count]

    Notes:

        The ADCP outputs the raw echo intensity as a 1-byte integer, so the ADCP_FILLVALUE
        cannot apply (requires 2 bytes).

    References:

        OOI (2012). Data Product Specification for Velocity Profile and Echo
            Intensity. Document Control Number 1341-00750.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00050_Data_Product_SPEC_VELPROF_OOI.pdf)
    """
    if np.isscalar(sfactor) is False:
        sfactor = sfactor.reshape(sfactor.shape[0], 1)

    # check raw for the presence of system fill values
    raw = replace_fill_with_nan(None, raw)

    dB = raw * sfactor
    return dB


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
        2015-06-24: Russell Desiderio. Incorporated int fillvalue -> Nan.

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
    b1 = np.atleast_2d(b1)
    b2 = np.atleast_2d(b2)
    b3 = np.atleast_2d(b3)
    b4 = np.atleast_2d(b4)

    b1, b2, b3, b4 = replace_fill_with_nan(ADCP_FILLVALUE, b1, b2, b3, b4)

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
        2015-06-24: Russell Desiderio. Changed implementation of 'vertical' in the roll
                    calculation so that if these values are equal to the CI fill value
                    (-999999999), when these fill values are replaced with nans, the nans
                    will propagate through to the data product output.
        2015-06-24: Russell Desiderio. Incorporated int fillvalue -> Nan.

    Usage:

        uu, vu, ww = adcp_ins2earth(u, v, w, heading, pitch, roll, vertical)

            where

        uu = "east" velocity profiles in earth coordinates [mm s-1]
        vv = "north" velocity profiles in earth coordinates [mm s-1]
        ww = "vertical" velocity profiles in earth coordinates [mm s-1]

        u = east velocity profiles in instrument coordinates [mm s-1]
        v = north velocity profiles in instrument coordinates [mm s-1]
        w = vertical velocity profiles in instrument coordinates [mm s-1]
        heading = instrument's uncorrected magnetic heading [centidegrees]
        pitch = instrument pitch [centidegrees]
        roll = instrument roll [centidegrees]
        vertical = instrument's vertical orientation (0 = downward looking and
            1 = upward looking)

    References:

        OOI (2012). Data Product Specification for Velocity Profile and Echo
            Intensity. Document Control Number 1341-00750.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00050_Data_Product_SPEC_VELPROF_OOI.pdf)
    """
    ### the input beam data for adcp_ins2earth are always called using the output
    ### of adcp_beam2ins, so the following lines are not needed.
    # insure we are dealing with array inputs
    #u = np.atleast_2d(u)
    #v = np.atleast_2d(v)
    #w = np.atleast_2d(w)

    # check for CI fill values before changing units.
    # this function 'conditions' (np.atleast_1d) its inputs.
    # TRDI does not apply its ADCP fill/bad value sentinels to compass data.
    heading, pitch, roll, vertical = replace_fill_with_nan(None, heading, pitch, roll, vertical)

    # change units from centidegrees to degrees
    heading = heading / 100.0
    pitch = pitch / 100.0
    roll = roll / 100.0

    # better way to calculate roll from the vertical orientation toggle;
    # this will propagate R as nans if the vertical variable is missing from the data.
    R = roll + vertical * 180.0

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

    # pack the coordinates (velocities) to be transformed into the appropriate
    # slices.
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
        2015-04-10: Russell Desiderio. Corrected a typo:
                    uv = np.atleast_2d(u)  ->  u = np.atleast_2d(u)

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
    u = np.atleast_2d(u)
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


# Calculates bin depths tRDI ADCPs configured to output data using the PD0 and PD12 formats
def adcp_bin_depths(dist_first_bin, bin_size, num_bins, pressure, adcp_orientation, latitude):
    """
    Description:

        Calculates the center bin depths for PD0 and PD12 ADCP data. As defined
        in the Data Product Specification for Velocity Profile and Echo
        Intensity - DCN 1341-00750.

    Implemented by:

        2015-01-29: Craig Risien. Initial code.
        2015-06-26: Russell Desiderio. Fixed the handling of the pressure variables.
                                       Time-vectorized the code by finessing the conditional.
        2015-06-30: Russell Desiderio. Incorporated int fillvalue -> Nan.


    Usage:

        bin_depths = adcp_bin_depths(dist_first_bin, bin_size, num_bins, pressure,
                                    adcp_orientation, latitude)

            where

        bin_depths =  [meters]

        dist_first_bin = distance to the first ADCP bin [centimeters]
        bin_size = depth of each ADCP bin [centimeters]
        num_bins = number of ADCP bins [unitless]
        pressure = pressure at the sensor head [dPa]
        adcp_orientation = 1=upward looking or 0=downward looking [unitless]
        latitude = latitude of the instrument [degrees]

    References:

        OOI (2012). Data Product Specification for Velocity Profile and Echo
            Intensity. Document Control Number 1341-00750.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00050_Data_Product_SPEC_VELPROF_OOI.pdf)
    """
    # check for CI fill values.
    dist_first_bin, bin_size, num_bins, pressure, adcp_orientation = replace_fill_with_nan(
        None, dist_first_bin, bin_size, num_bins, pressure, adcp_orientation)

    # note, there is a CI problem not yet addressed if the time-vectorized values
    # in num_bins are not all the same!! For now, assume they are all the same:
    num_bins_constant = num_bins[0]
    # make bin_numbers a row vector
    bin_numbers = np.array([np.arange(num_bins_constant)])

    # Convert from cm to meters
    dist_first_bin = dist_first_bin / 100.0
    bin_size = bin_size / 100.0

    # Convert pressure from decaPascal to decibar
    pressure_dbar = pressure / 1000.0

    # Calculate sensor depth using TEOS-10 toolbox z_from_p function
    # note change of sign to make the sensor_depth variable positive
    sensor_depth = -z_from_p(pressure_dbar, latitude)

    # For the PD0 convention:
    #     adcp_orientation = 0 is downward looking, bindepths are added to sensor depth
    #                      = 1 is upward looking, bindepths are subtracted from sensor depth
    z_sign = 1.0 - 2.0 * adcp_orientation

    # to broadcast the vertical time dimension correctly with the horizontal bin_numbers dimension,
    # make all the 1D time arrays into column vectors to be processed with the bin_numbers row vector.
    sensor_depth = sensor_depth.reshape(-1, 1)
    z_sign = z_sign.reshape(-1, 1)
    dist_first_bin = dist_first_bin.reshape(-1, 1)
    bin_size = bin_size.reshape(-1, 1)

    # Calculate bin depths
    bin_depths = sensor_depth + z_sign * (dist_first_bin + bin_size * bin_numbers)

    return bin_depths


def z_from_p(p, lat, geo_strf_dyn_height=0, sea_surface_geopotential=0):
    """Calculates height from sea pressure using the computationally-efficient
    75-term expression for density in terms of SA, CT and p (Roquet et al.,
    2015). Dynamic height anomaly, geo_strf_dyn_height, if provided, must be
    computed with its p_ref=0 (the surface). Also if provided, sea_surface_geopotental
    is the geopotential at zero sea pressure.

    Calls a function which calculates enthalpy assuming standard ocean salinity
    and 0 degrees celsius.

    Parameters
    ----------
    p : pressure [dbar]
    lat : latitude in decimal degrees north [-90..+90]
    geo_strf_dyn_height : dynamic height anomaly [m^2/s^2]
    sea_surface_geopotential : geopotential at zero sea pressure  [ m^2/s^2 ]

    Returns
    -------
    z : TEOS-10 height [m] : height is returned as a negative number; its
                             absolute value is the depth below the sea surface.

    #################################################################
    #  Check values from TEOS-10 version 3.05 (matlab code):        #
    #  from http://www.teos-10.org/pubs/gsw/html/gsw_z_from_p.html  #
    #################################################################

    p = [10, 50, 125, 250, 600, 1000]
    lat = 4

    z_from_p(p, lat) =
    [  -9.9445834469453,  -49.7180897012550, -124.2726219409978,
     -248.4700576548589, -595.8253480356214, -992.0919060719987]

    Notes
    -----
    At sea level z = 0, and since z (HEIGHT) is defined to be positive upwards,
    it follows that while z is positive in the atmosphere, it is NEGATIVE in
    the ocean.

    References
    ----------
    IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
     seawater - 2010: Calculation and use of thermodynamic properties.
     Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
     UNESCO (English), 196 pp.  Available from the TEOS-10 web site.

    McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003: 
     Accurate and computationally efficient algorithms for potential 
     temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
     pp. 730-741.

    Moritz, 2000: Goedetic reference system 1980. J. Geodesy, 74, 128-133.

    Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
     polynomial expressions for the density and specifc volume of seawater
     using the TEOS-10 standard. Ocean Modelling.

    Saunders, P. M., 1981: Practical conversion of pressure to depth.
     Journal of Physical Oceanography, 11, 573-574.

    IMPLEMENTATION NOTES:

        Russell Desiderio. 2015_07_01
            versions 3.04 and 3.05 of the main function z_from_p are identical.

            z_from_p calls the subroutine enthalpy_SSO_0_p; this subroutine
            has been updated from ver 3.04 to 3.05.

            the check values above for z_from_p have been updated to incorporate
            this change using enthalpy_SSO_0_p ver 3.05.

    """

    X = np.sin(np.deg2rad(lat))
    sin2 = X ** 2
    B = 9.780327 * (1.0 + (5.2792e-3 + (2.32e-5 * sin2)) * sin2)
    gamma = 2.26e-07
    A = -0.5 * gamma * B
    C = enthalpy_SSO_0_p(p) - geo_strf_dyn_height

    return -2 * C / (B + np.sqrt(B ** 2 - 4 * A * C))


def enthalpy_SSO_0_p(p):
    """
    This documentation and code is copy\pasted from the matlab coding of this function.

    %==========================================================================
    %  This function calculates enthalpy at the Standard Ocean Salinity, SSO,
    %  and at a Conservative Temperature of zero degrees C, as a function of
    %  pressure, p, in dbar, using a streamlined version of the 76-term
    %  computationally-efficient expression for specific volume, that is, a
    %  streamlined version of the code "gsw_enthalpy(SA,CT,p)".
    %
    % VERSION NUMBER: 3.05 (27th January 2015)
    %
    % REFERENCES:
    %  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
    %   polynomial expressions for the density and specifc volume of seawater
    %   using the TEOS-10 standard. Ocean Modelling.
    %
    %==========================================================================

    IMPLEMENTATION NOTES:

        Russell Desiderio. 2015_07_01. this subroutine has been updated
                                       from ver 3.04 to 3.05.
    """
    z = p * 1e-4

    h006 = -2.1078768810e-9
    h007 = 2.8019291329e-10

    dynamic_enthalpy_SSO_0_p = z * (9.726613854843870e-4 + z * (-2.252956605630465e-5 + z * (
        2.376909655387404e-6 + z * (-1.664294869986011e-7 + z * (
            -5.988108894465758e-9 + z * (h006 + h007 * z))))))

    enthalpy_SSO_0 = dynamic_enthalpy_SSO_0_p * 1.e8  # Note. 1e8 = db2Pa*1e4

    return enthalpy_SSO_0


# Calculates bin depths tRDI ADCPs configured to output data using the PD8 format
def adcp_bin_depths_pd8(dist_first_bin, bin_size, num_bins, sensor_depth, adcp_orientation):
    """
    Description:

        Calculates the center bin depths for PD8 ADCP data. As defined
        in the Data Product Specification for Velocity Profile and Echo
        Intensity - DCN 1341-00750.

    Implemented by:

        2015-01-30: Craig Risien. Initial code.
        2015-06-26: Russell Desiderio. Time-vectorized the code by finessing the conditionals.
        2015-06-30: Russell Desiderio. Incorporated int fillvalue -> Nan.

    Usage:

        bin_depths_pd8 = adcp_bin_depths(dist_first_bin, bin_size, num_bins, sensor_depth,
                                    adcp_orientation)

            where

        bin_depths_pd8 =  [meters]

        dist_first_bin = distance to the first ADCP bin [centimeters]
        bin_size = depth of each ADCP bin [centimeters]
        num_bins = number of ADCP bins [unitless]
        sensor_depth = estimated depth at the sensor head [meters]
        adcp_orientation = 1=upward looking or 0=downward looking [unitless]

    Notes:

        The PD8 output format is a very sparse format. Other than num_bins, it does *not* record
        any of the other input variables required by this DPA. Those must somehow be supplied "by
        hand".

    """
    # check for CI fill values.
    #
    # Note that these input parameters will not come from an IDD driver (except for possibly
    # (num_bins) because the PD8 output format does not output them. Therefore, I don't know
    # if they will be of type integer or not. However, ndarrays composed of float types are
    # passed through the check-code unchanged, so run the inputs through in case they are of
    # type int and in case -999999999 fill values are somehow present.
    dist_first_bin, bin_size, num_bins, sensor_depth, adcp_orientation = replace_fill_with_nan(
        None, dist_first_bin, bin_size, num_bins, sensor_depth, adcp_orientation)

    # note, there is a CI problem not yet addressed if the time-vectorized values
    # in num_bins are not all the same!! For now, assume they are all the same:
    num_bins_constant = num_bins[0]
    # make bin_numbers a row vector
    bin_numbers = np.array([np.arange(num_bins_constant)])

    # Convert from cm to meters
    # the input variables are type integer, so divide by a real number
    # to avoid truncation errors.
    dist_first_bin = dist_first_bin / 100.0
    bin_size = bin_size / 100.0

    # make sure sensor depth is positive
    sensor_depth = np.fabs(sensor_depth)

    # Following the PD0 convention where
    #     adcp_orientation = 0 is downward looking, bindepths are added to sensor depth
    #                      = 1 is upward looking, bindepths are subtracted from sensor depth
    z_sign = 1.0 - 2.0 * adcp_orientation

    # to broadcast the vertical time dimension correctly with the horizontal bin_numbers dimension,
    # make all the 1D time arrays into column vectors to be processed with the bin_numbers row vector.
    sensor_depth = sensor_depth.reshape(-1, 1)
    z_sign = z_sign.reshape(-1, 1)
    dist_first_bin = dist_first_bin.reshape(-1, 1)
    bin_size = bin_size.reshape(-1, 1)

    # Calculate bin depths
    bin_depths_pd8 = sensor_depth + z_sign * (dist_first_bin + bin_size * bin_numbers)

    return bin_depths_pd8
