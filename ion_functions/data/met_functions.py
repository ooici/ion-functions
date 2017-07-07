#!/usr/bin/env python
"""
@package ion_functions.data.met_functions
@file ion_functions/data/met_functions.py
@author Russell Desiderio
@brief Module containing functions for the met family of instruments
"""

import numpy as np
import numexpr as ne
from pygsw import vectors as gsw

from ion_functions.data.generic_functions import magnetic_declination, magnetic_correction

### July 2015
""" use_velptmn_with_metbk """
# Until VELPT current meters are
#     (1) co-located with the METBK instrumentation to measure surface currents and
#     (2) without ferrous interference aliasing the VELPT compass measurements,
# VELPT current measurements should not be used to calculate relative wind speed (with
# respect to water), which is the fundamental windspeed variable used in the METBK
# calculations. Almost all of the METBK L2 data products require the relative wind speed
# as a calling argument.
#
# The DPA to calculate relative wind speed over water is currently set to return actual
# wind speed (as if the current velocities were measured to be 0) for all cases of input
# current velocity values (absent, present, nan).
#
# A subset of the Endurance moorings are the only ones that have VELPT instruments
# mounted to measure surface currents. However, their compass readings are inaccurate due
# to the mounted instruments' proximity to iron ballast in the mooring.
#
# It is anticipated that these moorings will be modified to eliminate this magnetic
# interference. To use the velptmn measurements for METBK calculations on these moorings,
# a 5th variable, use_velptmn_with_metbk, has been added to the argument list of the
# function met_relwind_speed. Implementation of the use_velptmn_with_metbk variable will
# require that it be treated as a "platform/instrument instance specific metadata parameter
# changeable in time". It has been coded to accept time-vectorized input.
#
### Further documentation is contained in the Notes to function met_relwind_speed in this module.

"""
    METBK SENSOR HEIGHTS

    Note that the sensor heights may depend on the type of mooring:
"""
# these 4 sensor height variables were time-vectorized in the July 2015 revision.
#     zwindsp = height of the wind measurement [m]
#     ztmpair = height of air temperature measurement [m]
#     zhumair = height of air humidity measurement [m]
#     ztmpwat = depth of bulk sea surface water measurements [m]

#     zvelptm = depth of surface current measurement [m]:
#         this parameter is specified as metadata in the DPS;
#         however, it is not used in the code.

#     zinvpbl = planetary boundary layer/inversion height: this is
#               set to a default value of 600m in the code, as is
#               used in the DPS. this variable was written to accept
#               time-vectorized input in the initial METBK code.

"""
    Set algorithm switches used in METBK bulk flux calculations
"""
# The jcool and jwarm switches should be set to 1, always!
JCOOLFL = 1      # 1=do coolskin calc
JWARMFL = 1      # 1=do warmlayer calc
#JWAVEFL         # only the windspeed parametrization of the charnok
                 # variable is coded; therefore this switch is not used.

"""
    LISTING OF SUBROUTINES BY ORDER IN THIS MODULE
        Grouped by sections; alphabetical within each section.

        The functions which directly calculate data products, both formal and meta,
        are listed here as the actual data product name in upper case, rather than
        by the name of the function; the functions themselves are named as
        "met_prdname" except as noted.

        All other functions are listed by function name.
#...................................................................................
    Functions to compute the L1 BULKMET (METBK) data products:
    these do not require the 'warmlayer/coolskin' iteration algorithm:
        BARPRES
        WINDAVG-VLE  (function name: met_windavg_mag_corr_east)
        WINDAVG-VLN  (function name: met_windavg_mag_corr_north)
    These products are calculated at the native temporal resolution of the
    instrument suite (roughly each minute).
#...................................................................................
#...................................................................................
    Functions to compute the (simpler) metadata products that do not require
    the 'warmlayer/coolskin' iteration algorithms:
        CURRENT_DIR
        CURRENT_SPD
        RELWIND_DIR-AUX
        RELWIND_SPD-AUX
        TIMEFLX-AUX
    These products are calculated at the native temporal resolution of the
    instrument suite (roughly each minute), EXCEPT for TIMEFLX-AUX (hourly).
#...................................................................................
#...................................................................................
    Functions to compute the (simpler) L2 METBK data products that do not require
    the 'warmlayer/coolskin' iteration algorithm:
        NETSIRR (this may operationally be an L1 product)
        NETSIRR_HOURLY (on an hourly time base)
        RAINRTE
        SALSURF
        SPECHUM
    These products are calculated at the native temporal resolution of the
    instrument suite (roughly each minute) except as noted.
#...................................................................................
#...................................................................................
    Functions to compute the L2 METBK data products that do require
    the 'warmlayer/coolskin' iteration algorithm:
        BUOYFLS:  added DPA to match FDCHP, not in original DPS
        BUOYFLX:  added DPA to match FDCHP, not in original DPS
        FRSHFLX
        HEATFLX
        LATNFLX
        MOMMFLX
        NETLIRR
        RAINFLX
        SENSFLX
        SPHUM2M
        STABLTY:  metadata
        TEMPA2M
        TEMPSKN:  metadata
        WIND10M
    These products are calculated on hourly averages.
#...................................................................................
#...................................................................................
    Functions to compute the L2 METBK data products (that do require
    the 'warmlayer/coolskin' iteration algorithm) at the native temporal
    resolution of the instrument suite (roughly per minute).

    These functions were not specified in the DPS.
        HEATFLX_MINUTE
        LATNFLX_MINUTE
        NETLIRR_MINUTE
        SENSFLX_MINUTE
#...................................................................................
#...................................................................................
    Simple subroutines used in the routines in the sections above.
        air_density
        airtemp_at_refheight
        calc_rain_rate
        gravity
        latent_heat_vaporization_pure_water
        net_longwave_up
        psit_26
        psiu_26
        rain_heat_flux
        sea_spechum
        spechum_at_refheight
        water_thermal_expansion
        windspeed_at_refheight
#...................................................................................
#...................................................................................
    seasurface_skintemp_correct  (wrapper; calls warmlayer and coare35vn)
#...................................................................................
#...................................................................................
    warmlayer ('warmlayer' toga-coare routine)
#...................................................................................
#...................................................................................
    coare35vn (bulk calculation + 'coolskin' toga-coare routines; plus subroutines)
        charnock_wind
        coolskin_parameters
        effective_relwind
        obukhov_for_init
        obukhov_length_scale
        roughness_lengths
        roughness_lengths_for_init
        scaling_parameters
#...................................................................................
#...................................................................................
    Data conditioning and averaging routines
        vet_velptmn_data
        condition_data
        make_hourly_data
        warmlayer_time_keys
#...................................................................................

"""

####################################################################################
####################################################################################
####################################################################################

"""
#...................................................................................
#...................................................................................
    Functions to compute the L1 BULKMET (METBK) data products:
    do not require the 'warmlayer/coolskin' iteration algorithm:

        BARPRES
        WINDAVG-VLE  (met_windavg_mag_corr_east)
        WINDAVG-VLN  (met_windavg_mag_corr_north)

    These products are calculated at the native temporal resolution of the
    instrument suite (roughly each minute).
#...................................................................................
#...................................................................................

"""


def met_barpres(mbar):
    """
    Description:

        OOI Level 1 Air Pressure core data product, which is calculated 
        by scaling the measured barometer pressure from mbar to Pascals.

    Implemented by:

        2014-06-25: Christopher Wingard. Initial code.

    Usage:

        Pa = met_barpres(mbar)

            where

        Pa = Air pressure (BARPRES_L1) [Pa]
        mbar = Barometer pressure (BARPRES_L0) [mbar]

    References:

        OOI (2012). Data Product Specification for L1 Bulk Meterological Data
            Products. Document Control Number 1341-00360.
            https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00360_Data_Product_SPEC_BULKMET_OOI.pdf)
    """
    Pa = mbar * 100.
    return Pa


def met_windavg_mag_corr_east(uu, vv, lat, lon, timestamp, zwindsp=0.0):
    """
    Description:

        Calculates WINDAVG-VLE_L1, the OOI Level 1 core data product for windspeed in the
        true eastward direction, for the METBK instrument by correcting for magnetic declination.

    Implemented by:

        2014-06-25: Christopher Wingard. Initial code.
        2014-08-26: Russell Desiderio. Added documentation.

    Usage:

        uu_cor = windavg_mag_corr_east(uu, vv, lat, lon, timestamp[, zwindsp])

            where

        uu_cor = WINDAVG-VLE_L1 [m/s], METBK eastward windspeed corrected for magnetic declination.
        uu     = WINDAVG-VLE_L0 [m/s], METBK eastward windspeed, uncorrected.
        vv     = WINDAVG-VLN_L0 [m/s], METBK northward windspeed, uncorrected.
        lat    = instrument's deployment latitude [decimal degrees]
        lon    = instrument's deployment longitude [decimal degrees]
        timestamp = sample date and time value [seconds since 1900-01-01]
        zwindsp  = [optional] height of windspeed sensor above sealevel [m].

    References:

        OOI (2012). Data Product Specification for L1 Bulk Meterological Data
            Products. Document Control Number 1341-00360.
            https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00360_Data_Product_SPEC_BULKMET_OOI.pdf)

        The magnetic_declination function and the magnetic_correction function are
        called from the ion_functions.data.generic_function module.
            magnetic_declination calculates the value(s) for the declination;
            magnetic_correction rotates the velocity vectors from the magnetic
                compass headings to true compass headings.
    """
    # calculate the magnetic declination using the WMM model
    zflag = 1  # denotes that z is a height above sealevel.
    mag_dec = magnetic_declination(lat, lon, timestamp, zwindsp, zflag)

    # rotate the vectors from the magnetic to the true compass frame
    magvar = np.vectorize(magnetic_correction)
    uu_cor, vv_cor = magvar(mag_dec, uu, vv)

    return uu_cor


def met_windavg_mag_corr_north(uu, vv, lat, lon, timestamp, zwindsp=0.0):
    """
    Description:

        Calculates WINDAVG-VLN_L1, the OOI Level 1 core data product for windspeed in the
        true northward direction, for the METBK instrument by correcting for magnetic declination.

    Implemented by:

        2014-06-25: Christopher Wingard. Initial code.
        2014-08-26: Russell Desiderio. Added documentation.

    Usage:

        vv_cor = windavg_mag_corr_north(uu, vv, lat, lon, timestamp[, zwindsp])

            where

        vv_cor = WINDAVG-VLN_L1 [m/s], METBK northward windspeed corrected for magnetic declination.
        uu     = WINDAVG-VLE_L0 [m/s], METBK eastward windspeed, uncorrected.
        vv     = WINDAVG-VLN_L0 [m/s], METBK northward windspeed, uncorrected.
        lat    = instrument's deployment latitude [decimal degrees]
        lon    = instrument's deployment longitude [decimal degrees]
        timestamp = sample date and time value [seconds since 1900-01-01]
        zwindsp = [optional] height of windspeed sensor above sealevel [m].

    References:

        OOI (2012). Data Product Specification for L1 Bulk Meterological Data
            Products. Document Control Number 1341-00360.
            https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00360_Data_Product_SPEC_BULKMET_OOI.pdf)

        The magnetic_declination function and the magnetic_correction function are
        called from the ion_functions.data.generic_function module.
            magnetic_declination calculates the value(s) for the declination;
            magnetic_correction rotates the velocity vectors from the magnetic
                compass headings to true compass headings.
    """
    # calculate the magnetic declination using the WMM model
    zflag = 1  # denotes that z is a height above sealevel.
    mag_dec = magnetic_declination(lat, lon, timestamp, zwindsp, zflag)

    # rotate the vectors from the magnetic to the true compass frame
    magvar = np.vectorize(magnetic_correction)
    uu_cor, vv_cor = magvar(mag_dec, uu, vv)

    return vv_cor


"""
#...................................................................................
#...................................................................................
    Functions to compute the (simpler) metadata products that do not require
    the warmlayer/coolskin iteration algorithms:
        CURRENT_DIR
        CURRENT_SPD
        RELWIND_DIR-AUX
        RELWIND_SPD-AUX
        TIMEFLX-AUX

    These products are calculated at the native temporal resolution of the
    instrument suite (roughly each minute), EXCEPT for TIMEFLX-AUX (hourly).
#...................................................................................
#...................................................................................

"""


def met_current_direction(vle_water, vln_water, use_velptmn_with_metbk=0):
    """
    Description:

        Calculates the direction of the surface current using the eastward and northward
        velocity components from the VELPT mounted on the surface buoy.

    Implemented by:

        2014-08-27: Russell Desiderio. Initial Code
        2015-07-10: Russell Desiderio. Added data quality flags (use_velptmn_with_metbk)
                    to argument list. See Notes to the function met_relwind_speed.

    Usage:

        current_dir = met_current_direction(vle_water, vln_water[, use_velptmn_with_metbk])

            where

        current_dir = direction of the surface current (CURRENT_DIR) [0 360) degrees
        vle_water = eastward surface current (VELPTMN-VLE_L1) [m/s]
        vln_water = northward surface current (VELPTMN-VLN_L1) [m/s]
        use_velptmn_with_metbk = time-vectorized data quality flag:
                                 0 -> bad  velptmn current data
                                 1 -> good velptmn current data

    Notes:

        The auxiliary data product calculated by this function is not used in any of
        the METBK functions in this module. It is calculated because it is of scientific
        interest. See also the Notes to the function met_relwind_speed.

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
    # replace aliased current values with nans.
    vle_water, vln_water = vet_velptmn_data(vle_water, vln_water, use_velptmn_with_metbk)

    # use arctan2, which will properly handle arguments in all 4 cartesian quadrants,
    # and gives answers [-180 180] after conversion to degrees.
    cartesian_dir = np.degrees(np.arctan2(vln_water, vle_water))

    # the angle above is defined for a cartesian (x,y) coordinate system:
    # an angle of 0 points along the positive x-axis ("east" instead of
    # "north"), and increasingly positive angles indicate a ccw rotation
    # (instead of clockwise, which is the compass heading convention).

    # to convert to the compass convention, flip the sign and then add 90 degrees.
    # to change the range of values to [0 360), add 360 and 'mod' the result.
    current_dir = np.mod(450 - cartesian_dir, 360)

    return current_dir


def met_current_speed(vle_water, vln_water, use_velptmn_with_metbk=0):
    """
    Description:

        Estimate the magnitude of the surface current using the eastward and northward
        velocity components from the VELPT mounted on the surface buoy. This is the
        meta-data "product" CURRENT specified by the DPS referenced below (Section 4.3,
        step 4). This product is not used in the METBK code; rather, the magnitude of
        the vector difference of the wind and current vectors is the fundamental
        variable used in the METBK calculations (see RELWIND_SPD-AUX).

        Because the direction of the current will also be calculated so as to be
        made available, the CURRENT metadata product is sub-divided into:
        CURRENT_SPD (calculated by this code) and
        CURRENT_DIR.

    Implemented by:

        2014-06-26: Chris Wingard. Initial Code
        2014-08-27: Russell Desiderio. Added documentation, changed variable names.
        2015-07-10: Russell Desiderio. Added data quality flags (use_velptmn_with_metbk)
                    to argument list. See Notes to the function met_relwind_speed.

    Usage:

        current_spd = met_current_speed(vle_water, vln_water[, use_velptmn_with_metbk])

            where

        current_spd = magnitude (speed) of the surface current (CURRENT_SPD) [m/s]
        vle_water = eastward surface current (VELPTMN-VLE_L1) [m/s]
        vln_water = northward surface current (VELPTMN-VLN_L1) [m/s]
        use_velptmn_with_metbk = time-vectorized data quality flag:
                                 0 -> bad  velptmn current data
                                 1 -> good velptmn current data

    Notes:

        The auxiliary data product calculated by this function is not used in any of
        the METBK functions in this module. It is calculated because it is of scientific
        interest and because the DPS specified it. See also the Notes to the function
        met_relwind_speed.

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
    # replace aliased current values with nans.
    vle_water, vln_water = vet_velptmn_data(vle_water, vln_water, use_velptmn_with_metbk)

    current_spd = ne.evaluate("sqrt(vle_water**2 + vln_water**2)")
    return current_spd


def met_relwind_direction(vle_wind, vln_wind, vle_water=None, vln_water=None, use_velptmn_with_metbk=0):
    """
    Description:

        Calculates RELWIND_DIR-AUX, the direction of the vector difference of wind velocity
        (from METBK measurements) and surface current (from VELPT measurements).

        It is anticipated that the wind measurements will be roughly each minute and that the
        current measurements will be broadcast to that resolution.

    Implemented by:

        2014-08-26: Russell Desiderio. Initial Code.
        2015-07-10: Russell Desiderio. Set default calling water velocity values and implemented
                                       use_velptmn_with_metbk switch.

    Usage:

        u_dir = met_relwind_direction(vle_wind, vln_wind[, vle_water, vln_water[, use_velptmn_with_metbk]])

            where

        u_dir = direction of relative wind (RELWIND_DIR-AUX) [0 360) degrees
        vle_wind  = eastward wind speed (WINDAVG-VLE_L1) [m/s]
        vln_wind  = northward wind speed (WINDAVG-VLN_L1) [m/s]
        vle_water = eastward surface current (VELPTMN-VLE_L1) [m/s]
        vln_water = northward surface current (VELPTMN-VLN_L1) [m/s]
        use_velptmn_with_metbk = time-vectorized data quality flag:
                                 0 -> bad  velptmn current data
                                 1 -> good velptmn current data

    Notes:

        The auxiliary data product calculated by this function is not used in any of
        the METBK functions in this module. It is calculated because it is of scientific
        interest. See also the Notes to the function met_relwind_speed.

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
    # if this function is called without using surface current data, return nan
    if vle_water is None or vln_water is None:
        u_dir = vle_wind * np.nan
        return u_dir

    # replace aliased current values with nans.
    vle_water, vln_water = vet_velptmn_data(vle_water, vln_water, use_velptmn_with_metbk)

    # use arctan2, which will properly handle arguments in all 4 cartesian quadrants,
    # and gives answers [-180 180] after conversion to degrees.
    cartesian_dir = np.degrees(np.arctan2(vln_wind - vln_water, vle_wind - vle_water))

    # the angle above is defined for a cartesian (x,y) coordinate system:
    # an angle of 0 points along the positive x-axis ("east" instead of
    # "north"), and increasingly positive angles indicate a ccw rotation
    # (instead of clockwise, which is the compass heading convention).

    # to convert to the compass convention, flip the sign and then add 90 degrees.
    # to change the range of values to [0 360), add 360 and 'mod' the result.
    u_dir = np.mod(450 - cartesian_dir, 360)

    return u_dir


def met_relwind_speed(vle_wind, vln_wind, vle_water=None, vln_water=None, use_velptmn_with_metbk=0):
    """
    Description:

        Calculates RELWIND_SPD-AUX, the relative windspeed over water, calculated as the
        magnitude of the vector difference of surface current velocity (from VELPT
        measurements) subtracted from wind velocity (from METBK measurements).

        This is the fundamental windspeed variable used in the METBK toga-coare algorithms.

        It is anticipated that the wind measurements will be roughly each minute and that the
        current measurements will be broadcast to that resolution.

    Implemented by:

        2014-08-26: Russell Desiderio. Initial Code.
        2015-07-10: Russell Desiderio. Set default calling water velocity values.
                    Added the switch use_velptmn_with_metbk and code to vet surface current data.

    Usage:

        u_rel = met_relwind_speed(vle_wind, vln_wind[, vle_water, vln_water[, use_velptmn_with_metbk]])

            where

        u_rel = magnitude of windspeed relative to the ocean (RELWIND_SPD-AUX) [m/s]
        vle_wind  = eastward wind speed (WINDAVG-VLE_L1) [m/s]
        vln_wind  = northward wind speed (WINDAVG-VLN_L1) [m/s]
        vle_water = eastward surface current (VELPTMN-VLE_L1) [m/s]
        vln_water = northward surface current (VELPTMN-VLN_L1) [m/s]
        use_velptmn_with_metbk = time-vectorized data quality flag:
                                 0 -> bad  velptmn current data
                                 1 -> good velptmn current data

    Notes:

        Previous fortran and matlab implementations of the toga-coare algorithms (not associated with
        OOI) often used windspeed rather than windspeed relative to water, presumably because surface
        current data were not available. And, the only surface moorings which do have surface current
        meters are 4 Endurance moorings; however, at this time (July 2015) their compass readings are
        aliased because of the proximity of their VELPT current meters to iron mooring ballast.

        Therefore if current data are aliased or missing, met_relwind_speed will set current=0 and use
        the actual windspeed in place of the windspeed relative to water (as the DPS also specifies).
        This is only done in the routine met_relwind_speed, and not in met_relwind_direction,
        met_current_direction, nor met_current_speed. In these latter 3 routines, bad or missing
        current values are set to nan.

        This treatment of current data allows the following two relative wind calculation cases to
        be differentiated. Both cases calculate relative wind using the surface current values
        vle_water=0 and vln_water=0:

            (1) These 0 values came from actual VELPT data. In this case CURRENT_SPD calculated by
                met_current_speed will be 0.
            (2) These 0 values for current were used because the VELPT data were bad or missing. In
                this case CURRENT_SPD calculated by met_current_speed will be a nan.

        Implementation of the use_velptmn_with_metbk variable will require that it be treated as a
        "platform/instrument instance specific metadata parameter changeable in time". It has been
        coded to accept time-vectorized input.

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
    # If the surface current velocities are missing or invalid, the actual windspeed
    # will be used in place of the relative windspeed over water to calculate the METBK
    # data products.
    #
    # if this function is called without using surface current data, set current data to 0.
    if vle_water is None or vln_water is None:
        vle_water = np.zeros(vle_wind.shape[0])
        vln_water = np.zeros(vle_wind.shape[0])

    # find nans in the current record and if found replace both the east
    # and north components with values of 0.
    nanmask = np.isnan(vle_water * vln_water)
    vle_water[nanmask] = 0.0
    vln_water[nanmask] = 0.0

    # expand use_velptmn_with_metbk if it is called as a scalar
    if np.atleast_1d(use_velptmn_with_metbk).shape[0] == 1:
        use_velptmn_with_metbk = np.tile(use_velptmn_with_metbk, vle_wind.shape[0])

    # vet the surface current data - but don't use Nans
    #   when use_velptmn_with_metbk=0, set surface current velocities to 0.
    #   when use_velptmn_with_metbk=1, use surface current velocities as received.
    vle_water = vle_water * use_velptmn_with_metbk
    vln_water = vln_water * use_velptmn_with_metbk

    u_rel = np.sqrt((vle_water - vle_wind)**2 + (vln_water - vln_wind)**2)
    return u_rel


def met_timeflx(timestamp):
    """
    Description:

        Calculates TIMEFLX-AUX, the UTC timestamps corresponding to the hourly averaged
        METBK data products. The units of the timestamps are seconds since 01-01-1900.

        The timestamp values are selected to be at the midpoint of the bin intervals,
        starting half an hour after the timestamp of the first data record to be
        processed. For example, if the first data record for 30 days of data is at
        4:45 AM on a given day, the first timeflx stamp will be at 5:15 AM on that
        day, and all succeeding timestamps for the rest of the data will all be at
        15 minutes past the hour.

    Implemented by:

        2014-10-22: Russell Desiderio. Initial Code.

    Usage:

        fluxtime_hourly = met_timeflx(timestamp)

            where

        fluxtime_hourly = UTC timestamp for hourly data [seconds since 01-01-1900]
        timestamp = seconds since 01-01-1900 [UTC]

    References:

        See documentation in this module for the function make_hourly_data.
    """
    # here, the output of make_hourly_data is a list,
    # the only element of which is the desired rank 1 np.array
    fluxtime_hourly = make_hourly_data(timestamp)[0]
    return fluxtime_hourly


"""
#...................................................................................
#...................................................................................
    Functions to compute the (simpler) L2 METBK data products that do not require
    the 'warmlayer/coolskin' iteration algorithm:

        NETSIRR (this may operationally be an L1 product)
        NETSIRR_HOURLY (hourly time base)
        RAINRTE
        SALSURF
        SPECHUM

    These products are calculated at the native temporal resolution of the
    instrument suite (roughly each minute) except as noted.
#...................................................................................
#...................................................................................

"""


def met_netsirr(shortwave_down):
    """
    Description:

        Calculates NETSIRR_L2, the OOI core data product net shortwave radiation
        (wavelengths between 0.3 and 3.0 um) in the downward direction, for the METBK
        instrument. This data product may have been misclassified (it looks like L1).

        This data product is calculated on its native timebase (roughly per minute)

    Implemented by:

        2014-08-27: Russell Desiderio. Initial code.
        2017-02-03: Russell Desiderio. Added timebase documentation.

    Usage:

        net_shortwave_down = met_netsirr(shortwave_down)

            where

        net_shortwave_down = net shortwave radiation in the downward direction
                             (NETSIRR_L2) [W/m^2]
        shortwave_down = measured downward shortwave radiation (SHRTIRR_L1) [W/m^2]

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
    # net down = total down - reflected up
    albedo = 0.055    # value for reflection coefficient used in toga-coare code
    net_shortwave_down = (1.0 - albedo) * shortwave_down

    return net_shortwave_down


def met_netsirr_hourly(shortwave_down, timestamp):
    """
    Description:

        Calculates NETSIRR_HOURLY_L2, the OOI core data product net shortwave radiation
        (wavelengths between 0.3 and 3.0 um) in the downward direction, for the METBK
        instrument on an hourly time base. This data product does not require the
        coolskin\warmlayer algorithms.

    Implemented by:

        2017-02-03: Russell Desiderio. Initial code.

    Usage:

        net_shortwave_down_hourly = met_netsirr_hourly(shortwave_down, timestamp)

            where

        net_shortwave_down_hourly = net shortwave radiation in the downward direction
                                    on an hourly time base (NETSIRR_HOURLY_L2) [W/m^2]
        shortwave_down = measured downward shortwave radiation (SHRTIRR_L1) [W/m^2]
        timestamp = sample date and time value [seconds since 1900-01-01]

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
    shortwave_down, timestamp = condition_data(shortwave_down, timestamp)

    (shortwave_down_hourly, _) = make_hourly_data(shortwave_down, timestamp)

    net_shortwave_down_hourly = met_netsirr(shortwave_down_hourly)

    return net_shortwave_down_hourly


def met_rainrte(cumulative_precipitation, timestamp):
    """
    Description:

        Calculates RAINRTE_L2 (probably really an L1 product), the OOI core data
        product rain rate, for the METBK instrument. The DPS requires that the
        output data be hourly; METBK is set up to give roughly one data record
        per minute for the data needed to calculate RAINRTE.

    Implemented by:

        2014-08-27: Russell Desiderio. Initial code.
        2014-09-19: Russell Desiderio. Added front end to convert eachminute data to hourly.

    Usage:

        rainrte = met_rainrte(cumulative_precipitation, timestamp)

            where

        rainrte = rain rate (RAINRTE_L2) [mm/hr]
        cumulative_precipitation = measured rain level (PRECIPM_L1) [mm]
        timestamp = sample date and time value [seconds since 1900-01-01]

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
    cumulative_precipitation, timestamp = condition_data(cumulative_precipitation, timestamp)

    # trap out scalar case; return a value of 0 as does the DPS code.
    if cumulative_precipitation.size == 1:
        return 0.0

    cumu_prcp_hrly, time_hrly = make_hourly_data(cumulative_precipitation, timestamp)

    rainrte = calc_rain_rate(cumu_prcp_hrly, time_hrly)

    return rainrte


def met_salsurf(cond, tC_sea, ztmpwat):
    """
    Description:

        OOI Level 2 Sea Surface Salinity core data product, which is calculated
        using the Thermodynamic Equations of Seawater - 2010 (TEOS-10) Version
        3.0, with data from the conductivity, temperature and depth (CTD)
        family of instruments.

    Implemented by:

        2014-06-25: Christopher Wingard. Initial code.
        2014-08-26: Russell Desiderio. Changed variable names.

    Usage:

        SP = met_salsurf(cond, tC_sea, p_sea)

            where

        SP = practical salinity, PSS-78, (SALSURF_L2) [unitless]
        cond = sea surface conductivity (CONDSRF_L1) [S m-1]
        tC_sea = sea surface temperature (TEMPSRF_L1) [deg_C]
        ztmpwat = depth of conductivity and temperature measurements [m].
                the depth is near the surface and serves as a proxy for
                pressure [db].

    References:

        OOI (2012). Data Product Specification for Salinity. Document Control
            Number 1341-00040. https://alfresco.oceanobservatories.org/ (See: 
            Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00040_Data_Product_SPEC_PRACSAL_OOI.pdf)
    """
    # Convert L1 Conductivity from S/m to mS/cm
    C10 = cond * 10.0

    # Calculate the Practical Salinity (PSS-78) [unitless]
    SP = gsw.sp_from_c(C10, tC_sea, ztmpwat)
    return SP


def met_spechum(tC_air, pr_air, relhum):
    """
    Description:

        Calculates SPECHUM_L2, the OOI air specific humidity core data
        product, for the METBK instrument. Not to be confused with the
        SPHUM2M_L2 data product.

    Implemented by:

        2014-08-26: Russell Desiderio. Initial code.

    Usage:

        q_air = met_spechum(tC_air, pr_air, relhum)

            where

        q_air = air specific humidity (SPECHUM_L2) [g/kg]
        tC_air = air temperature (TEMPAIR_L1) [deg_C]
        pr_air = air pressure (BARPRES_L0) [mbar] {*NOT* BARPRES_L1 [Pa]}
        relhum = relative humidity (RELHUMI_L1) [%]

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
    # calculate saturated vapor pressure es in mbar
    es = 6.1121 * np.exp(17.502 * tC_air/(tC_air+240.97)) * (1.0007 + 3.46e-6 * pr_air)
    # calculate vapor pressure em from definition of relative humidity
    em = 0.01 * relhum * es
    # specific humidity q [g/kg] is then:
    q_air = 621.97 * em/(pr_air - 0.378 * em)
    return q_air


"""
#...................................................................................
#...................................................................................
    Functions to compute the L2 METBK data products that do require
    the 'warmlayer/coolskin' iteration algorithm:

        BUOYFLS:  added DPA to match FDCHP, not in original DPS
        BUOYFLX:  added DPA to match FDCHP, not in original DPS
        FRSHFLX
        HEATFLX
        LATNFLX
        MOMMFLX
        NETLIRR
        RAINFLX
        SENSFLX
        SPHUM2M
        STABLTY:  metadata
        TEMPA2M
        TEMPSKN:  metadata
        WIND10M

    These products are calculated on hourly averages.
#...................................................................................
#...................................................................................

"""


def met_buoyfls(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    Description:

        Calculates the sonic buoyancy flux (proposed) data product BUOYFLS_L2 using the
        sonic temperature instead of the virtual temperature. The FDCHP instrument
        calculates the analogous FLUXHOT_L2 data product also as a buoyancy flux using
        the sonic temperature. In contrast, the (METBK) proposed data product BUOYFLX_L2
        uses the virtual temperature in its calculation of buoyancy flux.

    Implemented by:

        2014-09-01: Russell Desiderio. Initial Code
        2014-09-19: Russell Desiderio. Added front end to convert eachminute data to hourly.

    Usage:

        Normally this routine will be called with all input arguments populated except
        for the last 3: zinvpbl is not a sensor height, and the jwarm and jcool switches
        should always be globally set to 1.

        The values for ztwmpwat, zwindsp, ztmpair, and zhumair
        may be dependent on mooring type.

        buoyfls = met_buoyfls(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                           zwindsp, ztmpair, zhumair, lat, pr_air,
                           Rshort_down, Rlong_down, cumu_prcp]

            where

        buoyfls = sonic buoyancy flux BUOYFLS_L2[W/m^2]
        tC_sea = sea temperature TEMPSRF_L1 [degC]
        wnd = windspeed relative to current RELWIND_SPD-AUX [m/s]
        tC_air = air temperature TEMPAIR [degC]
        relhum = relative humidity RELHUMI [%]
        timestamp = seconds since 01-01-1900
        lon = longitude of METBK instrument. East, positive; West, negative. [deg]
        ztmpwat = depth of sea temperature measurement TEMPSRF [m]
        zwindsp = height of windspeed measurement WINDAVG_L0 [m]
        ztmpair = height of air temperature measurement TEMPAIR [m]
        zhumair = height of air humidity measurement RELHUMI [m]
        lat = latitude of METBK instrument [deg]
        pr_air = air pressure BARPRES_L0 [mb] (mb, not pascal)
        Rshort_down = downwelling shortwave irradiation SHRTIRR_L1 [W/m^2]
        Rlong_down = downwelling longwave irradiation LONGIRR_L1 [W/m^2]
        cumu_prcp = cumulative precipitation PRECIPM_L1 [mm]

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args = make_hourly_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    (usr, tsr, qsr, _, _, _, _, _, _, _, _, _, _, _) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, _, _, _, tC_air, _, relhum, _, pr_air, _, _, _, _, _, _) = args

    rhoa = air_density(tC_air, pr_air, relhum)

    c2k = 273.15   # celsius to kelvin temperature constant
    tssr = tsr + 0.51 * (tC_air + c2k) * qsr
    cpa = 1004.67  # specific heat capacity of (dry) air [J/kg/K]
    # sonic buoyancy flux
    hsbb = -rhoa * cpa * usr * tssr

    return hsbb


def met_buoyflx(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    Description:

        Calculates the buoyancy flux (proposed) data product BUOYFLX_L2 using
        the virtual temperature. This is the more fundamental quantity for
        buoyancy flux (rather than using the sonic temperature).

    Implemented by:

        2014-09-01: Russell Desiderio. Initial Code
        2014-09-19: Russell Desiderio. Added front end to convert eachminute data to hourly.

    Usage:

        Normally this routine will be called with all input arguments populated except
        for the last 3: zinvpbl is not a sensor height, and the jwarm and jcool switches
        should always be globally set to 1.

        The values for ztwmpwat, zwindsp, ztmpair, and zhumair
        may be dependent on mooring type.

        buoyflx = met_buoyflx(tC_sea, wnd, tC_air, relhum, timestamp, lon,
                          ztmpwat, zwindsp, ztmpair, zhumair, lat,
                          pr_air, Rshort_down, Rlong_down, cumu_prcp]

            where

        buoyflx = buoyancy flux BUOYFLX_L2 [W/m^2]
        tC_sea = sea temperature TEMPSRF_L1 [degC]
        wnd = windspeed relative to current RELWIND_SPD-AUX [m/s]
        tC_air = air temperature TEMPAIR [degC]
        relhum = relative humidity RELHUMI [%]
        timestamp = seconds since 01-01-1900
        lon = longitude of METBK instrument. East, positive; West, negative. [deg]
        ztmpwat = depth of sea temperature measurement TEMPSRF [m]
        zwindsp = height of windspeed measurement WINDAVG_L0 [m]
        ztmpair = height of air temperature measurement TEMPAIR [m]
        zhumair = height of air humidity measurement RELHUMI [m]
        lat = latitude of METBK instrument [deg]
        pr_air = air pressure BARPRES_L0 [mb] (mb, not pascal)
        Rshort_down = downwelling shortwave irradiation SHRTIRR_L1 [W/m^2]
        Rlong_down = downwelling longwave irradiation LONGIRR_L1 [W/m^2]
        cumu_prcp = cumulative precipitation PRECIPM_L1 [mm]

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args = make_hourly_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    (usr, tsr, qsr, _, _, _, _, _, _, _, _, _, _, _) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, _, _, _, tC_air, _, relhum, _, pr_air, _, _, _, _, _, _) = args

    rhoa = air_density(tC_air, pr_air, relhum)

    c2k = 273.15   # celsius to kelvin temperature constant
    tvsr = tsr + 0.61 * (tC_air + c2k) * qsr
    cpa = 1004.67  # specific heat capacity of (dry) air [J/kg/K]
    # buoyancy flux
    hbb = -rhoa * cpa * usr * tvsr

    return hbb


def met_frshflx(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    Description:

        Calculates the freshwater upward flux data product FRSHFLX_L2.

    Implemented by:

        2014-09-01: Russell Desiderio. Initial Code
        2014-09-19: Russell Desiderio. Added front end to convert eachminute data to hourly.

    Usage:

        Normally this routine will be called with all input arguments populated except
        for the last 3: zinvpbl is not a sensor height, and the jwarm and jcool switches
        should always be globally set to 1.

        The values for ztwmpwat, zwindsp, ztmpair, and zhumair
        may be dependent on mooring type.

        frshflx = met_frshflx(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                           zwindsp, ztmpair, zhumair, lat, pr_air,
                           Rshort_down, Rlong_down, cumu_prcp]

            where

        frshflx = upward freshwater flux FRSHFLX_L2 [mm/hr]
        tC_sea = sea temperature TEMPSRF_L1 [degC]
        wnd = windspeed relative to current RELWIND_SPD-AUX [m/s]
        tC_air = air temperature TEMPAIR [degC]
        relhum = relative humidity RELHUMI [%]
        timestamp = seconds since 01-01-1900
        lon = longitude of METBK instrument. East, positive; West, negative. [deg]
        ztmpwat = depth of sea temperature measurement TEMPSRF [m]
        zwindsp = height of windspeed measurement WINDAVG_L0 [m]
        ztmpair = height of air temperature measurement TEMPAIR [m]
        zhumair = height of air humidity measurement RELHUMI [m]
        lat = latitude of METBK instrument [deg]
        pr_air = air pressure BARPRES_L0 [mb] (mb, not pascal)
        Rshort_down = downwelling shortwave irradiation SHRTIRR_L1 [W/m^2]
        Rlong_down = downwelling longwave irradiation LONGIRR_L1 [W/m^2]
        cumu_prcp = cumulative precipitation PRECIPM_L1 [mm]

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args = make_hourly_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    (usr, _, qsr, _, _, _, _, _, _, _, _, _, _, _) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (rain_rate, _, _, _, _, _, _, tC_air, _, relhum, _, pr_air, _, _, _, _, _, _) = args

    rhoa = air_density(tC_air, pr_air, relhum)

    # jim edson uses freshwater density; whoi dps uses seawater density;
    # perhaps the w/v concentration of pure water in seawater should be used
    # (which would be < 1000 kg/m^3).
    rho_purewater = 1000.0  # kg/m^3
    # the factor of 1000 converts from m -> mm; 3600, per sec -> per hr.
    evap = -rhoa * usr * qsr / rho_purewater * 1000.0 * 3600.0    # [mm/hr]
    frshflx = evap - rain_rate

    return frshflx


def met_heatflx(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    Description:

        Calculates the total net upward heat flux data product HEATFLX_L2.

    Implemented by:

        2014-09-01: Russell Desiderio. Initial Code
        2014-09-19: Russell Desiderio. Added front end to convert eachminute data to hourly.

    Usage:

        Normally this routine will be called with all input arguments populated except
        for the last 3: zinvpbl is not a sensor height, and the jwarm and jcool switches
        should always be globally set to 1.

        The values for ztwmpwat, zwindsp, ztmpair, and zhumair
        may be dependent on mooring type.

        heatflx = met_heatflx(tC_sea, wnd, tC_air, relhum, timestamp, lon,
                              ztmpwat, zwindsp, ztmpair, zhumair, lat,
                              pr_air, Rshort_down, Rlong_down, cumu_prcp]

            where

        heatflx = total net upward heat flux HEATFLX_L2 [W/m^2]
        tC_sea = sea temperature TEMPSRF_L1 [degC]
        wnd = windspeed relative to current RELWIND_SPD-AUX [m/s]
        tC_air = air temperature TEMPAIR [degC]
        relhum = relative humidity RELHUMI [%]
        timestamp = seconds since 01-01-1900
        lon = longitude of METBK instrument. East, positive; West, negative. [deg]
        ztmpwat = depth of sea temperature measurement TEMPSRF [m]
        zwindsp = height of windspeed measurement WINDAVG_L0 [m]
        ztmpair = height of air temperature measurement TEMPAIR [m]
        zhumair = height of air humidity measurement RELHUMI [m]
        lat = latitude of METBK instrument [deg]
        pr_air = air pressure BARPRES_L0 [mb] (mb, not pascal)
        Rshort_down = downwelling shortwave irradiation SHRTIRR_L1 [W/m^2]
        Rlong_down = downwelling longwave irradiation LONGIRR_L1 [W/m^2]
        cumu_prcp = cumulative precipitation PRECIPM_L1 [mm]

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args = make_hourly_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    (usr, tsr, qsr, _, dter, dqer, _, _, _, _, _, _, _, dsea) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (rain_rate, _, _, _, tC_sea, _, _, tC_air, _, relhum, _, pr_air, Rshort_down,
        Rlong_down, _, _, _, _) = args

    cpa = 1004.67  # specific heat capacity of (dry) air [J/kg/K]
    rhoa = air_density(tC_air, pr_air, relhum)
    Le = latent_heat_vaporization_pure_water(tC_sea + dsea)

    hlb = -rhoa * Le * usr * qsr                                              # positive up
    hsb = -rhoa * cpa * usr * tsr                                             # positive up
    Rns_down = met_netsirr(Rshort_down)                                       # positive down
    Rnl_up = net_longwave_up(tC_sea + dsea - dter, Rlong_down)                # positive up
    rainflx = rain_heat_flux(rain_rate, tC_sea+dsea, tC_air, relhum, pr_air)  # positive up

    heatflx = hlb + hsb - Rns_down + Rnl_up + rainflx

    return heatflx


def met_latnflx(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    Description:

        Calculates the upward latent heat flux data product LATNFLX_L2.

    Implemented by:

        2014-09-01: Russell Desiderio. Initial Code
        2014-09-19: Russell Desiderio. Added front end to convert eachminute data to hourly.

    Usage:

        Normally this routine will be called with all input arguments populated except
        for the last 3: zinvpbl is not a sensor height, and the jwarm and jcool switches
        should always be globally set to 1.

        The values for ztwmpwat, zwindsp, ztmpair, and zhumair
        may be dependent on mooring type.

        latnflx = met_latnflx(tC_sea, wnd, tC_air, relhum, timestamp, lon,
                              ztmpwat, zwindsp, ztmpair, zhumair, lat,
                              pr_air, Rshort_down, Rlong_down, cumu_prcp]

            where

        latnflx = upward latent heat flux LATNFLX_L2 [W/m^2]
        tC_sea = sea temperature TEMPSRF_L1 [degC]
        wnd = windspeed relative to current RELWIND_SPD-AUX [m/s]
        tC_air = air temperature TEMPAIR [degC]
        relhum = relative humidity RELHUMI [%]
        timestamp = seconds since 01-01-1900
        lon = longitude of METBK instrument. East, positive; West, negative. [deg]
        ztmpwat = depth of sea temperature measurement TEMPSRF [m]
        zwindsp = height of windspeed measurement WINDAVG_L0 [m]
        ztmpair = height of air temperature measurement TEMPAIR [m]
        zhumair = height of air humidity measurement RELHUMI [m]
        lat = latitude of METBK instrument [deg]
        pr_air = air pressure BARPRES_L0 [mb] (mb, not pascal)
        Rshort_down = downwelling shortwave irradiation SHRTIRR_L1 [W/m^2]
        Rlong_down = downwelling longwave irradiation LONGIRR_L1 [W/m^2]
        cumu_prcp = cumulative precipitation PRECIPM_L1 [mm]

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args = make_hourly_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    # dsea is the warmlayer correction to the sea surface temperature
    (usr, _, qsr, _, _, _, _, _, _, _, _, _, _, dsea) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, tC_sea, _, _, tC_air, _, relhum, _, pr_air, _, _, _, _, _, _) = args

    rhoa = air_density(tC_air, pr_air, relhum)

    # note that the original (coare ver. 3.5) code:
    #    (a) uses Le for pure water, not seawater.
    #    (b) does not include the coolskin correction to sea surface temperature.
    Le = latent_heat_vaporization_pure_water(tC_sea + dsea)

    hlb = -rhoa * Le * usr * qsr

    return hlb


def met_mommflx(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    Description:

        Calculates (the absolute value of) the momentum flux MOMMFLX_L2,
        also called the wind stress tau.

    Implemented by:

        2014-09-01: Russell Desiderio. Initial Code
        2014-09-19: Russell Desiderio. Added front end to convert eachminute data to hourly.

    Usage:

        Normally this routine will be called with all input arguments populated except
        for the last 3: zinvpbl is not a sensor input, and the jwarm and jcool switches
        should always be globally set to 1.

        The values for ztwmpwat, zwindsp, ztmpair, and zhumair
        may be dependent on mooring type.

        mommflx = met_mommflx(tC_sea, wnd, tC_air, relhum, timestamp, lon,
                              ztmpwat, zwindsp, ztmpair, zhumair, lat,
                              pr_air, Rshort_down, Rlong_down, cumu_prcp]

            where

        mommflx = the momentum flux MOMMFLX_L2 [N/m^2]
        tC_sea = sea temperature TEMPSRF_L1 [degC]
        wnd = windspeed relative to current RELWIND_SPD-AUX [m/s]
        tC_air = air temperature TEMPAIR [degC]
        relhum = relative humidity RELHUMI [%]
        timestamp = seconds since 01-01-1900
        lon = longitude of METBK instrument. East, positive; West, negative. [deg]
        ztmpwat = depth of sea temperature measurement TEMPSRF [m]
        zwindsp = height of windspeed measurement WINDAVG_L0 [m]
        ztmpair = height of air temperature measurement TEMPAIR [m]
        zhumair = height of air humidity measurement RELHUMI [m]
        lat = latitude of METBK instrument [deg]
        pr_air = air pressure BARPRES_L0 [mb] (mb, not pascal)
        Rshort_down = downwelling shortwave irradiation SHRTIRR_L1 [W/m^2]
        Rlong_down = downwelling longwave irradiation LONGIRR_L1 [W/m^2]
        cumu_prcp = cumulative precipitation PRECIPM_L1 [mm]

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args = make_hourly_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    (usr, _, _, ut, _, _, _, _, _, _, _, _, _, _) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, _, wnd, _, tC_air, _, relhum, _, pr_air, _, _, _, _, _, _) = args

    rhoa = air_density(tC_air, pr_air, relhum)

    # the wind stress tau is the magnitude of the momentum flux.
    tau = rhoa * usr * usr * wnd / ut

    return tau


def met_netlirr(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    Description:

        Calculates the net upward longwave irradiance NETLIRR_L2.

    Implemented by:

        2014-09-01: Russell Desiderio. Initial Code
        2014-09-19: Russell Desiderio. Added front end to convert eachminute data to hourly.

    Usage:

        Normally this routine will be called with all input arguments populated except
        for the last 3: zinvpbl is not a sensor input, and the jwarm and jcool switches
        should always be globally set to 1.

        The values for ztwmpwat, zwindsp, ztmpair, and zhumair
        may be dependent on mooring type.

        netlirr = met_netlirr(tC_sea, wnd, tC_air, relhum, timestamp, lon,
                              ztmpwat, zwindsp, ztmpair, zhumair, lat,
                              pr_air, Rshort_down, Rlong_down, cumu_prcp]

            where

        netlirr = net upward longwave irradiance NETLIRR_L2 [W/m^2]
        tC_sea = sea temperature TEMPSRF_L1 [degC]
        wnd = windspeed relative to current RELWIND_SPD-AUX [m/s]
        tC_air = air temperature TEMPAIR [degC]
        relhum = relative humidity RELHUMI [%]
        timestamp = seconds since 01-01-1900
        lon = longitude of METBK instrument. East, positive; West, negative. [deg]
        ztmpwat = depth of sea temperature measurement TEMPSRF [m]
        zwindsp = height of windspeed measurement WINDAVG_L0 [m]
        ztmpair = height of air temperature measurement TEMPAIR [m]
        zhumair = height of air humidity measurement RELHUMI [m]
        lat = latitude of METBK instrument [deg]
        pr_air = air pressure BARPRES_L0 [mb] (mb, not pascal)
        Rshort_down = downwelling shortwave irradiation SHRTIRR_L1 [W/m^2]
        Rlong_down = downwelling longwave irradiation LONGIRR_L1 [W/m^2]
        cumu_prcp = cumulative precipitation PRECIPM_L1 [mm]

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args = make_hourly_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    # dter is the coolskin temperature depression [degC]
    # dsea is the warmlayer correction to the sea surface temperature [degC]
    (_, _, _, _, dter, _, _, _, _, _, _, _, _, dsea) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, tC_sea, _, _, _, _, _, _, _, _, Rlong_down, _, _, _, _) = args

    Rnl = net_longwave_up(tC_sea + dsea - dter, Rlong_down)

    return Rnl


def met_rainflx(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    Description:

        Calculates the net upward rain heat flux RAINFLX_L2.

        A new derivation for rain heat flux is used, as calculated in the
        subroutine rain_heat_flux.

    Implemented by:

        2014-09-01: Russell Desiderio. Initial Code
        2014-09-19: Russell Desiderio. Added front end to convert eachminute data to hourly.
        2014-10-28: Russell Desiderio. Incorporated new subroutine for rain heat flux.

    Usage:

        Normally this routine will be called with all input arguments populated except
        for the last 3: zinvpbl is not a sensor input, and the jwarm and jcool switches
        should always be globally set to 1.

        The values for ztwmpwat, zwindsp, ztmpair, and zhumair
        may be dependent on mooring type.

        rainflx = met_rainflx(tC_sea, wnd, tC_air, relhum, timestamp, lon,
                              ztmpwat, zwindsp, ztmpair, zhumair, lat,
                              pr_air, Rshort_down, Rlong_down, cumu_prcp]

            where

        rainflx = net upward rain heat flux RAINFLX_L2 [W/m^2]
        tC_sea = sea temperature TEMPSRF_L1 [degC]
        wnd = windspeed relative to current RELWIND_SPD-AUX [m/s]
        tC_air = air temperature TEMPAIR [degC]
        relhum = relative humidity RELHUMI [%]
        timestamp = seconds since 01-01-1900
        lon = longitude of METBK instrument. East, positive; West, negative. [deg]
        ztmpwat = depth of sea temperature measurement TEMPSRF [m]
        zwindsp = height of windspeed measurement WINDAVG_L0 [m]
        ztmpair = height of air temperature measurement TEMPAIR [m]
        zhumair = height of air humidity measurement RELHUMI [m]
        lat = latitude of METBK instrument [deg]
        pr_air = air pressure BARPRES_L0 [mb] (mb, not pascal)
        Rshort_down = downwelling shortwave irradiation SHRTIRR_L1 [W/m^2]
        Rlong_down = downwelling longwave irradiation LONGIRR_L1 [W/m^2]
        cumu_prcp = cumulative precipitation PRECIPM_L1 [mm]

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args = make_hourly_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    # dsea is the warmlayer correction to the sea surface temperature [degC]
    (_, _, _, _, _, _, _, _, _, _, _, _, _, dsea) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (rain_rate, _, _, _, tC_sea, _, _, tC_air, _, relhum, _, pr_air,
        _, _, _, _, _, _) = args

    # the raindrops penetrate the sea surface on the order of cm, which is where the
    # heat is 'exchanged'. therefore, use the warmlayer correction but not the
    # coolskin correction (which is order microns (?) thick) to the sea temperature.
    rainflx = rain_heat_flux(rain_rate, tC_sea+dsea, tC_air, relhum, pr_air)

    return rainflx


def met_sensflx(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    Description:

        Calculates the net upward sensible heat flux SENSFLX_L2.

    Implemented by:

        2014-09-01: Russell Desiderio. Initial Code
        2014-09-19: Russell Desiderio. Added front end to convert eachminute data to hourly.

    Usage:

        Normally this routine will be called with all input arguments populated except
        for the last 3: zinvpbl is not a sensor input, and the jwarm and jcool switches
        should always be globally set to 1.

        The values for ztwmpwat, zwindsp, ztmpair, and zhumair
        may be dependent on mooring type.

        sensflx = met_sensflx(tC_sea, wnd, tC_air, relhum, timestamp, lon,
                              ztmpwat, zwindsp, ztmpair, zhumair, lat,
                              pr_air, Rshort_down, Rlong_down, cumu_prcp]

            where

        sensflx = net upward sensible heat flux SENSFLX_L2 [W/m^2]
        tC_sea = sea temperature TEMPSRF_L1 [degC]
        wnd = windspeed relative to current RELWIND_SPD-AUX [m/s]
        tC_air = air temperature TEMPAIR [degC]
        relhum = relative humidity RELHUMI [%]
        timestamp = seconds since 01-01-1900
        lon = longitude of METBK instrument. East, positive; West, negative. [deg]
        ztmpwat = depth of sea temperature measurement TEMPSRF [m]
        zwindsp = height of windspeed measurement WINDAVG_L0 [m]
        ztmpair = height of air temperature measurement TEMPAIR [m]
        zhumair = height of air humidity measurement RELHUMI [m]
        lat = latitude of METBK instrument [deg]
        pr_air = air pressure BARPRES_L0 [mb] (mb, not pascal)
        Rshort_down = downwelling shortwave irradiation SHRTIRR_L1 [W/m^2]
        Rlong_down = downwelling longwave irradiation LONGIRR_L1 [W/m^2]
        cumu_prcp = cumulative precipitation PRECIPM_L1 [mm]

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args = make_hourly_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    (usr, tsr, _, _, _, _, _, _, _, _, _, _, _, _) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, _, _, _, tC_air, _, relhum, _, pr_air, _, _, _, _, _, _) = args

    rhoa = air_density(tC_air, pr_air, relhum)

    cpa = 1004.67  # specific heat capacity of (dry) air [J/kg/K]
    hsb = -rhoa * cpa * usr * tsr

    return hsb


def met_sphum2m(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    Description:

        Calculates the modelled specific humidity at a reference height of 2m SPHUM2M_L2.

    Implemented by:

        2014-09-01: Russell Desiderio. Initial Code
        2014-09-19: Russell Desiderio. Added front end to convert eachminute data to hourly.

    Usage:

        Normally this routine will be called with all input arguments populated except
        for the last 3: zinvpbl is not a sensor input, and the jwarm and jcool switches
        should always be globally set to 1.

        The values for ztwmpwat, zwindsp, ztmpair, and zhumair
        may be dependent on mooring type.

        sphum2m = met_sphum2m(tC_sea, wnd, tC_air, relhum, timestamp, lon,
                              ztmpwat, zwindsp, ztmpair, zhumair, lat,
                              pr_air, Rshort_down, Rlong_down, cumu_prcp]

            where

        sphum2m = modelled specific humidity at 2m SPHUM2M_L2 [g/kg]
        tC_sea = sea temperature TEMPSRF_L1 [degC]
        wnd = windspeed relative to current RELWIND_SPD-AUX [m/s]
        tC_air = air temperature TEMPAIR [degC]
        relhum = relative humidity RELHUMI [%]
        timestamp = seconds since 01-01-1900
        lon = longitude of METBK instrument. East, positive; West, negative. [deg]
        ztmpwat = depth of sea temperature measurement TEMPSRF [m]
        zwindsp = height of windspeed measurement WINDAVG_L0 [m]
        ztmpair = height of air temperature measurement TEMPAIR [m]
        zhumair = height of air humidity measurement RELHUMI [m]
        lat = latitude of METBK instrument [deg]
        pr_air = air pressure BARPRES_L0 [mb] (mb, not pascal)
        Rshort_down = downwelling shortwave irradiation SHRTIRR_L1 [W/m^2]
        Rlong_down = downwelling longwave irradiation LONGIRR_L1 [W/m^2]
        cumu_prcp = cumulative precipitation PRECIPM_L1 [mm]

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
    zrefht = 2.0  # [m]

    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args = make_hourly_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    # L is the Obukhov length scale [m]
    (_, _, qsr, _, _, _, _, L, _, _, _, _, _, _) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, _, _, _, tC_air, _, relhum, zhumair, pr_air, _, _, _, _, _, _) = args

    sphum2m = spechum_at_refheight(tC_air, pr_air, relhum, qsr, zrefht, zhumair, L)

    return sphum2m


def met_stablty(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    Description:

        Calculates the Monin-Obukhov stability parameter metadata product STABLTY_L2.

    Implemented by:

        2014-09-01: Russell Desiderio. Initial Code
        2014-09-19: Russell Desiderio. Added front end to convert eachminute data to hourly.

    Usage:

        Normally this routine will be called with all input arguments populated except
        for the last 3: zinvpbl is not a sensor input, and the jwarm and jcool switches
        should always be globally set to 1.

        The values for ztwmpwat, zwindsp, ztmpair, and zhumair
        may be dependent on mooring type.

        stablty = met_stablty(tC_sea, wnd, tC_air, relhum, timestamp, lon,
                              ztmpwat, zwindsp, ztmpair, zhumair, lat,
                              pr_air, Rshort_down, Rlong_down, cumu_prcp]

            where

        stablty = Monin-Obukhov stability parameter STABLTY_L2 [unitless]
        tC_sea = sea temperature TEMPSRF_L1 [degC]
        wnd = windspeed relative to current RELWIND_SPD-AUX [m/s]
        tC_air = air temperature TEMPAIR [degC]
        relhum = relative humidity RELHUMI [%]
        timestamp = seconds since 01-01-1900
        lon = longitude of METBK instrument. East, positive; West, negative. [deg]
        ztmpwat = depth of sea temperature measurement TEMPSRF [m]
        zwindsp = height of windspeed measurement WINDAVG_L0 [m]
        ztmpair = height of air temperature measurement TEMPAIR [m]
        zhumair = height of air humidity measurement RELHUMI [m]
        lat = latitude of METBK instrument [deg]
        pr_air = air pressure BARPRES_L0 [mb] (mb, not pascal)
        Rshort_down = downwelling shortwave irradiation SHRTIRR_L1 [W/m^2]
        Rlong_down = downwelling longwave irradiation LONGIRR_L1 [W/m^2]
        cumu_prcp = cumulative precipitation PRECIPM_L1 [mm]

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args = make_hourly_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    # L is the Obukhov length scale [m]
    (_, _, _, _, _, _, _, L, _, _, _, _, _, _) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, _, _, zwindsp, _, _, _, _, _, _, _, _, _, _, _) = args

    return zwindsp / L


def met_tempa2m(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    Description:

        Calculates the modelled air temperature at a reference height of 2m TEMPA2M_L2.

    Implemented by:

        2014-09-01: Russell Desiderio. Initial Code
        2014-09-19: Russell Desiderio. Added front end to convert eachminute data to hourly.

    Usage:

        Normally this routine will be called with all input arguments populated except
        for the last 3: zinvpbl is not a sensor input, and the jwarm and jcool switches
        should always be globally set to 1.

        The values for ztwmpwat, zwindsp, ztmpair, and zhumair
        may be dependent on mooring type.

        tempa2m = met_tempa2m(tC_sea, wnd, tC_air, relhum, timestamp, lon,
                              ztmpwat, zwindsp, ztmpair, zhumair, lat,
                              pr_air, Rshort_down, Rlong_down, cumu_prcp]

            where

        tempa2m = modelled air temperature at 2m TEMPA2M_L2 [degC]
        tC_sea = sea temperature TEMPSRF_L1 [degC]
        wnd = windspeed relative to current RELWIND_SPD-AUX [m/s]
        tC_air = air temperature TEMPAIR [degC]
        relhum = relative humidity RELHUMI [%]
        timestamp = seconds since 01-01-1900
        lon = longitude of METBK instrument. East, positive; West, negative. [deg]
        ztmpwat = depth of sea temperature measurement TEMPSRF [m]
        zwindsp = height of windspeed measurement WINDAVG_L0 [m]
        ztmpair = height of air temperature measurement TEMPAIR [m]
        zhumair = height of air humidity measurement RELHUMI [m]
        lat = latitude of METBK instrument [deg]
        pr_air = air pressure BARPRES_L0 [mb] (mb, not pascal)
        Rshort_down = downwelling shortwave irradiation SHRTIRR_L1 [W/m^2]
        Rlong_down = downwelling longwave irradiation LONGIRR_L1 [W/m^2]
        cumu_prcp = cumulative precipitation PRECIPM_L1 [mm]

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
    zrefht = 2.0  # [m]

    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args = make_hourly_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    # L is the Obukhov length scale [m]
    (_, tsr, _, _, _, _, _, L, _, _, _, _, _, _) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, _, _, _, tC_air, ztmpair, _, _, _, _, _, lat, _, _, _) = args

    tempa2m = airtemp_at_refheight(tC_air, tsr, zrefht, ztmpair, L, lat)

    return tempa2m


def met_tempskn(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    Description:

        Calculates the skin sea temperature based on the warmlayer and coolskin (coare35vn)
        model: metadata product TEMPSKN_L2.

    Implemented by:

        2014-09-01: Russell Desiderio. Initial Code
        2014-09-19: Russell Desiderio. Added front end to convert eachminute data to hourly.

    Usage:

        Normally this routine will be called with all input arguments populated except
        for the last 3: zinvpbl is not a sensor input, and the jwarm and jcool switches
        should always be globally set to 1.

        The values for ztwmpwat, zwindsp, ztmpair, and zhumair
        may be dependent on mooring type.

        tempskn = met_tempskn(tC_sea, wnd, tC_air, relhum, timestamp, lon,
                              ztmpwat, zwindsp, ztmpair, zhumair, lat,
                              pr_air, Rshort_down, Rlong_down, cumu_prcp]

            where

        tempskn = skin seasurface temperature TEMPSKN_L2 [degC]
        tC_sea = sea temperature TEMPSRF_L1 [degC]
        wnd = windspeed relative to current RELWIND_SPD-AUX [m/s]
        tC_air = air temperature TEMPAIR [degC]
        relhum = relative humidity RELHUMI [%]
        timestamp = seconds since 01-01-1900
        lon = longitude of METBK instrument. East, positive; West, negative. [deg]
        ztmpwat = depth of sea temperature measurement TEMPSRF [m]
        zwindsp = height of windspeed measurement WINDAVG_L0 [m]
        ztmpair = height of air temperature measurement TEMPAIR [m]
        zhumair = height of air humidity measurement RELHUMI [m]
        lat = latitude of METBK instrument [deg]
        pr_air = air pressure BARPRES_L0 [mb] (mb, not pascal)
        Rshort_down = downwelling shortwave irradiation SHRTIRR_L1 [W/m^2]
        Rlong_down = downwelling longwave irradiation LONGIRR_L1 [W/m^2]
        cumu_prcp = cumulative precipitation PRECIPM_L1 [mm]

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args = make_hourly_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    # dter is the coolskin temperature depression [degC]
    # dsea is the warmlayer correction to the sea surface temperature [degC]
    (_, _, _, _, dter, _, _, _, _, _, _, _, _, dsea) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, tC_sea, _, _, _, _, _, _, _, _, _, _, _, _, _) = args

    # warmlayer corrections are added; coolskin corrections are subtracted
    tempskn = tC_sea + dsea - dter

    return tempskn


def met_wind10m(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    Description:

        Calculates the modelled windspeed at a reference height of 10m WIND10M_L2.

    Implemented by:

        2014-09-01: Russell Desiderio. Initial Code
        2014-09-19: Russell Desiderio. Added front end to convert eachminute data to hourly.

    Usage:

        Normally this routine will be called with all input arguments populated except
        for the last 3: zinvpbl is not a sensor input, and the jwarm and jcool switches
        should always be globally set to 1.

        The values for ztwmpwat, zwindsp, ztmpair, and zhumair
        may be dependent on mooring type.

        wind10m = met_wind10m(tC_sea, wnd, tC_air, relhum, timestamp, lon,
                              ztmpwat, zwindsp, ztmpair, zhumair, lat,
                              pr_air, Rshort_down, Rlong_down, cumu_prcp]

            where

        wind10m = modelled windspeed at 10m WIND10M_L2 [degC]
        tC_sea = sea temperature TEMPSRF_L1 [degC]
        wnd = windspeed relative to current RELWIND_SPD-AUX [m/s]
        tC_air = air temperature TEMPAIR [degC]
        relhum = relative humidity RELHUMI [%]
        timestamp = seconds since 01-01-1900
        lon = longitude of METBK instrument. East, positive; West, negative. [deg]
        ztmpwat = depth of sea temperature measurement TEMPSRF [m]
        zwindsp = height of windspeed measurement WINDAVG_L0 [m]
        ztmpair = height of air temperature measurement TEMPAIR [m]
        zhumair = height of air humidity measurement RELHUMI [m]
        lat = latitude of METBK instrument [deg]
        pr_air = air pressure BARPRES_L0 [mb] (mb, not pascal)
        Rshort_down = downwelling shortwave irradiation SHRTIRR_L1 [W/m^2]
        Rlong_down = downwelling longwave irradiation LONGIRR_L1 [W/m^2]
        cumu_prcp = cumulative precipitation PRECIPM_L1 [mm]

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
    zrefht = 10.0  # [m]

    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args = make_hourly_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    # L is the Obukhov length scale [m]
    (usr, _, _, ut, _, _, _, L, _, _, _, _, _, _) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, _, wnd, zwindsp, _, _, _, _, _, _, _, _, _, _, _) = args

    wind10m = windspeed_at_refheight(wnd, usr, ut, zrefht, zwindsp, L)

    return wind10m

"""
#...................................................................................
#...................................................................................
    Functions to compute the L2 METBK data products that do require
    the 'warmlayer/coolskin' iteration algorithm, but at the native
    temporal resolution of the instrument suite (roughly per minute).

    These functions were not specified in the DPS.
        HEATFLX_MINUTE
        LATNFLX_MINUTE
        NETLIRR_MINUTE
        SENSFLX_MINUTE
#...................................................................................
#...................................................................................
"""


def met_heatflx_minute(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                       zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                       Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                       zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    Description:

        Calculates the total net upward heat flux data product HEATFLX_MINUTE_L2,
        calculated on the native METBK time base (per minute).

    Implemented by:

        2017-02-13: Russell Desiderio. Initial code.

    Usage:

        Normally this routine will be called with all input arguments populated except
        for the last 3: zinvpbl is not a sensor height, and the jwarm and jcool switches
        should always be globally set to 1.

        The values for ztwmpwat, zwindsp, ztmpair, and zhumair
        may be dependent on mooring type.

        heatflx_minute = met_heatflx_minute(tC_sea, wnd, tC_air, relhum, timestamp, lon,
                                            ztmpwat, zwindsp, ztmpair, zhumair, lat,
                                            pr_air, Rshort_down, Rlong_down, cumu_prcp]

            where

        heatflx_minute = total net upward heat flux HEATFLX_MINUTE_L2 [W/m^2]
        tC_sea = sea temperature TEMPSRF_L1 [degC]
        wnd = windspeed relative to current RELWIND_SPD-AUX [m/s]
        tC_air = air temperature TEMPAIR [degC]
        relhum = relative humidity RELHUMI [%]
        timestamp = seconds since 01-01-1900
        lon = longitude of METBK instrument. East, positive; West, negative. [deg]
        ztmpwat = depth of sea temperature measurement TEMPSRF [m]
        zwindsp = height of windspeed measurement WINDAVG_L0 [m]
        ztmpair = height of air temperature measurement TEMPAIR [m]
        zhumair = height of air humidity measurement RELHUMI [m]
        lat = latitude of METBK instrument [deg]
        pr_air = air pressure BARPRES_L0 [mb] (mb, not pascal)
        Rshort_down = downwelling shortwave irradiation SHRTIRR_L1 [W/m^2]
        Rlong_down = downwelling longwave irradiation LONGIRR_L1 [W/m^2]
        cumu_prcp = cumulative precipitation PRECIPM_L1 [mm]

    Notes:

        This data product is not specifically included in the DPS.

        This routine differs from met_heatflx only in that the function
        make_hourly_data is not called to bin the input data into hourly bins.

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    (usr, tsr, qsr, _, dter, dqer, _, _, _, _, _, _, _, dsea) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (rain_rate, _, _, _, tC_sea, _, _, tC_air, _, relhum, _, pr_air, Rshort_down,
        Rlong_down, _, _, _, _) = args

    cpa = 1004.67  # specific heat capacity of (dry) air [J/kg/K]
    rhoa = air_density(tC_air, pr_air, relhum)
    Le = latent_heat_vaporization_pure_water(tC_sea + dsea)

    hlb = -rhoa * Le * usr * qsr                                              # positive up
    hsb = -rhoa * cpa * usr * tsr                                             # positive up
    Rns_down = met_netsirr(Rshort_down)                                       # positive down
    Rnl_up = net_longwave_up(tC_sea + dsea - dter, Rlong_down)                # positive up
    rainflx = rain_heat_flux(rain_rate, tC_sea+dsea, tC_air, relhum, pr_air)  # positive up

    heatflx = hlb + hsb - Rns_down + Rnl_up + rainflx

    return heatflx


def met_latnflx_minute(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                       zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                       Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                       zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    Description:

        Calculates the upward latent heat flux data product LATNFLX_MINUTE_L2,
        calculated on the native METBK time base (per minute).

    Implemented by:

        2017-02-13: Russell Desiderio. Initial code.

    Usage:

        Normally this routine will be called with all input arguments populated except
        for the last 3: zinvpbl is not a sensor height, and the jwarm and jcool switches
        should always be globally set to 1.

        The values for ztwmpwat, zwindsp, ztmpair, and zhumair
        may be dependent on mooring type.

        latnflx_minute = met_latnflx_minute(tC_sea, wnd, tC_air, relhum, timestamp, lon,
                                            ztmpwat, zwindsp, ztmpair, zhumair, lat,
                                            pr_air, Rshort_down, Rlong_down, cumu_prcp]

            where

        latnflx_minute = upward latent heat flux LATNFLX_MINUTE_L2 [W/m^2]
        tC_sea = sea temperature TEMPSRF_L1 [degC]
        wnd = windspeed relative to current RELWIND_SPD-AUX [m/s]
        tC_air = air temperature TEMPAIR [degC]
        relhum = relative humidity RELHUMI [%]
        timestamp = seconds since 01-01-1900
        lon = longitude of METBK instrument. East, positive; West, negative. [deg]
        ztmpwat = depth of sea temperature measurement TEMPSRF [m]
        zwindsp = height of windspeed measurement WINDAVG_L0 [m]
        ztmpair = height of air temperature measurement TEMPAIR [m]
        zhumair = height of air humidity measurement RELHUMI [m]
        lat = latitude of METBK instrument [deg]
        pr_air = air pressure BARPRES_L0 [mb] (mb, not pascal)
        Rshort_down = downwelling shortwave irradiation SHRTIRR_L1 [W/m^2]
        Rlong_down = downwelling longwave irradiation LONGIRR_L1 [W/m^2]
        cumu_prcp = cumulative precipitation PRECIPM_L1 [mm]

    Notes:

        This data product is not specifically included in the DPS.

        This routine differs from met_latnflx only in that the function
        make_hourly_data is not called to bin the input data into hourly bins.

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    # dsea is the warmlayer correction to the sea surface temperature
    (usr, _, qsr, _, _, _, _, _, _, _, _, _, _, dsea) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, tC_sea, _, _, tC_air, _, relhum, _, pr_air, _, _, _, _, _, _) = args

    rhoa = air_density(tC_air, pr_air, relhum)

    # note that the original (coare ver. 3.5) code:
    #    (a) uses Le for pure water, not seawater.
    #    (b) does not include the coolskin correction to sea surface temperature.
    Le = latent_heat_vaporization_pure_water(tC_sea + dsea)

    hlb = -rhoa * Le * usr * qsr

    return hlb


def met_netlirr_minute(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                       zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                       Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                       zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    Description:

        Calculates the data product NETLIRR_MINUTE_L2, the net upward longwave irradiance
        calculated on the native METBK time base (per minute).

    Implemented by:

        2017-02-13: Russell Desiderio. Initial code.

    Usage:

        Normally this routine will be called with all input arguments populated except
        for the last 3: zinvpbl is not a sensor input, and the jwarm and jcool switches
        should always be globally set to 1.

        The values for ztwmpwat, zwindsp, ztmpair, and zhumair
        may be dependent on mooring type.

        netlirr_minute = met_netlirr_minute(tC_sea, wnd, tC_air, relhum, timestamp, lon,
                                            ztmpwat, zwindsp, ztmpair, zhumair, lat,
                                            pr_air, Rshort_down, Rlong_down, cumu_prcp]

            where

        netlirr_minute = net upward longwave irradiance NETLIRR_MINUTE_L2 [W/m^2]
        tC_sea = sea temperature TEMPSRF_L1 [degC]
        wnd = windspeed relative to current RELWIND_SPD-AUX [m/s]
        tC_air = air temperature TEMPAIR [degC]
        relhum = relative humidity RELHUMI [%]
        timestamp = seconds since 01-01-1900
        lon = longitude of METBK instrument. East, positive; West, negative. [deg]
        ztmpwat = depth of sea temperature measurement TEMPSRF [m]
        zwindsp = height of windspeed measurement WINDAVG_L0 [m]
        ztmpair = height of air temperature measurement TEMPAIR [m]
        zhumair = height of air humidity measurement RELHUMI [m]
        lat = latitude of METBK instrument [deg]
        pr_air = air pressure BARPRES_L0 [mb] (mb, not pascal)
        Rshort_down = downwelling shortwave irradiation SHRTIRR_L1 [W/m^2]
        Rlong_down = downwelling longwave irradiation LONGIRR_L1 [W/m^2]
        cumu_prcp = cumulative precipitation PRECIPM_L1 [mm]

    Notes:

        This data product is not specifically included in the DPS.

        This routine differs from met_netlirr only in that the function
        make_hourly_data is not called to bin the input data into hourly bins.

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    # dter is the coolskin temperature depression [degC]
    # dsea is the warmlayer correction to the sea surface temperature [degC]
    (_, _, _, _, dter, _, _, _, _, _, _, _, _, dsea) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, tC_sea, _, _, _, _, _, _, _, _, Rlong_down, _, _, _, _) = args

    Rnl_native = net_longwave_up(tC_sea + dsea - dter, Rlong_down)

    return Rnl_native


def met_sensflx_minute(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                       zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                       Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                       zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    Description:

        Calculates the net upward sensible heat flux SENSFLX_MINUTE_L2.
        calculated on the native METBK time base (per minute).

    Implemented by:

        2017-02-13: Russell Desiderio. Initial code.

    Usage:

        Normally this routine will be called with all input arguments populated except
        for the last 3: zinvpbl is not a sensor input, and the jwarm and jcool switches
        should always be globally set to 1.

        The values for ztwmpwat, zwindsp, ztmpair, and zhumair
        may be dependent on mooring type.

        sensflx_minute = met_sensflx_minute(tC_sea, wnd, tC_air, relhum, timestamp, lon,
                                            ztmpwat, zwindsp, ztmpair, zhumair, lat,
                                            pr_air, Rshort_down, Rlong_down, cumu_prcp]

            where

        sensflx_minute = net upward sensible heat flux SENSFLX_MINUTE_L2 [W/m^2]
        tC_sea = sea temperature TEMPSRF_L1 [degC]
        wnd = windspeed relative to current RELWIND_SPD-AUX [m/s]
        tC_air = air temperature TEMPAIR [degC]
        relhum = relative humidity RELHUMI [%]
        timestamp = seconds since 01-01-1900
        lon = longitude of METBK instrument. East, positive; West, negative. [deg]
        ztmpwat = depth of sea temperature measurement TEMPSRF [m]
        zwindsp = height of windspeed measurement WINDAVG_L0 [m]
        ztmpair = height of air temperature measurement TEMPAIR [m]
        zhumair = height of air humidity measurement RELHUMI [m]
        lat = latitude of METBK instrument [deg]
        pr_air = air pressure BARPRES_L0 [mb] (mb, not pascal)
        Rshort_down = downwelling shortwave irradiation SHRTIRR_L1 [W/m^2]
        Rlong_down = downwelling longwave irradiation LONGIRR_L1 [W/m^2]
        cumu_prcp = cumulative precipitation PRECIPM_L1 [mm]

    Notes:

        This data product is not specifically included in the DPS.

        This routine differs from met_sensflx only in that the function
        make_hourly_data is not called to bin the input data into hourly bins.

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    (usr, tsr, _, _, _, _, _, _, _, _, _, _, _, _) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, _, _, _, tC_air, _, relhum, _, pr_air, _, _, _, _, _, _) = args

    rhoa = air_density(tC_air, pr_air, relhum)

    cpa = 1004.67  # specific heat capacity of (dry) air [J/kg/K]
    hsb = -rhoa * cpa * usr * tsr

    return hsb


"""
#...................................................................................
#...................................................................................
    Simple subroutines used in the routines in the sections above.
        Does NOT include:
            the routines condition_data and make_hourly_data
            warmlayer and coolskin (coare35vn) routines

    air_density
    airtemp_at_refheight
    calc_rain_rate
    gravity
    latent_heat_vaporization_pure_water
    net_longwave_up
    psit_26
    psiu_26
    rain_heat_flux
    rain_heat_flux_FLAWED (DPS code; not used in calculations)
    sea_spechum
    spechum_at_refheight
    water_thermal_expansion
    windspeed_at_refheight
#...................................................................................
#...................................................................................

"""


def air_density(tC_air, air_pressure_mbar, relhum):
    """
    Description

            Returns the density of air.

   Implemented by:

        2014-08-29: Russell Desiderio. Initial Code.

    Usage:

        rho_air = air_density(tC_air, air_pressure_mbar, relhum)

            where

        rho_air = the density of air [kg/m^3]
        tC_air = air temperature (TEMPAIR_L1) [deg_C].
        air_pressure_mbar = barometer pressure (BARPRES_L0) [mbar]
        relhum = relative humidity of air (RELHUMI_L1) [%]
    """
    Rgas = 287.05  # gas constant [J/kg/K] for dry(!) air
    c2k = 273.15  # celsius to kelvin temperature constant
    sp_hum_air = met_spechum(tC_air, air_pressure_mbar, relhum)/1000.0   # units kg/kg
    # the factor of 100 converts pressure from mbar to pascal [N/m^2]
    rho_air = air_pressure_mbar * 100.0 / (1.0 + 0.61 * sp_hum_air) / (Rgas * (tC_air + c2k))

    return rho_air


def airtemp_at_refheight(tC_air, tsr, zrefht, ztmpair, L, lat):
    """
        Formulation from J. Edson's coare35vn version 3.5 coded here is
        equivalent to code in section A.3.3 in the DPS, except that the
        Edson code includes the "lapse" term.
    """
    von = 0.4      # von Karman constant
    cpa = 1004.67  # heat capacity of dry air [J/kg/K]

    lapse = gravity(lat)/cpa
    temp = tC_air + tsr / von * (np.log(zrefht/ztmpair) - psit_26(zrefht/L) +
                                 psit_26(ztmpair/L)) + lapse * (ztmpair-zrefht)

    return temp


def calc_rain_rate(cumulative_precipitation, timestamp):
    """
    Description:

        Calculates rain rate [mm/hr]. This routine does not calculate the formal
        data product RAINRTE, but is used in the routine that does (and in several
        others).

    Implemented by:

        2014-09-19: Russell Desiderio. Initial code.

    Usage:

        rain_rate = calc_rain_rate(cumulative_precipitation, timestamp)

            where

        rain_rate = rain rate [mm/hr]
        cumulative_precipitation = measured rain level (PRECIPM_L1) [mm]
        timestamp = sample date and time value [seconds since 1900-01-01]

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
    # trap out scalar case; return a value of 0
    if cumulative_precipitation.size == 1:
        return 0.0

    # calculate the amount of rain fallen in each time interval
    rainfall = np.diff(cumulative_precipitation)   # [mm]
    # calculate each time interval and convert to hours
    delta_time = np.diff(timestamp) / 3600.0
    # calculate rainrate
    rain_rate = rainfall / delta_time
    # prepend a 0 to represent the first (unknown) rainrate value.
    rain_rate = np.hstack((0.0, rain_rate))
    # follow the DPS and set all negative values to 0.
    #     values of 0 could arise when (1) the sensor resets its rainlevel
    #     to 0 by draining the accumulated rainwater, or when (2) the
    #     evaporation is greater than the precipitation.
    rain_rate[np.less(rain_rate, 0.0)] = 0.0

    return rain_rate


def gravity(lat):
    """
    Description

            Returns acceleration due to gravity as a function of latitude.

    Implemented by:

        2014-06-26: Chris Wingard. Initial Code.
        2014-08-26: Russell Desiderio. Optimized.

    Usage:

        g = gravity(lat)

            where

        g = acceleration due to gravity [m/s/s]
        lat = latitude [degrees]

    from grv.m in the TOGA COARE 3.0 Matlab Toolbox
    """
    gamma = 9.7803267715
    c1 = 0.0052790414
    c2 = 0.0000232718
    c3 = 0.0000001262
    c4 = 0.0000000007

    phi = np.radians(lat)
    x = np.sin(phi)
    xsq = x * x
    g = gamma * (1.0 + xsq * (c1 + xsq * (c2 + xsq * (c3 + xsq * c4))))

    return g


def latent_heat_vaporization_pure_water(tC_water):
    """
    Description

            Returns the latent heat of vaporization of pure water.

   Implemented by:

        2014-08-29: Russell Desiderio. Initial Code.

    Usage:

        Le_water = latent_heat_vaporization_pure_water(tC_water)

            where

        Le_water is the latent heat of vaporization of pure water [J/kg]
        tC_water is the water temperature [degC].

    """
    return (2500.8 - 2.37 * tC_water) * 1000.0


def net_longwave_up(tC_water, total_longwave_down):
    """
    Description:

        Calculates the net upward longwave radiation flux. Note, this routine
        by itself does not calculate the related L2 data product NETLIRR, which
        specifically requires input sea surface skin temperatures corrected for
        warmlayer and coolskin effects.

    Implemented by:

        2014-09-01: Russell Desiderio. Initial code.

    Usage:

        Rnl = net_longwave_up(tC_water, total_longwave_down)

            where

        Rnl = net longwave radiation [W/m^2]; positive in the upward direction
        tC_water = water temperature [deg_C]
        total_longwave_down = measured downward longwave radiation (LONGIRR_L1)
                              [W/m^2]; positive in the downward direction.


        eps(ilon) is the blackbody emissivity; a loose application of Kirchoff's
        thermal radiation law sets the IR reflection coefficient as (1 - eps).

        total longwave radiation down                 = IR
        reflected longwave radiation up               = (1-eps) * IR
        blackbody energy radiated up from sea surface = eps * sigma * Tsea^4

        Rnl(net upward) = IR_up - total_IR_down
                        = eps*sigma*T^4 + (1-eps)*IR - IR
                        = eps * (sigma*T^4 - IR)

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
        """
    sigma = 5.67e-8  # Stefan-Boltzmann constant [W/(m^2 K^4)]
    eps = 0.97
    c2k = 273.15  # degC to kelvin conversion constant
    Rnl = eps * (sigma * (tC_water + c2k) ** 4 - total_longwave_down)

    return Rnl


def psit_26(zet):
    """
    Description

        Computes the temperature structure function and is used to calculate
        air temperature and specific humidity at reference heights above sea
        level.

    Implemented by:

        2014-06-26: Chris Wingard. Initial Code
        2014-09-02: Russell Desiderio. Prevented the possibility of raising
                    a negative number to a fractional power.

    Usage:

        psit = psit_26(zet)

            where

        psit = temperature structure function value at zet
        zet = Monin-Obukhov stability parameter (zu/L) [dimensionless]

    References:

        Edson, J.B., V. Jampana, R.A. Weller, S.P. Bigorre, A.J. Plueddemann,
            C.W. Fairall, S.D. Miller, L. Mahrt, D. Vickers, and H. Hersbach,
            (2013). On the exchange of momentum over the open ocean. Journal of
            Physical Oceanography, 43, 1589-1610.

        Fairall, C.W. (2013). Vectorized version of COARE 3 code with
            modification based on the CLIMODE, MBL and CBLAST experiments.
            ftp://ftp.etl.noaa.gov/users/cfairall/bulkalg/cor3_5/coare35vn.m
    """
    # force the shape of the input to a 1D array.
    zet = np.atleast_1d(zet)

    # stable case: zet > 0.
    # calculate for all zet, and overwrite zet < 0 cases in last section, as
    # done in the original code.
    #
    # because negative zet cases will be overwritten, can take the absolute value
    # of zet underneath the fractional exponent (1.5) without affecting the end
    # result. this is done to avoid the possibility of raising a negative number
    # to a fractional power which will result in complex values.
    dzet = np.minimum(50.0, 0.35 * zet)
    psit = -((1.0 + 0.6667 * np.abs(zet))**1.5 +
             0.6667 * (zet - 14.28) * np.exp(-dzet) + 8.525)

    # overwrite psit for zet < 0 values (unstable case).
    # trap out nans; psit already has nans where we want them
    zet[np.isnan(zet)] = 1.0
    k = zet < 0.0
    x = (1.0 - 16.0 * zet[k])**0.5
    psik = 2.0 * np.log((1.0 + x) / 2.0)
    x = (1.0 - 34.15 * zet[k])**0.3333
    psic = (1.5 * np.log((1.0 + x + x**2) / 3.0) -
            np.sqrt(3) * np.arctan((1.0 + 2.0 * x) / np.sqrt(3)) +
            4.0 * np.arctan(1) / np.sqrt(3))
    f = zet[k]**2 / (1.0 + zet[k]**2)
    psit[k] = (1.0 - f) * psik + f * psic

    return psit


def psiu_26(zet):
    """
    Description

        Computes the velocity structure function and is used to calculate
        wind speed at reference heights above sea level.

    Implemented by:

        2014-09-03: Russell Desiderio. Initial Code

    Usage:

        psiu = psiu_26(zet)

            where

        psiu = velocity structure function value at zet
        zet = Monin-Obukhov stability parameter (zu/L) [dimensionless]

    References:

        Edson, J.B., V. Jampana, R.A. Weller, S.P. Bigorre, A.J. Plueddemann,
            C.W. Fairall, S.D. Miller, L. Mahrt, D. Vickers, and H. Hersbach,
            (2013). On the exchange of momentum over the open ocean. Journal of
            Physical Oceanography, 43, 1589-1610.

        Fairall, C.W. (2013). Vectorized version of COARE 3 code with
            modification based on the CLIMODE, MBL and CBLAST experiments.
            ftp://ftp.etl.noaa.gov/users/cfairall/bulkalg/cor3_5/coare35vn.m
    """
    # force the shape of the input to a 1D array.
    zet = np.atleast_1d(zet)

    # stable case: zet > 0.
    # calculate for all zet, and overwrite zet < 0 cases in last section, as
    # done in the original code.
    dzet = np.minimum(50.0, 0.35 * zet)
    (a, b, c, d) = (0.7, 0.75, 5.0, 0.35)
    psiu = -(a * zet + b * (zet - c / d) * np.exp(-dzet) + b * c / d)

    # overwrite psiu for zet < 0 values (unstable case).
    # trap out nans; psiu already has nans where we want them
    zet[np.isnan(zet)] = 1.0
    k = zet < 0.0
    x = (1.0 - 16.0 * zet[k])**0.25
    psik = (2.0 * np.log((1.0 + x) / 2.0) + np.log((1.0 + x**2) / 2.0) -
            2.0 * np.arctan(x) + 2.0 * np.arctan(1))
    x = (1.0 - 10.15 * zet[k])**0.3333
    psic = (1.5 * np.log((1.0 + x + x**2) / 3.0) -
            np.sqrt(3) * np.arctan((1.0 + 2.0 * x) / np.sqrt(3)) +
            4.0 * np.arctan(1) / np.sqrt(3))
    f = zet[k]**2 / (1.0 + zet[k]**2)
    psiu[k] = (1.0 - f) * psik + f * psic

    return psiu


def rain_heat_flux(rainrate, Tsea, Tair, relhum, pr_air):
    """
    Description:

        Calculates the rain heat flux, which is the heat flux due to rain falling
        into the ocean. Positive heat flux means that the heat is flowing from
        the ocean to the atmosphere (the rain cools the ocean).

        This code represents a new derivation for the calculation of rain heat flux.
        See the documentation for the function rain_heat_flux_FLAWED following,
        for details of the problems with implementation of the old calculation.

    Implemented by:

        2014-10-28: Russell Desiderio. New derivation. Initial code.

    Usage:

        RHF = rain_heat_flux(rainrate, Tsea, Tair, relhum, pr_air)

            where

        RHF = rain heat flux [W/m^2]
        rainrate = rainfall [mm/hr]
        Tsea = bulk sea surface temperature corrected for warmlayer only [degC]
        Tair = air temperature [degC]
        relhum = relative humidity [%]
        pr_air = air pressure [mb]

    References:

        Gosnell, R., C.W. Fairall, and P.J. Webster. The sensible heat of rainfall in
        the tropical ocean. JGR (100)(C9) pp 18437-18442. 1995.

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)

    Derivation:

        This derivation follows that of the Gosnell et al reference cited above. The
        difference is that the wetbulb temperature (assumed to be that of the raindrop)
        is approximated by a Taylor series expansion at the air temperature, as suggested
        by Simon de Szoeke, instead of at the sea temperature, as done in the reference.
        This will be a more accurate calculation of the rain heat flux for normal (rainfall)
        conditions; as Simon pointed out, when Twetbulb < Tair < Tsea, expanding about Tair
        rather than Tsea will give a value closer to the actual wetbulb temperature.

        Equation numbers are from Gosnell et al, cited above. Note that their equations
        (9) and (13) are missing factors of the ratio of the diffusivities.


        The sensible rain heat flux is given by (eqn (5)):

            RHF = cp_rain * R * (Tsea - Train) * rho_rain

                RHF = rain heat flux [W/m^2]
                cp_rain = specific heat of fresh water [J/kg/K]
                R  = rain rate [mm/hr]
                Tsea = bulk sea temperature [degC]
                Train = raindrop temperature, approximated by the wetbulb temperature [degC]
                rho_rain = density of rainwater [kg/m^3] (needed to get RHF in units of W/m^2)

                Units check:

                RHF = J/kg/K * mm/hr * K * kg/m^3 = J/(3600*sec) /(1000*m^2)
                    = W/m^2/3600/1000  (rho_rain=1000, which will cancel)

        To find Train = Twet, use the psychrometric equation (eqn (6)):

            Train = Tair - psi * (qsat_air(Train) - qair)

                Tair = air temperature
                qsat_air = saturation humidity of air as f(temperature)
                qair = ambient specific humidity of air
                psi = (Lv/cp_air)*(dv/dh)
                    Lv = latent heat of vaporization of fresh water at Tair
                    cp_air = specific heat of air
                    dv = water vapor diffusivity
                    dh = heat diffusivity

        Expand qsat_air(Train) as a Taylor series expansion at T = Tair:

            qsat_air(Train) = qsat_air(Tair) + qsat_air'(Tair) * (Train - Tair)

        All terms of the expansion are known; the qsat_air derivative qsat_air'
        is the Clausius-Clapeyron (CC) equation to be evaluated at T=Tair.

        Solve for Train:

            Train = Tair - psi * (qsat_air(Tair) - qair) / (1 + psi * CC(Tair))

        Use this expression for Train in eqn (5) above to calculate RHF.

    """
    c2k = 273.15       # celsius to kelvin temperature constant
    Rgas_air = 287.05  # gas constant [J/kg/K] for air
    Rgas_wtr = Rgas_air / 0.622  # gas constant for water wapor [J/kg/K]
    cp_air = 1004.67   # specific heat capacity of dry air [J/kg/K]
    cp_rain = 4186.0   # specific heat capacity of freshwater at T=25 degC [J/kg/K]
    rho_rain = 1000.0  # density of freshwater [kg/m^3]

    rho_air = air_density(Tair, pr_air, relhum)               # units kg/m^3
    Lv_at_Tair = latent_heat_vaporization_pure_water(Tair)    # units J/kg
    # specific humidity of air
    qair = met_spechum(Tair, pr_air, relhum)/1000.0           # units kg/kg
    # saturation humidity of air at a temperature of Tair
    qsat_at_Tair = met_spechum(Tair, pr_air, 100.0)/1000.0    # units kg/kg

    # diffusivity expressions are taken from original code;
    # these equations were not checked.
    vapor_diffusivity = 2.11e-5 * ((Tair + c2k) / c2k) ** 1.94         # water vapour diffusivity
    heat_diffusivity = (1.0 + 3.309e-3 * Tair -
                        1.44e-6 * Tair * Tair) * 0.02411 / (rho_air * cp_air)  # heat diffusivity

    psi = Lv_at_Tair / cp_air * vapor_diffusivity / heat_diffusivity

    # Clausius-Clapeyron; the factor of 0.622 is included in Rgas_wtr
    CC_at_Tair = Lv_at_Tair * qsat_at_Tair / (Rgas_wtr * (Tair + c2k)**2)

    # rain temperature is assumed to be at the wet bulb temperature.
    Train = Tair - psi * (qsat_at_Tair - qair) / (1.0 + psi * CC_at_Tair)

    # most code versions omit rho_rain and the factor of 1/1000.
    RHF = rainrate * rho_rain * cp_rain * (Tsea - Train) / 1000.0 / 3600.0

    return RHF


def rain_heat_flux_FLAWED(rain_rate, tC_sea, tC_air, relhum, pr_air, dter, dqer, dsea):
    """
        2014-Oct-28.
        THIS CALCULATION OF RAIN HEAT FLUX IS INCORRECT AND NOT USED IN THIS MODULE.

            (1) The original derivation (Gosnell, Fairall, and Webster, JGR(100)(C9)
                18437-18442 1995) dropped a factor of the ratio of the heat and vapor
                diffusivities from eqns (9) and (13). (Verified, personal communication,
                C. Fairall, 2014). This paper's equations were used as they appeared in
                the text when the original coare algorithms were coded so that this error
                has propagated through all versions of the code up to this point.

            (2) The original derivation used a Taylor series expansion at Tsea. Therefore
                the whoi expression for the Clausius-Clapeyron equation is more correct
                than the Edson expression.

            (3) The expansion involves the saturation humidity Qsat of air, the environment
                through which the raindrops are falling, not the saturation humidity of the
                'sea' (of the air directly above the sea surface). Occurrences of the variables
                Qsea and dqer in the original code are incorrect. (If the original derivation
                is to be followed, Qsea should be replaced by Qsat of air at 100% relative
                humidity evaluated at the temperature T=Tsea; dqer would not appear at all).

            (4) The versions of the old code also used the specific heat of seawater instead
                of fresh rain water (almost 5% error).

            (5) The new derivation coded in a separate function follows the suggestion made
                by Simon de Szoeke and features a Taylor series expansion at the air
                temperature instead of at the sea temperature. This is more accurate and
                as an added benefit will avoid the confusion that arose when the sea temperature
                was used as the expansion temperature.

        Incorrectly calculates net upward rain heat flux.

    """
    Rgas = 287.05  # gas constant [J/kg/K] for dry(!) air
    c2k = 273.15   # celsius to kelvin temperature constant
    cpa = 1004.67  # specific heat capacity of (dry) air [J/kg/K]
    cpw = 4000.0   # specific heat capacity of sw at T=20 degC, S=35 [J/kg/K]

    rhoa = air_density(tC_air, pr_air, relhum)               # units kg/m^3
    Le = latent_heat_vaporization_pure_water(tC_sea + dsea)  # units J/kg
    Qair = met_spechum(tC_air, pr_air, relhum)/1000.0        # units kg/kg
    Qsea = sea_spechum(tC_sea + dsea, pr_air)/1000.0         # units kg/kg

    dwat = 2.11e-5 * ((tC_air + c2k) / c2k) ** 1.94        # water vapour diffusivity

    dtmp = (1.0 + 3.309e-3 * tC_air -
            1.44e-6 * tC_air * tC_air) * 0.02411 / (rhoa * cpa)  # heat diffusivity

    # in Clausius-Clayperon eqn, whoi (DPS) and fortran uses tC_sea (+ dsea);
    # jim edson and archived whoi rain_flux.m use tC_air;
    # some versions include dter in this expression

    # this expression affects *most* of the unit tests!
    # this is because rain_heat_flux is used in the warmlayer algorithm, so
    # any data product which uses warmmlayer will fail if a switch is made.

    #dqs_dt = Qair * Le / (Rgas * (tC_air + c2k) ** 2)      # Clausius-Clapeyron
    dqs_dt = Qsea * Le / (Rgas * (tC_sea + dsea + c2k) ** 2)      # Clausius-Clapeyron

    alfac = 1.0 / (1.0 + 0.622 * (dqs_dt * Le * dwat) / (cpa * dtmp))  # wet bulb factor

    factor = rain_rate * alfac * cpw / 3600.0
    rain_heat_flux = factor * ((tC_sea + dsea - tC_air - dter) +
                               (Qsea - Qair - dqer) * Le / cpa)

    return rain_heat_flux


def sea_spechum(tC_sea, p_air):
    """
    Description:

        Calculates the surface value of specific humidity Qsea. Not to be confused
        with any of the specific humidity in air variables (for example, Qair). In
        the original edson coare35vn version 3.5 matlab code, the equivalent function
        is named qsat26sea; in the v3.0 DPS code, qsee.

    Implemented by:

        2014-09-01: Russell Desiderio. Initial code.

    Usage:

        q_sea = sea_spechum(tC_sea, p_air)

            where

        q_sea = sea surface specific humidity [g/kg]
        tC_sea = seawater temperature [deg_C]
        p_air = air pressure (BARPRES_L0) [mbar]

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
    # calculate saturated vapor pressure es in mbar
    es = 6.1121 * np.exp(17.502 * tC_sea/(tC_sea + 240.97)) * (1.0007 + 3.46e-6 * p_air)
    # the factor of 0.98 arises from the fact that the saturation vapor pressure above a
    # salt solution is less than that for pure water; Raoult's law (applicable for the
    # case of an ideal solution) specifies this factor as the mole fraction of water,
    # which for a salinity of 30 is 0.99. Empirically the correct value to use is 0.98.
    esat_sea = 0.98 * es
    # specific humidity q_sea [g/kg] is then:
    q_sea = 621.97 * esat_sea/(p_air - 0.378 * esat_sea)
    return q_sea


def spechum_at_refheight(tC_air, pr_air, relhum, qsr, zrefht, zhumair, L):
    """
        Formulation from J. Edson's coare35vn version 3.5 coded here is
        equivalent to code in section A.3.3 in the DPS.
    """

    # qsr coming in to routine has units of kg/kg
    qsr = qsr * 1000.0  # change units to g/kg for this calculation
    von = 0.4      # von Karman constant
    Q_air = met_spechum(tC_air, pr_air, relhum)
    spechum = Q_air + qsr / von * (np.log(zrefht/zhumair) - psit_26(zrefht/L) +
                                   psit_26(zhumair/L))

    return spechum


def water_thermal_expansion(tC_water):
    """
    Description

            Returns the thermal expansion of water. Used in both the warmlayer
            and coolskin algorithms.

   Implemented by:

        2014-08-29: Russell Desiderio. Initial Code.

    Usage:

        Al = water_thermal_expansion(tC_water)

            where

        Al = water thermal expansion coefficient; no documentation as to units.
        tC_water = water temperature [degC].
    """

    return 2.1e-5 * (tC_water + 3.2)**0.79


def windspeed_at_refheight(wnd, usr, ut, zrefht, zwindsp, L):
    """
        Formulation from J. Edson's coare35vn version 3.5 coded here is
        equivalent to code in section A.3.3 in the DPS.
    """
    von = 0.4      # von Karman constant

    windspeed = wnd + usr / von / ut * wnd * (np.log(zrefht/zwindsp) -
                                              psiu_26(zrefht/L) +
                                              psiu_26(zwindsp/L))

    return windspeed


"""
#...................................................................................
#...................................................................................

    Wrapper function which calls the warmlayer and coolskin (coare35vn) routines:

        seasurface_skintemp_correct

#...................................................................................
#...................................................................................
"""


def seasurface_skintemp_correct(*args):
    """
    Description:

        Wrapper function which by OOI default applies both of the METBK seasurface
        skin temperature correction algorithms (warmlayer, coolskin in coare35vn).
        This behavior is set by the global switches JWARMFL=1 and JCOOLFL=1. The
        switch construction is retained for generality.

        Most of the METBK L2 data products and 2 of the metadata products require
        the skin corrections to be applied before their values can be calculated.

        Warmlayer corrections dsea are added.
        Coolskin corrections dter and dqer are subtracted.

    Implemented by:

        2014-09-01: Russell Desiderio. Initial code.

    Usage (command line spaced out for clarity):

        (usr, tsr, qsr, ut, dter, dqer, tkt, L, zou, zot, zoq,     # coare35vn output
        dt_wrm, tk_pwp, dsea) =                                    # warmlayer output

        seasurface_skintemp_correct

        (rain_rate, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm)


            where

        OUTPUTS (documentation from coare35vn matlab code):

            usr = friction veclocity that includes gustiness [m/s]
            tsr = temperature scaling parameter [K]
            qsr = specific humidity scaling parameter [g/g, I changed this from Edson code]
            ut = not an output of the original code
            dter = coolskin temperature depression [degC]
            dqer = coolskin humidity depression [kg/kg]
            tkt = coolskin thickness [m]
            L = Obukhov length scale [m]
            zou = wind roughness length [m]
            zot = thermal roughness length [m]
            zoq = moisture roughness length [m]

        OUTPUTS (documentation from coare35vnWarm matlab code):

            dt_wrm = warming across entire warmlayer [degC]
            tk_pwp = warmlayer thickness [m]
            dsea = additive warmlayer temperature correction [degC];
                   (this is warmlayer's key output)

        INPUTS:

            rain_rate = rainfall [mm/hr]
            timestamp = seconds since 01-01-1900
            lon = longitude [deg]
            ztmpwat = depth of bulk sea temperature measurement [m]
            tC_sea = bulk sea surface temperature [degC]
            wnd = windspeed relative to current [m/s]
            zwindsp = height of windspeed measurement[m]
            tC_air = air temperature [degC]
            ztmpair = height of air temperature measurement [m]
            relhum = relative humidity [%]
            zhumair = height of air humidity measurement [m]
            pr_air = air pressure [mb]
            Rshort_down = downwelling shortwave irradiation [W/m^2]
            Rlong_down = downwelling longwave irradiation [W/m^2]
            lat = latitude [deg]
            zinvpbl = inversion height; default is 600m [m]
            jcool = switch to activate coolskin algorithm (hardwired to 1 = true)
            jwarm = switch to activate warmlayer algoritgm (hardwired to 1 = true)

    References:

        Fairall, C.W., E.F. Bradley, J.S. Godfrey, G.A. Wick, J.B. Edson, and G.S. Young
        (1996) Cool-skin and warm-layer effects on sea surface temperature. JGR, Vol. 101,
        No. C1, 1295-1308, 1996.

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)

        OOI (2014). 1341-00370_BULKFLX Artifacts. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> REFERENCE >> Data Product Specification Artifacts
            >> 1341-00370_BULKFLX  (Original matlab code).

    Notes:

        (1) the jwarm switch selects whether or not the warmlayer code is run.
            the jcool 'switch' is itself a variable within the (original)
            coare35vn code; it was used as a multiplicative factor when
            calculating coolskin corrections, so that when jcool=0, the
            corrections are set to 0.
        (2) for OOI jwarm and jcool are always 1, because all of the OOI
            sea temperature measurements are bulk, not skin, measurements.
        (3) in the more general case, jwarm = jcool always, because:
            (a) jcool = 1 indicates that the input sea temperature values are
                bulk measurements, not surface skin measurements made with
                an infrared thermometer. in this bulk measurement case, both
                coolskin and warmlayer corrections to the bulk temperature are
                required to model the skin temperature (jcool = jwarm = 1).
            (b) jcool = 0 indicates that the input sea temperature values are
                surface skin temperatures directly measured with an infrared
                thermometer, and therefore both the coolskin and warmlayer
                corrections are not to be applied (jcool = jwarm = 0).
        (4) however, both switches are retained for generality in case this
            open source code is appropriated and adapted. (plus, the DPS
            specified archiving the jwarm and jcool switches as metadata).
        (5) the OOI cyberinfrastructure model originally required that each data
            product be specifically calculated by one function. This is the main
            reason that the wrapper function construct is used. In addition, I've
            chosen to explicitly write out its output tuple arguments for each
            data product call, so that the dependence of the various data products
            on these tuple arguments is obvious (underscores are used as placeholders
            for those arguments not used in any particular function call). In
            particular, this construct specifically identifies when coolskin and
            warmlayer temperature corrections have been applied to various variables
            in the original code. (For example - the latent heat of vaporization for
            water depends on water temperature, but only the warmlayer correction is
            used calculate it).
    """
    jwarm = args[-1]    # jwarm (and jcool) are scalars
    if jwarm:
        (dt_wrm, tk_pwp, dsea) = warmlayer(*args[0:-1])  # does not pass jwarm
    else:
        # the tk_pwp parameter is often used as a divisor in warmlayer calculations to
        # compare the warmlayer depth with the depth of the bulk temperature sensor.
        # when the warmlayer code is not run, the desired results will be obtained if
        # dt_warm and dsea are set to 0 where tk_pwp is nonzero so that a divide by
        # zero error does not result. the value chosen is the default value specified
        # in the warmlayer code itself.
        (dt_wrm, tk_pwp, dsea) = (0.0, 19.0, 0.0)

    # construct tuple containing coolskin input arguments;
    # add the warmlayer temperature correction to the msrd bulk sea temp.
    coolskin_args = (args[4]+dsea,) + args[5:-1]    # does not pass jwarm
    # append results of warmlayer calculation to output,
    # as is also done in original coare35vn warmlayer matlab code.
    return coare35vn(*coolskin_args) + (dt_wrm, tk_pwp, dsea)

"""
#...................................................................................
#...................................................................................

    warmlayer

#...................................................................................
#...................................................................................
"""


def warmlayer(rain_rate, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp, tC_air, ztmpair, relhum,
              zhumair, pr_air, Rshort_down, Rlong_down, lat, zinvpbl, jcool):
    """
    Description:

        Accurate parameterization of air-sea fluxes requires an accurate estimation of the
        interfacial sea temperature. The warmlayer algorithm accounts for solar heating to
        correct sea temperature readings, made by a temperature sensor below the sea surface,
        to the air-sea interface.

        Warmlayer code refactored from coare35vnWarm.m (see DPS artifacts in References).

    Implemented by:

        2014-09-01: Russell Desiderio. Initial code.
        2015-07-08: Russell Desiderio. Added array subscripts to sensor height arrays ztmpwat,
                    zwindsp, ztmpair, and zhumair so that now these parameters on input can
                    either be 1-element 1D arrays or time-vectorized 1D arrays.

    Usage :

        (dt_wrm, tk_pwp, dsea) = warmlayer(rain_rate, timestamp, lon, ztmpwat, tC_sea,
                                           wnd, zwindsp, tC_air, ztmpair, relhum,
                                           zhumair, pr_air, Rshort_down, Rlong_down,
                                           lat, zinvpbl, jcool)


            where

        OUTPUTS (documentation from coare35vnWarm matlab code):

            dt_wrm = warming across entire warmlayer [degC]
            tk_pwp = warmlayer thickness [m]
            dsea = additive warmlayer temperature correction [degC];
                   (this is warmlayer's key output)

        INPUTS:

            rain_rate = rainfall [mm/hr]
            timestamp = seconds since 01-01-1900; hourly
            lon = longitude [deg]
            ztmpwat = depth of bulk sea temperature measurement [m]
            tC_sea = bulk sea surface temperature [degC]
            wnd = windspeed relative to current [m/s]
            zwindsp = height of windspeed measurement[m]
            tC_air = air temperature [degC]
            ztmpair = height of air temperature measurement [m]
            relhum = relative humidity [%]
            zhumair = height of air humidity measurement [m]
            pr_air = air pressure [mb]
            Rshort_down = downwelling shortwave irradiation [W/m^2]
            Rlong_down = downwelling longwave irradiation [W/m^2]
            lat = latitude [deg]
            zinvpbl = inversion height; default is 600m [m]
            jcool = switch to activate coolskin algorithm (hardwired to 1 = true)

    References:

        Fairall, C.W., E.F. Bradley, J.S. Godfrey, G.A. Wick, J.B. Edson, and G.S. Young
        (1996) Cool-skin and warm-layer effects on sea surface temperature. JGR, Vol. 101,
        No. C1, 1295-1308, 1996.

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)

        OOI (2014). 1341-00370_BULKFLX Artifacts. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> REFERENCE >> Data Product Specification Artifacts
            >> 1341-00370_BULKFLX  (Original matlab code).

    Notes:

        Original code was refactored to produce this version. In-code documentation is a
        mixture of the original code's documentation and my additions.

        Two noteworthy changes made: local time calculation from longitude no longer can give
        negative local times; a mystery addend of 7.5 in the local time calculation, resulting
        in increasing local time by half an hour, was deleted.

    """
    # set constants
    c2k = 273.15        # Converts degC to Kelvin
    cpw = 4000.0        # specific heat capacity of sw at T=20 degC, S=35 [J/kg/K]
    rhow = 1022.0       # density of seawater at T=20C, S=31 kg/m^3.
    cpa = 1004.67       # specific heat capacity of (dry) air [J/kg/K]

    #.. hardcoded warmlayer parameters
    rich = 0.65         # critical Richardson number

    #.. initialize warmlayer variables
    fxp = 0.5           # initial value of solar flux absorption
    max_pwp = 19.0      # maximum depth of warm layer (adjustable)
    jamset = 0          # warmlayer threshold indicator
    qcol_ac = 0.0       # accumulates heat from integral
    tau_ac = 0.0        # accumulates stress from integral

    #.. vector calculation of variables used in loop.
    rhoa = air_density(tC_air, pr_air, relhum)
    Rns = met_netsirr(Rshort_down)
    Al = water_thermal_expansion(tC_sea)   # original code does not use dter nor dsea
    grav = gravity(lat)
    ctd1 = np.sqrt(2.0 * rich * cpw / (Al * grav * rhow))         # mess-o-constants 1
    ctd2 = np.sqrt(2.0 * Al * grav / (rich * rhow)) / (cpw**1.5)  # mess-o-constants 2

    #**********************************************************
    nx = timestamp.size        # number of lines of data

    #.. initialize warmlayer products with the default values applicable to
    #.. the case where there is no warmlayer correction to the coare35vn
    #.. coolskin calculation.
    dt_wrm = np.zeros(nx)            # warming across entire warm layer deg.C
    tk_pwp = np.zeros(nx) + max_pwp  # warm layer thickness m
    dsea = np.zeros(nx)              # correction to get to interface

    # local solar time adjustment is a function of longitude:
    #..  360 degrees = 24 *3600 seconds,
    #..  so each degree is worth 240 seconds of time.
    # the OOI timestamp is seconds since midnight 01-jan-1900; therefore
    # local time will still be positive for the case of lon = -180deg.
    local_date_time = timestamp + lon * 240.0

    # and calculate all the delta times [sec] for the integrals' abscissae.
    #.. prepend a zero to line up delta_time with iteration number.
    #.. values at newday records are not used in the calculations.
    delta_time = np.hstack((0.0, np.diff(local_date_time)))

    # determine:
    #    idx_warm: indices of data for days that have data before 6AM
    #    newday_bool: boolean mask, true for the first record of each day.
    #    nanmask: boolean, true for records of days that do not have data before 6AM.
    #             the output at these indices will be changed from initialized to nan.
    idx_warm, newday_bool, nanmask = warmlayer_time_keys(local_date_time)

    #.. the original code has been changed to show the explicit dependence
    #.. of the variables upon iteration count (data record number).
    for ii in idx_warm:   # step through each timepoint

        #.. warmlayer values for the following case are just the initialized
        #.. values. so, instead of using if-then-else, simplify indentation
        #.. by using 'if' only, reset variables, and jump to next iteration.
        if newday_bool[ii]:  # re-zero when starting a new day
            # dt_wrm[ii] = 0.0;
            # tk_pwp[ii] = max_pwp;
            # dsea[ii]   = 0.0;
            jamset = 0
            fxp = 0.5
            tau_ac = 0.0
            qcol_ac = 0.0
            continue  # go to next time (data) record
            # end midnight reset

        #*****  dependent variables for the [ii]th warm layer calculation
        #*****  of dsea are fluxes, coolskin correction dter, and dsea itself,
        #*****  which are derived from the previous ([ii-1]th) data record.
        #
        # because of the dependence on the previous value of dsea, this calculation
        # cannot be vectorized.
        tsea_corr = tC_sea[ii-1] + dsea[ii-1]

        # slicing 1D arrays with [ii-1:ii] returns a 1-element nd.array variable which
        # can be indexed, whereas slicing with [ii-1] returns a variable which cannot be
        # indexed. [ii-1:ii] slicing is used so that coare35vn can be run with both
        # 'scalar' and 'vector' input.
        args = (tsea_corr, wnd[ii-1:ii], zwindsp[ii-1:ii], tC_air[ii-1:ii], ztmpair[ii-1:ii], relhum[ii-1:ii], zhumair[ii-1:ii],
                pr_air[ii-1:ii], Rshort_down[ii-1:ii], Rlong_down[ii-1:ii], lat[ii-1:ii], zinvpbl[ii-1:ii],
                jcool)

        (usr, tsr, qsr, ut, dter, dqer, _, _, _, _, _) = coare35vn(*args)

        # in the original matlab code, Le was calculated inside of the coare35vn
        # subroutine, which was called using tC_sea+dsea for seawater temperature:
        Le = latent_heat_vaporization_pure_water(tsea_corr)
        tau_old = rhoa[ii-1] * usr * usr * wnd[ii-1] / ut  # stress
        hs_old = -rhoa[ii-1] * cpa * usr * tsr              # sensible heat flux
        hl_old = -rhoa[ii-1] * Le * usr * qsr                      # latent heat flux

        # note:
        #     the source v3.5 matlab code is followed here: it does not use dsea
        #     in the Rnl expression used in the warmlayer calculation, although
        #     dsea is used in the expression for RF_old.
        Rnl = net_longwave_up(tC_sea[ii]-dter, Rlong_down[ii])
        RF_old = rain_heat_flux(rain_rate[ii-1], tC_sea[ii-1]+dsea[ii-1], tC_air[ii-1],
                                relhum[ii-1], pr_air[ii-1])

        #********************************************************
        #****  Compute warm layer correction *******************
        #********************************************************
        qr_out = Rnl + hs_old + hl_old + RF_old  # total cooling at surface
        q_pwp = fxp * Rns[ii] - qr_out          # tot heat abs in warm layer

        # calculate dt_wrm and tk_pwp for this iteration.
        if q_pwp >= 50.0 or jamset == 1:         # Check for threshold
            jamset = 1			         # indicates threshold crossed
            tau_ac = tau_ac + np.maximum(.002, tau_old) * delta_time[ii]  # momentum integral

            # check threshold for warm layer existence
            if qcol_ac + q_pwp * delta_time[ii] > 0.0:
                #******************************************
                # Compute the absorption profile
                #******************************************
                #.. tk_pwp can iteratively change value in the following loop,
                #.. requiring the creation of the variable tkpwp.
                tkpwp = tk_pwp[ii-1]
                for i in range(5):               # loop 5 times for fxp
                    fxp = 1.0 - (0.28 * 0.014 * (1.0 - np.exp(-tkpwp / 0.014)) +
                                 0.27 * 0.357 * (1.0 - np.exp(-tkpwp / 0.357)) +
                                 0.45 * 12.82 * (1.0 - np.exp(-tkpwp / 12.82))) / tkpwp
                    qjoule = (fxp * Rns[ii] - qr_out) * delta_time[ii]
                    if qcol_ac + qjoule > 0.0:   # Compute warm-layer depth
                        tkpwp = np.minimum(max_pwp,
                                           ctd1[ii] * tau_ac / np.sqrt(qcol_ac + qjoule))
                tk_pwp[ii] = tkpwp
            else:                                # warm layer wiped out
                fxp = 0.75
                tk_pwp[ii] = max_pwp
                qjoule = (fxp * Rns[ii] - qr_out) * delta_time[ii]

            qcol_ac = qcol_ac + qjoule           # heat integral

            #*******  compute dt_warm  ******
            if qcol_ac > 0.0:
                dt_wrm[ii] = ctd2[ii] * (qcol_ac)**1.5 / tau_ac
            else:
                dt_wrm[ii] = 0.0

        else:   # propagate dt_wrm and tk_pwp values
            dt_wrm[ii] = dt_wrm[ii-1]
            tk_pwp[ii] = tk_pwp[ii-1]

        # Compute warm layer correction dsea
        if tk_pwp[ii] < ztmpwat[ii]:
            dsea[ii] = dt_wrm[ii]
        else:
            dsea[ii] = dt_wrm[ii] * ztmpwat[ii] / tk_pwp[ii]

    # for all days that did not begin before 6AM, return NaNs
    dt_wrm[nanmask] = np.nan
    tk_pwp[nanmask] = np.nan
    dsea[nanmask] = np.nan

    return dt_wrm, tk_pwp, dsea

"""
#...................................................................................
#...................................................................................

    coare35vn coolskin code and the subroutines unique to it.
        coare35vn also calls subroutines located elsewhere in this module.

    coare35vn

        charnock_wind
        coolskin_parameters
        effective_relwind
        obukhov_for_init
        obukhov_length_scale
        roughness_lengths
        roughness_lengths_for_init
        scaling_parameters
#...................................................................................
#...................................................................................
"""


def coare35vn(tC_sea, wnd, zwindsp, tC_air, ztmpair, relhum, zhumair, pr_air,
              Rshort_down, Rlong_down, lat, zinvpbl, jcool):
    """
    Description:

        Iteratively calculates fundamental bulk parameters, with coolskin correction,
        from which air-sea fluxes can be calculated. Transliterated from version 3.5
        of coare35vn.m (see DPS artifacts in References). In contrast to the original
        fortran and subsequent matlab versions, this python version does not directly
        calculate data products.

        Accurate parameterization of air-sea fluxes requires an accurate estimation of
        the interfacial sea temperature. The coolskin algorithm accounts for interfacial
        cooling due to longwave blackbody radiation, sensible heat flux, and latent heat
        flux.

    Implemented by:

        2014-09-01: Russell Desiderio. Initial code.

    Usage (command line spaced out for clarity):

        (usr, tsr, qsr, ut, dter, dqer, tkt, L, zou, zot, zoq) = coare35vn(
                            tC_sea, wnd, zwindsp, tC_air, ztmpair, relhum, zhumair,
                            pr_air, Rshort_down, Rlong_down, lat, zinvpbl, jcool)

            where

        OUTPUTS (documentation from coare35vn matlab code):

            usr = friction veclocity that includes gustiness [m/s]
            tsr = temperature scaling parameter [K]
            qsr = specific humidity scaling parameter [g/g, I changed this from Edson code]
            ut = not an output of the original code
            dter = coolskin temperature depression [degC]
            dqer = coolskin humidity depression [kg/kg]
            tkt = coolskin thickness [m]
            L = Obukhov length scale [m]
            zou = wind roughness length [m]
            zot = thermal roughness length [m]
            zoq = moisture roughness length [m]

        INPUTS:

            tC_sea = bulk sea surface temperature [degC]
            wnd = windspeed relative to current [m/s]
            zwindsp = height of windspeed measurement[m]
            tC_air = air temperature [degC]
            ztmpair = height of air temperature measurement [m]
            relhum = relative humidity [%]
            zhumair = height of air humidity measurement [m]
            pr_air = air pressure [mb]
            Rshort_down = downwelling shortwave irradiation [W/m^2]
            Rlong_down = downwelling longwave irradiation [W/m^2]
            lat = latitude [deg]
            zinvpbl = inversion height; default is 600m [m]
            jcool = switch to activate coolskin algorithm (hardwired to 1 = true)

    References:

        Fairall, C.W., E.F. Bradley, J.S. Godfrey, G.A. Wick, J.B. Edson, and G.S. Young
        (1996) Cool-skin and warm-layer effects on sea surface temperature. JGR, Vol. 101,
        No. C1, 1295-1308, 1996.

        Fairall, C. W., E.F. Bradley, D.P. Rogers, J.B. Edson, and G.S. Young (1996),
        Bulk Parameterization of air-sea fluxes for Tropical Ocean-Global Atmosphere
        Coupled-Ocean Atmosphere Response Experiment. JGR Vol. 101 No. C2, 3747-3764.

        Fairall, C.W., E.F. Bradley, J.E. Hare, A.A. Grachev, and J.B. Edson (2003),
        Bulk parameterization of air sea fluxes: updates and verification for the
        COARE algorithm. J. Climate, 16, 571-590.

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)

        OOI (2014). 1341-00370_BULKFLX Artifacts. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> REFERENCE >> Data Product Specification Artifacts
            >> 1341-00370_BULKFLX  (Original matlab code).

    Notes:

        The code within the iteration loop was segregated into subroutines to clarify which
        of the variables were key in that they were recalculated during each iteration, and
        afterwards used to calculate data products. These subroutine names are somewhat
        arbitrary in that there was essentially no documentation in the original code for
        what was being calculated in the loop.

        These key variables are the output of this stripped down version of the original code.
        This python code is not meant to directly calculate data products; rather, it is to
        calculate the key variables from which all associated data products can be calculated.

        This approach greatly increases the number of subroutines required for the calculations.
        However, it also makes explicit which parameters are necessary for the calculation of
        each data product, and in some cases may expose an inconsistent application of either
        the coolskin or warmlayer corrections to the bulk seasurface temperature when these
        data products are calculated.
    """
    # convert relative humidity to specific humidity [kg/kg]
    Qsea = sea_spechum(tC_sea, pr_air) / 1000.0          # surface water specific humidity
    Qair = met_spechum(tC_air, pr_air, relhum) / 1000.0  # specific humidity of air

    #***********  set constants **********************************************
    Beta = 1.2
    von = 0.4     # von karman constant
    fdg = 1.00    # Turbulent Prandtl number
    c2k = 273.15  # degC to kelvin
    grav = gravity(lat)

    #***********  air constants **********************************************
    Rgas = 287.05
    Le = latent_heat_vaporization_pure_water(tC_sea)
    cpa = 1004.67   # specific heat capacity of dry air, J/kg/K
    rhoa = air_density(tC_air, pr_air, relhum)
    visa = 1.326e-5 * (1.0 + 6.542e-3 * tC_air + 8.301e-6 * tC_air**2 -
                       4.84e-9 * tC_air**3)
    lapse = grav / cpa

    #***********  cool skin constants  ***************************************
    Al = water_thermal_expansion(tC_sea)
    be = 0.026     # salinity expansion coeff.
    cpw = 4000.0   # specific heat of seawater at T=20C, S=35 J/kg/K.
    rhow = 1022.0  # density of seawater at T=20C, S=31 kg/m^3.
    visw = 1.0e-6
    tcw = 0.6
    bigc = 16.0 * grav * cpw * (rhow * visw)**3 / (tcw**2 * rhoa**2)

    #.. derived variables unchanged by loop
    dq = Qsea - Qair
    dt = tC_sea - tC_air - lapse * ztmpair
    tK_air = tC_air + c2k
    tv = tK_air * (1.0 + 0.61 * Qair)   # virtual temperature
    ug = 0.5

    #*********** initialization **********************************************

    #.. coolskin parameters (changed inside loop and used after loop)
    #
    #.. if jcool is set to 0, then all dter and dqer values in this code will
    #.. also be 0 (and therefore not require a multiplicative factor of jcool
    #.. as in the original matlab code) *except* for when they are directly
    #.. calculated inside the loop.
    dter = 0.3 * jcool
    wetc = 0.622 * Le * Qsea / (Rgas * (tC_sea + c2k)**2)
    dqer = dter * wetc
    tkt = 0.001 + np.zeros(wnd.size)
    ut = np.sqrt(wnd**2 + ug**2)
    # for initialization of usr, tsr, qsr, charnC
    # original code does use (10.0/1e-4)
    u10 = ut * np.log(10.0/1.0e-4) / np.log(zwindsp/1.0e-4)

    #.. scaling parameters usr,qsr,tsr
    zo10, zot10 = roughness_lengths_for_init(0.035*u10, grav, visa, von)
    #.. k50 is used inside the loop to save loop variables for ii=0
    L10, k50 = obukhov_for_init(von, grav, tK_air, dt, dter, dq, ut,
                                zwindsp, ztmpair, zo10, zot10, zinvpbl, Beta)
    usr, qsr, tsr = scaling_parameters(dter, dqer, von, fdg, zwindsp,
                                       zhumair, ztmpair, zo10, zot10, zot10,
                                       L10, ut, dq, dt)

    #***********  net radiation fluxes ***************************************
    Rns = met_netsirr(Rshort_down)                  # net shortwave radiation DOWN
    Rnl = net_longwave_up(tC_sea-dter, Rlong_down)  # net longwave radiation UP

    #**********************************************************
    #  The following gives the new formulation for the
    #  Charnock variable in COARE 3.5
    #**********************************************************

    charnC = charnock_wind(u10)

    nits = 6  # hardwired number of iterations

    #**************  bulk loop ***********************************************
    for ii in range(nits):

        L = obukhov_length_scale(von, grav, tK_air, Qair, usr, tsr, qsr)
        zou, zoq, zot = roughness_lengths(charnC, usr, grav, visa)

        usr, qsr, tsr = scaling_parameters(dter, dqer, von, fdg, zwindsp,
                                           zhumair, ztmpair, zou, zoq, zot, L, ut,
                                           dq, dt)

        dter, dqer, tkt = coolskin_parameters(usr, qsr, tsr, Rnl, Rns, rhoa, cpa,
                                              Le, tkt, Al, be, cpw, visw, rhow,
                                              bigc, tcw, tC_sea, Qsea, pr_air)

        # these coolskin parameters must be reset to 0 if coolskin is off, so:
        dter = dter * jcool
        dqer = dqer * jcool

        Rnl = net_longwave_up(tC_sea-dter, Rlong_down)

        ut = effective_relwind(tsr, tK_air, qsr, grav, tv, usr, Beta, zinvpbl, wnd)

        #.. update charnock variable.
        u10N_for_charnC = usr / von / ut * wnd * np.log(10.0/zou)
        charnC = charnock_wind(u10N_for_charnC)

        # for stable cases designated by k50 save the results from the
        # first iteration as the algorithm output. this construction also
        # works as desired if k50 is empty.
        if ii == 0:
            stable_cases = (usr[k50], tsr[k50], qsr[k50], ut[k50],
                            dter[k50], dqer[k50], tkt[k50], L[k50],
                            zou[k50], zot[k50], zoq[k50])

    # loop is finished:
    # insert first iteration solution for stable cases.
    (usr[k50], tsr[k50], qsr[k50], ut[k50], dter[k50], dqer[k50],
        tkt[k50], L[k50], zou[k50], zot[k50], zoq[k50]) = stable_cases

    # the whoi v3.0 and jbe v3.5 qsr units differ.
    #    whoi (DPS): [kg/kg] (same as [g/g])
    #    jbe - output units are [g/kg]
    #
    #    because the calculations of fluxes are natively done with qsr in units
    #    of [kg/kg], this is what will be used.
    #
    #    HOWEVER - when calculating specific humidity at a reference height,
    #              qsr must have units of g/kg:
    #              qsr = qsr * 1000    # changes units from kg/kg to g/kg

    # note: original DPS code returns nonzero values for dter and dqer even
    # when jcool = 0; the edson v3.5 code returns a nonzero dter value, but
    # a zero value for dqer.
    #
    # this code gives dter = dqer = 0 when jcool = 0. this also obviates the
    # need for multiplicative factors of jcool when this output is used in
    # subsequent routines.

    return (usr, tsr, qsr, ut, dter, dqer, tkt, L, zou, zot, zoq)


#-------------------------------------------------------------------------
#---------------- subroutines unique to coare35vn-------------------------
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
def charnock_wind(u10N):
    """
        section essentially verbatim from original code.
    """
    umax = 19
    a1 = 0.0017
    a2 = -0.0050
    charnC = a1 * u10N + a2
    # trap out nans after nans have propagated into charnC
    u10N[np.isnan(u10N)] = -1.0  # anything less than umax
    charnC[u10N > umax] = a1 * umax + a2
    return charnC


#-------------------------------------------------------------------------
def coolskin_parameters(usr, qsr, tsr, Rnl, Rns, rhoa, cpa, Le, tkt, Al, be,
                        cpw, visw, rhow, bigc, tcw, tC_sea, Qsea, pr_air):
    """
        section essentially verbatim from original code.
    """
    #.. dter: coolskin temperature depression
    #.. dqer: coolskin humidity depression
    #.. tkt: coolskin thickness
    N = usr.shape[0]
    hsb = -rhoa * cpa * usr * tsr
    hlb = -rhoa * Le * usr * qsr
    qout = Rnl + hsb + hlb
    dels = Rns * (0.065 + 11.0 * tkt - 6.6e-5 / tkt *
                  (1.0 - np.exp(-tkt / 8.0e-4)))
    qcol = qout - dels
    alq = Al * qcol + be * hlb * cpw / Le
    xlamx = 6.0 * np.ones(N)
    tkt = np.minimum(0.01, xlamx * visw / (np.sqrt(rhoa / rhow) * usr))
    # trap out nans; tkt already has nans where we want them
    alq[np.isnan(alq)] = -1.0
    k = np.where(alq > 0)
    xlamx[k] = 6.0 / (1.0 + (bigc[k] * alq[k] / usr[k]**4)**0.75)**0.333
    tkt[k] = xlamx[k] * visw / (np.sqrt(rhoa[k] / rhow) * usr[k])
    dter = qcol * tkt / tcw
    # formerly, dqer = wetc * dter
    dqer = Qsea - sea_spechum(tC_sea - dter, pr_air) / 1000.0
    return dter, dqer, tkt


#-------------------------------------------------------------------------
def effective_relwind(tsr, tK_air, qsr, grav, tv, usr, Beta, zinvpbl, wnd):
    """
        section essentially verbatim from original code.
    """
    N = wnd.shape[0]
    tvsr = tsr + 0.61 * tK_air * qsr
    Bf = -grav / tv * usr * tvsr
    ug = 0.2 * np.ones(N)
    nanmask = np.isnan(Bf)
    Bf[nanmask] = -1.0
    ug[Bf > 0] = np.maximum(0.2,
                            Beta * (Bf[Bf > 0] * zinvpbl[Bf > 0]) ** 0.333)
    ug[nanmask] = np.nan
    ut = np.sqrt(wnd**2 + ug**2)
    return ut


#-------------------------------------------------------------------------
def obukhov_for_init(von, grav, tK_air, dt, dter, dq, ut, zwindsp,
                     ztmpair, zo10, zot10, zinvpbl, Beta):
    """
        section essentially verbatim from original code.
    """
    #.. calculates an obukhov length scale for loop variable initializations.
    #.. also finds indices (k50) of zetu stability values with very thin
    #..    Monin-Obukhov lengths relative to zwindsp.
    Ribu = -grav * zwindsp / tK_air * (dt - dter + 0.61 * tK_air * dq) / ut**2
    #.. CC calculation is left unchanged from earlier code.
    Cd = (von / np.log(zwindsp/zo10))**2
    Ct = von / np.log(ztmpair/zot10)
    CC = von * Ct / Cd
    zetu = CC * Ribu * (1.0 + 27.0 / 9.0 * Ribu / CC)
    # trap out nans
    nanmask = np.isnan(zetu)
    zetu[nanmask] = 0.0
    k50 = np.where(zetu > 50)  # stable with very thin M-O length relative to zwindsp
    # restore nans to zetu
    zetu[nanmask] = np.nan

    Ribcu = -zwindsp / zinvpbl / 0.004 / Beta**3
    # trap out nans; Ribu>0 values are not further processed, so:
    Ribu[np.isnan(Ribu)] = 1.0
    k = Ribu < 0
    zetu[k] = CC[k] * Ribu[k] / (1.0 + Ribu[k] / Ribcu[k])
    L_init = zwindsp / zetu
    return L_init, k50


#-------------------------------------------------------------------------
def obukhov_length_scale(von, grav, tK_air, Qair, usr, tsr, qsr):
    """
        DPS documentation:
        TO EVALUATE OBUKHOV LENGTH FROM AVERAGE TEMP T IN (kelvin)
        AVERAGE HUMIDITY Q (in air) IN GM/GM,
        AND FRICTIONAL VEL,TEMP.,HUM. IN MKS UNITS
        SEE LIU ET AL. (1979)
    """
    tv = tK_air * (1.0 + 0.61 * Qair)
    tvsr = tsr * (1.0 + 0.61 * Qair) + 0.61 * tK_air * qsr
    # add tvsr adjustment to avoid program failure when tvsr very small
    #.. (and evidently, negative).
    #.. note also that tvsr=0 is not trapped out in original code.
    nanmask = np.logical_or(np.isnan(tvsr), tvsr == 0)
    tvsr[nanmask] = 100.0
    mask = np.abs(tvsr) < 1.e-3
    tvsr[mask] = np.abs(tvsr[mask])
    tvsr[nanmask] = np.nan
    L = tv * usr * usr / (grav * von * tvsr)

    return L


#-------------------------------------------------------------------------
def roughness_lengths(charn, usr, grav, visa):
    """
        section essentially verbatim from original code.
    """
    zo = charn * usr**2 / grav + 0.11 * visa / usr  # surface roughness
    rr = zo * usr / visa
    # These thermal roughness lengths give Stanton and
    # Dalton numbers that closely approximate COARE 3.0
    zoq = np.minimum(1.6e-4, 5.8e-5 / rr**0.72)
    zot = zoq
    return zo, zoq, zot


#-------------------------------------------------------------------------
def roughness_lengths_for_init(usr, grav, visa, von):
    """
        section essentially verbatim from original code.
    """
    charn_standin = 0.011
    zo10 = charn_standin * usr**2 / grav + 0.11 * visa / usr
    #.. the following code is unchanged from the original version,
    #.. in case there's ever a question of where zot10 came from.
    Cd10 = (von / np.log(10./zo10))**2
    Ch10 = 0.00115
    Ct10 = Ch10 / np.sqrt(Cd10)
    zot10 = 10. / np.exp(von/Ct10)
    return zo10, zot10


#-------------------------------------------------------------------------
def scaling_parameters(dter, dqer, von, fdg, zwindsp, zhumair,
                       ztmpair, zou, zoq, zot, L, ut, dq, dt):
    """
        section essentially verbatim from original code.
    """
    cdhf = von / (np.log(zwindsp / zou) - psiu_26(zwindsp / L))
    cqhf = von * fdg / (np.log(zhumair / zoq) - psit_26(zhumair / L))
    cthf = von * fdg / (np.log(ztmpair / zot) - psit_26(ztmpair / L))
    usr = ut * cdhf
    qsr = -(dq - dqer) * cqhf
    tsr = -(dt - dter) * cthf
    return usr, qsr, tsr


"""
#...................................................................................
#...................................................................................
    Data conditioning and averaging routines

        vet_velptmn_data
        condition_data
        make_hourly_data
        warmlayer_time_keys
#...................................................................................
#...................................................................................
"""


def vet_velptmn_data(vle, vln, use_velptmn):
    """
    Description:

        When an element of the use_velptmn switch is 0, the corresponding elements in
        vle_in and vln_in are set to Nan.

        This function is not to be used with met_relwind_speed, where suspect velptmn
        values are replaced with 0.

    Implemented by:

        2015-07-10: Russell Desiderio. Initial Code

    Usage:

        vle_out, vln_out = vet_velptmn_data(vle, vln, use_velptmn)

            where

        vle_out, vln_out = eastward and northward current speeds with suspect values
                           returned as nan.
        vle, vln = input eastward and northward current speeds
        use_velptmn = time-vectorized data quality flag:
                      0 -> bad  current data
                      1 -> good current data
    """
    # expand use_velptmn_with_metbk if it is called as a scalar
    if np.atleast_1d(use_velptmn).shape[0] == 1:
        use_velptmn = np.tile(use_velptmn, vle.shape[0])

    # to prevent what turned out to be very much unanticipated "call-by-reference"-ish
    # ramifications in the unit tests
    vle_out = np.copy(vle)
    vln_out = np.copy(vln)
    # replace aliased current values with nans.
    nanmask = use_velptmn == 0
    vle_out[nanmask] = np.nan
    vln_out[nanmask] = np.nan

    return vle_out, vln_out


def condition_data(*args):
    """
    Description:

        (1) Makes sure that all relevant variables are at least 1D np.arrays.

        (2) For missing input variables with default scalar values set in the
            input argument list, this routine expands the size of those variable
            arrays from 1 element to the size of that of the other variables.

        (3) Also expands (time-vectorizes) sensor heights ztmpwat, zwindsp, ztmpair,
            and zhumair to maintain compatibility with scalar sensor height inputs.

        (4) Conditions the jcool and jwarm input variables; see code documentation.

    Implemented by:

        2014-09-19: Russell Desiderio. Initial code.
        2015-07-08: Russell Desiderio. Added sensor height indices [3, 6, 8, 10]
                    for [ztmpwat, zwindsp, ztmpair, zhumair] to idx_of_args_to_expand
                    so that now these parameters on input can either be 1-element 1D
                    arrays or time-vectorized 1D arrays. Also conditioned the jcool
                    and jwarm variables to convert them into 1-element switches if
                    they are input as time-vectorized variables.

    Usage:

        args_out = condition_data(*args_in)

            where

        args_out = argument list of conditioned output data
        args_in = argument list of input data

    Notes:

        This routine will be called at the front end of each DPA that requires
        the warmlayer/coolskin algorithm. It may also be called in other instances
        (for example, for rainrte) that may not require variable expansion, in
        which case the number of arguments will be less than number_of_bulk_vars.

    """
    # to enable modification of the input arguments in place,
    # make sure that the args are in a list
    args = list(args)

    number_of_bulk_vars = 18

    # zinvpbl [15] must always be expanded
    idx_of_args_to_expand = [0, 3, 6, 8, 10, 11, 12, 13, 14, 15]

    nargs = len(args)
    for ii in range(nargs):
        args[ii] = np.atleast_1d(args[ii])

    # for rainrte and testing
    if nargs < number_of_bulk_vars:
        return args

    # condition the jcool and jwarm switches.
        # (1) these should be type integer, so there is no need to trap out Nans.
        # (2) these switches should not be time-vectorized in the code. therefore,
        #     the code is made compatible with time-vectorized inputs of these
        #     switches by using only the 1st value; and, if it is not zero, change
        #     it to 1. This will also trap out -99999999 system fillvalues.
    for ii in [16, 17]:
        if args[ii][0] != 0:
            args[ii][0] = 1
        args[ii] = args[ii][0:1]  # return 1 element as an ndarray

    # expand if only one item in argument
    n_records = args[1].size
    for ii in idx_of_args_to_expand:
        if args[ii].size == 1:
            args[ii] = np.zeros(n_records) + args[ii]

    return args


def make_hourly_data(*args):
    """
    Description:

        Calculates hourly averaged data for the variables passed in the
        argument list based on timestamp. This routine will work for
        sporadically spaced data and for data with time gaps; in the latter
        case, no records are produced for the missing time bins.

        The timestamp for each hour of data is calculated as the midpoint
        of the bin interval. If this is changed to either the beginning or
        end of the interval (the latter may be more appropriate given that
        the warmlayer routine numerically calculates integrals over time)
        then the check values in the unit tests may also need to be changed.

        Note also that the first hourly timestamp will be one half hour after
        the timestamp of the first datum in, so that in the general case the
        hourly timestamps will not correspond to times at the top of the hour
        (ie, the hourly timestamps will not correspond to 6:00, 7:00, 8:00, etc).

    Implemented by:

        2014-09-18: Russell Desiderio. Initial code.
        2014-10-22: Russell Desiderio. Added capability to process timestamps if
                    it is the only input argument. This was added so that a
                    simple call could be made to this function for the purposes
                    of calculating an hourly timestamp metadata product not
                    included in the DPS but which I anticipate will be needed.
        2015-07-08: Russell Desiderio. Deleted sensor height indices [3, 6, 8, 10]
                    for [ztmpwat, zwindsp, ztmpair, zhumair] from idx_to_skip
                    so that now these parameters on input can either be 1-element 1D
                    arrays or time-vectorized 1D arrays.

    Usage

        args_out = make_hourly_data(*args_in)

            where

        args_out = argument list of hourly data
        args_in = argument list of each minute data


        Timestamps in units of seconds must be the second element in args_in, unless
        it is the only input in the argument list; this was dictated by convenience for
        the ordering of the input arguments into the seasurface_skintemp_correct routine.

        If the number of arguments in is 18, then the inputs are data to be
        prepared for running in the warmlayer and coolskin routines, and so not
        all of the input arguments are to be processed (for example, switches
        and the sensor heights). The variables to be processed is determined by
        a hard-wired list argument index number idx. For any other args_in sized
        input, all inputs are processed.

        All arguments must be 1D arrays of the same length.

    Notes:

        The np.bincount routine is used in much the same way accumarray in matlab
        would be used to construct the hourly data. np.histogram could also have
        been used.

        The key to the routine is to convert the timestamps into elapsed hours
        [0.0, 0.1, 0.7, 1.3, 1.7, 3.1, 3.3] so that when these values are floored
        [ 0 ,  0 ,  0 ,  1 ,  1 ,  3 ,  3 ] the entry at an index represents the
        bin number into which the various variables with that index will belong.

        The summing is carried out by using the weighting feature of the np.bincount
        function, as described in the example in the numpy.bincount documentation at:
        http://docs.scipy.org/doc/numpy-1.8.1/reference/generated/numpy.bincount.html.

        The same lines of code are executed when calculating the hourly timestamps
        regardless of the number of input arguments, so that if the current method
        is changed only the code at the end of this routine needs to be modified.
    """
    args = list(args)

    number_of_bulk_vars = 18

    # prep all variables for rainrte (nargs<18).
    # for nargs >= 18, skip arguments that are constants (except for zinvpbl)
    # and switches.
    idx_to_skip = [17, 16]

    nargs = len(args)
    idx = range(nargs)
    if nargs >= number_of_bulk_vars:
        for ii in idx_to_skip:
            del idx[ii]

    # timestamps must be the 2nd variable in the input argument list,
    # unless there is only 1 variable.
    index_timedata = np.sign(nargs-1)
    time_sec = args[index_timedata]
    time_elapsed_hr = (time_sec - time_sec[0])/3600.0

    # assign each timestamp a bin number index based on its elapsed time in hrs.
    bin_number = np.floor(time_elapsed_hr).astype(int)
    # the number of elements in each hourly bin is given by
    bin_count = np.bincount(bin_number).astype(float)
    # create a logical mask of non-zero bin_count values
    mask = (bin_count != 0)
    # and keep only the non-zero values for calculating the average
    bin_count_no_zeros = bin_count[mask]

    # average the values in each hourly bin for each input variable;
    # the np.bincount function only works on 1D arrays, so
    for ivar in idx:
        # sum the values in each hourly bin for the ivar[th] variable
        args[ivar] = np.bincount(bin_number, args[ivar])
        # discard trivial bin sums of 0 where there were no bin elements
        args[ivar] = args[ivar][mask]
        # divide the bin sums by the number of elements in each bin
        args[ivar] = args[ivar] / bin_count_no_zeros

    # hourly timestamp calculation:
    #     note that the midpoint of the data interval is used, not the timestamp
    #     of the first nor last point.
    #
    #     use the midpoint of the bins as the timestamp, instead of the average
    #     of the timestamps within the bin as calculated in the above loop; this
    #     would give significantly different values only if there are missing data.
    bin_time_sec = time_sec[0] + 1800.0 + 3600.0 * np.array(range(len(mask)))
    # delete bins with no entries as before
    bin_time_sec = bin_time_sec[mask]

    args[index_timedata] = bin_time_sec

    #for ii in idx:
    #    print args[ii]

    return args


def warmlayer_time_keys(localdate):
    """
    Description:

        Calculates time variables for the warmlayer routine. See the
        documentation in the Usage and Notes sections below.

    Implemented by:

        2014-10-22: Russell Desiderio. Initial code.

    Usage

        idx_warm, newday, nanmask = warmlayer_time_keys(localdate)

            where

        idx_warm = indices of data records to be processed by the warmlayer routine;
                   these are data for days for which there are data before a threshold
                   time value early in the morning (set to 6AM).
        newday = boolean array: true for the first record of a day, false otherwise.
        nanmask = boolean array: true for indices of data records not to be processed
                  by the warmlayer routine.
        localdate = local (not UTC) date and time [sec since 01-01-1900].

    Notes

        Although this routine was written to process hourly data, it will work with
        data on any time base.

        This function provides a measure of quality control for the automated calculation
        of metbk data products. One case: if there are data during one day from 5AM to 10AM,
        then an absence of data so that the next data immediately following are from the
        next day from 11AM to 4PM, the matlab code will treat these data as all coming from
        the same day. Using warmlayer_time_keys will result in the first day's data products
        being calculated as expected and will result in Nans for the second day's data as
        desired.

        The original fortran and matlab versions of the warmlayer function hard-coded
        a 6 AM threshold (equatorial sunrise). The OOI routines coded here construct
        hourly averages from each-minute data and assign a timestamp at the midpoint of
        the binning interval. So, if a day's first each-minute data record timestamp is
        at 5:45 AM local, then the the first hourly timestamp will be 6:15 AM, in which
        case the warmlayer routine would not be run on that day's data.

        The threshold is set here at 6:00 AM local just as it is in the original code. The
        rationale is that the first day's data should start before sunrise at the beginning
        of the daily heating cycle. At some point it may be desirable to replace this value
        with one calculated as a function of time of year, latitude, and perhaps longitude.

        If the timestamp assigned to the hourly intervals is changed to either the first or
        last time of the binning interval, it may be desirable to also change the warmlayer
        threshold value.

    """
    warmlayer_threshold_OOI = 21600.0  # equatorial sunrise: 6AM

    # in the original code, no checks were included to trap out the kind
    # of situation in which there are data for a given day from local
    # times 0500-1200 immediately followed by data for the following day
    # from 1300-1800.
    #
    # finding the start of each day when the timestamps have units of seconds
    # since 01-jan-1900 is straightforward:
    newday = np.diff(np.floor(localdate/86400.0)) > 0
    # prepend a True value to get the index count correct and to start data at a new day.
    newday = np.hstack((True, newday))
    # find the indices of the start of newdays
    idx_nd = np.nonzero(newday)[0]
    # append a bracketing end index for the for loop to follow
    idx_nd = np.hstack((idx_nd, newday.size))

    # the warmlayer routine is to be run only on days which start earlier than threshold.
    time_of_day = np.mod(localdate, 86400.0)
    earlier_than_threshold = time_of_day <= warmlayer_threshold_OOI

    # initialize warmmask to all False (do not process any data records with warmlayer)
    warmmask = np.zeros(time_of_day.size, dtype=bool)
    # every warmmask value will be overwritten by the loop.
    #..  each iteration processes one day's worth of data records whose indices are
    #..     [idx_nd[ii]:idx_nd[ii+1]] (python does not use last index in its range).
    #..  the values for that day's records are determined by whether the first local
    #..     time for that day is earlier than the threshold (T) or not (F).
    for ii in range(idx_nd.size-1):
        warmmask[idx_nd[ii]:idx_nd[ii+1]] = earlier_than_threshold[idx_nd[ii]]

    idx_warm = np.nonzero(warmmask)[0]
    nanmask = ~warmmask

    return idx_warm, newday, nanmask
