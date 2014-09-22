#!/usr/bin/env python
"""
@package ion_functions.data.met_functions
@file ion_functions/data/met_functions.py
@author Russell Desiderio, Chris Wingard, Stuart Pearce
@brief Module containing functions for the met family of instruments
"""

import numpy as np
import numexpr as ne
from pygsw import vectors as gsw

from ion_functions.data.generic_functions import magnetic_declination, magnetic_correction


# Set global switches used in METBK bulk flux calculations
# These should be set to 1, always!

JCOOLFL = 1      # 1=do coolskin calc
JWARMFL = 1      # 1=do warmlayer calc
#JWAVEFL         # only the windspeed parametrization of the charnok
                 # variable is coded; this switch is not used.

"""
    LISTING OF SUBROUTINES BY ORDER IN THIS MODULE
        Grouped by sections; alphabetical within each section.
        Data product names, formal and meta, are named as "met_prdname".

#...................................................................................
    Functions to compute the L1 BULKMET (METBK) data products:
    these do not require the 'warmlayer/coolskin' iteration algorithm:
        BARPRES
        WINDAVG-VLE  (met_windavg_mag_corr_east)
        WINDAVG-VLN  (met_windavg_mag_corr_north)
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
    These products are calculated at the native temporal resolution of the
    instrument suite (roughly each minute).
#...................................................................................
#...................................................................................
    Functions to compute the (simpler) L2 METBK data products that do not require
    the 'warmlayer/coolskin' iteration algorithm:
        NETSIRR (this may operationally be an L1 product)
        RAINRTE
        SALSURF
        SPECHUM
    These products are calculated at the native temporal resolution of the
    instrument suite (roughly each minute).
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
    Data conditioning and averaging routines
        condition_data
        make_hourly_data
#...................................................................................
#...................................................................................
    seasurface_temp_correct  (wrapper; calls warmlayer and coare35vn)
#...................................................................................
#...................................................................................
    warmlayer ('warmlayer' code)
#...................................................................................
#...................................................................................
    coare35vn ('coolskin' main routine)
        charnock_wind
        coolskin_parameters
        effective_relwind
        obukhov_for_init
        obukhov_length_scale
        roughness_lengths
        roughness_lengths_for_init
        scaling_parameters
#...................................................................................

"""

# METBK SENSOR HEIGHTS:

"""
    Note that the sensor heights may depend on the type of mooring:
"""
#     zwindsp = height of the wind measurement [m]
#     ztmpair = height of air temperature measurement [m]
#     zhumair = height of air humidity measurement [m]
#     ztmpwat = depth of bulk sea surface water measurements [m]

#     zvelptm = depth of surface current measurement [m]:
#         this parameter is specified as metadata in the DPS;
#         however, it is not used in the code.

#     zinvpbl = planetary boundary layer/inversion height: this is
#               set to a default value of 600m in the code, as is
#               used in the DPS.

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

        OOI Level 1 Barometric Pressure core data product, which is calculated
        by scaling the measured barometric pressure from mbar to Pascals.

    Implemented by:

        2014-06-25: Christopher Wingard. Initial code.

    Usage:

        Pa = met_barpres(mbar)

            where

        Pa = Barometric pressure (BARPRES_L1) [Pa]
        mbar = Barometric pressure (BARPRES_L0) [mbar]

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

    These products are calculated at the native temporal resolution of the
    instrument suite (roughly each minute).
#...................................................................................
#...................................................................................

"""


def met_current_direction(vle_water, vln_water):
    """
    Description:

        Calculates the direction of the surface current using the eastward and northward
        velocity components from the VELPT mounted on the surface buoy.

    Implemented by:

        2014-08-27: Russell Desiderio. Initial Code

    Usage:

        current_dir = met_current_direction(vle_water, vln_water)

            where

        current_dir = direction of the surface current (CURRENT_DIR) [0 360) degrees
        vle_water = eastward surface current (VELPTMN-VLE_L1) [m/s]
        vln_water = northward surface current (VELPTMN-VLN_L1) [m/s]

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
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


def met_current_speed(vle_water, vln_water):
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

    Usage:

        current_spd = met_current_speed(vle_water, vln_water)

            where

        current_spd = magnitude (speed) of the surface current (CURRENT_SPD) [m/s]
        vle_water = eastward surface current (VELPTMN-VLE_L1) [m/s]
        vln_water = northward surface current (VELPTMN-VLN_L1) [m/s]

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
    current_spd = ne.evaluate("sqrt(vle_water**2 + vln_water**2)")
    return current_spd


def met_relwind_direction(vle_wind, vln_wind, vle_water, vln_water):
    """
    Description:

        Calculates RELWIND_DIR-AUX, the direction of the vector difference of wind velocity
        (from METBK measurements) and surface current (from VELPT measurements).

        It is anticipated that the wind measurements will be roughly each minute and that the
        current measurements will be broadcast to that resolution.

    Implemented by:

        2014-08-26: Russell Desiderio. Initial Code.

    Usage:

        u_dir = met_relwind_direction(vle_wind, vln_wind, vle_water, vln_water)

            where

        u_dir = direction of relative wind (RELWIND_DIR-AUX) [0 360) degrees
        vle_wind  = eastward wind speed (WINDAVG-VLE_L1) [m/s]
        vln_wind  = northward wind speed (WINDAVG-VLN_L1) [m/s]
        vle_water = eastward surface current (VELPTMN-VLE_L1) [m/s]
        vln_water = northward surface current (VELPTMN-VLN_L1) [m/s]

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
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


def met_relwind_speed(vle_wind, vln_wind, vle_water, vln_water):
    """
    Description:

        Calculates RELWIND_SPD-AUX, the magnitude of the vector difference of surface current
        (from VELPT measurements) from wind velocity (from METBK measurements). This is the
        fundamental windspeed variable used in the METBK toga-coare algorithms.

        It is anticipated that the wind measurements will be roughly each minute and that the
        current measurements will be broadcast to that resolution.

    Implemented by:

        2014-08-26: Russell Desiderio. Initial Code.

    Usage:

        u_rel = met_relwind_speed(vle_wind, vln_wind, vle_water, vln_water)

            where

        u_rel = magnitude of windspeed relative to the ocean (RELWIND_SPD-AUX) [m/s]
        vle_wind  = eastward wind speed (WINDAVG-VLE_L1) [m/s]
        vln_wind  = northward wind speed (WINDAVG-VLN_L1) [m/s]
        vle_water = eastward surface current (VELPTMN-VLE_L1) [m/s]
        vln_water = northward surface current (VELPTMN-VLN_L1) [m/s]

    References:

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)
    """
    u_rel = np.sqrt((vle_water - vle_wind)**2 + (vln_water - vln_wind)**2)
    return u_rel


"""
#...................................................................................
#...................................................................................
    Functions to compute the (simpler) L2 METBK data products that do not require
    the 'warmlayer/coolskin' iteration algorithm:

        NETSIRR (this may operationally be an L1 product)
        RAINRTE
        SALSURF
        SPECHUM

    These products are calculated at the native temporal resolution of the
    instrument suite (roughly each minute).
#...................................................................................
#...................................................................................

"""


def met_netsirr(shortwave_down):
    """
    Description:

        Calculates NETSIRR_L2, the OOI core data product net shortwave radiation
        in the downward direction, for the METBK instrument. This data product may
        have been misclassified (it looks like L1).

    Implemented by:

        2014-08-27: Russell Desiderio. Initial code.

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

        Calculates the sonic buoyancy flux (proposed) data product BUOYFLS_L2. The FDCHP
        instrument provides the somewhat analogous FLUXHOT_L2 data product, and may
        provide the corresponding sonic buoyancy flux product derived from its
        covariance data FLUSHOT_L2.

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

    (usr, tsr, qsr, _, _, _, _, _, _, _, _, _, _, _) = seasurface_temp_correct(*args)

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

        Calculates the buoyancy flux (proposed) data product BUOYFLX_L2. The FDCHP
        instrument provides the corresponding FLUXHOT_L2 data product.

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

    (usr, tsr, qsr, _, _, _, _, _, _, _, _, _, _, _) = seasurface_temp_correct(*args)

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

    (usr, _, qsr, _, _, _, _, _, _, _, _, _, _, _) = seasurface_temp_correct(*args)

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

    (usr, tsr, qsr, _, dter, dqer, _, _, _, _, _, _, _, dsea) = seasurface_temp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (rain_rate, _, _, _, tC_sea, _, _, tC_air, _, relhum, _, pr_air, Rshort_down,
        Rlong_down, _, _, _, _) = args

    cpa = 1004.67  # specific heat capacity of (dry) air [J/kg/K]
    rhoa = air_density(tC_air, pr_air, relhum)
    Le = latent_heat_vaporization_pure_water(tC_sea + dsea)

    hlb = -rhoa * Le * usr * qsr                                         # positive up
    hsb = -rhoa * cpa * usr * tsr                                        # positive up
    Rns_down = met_netsirr(Rshort_down)                                  # positive down
    Rnl_up = net_longwave_up(tC_sea + dsea - dter, Rlong_down)           # positive up
    rainflx = rain_heat_flux(rain_rate, tC_sea, tC_air, relhum, pr_air,  # positive up
                             dter, dqer, dsea)

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
    (usr, _, qsr, _, _, _, _, _, _, _, _, _, _, dsea) = seasurface_temp_correct(*args)

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

    (usr, _, _, ut, _, _, _, _, _, _, _, _, _, _) = seasurface_temp_correct(*args)

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
    (_, _, _, _, dter, _, _, _, _, _, _, _, _, dsea) = seasurface_temp_correct(*args)

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

        Note: As of 22-Sep-2014, there is an unresolved issue with respect to
        how this product is correctly calculated (see code in the subroutine
        rain_heat_flux which is called by met_rainflx). This will affect all
        calculations of data products in this section, because rain heat flux
        is used in the warmlayer calculation.

    Implemented by:

        2014-09-01: Russell Desiderio. Initial Code
        2014-09-19: Russell Desiderio. Added front end to convert eachminute data to hourly.

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

    # dter is the coolskin temperature depression [degC]
    # dqer is the coolskin humidity depression
    # dsea is the warmlayer correction to the sea surface temperature [degC]
    (_, _, _, _, dter, dqer, _, _, _, _, _, _, _, dsea) = seasurface_temp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (rain_rate, _, _, _, tC_sea, _, _, tC_air, _, relhum, _, pr_air,
        _, _, _, _, _, _) = args

    rainflx = rain_heat_flux(rain_rate, tC_sea, tC_air, relhum, pr_air,
                             dter, dqer, dsea)

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

    (usr, tsr, _, _, _, _, _, _, _, _, _, _, _, _) = seasurface_temp_correct(*args)

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
    (_, _, qsr, _, _, _, _, L, _, _, _, _, _, _) = seasurface_temp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, _, _, _, tC_air, _, relhum, _, pr_air, _, _, _, _, _, _) = args

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
    (_, _, _, _, _, _, _, L, _, _, _, _, _, _) = seasurface_temp_correct(*args)

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
    (_, tsr, _, _, _, _, _, L, _, _, _, _, _, _) = seasurface_temp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, _, _, _, tC_air, _, _, _, _, _, _, lat, _, _, _) = args

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
    (_, _, _, _, dter, _, _, _, _, _, _, _, _, dsea) = seasurface_temp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, tC_sea, _, _, _, _, _, _, _, _, _, _, _, _, _) = args

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
    (usr, _, _, ut, _, _, _, L, _, _, _, _, _, _) = seasurface_temp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, _, wnd, _, _, _, _, _, _, _, _, _, _, _, _) = args

    wind10m = windspeed_at_refheight(wnd, usr, ut, zrefht, zwindsp, L)

    return wind10m


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
        air_pressure_mbar = barometric pressure (BARPRES_L0) [mbar]
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

        blackbody energy radiated up from sea surface = eps * sigma * Tsea^4
        reflected longwave radiation up               = (1-eps) * IR
        total longwave radiation down                 = IR

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
    # to a fractional power.
    dzet = np.minimum(50.0, 0.35 * zet)
    psit = -((1.0 + 0.6667 * np.abs(zet))**1.5 +
             0.6667 * (zet - 14.28) * np.exp(-dzet) + 8.525)

    # overwrite psit for zet < 0 values (unstable case).
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


def rain_heat_flux(rain_rate, tC_sea, tC_air, relhum, pr_air, dter, dqer, dsea):
    """

        Calculates net upward rain heat flux.

        Note: As of 22-Sep-2014, there is an unresolved issue with respect to
        how this product is correctly calculated (see code below). This will
        affect all calculations of data products in this section, because
        rain heat flux is used in the warmlayer calculation.

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

    # in Clausius-Clayperon eqn, whoi (DPS) uses tC_sea (+ dsea);
    # jim edson, archived whoi rain_flux.m, and fortran code all use tC_air;
    # neither include dter in this expression

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
    Data conditioning and averaging routines

        condition_data
        make_hourly_data
#...................................................................................
#...................................................................................
"""


def condition_data(*args):
    """
    Description:

        (1) Makes sure that all relevant variables are at least 1D np.arrays.

        (2) For missing input variables with default scalar values set in the
            input argument list, this routine expands the size of those variable
            arrays from 1 element to the size of that of the other variables.

    Implemented by:

        2014-09-19: Russell Desiderio. Initial code.

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
    idx_of_args_to_expand = [0, 11, 12, 13, 14, 15]

    nargs = len(args)
    for ii in range(nargs):
        args[ii] = np.atleast_1d(args[ii])

    # for rainrte and testing
    if nargs < number_of_bulk_vars:
        return args

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
        sporadicly spaced data and for data with time gaps; in the latter
        case, no records are produced for the missing time bins.

    Implemented by:

        2014-09-18: Russell Desiderio. Initial code.

    Usage

        args_out = make_hourly_data(*args_in)

            where

        args_out = argument list of hourly data
        args_in = argument list of each minute data


        If the number of arguments in is 18, then the input are data to be
        prepared for running in the warmlayer and coolskin routines, and so
        a hard-wired list argument index number idx of the variables to be
        processed is executed. For any other args_in sized input, all inputs
        are processed.

        Timestamps in units of seconds must be the second element in args_in;
        this was dictated by convenience for the ordering of the input arguments
        into the seasurface_temp_correct routine.

        All arguments must be 1D arrays of the same length.

    Notes:

        The np.bincount routine is used in much the same way accumarray in matlab
        would be used to construct the hourly data. np.histogram could also be
        used.

        The key to the routine is to convert the timestamps into elapsed hours
        [0.0, 0.1, 0.7, 1.3, 1.7, 3.1, 3.3] so that when these values are floored
        [ 0 ,  0 ,  0 ,  1 ,  1 ,  3 ,  3 ] the entry at an index represents the
        bin number into which the various variables with that index will belong.

        The summing is carried out by using the weighting feature of the np.bincount
        function, as described in the example in the numpy.bincount documentation at:
        http://docs.scipy.org/doc/numpy-1.8.1/reference/generated/numpy.bincount.html.
    """
    args = list(args)

    number_of_bulk_vars = 18

    # timestamps must be the 2nd variable (counting from 0)
    index_timedata = 1

    # prep all variables for rainrte (nargs<18).
    # for nargs >= 18, skip arguments that are constants (except for zinvpbl)
    # and switches.
    idx_to_skip = [17, 16, 10, 8, 6, 3]

    nargs = len(args)
    idx = range(nargs)
    if nargs >= number_of_bulk_vars:
        for ii in idx_to_skip:
            del idx[ii]

    # timestamps are in this element
    time_sec = args[index_timedata]
    time_elapsed_hr = (time_sec - time_sec[0])/3600.0

    # assign each timestamp a bin number based on its elapsed time in hrs.
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

    # a choice:
    # use the midpoint of the bins as the timestamp, instead of the average of
    # the timestamps within the bin; this matters only if there are missing data.
    bin_time_sec = time_sec[0] + 1800.0 + 3600.0 * np.array(range(len(mask)))
    # delete bins with no entries as before
    bin_time_sec = bin_time_sec[mask]

    args[index_timedata] = bin_time_sec

    #for ii in idx:
    #    print args[ii]

    return args


"""
#...................................................................................
#...................................................................................

    Wrapper function which calls the warmlayer and coolskin (coare35vn) routines:

        seasurface_temp_correct

#...................................................................................
#...................................................................................
"""


def seasurface_temp_correct(*args):
    """
    Description:

        Wrapper function which by OOI default applies both of the METBK seasurface
        temperature correction algorithms (warmlayer, coolskin in coare35vn). This
        behavior is set by the global switches JWARMFL=1 and JCOOLFL=1. The switch
        construction is retained for generality.

        Most of the METBK L2 data products and 2 of the metadata products require
        the skin corrections to be applied before their values can be calculated.

    Implemented by:

        2014-09-01: Russell Desiderio. Initial code.

    Usage (command line spaced out for clarity):

        (usr, tsr, qsr, ut, dter, dqer, tkt, L, zou, zot, zoq,     # coare35vn output
        dt_wrm, tk_pwp, dsea) =                                    # warmlayer output

        seasurface_temp_correct

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
            dqer = coolskin humidity depression [degC, but this must be wrong]
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
        # the warmlayer code requires at least 2 time records to run.
        # this is checked in the warmlayer code itself.
        (dt_wrm, tk_pwp, dsea) = warmlayer(*args[0:-1])
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
    coolskin_args = (args[4]+dsea,) + args[5:-1]
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

        Warmlayer code refactored from coare35vnWarm.m (see DPS artifacts in References).

    Implemented by:

        2014-09-01: Russell Desiderio. Initial code.

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

    References:

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
    # this routine requires at least 2 data records:
    if wnd.size < 2:
        raise ValueError('There must be at least 2 data records to run warmlayer code.')

    # set constants
    c2k = 273.15        # Converts degC to Kelvin
    cpw = 4000.0        # Specific heat of water
    rhow = 1022.0       # Density of seawater
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

    # in the original code, no checks were included to trap out the kind
    # of situation in which there are data for a given day from local
    # times 0500-1200 immediately followed by data for the following day
    # from 1300-1800.

    # local solar time adjustment is a function of longitude:
    #.. 360 degrees = 24 *3600 seconds,
    #.. so each degree is worth 240 seconds of time.

    # the OOI timestamp is seconds since midnight 01-jan-1900; therefore
    # local time will still be positive for the case of lon = -180deg.
    localtime_sec = np.mod(timestamp + lon * 240.0, 86400.0)

    #.. and calculate all the delta times;
    #.. prepend a zero to line up delta_time with iteration number.
    delta_time = np.hstack((0.0, np.diff(localtime_sec)))

    #.. the original code used a 'jump' variable, the action of which
    #.. was to do the first loop calculation at the first encounter of
    #.. localtime <= 6:00 AM, then continue processing all records.
    #.. replace this with directly finding first local time <= 6:00 AM:
    ibg_begin = np.where(localtime_sec <= 21600.0)[0]

    # check to see if ibg_begin is empty:
    if ibg_begin.size == 0:
        raise ValueError('There must be at least one local time earlier than '
                         '6:00 AM to run warmlayer code.')

    #.. original code initializes with the first data record, then
    #.. starts checking whether to start processing with the 2nd.
    #.. refactored matlab code initializes the loop variables so that the
    #.. 1st iteration processes the 2nd data record (with loop index=2).
    #
    #.. python starts index numbering at 0, so:
    ibg_begin = np.maximum(ibg_begin[0], 1)

    #.. the original code has been changed to show the explicit dependence
    #.. of the variables upon iteration count (data record number).
    for ibg in range(ibg_begin, nx):   # step through each timepoint

        #.. warmlayer values for the following case are just the initialized
        #.. values. so, instead of using if-then-else, simplify indentation
        #.. by using 'if' only, reset variables, and jump to next iteration.
        if delta_time[ibg] < 0.0:  # re-zero at midnight
            # dt_wrm[ibg] = 0.0;
            # tk_pwp[ibg] = max_pwp;
            # dsea[ibg]   = 0.0;
            jamset = 0
            fxp = 0.5
            tau_ac = 0.0
            qcol_ac = 0.0
            continue  # go to next time (data) record
            # end midnight reset

        #*****  dependent variables for the [ibg]th warm layer calculation
        #*****  of dsea are fluxes, coolskin correction dter, and dsea itself,
        #*****  which are derived from the previous ([ibg-1]th) data record.
        #
        # because of the dependence on the previous value of dsea, this calculation
        # cannot be vectorized.
        tsea_corr = tC_sea[ibg-1] + dsea[ibg-1]

        # slicing 1D arrays with [ibg-1:ibg] returns a 1-element nd.array variable which
        # can be indexed, whereas slicing with [ibg-1] returns a variable which cannot be
        # indexed. [ibg-1:ibg] slicing is used so that coare35vn can be run with both
        # 'scalar' and 'vector' input.
        args = (tsea_corr, wnd[ibg-1:ibg], zwindsp, tC_air[ibg-1:ibg], ztmpair, relhum[ibg-1:ibg], zhumair,
                pr_air[ibg-1:ibg], Rshort_down[ibg-1:ibg], Rlong_down[ibg-1:ibg], lat[ibg-1:ibg], zinvpbl[ibg-1:ibg],
                jcool)

        (usr, tsr, qsr, ut, dter, dqer, _, _, _, _, _) = coare35vn(*args)

        # in the original matlab code, Le was calculated inside of the coare35vn
        # subroutine, which was called using tC_sea+dsea for seawater temperature:
        Le = latent_heat_vaporization_pure_water(tsea_corr)
        tau_old = rhoa[ibg-1] * usr * usr * wnd[ibg-1] / ut  # stress
        hs_old = -rhoa[ibg-1] * cpa * usr * tsr              # sensible heat flux
        hl_old = -rhoa[ibg-1] * Le * usr * qsr                      # latent heat flux

        # note:
        #     the original matlab code is followed here: it does not use dsea
        #     in the Rnl expression used in the warmlayer calculation, although
        #     dsea is used in the expression for RF_old.
        Rnl = net_longwave_up(tC_sea[ibg]-dter, Rlong_down[ibg])
        RF_old = rain_heat_flux(rain_rate[ibg-1], tC_sea[ibg-1], tC_air[ibg-1], relhum[ibg-1],
                                pr_air[ibg-1], dter, dqer, dsea[ibg-1])

        #********************************************************
        #****  Compute warm layer correction *******************
        #********************************************************
        qr_out = Rnl + hs_old + hl_old + RF_old  # total cooling at surface
        q_pwp = fxp * Rns[ibg] - qr_out          # tot heat abs in warm layer

        # calculate dt_wrm and tk_pwp for this iteration.
        if q_pwp >= 50.0 or jamset == 1:         # Check for threshold
            jamset = 1			         # indicates threshold crossed
            tau_ac = tau_ac + np.maximum(.002, tau_old) * delta_time[ibg]  # momentum integral

            # check threshold for warm layer existence
            if qcol_ac + q_pwp * delta_time[ibg] > 0.0:
                #******************************************
                # Compute the absorption profile
                #******************************************
                #.. tk_pwp can iteratively change value in the following loop,
                #.. requiring the creation of the variable tkpwp.
                tkpwp = tk_pwp[ibg-1]
                for i in range(5):               # loop 5 times for fxp
                    fxp = 1.0 - (0.28 * 0.014 * (1.0 - np.exp(-tkpwp / 0.014)) +
                                 0.27 * 0.357 * (1.0 - np.exp(-tkpwp / 0.357)) +
                                 0.45 * 12.82 * (1.0 - np.exp(-tkpwp / 12.82))) / tkpwp
                    qjoule = (fxp * Rns[ibg] - qr_out) * delta_time[ibg]
                    if qcol_ac + qjoule > 0.0:   # Compute warm-layer depth
                        tkpwp = np.minimum(max_pwp,
                                           ctd1[ibg] * tau_ac / np.sqrt(qcol_ac + qjoule))
                tk_pwp[ibg] = tkpwp
            else:                                # warm layer wiped out
                fxp = 0.75
                tk_pwp[ibg] = max_pwp
                qjoule = (fxp * Rns[ibg] - qr_out) * delta_time[ibg]
            qcol_ac = qcol_ac + qjoule           # heat integral

            #*******  compute dt_warm  ******
            if qcol_ac > 0.0:
                dt_wrm[ibg] = ctd2[ibg] * (qcol_ac)**1.5 / tau_ac
            else:
                dt_wrm[ibg] = 0.0

        else:   # propagate dt_wrm and tk_pwp values
            dt_wrm[ibg] = dt_wrm[ibg-1]
            tk_pwp[ibg] = tk_pwp[ibg-1]

        # Compute warm layer correction dsea
        if tk_pwp[ibg] < ztmpwat:
            dsea[ibg] = dt_wrm[ibg]
        else:
            dsea[ibg] = dt_wrm[ibg] * ztmpwat / tk_pwp[ibg]

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
        of coare35vn.m (see DPS artifacts in References).

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
            dqer = coolskin humidity depression [degC, but this must be wrong]
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

        OOI (2014). Data Product Specification for L2 BULKFLX Data Products.
            Document Control Number 1341-00370.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00370_Data_Product_Spec_BULKFLX_OOI.pdf)

        OOI (2014). 1341-00370_BULKFLX Artifacts. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> REFERENCE >> Data Product Specification Artifacts
            >> 1341-00370_BULKFLX  (Original matlab code).

        Fairall, C.W., E.F. Bradley, J.E. Hare, A.A. Grachev, and J.B. Edson (2003),
        Bulk parameterization of air sea fluxes: updates and verification for the
        COARE algorithm, J. Climate, 16, 571-590.

    Notes:

        The code within the iteration loop was segregated into subroutines to clarify which
        of the variables were key in that they were recalculated during each iteration, and
        afterwards used to calculate data products.

        These key variables are the output of this stripped down version of the original code.
        This code is not meant to calculate data products; rather, it is to provide the key
        variables from which all bulk parameters can be calculated.
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
    #              qsr must have units of g/kg
    # qsr = qsr * 1000    # changes units from kg/kg to g/kg

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
    ug[Bf > 0] = np.maximum(0.2,
                            Beta * (Bf[Bf > 0] * zinvpbl[Bf > 0]) ** 0.333)
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
    k50 = np.where(zetu > 50)  # stable with very thin M-O length relative to zwindsp
    Ribcu = -zwindsp / zinvpbl / 0.004 / Beta**3
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
    mask = np.abs(tvsr) < 1.e-3
    tvsr[mask] = np.abs(tvsr[mask])
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
