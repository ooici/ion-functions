#!/usr/bin/env python
"""
@package ion_functions.data.fdc_functions
@file ion_functions/data/fdc_functions.py
@author Russell Desiderio
@brief Module containing FDC related data-calculations.
"""

import numpy as np
import scipy as sp
from scipy import integrate
from scipy import interpolate
from scipy import signal


"""
#...................................................................................
#...................................................................................
    The FDCHP instrument outputs a dataset of 10 Hz data for approximately 20 minutes
        every hour, resulting in around 12000 data records for each dataset.
    The FDCHP data product algorithms parse incoming 1D array data into separate
        datasets and operate on these datasets individually.
    Each of the four L1 products (wind and temperature) consists of 11400 values for
        each 20 minute dataset.
    Each of the three L2 flux products consists of one value per 20 minute dataset.

    Two additional auxiliary data products, not specified in the DPS, have been coded
        to provide time bases for the L1 and L2 products.
#...................................................................................
#...................................................................................

    LISTING OF SUBROUTINES BY ORDER IN THIS MODULE
        Grouped by sections; alphabetical within each section.

#...................................................................................
#...................................................................................
    Functions to compute the L1 FDCHP data products:

        fdc_tmpatur:        TMPATUR
        fdc_windtur_north:  WINDTUR-VLN
        fdc_windtur_up:     WINDTUR-VLU
        fdc_windtur_west:   WINDTUR-VLW
#...................................................................................
#...................................................................................
    Functions to compute the L2 FDCHP data products:

        fdc_fluxhot:            FLUXHOT
        fdc_fluxmom_alongwind:  FLUXMOM-U
        fdc_fluxmom_crosswind:  FLUXMOM-V
#...................................................................................
#...................................................................................
    Functions to compute the auxiliary time base data products:

        fdc_time_L1:  TIME_L1-AUX
        fdc_time_L2:  TIME_L2-AUX
#...................................................................................
#...................................................................................
    Primary routine to directly compute L1 wind products and L2 flux products:

        fdc_flux_and_wind
#...................................................................................
#...................................................................................
    Subroutines called by the primary routine fdc_flux_and_wind and its subroutines:

        fdc_accelsclimode
        fdc_alignwind
        fdc_anglesclimodeyaw
        fdc_despikesimple
        fdc_detrend
        fdc_filtcoef
        fdc_grv
        fdc_process_compass_data
        fdc_quantize_data
        fdc_sonic
        fdc_trans
        fdc_update
#...................................................................................
#...................................................................................

"""
####################################################################################
####################################################################################
####################################################################################
"""
#...................................................................................
#...................................................................................
    Functions to compute the L1 FDCHP data products:

        fdc_tmpatur:        TMPATUR
        fdc_windtur_north:  WINDTUR-VLN
        fdc_windtur_up:     WINDTUR-VLU
        fdc_windtur_west:   WINDTUR-VLW
#...................................................................................
#...................................................................................
"""


def fdc_tmpatur(timestamp, sonicT):
    """
    Description:

        Calculates the L1 temperature data product TMPATUR_L1 from the FDCHP
        instrument, which collects 20 minutes of data every hour. The L1 data
        consists of these values less 30 seconds from both the beginning and
        end of each 12000 point dataset.

        The L1 temperature data product is also separately calculated in the
        routine fdc_flux_and_wind which directly calculates the fluxhot product.

    Implemented by:

        2014-11-17: Russell Desiderio. Initial Code

    Usage:

        Ts = fdc_tmpatur(timestamp, sonicT)

            where

        Ts = sonic temperature [K]
        timestamp = data date and time values [seconds since 1900-01-01]
        sonicT = TMPATUR_L0 [counts]; speed of sound measured by the sonic anemometer

    References:

        OOI (2014). Data Product Specification for FDCHP Data Products. Document
            Control Number 1341-00280. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00280_Data_Product_Spec_FDCHP_OOI.pdf)
    """
    # condition data and parse it into discrete datasets.
    data = fdc_quantize_data(timestamp, sonicT)
    # the shape of the data array is [n_var, n_packets, pts per packet].

    # there is no despiking or filtering to be done on these data, so truncate now.
    # number of seconds of data to remove from the beginning and end of the processed
    # data before calculating the mean of the products of the elements of two vectors.
    edge_sec = 30
    # sampling frequency
    fs = 10
    # number of edge data values to remove, based on sampling frequency fs in Hz
    edge = fs * edge_sec
    data = data[:, :, edge:-edge]

    # for clarity in following the DPS code, unpack the data into its constituent
    # variables as 2D arrays so that the index of the lead dimension of each
    # variable array indicates the dataset number; sonicT[0, :] will be a vector
    # containing the L0 temperature data for the first dataset packet.
    sonicT = data[1, :, :]

    # process L0 temperature data
    Ts = 0.01 * sonicT
    Ts = Ts * Ts / 403.0
    Ts = Ts.flatten()

    return Ts


def fdc_windtur_north(timestamp, sonicU, sonicV, sonicW, sonicT, heading,
                      rateX, rateY, rateZ, accX, accY, accZ, lat):
    """
    Description:

        Calculates the L1 windspeed data product WINDTUR-VLN_L1 from the FDCHP
        instrument, which collects 20 minutes of data every hour. The L1 data
        consists of these values less 30 seconds from both the beginning and
        end of each 12000 point dataset.

    Implemented by:

        2014-11-17: Russell Desiderio. Initial Code

    Usage:

        wind_north = fdc_windtur_north(timestamp, sonicU, sonicV, sonicW, sonicT, heading,
                                       rateX, rateY, rateZ, accX, accY, accZ, lat)

            where

        wind_north = windspeed North WINDTUR-VLN_L1 [m/s], UNcorrected for magnetic variation
        timestamp = data date and time values [seconds since 1900-01-01]
        sonicU = WINDTUR-U_L0 [cm/s]; u-component of windspeed measured in the buoy
                 frame of reference
        sonicV = WINDTUR-V_L0 [cm/s]; v-component of windspeed measured in the buoy
                 frame of reference
        sonicW = WINDTUR-W_L0 [cm/s]; w-component of windspeed measured in the buoy
                 frame of reference
        sonicT = TMPATUR_L0 [counts]; speed of sound measured by the sonic anemometer
        heading = MOTFLUX-YAW_L0 [radians] measured by the magnetometer (NOT msrd by the gyro).
        ***NOT USED*** roll: MOTFLUX-ROLL_L0 [radians] ***NOT USED***
        ***NOT USED*** pitch: MOTFLUX-PITCH_L0 [radians] ***NOT USED***
        rateX = MOTFLUX-ROLL_RATE_L0 [radians/s] measured by the gyro
        rateY = MOTFLUX-PITCH_RATE_L0 [radians/s] measured by the gyro
        rateZ = MOTFLUX-YAW_RATE_L0 [radians/s] measured by the gyro
        accX = MOTFLUX-ACX_L0 [9.80665 m^2/s^2] x-component of platform linear acceleration
        accY = MOTFLUX-ACY_L0 [9.80665 m^2/s^2] y-component of platform linear acceleration
        accZ = MOTFLUX-ACZ_L0 [9.80665 m^2/s^2] z-component of platform linear acceleration
        lat = latitude of instrument in decimal degrees

    References:

        OOI (2014). Data Product Specification for FDCHP Data Products. Document
            Control Number 1341-00280. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00280_Data_Product_Spec_FDCHP_OOI.pdf)
    """
    _, windspeeds = fdc_flux_and_wind(timestamp, sonicU, sonicV, sonicW, sonicT,
                                      heading, rateX, rateY, rateZ, accX, accY,
                                      accZ, lat)

    wind_north = np.asarray(windspeeds[0]).flatten()

    return wind_north


def fdc_windtur_up(timestamp, sonicU, sonicV, sonicW, sonicT, heading,
                   rateX, rateY, rateZ, accX, accY, accZ, lat):
    """
    Description:

        Calculates the L1 windspeed data product WINDTUR-VLU_L1 from the FDCHP
        instrument, which collects 20 minutes of data every hour. The L1 data
        consists of these values less 30 seconds from both the beginning and
        end of each 12000 point dataset.

    Implemented by:

        2014-11-17: Russell Desiderio. Initial Code

    Usage:

        wind_up = fdc_windtur_up(timestamp, sonicU, sonicV, sonicW, sonicT, heading,
                                 rateX, rateY, rateZ, accX, accY, accZ, lat)

            where

        wind_up = windspeed UP WINDTUR-VLU_L1 [m/s]
        timestamp = data date and time values [seconds since 1900-01-01]
        sonicU = WINDTUR-U_L0 [cm/s]; u-component of windspeed measured in the buoy
                 frame of reference
        sonicV = WINDTUR-V_L0 [cm/s]; v-component of windspeed measured in the buoy
                 frame of reference
        sonicW = WINDTUR-W_L0 [cm/s]; w-component of windspeed measured in the buoy
                 frame of reference
        sonicT = TMPATUR_L0 [counts]; speed of sound measured by the sonic anemometer
        heading = MOTFLUX-YAW_L0 [radians] measured by the magnetometer (NOT msrd by the gyro).
        ***NOT USED*** roll: MOTFLUX-ROLL_L0 [radians] ***NOT USED***
        ***NOT USED*** pitch: MOTFLUX-PITCH_L0 [radians] ***NOT USED***
        rateX = MOTFLUX-ROLL_RATE_L0 [radians/s] measured by the gyro
        rateY = MOTFLUX-PITCH_RATE_L0 [radians/s] measured by the gyro
        rateZ = MOTFLUX-YAW_RATE_L0 [radians/s] measured by the gyro
        accX = MOTFLUX-ACX_L0 [9.80665 m^2/s^2] x-component of platform linear acceleration
        accY = MOTFLUX-ACY_L0 [9.80665 m^2/s^2] y-component of platform linear acceleration
        accZ = MOTFLUX-ACZ_L0 [9.80665 m^2/s^2] z-component of platform linear acceleration
        lat = latitude of instrument in decimal degrees

    References:

        OOI (2014). Data Product Specification for FDCHP Data Products. Document
            Control Number 1341-00280. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00280_Data_Product_Spec_FDCHP_OOI.pdf)
    """
    _, windspeeds = fdc_flux_and_wind(timestamp, sonicU, sonicV, sonicW, sonicT,
                                      heading, rateX, rateY, rateZ, accX, accY,
                                      accZ, lat)

    wind_up = np.asarray(windspeeds[2]).flatten()

    return wind_up


def fdc_windtur_west(timestamp, sonicU, sonicV, sonicW, sonicT, heading,
                     rateX, rateY, rateZ, accX, accY, accZ, lat):
    """
    Description:

        Calculates the L1 windspeed data product WINDTUR-VLW_L1 from the FDCHP
        instrument, which collects 20 minutes of data every hour. The L1 data
        consists of these values less 30 seconds from both the beginning and
        end of each 12000 point dataset.

    Implemented by:

        2014-11-17: Russell Desiderio. Initial Code

    Usage:

        wind_west = fdc_windtur_west(timestamp, sonicU, sonicV, sonicW, sonicT, heading,
                                     rateX, rateY, rateZ, accX, accY, accZ, lat)

            where

        wind_west = windspeed West WINDTUR-VLW_L1 [m/s], UNcorrected for magnetic variation
        timestamp = data date and time values [seconds since 1900-01-01]
        sonicU = WINDTUR-U_L0 [cm/s]; u-component of windspeed measured in the buoy
                 frame of reference
        sonicV = WINDTUR-V_L0 [cm/s]; v-component of windspeed measured in the buoy
                 frame of reference
        sonicW = WINDTUR-W_L0 [cm/s]; w-component of windspeed measured in the buoy
                 frame of reference
        sonicT = TMPATUR_L0 [counts]; speed of sound measured by the sonic anemometer
        heading = MOTFLUX-YAW_L0 [radians] measured by the magnetometer (NOT msrd by the gyro).
        ***NOT USED*** roll: MOTFLUX-ROLL_L0 [radians] ***NOT USED***
        ***NOT USED*** pitch: MOTFLUX-PITCH_L0 [radians] ***NOT USED***
        rateX = MOTFLUX-ROLL_RATE_L0 [radians/s] measured by the gyro
        rateY = MOTFLUX-PITCH_RATE_L0 [radians/s] measured by the gyro
        rateZ = MOTFLUX-YAW_RATE_L0 [radians/s] measured by the gyro
        accX = MOTFLUX-ACX_L0 [9.80665 m^2/s^2] x-component of platform linear acceleration
        accY = MOTFLUX-ACY_L0 [9.80665 m^2/s^2] y-component of platform linear acceleration
        accZ = MOTFLUX-ACZ_L0 [9.80665 m^2/s^2] z-component of platform linear acceleration
        lat = latitude of instrument in decimal degrees

    References:

        OOI (2014). Data Product Specification for FDCHP Data Products. Document
            Control Number 1341-00280. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00280_Data_Product_Spec_FDCHP_OOI.pdf)
    """
    _, windspeeds = fdc_flux_and_wind(timestamp, sonicU, sonicV, sonicW, sonicT,
                                      heading, rateX, rateY, rateZ, accX, accY,
                                      accZ, lat)

    wind_west = np.asarray(windspeeds[1]).flatten()

    return wind_west


"""
#...................................................................................
#...................................................................................
    Functions to compute the L2 FDCHP data products:

        fdc_fluxhot:            FLUXHOT
        fdc_fluxmom_alongwind:  FLUXMOM-U
        fdc_fluxmom_crosswind:  FLUXMOM-V
#...................................................................................
#...................................................................................
"""


def fdc_fluxhot(timestamp, sonicU, sonicV, sonicW, sonicT, heading,
                rateX, rateY, rateZ, accX, accY, accZ, lat):
    """
    Description:

        Calculates FLUXHOT_L2, the sonic buoyancy flux, from the FDCHP instrument, which
        collects 20 minutes of data every hour. There is one FLUXHOT value calculated
        for each 20 minute dataset of 12000 data records.

    Implemented by:

        2014-11-17: Russell Desiderio. Initial Code

    Usage:

        fluxhot = fdc_fluxhot(timestamp, sonicU, sonicV, sonicW, sonicT,
                              heading, rateX, rateY, rateZ, accX, accY,
                              accZ, lat)

            where

        fluxhot = FLUXHOT_L2, the sonic buoyancy flux [m/s * K]
        timestamp = data date and time values [seconds since 1900-01-01]
        sonicU = WINDTUR-U_L0 [cm/s]; u-component of windspeed measured in the buoy
                 frame of reference
        sonicV = WINDTUR-V_L0 [cm/s]; v-component of windspeed measured in the buoy
                 frame of reference
        sonicW = WINDTUR-W_L0 [cm/s]; w-component of windspeed measured in the buoy
                 frame of reference
        sonicT = TMPATUR_L0 [counts]; speed of sound measured by the sonic anemometer
        heading = MOTFLUX-YAW_L0 [radians] measured by the magnetometer (NOT msrd by the gyro).
        ***NOT USED*** roll: MOTFLUX-ROLL_L0 [radians] ***NOT USED***
        ***NOT USED*** pitch: MOTFLUX-PITCH_L0 [radians] ***NOT USED***
        rateX = MOTFLUX-ROLL_RATE_L0 [radians/s] measured by the gyro
        rateY = MOTFLUX-PITCH_RATE_L0 [radians/s] measured by the gyro
        rateZ = MOTFLUX-YAW_RATE_L0 [radians/s] measured by the gyro
        accX = MOTFLUX-ACX_L0 [9.80665 m^2/s^2] x-component of platform linear acceleration
        accY = MOTFLUX-ACY_L0 [9.80665 m^2/s^2] y-component of platform linear acceleration
        accZ = MOTFLUX-ACZ_L0 [9.80665 m^2/s^2] z-component of platform linear acceleration
        lat = latitude of instrument in decimal degrees

    References:

        OOI (2014). Data Product Specification for FDCHP Data Products. Document
            Control Number 1341-00280. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00280_Data_Product_Spec_FDCHP_OOI.pdf)
    """
    fluxes, _ = fdc_flux_and_wind(timestamp, sonicU, sonicV, sonicW, sonicT,
                                  heading, rateX, rateY, rateZ, accX, accY,
                                  accZ, lat)

    fluxhot = fluxes[2]

    return fluxhot


def fdc_fluxmom_alongwind(timestamp, sonicU, sonicV, sonicW, sonicT, heading,
                          rateX, rateY, rateZ, accX, accY, accZ, lat):
    """
    Description:

        Calculates FLUXMOM-U_L2, the along-wind component of the momentum flux, from
        the FDCHP instrument, which collects 20 minutes of data every hour. There is one
        FLUXMOM-U_L2 value calculated for each 20 minute dataset of 12000 data records.

    Implemented by:

        2014-11-17: Russell Desiderio. Initial Code

    Usage:

        fluxmom_along = fdc_fluxmom_alongwind(timestamp, sonicU, sonicV, sonicW, sonicT,
                                              heading, rateX, rateY, rateZ, accX, accY,
                                              accZ, lat)

            where

        fluxmom_along = FLUXMOM-U_L2, the along-wind component of the momentum flux [m^2/s^2]
        timestamp = data date and time values [seconds since 1900-01-01]
        sonicU = WINDTUR-U_L0 [cm/s]; u-component of windspeed measured in the buoy
                 frame of reference
        sonicV = WINDTUR-V_L0 [cm/s]; v-component of windspeed measured in the buoy
                 frame of reference
        sonicW = WINDTUR-W_L0 [cm/s]; w-component of windspeed measured in the buoy
                 frame of reference
        sonicT = TMPATUR_L0 [counts]; speed of sound measured by the sonic anemometer
        heading = MOTFLUX-YAW_L0 [radians] measured by the magnetometer (NOT msrd by the gyro).
        ***NOT USED*** roll: MOTFLUX-ROLL_L0 [radians] ***NOT USED***
        ***NOT USED*** pitch: MOTFLUX-PITCH_L0 [radians] ***NOT USED***
        rateX = MOTFLUX-ROLL_RATE_L0 [radians/s] measured by the gyro
        rateY = MOTFLUX-PITCH_RATE_L0 [radians/s] measured by the gyro
        rateZ = MOTFLUX-YAW_RATE_L0 [radians/s] measured by the gyro
        accX = MOTFLUX-ACX_L0 [9.80665 m^2/s^2] x-component of platform linear acceleration
        accY = MOTFLUX-ACY_L0 [9.80665 m^2/s^2] y-component of platform linear acceleration
        accZ = MOTFLUX-ACZ_L0 [9.80665 m^2/s^2] z-component of platform linear acceleration
        lat = latitude of instrument in decimal degrees

    References:

        OOI (2014). Data Product Specification for FDCHP Data Products. Document
            Control Number 1341-00280. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00280_Data_Product_Spec_FDCHP_OOI.pdf)
    """
    fluxes, _ = fdc_flux_and_wind(timestamp, sonicU, sonicV, sonicW, sonicT,
                                  heading, rateX, rateY, rateZ, accX, accY,
                                  accZ, lat)

    fluxmom_along = fluxes[0]

    return fluxmom_along


def fdc_fluxmom_crosswind(timestamp, sonicU, sonicV, sonicW, sonicT, heading,
                          rateX, rateY, rateZ, accX, accY, accZ, lat):
    """
    Description:

        Calculates FLUXMOM-V_L2, the cross-wind component of the momentum flux, from
        the FDCHP instrument, which collects 20 minutes of data every hour. There is one
        FLUXMOM-V_L2 value calculated for each 20 minute dataset of 12000 data records.

    Implemented by:

        2014-11-17: Russell Desiderio. Initial Code

    Usage:

        fluxmom_cross = fdc_fluxmom_crosswind(timestamp, sonicU, sonicV, sonicW, sonicT,
                                              heading, rateX, rateY, rateZ, accX, accY,
                                              accZ, lat)

            where

        fluxmom_cross = FLUXMOM-V_L2, the cross-wind component of the momentum flux [m^2/s^2]
        timestamp = data date and time values [seconds since 1900-01-01]
        sonicU = WINDTUR-U_L0 [cm/s]; u-component of windspeed measured in the buoy
                 frame of reference
        sonicV = WINDTUR-V_L0 [cm/s]; v-component of windspeed measured in the buoy
                 frame of reference
        sonicW = WINDTUR-W_L0 [cm/s]; w-component of windspeed measured in the buoy
                 frame of reference
        sonicT = TMPATUR_L0 [counts]; speed of sound measured by the sonic anemometer
        heading = MOTFLUX-YAW_L0 [radians] measured by the magnetometer (NOT msrd by the gyro).
        ***NOT USED*** roll: MOTFLUX-ROLL_L0 [radians] ***NOT USED***
        ***NOT USED*** pitch: MOTFLUX-PITCH_L0 [radians] ***NOT USED***
        rateX = MOTFLUX-ROLL_RATE_L0 [radians/s] measured by the gyro
        rateY = MOTFLUX-PITCH_RATE_L0 [radians/s] measured by the gyro
        rateZ = MOTFLUX-YAW_RATE_L0 [radians/s] measured by the gyro
        accX = MOTFLUX-ACX_L0 [9.80665 m^2/s^2] x-component of platform linear acceleration
        accY = MOTFLUX-ACY_L0 [9.80665 m^2/s^2] y-component of platform linear acceleration
        accZ = MOTFLUX-ACZ_L0 [9.80665 m^2/s^2] z-component of platform linear acceleration
        lat = latitude of instrument in decimal degrees

    References:

        OOI (2014). Data Product Specification for FDCHP Data Products. Document
            Control Number 1341-00280. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00280_Data_Product_Spec_FDCHP_OOI.pdf)
    """
    fluxes, _ = fdc_flux_and_wind(timestamp, sonicU, sonicV, sonicW, sonicT,
                                  heading, rateX, rateY, rateZ, accX, accY,
                                  accZ, lat)

    fluxmom_cross = fluxes[1]

    return fluxmom_cross


"""
#...................................................................................
#...................................................................................
    Functions to compute the auxiliary time base data products:

        fdc_time_L1:  TIME_L1-AUX
        fdc_time_L2:  TIME_L2-AUX
#...................................................................................
#...................................................................................
"""


def fdc_time_L1(timestamp):
    """
    Description:

        Calculates the time metadata product TIME_L1-AUX associated with the L1 wind
        and temperature data products from the FDCHP instrument, which collects 20
        minutes of data every hour. For each of the L1 data products, 30 seconds
        of data is stripped out from the beginning and end of each 20 minute, 12000
        record dataset. The TIME_L1-AUX values are the remaining timestamps for each
        dataset.

        There were no time metadata products specified in the DPS; however, it is
        currently not clear to me how in the current CI times would be associated
        with the L1 data products without this metaproduct.

    Implemented by:

        2014-11-17: Russell Desiderio. Initial Code

    Usage:

        time_L1 = fdc_time_L1(timestamp)

            where

        time_L1 = timestamps associated with the L1 data products [seconds since 1900-01-01]
        timestamp = input data date and time values [seconds since 1900-01-01]

    References:

        OOI (2014). Data Product Specification for FDCHP Data Products. Document
            Control Number 1341-00280. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00280_Data_Product_Spec_FDCHP_OOI.pdf)
    """
    # condition data and parse it into discrete datasets.
    data = fdc_quantize_data(timestamp)
    # the shape of the data array is [n_var, n_packets, pts per packet].

    # number of seconds of data to remove from the beginning and end of the processed data
    edge_sec = 30
    # sampling frequency
    fs = 10
    # number of edge data values to remove, based on sampling frequency fs in Hz
    edge = fs * edge_sec
    data = data[:, :, edge:-edge]

    # for clarity in following the DPS code, unpack the data into its constituent
    # variables as 2D arrays so that the index of the lead dimension of each
    # variable array indicates the dataset number.
    tmstmp = data[0, :, :]

    time_L1 = tmstmp.flatten()

    return time_L1


def fdc_time_L2(timestamp):
    """
    Description:

        Calculates the time metadata product TIME_L2-AUX associated with the L2
        flux data products from the FDCHP instrument. FDCHP collects 20 minutes
        of data every hour; for each of the L2 flux data products, one data value
        is calculated for each 20 minute dataset. The TIME_L2-AUX values are the
        median timestamps for each dataset.

        There were no time metadata products specified in the DPS; however, it is
        currently not clear to me how in the current CI times would be associated
        with the L2 data products without this metaproduct.

    Implemented by:

        2014-11-17: Russell Desiderio. Initial Code

    Usage:

        time_L2 = fdc_time_L2(timestamp)

            where

        time_L2 = timestamps associated with the L2 data products [seconds since 1900-01-01]
                  time_L2 is a numpy array, one timestamp for each dataset packet.
        timestamp = data date and time values [seconds since 1900-01-01]

    References:

        OOI (2014). Data Product Specification for FDCHP Data Products. Document
            Control Number 1341-00280. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00280_Data_Product_Spec_FDCHP_OOI.pdf)
    """
    # condition data and parse it into discrete datasets.
    data = fdc_quantize_data(timestamp)
    # the shape of the data array is [n_var, n_packets, pts per packet].

    # number of seconds of data to remove from the beginning and end of the processed
    # data before calculating the mean of the products of the elements of two vectors.
    edge_sec = 30
    # sampling frequency
    fs = 10
    # number of edge data values to remove, based on sampling frequency fs in Hz
    edge = fs * edge_sec
    data = data[:, :, edge:-edge]

    # for clarity in following the DPS code, unpack the data into its constituent
    # variables as 2D arrays so that the index of the lead dimension of each
    # variable array indicates the dataset number.
    tmstmp = data[0, :, :]

    # round the median L1 time values to get the L2 flux timestamps
    time_L2 = np.around(np.median(tmstmp, axis=-1))

    return time_L2


"""
#...................................................................................
#...................................................................................
    Primary routine to directly compute L1 wind products and L2 flux products:

        fdc_flux_and_wind
#...................................................................................
#...................................................................................
"""


def fdc_flux_and_wind(timestamp, sonicU, sonicV, sonicW, sonicT, heading,
                      rateX, rateY, rateZ, accX, accY, accZ, lat):
    """
    Description:

        Calculates the 3 L2 flux data products and the 3 L1 wind direction data products
        from the FDCHP instrument. It is anticipated that wrapper functions will be written
        that will call this routine in order to furnish discrete data products as originally
        (and possibly still) required by OOI CyberInfrastructure. This python code is derived
        from Matlab code updated by the DPS author from the code in the original DPS.

    Implemented by:

        2014-05-20: Russell Desiderio. Initial Code
        2014-11-06: Russell Desiderio. Incorporated fcd_quantize_data routine.

    Usage:

        fluxes, windspeeds = fdc_flux_and_wind(timestamp, sonicU, sonicV, sonicW, sonicT,
                                               heading, rateX, rateY, rateZ, accX, accY,
                                               accZ, lat)

            where

        fluxes = 3 element tuple of L2 numpy array flux products:
            fluxmom_u = along-wind component of momentum flux FLUXMOM-U_L2 [m^2/s^2]
            fluxmom_v = cross-wind component of momentum flux FLUXMOM-V_L2 [m^2/s^2]
            fluxhot = sonic buoyancy flux FLUXHOT_L2 [m/s * K]
        windspeeds = 3 element tuple of L1 windspeed product lists:
            windtur_vln = windspeed North WINDTUR-VLN_L1 [m/s]
            windtur_vlw = windspeed West WINDTUR-VLW_L1 [m/s]
            windtur_vlu = windspeed Up WINDTUR-VLU_L1 [m/s]
        timestamp = data date and time values [seconds since 1900-01-01]
        sonicU, sonicV, sonicW = L0 wind velocities from the sonic anemometer
        sonicT = L0 sonic temperature from the sonic anemometer
        heading = L0 variable from the magnetometer (not gyro)
        rateX, rateY, rateZ = L0 angular rates
        accX, accY, accZ = L0 linear accelerations
        lat = latitude of instrument in decimal degrees

    Notes:

        This routine is directly called by the functions calculating the flux and wind
        data products. It is anticipated that all of the input variables with the possible
        exception of latitude will be 1D arrays. The function fdc_quantize_data parses
        the data in these 1D arrays into the discrete dataset packets necessary for
        calculation of the L1 and L2 data products.

        The pitch and roll variables are not used in any of the FDCHP calculations.
        However, the machinery to process these L0 variables is kept (1) for future
        use (2) for users downloading this code.

    References:

        OOI (2014). Data Product Specification for FDCHP Data Products. Document
            Control Number 1341-00280. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00280_Data_Product_Spec_FDCHP_OOI.pdf)
    """
    # uncertainty in how latitude will be broadcasted. so.
    lat = np.atleast_1d(lat)
    if lat.size == 1:
        lat = np.repeat(lat, sonicT.size)

    # condition data and parse it into discrete datasets.
    # the heading data is passed 2 extra times, once to take the place of pitch
    # data, once for roll, which in the original matlab code are processed but
    # not used to calculate any data products.
    data = fdc_quantize_data(timestamp, sonicU, sonicV, sonicW, sonicT, heading,
                             heading, heading, rateX, rateY, rateZ, accX, accY,
                             accZ, lat)

    # the shape of the data array is [n_var, n_packets, pts per packet].
    #print data.shape
    n_pack = data.shape[1]

    # for clarity in following the DPS code, unpack the data into its constituent
    # variables as 2D arrays so that the index of the lead dimension of each
    # variable array indicates the dataset number; sonicU[0, :] will be a vector
    # containing the U wind data for the first dataset packet.
    sonicU = data[1, :, :]
    sonicV = data[2, :, :]
    sonicW = data[3, :, :]
    sonicT = data[4, :, :]
    heading = data[5, :, :]
    roll = data[6, :, :]      # not used to calculate any data products
    pitch = data[7, :, :]     # not used to calculate any data products
    rateX = data[8, :, :]
    rateY = data[9, :, :]
    rateZ = data[10, :, :]
    accX = data[11, :, :]
    accY = data[12, :, :]
    accZ = data[13, :, :]
    lat = data[14, :, :]

    # pitch and roll aren't currently used. to emphasize this:
    roll = np.nan
    pitch = np.nan

    # calculate the gravitational acceleration for each dataset
    gv = fdc_grv(np.median(lat, axis=-1))

    # process L0 data
    sonicU = 0.01 * sonicU
    sonicV = 0.01 * sonicV
    sonicW = 0.01 * sonicW
    sonicT = 0.01 * sonicT
    sonicT = sonicT * sonicT / 403.0

    # convert IMU from N,E,Down to match Sonic N,W,Up coordinate system
    rateY = -rateY
    rateZ = -rateZ
    accY = -accY
    accZ = -accZ
    pitch = -pitch
    heading = -heading

    # hardcoded variables in the DPS code
    G = 9.80665  # units of accelerometer values
    roffset = 0.0
    poffset = 0.0
    # z distance between IMU and sonic sampling volume: ok to hardcode
    z_imu_2_smplvol = 0.753
    # distance vector between IMU and sonic sampling volume
    Rvec = np.array([0.0, 0.0, z_imu_2_smplvol])

    fs = 10.0                      # sampling frequency, Hz
    #fltr_cutoff_freq = 10.0        # cutoff frequency to generate filter coeffs
    # 16-oct-2014 e-mail from Jim Edson (DPS author): use
    fltr_cutoff_freq = 12.0

    bhi, ahi = fdc_filtcoef(fs, 1.0/fltr_cutoff_freq)

    # gyro is the processed compass data, and,
    # goodcompass is a switch signifying whether these data are good.
    # this subroutine is vectorized.
    gyro, goodcompass = fdc_process_compass_data(heading)

    # number of seconds of data to remove from the beginning and end of the processed
    # data before calculating the mean of the products of the elements of two vectors.
    edge_sec = 30
    # number of edge data values to remove, based on sampling frequency
    edge = fs * edge_sec
    # set up sonic temperature for buoyancy flux calculation;
    # the temperature processing can be vectorized outside the loop
    Ts_L1 = sonicT[:, edge:-edge]
    Ts = fdc_detrend(Ts_L1, -1, 'linear')

    # initialize L2 dataproduct arrays
    fluxmom_u = np.zeros(n_pack)
    fluxmom_v = np.zeros(n_pack)
    fluxhot = np.zeros(n_pack)

    # and lists to contain the L1 dataproducts
    vln = [None] * n_pack
    vlw = [None] * n_pack
    vlu = [None] * n_pack
    Tmp = [None] * n_pack

    # process one datapacket at a time
    for ii in range(n_pack):
        # wind speeds
        sonics = np.vstack((sonicU[ii, :], sonicV[ii, :], sonicW[ii, :]))

        # process angular rate data; already in radians
        deg_rate = np.vstack((rateX[ii, :], rateY[ii, :], rateZ[ii, :]))
        deg_rate = fdc_despikesimple(deg_rate)

        # process the linear accelerometer data:
        platform = np.vstack((accX[ii, :], accY[ii, :], accZ[ii, :])) * G
        platform = fdc_despikesimple(platform)
        gcomp = np.mean(platform, axis=-1)
        g = np.array([np.sqrt(np.sum(gcomp*gcomp))])
        platform = platform * gv[ii]/g

        platform[0, :] = platform[0, :] + poffset
        platform[1, :] = platform[1, :] + roffset

        gcomp = np.mean(platform, axis=-1)
        g = np.array([np.sqrt(np.sum(gcomp*gcomp))])
        platform = platform * gv[ii] / g

        euler, dr = fdc_anglesclimodeyaw(ahi, bhi, fs, platform, deg_rate,
                                         gyro[ii, :], goodcompass[ii, 0])

        # euler angles are right-handed
        _, uvwplat, _ = fdc_accelsclimode(bhi, ahi, fs, platform, euler)

        uvw, _, _ = fdc_sonic(sonics, dr, euler, uvwplat, Rvec)

        UVW_L1 = uvw[:, edge:-edge]

        # rotate wind velocity components into windstream
        u = fdc_alignwind(UVW_L1)

        u = fdc_detrend(u, -1, 'linear')

        # calculate flux products
        fluxmom_u[ii] = np.mean(u[2, :] * u[0, :])
        fluxmom_v[ii] = np.mean(u[2, :] * u[1, :])
        fluxhot[ii] = np.mean(u[2, :] * Ts[ii])

        # save the L1 wind data products
        (vln[ii], vlw[ii], vlu[ii]) = (UVW_L1[0, :], UVW_L1[1, :], UVW_L1[2, :])

    fluxes = (fluxmom_u, fluxmom_v, fluxhot)

    windspeeds = (vln, vlw, vlu)

    return fluxes, windspeeds


"""
#...................................................................................
#...................................................................................
    Subroutines called by the primary routine fdc_flux_and_wind and its subroutines:

        fdc_accelsclimode
        fdc_alignwind
        fdc_anglesclimodeyaw
        fdc_despikesimple
        fdc_detrend
        fdc_filtcoef
        fdc_grv
        fdc_process_compass_data
        fdc_quantize_data
        fdc_sonic
        fdc_trans
        fdc_update
#...................................................................................
#...................................................................................
"""


def fdc_accelsclimode(bhi, ahi, sf, accm, euler):
    """
    Description:

        Rotates linear accelerations measured on the FDCHP platform into an earth reference
        system, then integrates to get platform velocity and displacement. This code is
        straight from the FDCHP DPS.

    Implemented by:

        2014-05-19: Russell Desiderio. Initial Code

    Usage:

        acc, uvwplat, xyzplat = fdc_accelsclimode(bhi, ahi, sf, accm, euler)

            where

        acc = (3xn_rec) linear accelerations in "FLIP/Earth" reference (NOT USED IN DPA)
        uvwplat = (3xn_rec) linear velocities at the point of measurement
        xyzplat = (3xn_rec) platform displacements from mean position (NOT USED IN DPA)
        bhi = numerator coefficients for high pass filter
        ahi = denominator coefficients for high pass filter
        sf = sampling frequency
        accm = measured platform linear accelerations
        euler = (3Xn_rec) euler angles phi, theta, psi

    References:

        OOI (2014). Data Product Specification for FDCHP Data Products. Document
            Control Number 1341-00280. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00280_Data_Product_Spec_FDCHP_OOI.pdf)
    """
    # keep DPS variable names and comments
    gravxyz = np.mean(accm, axis=-1)
    gravity = np.sqrt(gravxyz.dot(gravxyz))

    # rotate measured accelerations into earth frame
    acc = fdc_trans(accm, euler)
    acc[2, :] = acc[2, :] - gravity

    # integrate accelerations to get velocities
    uvwplat = integrate.cumtrapz(acc, axis=-1, initial=0.0) / sf
    uvwplat = signal.filtfilt(bhi, ahi, uvwplat, axis=-1, padtype='odd', padlen=15)

    # integrate again to get displacements
    xyzplat = integrate.cumtrapz(uvwplat, axis=-1, initial=0.0) / sf
    xyzplat = signal.filtfilt(bhi, ahi, xyzplat, axis=-1, padtype='odd', padlen=15)

    return acc, uvwplat, xyzplat


def fdc_alignwind(u):
    """
    Description:

        Rotates wind velocity components into the streamwise wind. This code is straight from
        the FDCHP DPS.

    Implemented by:

        2014-05-19: Russell Desiderio. Initial Code; converted arithmetic to matrix multiplication.
        2014-05-29: Russell Desiderio. Programmed original arithmetic, faster than matrix operations.

    Usage:

        u_rot = fdc_alignwind(u)

            where

        u_rot = (3xn_rec) wind velocity rotated into the "streamwise" coordinate system
        u = (3Xn_rec) wind velocity in the earth frame (uncorrected for magnetic declination)

    References:

        OOI (2014). Data Product Specification for FDCHP Data Products. Document
            Control Number 1341-00280. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00280_Data_Product_Spec_FDCHP_OOI.pdf)
    """
    # mean wind velocity components
    u_mean = np.mean(u, axis=-1)

    # calculate angles for coordinate rotation
    u_hor = np.sqrt(u_mean[0] * u_mean[0] + u_mean[1] * u_mean[1])
    beta = np.arctan2(u_mean[2], u_hor)
    alpha = np.arctan2(u_mean[1], u_mean[0])

    # populate rotation matrix
    sin_a = np.sin(alpha)
    cos_a = np.cos(alpha)
    sin_b = np.sin(beta)
    cos_b = np.cos(beta)

    ur = u[0, :] * cos_a * cos_b + u[1, :] * sin_a * cos_b + u[2, :] * sin_b
    vr = -u[0, :] * sin_a + u[1, :] * cos_a
    wr = -u[0, :] * cos_a * sin_b - u[1, :] * sin_a * sin_b + u[2, :] * cos_b

    u_rot = np.vstack((ur, vr, wr))

    return u_rot


def fdc_anglesclimodeyaw(ahi, bhi, sf, accm, ratem, gyro, goodcompass):
    """
    Description:

        Calculates the euler angles for the FDCHP instrument platform. This code is straight from
        the FDCHP DPS.

    Implemented by:

        2014-05-19: Russell Desiderio. Initial Code
        2014-09-25: Russell Desiderio. Incorporated latest changes to pitch and roll calculation.

    Usage:

        euler, dr = fdc_anglesclimodeyaw(ahi, bhi, sf, accm, ratem, gyro, goodcompass)

            where

        euler = (3xn_rec) array of euler angles (phi, theta, psi) in radians
        dr = corrected angular rate velocities
        ahi = denominator coefficients for high pass filter
        bhi = numerator coefficients for high pass filter
        sf = sampling frequency
        accm = (3xn_rec) array of recalibrated platform linear accelerations
        ratem = (3xn_rec) array of recalibrated angular rates
        gyro = (1xn_rec) array of gyro signal (= heading = compass)
        goodcompass = boolean signifying whether gyro measurements are to be used.

    References:

        OOI (2014). Data Product Specification for FDCHP Data Products. Document
            Control Number 1341-00280. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00280_Data_Product_Spec_FDCHP_OOI.pdf)
    """
    # hardcoded number of iterations for integrations
    n_iterations = 5

    # keep DPS variable names and comments
    gravxyz = np.mean(accm, axis=-1)
    gravity = np.sqrt(np.sum(gravxyz*gravxyz))

    # unwrap compass
    gyro = np.unwrap(gyro)

    # DPS comment: "remove mean from rate sensors".
    # However, the DPS code removes the linear trend.
    # note that the matlab function detrend used with the 'linear'
    # option specified subtracts piecewise linear values, whereas
    # the scipy version acts to subtract least fit values (as does
    # the matlab version with no options specified).
    ratem = fdc_detrend(ratem, -1, 'linear')

    ### documentation from DPS code verbatim:
    # low frequency angles from accelerometers and gyro
    # slow roll from gravity effects on horizontal accelerations. low pass
    # filter since high frequency horizontal accelerations may be 'real'

    ## avoid imaginary numbers by making sure that arcsin operates on [-1 1].
    # RAD version of revised DPS pitch
    # also, trap out runtime warnings when theta has nans
    theta = np.array(-accm[0, :] / gravity, ndmin=2)
    nanmask = np.isnan(theta)
    theta[nanmask] = 0.0
    theta[theta <= -1.0] = -np.pi/2.0
    theta[theta >= 1.0] = np.pi/2.0
    mask = np.absolute(theta) < 1.0
    theta[mask] = np.arcsin(theta[mask])
    theta[nanmask] = np.nan
    theta_slow = theta - signal.filtfilt(bhi, ahi, theta, axis=-1, padtype='odd', padlen=15)

    # RAD version of revised DPS roll
    # nanmask used for the same reason as above
    phi = np.array(accm[1, :] / gravity, ndmin=2) / np.cos(theta_slow)
    nanmask = np.isnan(phi)
    phi[nanmask] = 0.0
    phi[phi <= -1.0] = -np.pi/2.0
    phi[phi >= 1.0] = np.pi/2.0
    mask = np.absolute(phi) < 1.0
    phi[mask] = np.arcsin(phi[mask])
    phi[nanmask] = np.nan
    phi_slow = phi - signal.filtfilt(bhi, ahi, phi, axis=-1, padtype='odd', padlen=15)

    ### documentation from DPS code verbatim:
    # yaw
    # here, we estimate the slow heading. the 'fast heading' is not needed
    # for the euler angle update matrix. the negative sign puts the gyro
    # signal into a right handed system.

    # these are from revised DPS code: fs = 10, cutoff_freq = 1/240
    #ahi2 = np.array([1.0, -4.989457431527359, 9.957885277614746, -9.936911062858101,
    #                 4.957996019344925, -0.989512802573847])
    #bhi2 = np.array([0.994742581059968, -4.973712905299841, 9.947425810599682, -9.947425810599682,
    #                 4.973712905299841, -0.994742581059968])

    # Using 1.0/240.0 as the second argument results in matrices approaching singularity
    # (computation on the edge of robustness).
    bhi2, ahi2 = fdc_filtcoef(sf, 1.0/240.0)

    if goodcompass:
        psi_slow = -gyro - signal.filtfilt(bhi2, ahi2, -gyro, axis=-1, padtype='odd', padlen=15)
    else:
        psi_slow = -np.median(gyro)*np.ones(phi.shape)

    # use slow angles as first guess
    euler = np.vstack((phi_slow, theta_slow, psi_slow))
    rates = fdc_update(ratem, euler)

    # "i will use this filter with a lower cutoff for yaw
    #  since the compass is having issues"

    # integrate and filter angle rates, and add to slow angles
    for ii in range(n_iterations):
        phi_int = integrate.cumtrapz(rates[0, :], axis=-1, initial=0.0) / sf
        phi = phi_slow + signal.filtfilt(bhi, ahi, phi_int, axis=-1, padtype='odd', padlen=15)
        theta_int = integrate.cumtrapz(rates[1, :], axis=-1, initial=0.0) / sf
        theta = theta_slow + signal.filtfilt(bhi, ahi, theta_int, axis=-1, padtype='odd', padlen=15)
        psi_int = integrate.cumtrapz(rates[2, :], axis=-1, initial=0.0) / sf
        # rad: note that psi_slow values are also a function of the goodcompass value
        if goodcompass:
            psi = psi_slow + signal.filtfilt(bhi2, ahi2, psi_int, axis=-1, padtype='odd', padlen=15)
        else:
            psi = psi_slow + psi_int

        euler = np.vstack((phi, theta, psi))
        rates = fdc_update(ratem, euler)
        rates = fdc_detrend(rates, -1, 'constant')

    dr = ratem

    return euler, dr


def fdc_despikesimple(data):
    """
    Description:

        Function to remove outliers.
        This function is a weak point in the FDCHP programming.

    Implemented by:

        2014-05-19: Craig Risien. Initial Code
        2014-05-30: Russell Desiderio. Vectorized some of the code.
        2014-11-19: Russell Desiderio. Made code more robust, so that out of range
                    interpolations are returned as nans to avoid execution runtime
                    errors.

    Usage:

        [data] = fdc_despikesimple(data)

            where

        data = (3xN) array of data values

    Notes:

        By default, nearest neighbor interpolation causes a ValueError to be raised
        in instances where extrapolation would be necessary to replace outliers. This
        is trapped out by setting bounds_error to False; then, to match the testcode
        (Matlab) output, a fill value of np.nan (which is also the default) is used.
        This modification also requires the use of nanmedian and nanstd functions,
        else it is possible that the input arrays to interpolate.interp1d be empty.

        A more robust method to find outliers is needed, one that will recognize and
        not include the low frequency trend/variability of the data when specifying data
        points as spikes; some form of a Tukey-style implementation may be appropriate.

        When the data abscissae are equidistant,
            Python Scipy nearest neighbor interpolation replaces from the left.
            Matlab interp1 nearest neighbor interpolation replaces from the right.

    References:

        OOI (2014). Data Product Specification for FDCHP Data Products. Document
            Control Number 1341-00280. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00280_Data_Product_Spec_FDCHP_OOI.pdf)
    """
    # at first the matlab test code flux results for the goodcompass=0 case differed
    # in the 4th decimal place compared to the python results; this agreement extends
    # out to the 11th decimal place if the data is flipped before and after this routine
    # is run. (the rest of the discrepancy results because filtfilt.m and scipy.filtfilt
    # do not give exactly the same results).

    # the reason flipping the python data increases agreement is because matlab
    # nearest neighbor interpolation replaces from the right while scipy replaces
    # from the left.
    data = np.fliplr(data)  # to match the matlab results

    # number of times to run each vector stream through the despiking routine
    n_iterations = 3
    # standard deviation span; this was 6 in the original DPS and revised code;
    # Jim Edson (DPS author) says to set this at 4
    n_std = 4
    # interpolation method; 'nearest' in test code
    ntrpmeth = 'nearest'
    #print "interpolation method", ntrpmeth

    array_size = np.shape(data)
    t = np.arange(0, array_size[1])

    for jj in range(n_iterations):

        # vectorize the median, stdev, and masking operations outside of the inner
        # loop, which should help program efficiency when processing large numbers
        # of datasets.

        # calculate the median and stdev as column vectors for broadcasting
        M = np.atleast_2d(sp.stats.nanmedian(data, axis=-1)).T
        # original code used matlab std function, which has a "N-1" in denominator -
        # so, ddof=1.
        Sn = np.nanstd(data, axis=-1, ddof=1, keepdims=True) * n_std
        mask = np.logical_and(data < M + Sn, data > M - Sn)
        # the interp1d function is vectorized for the second argument 2D array ONLY IF
        # the first argument 1D array is unchanging - which it's not.
        # therefore, here a for loop is required.
        for ii in range(array_size[0]):
            f = interpolate.interp1d(t[mask[ii, :]], data[ii, mask[ii, :]], kind=ntrpmeth,
                                     bounds_error=False, fill_value=np.nan)
            data[ii, :] = f(t)

        ## as coded in DPS
        #for tot in range(array_size[0]):
        #    M = np.median(data[tot, :])
        #    # original code used matlab std function, which has a "N-1" in denom, so ddof=1
        #    stan = np.std(data[tot, :], ddof=1)
        #    Sn = stan * n_std
        #    mask = np.logical_and(data[tot, :] < M + Sn, data[tot, :] > M - Sn)
        #    f = interpolate.interp1d(t[mask], data[tot, mask], kind=ntrpmeth)
        #    data[tot, :] = f(t)

    # re-orient the data to the way it came into the routine
    data = np.fliplr(data)

    return data


def fdc_detrend(data, axis=-1, type='linear'):
    """
    Description:

        Calls scipy.signal.detrend to detrend data. fdc_detrend was written to trap
        out nan values because scipy.signal.detrend throws runtime execution errors
        if nans are encountered. In these cases the function fdc_detrend sets *all*
        output values to nan, regardless of the dimensionality of data and what axis
        is designated.

    Implemented by:

        2014-11-19: Russell Desiderio. Initial Code

    Usage:

        detrended_data = fdc_detrend(data, axis, type)

            where

        axis = data axis along which detrend will be applied; default is last axis (-1)
        type = detrend method, either 'linear' or 'constant'.

    References:

        Scipy documentation for scipy.signal.detrend.
    """
    # trap out nans and infs
    if np.any(~np.isfinite(data)):
        detrended = np.zeros(data.shape) + np.nan
        return detrended
    # else run the data through the scipy function
    detrended = signal.detrend(data, axis=axis, type=type)
    return detrended


def fdc_filtcoef(fs, fc):
    """
    Description:

        Function to calculate Butterworth filter coefficients as a function of
        sampling frequency and cutoff frequency.

    Implemented by:

        2014-09-25: Russell Desiderio. Initial code.

    Usage:

        bhigh, ahigh = fdc_filtcoef(fs, fc)

            where

        bhigh, ahigh = Butterworth filter coefficients (used in function filtfilt)
        fs = sampling frequency in Hz
        fc = either the cutoff frequency or its inverse (?)

    Notes:

        Edson documentation:

        DEFINE A HIGH PASS FILTER WHICH RETAINS REAL ACCELERATION BUT REMOVES DRIFT
        3/2/98 - MODIFIED FILTER DESIGN.  THE OLD FILTER WAS

          wp=1/40/(sf/2);ws=1/20/(sf/2);
          [n,wn]=buttord(wp,ws,1,14);

        THE NEW FILTER WAS CHOSEN SO THAT THE SPECTRA OF THE DOUBLY INTEGRATED ACCELERATION
        MATCHED THE FREQUENCY DOMAIN INTEGRATED POWER SPECTRUM IN THE PASS BAND. THE
        COMPARISON WAS MOST SENSITIVE TO THE TRANSITION WIDTH

        THE NEW FILTER WAS CHOSEN TO HAVE A TRANSITION REGION WHICH LIES
        IN THE OVERLAP REGION BETWEEN THE INTEGRATED ANGLE RATE AND THE
        ACCELEROMETER BASED ANGLE ESTIMATES. THE NEW FILTER IS SHIFT
        TO HIGHER FREQUENCIES BY ABOUT 1/4 TO 1/2 DECADE

    References:

        Scipy.signal documentation for the functions buttord and butter.
    """
    nfreq = fs / 2.0
    wp = fc / nfreq
    ws = .7 * wp
    n, wn = signal.buttord(wp, ws, 10.0, 25.0)
    bhigh, ahigh = signal.butter(n, wn, 'high')
    return bhigh, ahigh


def fdc_grv(lat):
    """
    Description:

        Calculates gravity (acceleration due to earth's gravitational field)
        as a function of latitude. This code is from the FDCHP DPS.

    Implemented by:

        2014-05-15: Russell Desiderio. Initial Code

    Usage:

        g = fdc_grv(lat)

            where

        g = acceleration due to earth's gravitational field [m/s/s]
        lat = latitude of instrument in decimal degrees.

    References:

        OOI (2014). Data Product Specification for FDCHP Data Products. Document
            Control Number 1341-00280. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00280_Data_Product_Spec_FDCHP_OOI.pdf)
    """
    # constants from the DPS:
    # equatorial value for 'g'
    gamma = 9.7803267715
    # coefficients of polynomial in (sin(lat))^2
    c1 = 0.0052790414
    c2 = 0.0000232718
    c3 = 0.0000001262
    c4 = 0.0000000007

    x = sp.sin(np.radians(lat))
    xsq = x * x

    # Horner's method for calculating polynomials
    g = gamma * (1.0 + xsq * (c1 + xsq * (c2 + xsq * (c3 + xsq * c4))))

    ## straightforward powers method
    #g=gamma*(1.0+c1*x**2+c2*x**4+c3*x**6+c4*x**8)

    return g


def fdc_process_compass_data(heading):
    """
    Description:

        Vectorized routine which returns processed compass data from the FDCHP
        instrument and a switch denoting how the data is to be processed.

    Implemented by:

        2014-09-25: Russell Desiderio. Separated out from main code.

    Usage:

        gyro, goodcompass = fdc_process_compass_data(heading)

            where

        gyro = processed compass data in an array with the same shape as the input variable.
        goodcompass = switch denoting reliability of heading data: (0,1) = (bad,good).
        heading = heading in a (N,W,Up) coordinate system; this must be a 2D array
                  shaped as [n_packets, n_pts].
    """
    # this routine is vectorized as noted above.
    # number of values on either edge of the compass readings to overwrite:
    edge_compass = 10

    # gyro = 'heading' = compass (may be in N,W,Up cpoordinate sytem); already in radians
    gyro = heading
    # overwrite edge values
    gyro[:, 0:edge_compass] = gyro[:, [edge_compass]]     # right side is a column vector
    gyro[:, -edge_compass:] = gyro[:, [-edge_compass-1]]  # right side is a column vector

    # process gyro values
    gx = np.cos(gyro)
    gy = np.sin(gyro)
    gx = fdc_despikesimple(gx)
    gy = fdc_despikesimple(gy)
    gyro = np.arctan2(gy, gx)
    gyro[gyro < 0] = gyro[gyro < 0] + 2.0 * np.pi

    # determine whether gyro data is good
    gchk = np.unwrap(gyro)
    # matlab std uses (N-1) in the denominator, so set ddof=1
    stdhdg = np.std(gchk, axis=-1, ddof=1, keepdims=True)
    hdg_range = np.amax(gchk, axis=-1, keepdims=True) - np.amin(gchk, axis=-1, keepdims=True)

    # set the goodcompass vector:
    #    if ( hdg_range>(120/180*pi) or stdhdg>(45/180*pi) )
    #        goodcompass = 0
    #    else
    #        goodcompass = 1
    #    end
    #
    # do the same thing without conditional
    goodcompass = np.logical_not(np.logical_or(
        hdg_range > (120.0/180.0*np.pi), stdhdg > (45.0/180.0*np.pi)))

    #print 'goodcompass value: ', goodcompass
    return gyro, goodcompass


def fdc_quantize_data(*args):
    """
    Description:

        Groups data from the FDCHP instrument into discrete packets of npts=12000 records
        each by parsing the records' timestamps. FDCHP is set up to acquire data at 10
        Hz for 20 minutes every hour, so that there will be about 40 minutes in between
        the first record of a dataset packet and the last record of the preceding dataset
        packet.

        Instructions solicited from the DPS author Jim Edson: If more than 12000 records
        are collected (as was true in the testset), truncate to 12000. If less than 12000
        records are collected, pad to 12000 using the last data record.

    Implemented by:

        2014-11-06: Russell Desiderio. Initial code.

    Usage:

        data = fdc_quantize_data(args)

            where

        data = 3D array of input argument variables grouped into datasets. The shape
               of the array is [n_variables, n_packets, n_pts_per_packet].
        args = argument list (tuple) of input data. The number of elements in the list
                  is n_variables; the first element must be (OOI CI) timestamps. Each
                  element must be a 1D array with the same number of components, which
                  will be on the order of npts * n_packets.

    Notes:

        This routine will also work if the only input variable is a 1D array of timestamps.

        The DPS code was written to process 1 dataset of 12000 records. Therefore, given the
        anticipated form of the L0 data inputs into the DPAs (as 1D arrays not divided up into
        data packets), it was necessary to write this function. It might be beneficial to build
        into it a feature deleting data chunks that are less than some fixed number of records
        as not containing enough data to reliably calculate data product values.

    """
    # target number of datapoints per dataset: 10Hz for 20 minutes.
    npts = 12000
    # time in between dataset chunks (last of n_th and first of n+1_th)
    # is expected to be 40 minutes = 2400 seconds; use a lower value as
    # the time discriminant.
    time_gap = 1800.0

    data = np.atleast_2d(args)  # data is a 2D array

    # parse into discrete datasets by finding the number of datapoints in each chunk.
    # data are expected to come in chunks of 20 minutes duration, separated by 40 min.
    # the first row of data must be the timestamps.
    # first find the indices at these timegaps
    idx_at_gap = np.where(np.diff(data[0, :]) > time_gap)[0]

    # prepend and append values to get accurate counts for 1st and last dataset
    idx_at_gap = np.hstack((-1, idx_at_gap, data.shape[1]-1))

    # difference to get answer
    chunklengths = np.diff(idx_at_gap)
    #print chunklengths

    # process one dataset chunk at a time.
    # process all data streams for each dataset at the same time,
    # inspecting each chunk to see if it has (less than), (equal to),
    # or (greater than) npts and then processing it accordingly.
    for ii in range(chunklengths.size):
        if chunklengths[ii] < npts:
            # pad dataset with last set of datapoints.
            # here idx is the insertion point index (1st element = 1);
            # pad data and insert after this point and before start of next dataset
            idx = ii * npts + chunklengths[ii]
            # number of rows to insert
            n_nsrt = npts-chunklengths[ii]
            # now go to python array indexing conventions;
            # tile the column vector data[:,idx-1:idx]
            filldata = np.tile(data[:, idx-1:idx], (1, n_nsrt))
            # correct the timestamps of the filldata.
            delta_time = np.median(np.diff(data[0, idx-chunklengths[ii]:idx]))
            filldata[0, :] = filldata[0, :] + np.arange(1.0, n_nsrt+1) * delta_time
            # and insert filldata into the data array
            data = np.hstack((data[:, 0:idx], filldata, data[:, idx:]))
        elif chunklengths[ii] == npts:
            continue  # no action needed
        else:  # chunklengths[ii] > npts
            # delete data after npts in this dataset and before next dataset.
            idx_del_beg = (ii + 1) * npts
            idx_del_end = ii * npts + chunklengths[ii]
            data = np.delete(data, np.s_[idx_del_beg:idx_del_end], 1)

    # convert data to a 3D array so that calling program can parse its shape
    # to figure out dataset dimensions (n_var, n_dataset_packets, npts per dataset)
    data = np.reshape(data, (data.shape[0], -1, npts))

    return data


def fdc_sonic(sonics, omegam, euler, uvwplat, dist_vec):
    """
    Description:

        This function, which comes from the EDDYCORR toolbox,
        corrects the sonic anemometer components for platform
        motion and orientation.

    Implemented by:

        2014-05-16: Craig Risien. Initial Code

    Usage:

        [uvw, uvwr, uvwrot] = fdc_sonic(sonics, omegam, euler, uvwplat, dist_vec)

            where

        uvw = (MxN) array of corrected sonic anemometer components, in the fixed
               earth frame (North_West-Up)
        uvwr = (NOT USED IN DPA)
        uvwrot = (NOT USED IN DPA)
        sonics = row of integers corre to sonic numbers which are to be corrected
        omegam = (3xN) measured angular rate 'vector' in platform frame
        euler = (3xN) array of euler angles (phi, theta, psi)
        uvwplat = (3xN) array of platform velocities
        dist_vec = (3x1) distance vector between IMU and Sonic sampling volume

    References:

        OOI (2014). Data Product Specification for FDCHP Data Products. Document
            Control Number 1341-00280. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00280_Data_Product_Spec_FDCHP_OOI.pdf)
    """

    n_rec = euler.shape[-1]
    Rvec = np.transpose(np.tile(dist_vec, (n_rec, 1)))
    # override default cross product vector axis definition, which is -1
    uvwrot = np.cross(omegam, Rvec, axis=-2)

    uvwr = fdc_trans(sonics + uvwrot, euler)
    uvw = uvwr + uvwplat

    return uvw, uvwr, uvwrot


def fdc_trans(ang_rates, angles):
    """
    Description:

        It is likely this function serves to transorm the input arguments into the
        earth frame of reference.

    Implemented by:

        2014-05-16: Craig Risien. Initial Code.
        2014-05-30: Russell Desiderio. Removed conditional iflag (always True).

    Usage:

        values = fdc_trans(ang_rates, angles)

            where

        values = output matrix of values
        ang_rates = (3xN) array of angular rates.
        angles = (3xN) array of Euler angles phi,theta,psi.

    References:

        OOI (2014). Data Product Specification for FDCHP Data Products. Document
            Control Number 1341-00280. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00280_Data_Product_Spec_FDCHP_OOI.pdf)
    """

    p = angles[0, :]
    t = angles[1, :]
    ps = angles[2, :]

    up = ang_rates[0, :]
    vp = ang_rates[1, :]
    wp = ang_rates[2, :]

    u = (up * sp.cos(t) * sp.cos(ps) + vp * (sp.sin(p) * sp.sin(t) *
         sp.cos(ps) - sp.cos(p) * sp.sin(ps)) + wp * (sp.cos(p) *
         sp.sin(t) * sp.cos(ps) + sp.sin(p) * sp.sin(ps)))
    v = (up * sp.cos(t) * sp.sin(ps) + vp * (sp.sin(p) * sp.sin(t) *
         sp.sin(ps) + sp.cos(p) * sp.cos(ps)) + wp * (sp.cos(p) *
         sp.sin(t) * sp.sin(ps) - sp.sin(p) * sp.cos(ps)))
    w = (up * (-sp.sin(t)) + vp * (sp.cos(t) * sp.sin(p)) + wp *
         (sp.cos(t) * sp.cos(p)))

    values = np.vstack((u, v, w))

    return values


def fdc_update(ang_rates, angles):
    """
    Description:

        This function is derived from the EDDYCORR toolbox.
        It computes the angular update matrix as described in
        Edson et al. (1998) and Thwaites (1995) page 50.

    Implemented by:

        2014-05-16: Craig Risien. Initial Code

    Usage:

        values = fdc_update(ang_rates, angles)

            where

        values = output matrix of values
        ang_rates = (3xN) array of angular rates.
        angles = (3xN) array of Euler angles phi,theta,psi.

    References:

        OOI (2014). Data Product Specification for FDCHP Data Products. Document
            Control Number 1341-00280. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00280_Data_Product_Spec_FDCHP_OOI.pdf)

        Edson et al, 1998: "Direct covariance flux estimates from mobile platforms
            at sea", J. Atmos. Oceanic Tech., 15, 547-562.

    """

    p = angles[0, :]
    t = angles[1, :]
    ps = angles[2, :]

    up = ang_rates[0, :]
    vp = ang_rates[1, :]
    wp = ang_rates[2, :]

    u = up + vp * sp.sin(p) * np.tan(t) + wp * sp.cos(p) * np.tan(t)
    v = 0 + vp * sp.cos(p) - wp * sp.sin(p)
    w = 0 + vp * sp.sin(p) / sp.cos(t) + wp * sp.cos(p) / sp.cos(t)

    values = np.vstack((u, v, w))
    return values
