#!/usr/bin/env python
"""
@package ion_functions.data.nit_functions
@file ion_functions/data/nit_functions.py
@author Craig Risien
@brief Module containing NIT related data-calculations.
"""

import numpy as np


def ts_corrected_nitrate(cal_temp, wl, eno3, eswa, di, dark_value, ctd_t,
                         ctd_sp, data_in, frame_type, wllower=217, wlupper=240):
    """
    Description:

        This Python code is based on Matlab code
        (NUTNR_Example_MATLAB_Code_20140521_ver_1_00.m) that was
        developed by O.E. Kawka (UW/RSN).

        The code below calculates the Dissolved Nitrate Concentration
        with the Sakamoto et. al. (2009) algorithm that uses the observed
        sample salinity and temperature to subtract the bromide component
        of the overall seawater UV absorption spectrum before solving for
        the nitrate concentration.

        The output represents the OOI L2 Dissolved Nitrate Concentration,
        Temperature and Salinity Corrected (NITRTSC).

    Implemented by:

        2014-05-22: Craig Risien. Initial Code
        2014-05-27: Craig Risien. This function now looks for the light vs
                    dark frame measurements and only calculates nitrate
                    concentration based on the light frame measurements.

    Usage:

        NO3_conc = ts_corrected_nitrate(wllower, wlupper, cal_temp, wl,
                        eno3, eswa, di, dark_value, ctd_t, ctd_sp,
                        data_in,frame_type)

            where

        wllower = Lower wavelength limit for spectra fit.
                  From DPS: 217 nm (1-cm pathlength probe tip) or
                            220 nm (4-cm pathlength probe tip)
        wlupper = Upper wavelength limit for spectra fit.
                  From DPS: 240 nm (1-cm pathlength probe tip) or
                            245 nm (4-cm pathlength probe tip)
        cal_temp = Calibration water temperature value
        wl = (1 x 256) array of wavelength bins
        eno3 = (1 x 256) array of wavelength-dependent nitrate
                extinction coefficients
        eswa = (1 x 256) array of seawater extinction coefficients
        di = (1 x 256) array of deionized water reference spectrum
        dark_value = Dark average scalar value
        ctd_t = (1 x N) array of water temperature values from
                colocated CTD [deg C].
                (see 1341-00010_Data_Product_Spec_TEMPWAT)
        ctd_sp = (1 x N) array of practical salinity values from
                colocated CTD [unitless].
                (see 1341-00040_Data_Product_Spec_PRACSAL)
        data_in = (N x 256) array of nitrate measurement values
                from the UV absorption spectrum data product
                (L0 NITROPT) [unitless]
        NO3_conc = L2 Dissolved Nitrate Concentration, Temperature and
                Corrected (NITRTSC) [uM]
        frame_type = Frame type, either a light or dark measurement.
                This function only uses the light frame measurements.

    References:

        OOI (2014). Data Product Specification for NUTNR Data Products.
            Document Control Number 1341-00620.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00620_Data_Product_Spec_NUTNR_OOI.pdf)
        Johnson, K. S., and L. J. Coletti. 2002. In situ ultraviolet
            spectrophotometry for high resolution and long-term monitoring
            of nitrate, bromide and bisulfide in the ocean. Deep-Sea Res.
            I 49:1291-1305
        Sakamoto, C.M., K.S. Johnson, and L.J. Coletti (2009). Improved
            algorithm for the computation of nitrate concentrations in
            seawater using an in situ ultraviolet spectrophotometer.
            Limnology and Oceanography: Methods 7: 132-143
    """

    data_array_size = np.shape(data_in)
    NO3_conc = np.ones(data_array_size[0])

    # Find wavelength bins that fall between the upper and lower
    # limits for spectra fit
    useindex = np.logical_and(wllower <= wl, wl <= wlupper)

    # subset data so that we only use wavelengths between wllower & wlupper
    WL = wl[useindex]
    ENO3 = eno3[useindex]
    ESWA = eswa[useindex]
    DI = np.array(di[useindex], dtype='float64')  # convert to float64

    # coefficients to equation 4 of Sakamoto et al 2009 that give the
    # absorbance of seasalt at 35 salinity versus temperature
    Asak = 1.1500276
    Bsak = 0.02840
    Csak = -0.3101349
    Dsak = 0.001222

    # ENO3 plus a linear baseline
    subset_array_size = np.shape(ENO3)
    # for the constant in the linear baseline
    Ones = np.ones((subset_array_size[0],), dtype='float64') / 100
    M = np.vstack((ENO3, Ones, WL / 1000)).T

    for i in range(0, data_array_size[0]):

        if frame_type[i] == 'SDB' or frame_type[i] == 'SDF' or frame_type[i] == "NDF":

            # Ignore and fill dark frame measurements
            NO3_conc[i] = -9999999.0

        else:

            SW = np.array(data_in[i, useindex], dtype='float64')  # convert to float64

            # correct each SW intensity for dark current
            SWcorr = SW - dark_value[i]

            # calculate absorbance
            Absorbance = np.log10(DI / SWcorr)

            # now estimate molar absorptivity of seasalt at in situ temperature
            # use Satlantic calibration and correct as in Sakamoto et al. 2009.
            SWA_Ext_at_T = (ESWA * ((Asak + Bsak * ctd_t[i]) / (Asak + Bsak * cal_temp))
                            * np.exp(Dsak * (ctd_t[i] - cal_temp) * (WL - 210.0)))

            # absorbance due to seasalt
            A_SWA = ctd_sp[i] * SWA_Ext_at_T
            # subtract seasalt absorbance from measured absorbance
            Acomp = np.array(Absorbance - A_SWA, ndmin=2).T

            # C has NO3, baseline constant, and slope (vs. WL)
            C = np.dot(np.linalg.pinv(M), Acomp)

            NO3_conc[i] = C[0, 0]

    return NO3_conc
