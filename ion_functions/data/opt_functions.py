#!/usr/bin/env python
"""
@package ion_functions.data.opt_functions
@file ion_functions/data/opt_functions.py
@author Christopher Wingard
@brief Module containing OPTAA and PAR data product algorithms.
"""

import numpy as np
import numexpr as ne


# wrapper function to calculate the beam attenuation coefficients (OPTATTN_L2)
#   from the WET Labs, Inc. ACS (OPTAA).
def opt_beam_attenuation(cref, csig, traw, cwl, coff, tcal, tbins, tc_arr,
                         T, PS):

    """
    Description:

        Wrapper function to calculate the L2 beam attenuation coefficients OPTATTN
        from the WET Labs, Inc. ACS instrument.

    Implemented by:

        2013-04-25: Christopher Wingard. Initial implementation.
        2014-03-06: Russell Desiderio. Reset dimensions of arguments and implemented
                    for loop to handle vectorized input.
        2014-03-07: Russell Desiderio. Added Usage documentation.

    Usage:

        cpd_ts = opt_beam_attenaution(cref, csig, traw, cwl, coff, tcal, tbins,
                                      tc_arr, T, PS)

            where

        cpd_ts = beam attenuation coefficients corrected for temperature and salinity
            (OPTATTN_L2) [m-1]
        cref = raw reference light measurements (OPTCREF_L0) [counts]
        csig = raw signal light transmission measurements (OPTCSIG_L0) [counts]
        traw = raw internal instrument temperature (OPTTEMP_L0) [counts]
        cwl = wavelengths at which the beam attenuation measurements were made [nm].
        coff = pure water offsets for the beam attenuation channels from ACS device
            (calibration) file [m-1].
        tcal = factory calibration reference (pure water) temperature [deg_C].
            supplied by the instrument manufacturer (WETLabs).
        tbins = instrument specific internal temperature calibration bin values from
            ACS device (calibration) file [deg_C].
        tc_arr = instrument, wavelength, and channel ('c') specific internal
            temperature calibration correction coefficients from ACS device
            (calibration) file [m-1].
        T  = TEMPWAT(L1): In situ temperature from co-located CTD [deg_C]
        PS = PRACSAL(L2): In situ practical salinity from co-located CTD [unitless]

    References:

        OOI (2013). Data Product Specification for Optical Beam Attenuation
            Coefficient. Document Control Number 1341-00690.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00690_Data_Product_SPEC_OPTATTN_OOI.pdf)

        OOI (2013). Data Product Specification for Optical Absorption
            Coefficient. Document Control Number 1341-00700.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00700_Data_Product_SPEC_OPTABSN_OOI.pdf)

        OOI (2014). OPTAA Unit Test. 1341-00700_OPTABSN Artifact.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI >>
            >> REFERENCE >> Data Product Specification Artifacts >> 1341-00700_OPTABSN >>
            OPTAA_unit_test.xlsx)

    Notes:

        Variables in the argument list, regardless of the number of dimensions, are assumed
        to be vectorized to contain multiple data packets such that the first dimension
        iterates over data packet number.
    """
    # reset shapes of input arguments
    #    using np.array([], ndmin=#) seems faster than using np.atleast_#d
    cref = np.array(cref, ndmin=2)
    csig = np.array(csig, ndmin=2)
    traw = np.array(traw, ndmin=1)
    cwl = np.array(cwl, ndmin=2)
    coff = np.array(coff, ndmin=2)
    tcal = np.array(tcal, ndmin=1)
    tbins = np.array(tbins, ndmin=2)
    T = np.array(T, ndmin=1)
    PS = np.array(PS, ndmin=1)
    # note, np.atleast_3d appends the extra dimension;
    # np.array using ndmin prepends the extra dimension.
    tc_arr = np.array(tc_arr, ndmin=3)

    # size up inputs
    npackets = cwl.shape[0]
    nwavelengths = cwl.shape[1]
    # initialize output array
    cpd_ts = np.zeros([npackets, nwavelengths])

    for ii in range(npackets):

        # calculate the internal instrument temperature [deg_C]
        tintrn = opt_internal_temp(traw[ii])

        # calculate the uncorrected beam attenuation coefficient [m^-1]
        cpd, _ = opt_pd_calc(cref[ii, :], csig[ii, :], coff[ii, :], tintrn,
                             tbins[ii, :], tc_arr[ii, :, :])

        # correct the beam attenuation coefficient for temperature and salinity.
        cpd_ts_row = opt_tempsal_corr('c', cpd, cwl[ii, :], tcal[ii], T[ii], PS[ii])
        cpd_ts[ii, :] = cpd_ts_row

    # return the temperature and salinity corrected beam attenuation
    # coefficient OPTATTN_L2 [m^-1]
    return cpd_ts


# wrapper function to calculate the optical absorption coefficients (OPTABSN_L2)
#   from the WET Labs, Inc. ACS (OPTAA).
def opt_optical_absorption(aref, asig, traw, awl, aoff, tcal, tbins, ta_arr,
                           cpd_ts, cwl, T, PS, rwlngth=715.):
    """
    Wrapper function to calculate the L2 optical absorption coefficient OPTABSN
    from the WET Labs, Inc. ACS instrument.

        2013-04-25: Christopher Wingard. Initial implementation.
        2014-02-19: Russell Desiderio. Added rwlngth to argument lists, so that
                    a non-default scatter correction wavelength could be passed
                    to function opt_scatter_corr.
        2014-03-06: Russell Desiderio. Reset dimensions of arguments and implemented
                    for loop to handle vectorized input.
        2014-03-07: Russell Desiderio. Added Usage documentation.

    Usage:

        apd_ts_s = opt_optical_absorption(aref, asig, traw, awl, aoff, tcal, tbins,
                                      ta_arr, cpd_ts, cwl, T, PS[, rwlngth])

            where

        apd_ts_s = optical absorption coefficients corrected for temperature, salinity,
            and measurement error due to scattering (OPTABSN_L2) [m-1]
        aref = raw reference light measurements (OPTAREF_L0) [counts]
        asig = raw signal light transmission measurements (OPTASIG_L0) [counts]
        traw = raw internal instrument temperature (OPTTEMP_L0) [counts]
        awl = wavelengths at which the absorption measurements were made [nm].
        aoff = pure water offsets for the absorption channels from ACS device
            (calibration) file [m-1].
        tcal = factory calibration reference (pure water) temperature [deg_C].
            supplied by the instrument manufacturer (WETLabs).
        tbins = instrument specific internal temperature calibration bin values from
            ACS device (calibration) file [deg_C].
        ta_arr = instrument, wavelength, and channel ('a') specific internal
            temperature calibration correction coefficients from ACS device
            (calibration) file [m-1].
        cpd_ts = beam attenuation coefficient corrected for temperature and
            salinity effects (OPTATTN_L2) [m-1], from function opt_beam_attenuation.
        cwl = attenuation channel wavelengths [nm], from ACS device (calibration) file.
        T  = TEMPWAT(L1): In situ temperature from co-located CTD [deg_C]
        PS = PRACSAL(L2): In situ practical salinity from co-located CTD [unitless]
        rwlngth = [optional] user selected scattering correction reference wavelength
            (default = 715) [nm]

    References:

        OOI (2013). Data Product Specification for Optical Beam Attenuation
            Coefficient. Document Control Number 1341-00690.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00690_Data_Product_SPEC_OPTATTN_OOI.pdf)

        OOI (2013). Data Product Specification for Optical Absorption
            Coefficient. Document Control Number 1341-00700.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00700_Data_Product_SPEC_OPTABSN_OOI.pdf)

        OOI (2014). OPTAA Unit Test. 1341-00700_OPTABSN Artifact.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI >>
            >> REFERENCE >> Data Product Specification Artifacts >> 1341-00700_OPTABSN >>
            OPTAA_unit_test.xlsx)

    Notes:

        Variables in the argument list, regardless of the number of dimensions, are assumed
        to be vectorized to contain multiple data packets such that the first dimension
        iterates over data packet number.
    """
    # reset shapes of input arguments
    #    using np.array ndmin=# seems faster than using np.atleast_#d
    aref = np.array(aref, ndmin=2)
    asig = np.array(asig, ndmin=2)
    traw = np.array(traw, ndmin=1)
    awl = np.array(awl, ndmin=2)
    aoff = np.array(aoff, ndmin=2)
    tcal = np.array(tcal, ndmin=1)
    tbins = np.array(tbins, ndmin=2)
    # note, np.atleast_3d appends the extra dimension;
    # np.array using ndmin prepends the extra dimension.
    ta_arr = np.array(ta_arr, ndmin=3)
    cpd_ts = np.array(cpd_ts, ndmin=2)
    cwl = np.array(cwl, ndmin=2)
    T = np.array(T, ndmin=1)
    PS = np.array(PS, ndmin=1)

    # size up inputs
    npackets = awl.shape[0]
    nwavelengths = awl.shape[1]
    # initialize output array
    apd_ts_s = np.zeros([npackets, nwavelengths])

    for ii in range(npackets):

        # calculate the internal instrument temperature [deg_C]
        tintrn = opt_internal_temp(traw[ii])

        # calculate the uncorrected optical absorption coefficient [m^-1]
        apd, _ = opt_pd_calc(aref[ii, :], asig[ii, :], aoff[ii, :], tintrn,
                             tbins[ii, :], ta_arr[ii, :, :])

        # correct the optical absorption coefficient for temperature and salinty.
        apd_ts = opt_tempsal_corr('a', apd, awl[ii, :], tcal[ii], T[ii], PS[ii])

        # correct the optical absorption coefficient for scattering effects
        apd_ts_s_row = opt_scatter_corr(apd_ts, awl[ii, :], cpd_ts[ii, :], cwl[ii, :], rwlngth)
        apd_ts_s[ii, :] = apd_ts_s_row

    # return the temperature, salinity and scattering corrected optical
    # absorption coefficient OPTABSN_L2 [m^-1]
    return apd_ts_s


# Functions used in calculating optical absorption and beam attenuation
# coefficients from the OPTAA family of instruments.
def opt_internal_temp(traw):
    """
    Description:

        Calculates the internal instrument temperature. Used in subsequent
        OPTAA calculations.

    Implemented by:

        2013-04-25: Christopher Wingard. Initial implementation.
        2014-03-07: Russell Desiderio. Reduced calls to np.log.

    Usage:

        tintrn = opt_internal_temp(traw)

            where

        tintrn = calculated internal instrument temperature [deg_C]
        traw = raw internal instrument temperature (OPTTEMP_L0) [counts]

    References:

        OOI (2013). Data Product Specification for Optical Beam Attenuation
            Coefficient. Document Control Number 1341-00690.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00690_Data_Product_SPEC_OPTATTN_OOI.pdf)

        OOI (2013). Data Product Specification for Optical Absorption
            Coefficient. Document Control Number 1341-00700.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00700_Data_Product_SPEC_OPTABSN_OOI.pdf)
    """
    # convert counts to volts
    volts = 5. * traw / 65535.

    # calculate the resistance of the thermistor
    res = 10000. * volts / (4.516 - volts)

    # convert resistance to temperature
    a = 0.00093135
    b = 0.000221631
    c = 0.000000125741

    log_res = np.log(res)
    degC = (1. / (a + b * log_res + c * log_res**3)) - 273.15
    return degC


def opt_pd_calc(ref, sig, offset, tintrn, tbins, tarray):
    """
    Description:

        Convert raw reference and signal measurements to scientific units.

        The calculations for the beam attenuation ('c') and absortion ('a') cases
        are isomorphic; they differ in just the values used for the input arguments.

        The returned values are not final data products.

    Implemented by:

        2013-04-25: Christopher Wingard. Initial implementation.
        2014-02-19: Russell Desiderio. Expanded Usage documentation.

    Usage:

        pd, deltaT = opt_pd_calc(ref, sig, offset, tintrn, tbins, tarray)

            where

        pd = uncorrected beam attenuation or optical absorption coefficients [m-1]
        deltaT = correction due to instrument internal temperature [m-1]
            (this value is returned so that it can be checked in unit tests. it
            is not used in subsequent processing).
        ref = raw reference light measurements (OPTCREF_L0 or OPTAREF_L0, as
            appropriate) [counts]
        sig = raw signal light measurements (OPTCSIG_L0 or OPTASIG_L0, as
            appropriate) [counts]
        offset = clear water offsets from ACS device (calibration) file [m-1].
            use 'c' offsets or 'a' offsets, as appropriate.
        tintrn = internal instrument temperature [deg_C]; output from function
            opt_internal_temp
        tbins = instrument specific internal temperature calibration bin values from
            ACS device (calibration) file [deg_C].
        tarray = instrument, wavelength and channel ('c' or 'a') specific internal
            temperature calibration correction coefficients from ACS device
            (calibration) file [m-1]. use 'c' or 'a' coefficients as appropriate.

    References:

        OOI (2013). Data Product Specification for Optical Beam Attenuation
            Coefficient. Document Control Number 1341-00690.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00690_Data_Product_SPEC_OPTATTN_OOI.pdf)

        OOI (2013). Data Product Specification for Optical Absorption
            Coefficient. Document Control Number 1341-00700.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00700_Data_Product_SPEC_OPTABSN_OOI.pdf)

        OOI (2014). OPTAA Unit Test. 1341-00700_OPTABSN Artifact.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI >>
            >> REFERENCE >> Data Product Specification Artifacts >> 1341-00700_OPTABSN >>
            OPTAA_unit_test.xlsx)
    """
    # Raw reference and signal values are imported as 1D arrays. They must be
    # the same length.
    ref = np.atleast_1d(ref).astype(np.float)
    sig = np.atleast_1d(sig).astype(np.float)
    lFlag = len(ref) != len(sig)
    if lFlag:
        raise ValueError('Reference and Signal arrays must be the same length')

    nValues = len(sig)

    # The offsets are imported as a 1D array. They must be the same length as
    # ref and sig.
    offset = np.atleast_1d(offset)
    lFlag = len(offset) != nValues
    if lFlag:
        raise ValueError('The number of offsets must match the number of ',
                         'Signal and Reference values.')

    # The temperature bins are imported as a 1D array
    tbins = np.atleast_1d(tbins)
    tValues = np.size(tbins)

    # The temperature array tarray is a 2D array. The # of "columns" must equal
    # the length of temperature bins. The number of "rows" must equal the number
    # of wavelengths.
    tarray = np.atleast_2d(tarray)
    r, c = tarray.shape

    if r != nValues:
        raise ValueError('The number of rows in the temperature array must ',
                         'match the number of Signal and Reference values.')

    if c != tValues:
        raise ValueError('The number of columns in the temperature array must ',
                         'match the number of temperature bin values.')

    # find the indexes in the temperature bins corresponding to the values
    # bracketing the internal temperature.
    ind1 = np.nonzero(tbins-tintrn < 0)[0][-1]
    ind2 = np.nonzero(tintrn-tbins < 0)[0][0]
    T0 = tbins[ind1]    # set first bracketing temperature
    T1 = tbins[ind2]    # set second bracketing temperaure

    # Calculate the linear temperature correction.
    dT0 = tarray[:, ind1]
    dT1 = tarray[:, ind2]
    deltaT = dT0 + ((tintrn - T0) / (T1 - T0)) * (dT1 - dT0)

    # Calculate the uncorrected signal [m-1]; the pathlength is 0.25m.
    # Apply the corrections for the clean water offsets (offset) and
    # the instrument's internal temperature (deltaT).
    pd = (offset - (1./0.25) * np.log(sig/ref)) - deltaT

    return pd, deltaT


def opt_tempsal_corr(channel, pd, wlngth, tcal, T, PS):
    """
    Description:

        Apply the wavelength and optical channel temperature and salinity corrections.

    Implemented by:

        2013-04-25: Christopher Wingard. Initial implementation.
        2014-02-19: Russell Desiderio. Expanded Usage documentation.
                    Deleted incorrect requirement that T, PS vector lengths
                        are the same as that of the number of wavelengths.

    Usage:

        pd_ts = opt_tempsal_corr(channel, pd, wlngth, tcal, T, PS)

            where

        pd_ts = temperature and salinity corrected data [m-1]
                case 'c': OPTATTN_L2
                case 'a': intermediate absorption product: will also
                          need to have the scattering correction applied.
        channel = which measurement channel is this? 'c' or 'a'
                'c' denotes beam attenuation [m-1]
                'a' denotes absorption [m-1]
        pd = uncorrected absorption or attenuation data [m-1]
                (from function opt_pd_calc)
        wlngth = wavelengths at which measurements were made [nm].
                from ACS device (calibration) file. use 'c' wavelengths or
                'a' wavelengths as appropriate.
        tcal = factory calibration reference (pure water) temperature [deg_C].
                supplied by the instrument manufacturer (WETLabs).
        T  = TEMPWAT(L1): In situ temperature from co-located CTD [deg_C]
        PS = PRACSAL(L2): In situ practical salinity from co-located CTD [unitless]

    References:

        OOI (2013). Data Product Specification for Optical Beam Attenuation
            Coefficient. Document Control Number 1341-00690.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00690_Data_Product_SPEC_OPTATTN_OOI.pdf)

        OOI (2013). Data Product Specification for Optical Absorption
            Coefficient. Document Control Number 1341-00700.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>

        OOI (2014). OPTAA Unit Test. 1341-00700_OPTABSN Artifact.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI >>
            >> REFERENCE >> Data Product Specification Artifacts >> 1341-00700_OPTABSN >>
            OPTAA_unit_test.xlsx)
            1341-00700_Data_Product_SPEC_OPTABSN_OOI.pdf)
    """
    # load the temperature and salinity correction coefficients table
    from ion_functions.data.opt_functions_tscor import tscor

    # Absorption/attenuation and the wavelength values are imported as 1D
    # arrays. They must be the same length.
    pd = np.atleast_1d(pd)
    wlngth = np.atleast_1d(wlngth)
    lFlag = len(pd) != len(wlngth)
    if lFlag:
        raise ValueError('pd and wavelength arrays must be the same length')

    nValues = np.size(pd)

    # apply the temperature and salinity corrections for each wavelength
    pd_ts = np.zeros(nValues)
    for i in range(nValues):
        # find the temperature and salinity correction coefficients
        psi_t = tscor[wlngth[i]][0]
        psi_sc = tscor[wlngth[i]][1]
        psi_sa = tscor[wlngth[i]][2]

        # apply based on optical channel
        if channel == 'a':
            pd_ts[i] = pd[i] - (psi_t * (T - tcal) + psi_sa * PS)
        elif channel == 'c':
            pd_ts[i] = pd[i] - (psi_t * (T - tcal) + psi_sc * PS)
        else:
            raise ValueError('Channel must be either "a" or "c"')

    return pd_ts


def opt_scatter_corr(apd_ts, awlngth, cpd_ts, cwlngth, rwlngth=715.):
    """
    Description:

        Apply the scattering correction to the temperature and salinity
        corrected optical absorption coefficient.

    Implemented by:

        2013-04-25: Christopher Wingard. Initial implementation.
        2014-02-19: Russell Desiderio. Trapped out potential problems in
                    scat_ratio calculation.

    Usage:

        apd_ts_s = opt_scatter_corr(apd_ts, awlngth, cpd_ts, cwlngth[, rwlngth])

            where

        apd_ts_s = optical absorption coefficient corrected for temperature,
            salinity, and light scattering effects (OPTABSN_L2) [m-1]
        apd_ts = optical absorption coefficient corrected for temperature and
            salinity effects [m-1], from function opt_tempsal_corr.
        awlngth = absorption channel wavelengths [nm] from ACS device (calibration) file.
        cpd_ts = beam attenuation coefficient corrected for temperature and
            salinity effects (OPTATTN_L2) [m-1], from function opt_tempsal_corr.
        cwlngth = attenuation channel wavelengths [nm], from ACS device (calibration) file.
        rwlngth = user selected scattering correction reference wavelength
            (default = 715) [nm]

    References:

        OOI (2013). Data Product Specification for Optical Beam Attenuation
            Coefficient. Document Control Number 1341-00690.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00690_Data_Product_SPEC_OPTATTN_OOI.pdf)

        OOI (2013). Data Product Specification for Optical Absorption
            Coefficient. Document Control Number 1341-00700.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00700_Data_Product_SPEC_OPTABSN_OOI.pdf)

        OOI (2014). OPTAA Unit Test. 1341-00700_OPTABSN Artifact.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI >>
            >> REFERENCE >> Data Product Specification Artifacts >> 1341-00700_OPTABSN >>
            OPTAA_unit_test.xlsx)
    """
    # Absorption and the absorption wavelength values are imported as 1D
    # arrays. They must be the same length.
    apd_ts = np.atleast_1d(apd_ts)
    awlngth = np.atleast_1d(awlngth)
    lFlag = len(apd_ts) != len(awlngth)
    if lFlag:
        raise ValueError('Absorption and absorption wavelength arrays must ',
                         'be the same length')

    # Attenuation and the attenuation wavelength values are imported as 1D
    # arrays. They must be the same length.
    cpd_ts = np.atleast_1d(cpd_ts)
    cwlngth = np.atleast_1d(cwlngth)
    lFlag = len(cpd_ts) != len(cwlngth)
    if lFlag:
        raise ValueError('Attenuation and attenuation wavelength arrays must ',
                         'be the same length')

    # find the the 'a' channel wavelength closest to the reference wavelength
    # for scattering and set the 'a' scattering reference value.
    idx = (np.abs(awlngth-rwlngth)).argmin()
    aref = apd_ts[idx]

    # interpolate the 'c' channel cpd_ts values to match the 'a' channel
    # wavelengths and set the 'c' scattering reference value. 
    cintrp = np.interp(awlngth, cwlngth, cpd_ts)
    cref = cintrp[idx]

    # trap out potential problems in scat_ratio calculation:
    # scat_ratio = aref / (cref - aref).
    # aref must be > 0 AND scat_ratio must be > 0; else, scat_ratio = 0.
    if aref <= 0.0:
        scat_ratio = 0.
    elif cref - aref <= 0.0:
        scat_ratio = 0.
    else:
        scat_ratio = aref / (cref - aref)

    # apply the scattering corrections
    apd_ts_s = apd_ts - scat_ratio * (cintrp - apd_ts)
    return apd_ts_s


# The next 2 functions are not used in calculating optical absorption and beam attenuation
# coefficients from the OPTAA family of instruments. However, some of these instruments
# may be outfitted with an auxiliary pressure sensor and/or external temperature sensor.
#
# opt_pressure is not used in the calculation of the final OPTAA data products.
def opt_pressure(praw, offset, sfactor):
    """
    Description:

        Calculates the pressure (depth) of the ACS, if the unit is equipped
        with an auxiliary pressure sensor.

    Implemented by:

        2013-04-25: Christopher Wingard. Initial implementation.

    Usage:

        depth = opt_pressure(praw, offset, sfactor)

            where

        depth = depth of the instrument [m]
        praw = raw pressure reading [counts]
        offset = depth offset from instrument device file [m]
        sfactor = scale factor from instrument device file [m counts-1]

    References:

        OOI (2013). Data Product Specification for Optical Beam Attenuation
            Coefficient. Document Control Number 1341-00690.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00690_Data_Product_SPEC_OPTATTN_OOI.pdf)

        OOI (2013). Data Product Specification for Optical Absorption
            Coefficient. Document Control Number 1341-00700.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00700_Data_Product_SPEC_OPTABSN_OOI.pdf)
    """
    depth = praw * sfactor + offset
    return depth


# opt_external_temp is not used in the calculation of the final OPTAA data products.
def opt_external_temp(traw):
    """
    Description:

        Calculates the external environmental temperature of the ACS, if the unit
        is equipped with an auxiliary temperature sensor.


    Implemented by:

        2013-04-25: Christopher Wingard. Initial implementation.

    Usage:

        textrn = opt_external_temp(traw)

            where

        textrn = calculated external environment temperature [deg_C]
        traw = raw external temperature [counts]

    References:

        OOI (2013). Data Product Specification for Optical Beam Attenuation
            Coefficient. Document Control Number 1341-00690.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00690_Data_Product_SPEC_OPTATTN_OOI.pdf)

        OOI (2013). Data Product Specification for Optical Absorption
            Coefficient. Document Control Number 1341-00700.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00700_Data_Product_SPEC_OPTABSN_OOI.pdf)
    """
    # convert counts to degrees Centigrade
    a = -7.1023317e-13
    b = 7.09341920e-08
    c = -3.87065673e-03
    d = 95.8241397

    degC = a * traw**3 + b * traw**2 + c * traw + d
    return degC


def opt_par_satlantic(counts_output, a0, a1, Im):
    """
    Description:

        The OOI Level 1 Photosynthetically Active Radiation (PAR)
        (OPTPARW) core data product is the spectral range
        (wavelength) of solar radiation from 400 to 700 nanometers
        that photosynthetic organisms are able to use in the process
        of photosynthesis.

    Implemented by:

        2014-01-31: Craig Risien. Initial Code

    Usage:

        OPTPARW_L1 = opt_par_satlantic(counts_output, a0, a1, Im):

        Calculates the L1 OPTPARW from the Satlantic instrument on the
        RSN Shallow Profiler:
        PAR (umol photons m^-2 s^-1) = Im * a1 (counts_output - a0)

            where

        counts_output is the OPTPARW L0 output [ADC counts]
        a0 is the voltage offset [counts]
        a1 is the scaling factor [umol photons per m^2 per second per count]
        Im = immersion coefficient

    References:

        OOI (2012). Data Product Specification for PHOTOSYNTHETICALLY
        ACTIVE RADIATION (PAR) FROM SATLANTIC INSTRUMENT ON RSN SHALLOW
        PROFILER Document Control Number 1341-00720.
        https://alfresco.oceanobservatories.org/
        (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
        1341-00720_Data_Product_SPEC_OPTPARW_Satl_OOI.pdf)
    """

    OPTPARW_L1 = ne.evaluate('Im * a1 * (counts_output - a0)')

    return OPTPARW_L1


def opt_par_biospherical_mobile(output, dark_offset, scale_wet):
    """
    Description:

        The OOI Level 1 Photosynthetically Active Radiation (PAR)
        (OPTPARW) core data product is the spectral range (wavelength)
        of solar radiation from 400 to 700 nanometers that photosynthetic
        organisms are able to use in the process of photosynthesis.

    Implemented by:

        2014-01-31: Craig Risien. Initial Code

    Usage:

        OPTPARW_L1 = opt_par_biospherical_mobile(output, dark_offset, scale_wet):

        Calculate the L1 OPTPARW from the Biospherical QSP-2100 series of scalar
        instruments.

        where

        OPTPARW_L1 is Level 1 Photosynthetically Active Radiation [umol photons m^-2 s^-1]
        output is the OPTPARW L0 output [volts]
        dark offset is the dark reading [volts]
        scale_wet is the wet calibration scale factor [volts per umol photons / m^2 s^1]

    References:

        OOI (2012). Data Product Specification for PHOTOSYNTHETICALLY
        ACTIVE RADIATION (PAR) FROM BIOSPHERICAL INSTRUMENT QSP 2100
        ON CGSN MOBILE ASSETS Document Control Number 1341-00721.
        https://alfresco.oceanobservatories.org/
        (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
        1341-00721_Data_Product_SPEC_OPTPARW_Bios_OOI.pdf)
    """

    OPTPARW_L1 = ne.evaluate('(output - dark_offset) / scale_wet')

    return OPTPARW_L1


def opt_par_biospherical_wfp(output, dark_offset, scale_wet):
    """
    Description:

        The OOI Level 1 Photosynthetically Active Radiation (PAR)
        (OPTPARW) core data product is the spectral range (wavelength)
        of solar radiation from 400 to 700 nanometers that photosynthetic
        organisms are able to use in the process of photosynthesis.

    Implemented by:

        2014-03-07: Craig Risien. Initial Code

    Usage:

        OPTPARW_L1 = opt_par_biospherical_wfp(output, dark_offset, scale_wet):

        Calculate the L1 OPTPARW from the Biospherical QSP-2200 series of scalar
        instruments.

        where

        OPTPARW_L1 is Level 1 Photosynthetically Active Radiation [umol photons m^-2 s^-1]
        output is the OPTPARW L0 output [millivolts (mV)]
        dark offset is the dark reading [millivolts (mV)]
        scale_wet is the wet calibration scale factor [volts / (quanta / cm^2 s^1)]

    References:

        OOI (2012). Data Product Specification for PHOTOSYNTHETICALLY
        ACTIVE RADIATION (PAR) FROM BIOSPHERICAL INSTRUMENT QSP 2200
        ON CGSN PROFILERS Document Control Number 1341-00721.
        https://alfresco.oceanobservatories.org/
        (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
        1341-00721_Data_Product_SPEC_OPTPARW_Bios_OOI.pdf)
    """

    #Convert output from mvolts to volts
    output_volts = output / 1000.

    #Convert dark_offset from mvolts to volts
    dark_offset_volts = dark_offset / 1000.

    #Convert scale_wet from Volts/(quanta/cm^2.s^1) to Volts/(umol photons/m^2.s^1)
    #1uE/sec/m^2 PAR= 1umole/sec/m^2 PAR = 6.02*10**13 quanta/sec/cm^2 PAR
    scale_wet_converted = scale_wet * (6.02 * 10**13)

    OPTPARW_L1 = ne.evaluate('(output_volts - dark_offset_volts) / scale_wet_converted')

    return OPTPARW_L1


def opt_ocr507_irradiance(counts, offset, scale, immersion_factor):
    """
    Description:

        This function calculates the OOI Level 1 Data Product SPECTIR_L1, which is the
        'vector' (cosine-weighted) downwelling irradiance (Ed) in units of uW/cm^2/nm.
        This data product is produced by the Satlantic OCR-507 multispectral radiometer
        in the SPKIR instrument class.

    Implemented by:

        2014-03-14: Russell Desiderio. Initial Code.
            At this time the DPS for this data product is not finished. The OCR-507
            has 7 wavelength channels, and each data packet coming out of the instrument
            contains one record from each channel. It has not been determined whether
            the raw input data to this function will be in arrays of 7 wavelengths or not.

    Usage:

        Ed = opt_ocr507_irradiance(counts, offset, scale, immersion_factor):

        where

        Ed = downwelling vector irradiance SPECTIR_L1 [uW/cm^2/nm].
        counts = raw downwelling vector irradiance signal SPECTIR_L0 [counts].
        offset = calibration coefficient supplied by the manufacturer.
        scale = calibration coefficient supplied by the manufacturer.
        immersion_factor = calibration coefficient supplied by the manufacturer.

    References:

        OOI (2014). Data Product Specification for Downwelling Spectral Irradiance.
            Document Control Number 1341-00730. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00730__???.pdf)
    """

    # Apply cal coeffs to raw data
    Ed = (counts - offset) * scale * immersion_factor
    return Ed
