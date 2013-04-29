#!/usr/bin/env python
"""
@package ion_functions.data.opt_functions
@file ion_functions/data/opt_functions.py
@author Christopher Wingard
@brief Module containing OPTAA related data-calculations for calculating beam
    attenuation and the optical absorption.
"""

# wrapper functions to calculate the beam attenuation and optical absorption
# from the WET Labs, Inc. ACS (OPTAA).
def opt_beam_attenuation(cref, csig, traw, cwl, coff, tcal, tbins, tc_arr,
                         T, PS):
    """
    Wrapper function to calculate the L2 beam attenuation coefficient from the
    WET Labs, Inc. ACS instrument.
    """ 
    # calculate the internal instrument temperature [deg_C]
    tintrn = opt_internal_temp(traw)
    
    # calculate the beam attenuation coefficient [m^-1]   
    cpd, deltaT = opt_pd_calc(cref, csig, coff, tintrn, tbins, tc_arr)
    
    # correct the beam attenuation coefficient for temperature and salinity.
    cpd_ts = opt_tempsal_corr('c', cpd, cwl, tcal, T, PS)
    
    # return the temperature and salinity corrected beam attenuation
    # coefficient [m^-1]
    return cpd_ts


def opt_optical_absorption(aref, asig, traw, awl, aoff, tcal, tbins, ta_arr,
                       cpd_ts, cwl, T, PS):
    """
    Wrapper function to calculate the L2 optical absorption coefficient from
    the WET Labs, Inc. ACS instrument.
    """ 
    # calculate the internal instrument temperature [deg_C]
    tintrn = opt_internal_temp(traw)
    
    # calculate the optical absorption coefficient [m^-1]   
    apd, deltaT = opt_pd_calc(aref, asig, aoff, tintrn, tbins, ta_arr)
    
    # correct the optical absorption coefficient for temperature and salinty.
    apd_ts = opt_tempsal_corr('a', apd, awl, tcal, T, PS)

    # correct the optical absorption coefficient for scattering effects
    apd_ts_s = opt_scatter_corr(apd_ts, awl, cpd_ts, cwl)
    
    # return the temperature, salinity and scattering corrected optical
    # absorption coefficient [m^-1] 
    return apd_ts_s


# Functions used in calculating optical absorption and beam attenuation
# coefficients from the OPTAA family of instruments.
def opt_pressure(praw, offset, sfactor):
    """
    Description:

        Calculates the pressure (depth) of the ACS, if the unit is equipped
        with a pressure sensor.
        
    Implemented by:

        2013-04-25: Christopher Wingard. Initial implementation.
        
    Usage:

        depth = opt_pressure(praw, offset, sfactor)
        
            where
            
        depth = depth of the instrument [m]
        offset = depth offest from instrument device file [m]
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


def opt_internal_temp(traw):
    """
    Description:

        Calculates the internal instrument temperature. Used in subsequent
        OPTAA calculations.
        
    Implemented by:

        2013-04-25: Christopher Wingard. Initial implementation.
        
    Usage:

        tintrn = opt_internal_temp(traw)
        
            where
        
        tintrn = calculated internal instrument temperature [deg_C]
        traw = raw internal instrument temperature [counts]
    
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
    import numpy as np
    
    # convert counts to volts
    volts = 5. * traw / 65535.
    
    # calculate the resistance of the thermistor
    res = 10000. * volts / (4.516 - volts)
    
    # convert resistance to temperature
    a = 0.00093135
    b = 0.000221631
    c = 0.000000125741
    
    degC = (1. / (a + b * np.log(res) + c * np.log(res)**3)) - 273.15 
    return degC


def opt_external_temp(traw):
    """
    Description:

        Calculates the external environmental temperature.
        
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


def opt_pd_calc(ref, sig, offset, tintrn, tbins, tarray):
    """
    Description:

        Convert raw reference and signal measurements to scientific units.
        
    Implemented by:

        2013-04-25: Christopher Wingard. Initial implementation.
        
    Usage:

        pd, deltaT = opt_pd_calc(ref, sig, offset, tintrn, tbins, tarray)
        
            where
        
        pd = beam attenuation or optical absorption coefficients [m-1]
        ref = raw reference light measurements [counts]
        sig = raw signal light measurements [counts]
        offset = clear water offsets [m-1]
        tintrn = internal instrument temperature [deg_C]
        tbins = instrument specific temperature calibration bins [deg_C]
        tarray = instrument, wavelength and channel specific temperature
            calibration correction coefficients [m-1] 
    
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
    import numpy as np
    
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
    
    # The temperature array, is a 2D array. The # of "columns" must equal the
    # length of temperature bins. The number of "rows" must equal the number of
    # wavelengths.
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
    dT0 = tarray[:,ind1]
    dT1 = tarray[:,ind2]
    deltaT = dT0 + ((tintrn - T0) / (T1 - T0)) * (dT1 - dT0)
    
    # Compute uncorrected egineering units (m^-1) from the Signal and Reference
    # values.
    eng = 4. * np.log(sig/ref)
    
    # Apply the clean water offsets
    pd = (offset - (1./0.25) * np.log(sig/ref)) - deltaT 

    return pd, deltaT


def opt_tempsal_corr(channel, pd, wlngth, tcal, T, PS):
    """
    Description:

        Apply the wavelength and channel temperature and salinity corrections.
        
    Implemented by:

        2013-04-25: Christopher Wingard. Initial implementation.
        
    Usage:

        pd_ts = opt_tempsal_corr(channel, pd, wlngth, tcal, T, PS)
        
            where
        
        pd_ts = temperature and salinity corrected data [m-1]
        channel = which measurement channel is this? 'a' or 'c'
        pd = uncorrected absoprtion/attenuation data [m-1]
        wlngth = wavelengths at which measurements were made [nm]
        tcal = factory calibration reference temperature [deg_C] 
        T = In situ temperature from co-located CTD [deg_C]
        PS = In situ practical salinity from co-located CTD [unitless]
    
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
    import numpy as np
   
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
       
    # Temperature and practical salinity values must be 1D arrays and must be
    # the same length as each other and pd.
    T = np.atleast_1d(T)
    PS = np.atleast_1d(PS)
    lFlag = len(T) != len(PS) and len(T) != len(pd)
    if lFlag:
        raise ValueError('pd, T and PS arrays must be the same length')

    # apply the temperature and salinity corrections to each wavelength
    pd_ts = np.zeros(nValues)
    for i in range(nValues):
        
        # find the temperature and salinity correction coefficients
        ind = np.nonzero(tscor[:,0]-wlngth[i] == 0)
        if np.atleast_1d(ind).size == 1:
            psi_t = tscor[ind[0][0],1]
            psi_sc = tscor[ind[0][0],2]
            psi_sa = tscor[ind[0][0],3]
        else:
            ind1 = np.nonzero(tscor[:,0]-wlngth[i] < 0)[0][0]    
            ind2 = np.nonzero(wlngth[i]-tscor[:,0] < 0)[0][-1]
            
            wv = wlngth[i]
            wv1 = tscor[ind1,0]; wv2 = tscor[ind2,0]
            pt1 = tscor[ind1,1]; pt2 = tscor[ind2,1]
            psc1 = tscor[ind1,2]; psc2 = tscor[ind2,2]
            psa1 = tscor[ind1,3]; psa2 = tscor[ind2,3]
            
            psi_t = pt1 + ((wv-wv1) / (wv2-wv1)) * (pt2 - pt1)
            psi_sa = psa1 + ((wv-wv1) / (wv2-wv1)) * (psa2 - psa1)
            psi_sc = psc1 + ((wv-wv1) / (wv2-wv1)) * (psc2 - psc1)
        
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
        
    Usage:

        apd_ts_s = opt_scatter_corr(apd_ts, wlngth, cpd_ts, cwlngth, rwlngth)
        
            where
            
        apd_ts_s = optical absorption coefficient corrected for temperature,
            salinity, and light scattering effects [m-1]
        apd_ts = optical absorption coefficient corrected for temperature and
            salinity effects [m-1]
        awlngth = absorption channel wavelengths [nm]
        cpd_ts = beam attenuation coefficient corrected for temperature and
            salinity effects [m-1]
        cwlngth = atteunation channel wavelengths [nm]
        rwlngth = scattering reference wavelength (default is 715) [nm]
        
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
    import numpy as np

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
    
    # apply the scattering corrections
    scat_ratio = aref / (cref - aref)
    apd_ts_s = apd_ts - scat_ratio * (cintrp - apd_ts)
    
    return apd_ts_s