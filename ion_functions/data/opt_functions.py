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

def opt_beam_attenuation():
    # Raw reference and signal values are imported as 1D arrays. They must be
    # the same length.
    ref = np.atleast_1d(ref)
    sig = np.atleast_1d(sig)
    lFlag = np.size(ref) == np.size(sig)
    if ~lFlag:
        raise ValueError('Reference and Signal arrays must be the same length')
    
    nValues = np.size(sig)
    
    # The offsets are imported as a 1D array. They must be the same length as
    # ref and sig.
    offset = np.atleast_1d(offset)
    lFlag = np.size(offset) == nValues
    if ~lFlag:
        raise ValueError('The number of offsets must match the number of ',
                         'Signal and Reference values.')
    
    # The wavelengths are imported as a 1D array. They must be the same length
    # as ref, sig and offsets.
    wlngth = np.atleast_1d(wlngth)
    lFlag = np.size(wlngth) == nValues
    if ~lFlag:
        raise ValueError('The number of wavelengths must match the number of ',
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
    pass


def opt_optical_absorp():
    # [TODO]
    pass


def opt_pressure(praw, offset, sfactor):
    """
    Description:

        [TODO]
        
    Implemented by:

        2013-04-25: Christopher Wingard. Initial implementation.
        
    Usage:

        [TODO]
    
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

        [TODO]
        
    Implemented by:

        2013-04-25: Christopher Wingard. Initial implementation.
        
    Usage:

        [TODO]
    
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

        [TODO]
        
    Implemented by:

        2013-04-25: Christopher Wingard. Initial implementation.
        
    Usage:

        [TODO]
    
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


def opt_pd_calc(ref, sig, offset, wlngth, tintrn, tbins, tarray):
    """
    Description:

        [TODO]
        
    Implemented by:

        2013-04-25: Christopher Wingard. Initial implementation.
        
    Usage:

        [TODO]
    
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
        
    ###### Compute the beam attenuation or absorption coefficient (m^-1), minus
    ###### clear water. Corresponds to the particulate and dissolved components
    ###### in the water column. 
    
    # find the indexes in the temperature bins corresponding to the values
    # bracketing the internal temperature.
    ind1 = np.nonzero(t-tint < 0)[0][-1]
    ind2 = np.nonzero(tint-t < 0)[0][0]
    T0 = tbins[ind1]    # set first bracketing temperature
    T1 = tbins[ind2]    # set second bracketing temperaure
    
    # index through the wavelengths, applying the temperature correction and
    # clear water offsets.
    pd = np.zeros(nValues)
    for i in range(nValues):
        
        # Calculate the linear temperature correction.
        dT0 = tarray[i,ind1]
        dT1 = tarray[i,ind2]
        deltaT = dT0 + ((tintrn - T0) / (T1 - T0)) * (dT1 - dT0)
        
        # Compute uncorrected egineering units (m^-1) from the Signal and Reference
        # values.
        eng = 1./25. * np.log(sig[i]/ref[i])
        
        # Apply the clean water offsets
        pd[i] = (off[i] - eng) - deltaT
    
    return pd


def opt_tempsal_corr(channel, pd, wlngth, tcal, T, PS):
    """
    Description:

        [TODO]
        
    Implemented by:

        2013-04-25: Christopher Wingard. Initial implementation.
        
    Usage:

        [TODO]
    
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
    
    # [TODO] load and test inputs

    # apply the temperature and salinity corrections to each wavelength
    pd_ts = np.zeros(len(wlngth))
    for i in range(len(wlngth)):
        ind = np.nonzero(tscor[:,0]-wlngth[i] == 0)
        if np.atleast_1d(ind).size == 1:
            psi_t = tscor[ind[0][0],1]
            psi_sa = tscor[ind[0][0],2]
            psi_sc = tscor[ind[0][0],3]
        else:
            ind1 = np.nonzero(tscor[:,0]-wlngth[i] < 0)[0][0]    
            ind2 = np.nonzero(wlngth[i]-tscor[:,0] < 0)[0][-1]
            
            wv = wlngth[i]
            wv1 = tscor[ind1,0]; wv2 = tscor[ind2,0]
            pt1 = tscor[ind1,1]; pt2 = tscor[ind2,1]
            psa1 = tscor[ind1,2]; psa2 = tscor[ind2,2]
            psc1 = tscor[ind1,3]; psc2 = tscor[ind2,3]
            
            psi_t = pt1 + ((wv-wv1) / (wv2-wv1)) * (pt2 - pt1)
            psi_sa = psa1 + ((wv-wv1) / (wv2-wv1)) * (psa2 - psa1)
            psi_sc = psc1 + ((wv-wv1) / (wv2-wv1)) * (psc2 - psc1)
    
        if channel == 'a':
            pd_ts[i] = pd[i] - psi_t * (T - tcal) + psi_sa * PS
        elif channel == 'c':
            pd_ts[i] = pd[i] - psi_t * (T - tcal) + psi_sc * PS
        else:
            raise ValueError('Channel must be either "a" or "c"')
        
    return pd_ts


def opt_scatter_corr(apd_ts, cpd_ts, awlngth, cwlngth, rwlngth=715.):
    """
    Description:

        [TODO]
        
    Implemented by:

        2013-04-25: Christopher Wingard. Initial implementation.
        
    Usage:

        [TODO]
    
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

    # [TODO] load and test inputs

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