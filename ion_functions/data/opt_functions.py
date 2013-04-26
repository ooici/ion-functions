#!/usr/bin/env python
"""
@package ion_functions.data.opt_functions
@file ion_functions/data/opt_functions.py
@author Christopher Wingard
@brief Module containing OPTAA related data-calculations.
"""

# Read in the temperature and salinity correction coefficients.
if 'tscor' not in locals():
    from ion_functions.data.opt_functions_tscor import tscor


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
    resistance = 10000. * volts / (4.516 - volts)
    
    # convert resistance to temperature
    a = 0.00093135
    b = 0.000221631
    c = 0.000000125741
    
    degC = (1. / (a + b * np.log(resistance) + c * np.log(resistance)**3)) - 273.15
    
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


def opt_pd_calc(ref, sig, offset, wvlngth, Tbins, Tarray):
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
    
    # raw reference and signal values are imported as 1D arrays, comprised of a
    # single row, reshaped to column. They must be the same length.
    ref = np.atleast_1d(ref).reshape(-1)
    sig = np.atleast_1d(sig)
    
    # the offsets are imported as a 1D array, need to make sure it is a "column
    # vector". They must be the same length as ref and sig.
    
    # the wavelengths are imported as a 1D array, need to make sure it is a
    # "column vector". they must be the same length as ref, sig and offsets.
    
    # Tint is a scalar
    
    # the Tbins are imported as a 1D array, need to make sure it is a row
    # vector.
    
    # the Tarray, is a 2D array. The # of columns must equal the length of
    # Tbins. The number of rows must equal the number of wvlngths.
    
    
    
    
