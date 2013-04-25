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

