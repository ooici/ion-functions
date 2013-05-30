#!/usr/bin/env python
"""
@package ion_functions.data.ctd_functions
@file ion_functions/data/ctd_functions.py
@author Christopher Wingard
@brief Module containing CTD related data-calculations.
"""

# Import Numpy and the TEOS-10 GSW libraries
import numpy as np
from pygsw import vectors as gsw
    
def ctd_sbe16plus_tempwat(t0, a0, a1, a2, a3):
    """
    Description:

        OOI Level 1 Water Temperature data product, which is calculated using
        data from the Sea-Bird Electronics conductivity, temperature and depth
        (CTD) family of instruments. This document is intended to be used by
        OOI programmers to construct appropriate processes to create the L1
        water temperature product. 

    Implemented by:

        2013-04-12: Luke Campbell. Initial Code
        2013-04-12: Christopher Wingard. Minor edits
        2013-05-10: Christopher Wingard. Minor edits to comments.
        
    Usage:

        t = ctd_sbe16plus_tempwat(t0, a0, a1, a2, a3)

            where

        t = sea water temperature (TEMPWAT_L1) [deg_C] 
        t0 = raw temperature (TEMPWAT_L0) [counts]        
        a0 = temperature calibration coefficients
        a1 = temperature calibration coefficients
        a2 = temperature calibration coefficients
        a3 = temperature calibration coefficients

    References:
    
        OOI (2012). Data Product Specification for Water Temperature. Document
            Control Number 1341-00010. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00010_Data_Product_SPEC_TEMPWAT_OOI.pdf)
    """
    mv = (t0 - 524288) / 1.6e7
    r = (mv * 2.9e9 + 1.024e8)/(2.048e4 - mv * 2.0e5)
    t = 1 / (a0 + a1 * np.log(r) + a2 * np.power(np.log(r),2)
           + a3 * np.power(np.log(r),3)) - 273.15
    return t


def ctd_sbe16plus_preswat(p0, t0, ptempa0, ptempa1, ptempa2,
                          ptca0, ptca1, ptca2, ptcb0, ptcb1, ptcb2,
                          pa0, pa1, pa2):
    """
    Description:

        OOI Level 1 Pressure (Depth) data product, which is calculated using
        data from the Sea-Bird Electronics conductivity, temperature and depth
        (CTD) family of instruments. This document is intended to be used by
        OOI programmers to construct appropriate processes to create the L1
        water temperature product.
        
        This data product is derived from SBE 16Plus instruments outfitted with
        a strain gauge pressure sensor. This is the default for most of the
        CTDBP instruments

    Implemented by:

        2013-04-12: Chris Wingard. Initial Code.
        2013-05-10: Christopher Wingard. Minor edits to comments.

    Usage:

        p = ctd_sbe16plus_preswat(p0, therm0, ptempa0, ptempa1, ptempa2,
                          ptca0, ptca1, ptca2, ptcb0, ptcb1, ptcb2,
                          pa0, pa1, pa2)

            where

        p = sea water pressure (PRESWAT_L1) [dbar]
        p0 = raw pressure (PRESWAT_L0) [counts]
        t0 = raw temperature (TEMPWAT_L0) [counts]
        ptempa0 = strain gauge pressure calibration coefficients
        ptempa1 = strain gauge pressure calibration coefficients
        ptempa2 = strain gauge pressure calibration coefficients
        ptca0 = strain gauge pressure calibration coefficients
        ptca1 = strain gauge pressure calibration coefficients
        ptca2 = strain gauge pressure calibration coefficients
        ptcb0 = strain gauge pressure calibration coefficients
        ptcb1 = strain gauge pressure calibration coefficients
        ptcb2 = strain gauge pressure calibration coefficients
        pa0 = strain gauge pressure calibration coefficients
        pa1 = strain gauge pressure calibration coefficients
        pa2 = strain gauge pressure calibration coefficients
        
    References:
    
        OOI (2012). Data Product Specification for Pressure (Depth). Document
            Control Number 1341-00020. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00020_Data_Product_SPEC_PRESWAT_OOI.pdf)
    """
    # compute calibration parameters
    tv = t0 / 13107.0
    t = ptempa0 + ptempa1 * tv + ptempa2 * tv**2
    x = p0 - ptca0 - ptca1 * t - ptca2 * t**2
    n = x * ptcb0 / (ptcb0 + ptcb1 * t + ptcb2 * t**2)
    
    # compute pressure in psi, rescale and compute in dbar and return
    p_psi = pa0 + pa1 * n + pa2 * n**2
    p_dbar = (p_psi * 0.689475729) - 10.1325
    return p_dbar


def ctd_sbe16digi_preswat(p0, t0, C1, C2, C3, D1, D2, T1, T2, T3, T4, T5):
    """
    Description:

        OOI Level 1 Pressure (Depth) data product, which is calculated using
        data from the Sea-Bird Electronics conductivity, temperature and depth
        (CTD) family of instruments. This document is intended to be used by
        OOI programmers to construct appropriate processes to create the L1
        water temperature product.
        
        This data product is derived from SBE 16Plus instruments outfitted with
        a digiquartz pressure sensor. This applies to the CTDBP-N,O instruments
        only.

    Implemented by:

        2013-05-10: Christopher Wingard. Initial Code.
        2013-05-10: Christopher Wingard. Minor edits to comments.

    Usage:

        p = ctd_sbe16digi_preswat(p0,t0,C1,C2,C3,D1,D2,T1,T2,T3,T4,T5)

            where

        p = sea water pressure (PRESWAT_L1) [dbar]
        p0 = raw pressure frequency (PRESWAT_L0) [counts]
        t0 = raw temperature (TEMPWAT_L0) [counts]
        C1 = digiquartz pressure calibration coefficients
        C2 = digiquartz pressure calibration coefficients
        C3 = digiquartz pressure calibration coefficients
        D1 = digiquartz pressure calibration coefficients
        D2 = digiquartz pressure calibration coefficients
        T1 = digiquartz pressure calibration coefficients
        T2 = digiquartz pressure calibration coefficients
        T3 = digiquartz pressure calibration coefficients
        T4 = digiquartz pressure calibration coefficients
        T5 = digiquartz pressure calibration coefficients
        
    References:
    
        OOI (2012). Data Product Specification for Pressure (Depth). Document
            Control Number 1341-00020. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00020_Data_Product_SPEC_PRESWAT_OOI.pdf)
    """
    # Convert raw temperature input to voltage
    tv = t0 / 13107.0

    # Calculate U (thermistor temp):
    U = (23.7 * (tv + 9.7917)) - 273.15

    # Calculate calibration parameters
    C = C1 + C2 * U + C3 * U**2	
    D = D1 + D2 * U
    T0 = T1 + T2 * U + T3 * U**2 + T4 * U**3 + T5 * U**4

    # Calculate T (pressure period, in microseconds):
    T = (1.0 / p0) * 1.0e6

    # compute pressure in psi, rescale and compute in dbar and return
    p_psi = C * (1.0 - T0**2 / T**2) * (1.0 - D * (1.0 - T0**2 / T**2))
    p_dbar = (p_psi * 0.689475729) - 10.1325
    return p_dbar
    

def ctd_sbe16plus_condwat(c0, t1, p1, g, h, i, j, cpcor, ctcor):
    """
    Description:

        OOI Level 1 Conductivity core data product, which is calculated using
        data from the Sea-Bird Electronics conductivity, temperature and depth
        (CTD) family of instruments. This document is intended to be used by
        OOI programmers to construct appropriate processes to create the L1
        water temperature product. 

    Implemented by:

        2013-04-12: Christopher Wingard. Initial Code
        2013-05-10: Christopher Wingard. Minor edits to comments.

    Usage:

        c = ctd_sbe16plus_condwat(c0, t1, p1, g, h, i, j, cpcor, ctcor)

            where

        c = sea water conductivity (CONDWAT_L1) [S m-1]
        t1 = sea water temperature (TEMPWAT_L1) [deg_C]
        p1 = sea water pressure (PRESWAT_L1) [deg_C]
        g = conductivity calibration coefficients
        h = conductivity calibration coefficients
        i = conductivity calibration coefficients
        j = conductivity calibration coefficients
        cpcor = conductivity calibration coefficients
        ctcor = conductivity calibration coefficients
        
    References:
    
        OOI (2012). Data Product Specification for Conductivity. Document
            Control Number 1341-00030. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00030_Data_Product_SPEC_CONDWAT_OOI.pdf)
    """
    # convert raw conductivty measurement to frequency
    f = (c0 / 256.0) / 1000.0
    
    # calculate conductivity [S m-1]
    c = (g + h * f**2 + i * f**3 + j * f**4) / (1 + ctcor * t1 + cpcor * p1)
    return c


def ctd_pracsal(c, t, p):
    """
    Description:

        OOI Level 2 Practical Salinity core data product, which is calculated
        using the Thermodynamic Equations of Seawater - 2010 (TEOS-10) Version
        3.0, with data from the conductivity, temperature and depth (CTD)
        family of instruments. This calculation is defined in the Data Product
        Specification for Salinty - DCN 1341-00040.

    Implemented by:

        2013-03-13: Christopher Wingard. Initial code.
        2013-05-10: Christopher Wingard. Minor edits to comments.

    Usage:

        SP = ctd_pracsal(c, t, p)

            where

        SP = Practical Salinity (seawater salinity, PSS-78) [unitless]
        c = conductivity (seawater conductivity) [S m-1], (see
            1341-00010_Data_Product_Spec_CONDWAT) 
        t = temperature (seawater temperature) [deg_C], (see
            1341-00030_Data_Product_Spec_TEMPWAT)
        p = pressure (sea pressure) [dbar], (see
            1341-00020_Data_Product_Spec_PRESWAT)

    References:
    
        OOI (2012). Data Product Specification for Salinty. Document Control
            Number 1341-00040. https://alfresco.oceanobservatories.org/ (See: 
            Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00040_Data_Product_SPEC_PRACSAL_OOI.pdf)
    """
    # Convert L1 Conductivity from S/m to mS/cm
    C10 = c * 10.0
    
    # Calculate the Practical Salinity (PSS-78) [unitless]
    SP = gsw.sp_from_c(C10, t, p)
    return SP


def ctd_density(SP, t, p, lat, lon):
    """
    Description:
    
        OOI Level 2 Density core data product, which is calculated using the
        Thermodynamic Equations of Seawater - 2010 (TEOS-10) Version 3.0, with
        data from the conductivity, temperature and depth (CTD) family of
        instruments. This calculation is defined in the Data Product
        Specification for Density - DCN 1341-00050.
        
    Implemented by:
    
        2013-03-11: Christopher Mueller. Initial code.
        2013-03-13: Christopher Wingard. Added commenting and moved to
            ctd_functions
        2013-05-10: Christopher Wingard. Minor edits to comments.

    Usage:
    
        rho = ctd_density(SP, t, p, lat, lon)
        
            where
    
        rho = Density (seawater density) [kg m^-3]
        SP = Practical Salinity (PSS-78) [unitless], (see
            1341-00040_Data_Product_Spec_PRACSAL)
        t = temperature (seawater temperature) [deg_C], (see
            1341-00010_Data_Product_Spec_TEMPWAT)
        p = pressure (sea pressure) [dbar], (see
            1341-00020_Data_Product_Spec_PRESWAT)
        lat = latitude where input data was collected [decimal degree]
        lon = longitude where input data was collected [decimal degree]
            
    References:
    
        OOI (2012). Data Product Specification for Density. Document Control
            Number 1341-00050. https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00050_Data_Product_SPEC_DENSITY_OOI.pdf)
    """
    # Calculate the Absolute Salinity (SA) from the Practical Salinity (SP)
    # [g kg^-1]
    rho = gsw.ctd_density(SP,t,p,lat,lon)
    # Calculate the Conservative Temperature (CT) [deg_C]
    return rho
