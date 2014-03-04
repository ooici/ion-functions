#!/usr/bin/env python
"""
@package ion_functions.data.sfl_functions
@file ion_functions/data/sfl_functions.py
@author Christopher Wingard
@brief Module containing Seafloor Properties related data-calculations.
"""
from ion_functions.utils import fill_value
def sfl_trhph_vfltemp(V_s, V_c, a, b, c, d, e):
    """
    Description:
    
        OOI Level 1 Vent Fluid Temperature from TRHPH (TRHPHTE) data product,
        which is calculated using data from the Temperature Resistivity Probe
        (TRHPH) instrument.
    
    Implemented by:
    
        2013-05-01: Christopher Wingard. Initial Code
    
    Usage:
    
        T = sfl_trhph_vftemp(V_s, V_c, a, b, c, d, e)
    
            where
    
        T = Vent fluid temperature from TRHPH [deg_C]
        V_c = Thermocouple voltage [volts]
        V_s = Thermistor voltage [volts]
        a = coefficient from 4th degree polynomial fit of laboratory
            calibration correction curve.
        b = coefficient from 4th degree polynomial fit of laboratory
            calibration correction curve.
        c = coefficient from 4th degree polynomial fit of laboratory
            calibration correction curve.
        d = coefficient from 4th degree polynomial fit of laboratory
            calibration correction curve.
        e = coefficient from 4th degree polynomial fit of laboratory
            calibration correction curve.
    
    References:
    
        OOI (2012). Data Product Specification for Vent Fluid Temperature from
            TRHPH. Document Control Number 1341-00150.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00150_Data_Product_SPEC_TRHPHTE_OOI.pdf)
    """
    # raw thermistor temperature
    T_s = 27.50133 - 17.2658 * V_s + 15.83424 / V_s
    
    # raw thermocouple temperature
    T_c = 244970. * V_c / 1000.
    
    # uncorrected total temperature
    T_u = T_s + T_c
    
    # correction based on laboratory calibration
    T_lc = a * T_u**4 + b * T_u**3 + c * T_u**2 + d * T_u + e 
    
    # final, corrected temperature at sensor tip    
    T = T_u + T_lc
    
    return T


def sfl_trhph_chlorconc(V_R1, V_R2, V_R3, T):
    """
    Wrapper function to vectorize the vent fluid chloride calculation defined
    below in sfl_trhph_chloride.
    """
    import numpy as np
    
    sfunc = np.vectorize(sfl_trhph_chloride)
    Cl = sfunc(V_R1, V_R2, V_R3, T)
    
    return Cl


def sfl_trhph_chloride(V_R1, V_R2, V_R3, T):
    """
    Description:
    
        OOI Level 1 Vent Fluid Chloride Concentration from TRHPH (TRHPHCC) data
        product, which is calculated using data from the Temperature
        Resistivity Probe (TRHPH) instrument.
    
    Implemented by:
    
        2013-05-01: Christopher Wingard. Initial Code
    
    Usage:
    
        Cl = sfl_trhph_chloride(V_R1, V_R2, V_R3, T)
    
            where
    
        Cl = Vent fluid chloride concentration from TRHPH [mmol kg-1]
        V_R1 = Resistivity voltage 1 (TRHPHR1) [volts]
        V_R2 = Resistivity voltage 1 (TRHPHR2) [volts]
        V_R3 = Resistivity voltage 1 (TRHPHR3) [volts]
        T = Vent fluid temperature from TRHPH (TRHPHTE) [deg_C]
    
    References:
    
        OOI (2012). Data Product Specification for Vent Fluid Temperature from
            TRHPH. Document Control Number 1341-00150.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00150_Data_Product_SPEC_TRHPHTE_OOI.pdf)
    """
    import numpy as np
    from scipy.interpolate import RectBivariateSpline
    
    # load sfl_functions_surface.py This loads the 3-dimensional calibration
    # surface of temperature, conductivity and chloride reproduced as numpy
    # arrays from Larson_2007surface.mat
    from ion_functions.data.sfl_functions_surface import tdat, sdat, cdat
    
    # select the optimal L0 Resistivity voltage
    V_R = V_R1 * 5.              # Option 1, default (V_R1 * 5)
    
    vflag = np.where(V_R2 < 0.75)         # Option 2 
    V_R[vflag] = V_R3[vflag] / 5.
    
    vflag = np.where((V_R2 >= 0.75) & (V_R2 < 3.90))    # Option 3
    V_R[vflag] = V_R2[vflag]
    
    # convert resistivity to conductivity
    C = 1. / V_R
    
    # extract curves of constant temperature out of the data surfaces
    Cl = np.zeros(len(C))
    Scurve = np.linspace(np.min(sdat), np.max(sdat), 100,
                         endpoint='True')
    f = RectBivariateSpline(tdat, sdat, cdat.T, kx=1, ky=1, s=0)
    for i in range(len(Cl)):
        Tcurve = np.zeros(len(Scurve)) + T[i]
        Ccurve = f(Tcurve, Scurve)
        if (np.all(np.isfinite(Ccurve)) and
            (T[i] >= min(tdat) and T[i] <= max(tdat))):
            #now interpolate onto the Scurve/Ccurve
            S = np.interp(C[i], Ccurve[0,:], Scurve, left=np.nan, right=np.nan)
            Cl[i] = np.round(S * 1000.)
        else:
            Cl[i] = np.nan
        
    return Cl

import numexpr as ne

def sfl_SFLPRES_L1(p_psia):
    """
    Description:

        The OOI Level 1 Seafloor Pressure core data products, SFLPRES and sub-parameters
        SFLPRES-RTIME, SFLPRES-TIDE, and SFLPRES-WAVE, are created from the Sea-Bird
        Electronics SBE 26plus member of the Pressure SF (PRESF) family of instruments by either a)
        polling, in real-time, for L0 ASCII text format data output and converting from psia to decibar units
        or b) converting, after instrument recovery, L0 hexadecimal pressure data into decimal format and
        the resulting tide and wave pressure data in psia to decibar units.

    Implemented by:

        2014-01-31: Craig Risien. Initial Code
        
    Usage:

        SFLPRES_L1 = sfl_SFLPRES_L1(p_psia):

        Scaling: To convert from psia to dbar, use the Sea-Bird-specified conversion:
        p_dbar = 0.689475728 *(p_psia)

            where

        p_dbar = pressure (SFLPRES_L1) [dbar]
        p_psia = pressure (SFLPRES_L0) [psi].

    References:
    
        OOI (2013). Data Product Specification for Seafloor Pressure from Sea-Bird SBE 26PLUS.
        Document Control Number 1341-00230. https://alfresco.oceanobservatories.org/
        (See: Company Home >> OOI >> Controlled >> 1000 System Level >> 1341-00230_Data_Product_SPEC_SFLPRES_OOI.pdf)
    """

    SFLPRES_L1 = ne.evaluate('p_psia * 0.689475728')
    
    return SFLPRES_L1
