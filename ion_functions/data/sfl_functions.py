#!/usr/bin/env python
"""
@package ion_functions.data.sfl_functions
@file ion_functions/data/sfl_functions.py
@author Christopher Wingard
@brief Module containing Seafloor Properties related data-calculations.
"""

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
    
    vflag = V_R2 < 0.75         # Option 2 
    V_R[vflag] = V_R3[vflag] / 5.
    
    vflag = (V_R2 >= 0.75) & (V_R2 < 3.90)    # Option 3
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
        
    # reset NaN values generated in interpolation functions above to system
    # default of -99999999
    Cl[np.isnan(Cl)] = -99999999.
    
    return Cl
