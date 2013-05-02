#!/usr/bin/env python
"""
@package ion_functions.data.sfl_functions
@file ion_functions/data/sfl_functions.py
@author Christopher Wingard
@brief Module containing Seafloor Properties related data-calculations.
"""

def sfl_trhph_vftemp(V_s, V_c, a, b, c, d, e):
    """
    Description:

        OOI Level 1 Vent Fluid Temperature from TRHPH (TRHPHTE) data product,
        which is calculated using data from the Temperature Resistivity Probe
        (TRHPH) instruments.

    Implemented by:

        2013-05-01: Christopher Wingard. Initial Code
        
    Usage:

        vftemp = sfl_trhph_vftemp(V_s, V_c, a, b, c, d, e)

            where

        vftemp = Vent fluid temperature from TRHPH [deg_C]
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
    vftemp = T_u + T_lc
    
    return vftemp


