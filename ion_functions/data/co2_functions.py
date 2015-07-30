#!/usr/bin/env python
"""
@package ion_functions.data.pco2_functions
@file ion_functions/data/pco2_functions.py
@author Christopher Wingard
@brief Module containing CO2 instrument family related functions
"""

import numpy as np
import numexpr as ne
import scipy as sp
from ion_functions.utils import fill_value


# wrapper functions to extract parameters from SAMI-II CO2 instruments (PCO2W)
# and process these extracted parameters to calculate pCO2
def pco2_abs434_ratio(light):
    """
    Description:

        Extract the absorbance ratio at 434 nm from the pCO2 instrument light
        measurements. This will extract the CO2ABS1_L0 data product from the
        instrument light measurement arrays.

    Implemented by:

        2013-04-20: Christopher Wingard. Initial code.
        2014-02-19: Christopher Wingard. Updated comments.

    Usage:

        a434ratio = pco2_abs434_ratio(light)

            where

        a434ratio = optical absorbance Ratio at 434 nm (CO2ABS1_L0) [unitless]
        light = array of light measurements

    References:

        OOI (2012). Data Product Specification for Partial Pressure of CO2 in
            Seawater. Document Control Number 1341-00510.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00490_Data_Product_SPEC_PCO2WAT_OOI.pdf)
    """
    light = np.atleast_2d(light)
    a434ratio = light[:, 6]
    return a434ratio


def pco2_abs620_ratio(light):
    """
    Description:

        Extract the absorbance ratio at 620 nm from the pCO2 instrument light
        measurements. This will extract the CO2ABS2_L0 data product from the
        instrument light measurement arrays.

    Implemented by:

        2013-04-20: Christopher Wingard. Initial code.
        2014-02-19: Christopher Wingard. Updated comments.

    Usage:

        a620ratio = pco2_abs620_ratio(light)

            where

        a620ratio = optical absorbance Ratio at 620 nm (CO2ABS2_L0) [unitless]
        light = array of light measurements

    References:

        OOI (2012). Data Product Specification for Partial Pressure of CO2 in
            Seawater. Document Control Number 1341-00510.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00490_Data_Product_SPEC_PCO2WAT_OOI.pdf)
    """
    light = np.atleast_2d(light)
    a620ratio = light[:, 7]
    return a620ratio


def pco2_blank(raw_blank):
    """
    Description:

        Calculates the absorbance blank at 434 or 620 nm from the SAMI2-pCO2
        instrument.

    Implemented by:

        2013-04-20: Christopher Wingard. Initial code.
        2014-02-19: Christopher Wingard. Updated comments.
        2014-02-28: Christopher Wingard. Updated to except raw blank values
                    from a sparse array.

    Usage:

        blank = pco2_blank(raw_blank)

            where

        blank = optical absorbance blank at 434 or 620 nm [unitless]
        raw_blank = raw optical absorbance blank at 434 or 620 nm [counts]

    References:

        OOI (2012). Data Product Specification for Partial Pressure of CO2 in
            Seawater. Document Control Number 1341-00510.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00490_Data_Product_SPEC_PCO2WAT_OOI.pdf)
    """
    #blank = -1. * sp.log10(raw_blank / 16384.)
    blank = -1. * sp.log10(raw_blank)

    return blank


def pco2_thermistor(traw):
    """
    Description:

        Convert the thermistor data from counts to degrees Centigrade from the
        pCO2 instrument.

    Implemented by:

        2013-04-20: Christopher Wingard. Initial code.
        2014-02-19: Christopher Wingard. Updated comments.

    Usage:

        therm = pco2_thermistor(traw)

            where

        therm = converted thermistor temperature [degC]
        traw = raw thermistor temperature (CO2THRM_L0) [counts]

    References:

        OOI (2012). Data Product Specification for Partial Pressure of CO2 in
            Seawater. Document Control Number 1341-00510.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00490_Data_Product_SPEC_PCO2WAT_OOI.pdf)
    """
    # convert raw thermistor readings from counts to degrees Centigrade
    Rt = ne.evaluate('log((traw / (4096. - traw)) * 17400.)')
    InvT = ne.evaluate('0.0010183 + 0.000241 * Rt + 0.00000015 * Rt**3')
    therm = ne.evaluate('(1 / InvT) - 273.15')
    return therm


def pco2_pco2wat(mtype, light, therm, ea434, eb434, ea620, eb620,
                 calt, cala, calb, calc, a434blank, a620blank):
    """
    Description:

        Function to calculate the L1 PCO2WAT core data from the pCO2 instrument
        if the measurement type is 4 (pCO2 measurement), otherwise it is a
        blank and return a fill value.

    Implemented by:

        2013-04-20: Christopher Wingard. Initial code.
        2014-02-19: Christopher Wingard. Updated comments.
        2014-03-19: Christopher Wingard. Optimized using feedback provided by
                    Chris Fortin.

    Usage:

        pco2 = pco2_pco2wat(mtype, light, therm, ea434, eb434, ea620, eb620,
                            calt, cala, calb, calc, a434blank, a620blank)

            where

        pco2 = measured pco2 in seawater (PCO2WAT_L1) [uatm]
        mtype = measurement type, where 4 == actual measurement and 5 == a
            blank measurement [unitless]
        light = array of light measurements
        therm = PCO2W thermistor temperature (CO2THRM_L0) [counts]
        ea434 = Reagent specific calibration coefficient
        eb434 = Reagent specific calibration coefficient
        ea620 = Reagent specific calibration coefficient
        eb620 = Reagent specific calibration coefficient
        calt = Instrument specific calibration coefficient for temperature
        cala = Instrument specific calibration coefficient for the pCO2 measurements
        calb = Instrument specific calibration coefficient for the pCO2 measurements
        calc = Instrument specific calibration coefficient for the pCO2 measurements
        a434blank = Blank measurements at 434 nm (CO2ABS1_L0) [counts]
        a620blank = Blank measurements to 620 nm (CO2ABS2_L0) [counts]

    References:

        OOI (2012). Data Product Specification for Partial Pressure of CO2 in
            Seawater. Document Control Number 1341-00510.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00490_Data_Product_SPEC_PCO2WAT_OOI.pdf)
    """
    # reset inputs to arrays
    ### measurements
    mtype = np.atleast_1d(mtype)
    light = np.atleast_2d(light)
    therm = np.atleast_1d(therm)
    ### calibration coefficients
    ea434 = np.atleast_1d(ea434)
    eb434 = np.atleast_1d(eb434)
    ea620 = np.atleast_1d(ea620)
    eb620 = np.atleast_1d(eb620)
    calt = np.atleast_1d(calt)
    cala = np.atleast_1d(cala)
    calb = np.atleast_1d(calb)
    calc = np.atleast_1d(calc)
    ### blank measurements
    a434blank = np.atleast_1d(a434blank)
    a620blank = np.atleast_1d(a620blank)

    # calculate the pco2 value
    pco2 = pco2_calc_pco2(light, therm, ea434, eb434, ea620, eb620,
                          calt, cala, calb, calc, a434blank, a620blank)

    # reset dark measurements to the fill value
    m = np.where(mtype == 5)[0]
    pco2[m] = fill_value

    return pco2


# L1a PCO2WAT calculation
def pco2_calc_pco2(light, therm, ea434, eb434, ea620, eb620,
                   calt, cala, calb, calc, a434blank, a620blank):
    """
    Description:

        OOI Level 1 Partial Pressure of CO2 (pCO2) in seawater core data
        product, which is calculated from the Sunburst SAMI-II CO2 instrument
        (PCO2W).

    Implemented by:

        20??-??-??: J. Newton (Sunburst Sensors, LLC). Original Matlab code.
        2013-04-20: Christopher Wingard. Initial python code.
        2014-02-19: Christopher Wingard. Updated comments.
        2014-03-19: Christopher Wingard. Optimized.

    Usage:

        pco2 = pco2_pco2wat(light, therm, ea434, eb434, ea620, eb620,
                            calt, cala, calb, calc, a434blank, a620blank)

            where

        pco2 = measured pco2 in seawater (PCO2WAT_L1) [uatm]
        light = array of light measurements
        therm = PCO2W thermistor temperature (CO2THRM_L0) [counts]
        ea434 = Reagent specific calibration coefficient
        eb434 = Reagent specific calibration coefficient
        ea620 = Reagent specific calibration coefficient
        eb620 = Reagent specific calibration coefficient
        calt = Instrument specific calibration coefficient for temperature
        cala = Instrument specific calibration coefficient for the pCO2 measurements
        calb = Instrument specific calibration coefficient for the pCO2 measurements
        calc = Instrument specific calibration coefficient for the pCO2 measurements
        a434blank = Blank measurements at 434 nm (CO2ABS1_L0) [counts]
        a620blank = Blank measurements to 620 nm (CO2ABS2_L0) [counts]

    References:

        OOI (2012). Data Product Specification for Partial Pressure of CO2 in
            Seawater. Document Control Number 1341-00510.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00490_Data_Product_SPEC_PCO2WAT_OOI.pdf)
    """
    # set constants
    ea434 = ea434 - 29.3 * calt
    eb620 = eb620 - 70.6 * calt
    e1 = ea620 / ea434
    e2 = eb620 / ea434
    e3 = eb434 / ea434

    # Extract variables from light array
    light = light.astype(np.float)
    #DRef1 = light[0]  # Dark Reference LED
    #DSig1 = light[1]  # Dark Signal LED
    #R434 = light[2]   # 434nm Reference LED intensity
    #S434 = light[3]   # 434nm Signal Signal LED intensity
    #R620 = light[4]   # 620nm Reference LED intensity
    #S620 = light[5]   # 434nm Signal Signal LED intensity
    Ratio434 = light[:, 6]     # 434nm Ratio
    Ratio620 = light[:, 7]     # 620nm Ratio

    # calculate absorbance ratio, correcting for blanks
    A434 = -1. * sp.log10(Ratio434 / a434blank)  # 434 absorbance
    A620 = -1. * sp.log10(Ratio620 / a620blank)  # 620 absorbance
    Ratio = A620 / A434      # Absorbance ratio

    # calculate pCO2
    V1 = Ratio - e1
    V2 = e2 - e3 * Ratio
    RCO21 = -1. * sp.log10(V1 / V2)
    RCO22 = (therm - calt) * 0.007 + RCO21
    Tcoeff = 0.0075778 - 0.0012389 * RCO22 - 0.00048757 * RCO22**2
    Tcor_RCO2 = RCO21 + Tcoeff * (therm - calt)
    pco2 = 10.**((-1. * calb + (calb**2 - (4. * cala * (calc - Tcor_RCO2)))**0.5) / (2. * cala))

    return np.real(pco2)


def pco2_ppressure(xco2, p, std=1013.25):
    """
    Description:

        OOI Level 1 Partial Pressure of CO2 in Air (PCO2ATM) or Surface
        Seawater (PCO2SSW) core date product is computed by using an
        equation that incorporates the Gas Stream Pressure (PRESAIR) and the
        CO2 Mole Fraction in Air or Surface Seawater (XCO2ATM or XCO2SSW,
        respectively). It is computed using data from the pCO2 air-sea (PCO2A)
        family of instruments.

    Implemented by:

        2014-10-27: Christopher Wingard. Initial python code.

    Usage:

        ppres = pco2_ppressure(xco2, p, std)

            where

        ppres = partial pressure of CO2 in air or surface seawater [uatm]
                (PCO2ATM_L1 or PCO2SSW_L1)
        xco2 = CO2 mole fraction in air or surface seawater [ppm] (XCO2ATM_LO
               or XCO2SSW_L0)
        p = gas stream pressure [mbar] (PRESAIR_L0)
        std = standard atmospheric pressure set to default of 1013.25 [mbar/atm]

    References:

        OOI (2012). Data Product Specification for Partial Pressure of CO2 in
            Air and Surface Seawater. Document Control Number 1341-00260.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00260_Data_Product_SPEC_PCO2ATM_PCO2SSW_OOI.pdf)
    """
    ppres = ne.evaluate('xco2 * p / std')
    return ppres


def pco2_co2flux(pco2w, pco2a, u10, t, s):
    """
    Description:

        OOI Level 2 core date product CO2FLUX is an estimate of the CO2 flux
        from the ocean to the atmosphere. It is computed using data from the
        pCO2 air-sea (PCO2A) and bulk meteorology (METBK) families of
        instruments.

    Implemented by:

        2012-03-28: Mathias Lankhorst. Original Matlab code.
        2013-04-20: Christopher Wingard. Initial python code.

    Usage:

        flux = pco2_co2flux(pco2w, pco2a, u10, t, s)

            where

        flux = estimated flux of CO2 from the ocean to atmosphere [mol m-2 s-1]
               (CO2FLUX_L2)
        pco2w = partial pressure of CO2 in sea water [uatm] (PCO2SSW_L1)
        pco2a = partial pressure of CO2 in air [uatm] (PCO2ATM_L1)
        u10 = normalized wind speed at 10 m height from METBK [m s-1] (WIND10M_L2)
        t = sea surface temperature from METBK [deg_C] (TEMPSRF_L1)
        s = sea surface salinity from METBK [psu] (SALSURF_L2)

    References:

        OOI (2012). Data Product Specification for Flux of CO2 into the
            Atmosphere. Document Control Number 1341-00270.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00270_Data_Product_SPEC_CO2FLUX_OOI.pdf)
    """
    # convert micro-atm to atm
    pco2a = pco2a / 1.0e6
    pco2w = pco2w / 1.0e6

    # Compute Schmidt number (after Wanninkhof, 1992, Table A1)
    Sc = 2073.1 - (125.62 * t) + (3.6276 * t**2) - (0.043219 * t**3)

    # Compute gas transfer velocity (after Sweeney et al. 2007, Fig. 3 and Table 1)
    k = 0.27 * u10**2 * sp.sqrt(660.0 / Sc)

    # convert cm h-1 to m s-1
    k = k / (100.0 * 3600.0)

    # Compute the absolute temperature
    T = t + 273.15

    # Compute solubility (after Weiss 1974, Eqn. 12 and Table I).
    # Note that there are two versions, one for units per volume and
    # one per mass. Here, the volume version is used.
    # mol atm-1 m-3
    T100 = T / 100
    K0 = 1000 * np.exp(-58.0931 + (90.5069 * (100/T)) + (22.2940 * np.log(T100)) +
                       s * (0.027766 - (0.025888 * T100) + (0.0050578 * T100**2)))

    # mol atm-1 kg-1
    #K0 = np.exp(-60.2409 + (93.4517 * (100/T)) + (23.3585 * np.log(T100)) +
    #            s * (0.023517 - (0.023656 * T100) + (0.0047036 * T100**2)))

    # Compute flux (after Wanninkhof, 1992, eqn. A2)
    flux = k * K0 * (pco2w - pco2a)
    return flux
