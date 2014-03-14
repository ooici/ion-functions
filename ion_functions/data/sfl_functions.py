#!/usr/bin/env python
"""
@package ion_functions.data.sfl_functions
@file ion_functions/data/sfl_functions.py
@author Christopher Wingard
@brief Module containing Seafloor Properties related data-calculations.
"""

import numpy as np
import numexpr as ne

from ion_functions.utils import fill_value


def sfl_trhph_vfltemp(V_s, V_c, a, b, c, d, e):
    """
    Description:

        OOI Level 1 Vent Fluid Temperature from TRHPH (TRHPHTE) data product,
        which is calculated using data from the Temperature Resistivity Probe
        (TRHPH) instrument.

    Implemented by:

        2013-05-01: Christopher Wingard. Initial Code
        2014-02-27: Russell Desiderio. Added documentation.
                    Implemented old-style polynomial calculation.

    Usage:

        T = sfl_trhph_vfltemp(V_s, V_c, a, b, c, d, e)

            where

        T = Vent fluid temperature from TRHPH (TRHPHTE_L1) [deg_C]
        V_s = Thermistor voltage (TRHPHVS_L0) [volts]
        V_c = Thermocouple voltage (TRHPHVC_L0) [volts]
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
            1341-00150_Data_Product_Spec_TRHPHTE_OOI.pdf)
    """
    # raw thermistor temperature
    T_s = 27.50133 - 17.2658 * V_s + 15.83424 / V_s

    # raw thermocouple temperature (PD431)
    T_c = 244970. * V_c / 1000.

    # uncorrected total temperature
    T_u = T_s + T_c

    # correction based on laboratory calibration
    #Tu_sq = T_u * T_u
    #T_lc = a * Tu_sq*Tu_sq + b * T_u*Tu_sq + c * Tu_sq + d * T_u + e

    # polynomial calculation optimization - this appears to be a little faster
    T_lc = e + T_u * (d + T_u * (c + T_u * (b + T_u * a)))

    # final, corrected temperature at sensor tip
    T = T_u + T_lc

    return T


def sfl_trhph_vfl_thermistor_temp(V_s):
    """
    Description:

        Calculates T_S, which is an intermediate data product (not a core data product)
        requested by the authors of the TRHPHTE DPS. It is the instrument's thermistor
        temperature, useful as an important instrument diagnostic. It is the same variable
        as T_s in the function sfl_trhph_vfltemp.

    Implemented by:

        2014-02-28: Russell Desiderio. Initial Code

    Usage:

        T_s = sfl_trhph_vfl_thermistor_temp(V_s)

            where

        T_s = Thermistor reference temperature [deg_C]
        V_s = Thermistor voltage (TRHPHVS_L0) [volts]

    References:

        OOI (2012). Data Product Specification for Vent Fluid Temperature from
            TRHPH. Document Control Number 1341-00150.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00150_Data_Product_Spec_TRHPHTE_OOI.pdf)
    """
    # thermistor temperature
    T_s = 27.50133 - 17.2658 * V_s + 15.83424 / V_s

    return T_s


def sfl_trhph_vflorp(V, offset, gain):
    """
    Description:

        OOI Level 1 Vent Fluid Oxidation-Reduction Potential (ORP) from TRHPH
        (TRHPHEH_L1) data product, which is calculated using data from the Resistivity-
        Temperature Probe (TRHPH) instrument.

    Implemented by:

        2014-02-28: Russell Desiderio. Initial Code

    Usage:

        ORP = sfl_trhph_vflorp(V)

            where

        ORP = Oxidation-Reduction Potential (TRHPHEH_L1) [mV]
        V = ORP sensor voltage (TRHPHVO_L0) [volts]
        offset = calibration coefficient; electronic offset [mV]
        gain = calibration coefficient; gain multiplier [unitless]

    References:

        OOI (2012). Data Product Specification for Vent Fluid Oxidation-
            Reduction Potential (ORP). Document Control Number 1341-00170.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00170_Data_Product_Spec_TRHPHEH_OOI.pdf)
    """
    # convert sensor voltage V from volts to mV;
    # subtract offset; undo gain multiplier.
    ORP = np.round((V * 1000.0 - offset)/gain)

    return ORP


def sfl_trhph_chloride(V_R1, V_R2, V_R3, T):
    """
    Description:

        OOI Level 2 Vent Fluid Chloride Concentration TRHPHCC data
        product, which is calculated using data from the Temperature
        Resistivity Probe (TRHPH) instrument.

    Implemented by:

        2013-05-01: Christopher Wingard. Initial Code
        2014-02-28: Russell Desiderio. Modified code to better handle nans and fill_values.
                                       Added more documentation to algorithm.
        2014-03-10: Russell Desiderio. Removed unnecessary np.vectorized wrapper function.
                                       Improved speed by removing a temperature conditional
                                           statement from inside the for loop and incorporating
                                           it into the range of the for loop.

    Usage:

        Cl = sfl_trhph_chloride(V_R1, V_R2, V_R3, T)

            where

        Cl = Vent fluid chloride concentration from TRHPH (TRHPHCC_L2) [mmol kg-1]
        V_R1 = Resistivity voltage 1 (TRHPHR1_L0) [volts]
        V_R2 = Resistivity voltage 1 (TRHPHR2_L0) [volts]
        V_R3 = Resistivity voltage 1 (TRHPHR3_L0) [volts]
        T = Vent fluid temperature from TRHPH (TRHPHTE_L1) [deg_C]

    References:

        OOI (2012). Data Product Specification for Vent Fluid Chloride
        Concentration. Document Control Number 1341-00160.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00160_Data_Product_Spec_TRHPHCC_OOI.pdf)
    """
    from scipy.interpolate import RectBivariateSpline

    # load sfl_functions_surface.py This loads the 3-dimensional calibration
    # surface of temperature, chloride, and conductivity reproduced as numpy
    # arrays from Larson_2007surface.mat
    from ion_functions.data.sfl_functions_surface import tdat, sdat, cdat

    # select the optimal L0 Resistivity voltage.
    V_R = V_R3 / 5.0

    vflag = np.where((V_R2 >= 0.75) & (V_R2 < 3.90))
    V_R[vflag] = V_R2[vflag]

    vflag = np.where(V_R2 >= 3.90)
    V_R[vflag] = V_R1[vflag] * 5.0

    # convert resistivity to conductivity
    C = 1. / V_R

    # initialize product array Cl [mmol/kg] values to nans
    Cl = np.zeros(len(C)) + np.nan
    # set up chloride ['S' in units of mol/kg] range
    Scurve = np.linspace(np.min(sdat), np.max(sdat), 100,
                         endpoint='True')
    # create bivariate spline for interpolation
    f = RectBivariateSpline(tdat, sdat, cdat.T, kx=1, ky=1, s=0)

    # Note that when T is out-of-range, the interpolation np.interp does not
    # always give nan values for Cl as is required. Since Cl has been initialized
    # to nan values, iterate only over good T values, which also improves speed.
    for ii in np.where(np.logical_and(T >= min(tdat), T <= max(tdat)))[0]:
        # form constant T vector (Tcurve) for a given T value
        Tcurve = np.zeros(len(Scurve)) + T[ii]
        # find conductivity curve Ccurve as f(T=constant, chloride).
        Ccurve = f(Tcurve, Scurve)
        # now interpolate measured conductivity C into (Ccurve,Scurve) to get Cl.
        # this conditional statement is in the DPS and therefore retained.
        if (np.all(np.isfinite(Ccurve))):
            Cl[ii] = np.interp(C[ii], Ccurve[0, :], Scurve, left=np.nan, right=np.nan)

    # change units to mmol/kg; round to required # of sigfigs as specified in the DPS
    Cl = np.round(Cl * 1000.)
    Cl[np.isnan(Cl)] = fill_value

    return Cl


def sfl_sflpres_l1(p_psia):
    """
    Description:

        The OOI Level 1 Seafloor Pressure core data products, SFLPRES
        and sub-parameters SFLPRES-RTIME, SFLPRES-TIDE, and SFLPRES-WAVE,
        are created from the Sea-Bird Electronics SBE 26plus member of
        the Pressure SF (PRESF) family of instruments by either a)
        polling, in real-time, for L0 ASCII text format data output and
        converting from psia to decibar units or b) converting, after
        instrument recovery, L0 hexadecimal pressure data into decimal
        format and the resulting tide and wave pressure data in psia to
        decibar units.

    Implemented by:

        2014-01-31: Craig Risien. Initial Code

    Usage:

        sflpres_l1 = sfl_sflpres_l1(p_psia):

        Scaling: To convert from psia to dbar, use the Sea-Bird-specified
        conversion: p_dbar = 0.689475728 *(p_psia)

            where

        p_dbar = pressure (sflpres_L1) [dbar]
        p_psia = pressure (sflpres_L0) [psi].

    References:

        OOI (2013). Data Product Specification for Seafloor Pressure from
        Sea-Bird SBE 26PLUS. Document Control Number 1341-00230.
        https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
        >> Controlled >> 1000 System Level >>
        1341-00230_Data_Product_SPEC_SFLPRES_OOI.pdf)
    """

    sflpres_l1 = ne.evaluate('p_psia * 0.689475728')

    return sflpres_l1
