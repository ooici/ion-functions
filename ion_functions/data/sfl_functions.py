#!/usr/bin/env python
"""
@package ion_functions.data.sfl_functions
@file ion_functions/data/sfl_functions.py
@author Christopher Wingard, Russell Desiderio
@brief Module containing Seafloor Properties related data-calculations.
"""

import numpy as np
import numexpr as ne
from scipy.interpolate import RectBivariateSpline

# used by def sfl_trhph_chloride
from ion_functions.data.sfl_functions_surface import tdat, sdat, cdat
from ion_functions.utils import fill_value


def sfl_thsph_temp_th(tc_rawdec_H, e2l_H, l2s_H, ts_rawdec_r, e2l_r, l2s_r, s2v_r):
    """
    Description:

        OOI Level 1 THSPH data product THSPHTE-TH (final temperature at position
        "H" near sample inlet), which is calculated using data from the Hydrothermal
        Vent Fluid In-situ Chemistry (THSPH) instrument, series A (one series for all
        instruments).

    Implemented by:

        2014-05-01: Russell Desiderio. Initial Code
        2014-06-30: Russell Desiderio. DPS modifications to cal equations implemented.

    Usage:

        T_H = sfl_thsph_temp_th(tc_rawdec_H, e2l_H, l2s_H, ts_rawdec_r, e2l_r, l2s_r, s2v_r)

            where

        T_H = final temperature "H" near sample inlet THSPHTE-TH_L1 [deg_C]
        #
        tc_rawdec_H = "H" thermocouple, decimal counts (THSPHTE-TCH_L0) [counts]
        e2l_H = array of calibration coefficients to convert the 'H' thermocouple
                engineering values to lab calibrated values.
        l2s_H = array of calibration coefficients to convert the 'H' thermocouple
                lab calibrated values to scientific values.
        ts_rawdec_r = reference thermistor, decimal counts (THSPHTE-REF_L0) [counts]
        e2l_r = array of calibration coefficients to convert the 'r' thermistor
                engineering values to lab calibrated values.
        l2s_r = array of calibration coefficients to convert the 'r' thermistor
                lab calibrated values to scientific values.
        s2v_r = array of calibration coefficients to convert the 'r' thermistor
                scientific values to thermocouple equivalent voltage [mV].

    References:

        OOI (2014). Data Product Specification for Vent Fluid Temperature from
            THSPH. Document Control Number 1341-00120.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00120_Data_Product_Specification_THSPHTE_OOI.pdf)
    """
    # calculate intermediate product V_tc_actual_H (= V_tc_labcal_H)
    V_tc_actual_H = sfl_thsph_temp_labcal_h(tc_rawdec_H, e2l_H)

    # calculate intermediate products T_ts_r, then V_ts_r (June 2014 DPS)
    T_ts_r = sfl_thsph_temp_ref(ts_rawdec_r, e2l_r, l2s_r)
    # the calculation of T_ts_r could result in non-finite (inf, nan) values, which
    # the sfl_thsph_temp_ref function above returns as fill values. so, replace
    # the fill values with np.nans so they can be tracked through to the final
    # data product variable.
    T_ts_r[np.equal(T_ts_r, fill_value)] = np.nan
    V_ts_r = eval_poly(T_ts_r, s2v_r)

    # Correct thermocouple temperature to account for offset from cold junction as
    # measured by the reference thermistor
    T_H = eval_poly((V_tc_actual_H + V_ts_r), l2s_H)

    # replace nans with fill values
    T_H[np.isnan(T_H)] = fill_value

    return T_H


def sfl_thsph_temp_tl(tc_rawdec_L, e2l_L, l2s_L, ts_rawdec_r, e2l_r, l2s_r, s2v_r):
    """
    Description:

        OOI Level 1 THSPH data product THSPHTE-TL (final temperature at position
        "L" near vent), which is calculated using data from the Hydrothermal Vent
        Fluid In-situ Chemistry (THSPH) instrument, series A (one series for all
        instruments).

    Implemented by:

        2014-05-01: Russell Desiderio. Initial Code
        2014-06-30: Russell Desiderio. DPS modifications to cal equations implemented.

    Usage:

        T_L = sfl_thsph_temp_tl(tc_rawdec_L, e2l_L, l2s_L, ts_rawdec_r, e2l_r, l2s_r, s2v_r)

            where

        T_L = final temperature "L" near vent THSPHTE-TL_L1 [deg_C]
        #
        tc_rawdec_L = "L" thermocouple, decimal counts (THSPHTE-TCL_L0) [counts]
        e2l_L = array of calibration coefficients to convert the 'L' thermocouple
                engineering values to lab calibrated values.
        l2s_L = array of calibration coefficients to convert the 'L' thermocouple
                lab calibrated values to scientific values.
        ts_rawdec_r = reference thermistor, decimal counts (THSPHTE-REF_L0) [counts]
        e2l_r = array of calibration coefficients to convert the 'r' thermistor
                engineering values to lab calibrated values.
        l2s_r = array of calibration coefficients to convert the 'r' thermistor
                lab calibrated values to scientific values.
        s2v_r = array of calibration coefficients to convert the 'r' thermistor
                scientific values to thermocouple equivalent voltage [mV].

    References:

        OOI (2014). Data Product Specification for Vent Fluid Temperature from
            THSPH. Document Control Number 1341-00120.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00120_Data_Product_Specification_THSPHTE_OOI.pdf)
    """
    # calculate intermediate product V_tc_actual_L (= V_tc_labcal_L)
    V_tc_actual_L = sfl_thsph_temp_labcal_l(tc_rawdec_L, e2l_L)

    # calculate intermediate products T_ts_r, then V_ts_r (June 2014 DPS)
    T_ts_r = sfl_thsph_temp_ref(ts_rawdec_r, e2l_r, l2s_r)
    # the calculation of T_ts_r could result in non-finite (inf, nan) values, which
    # the sfl_thsph_temp_ref function above returns as fill values. so, replace
    # the fill values with np.nans so they can be tracked through to the final
    # data product variable.
    T_ts_r[np.equal(T_ts_r, fill_value)] = np.nan
    V_ts_r = eval_poly(T_ts_r, s2v_r)

    # Correct thermocouple temperature to account for offset from cold junction as
    # measured by the reference thermistor
    T_L = eval_poly((V_tc_actual_L + V_ts_r), l2s_L)

    # replace nans with fill values
    T_L[np.isnan(T_L)] = fill_value

    return T_L


def sfl_thsph_temp_tch(tc_rawdec_H, e2l_H, l2s_H):
    """
    Description:

        OOI Level 1 THSPH data product THSPHTE-TCH (intermediate thermocouple
        temperature at position "H"), which is calculated using data from the
        Hydrothermal Vent Fluid In-situ Chemistry (THSPH) instrument, series A
        (one series for all instruments).

    Implemented by:

        2014-05-01: Russell Desiderio. Initial Code
        2014-06-30: Russell Desiderio. DPS modifications to cal equations implemented.

    Usage:

        T_tc_H = sfl_thsph_temp_tch(tc_rawdec_H, e2l_H, l2s_H)

            where

        T_tc_H = intermediate thermocouple temperature "H" THSPHTE-TCH_L1 [deg_C]
        tc_rawdec_H = "H" thermocouple, decimal counts (THSPHTE-TCH_L0) [counts]
        e2l_H = array of calibration coefficients to convert the 'H' thermocouple
                engineering values to lab calibrated values.
        l2s_H = array of calibration coefficients to convert the 'H' thermocouple
                lab calibrated values to scientific values.

    References:

        OOI (2014). Data Product Specification for Vent Fluid Temperature from
            THSPH. Document Control Number 1341-00120.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00120_Data_Product_Specification_THSPHTE_OOI.pdf)
    """
    # convert raw decimal output to lab calibrated values [mV]
    V_tc_actual_H = sfl_thsph_temp_labcal_h(tc_rawdec_H, e2l_H)

    # convert lab calibrated values to scientific values [degC]
    T_tc_H = eval_poly(V_tc_actual_H, l2s_H)

    return T_tc_H


def sfl_thsph_temp_tcl(tc_rawdec_L, e2l_L, l2s_L):
    """
    Description:

        OOI Level 1 THSPH data product THSPHTE-TCL (intermediate thermocouple
        temperature at position "L"), which is calculated using data from the
        Hydrothermal Vent Fluid In-situ Chemistry (THSPH) instrument, series A
        (one series for all instruments).

    Implemented by:

        2014-05-01: Russell Desiderio. Initial Code
        2014-06-30: Russell Desiderio. DPS modifications to cal equations implemented.

    Usage:

        T_tc_L = sfl_thsph_temp_tcl(tc_rawdec_L, e2l_L, l2s_L)

            where

        T_tc_L = intermediate thermocouple temperature "L" THSPHTE-TCL_L1 [deg_C]
        tc_rawdec_L = "L" thermocouple, decimal counts (THSPHTE-TCL_L0) [counts]
        e2l_L = array of calibration coefficients to convert the 'L' thermocouple
                engineering values to lab calibrated values.
        l2s_L = array of calibration coefficients to convert the 'L' thermocouple
                lab calibrated values to scientific values.

    References:

        OOI (2014). Data Product Specification for Vent Fluid Temperature from
            THSPH. Document Control Number 1341-00120.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00120_Data_Product_Specification_THSPHTE_OOI.pdf)
    """
    # convert raw decimal output to lab calibrated values [mV]
    V_tc_actual_L = sfl_thsph_temp_labcal_l(tc_rawdec_L, e2l_L)

    # convert lab calibrated values to scientific values [degC]
    T_tc_L = eval_poly(V_tc_actual_L, l2s_L)

    return T_tc_L


def sfl_thsph_temp_ref(ts_rawdec_r, e2l_r, l2s_r):
    """
    Description:

        OOI Level 1 THSPH data product THSPHTE-REF (reference thermistor
        temperature), which is calculated using data from the Hydrothermal
        Vent Fluid In-situ Chemistry (THSPH) instrument, series A (one series
        for all instruments).

    Implemented by:

        2014-05-01: Russell Desiderio. Initial Code
        2014-06-30: Russell Desiderio. DPS modifications to cal equations implemented.

    Usage:

        T_ts_r = sfl_thsph_temp_ref(ts_rawdec_r, e2l_r, l2s_r)

            where

        T_ts_r = reference thermistor temperature THSPHTE-REF_L1 [deg_C]
        ts_rawdec_r = reference thermistor, decimal counts (THSPHTE-REF_L0) [counts]
        e2l_r = array of calibration coefficients to convert the 'r' thermistor
                engineering values to lab calibrated values.
        l2s_r = array of calibration coefficients to convert the 'r' thermistor
                lab calibrated values to scientific values.

    References:

        OOI (2014). Data Product Specification for Vent Fluid Temperature from
            THSPH. Document Control Number 1341-00120.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00120_Data_Product_Specification_THSPHTE_OOI.pdf)
    """
    # reset exception handling so that when divide by zeros and np.log(x<=0) are
    # encountered, the warnings are suppressed.
    oldsettings = np.seterr(divide='ignore', invalid='ignore')

    # convert raw decimal output to engineering values [ohms]
    ts_rawdec_r_scaled = ts_rawdec_r * 0.125
    R_ts_eng_r = 10000.0 * ts_rawdec_r_scaled / (2048.0 - ts_rawdec_r_scaled)

    # convert engineering values to lab calibrated values [ohms]
    R_ts_actual_r = eval_poly(R_ts_eng_r, e2l_r)

    # convert lab calibrated values to scientific values [degC]
    pval = eval_poly(np.log(R_ts_actual_r), l2s_r)
    T_ts_r = 1.0 / pval - 273.15

    # restore default exception handling settings
    np.seterr(**oldsettings)

    # trap out possible occurrences of nans and infs due to log(val<=0) and divide by zero.
    # nans and infs will be returned only if the variables are elements of numpy arrays.
    # do not use np.isinf, it does not work as desired if its argument is np.nan.
    T_ts_r[~np.isfinite(T_ts_r)] = fill_value

    return T_ts_r


def sfl_thsph_temp_int(ts_rawdec_b, e2l_b, l2s_b):
    """
    Description:

        OOI Level 1 THSPH data product THSPHTE-INT (internal board thermistor
        temperature), which is calculated using data from the Hydrothermal
        Vent Fluid In-situ Chemistry (THSPH) instrument, series A (one series
        for all instruments).

    Implemented by:

        2014-05-01: Russell Desiderio. Initial Code
        2014-06-30: Russell Desiderio. DPS modifications to cal equations implemented.

    Usage:

        T_ts_b = sfl_thsph_temp_int(ts_rawdec_b, e2l_b, l2s_b)

            where

        T_ts_b = board thermistor temperature THSPHTE-INT_L1 [deg_C]
        ts_rawdec_b = board thermistor, decimal counts (THSPHTE-INT_L0) [counts]
        e2l_b = array of calibration coefficients to convert the 'b' thermistor
                engineering values to lab calibrated values.
        l2s_b = array of calibration coefficients to convert the 'b' thermistor
                lab calibrated values to scientific values.

    References:

        OOI (2014). Data Product Specification for Vent Fluid Temperature from
            THSPH. Document Control Number 1341-00120.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00120_Data_Product_Specification_THSPHTE_OOI.pdf)
    """
    # reset exception handling so that when divide by zeros and np.log(x<=0) are
    # encountered, the warnings are suppressed.
    oldsettings = np.seterr(divide='ignore', invalid='ignore')

    # convert raw decimal output to engineering values [ohms]
    ts_rawdec_b_scaled = ts_rawdec_b * 0.125
    R_ts_eng_b = 10000.0 * ts_rawdec_b_scaled / (2048.0 - ts_rawdec_b_scaled)

    # convert engineering values to lab calibrated values [ohms]
    R_ts_actual_b = eval_poly(R_ts_eng_b, e2l_b)

    # convert lab calibrated values to scientific values [degC]
    pval = eval_poly(np.log(R_ts_actual_b), l2s_b)

    T_ts_b = 1.0 / pval - 273.15

    # restore default exception handling settings
    np.seterr(**oldsettings)

    # trap out possible occurrences of nans and infs due to log(val<=0) and divide by zero.
    # nans and infs will be returned only if the variables are elements of numpy arrays.
    # do not use np.isinf, it does not work as desired if its argument is np.nan.
    T_ts_b[~np.isfinite(T_ts_b)] = fill_value

    return T_ts_b


def sfl_thsph_temp_labcal_h(tc_rawdec_H, e2l_H):
    """
    Description:

        OOI Level 1 THSPH data products THSPHTE-TCH and THSPHTE-TH require this subfunction,
        which calculates lab calibrated mV values for the 'H' thermistor.

    Implemented by:

        2014-06-30: Russell Desiderio. Initial Code

    Usage:

        V_tc_labcal_H = sfl_thsph_temp_tch(tc_rawdec_H, e2l_H)

            where

        V_tc_labcal_H = intermediate variable used in calculation of THSPHTE-TCH and THSPHTE-TH.
        tc_rawdec_H = "H" thermocouple, decimal counts (THSPHTE-TCH_L0) [counts]
        e2l_H = array of calibration coefficients to convert the 'H' thermocouple
                engineering values to lab calibrated values.

    References:

        OOI (2014). Data Product Specification for Vent Fluid Temperature from
            THSPH. Document Control Number 1341-00120.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00120_Data_Product_Specification_THSPHTE_OOI.pdf)
    """
    # convert raw decimal output to engineering values [mV]
    # leave constants as is for clarity
    V_tc_eng_H = (tc_rawdec_H * 0.25 - 1024.0) / 61.606

    # convert engineering values to lab calibrated values [mV]
    V_tc_labcal_H = eval_poly(V_tc_eng_H, e2l_H)

    return V_tc_labcal_H


def sfl_thsph_temp_labcal_l(tc_rawdec_L, e2l_L):
    """
    Description:

        OOI Level 1 THSPH data products THSPHTE-TCL and THSPHTE-TL require this subfunction,
        which calculates lab calibrated mV values for the 'L' thermistor.

    Implemented by:

        2014-06-30: Russell Desiderio. Initial Code

    Usage:

        V_tc_labcal_L = sfl_thsph_temp_tcl(tc_rawdec_L, e2l_L)

            where

        V_tc_labcal_L = intermediate variable used in calculation of THSPHTE-TCL and THSPHTE-TL.
        tc_rawdec_L = "L" thermocouple, decimal counts (THSPHTE-TCL_L0) [counts]
        e2l_L = array of calibration coefficients to convert the 'L' thermocouple
                engineering values to lab calibrated values.

    References:

        OOI (2014). Data Product Specification for Vent Fluid Temperature from
            THSPH. Document Control Number 1341-00120.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00120_Data_Product_Specification_THSPHTE_OOI.pdf)
    """
    # convert raw decimal output to engineering values [mV]
    # leave constants as is for clarity
    V_tc_eng_L = (tc_rawdec_L * 0.25 - 1024.0) / 61.606

    # convert engineering values to lab calibrated values [mV]
    V_tc_labcal_L = eval_poly(V_tc_eng_L, e2l_L)

    return V_tc_labcal_L


def eval_poly(x, c):
    """
    Description:

        Calculates polynomial values for use with THSPH data products using the
        Horner algorithm. All coefficient sets are 5th order (6 terms), so that
        this function is written to be "vectorized" for speed for multiple data
        sets.

        The documentation for the numpy v1.7 function to evaluate polynomials
        was only available in draft form; plus, it won't handle "vectorized"
        calibration coefficients (2D arrays, in which each row is a separate
        set of calibration coeffs).

        The standard convention of storing polynomial coefficients in an array
        is used, namely, the highest order coefficient is the first element and
        coefficients are stored in descending order.

    Implemented by:

        2014-05-01: Russell Desiderio. Initial Code (no arrays, no recursion).
        2014-07-02: Russell Desiderio. 2D calcoeff array implementation.

    Usage:

        value = eval_poly(x, c)

            where

        value  = c[:,0]*x^(5) + c[:,1]*x^(4) + ... c[:,4]*x + c[:,5]

        x = the argument(s) of the polynomial to be evaluated; can be a scalar or vector
        c = array containing the polynomial coefficients:
                 if x is a scalar, then c is a vector.
                 if x is a vector, then c is a 2D array, where each row j is a set of
                                   polynomial coefficients associated with x[j].
    """
    # the "c = np.atleast_2d(c)" statement is necessary so that both single and
    # "vectorized" (in the OOI CI sense) calls to the eval_poly subroutine work.
    c = np.atleast_2d(c)

    val = c[:, 5] + x * (c[:, 4] + x * (c[:, 3] + x * (c[:, 2] + x * (c[:, 1] + x * c[:, 0]))))

    return val


def sfl_trhph_vfltemp(V_ts, V_tc, tc_slope, ts_slope, c0=0.015, c1=0.0024, c2=7.00e-5, c3=-1.00e-6):
    """
    Description:

        OOI Level 1 Vent Fluid Temperature from TRHPH (TRHPHTE) data product,
        which is calculated using data from the Temperature Resistivity Probe
        (TRHPH) instrument placed in a high temperature hydrothermal vent.

    Implemented by:

        2013-05-01: Christopher Wingard. Initial Code
        2014-02-27: Russell Desiderio. Added documentation.
                    Implemented Horner's method for polynomial calculation.

    Usage:

        T = sfl_trhph_vfltemp(V_ts, V_tc, tc_slope, ts_slope, c0, c1, c2, c3)

            where

        T = Vent fluid temperature from TRHPH (TRHPHTE_L1) [deg_C]
        V_ts = Thermistor voltage (TRHPHVS_L0) [volts]
        V_tc = Thermocouple voltage (TRHPHVC_L0) [volts]
        tc_slope = thermocople slope laboratory calibration coefficients
        ts_slope = thermistor slope laboratory calibration coefficients
        c0 = coefficient from 3rd degree polynomial fit of laboratory
            calibration correction curve (not expected to change).
        c1 = coefficient from 3rd degree polynomial fit of laboratory
            calibration correction curve (not expected to change).
        c2 = coefficient from 3rd degree polynomial fit of laboratory
            calibration correction curve (not expected to change).
        c3 = coefficient from 3rd degree polynomial fit of laboratory
            calibration correction curve (not expected to change).

    References:

        OOI (2012). Data Product Specification for Vent Fluid Temperature from
            TRHPH. Document Control Number 1341-00150.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00150_Data_Product_Spec_TRHPHTE_OOI.pdf)
    """
    # Test if polynomial coefficients are scalars (set via defaults), set to
    # same size as other inputs if required. Assumes if 'a' is a default, they
    # all are.
    if np.isscalar(c0):
        c0 = np.tile(c0, (V_ts.shape))

    if np.isscalar(c1):
        c1 = np.tile(c1, (V_ts.shape))

    if np.isscalar(c2):
        c2 = np.tile(c2, (V_ts.shape))

    if np.isscalar(c3):
        c3 = np.tile(c3, (V_ts.shape))

    # raw thermistor temperature
    T_ts = 27.50133 - 17.2658 * V_ts + 15.83424 / V_ts

    # where V_tc is less than or equal to 0, T = T_ts, otherwise...
    T = T_ts

    # Adjust raw thermistor temperature when V_tc is greater than 0 and ...
    tFlag = (V_tc > 0) & (T_ts > 10)  # T_ts is greater than 10
    poly = (c3[tFlag] * T_ts[tFlag]**3 + c2[tFlag] * T_ts[tFlag]**2 +
            c1[tFlag] * T_ts[tFlag] + c0[tFlag])
    T[tFlag] = (V_tc[tFlag] + poly) * 244.97

    tFlag = (V_tc > 0) & ((T_ts > 0) & (T_ts <= 10))  # T_ts is greater than 0 and less than 10
    T[tFlag] = (V_tc[tFlag] + V_tc[tFlag] * 244.97 * tc_slope[tFlag] + T_ts[tFlag]
                * ts_slope[tFlag]) * 244.97

    return T


def sfl_trhph_vfl_thermistor_temp(V_ts):
    """
    Description:

        Calculates T_S, which is an intermediate data product (not a core data
        product) requested by the authors of the TRHPHTE DPS. It is the
        instrument's thermistor temperature, useful as an important instrument
        diagnostic. It is the same variable as T_ts in the function
        sfl_trhph_vfltemp.

    Implemented by:

        2014-02-28: Russell Desiderio. Initial Code

    Usage:

        T_ts = sfl_trhph_vfl_thermistor_temp(V_ts)

            where

        T_ts = Thermistor reference temperature [deg_C]
        V_ts = Thermistor voltage (TRHPHVS_L0) [volts]

    References:

        OOI (2012). Data Product Specification for Vent Fluid Temperature from
            TRHPH. Document Control Number 1341-00150.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00150_Data_Product_Spec_TRHPHTE_OOI.pdf)
    """
    # thermistor temperature
    T_ts = 27.50133 - 17.2658 * V_ts + 15.83424 / V_ts
    return T_ts


def sfl_trhph_vflorp(V, offset, gain):
    """
    Description:

        OOI Level 1 Vent Fluid Oxidation-Reduction Potential (ORP) from TRHPH
        (TRHPHEH_L1) data product, which is calculated using data from the
        Resistivity- Temperature Probe (TRHPH) instrument.

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
        2014-02-28: Russell Desiderio. Modified code to better handle nans and
                    fill_values. Added more documentation to algorithm.
        2014-03-10: Russell Desiderio. Removed unnecessary np.vectorized
                    wrapper function. Improved speed by removing a temperature
                    conditional statement from inside the for loop and
                    incorporating it into the range of the for loop.
        2014-03-26: Russell Desiderio. Incorporated optimization due to Chris
                    Fortin: calculate Ccurve using scalar T instead of a vector
                    of constant T values. Sped up execution by factor of 5.

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
    # load sfl_functions_surface.py This loads the 3-dimensional calibration
    # surface of temperature, chloride, and conductivity reproduced as numpy
    # arrays from Larson_2007surface.mat.
    #
    # imported at module top
    # from ion_functions.data.sfl_functions_surface import tdat, sdat, cdat

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
    # always give nan values for Cl as is required. Since Cl has been
    # initialized to nan values, iterate only over good T values, which also
    # improves speed.
    for ii in np.where(np.logical_and(T >= min(tdat), T <= max(tdat)))[0]:
        # find conductivity Ccurve as f(T=constant, chloride).
        Ccurve = f(T[ii], Scurve)
        # now interpolate measured conductivity C into (Ccurve,Scurve) to get
        # Cl. the conditional statement is in the DPS and therefore retained.
        if (np.all(np.isfinite(Ccurve))):
            Cl[ii] = np.interp(C[ii], Ccurve[0], Scurve, left=np.nan, right=np.nan)

    # change units to mmol/kg; round to required # of sigfigs as specified in
    # the DPS
    Cl = np.round(Cl * 1000.)
    Cl[np.isnan(Cl)] = fill_value

    return Cl


def sfl_sflpres_l1(p_psia):
    """
    Description:

        The OOI Level 1 Seafloor Pressure core data products, SFLPRES and
        sub-parameters SFLPRES-RTIME, SFLPRES-TIDE, and SFLPRES-WAVE, are
        created from the Sea-Bird Electronics SBE 26plus member of the Pressure
        SF (PRESF) family of instruments by either a) polling, in real-time,
        for L0 ASCII text format data output and converting from psia to
        decibar units or b) converting, after instrument recovery, L0
        hexadecimal pressure data into decimal format and the resulting tide
        and wave pressure data in psia to decibar units.

    Implemented by:

        2014-01-31: Craig Risien. Initial Code

    Usage:

        sflpres_l1 = sfl_sflpres_l1(p_psia):

        Scaling: To convert from psia to dbar, use the Sea-Bird-specified
        conversion: p_dbar = 0.689475728 * (p_psia)

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
