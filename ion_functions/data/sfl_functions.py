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


def sfl_thsph_temp_th(tc_rawdec_H,
                      c0_e2l_H, c1_e2l_H, c2_e2l_H, c3_e2l_H, c4_e2l_H,
                      c0_l2s_H, c1_l2s_H, c2_l2s_H, c3_l2s_H, c4_l2s_H, c5_l2s_H,
                      c0_s2f_H, c1_s2f_H,
                      ts_rawdec_r,
                      c0_e2l_r, c1_e2l_r, c2_e2l_r, c3_e2l_r, c4_e2l_r,
                      c0_l2s_r, c1_l2s_r, c2_l2s_r, c3_l2s_r, c4_l2s_r):
    """
    Description:

        OOI Level 1 THSPH data product THSPHTE-TH (final temperature at position
        "H" near sample inlet), which is calculated using data from the Hydrothermal
        Vent Fluid In-situ Chemistry (THSPH) instrument, series A (one series for all
        instruments).

    Implemented by:

        2014-05-01: Russell Desiderio. Initial Code

    Usage:

        T_H = sfl_thsph_temp_tcl(tc_rawdec_H,
                       c0_e2l_H, c1_e2l_H, c2_e2l_H, c3_e2l_H, c4_e2l_H,
                       c0_l2s_H, c1_l2s_H, c2_l2s_H, c3_l2s_H, c4_l2s_H, c5_l2s_H,
                       c0_s2f_H, c1_s2f_H,
                                 ts_rawdec_r,
                       c0_e2l_r, c1_e2l_r, c2_e2l_r, c3_e2l_r, c4_e2l_r,
                       c0_l2s_r, c1_l2s_r, c2_l2s_r, c3_l2s_r, c4_l2s_r)

            where

        T_H = final temperature "H" near sample inlet THSPHTE-TH_L1 [deg_C]
        #
        tc_rawdec_H = "H" thermocouple, decimal counts (THSPHTE-TCH_L0) [counts]
        ### the e2l_H series of calibration coefficients convert the 'H' thermocouple
        ### engineering values to lab calibrated values using a 4th degree polynomial.
        c0_e2l_H = constant  coefficient for e2l polynomial for 'H' thermocouple
        c1_e2l_H = linear    coefficient for e2l polynomial for 'H' thermocouple
        c2_e2l_H = quadratic coefficient for e2l polynomial for 'H' thermocouple
        c3_e2l_H = cubic     coefficient for e2l polynomial for 'H' thermocouple
        c4_e2l_H = quartic   coefficient for e2l polynomial for 'H' thermocouple
        ### the l2s_H series of calibration coefficients convert the 'H' thermocouple
        ### lab calibrated values to scientific values using a 5th degree polynomial.
        c0_l2s_H = constant  coefficient for l2s polynomial for 'H' thermocouple
        c1_l2s_H = linear    coefficient for l2s polynomial for 'H' thermocouple
        c2_l2s_H = quadratic coefficient for l2s polynomial for 'H' thermocouple
        c3_l2s_H = cubic     coefficient for l2s polynomial for 'H' thermocouple
        c4_l2s_H = quartic   coefficient for l2s polynomial for 'H' thermocouple
        c5_l2s_H = quintic   coefficient for l2s polynomial for 'H' thermocouple
        ###
        c0_s2f_H = offset of final "H" temperature linear calibration
        c1_s2f_H = slope  of final "H" temperature linear calibration
        #
        ts_rawdec_r = reference thermistor, decimal counts (THSPHTE-REF_L0) [counts]
        ### the e2l_r series of calibration coefficients convert the 'r' thermistor
        ### engineering values to lab calibrated values using a 4th degree polynomial.
        c0_e2l_r = constant  coefficient for e2l polynomial for 'r' thermistor
        c1_e2l_r = linear    coefficient for e2l polynomial for 'r' thermistor
        c2_e2l_r = quadratic coefficient for e2l polynomial for 'r' thermistor
        c3_e2l_r = cubic     coefficient for e2l polynomial for 'r' thermistor
        c4_e2l_r = quartic   coefficient for e2l polynomial for 'r' thermistor
        ### the l2s_r series of calibration coefficients convert the 'r' thermistor
        ### lab calibrated values to scientific values using a 4th degree polynomial.
        c0_l2s_r = constant  coefficient for l2s polynomial for 'r' thermistor
        c1_l2s_r = linear    coefficient for l2s polynomial for 'r' thermistor
        c2_l2s_r = quadratic coefficient for l2s polynomial for 'r' thermistor
        c3_l2s_r = cubic     coefficient for l2s polynomial for 'r' thermistor
        c4_l2s_r = quartic   coefficient for l2s polynomial for 'r' thermistor

    References:

        OOI (2014). Data Product Specification for Vent Fluid Temperature from
            THSPH. Document Control Number 1341-00120.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00120_Data_Product_Specification_THSPHTE_OOI.pdf)
    """
    # calculate required intermediate products
    T_tc_H = sfl_thsph_temp_tcl(tc_rawdec_H,
                                c0_e2l_H, c1_e2l_H, c2_e2l_H, c3_e2l_H, c4_e2l_H,
                                c0_l2s_H, c1_l2s_H, c2_l2s_H, c3_l2s_H, c4_l2s_H, c5_l2s_H)
    T_ts_r = sfl_thsph_temp_ref(ts_rawdec_r,
                                c0_e2l_r, c1_e2l_r, c2_e2l_r, c3_e2l_r, c4_e2l_r,
                                c0_l2s_r, c1_l2s_r, c2_l2s_r, c3_l2s_r, c4_l2s_r)

    # apply final linear calibration
    T_H = c1_s2f_H * (T_ts_r + T_tc_H) + c0_s2f_H

    return T_H


def sfl_thsph_temp_tl(tc_rawdec_L,
                      c0_e2l_L, c1_e2l_L, c2_e2l_L, c3_e2l_L, c4_e2l_L,
                      c0_l2s_L, c1_l2s_L, c2_l2s_L, c3_l2s_L, c4_l2s_L, c5_l2s_L,
                      c0_s2f_L, c1_s2f_L,
                      ts_rawdec_r,
                      c0_e2l_r, c1_e2l_r, c2_e2l_r, c3_e2l_r, c4_e2l_r,
                      c0_l2s_r, c1_l2s_r, c2_l2s_r, c3_l2s_r, c4_l2s_r):
    """
    Description:

        OOI Level 1 THSPH data product THSPHTE-TL (final temperature at position
        "L" near vent), which is calculated using data from the Hydrothermal Vent
        Fluid In-situ Chemistry (THSPH) instrument, series A (one series for all
        instruments).

    Implemented by:

        2014-05-01: Russell Desiderio. Initial Code

    Usage:

        T_L = sfl_thsph_temp_tcl(tc_rawdec_L,
                       c0_e2l_L, c1_e2l_L, c2_e2l_L, c3_e2l_L, c4_e2l_L,
                       c0_l2s_L, c1_l2s_L, c2_l2s_L, c3_l2s_L, c4_l2s_L, c5_l2s_L,
                       c0_s2f_L, c1_s2f_L,
                                 ts_rawdec_r,
                       c0_e2l_r, c1_e2l_r, c2_e2l_r, c3_e2l_r, c4_e2l_r,
                       c0_l2s_r, c1_l2s_r, c2_l2s_r, c3_l2s_r, c4_l2s_r)

            where

        T_L = final temperature "L" near vent THSPHTE-TL_L1 [deg_C]
        #
        tc_rawdec_L = "L" thermocouple, decimal counts (THSPHTE-TCL_L0) [counts]
        ### the e2l_L series of calibration coefficients convert the 'L' thermocouple
        ### engineering values to lab calibrated values using a 4th degree polynomial.
        c0_e2l_L = constant  coefficient for e2l polynomial for 'L' thermocouple
        c1_e2l_L = linear    coefficient for e2l polynomial for 'L' thermocouple
        c2_e2l_L = quadratic coefficient for e2l polynomial for 'L' thermocouple
        c3_e2l_L = cubic     coefficient for e2l polynomial for 'L' thermocouple
        c4_e2l_L = quartic   coefficient for e2l polynomial for 'L' thermocouple
        ### the l2s_L series of calibration coefficients convert the 'L' thermocouple
        ### lab calibrated values to scientific values using a 5th degree polynomial.
        c0_l2s_L = constant  coefficient for l2s polynomial for 'L' thermocouple
        c1_l2s_L = linear    coefficient for l2s polynomial for 'L' thermocouple
        c2_l2s_L = quadratic coefficient for l2s polynomial for 'L' thermocouple
        c3_l2s_L = cubic     coefficient for l2s polynomial for 'L' thermocouple
        c4_l2s_L = quartic   coefficient for l2s polynomial for 'L' thermocouple
        c5_l2s_L = quintic   coefficient for l2s polynomial for 'L' thermocouple
        ###
        c0_s2f_L = offset of final "L" temperature linear calibration
        c1_s2f_L = slope  of final "L" temperature linear calibration
        #
        ts_rawdec_r = reference thermistor, decimal counts (THSPHTE-REF_L0) [counts]
        ### the e2l_r series of calibration coefficients convert the 'r' thermistor
        ### engineering values to lab calibrated values using a 4th degree polynomial.
        c0_e2l_r = constant  coefficient for e2l polynomial for 'r' thermistor
        c1_e2l_r = linear    coefficient for e2l polynomial for 'r' thermistor
        c2_e2l_r = quadratic coefficient for e2l polynomial for 'r' thermistor
        c3_e2l_r = cubic     coefficient for e2l polynomial for 'r' thermistor
        c4_e2l_r = quartic   coefficient for e2l polynomial for 'r' thermistor
        ### the l2s_r series of calibration coefficients convert the 'r' thermistor
        ### lab calibrated values to scientific values using a 4th degree polynomial.
        c0_l2s_r = constant  coefficient for l2s polynomial for 'r' thermistor
        c1_l2s_r = linear    coefficient for l2s polynomial for 'r' thermistor
        c2_l2s_r = quadratic coefficient for l2s polynomial for 'r' thermistor
        c3_l2s_r = cubic     coefficient for l2s polynomial for 'r' thermistor
        c4_l2s_r = quartic   coefficient for l2s polynomial for 'r' thermistor

    References:

        OOI (2014). Data Product Specification for Vent Fluid Temperature from
            THSPH. Document Control Number 1341-00120.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00120_Data_Product_Specification_THSPHTE_OOI.pdf)
    """
    # calculate required intermediate products
    T_tc_L = sfl_thsph_temp_tcl(tc_rawdec_L,
                                c0_e2l_L, c1_e2l_L, c2_e2l_L, c3_e2l_L, c4_e2l_L,
                                c0_l2s_L, c1_l2s_L, c2_l2s_L, c3_l2s_L, c4_l2s_L, c5_l2s_L)
    T_ts_r = sfl_thsph_temp_ref(ts_rawdec_r,
                                c0_e2l_r, c1_e2l_r, c2_e2l_r, c3_e2l_r, c4_e2l_r,
                                c0_l2s_r, c1_l2s_r, c2_l2s_r, c3_l2s_r, c4_l2s_r)

    # apply final linear calibration
    T_L = c1_s2f_L * (T_ts_r + T_tc_L) + c0_s2f_L

    return T_L


def sfl_thsph_temp_tch(tc_rawdec_H,
                       c0_e2l_H, c1_e2l_H, c2_e2l_H, c3_e2l_H, c4_e2l_H,
                       c0_l2s_H, c1_l2s_H, c2_l2s_H, c3_l2s_H, c4_l2s_H, c5_l2s_H):
    """
    Description:

        OOI Level 1 THSPH data product THSPHTE-TCH (intermediate thermocouple
        temperature at position "H"), which is calculated using data from the
        Hydrothermal Vent Fluid In-situ Chemistry (THSPH) instrument, series A
        (one series for all instruments).

    Implemented by:

        2014-05-01: Russell Desiderio. Initial Code

    Usage:

        T_tc_H = sfl_thsph_temp_tch(tc_rawdec_H,
                       c0_e2l_H, c1_e2l_H, c2_e2l_H, c3_e2l_H, c4_e2l_H,
                       c0_l2s_H, c1_l2s_H, c2_l2s_H, c3_l2s_H, c4_l2s_H, c5_l2s_H)

            where

        T_tc_H = intermediate thermocouple temperature "H" THSPHTE-TCH_L1 [deg_C]
        tc_rawdec_H = "H" thermocouple, decimal counts (THSPHTE-TCH_L0) [counts]
        ### the e2l_H series of calibration coefficients convert the 'H' thermocouple
        ### engineering values to lab calibrated values using a 4th degree polynomial.
        c0_e2l_H = constant  coefficient for e2l polynomial for 'H' thermocouple
        c1_e2l_H = linear    coefficient for e2l polynomial for 'H' thermocouple
        c2_e2l_H = quadratic coefficient for e2l polynomial for 'H' thermocouple
        c3_e2l_H = cubic     coefficient for e2l polynomial for 'H' thermocouple
        c4_e2l_H = quartic   coefficient for e2l polynomial for 'H' thermocouple
        ### the l2s_H series of calibration coefficients convert the 'H' thermocouple
        ### lab calibrated values to scientific values using a 5th degree polynomial.
        c0_l2s_H = constant  coefficient for l2s polynomial for 'H' thermocouple
        c1_l2s_H = linear    coefficient for l2s polynomial for 'H' thermocouple
        c2_l2s_H = quadratic coefficient for l2s polynomial for 'H' thermocouple
        c3_l2s_H = cubic     coefficient for l2s polynomial for 'H' thermocouple
        c4_l2s_H = quartic   coefficient for l2s polynomial for 'H' thermocouple
        c5_l2s_H = quintic   coefficient for l2s polynomial for 'H' thermocouple

    References:

        OOI (2014). Data Product Specification for Vent Fluid Temperature from
            THSPH. Document Control Number 1341-00120.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00120_Data_Product_Specification_THSPHTE_OOI.pdf)
    """
    # convert raw decimal output to engineering values [V]
    # leave constants as is for clarity
    V_tc_eng_H = (tc_rawdec_H * 0.25 - 1024.0) / 61606.0

    # convert engineering values to lab calibrated values [V]
    V_tc_actual_H = eval_poly(V_tc_eng_H,
                              c0_e2l_H, c1_e2l_H, c2_e2l_H, c3_e2l_H, c4_e2l_H)

    # convert to mV
    mV_tc_actual_H = 1000.0 * V_tc_actual_H

    # convert lab calibrated values to scientific values [degC]
    T_tc_H = eval_poly(mV_tc_actual_H,
                       c0_l2s_H, c1_l2s_H, c2_l2s_H, c3_l2s_H, c4_l2s_H, c5_l2s_H)

    return T_tc_H


def sfl_thsph_temp_tcl(tc_rawdec_L,
                       c0_e2l_L, c1_e2l_L, c2_e2l_L, c3_e2l_L, c4_e2l_L,
                       c0_l2s_L, c1_l2s_L, c2_l2s_L, c3_l2s_L, c4_l2s_L, c5_l2s_L):
    """
    Description:

        OOI Level 1 THSPH data product THSPHTE-TCL (intermediate thermocouple
        temperature at position "L"), which is calculated using data from the
        Hydrothermal Vent Fluid In-situ Chemistry (THSPH) instrument, series A
        (one series for all instruments).

    Implemented by:

        2014-05-01: Russell Desiderio. Initial Code

    Usage:

        T_tc_L = sfl_thsph_temp_tcl(tc_rawdec_L,
                       c0_e2l_L, c1_e2l_L, c2_e2l_L, c3_e2l_L, c4_e2l_L,
                       c0_l2s_L, c1_l2s_L, c2_l2s_L, c3_l2s_L, c4_l2s_L, c5_l2s_L)

            where

        T_tc_L = intermediate thermocouple temperature "L" THSPHTE-TCL_L1 [deg_C]
        tc_rawdec_L = "L" thermocouple, decimal counts (THSPHTE-TCL_L0) [counts]
        ### the e2l_L series of calibration coefficients convert the 'L' thermocouple
        ### engineering values to lab calibrated values using a 4th degree polynomial.
        c0_e2l_L = constant  coefficient for e2l polynomial for 'L' thermocouple
        c1_e2l_L = linear    coefficient for e2l polynomial for 'L' thermocouple
        c2_e2l_L = quadratic coefficient for e2l polynomial for 'L' thermocouple
        c3_e2l_L = cubic     coefficient for e2l polynomial for 'L' thermocouple
        c4_e2l_L = quartic   coefficient for e2l polynomial for 'L' thermocouple
        ### the l2s_L series of calibration coefficients convert the 'L' thermocouple
        ### lab calibrated values to scientific values using a 5th degree polynomial.
        c0_l2s_L = constant  coefficient for l2s polynomial for 'L' thermocouple
        c1_l2s_L = linear    coefficient for l2s polynomial for 'L' thermocouple
        c2_l2s_L = quadratic coefficient for l2s polynomial for 'L' thermocouple
        c3_l2s_L = cubic     coefficient for l2s polynomial for 'L' thermocouple
        c4_l2s_L = quartic   coefficient for l2s polynomial for 'L' thermocouple
        c5_l2s_L = quintic   coefficient for l2s polynomial for 'L' thermocouple

    References:

        OOI (2014). Data Product Specification for Vent Fluid Temperature from
            THSPH. Document Control Number 1341-00120.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00120_Data_Product_Specification_THSPHTE_OOI.pdf)
    """
    # convert raw decimal output to engineering values [V]
    # leave constants as is for clarity
    V_tc_eng_L = (tc_rawdec_L * 0.25 - 1024.0) / 61606.0

    # convert engineering values to lab calibrated values [V]
    V_tc_actual_L = eval_poly(V_tc_eng_L,
                              c0_e2l_L, c1_e2l_L, c2_e2l_L, c3_e2l_L, c4_e2l_L)

    # convert to mV
    mV_tc_actual_L = 1000.0 * V_tc_actual_L

    # convert lab calibrated values to scientific values [degC]
    T_tc_L = eval_poly(mV_tc_actual_L,
                       c0_l2s_L, c1_l2s_L, c2_l2s_L, c3_l2s_L, c4_l2s_L, c5_l2s_L)

    return T_tc_L


def sfl_thsph_temp_ref(ts_rawdec_r,
                       c0_e2l_r, c1_e2l_r, c2_e2l_r, c3_e2l_r, c4_e2l_r,
                       c0_l2s_r, c1_l2s_r, c2_l2s_r, c3_l2s_r, c4_l2s_r):
    """
    Description:

        OOI Level 1 THSPH data product THSPHTE-REF (reference thermistor
        temperature), which is calculated using data from the Hydrothermal
        Vent Fluid In-situ Chemistry (THSPH) instrument, series A (one series
        for all instruments).

    Implemented by:

        2014-05-01: Russell Desiderio. Initial Code

    Usage:

        T_ts_r = sfl_thsph_temp_ref(ts_rawdec_r,
                       c0_e2l_r, c1_e2l_r, c2_e2l_r, c3_e2l_r, c4_e2l_r,
                       c0_l2s_r, c1_l2s_r, c2_l2s_r, c3_l2s_r, c4_l2s_r)

            where

        T_ts_r = reference thermistor temperature THSPHTE-REF_L1 [deg_C]
        ts_rawdec_r = reference thermistor, decimal counts (THSPHTE-REF_L0) [counts]
        ### the e2l_r series of calibration coefficients convert the 'r' thermistor
        ### engineering values to lab calibrated values using a 4th degree polynomial.
        c0_e2l_r = constant  coefficient for e2l polynomial for 'r' thermistor
        c1_e2l_r = linear    coefficient for e2l polynomial for 'r' thermistor
        c2_e2l_r = quadratic coefficient for e2l polynomial for 'r' thermistor
        c3_e2l_r = cubic     coefficient for e2l polynomial for 'r' thermistor
        c4_e2l_r = quartic   coefficient for e2l polynomial for 'r' thermistor
        ### the l2s_r series of calibration coefficients convert the 'r' thermistor
        ### lab calibrated values to scientific values using a 4th degree polynomial.
        c0_l2s_r = constant  coefficient for l2s polynomial for 'r' thermistor
        c1_l2s_r = linear    coefficient for l2s polynomial for 'r' thermistor
        c2_l2s_r = quadratic coefficient for l2s polynomial for 'r' thermistor
        c3_l2s_r = cubic     coefficient for l2s polynomial for 'r' thermistor
        c4_l2s_r = quartic   coefficient for l2s polynomial for 'r' thermistor

    References:

        OOI (2014). Data Product Specification for Vent Fluid Temperature from
            THSPH. Document Control Number 1341-00120.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00120_Data_Product_Specification_THSPHTE_OOI.pdf)
    """
    # convert raw decimal output to engineering values [ohms]
    ts_rawdec_r_scaled = ts_rawdec_r * 0.125
    R_ts_eng_r = 10000.0 * ts_rawdec_r_scaled / (2048.0 - ts_rawdec_r_scaled)

    # convert engineering values to lab calibrated values [ohms]
    R_ts_actual_r = eval_poly(R_ts_eng_r,
                              c0_e2l_r, c1_e2l_r, c2_e2l_r, c3_e2l_r, c4_e2l_r)

    # convert lab calibrated values to scientific values [degC]
    T_ts_r = eval_poly(R_ts_actual_r,
                       c0_l2s_r, c1_l2s_r, c2_l2s_r, c3_l2s_r, c4_l2s_r)

    return T_ts_r


def sfl_thsph_temp_int(ts_rawdec_b,
                       c0_e2l_b, c1_e2l_b, c2_e2l_b, c3_e2l_b, c4_e2l_b,
                       c0_l2s_b, c1_l2s_b, c2_l2s_b, c3_l2s_b, c4_l2s_b):
    """
    Description:

        OOI Level 1 THSPH data product THSPHTE-INT (internal board thermistor
        temperature), which is calculated using data from the Hydrothermal
        Vent Fluid In-situ Chemistry (THSPH) instrument, series A (one series
        for all instruments).

    Implemented by:

        2014-05-01: Russell Desiderio. Initial Code

    Usage:

        T_ts_b = sfl_thsph_temp_int(ts_rawdec_b,
                       c0_e2l_b, c1_e2l_b, c2_e2l_b, c3_e2l_b, c4_e2l_b,
                       c0_l2s_b, c1_l2s_b, c2_l2s_b, c3_l2s_b, c4_l2s_b)

            where

        T_ts_b = board thermistor temperature THSPHTE-INT_L1 [deg_C]
        ts_rawdec_b = board thermistor, decimal counts (THSPHTE-INT_L0) [counts]
        ### the e2l_b series of calibration coefficients convert the 'b' thermistor
        ### engineering values to lab calibrated values using a 4th degree polynomial.
        c0_e2l_b = constant  coefficient for e2l polynomial for 'b' thermistor
        c1_e2l_b = linear    coefficient for e2l polynomial for 'b' thermistor
        c2_e2l_b = quadratic coefficient for e2l polynomial for 'b' thermistor
        c3_e2l_b = cubic     coefficient for e2l polynomial for 'b' thermistor
        c4_e2l_b = quartic   coefficient for e2l polynomial for 'b' thermistor
        ### the l2s_b series of calibration coefficients convert the 'b' thermistor
        ### lab calibrated values to scientific values using a 4th degree polynomial.
        c0_l2s_b = constant  coefficient for l2s polynomial for 'b' thermistor
        c1_l2s_b = linear    coefficient for l2s polynomial for 'b' thermistor
        c2_l2s_b = quadratic coefficient for l2s polynomial for 'b' thermistor
        c3_l2s_b = cubic     coefficient for l2s polynomial for 'b' thermistor
        c4_l2s_b = quartic   coefficient for l2s polynomial for 'b' thermistor

    References:

        OOI (2014). Data Product Specification for Vent Fluid Temperature from
            THSPH. Document Control Number 1341-00120.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00120_Data_Product_Specification_THSPHTE_OOI.pdf)
    """
    # convert raw decimal output to engineering values [ohms]
    ts_rawdec_b_scaled = ts_rawdec_b * 0.125
    R_ts_eng_b = 10000.0 * ts_rawdec_b_scaled / (2048.0 - ts_rawdec_b_scaled)

    # convert engineering values to lab calibrated values [ohms]
    R_ts_actual_b = eval_poly(R_ts_eng_b,
                              c0_e2l_b, c1_e2l_b, c2_e2l_b, c3_e2l_b, c4_e2l_b)

    # convert lab calibrated values to scientific values [degC]
    T_ts_b = eval_poly(R_ts_actual_b,
                       c0_l2s_b, c1_l2s_b, c2_l2s_b, c3_l2s_b, c4_l2s_b)

    return T_ts_b


def eval_poly(x, c0, c1, c2, c3, c4, c5=0.0):
    """
    Description:

        Calculates polynomial values for use with THSPSTE DPAs.

        The documentation for the numpy v1.7 function to evaluate polynomials was
        only available in draft form. Also, it uses the convention that the 0th
        coefficient multiplies the highest degree term, whereas the DPA uses the
        convention that the 0th coefficient is the constant coefficient, the 1st
        coefficient is the linear coefficient, etc. For clarity the present routine
        uses the convention followed in the DPA.

        Horner's method is used to evaluate the polynomial.

    Implemented by:

        2014-05-01: Russell Desiderio. Initial Code

    Usage:

        value = eval_poly(x, c0, c1, c2, c3, c4[, c5])

            where

        value = P(x) = c0 + c1*x + c2*x^2 + c3*x^3 + c4*x^4 + c5*x^5
    """
    return c0 + x * (c1 + x * (c2 + x * (c3 + x * (c4 + x * c5))))


def sfl_trhph_vfltemp(V_s, V_c, a, b, c, d, e):
    """
    Description:

        OOI Level 1 Vent Fluid Temperature from TRHPH (TRHPHTE) data product,
        which is calculated using data from the Temperature Resistivity Probe
        (TRHPH) instrument.

    Implemented by:

        2013-05-01: Christopher Wingard. Initial Code
        2014-02-27: Russell Desiderio. Added documentation.
                    Implemented Horner's method for polynomial calculation.

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
        2014-03-26: Russell Desiderio. Incorporated optimization due to Chris Fortin: calculate
                                           Ccurve using scalar T instead of a vector of constant
                                           T values. Sped up execution by factor of 5.

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
    # always give nan values for Cl as is required. Since Cl has been initialized
    # to nan values, iterate only over good T values, which also improves speed.
    for ii in np.where(np.logical_and(T >= min(tdat), T <= max(tdat)))[0]:
        # find conductivity Ccurve as f(T=constant, chloride).
        Ccurve = f(T[ii], Scurve)
        # now interpolate measured conductivity C into (Ccurve,Scurve) to get Cl.
        # the conditional statement is in the DPS and therefore retained.
        if (np.all(np.isfinite(Ccurve))):
            Cl[ii] = np.interp(C[ii], Ccurve[0], Scurve, left=np.nan, right=np.nan)

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
