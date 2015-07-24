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

from ion_functions.data.generic_functions import replace_fill_with_nan
# used by def sfl_trhph_chloride
from ion_functions.data.sfl_functions_surface import tdat, sdat, cdat

# .............................................................................
# THSPH data products: THSPHHC, THSPHHS, THSPHPH (4 PH products) ..............
# .............................................................................


def sfl_thsph_ph(counts_ysz, counts_agcl, temperature, e2l_ysz, e2l_agcl,
                 arr_hgo, arr_agcl, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3, chl):
    """
    Description:

        Calculates the THSPHPH-PH_L2 data product, one of the 4 THSPHPH data
        products for the THSPH instruments. The PH data product algorithm
        calculates pH assuming good chloride data is available from
        TRHPH (TRHPHCC_L2) and a working AgCl reference electrode.

    Implemented by:

        2014-07-08: Russell Desiderio. Initial Code.
        2015-07-24: Russell Desiderio. Incorporated calculate_vent_pH function.

    Usage:

        pH = sfl_thsph_ph(counts_ysz, counts_agcl, temperature, e2l_ysz, e2l_agcl,
                          arr_hgo, arr_agcl, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3, chl)

            where

        pH = vent fluid pH: THSPHPH-PH_L2) [unitless]
        counts_ysz =  raw data recorded by ysz electrode THSPHPH-YSZ_L0 [counts]
        counts_agcl =  raw data recorded by AgCl electrode THSPHPH-AGCL_L0 [counts]
        temperature = temperature near sample inlet THSPHTE-TH_L1 [deg_C].
        e2l_ysz =     array of 5th degree polynomial coefficients to convert the
                      ysz electrode engineering values to lab calibrated values.
        e2l_agcl =    array of 5th degree polynomial coefficients to convert the
                      agcl electrode engineering values to lab calibrated values.
        arr_hgo =     array of 5th degree polynomial coefficients to calculate the
                      electrode material response to temperature.
        arr_agcl =    array of 5th degree polynomial coefficients to calculate the
                      AgCl electrode material response to temperature.
        arr_tac = array containing the 5th degree polynomial coefficients to calculate tac (=tbc0).
        arr_tbc1 = array containing the 5th degree polynomial coefficients to calculate tbc1.
        arr_tbc2 = array containing the 5th degree polynomial coefficients to calculate tbc2.
        arr_tbc3 = array containing the 5th degree polynomial coefficients to calculate tbc3.
        chl = vent fluid chloride concentration from TRHPHCC_L2 [mmol kg-1].

    References:

        OOI (2014). Data Product Specification for Vent Fluid pH. Document Control
            Number 1341-00190. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI>> Controlled >> 1000 System Level >>
            1341-00190_Data_Product_Spec_THSPHPH_OOI.pdf)
    """
    # calculate lab calibrated electrode response [V]
    v_labcal_ysz = v_labcal(counts_ysz, e2l_ysz)
    # AgCl reference electrode
    v_labcal_agcl = v_labcal(counts_agcl, e2l_agcl)

    # calculate chloride activity
    act_chl = chloride_activity(temperature, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3, chl)

    pH = calculate_vent_pH(v_labcal_ysz, v_labcal_agcl, temperature, arr_hgo, arr_agcl, act_chl)

    return pH


def sfl_thsph_ph_acl(counts_ysz, counts_agcl, temperature, e2l_ysz, e2l_agcl,
                     arr_hgo, arr_agcl, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3):
    """
    Description:

        Calculates the THSPHPH-PH-ACL_L2 data product, one of the 4 THSPHPH
        data products for the THSPH instruments. The PH-ACL data product
        algorithm calculates pH assuming no good chloride data available from
        TRHPH (TRHPHCC_L2) (assumes instead a pre-determined chloride concentration
        which is set in the chloride_activity function). The data from the AgCl
        reference electrode is also assumed to be good and used in this
        calculation.

    Implemented by:

        2014-07-08: Russell Desiderio. Initial Code.
        2015-07-24: Russell Desiderio. Incorporated calculate_vent_pH function.

    Usage:

        pH = sfl_thsph_ph_acl(counts_ysz, counts_agcl, temperature, e2l_ysz, e2l_agcl,
                              arr_hgo, arr_agcl, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3)

            where

        pH = vent fluid pH: THSPHPH-PH-ACL_L2) [unitless]
        counts_ysz =  raw data recorded by ysz electrode THSPHPH-YSZ_L0 [counts]
        counts_agcl =  raw data recorded by AgCl electrode THSPHPH-AGCL_L0 [counts]
        temperature = temperature near sample inlet THSPHTE-TH_L1 [deg_C].
        e2l_ysz =     array of 5th degree polynomial coefficients to convert the
                      ysz electrode engineering values to lab calibrated values.
        e2l_agcl =    array of 5th degree polynomial coefficients to convert the
                      agcl electrode engineering values to lab calibrated values.
        arr_hgo =     array of 5th degree polynomial coefficients to calculate the
                      electrode material response to temperature.
        arr_agcl =    array of 5th degree polynomial coefficients to calculate the
                      AgCl electrode material response to temperature.
        arr_tac = array containing the 5th degree polynomial coefficients to calculate tac (=tbc0).
        arr_tbc1 = array containing the 5th degree polynomial coefficients to calculate tbc1.
        arr_tbc2 = array containing the 5th degree polynomial coefficients to calculate tbc2.
        arr_tbc3 = array containing the 5th degree polynomial coefficients to calculate tbc3.

    References:

        OOI (2014). Data Product Specification for Vent Fluid pH. Document Control
            Number 1341-00190. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI>> Controlled >> 1000 System Level >>
            1341-00190_Data_Product_Spec_THSPHPH_OOI.pdf)
    """
    # calculate lab calibrated electrode response [V]
    v_labcal_ysz = v_labcal(counts_ysz, e2l_ysz)
    # AgCl reference electrode
    v_labcal_agcl = v_labcal(counts_agcl, e2l_agcl)

    # chloride activity assuming the default value for chloride concentration
    # set in the chloride_activity subroutine
    act_chl = chloride_activity(temperature, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3)

    pH = calculate_vent_pH(v_labcal_ysz, v_labcal_agcl, temperature, arr_hgo, arr_agcl, act_chl)

    return pH


def sfl_thsph_ph_noref(counts_ysz, temperature, arr_agclref, e2l_ysz, arr_hgo,
                       arr_agcl, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3, chl):
    """
    Description:

        Calculates the THSPHPH-PH-NOREF_L2 data product, one of the 4 THSPHPH
        data products for the THSPH instruments. The PH-NOREF data product
        algorithm calculates pH assuming no good reference (AgCl) electrode data
        (uses instead a theoretical value calculated from vent temperature) and
        also uses (presumably good) chloride data from TRHPH (TRHPHCC_L2).

    Implemented by:

        2014-07-08: Russell Desiderio. Initial Code.
        2015-07-24: Russell Desiderio. Incorporated calculate_vent_pH function.

    Usage:

        pH = sfl_thsph_ph_noref(counts_ysz, temperature, arr_agclref, e2l_ysz, arr_hgo,
                                arr_agcl, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3, chl)

            where

        pH = vent fluid pH: THSPHPH-PH-NOREF_L2) [unitless]
        counts_ysz =  raw data recorded by ysz electrode THSPHPH-YSZ_L0 [counts]
        temperature = temperature near sample inlet THSPHTE-TH_L1 [deg_C].
        arr_agclref = array of 5th degree polynomial coefficients to calculate the
                      theoretical reference electrode potential, replacing measured
                      reference AgCl electrode potential values.
        e2l_ysz =     array of 5th degree polynomial coefficients to convert the
                      ysz electrode engineering values to lab calibrated values.
        arr_hgo =     array of 5th degree polynomial coefficients to calculate the
                      electrode material response to temperature.
        arr_agcl =    array of 5th degree polynomial coefficients to calculate the
                      AgCl electrode material response to temperature.
        arr_tac = array containing the 5th degree polynomial coefficients to calculate tac (=tbc0).
        arr_tbc1 = array containing the 5th degree polynomial coefficients to calculate tbc1.
        arr_tbc2 = array containing the 5th degree polynomial coefficients to calculate tbc2.
        arr_tbc3 = array containing the 5th degree polynomial coefficients to calculate tbc3.
        chl = vent fluid chloride concentration from TRHPHCC_L2 [mmol kg-1].

    References:

        OOI (2014). Data Product Specification for Vent Fluid pH. Document Control
            Number 1341-00190. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI>> Controlled >> 1000 System Level >>
            1341-00190_Data_Product_Spec_THSPHPH_OOI.pdf)
    """
    # calculate lab calibrated electrode response [V]
    v_labcal_ysz = v_labcal(counts_ysz, e2l_ysz)

    # theoretical reference value calculated from vent temperature
    e_refcalc = eval_poly(temperature, arr_agclref)
    # calculate chloride activity
    act_chl = chloride_activity(temperature, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3, chl)

    pH = calculate_vent_pH(v_labcal_ysz, e_refcalc, temperature, arr_hgo, arr_agcl, act_chl)

    return pH


def sfl_thsph_ph_noref_acl(counts_ysz, temperature, arr_agclref, e2l_ysz, arr_hgo,
                           arr_agcl, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3):
    """
    Description:

        Calculates the THSPHPH-PH-NOREF-ACL_L2 data product, one of the 4 THSPHPH
        data products for the THSPH instruments. The PH-NOREF-ACL data product
        algorithm calculates pH assuming no good reference (AgCl) electrode data
        (uses instead a theoretical value calculated from vent temperature) and
        assuming no good chloride data from TRHPH (TRHPHCC_L2) (assumes instead a
        pre-determined chloride concentration which is set in the chloride_activity
        function).

    Implemented by:

        2014-07-08: Russell Desiderio. Initial Code.
        2015-07-24: Russell Desiderio. Incorporated calculate_vent_pH function.

    Usage:

        pH = sfl_thsph_ph_noref_acl(counts_ysz, temperature, arr_agclref, e2l_ysz, arr_hgo,
                                    arr_agcl, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3)

            where

        pH = vent fluid pH: THSPHPH-PH-NOREF-ACL_L2) [unitless]
        counts_ysz =  raw data recorded by ysz electrode THSPHPH-YSZ_L0 [counts]
        temperature = temperature near sample inlet THSPHTE-TH_L1 [deg_C].
        arr_agclref = array of 5th degree polynomial coefficients to calculate the
                      theoretical reference electrode potential, replacing measured
                      reference AgCl electrode potential values.
        e2l_ysz =     array of 5th degree polynomial coefficients to convert the
                      ysz electrode engineering values to lab calibrated values.
        arr_hgo =     array of 5th degree polynomial coefficients to calculate the
                      electrode material response to temperature.
        arr_agcl =    array of 5th degree polynomial coefficients to calculate the
                      AgCl electrode material response to temperature.
        arr_tac = array containing the 5th degree polynomial coefficients to calculate tac (=tbc0).
        arr_tbc1 = array containing the 5th degree polynomial coefficients to calculate tbc1.
        arr_tbc2 = array containing the 5th degree polynomial coefficients to calculate tbc2.
        arr_tbc3 = array containing the 5th degree polynomial coefficients to calculate tbc3.

    References:

        OOI (2014). Data Product Specification for Vent Fluid pH. Document Control
            Number 1341-00190. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI>> Controlled >> 1000 System Level >>
            1341-00190_Data_Product_Spec_THSPHPH_OOI.pdf)
    """
    # calculate lab calibrated electrode response [V]
    v_labcal_ysz = v_labcal(counts_ysz, e2l_ysz)

    # theoretical reference value calculated from vent temperature
    e_refcalc = eval_poly(temperature, arr_agclref)
    # chloride activity assuming the default value for chloride concentration
    # set in the subroutine
    act_chl = chloride_activity(temperature, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3)

    pH = calculate_vent_pH(v_labcal_ysz, e_refcalc, temperature, arr_hgo, arr_agcl, act_chl)

    return pH


def calculate_vent_pH(e_ph, e_ref, temperature, arr_hgo, arr_agcl, act_chl):
    """
    Description:

        Worker function to calculate the vent fluid pH for the THSPH instruments. This
        function is called by
            sfl_thsph_ph
            sfl_thsph_ph_acl
            sfl_thsph_ph_noref
            sfl_thsph_ph_noref_acl.

    Implemented by:

        2015-07-24: Russell Desiderio. Initial Code.

    Usage:

        pH = calculate_vent_pH(e_ph, e_ref, temperature, arr_hgo, arr_agcl, act_chl)

            where

        pH = vent fluid pH
        e_ph = intermediate pH potential uncorrected for reference
        e_ref = reference pH potential, either measured or calculated
        temperature = temperature near sample inlet THSPHTE-TH_L1 [deg_C].
        arr_hgo = array of 5th degree polynomial coefficients to calculate the
                  electrode material response to temperature.
        arr_agcl = array of 5th degree polynomial coefficients to calculate the
                  AgCl electrode material response to temperature.
        act_chl = calculated chloride activity.

    References:

        OOI (2014). Data Product Specification for Vent Fluid pH. Document Control
            Number 1341-00190. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI>> Controlled >> 1000 System Level >>
            1341-00190_Data_Product_Spec_THSPHPH_OOI.pdf)
    """
    # fill value local to this function to avoid python warnings when nans are encountered
    # in boolean expressions. the masking will convert values derived from this local fill
    # back to nans.
    unphysical_pH_fill_value = -99999.0

    # calculate intermediate quantities that depend upon temperature
    e_nernst = nernst(temperature)
    e_hgo = eval_poly(temperature, arr_hgo)
    e_agcl = eval_poly(temperature, arr_agcl)

    # calculate pH potential
    e_phcalc = e_ph - e_ref
    # check for unphysical values as specified in the DPS.
    # logical indexing with boolean arrays is faster than integer indexing using np.where.
    # ok to apply mask at end of calculation.
    e_phcalc[np.isnan(e_phcalc)] = unphysical_pH_fill_value
    bad_eph_mask = np.logical_or(np.less(e_phcalc, -0.7), np.greater(e_phcalc, 0.0))

    # final data product calculation
    act_chl[act_chl <= 0.0] = np.nan  # trap out python warning
    pH = (e_phcalc - e_agcl + e_hgo) / e_nernst + np.log10(act_chl)

    # second check for unphysical values, as specified in the DPS
    pH[np.isnan(pH)] = unphysical_pH_fill_value
    bad_ph_mask = np.logical_or(np.less(pH, 3.0), np.greater(pH, 7.0))

    # set all out-of-range values to fill values
    pH[np.logical_or(bad_eph_mask, bad_ph_mask)] = np.nan

    return pH


def sfl_thsph_sulfide(counts_hs, counts_ysz, temperature, e2l_hs, e2l_ysz, arr_hgo,
                      arr_logkfh2g, arr_eh2sg, arr_yh2sg):
    """
    Description:

        Calculates the THSPHHS_L2 data product (hydrogen sulfide concentration) for
        the THSPH instruments from vent temperature and from data from its sulfide
        and YSZ electrodes. Note that the chemical formula for hydrogen is H2, and
        that for hydrogen sulfide is H2S; this could lead to confusion in the
        variable and array names from the DPS if care is not taken. Note also that
        this hydrogen sulfide DPA does use an intermediate data product and its
        'calibration' coefficients (hydrogen fugacity) that are also used in the
        hydrogen concentration DPA.

    Implemented by:

        2014-07-08: Russell Desiderio. Initial Code.

    Usage:

        h2s = sfl_thsph_sulfide(counts_hs, counts_ysz, temperature, e2l_hs, e2l_ysz,
                                arr_hgo, arr_logkfh2g, arr_eh2sg, arr_yh2sg)

            where

        h2s = hydrogen sulfide concentration at the vent: THSPHHS_L2 [mmol kg-1]
        counts_hs = raw data recorded by sulfide electrode THSPHHS_L0 [counts]
        counts_ysz = raw data recorded by ysz electrode THSPHPH-YSZ_L0 [counts]
        temperature = temperature near sample inlet THSPHTE-TH_L1 [deg_C].
        e2l_hs  = array of 5th degree polynomial coefficients to convert the
                  sulfide electrode engineering values to lab calibrated values.
        e2l_ysz = array of 5th degree polynomial coefficients to convert the
                  ysz electrode engineering values to lab calibrated values.
        arr_hgo = array of 5th degree polynomial coefficients to calculate the
                  electrode material response to temperature.
        arr_logkfh2g = array of 5th degree polynomial coefficients to calculate the
                  equilibrium hydrogen fugacity as a function of temperature.
        arr_eh2sg  = array of 5th degree polynomial coefficients to calculate the
                  theoretical potential of gas phase hydrogen sulfide to temperature;
                  in the current DPS, highest degree term is first order, so pad
                  this array with entries of zero [0., 0., 0., 0., c1, c0].
        arr_yh2sg  = array of 5th degree polynomial coefficients to calculate the
                  fugacity/concentration quotient yh2sg from hydrogen fugacity.

    References:

        OOI (2014). Data Product Specification for Vent Fluid Hydrogen Sulfide
            Concentration. Document Control Number 1341-00200.
            https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI>> Controlled >> 1000 System Level >>
            1341-00200_Data_Product_Spec_THSPHHS_OOI.pdf)
    """
    # calculate lab calibrated electrode responses [V]
    v_labcal_hs = v_labcal(counts_hs, e2l_hs)
    v_labcal_ysz = v_labcal(counts_ysz, e2l_ysz)

    # calculate intermediate products that depend upon temperature
    e_nernst = nernst(temperature)
    e_hgo = eval_poly(temperature, arr_hgo)
    e_h2sg = eval_poly(temperature, arr_eh2sg)
    log_kfh2g = eval_poly(temperature, arr_logkfh2g)
    # y_h2sg depends on temperature because hydrogen fugacity depends on temperature
    y_h2sg = eval_poly(log_kfh2g, arr_yh2sg)

    # explicitly follow the DPS calculation for clarity:

    # measured potential of the sulfide electrode [V]
    e_h2s = v_labcal_ysz - v_labcal_hs

    # (common) log of measured hydrogen sulfide fugacity
    log_fh2sg = 2.0 * (e_h2s - e_hgo + e_h2sg) / e_nernst

    # final data product, hydrogen sulfide concentration, [mmol/kg]
    # in the DPS, this is 1000 * 10^( logfh2sg - log( yh2sg ) )
    h2s = 1000.0 * (10.0 ** (log_fh2sg)) / y_h2sg

    return h2s


def sfl_thsph_hydrogen(counts_h2, counts_ysz, temperature, e2l_h2, e2l_ysz, arr_hgo,
                       arr_logkfh2g):
    """
    Description:

        Calculates the THSPHHC_L2 data product (hydrogen concentration) for the THSPH
        instruments from vent temperature and from data from its hydrogen and YSZ
        electrodes.

    Implemented by:

        2014-07-08: Russell Desiderio. Initial Code.

    Usage:

        h2 = sfl_thsph_hydrogen(counts_h2, counts_ysz, temperature, e2l_h2, e2l_ysz, arr_hgo,
                                arr_logkfh2g)

            where

        h2 = hydrogen concentration at the vent: THSPHHC_L2 [mmol kg-1]
        counts_h2 = raw data recorded by hydrogen electrode THSPHHC_L0 [counts]
        counts_ysz = raw data recorded by ysz electrode THSPHPH-YSZ_L0 [counts]
        temperature = temperature near sample inlet THSPHTE-TH_L1 [deg_C].
        e2l_h2  = array of 5th degree polynomial coefficients to convert the
                  hydrogen electrode engineering values to lab calibrated values.
        e2l_ysz = array of 5th degree polynomial coefficients to convert the
                  ysz electrode engineering values to lab calibrated values.
        arr_hgo = array of 5th degree polynomial coefficients to calculate the
                  electrode material response to temperature.
        arr_logkfh2g = array of 5th degree polynomial coefficients to calculate the
                  equilibrium hydrogen fugacity as a function of temperature.

        OOI (2014). Data Product Specification for Vent Fluid Hydrogen Concentration.
            Document Control Number 1341-00210. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI>> Controlled >> 1000 System Level >>
            1341-00210_Data_Product_Spec_THSPHHC_OOI.pdf)
    """
    # calculate lab calibrated electrode responses [V]
    v_labcal_h2 = v_labcal(counts_h2, e2l_h2)
    v_labcal_ysz = v_labcal(counts_ysz, e2l_ysz)

    # calculate intermediate products that depend upon temperature
    e_nernst = nernst(temperature)
    e_hgo = eval_poly(temperature, arr_hgo)
    log_kfh2g = eval_poly(temperature, arr_logkfh2g)

    # explicitly follow the DPS calculation for clarity:

    # measured potential of the h2 electrode [V]
    e_h2 = v_labcal_ysz - v_labcal_h2

    # (common) log of measured hydrogen fugacity
    log_fh2g = 2.0 * (e_h2 - e_hgo) / e_nernst

    # final data product, hydrogen concentration, [mmol/kg]
    h2 = 1000.0 * (10.0 ** (log_fh2g - log_kfh2g))

    return h2


def chloride_activity(temperature, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3, chloride=250.0):
    """
    Description:

        Subfunction to calculate the chloride activity as a function of temperature
        needed by the THSPHPH_L2 data products for the THSPH instruments. The chloride
        value can either come from the TRHPHCC_L2 data product or the default value
        of 250.0 mmol/kg can be used.

    Implemented by:

        2014-07-08: Russell Desiderio. Initial Code.

    Usage:

        act_chl = chloride_activity(temperature, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3[, chloride])

            where

        act_chl = calculated chloride activity.
        temperature = temperature near sample inlet THSPHTE-TH_L1 [deg_C].
        arr_tac = array containing the 5th degree polynomial coefficients to calculate tac (=tbc0).
        arr_tbc1 = array containing the 5th degree polynomial coefficients to calculate tbc1.
        arr_tbc2 = array containing the 5th degree polynomial coefficients to calculate tbc2.
        arr_tbc3 = array containing the 5th degree polynomial coefficients to calculate tbc3.
        chloride [optional] = if specified, vent fluid chloride concentration from TRHPH
                              (TRHPHCC_L2) [mmol kg-1], else a value of 250.0 mmol/kg will be used.
    """
    # find number of data packets to be processed;
    # this also works if temperature is not an np.array.
    nvalues = np.array([temperature]).shape[-1]

    # change units of chloride from mmol/kg to mol/kg
    chloride = chloride/1000.0

    # if chloride is not given in the argument list,
    # replicate its default value into a vector with
    # the same number of elements as temperature;
    # do so without using a conditional
    nreps = nvalues / np.array([chloride]).shape[-1]
    chloride = np.tile(chloride, nreps)

    # calculate the 4 coefficients needed to calculate the chloride activity from temperature
    tbc0 = eval_poly(temperature, arr_tac)
    tbc1 = eval_poly(temperature, arr_tbc1)
    tbc2 = eval_poly(temperature, arr_tbc2)
    tbc3 = eval_poly(temperature, arr_tbc3)

    # form these coeffs into a 2D array for the eval_poly routine.
    # need to pad the first two columns with zeros
    zeros = np.array([np.tile(0.0, nvalues)]).T
    arr_chloride_coeff = np.hstack((zeros, zeros, tbc3[:, np.newaxis], tbc2[:, np.newaxis],
                                    tbc1[:, np.newaxis], tbc0[:, np.newaxis]))

    # evaluate the activity
    act_chl = eval_poly(chloride, arr_chloride_coeff)

    return act_chl


def v_labcal(counts, array_e2l_coeff):
    """
    Description:

        Calculates any one of the 4 "lab calibrated" voltages from electrode sensor
        (not thermistor nor thermocouple) raw data used in the THSPH instruments.
        For use with the THSPH L2 data products (THSPHHC, THSPHHS, THSPHPH).

    Implemented by:

        2014-07-08: Russell Desiderio. Initial Code.
        2015-07-22: Russell Desiderio. Added call to replace_fill_with_nan.

    Usage:

        v_labcal_electrode = v_labcal(counts, array_e2l_coeff)

            where

        v_labcal_electrode = lab calibrated value ("V_actual" in DPS) for the electrode [V].
        counts = L0 output of one of the 4 electrodes:
                 THSPHPH-YSZ_L0, THSPHPH-AGC_L0, THSPHHC_L0, THSPHHS_L0 [decimal counts].
        array_e2l_coeff = 6 element array containing the 5th degree polynomial calibration
                          coefficients for the electrode for which the lab cal values are
                          desired. The coefficients are assumed to be stored in descending
                          order.

    Notes:

        All the THSPH L2 data products call v_labcal to process raw count input data. The
        action of the replace_fill_with_nan call in this code therefore replaces system fill
        values with nans for all of these DPAs.

    """
    counts = replace_fill_with_nan(None, counts)

    # transform decimal counts to engineering values [volts]
    v_eng = (counts * 0.25 - 2048.0) / 1000.0

    # transform engineering values to lab calibrated values [volts]
    v_labcal_electrode = eval_poly(v_eng, array_e2l_coeff)

    # in the DPSs, these values are designated as "V_actual"
    return v_labcal_electrode


def nernst(temperature):
    """
    Description:

        Calculates the value of the temperature dependent term of the Nernst
        equation to provide the link between measured electrode potentials and
        concentration. For use with the THSPH L2 data products (THSPHHC, THSPHHS,
        THSPHPH), all of which use electrodes to provide the raw data. The
        temperature to be used is specified in the DPSs to be THSPHTE-TH.

    Implemented by:

        2014-07-08: Russell Desiderio. Initial Code.

    Usage:

        e_nernst = nernst(temperature)

            where

        e_nernst = value of the temperature dependent term of the Nernst equation [V]
        temperature = temperature near sample inlet THSPHTE-TH_L1 [deg_C]
    """
    # e_nernst = ln(10) * (gas constant) * (temperature, Kelvin)/(Faraday's constant)
    #          = 2.30259 * 8.31446 [J/mole/K] / 96485.3 [coulombs/mole] * (T + 273.15)
    return 1.9842e-4 * (temperature + 273.15)


# .............................................................................
# THSPH data products: THSPHTE -TH, -TL, -TCH, -TCL, -REF, -INT ...............
# .............................................................................


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
    V_ts_r = eval_poly(T_ts_r, s2v_r)

    # Correct thermocouple temperature to account for offset from cold junction as
    # measured by the reference thermistor
    T_H = eval_poly((V_tc_actual_H + V_ts_r), l2s_H)

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
    V_ts_r = eval_poly(T_ts_r, s2v_r)

    # Correct thermocouple temperature to account for offset from cold junction as
    # measured by the reference thermistor
    T_L = eval_poly((V_tc_actual_L + V_ts_r), l2s_L)

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
        2015-07-24: Russell Desiderio. Added call to replace_fill_with_nan.
                                       Cleaned up error-checking.

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
    ts_rawdec_r = replace_fill_with_nan(None, ts_rawdec_r)

    # convert raw decimal output to engineering values [ohms]
    ts_rawdec_r_scaled = ts_rawdec_r * 0.125
    denom = 2048.0 - ts_rawdec_r_scaled
    denom[denom == 0.0] = np.nan
    R_ts_eng_r = 10000.0 * ts_rawdec_r_scaled / denom

    # convert engineering values to lab calibrated values [ohms]
    R_ts_actual_r = eval_poly(R_ts_eng_r, e2l_r)

    # convert lab calibrated values to scientific values [degC]
    R_ts_actual_r[R_ts_actual_r <= 0.0] = np.nan
    pval = eval_poly(np.log(R_ts_actual_r), l2s_r)
    T_ts_r = 1.0 / pval - 273.15

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
        2015-07-24: Russell Desiderio. Added call to replace_fill_with_nan.
                                       Cleaned up error-checking.

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
    ts_rawdec_b = replace_fill_with_nan(None, ts_rawdec_b)

    # convert raw decimal output to engineering values [ohms]
    ts_rawdec_b_scaled = ts_rawdec_b * 0.125
    denom = 2048.0 - ts_rawdec_b_scaled
    denom[denom == 0.0] = np.nan
    R_ts_eng_b = 10000.0 * ts_rawdec_b_scaled / denom

    # convert engineering values to lab calibrated values [ohms]
    R_ts_actual_b = eval_poly(R_ts_eng_b, e2l_b)

    # convert lab calibrated values to scientific values [degC]
    R_ts_actual_b[R_ts_actual_b <= 0.0] = np.nan
    pval = eval_poly(np.log(R_ts_actual_b), l2s_b)

    T_ts_b = 1.0 / pval - 273.15

    return T_ts_b


def sfl_thsph_temp_labcal_h(tc_rawdec_H, e2l_H):
    """
    Description:

        OOI Level 1 THSPH data products THSPHTE-TCH and THSPHTE-TH require this subfunction,
        which calculates lab calibrated mV values for the 'H' thermistor.

    Implemented by:

        2014-06-30: Russell Desiderio. Initial Code
        2015-07-24: Russell Desiderio. Added call to replace_fill_with_nan.

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
    tc_rawdec_H = replace_fill_with_nan(None, tc_rawdec_H)

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
        2015-07-24: Russell Desiderio. Added call to replace_fill_with_nan.

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
    tc_rawdec_L = replace_fill_with_nan(None, tc_rawdec_L)

    # convert raw decimal output to engineering values [mV]
    # leave constants as is for clarity
    V_tc_eng_L = (tc_rawdec_L * 0.25 - 1024.0) / 61.606

    # convert engineering values to lab calibrated values [mV]
    V_tc_labcal_L = eval_poly(V_tc_eng_L, e2l_L)

    return V_tc_labcal_L


# .............................................................................
# THSPH data products: eval_poly, used for all THSPH data products ............
# .............................................................................


def eval_poly(x, c):
    """
    Description:

        Calculates polynomial values for use with THSPH data products using the
        Horner algorithm. All coefficient sets are 5th degree (6 terms), so that
        this function is written to be "vectorized" for speed for multiple data
        sets.

        The documentation for the numpy v1.7 function to evaluate polynomials
        was only available in draft form; plus, it won't handle "vectorized"
        calibration coefficients (2D arrays, in which each row is a separate
        set of calibration coeffs).

        The standard convention of storing polynomial coefficients in an array
        is used, namely, the highest degree coefficient is the first element and
        coefficients are stored in descending order.

    Implemented by:

        2014-05-01: Russell Desiderio. Initial Code (no arrays).
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

    # Horner's algorithm
    val = c[:, 5] + x * (c[:, 4] + x * (c[:, 3] + x * (c[:, 2] + x * (c[:, 1] + x * c[:, 0]))))

    return val

# .............................................................................
# TRHPH data products .........................................................
# .............................................................................


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
        tc_slope = thermocouple slope laboratory calibration coefficients
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

        Calculates TRHPHTE-T_TS-AUX, which is an auxiliary data product (not a
        core data product) requested by the authors of the TRHPHTE DPS. It is the
        instrument's thermistor temperature, useful as an important instrument
        diagnostic. It is the same variable as T_ts in the function
        sfl_trhph_vfltemp.

    Implemented by:

        2014-02-28: Russell Desiderio. Initial Code
        2015-01-06: Russell Desiderio. Documented this product as TRHPHTE-T_TS-AUX,
                    following convention established after initial coding.

    Usage:

        T_ts = sfl_trhph_vfl_thermistor_temp(V_ts)

            where

        T_ts = TRHPHTE-T_TS-AUX, thermistor reference temperature [deg_C]
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
        V_R2 = Resistivity voltage 2 (TRHPHR2_L0) [volts]
        V_R3 = Resistivity voltage 3 (TRHPHR3_L0) [volts]
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

    return Cl


# .............................................................................
# PRESF data products .........................................................
# .............................................................................
def sfl_sflpres_rtime(p_psia):
    """
    Description:

        The OOI Level 1 Seafloor Pressure core data products, SFLPRES and
        sub-parameters SFLPRES-RTIME, SFLPRES-TIDE, and SFLPRES-WAVE, are
        created from the Sea-Bird Electronics SBE 26plus member of the Seafloor
        Pressure (SFL, PRESF) family of instruments by either a) polling, in
        real-time, for L0 ASCII text format data output and converting from
        psia to decibar units or b) converting, after instrument recovery, L0
        hexadecimal pressure data into decimal format and the resulting tide
        and wave pressure data in psia to decibar units.

        This code creates the SFLPRES-RTIME data product.

    Implemented by:

        2014-01-31: Craig Risien. Initial Code
        2014-09-23: Christopher Wingard. Minor edits and adds code for -TIDE
                    and -WAVE.
        2015-07-22: Russell Desiderio. There are no type integer input arguments,
                                       don't need replace_fill_with_nan call.

    Usage:

        rtime = sfl_sflpres_rtime(p_psia):

            where

        rtime = real-time pressure (SFLPRES-RTIME_L1) (hydrostatic + atmospheric) [dbar]
        p_psia = pressure (SFLPRES-RTIME_L0) [psia].

    References:

        OOI (2013). Data Product Specification for Seafloor Pressure from
        Sea-Bird SBE 26PLUS. Document Control Number 1341-00230.
        https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
        >> Controlled >> 1000 System Level >>
        1341-00230_Data_Product_SPEC_SFLPRES_OOI.pdf)
    """
    rtime = ne.evaluate('p_psia * 0.689475728')
    return rtime


def sfl_sflpres_tide(p_dec_tide, b, m, slope=1.0, offset=0.0):
    """
    Description:

        The OOI Level 1 Seafloor Pressure core data products, SFLPRES and
        sub-parameters SFLPRES-RTIME, SFLPRES-TIDE, and SFLPRES-WAVE, are
        created from the Sea-Bird Electronics SBE 26plus member of the Seafloor
        Pressure (SFL, PRESF) family of instruments by either a) polling, in
        real-time, for L0 ASCII text format data output and converting from
        psia to decibar units or b) converting, after instrument recovery, L0
        hexadecimal pressure data into decimal format and the resulting tide
        and wave pressure data in psia to decibar units.

        This code creates the SFLPRES-TIDE data product.

    Implemented by:

        2014-09-23: Christopher Wingard. Initial code
        2015-07-22: Russell Desiderio. Added call to replace_fill_with_nan.

    Usage:

        tide = sfl_sflpres_tide(p_dec_tide, b, m, slope, offset):

            where

        tide = tidal pressure (SFLPRES-TIDE_L1) (hydrostatic + atmospheric) [dbar]
        p_dec_tide = tidal pressure (SFLPRES-TIDE_L0) [].
        b = calibration coefficient.
        m = calibration coefficient.
        slope = slope correction factor, 1.0 by default
        offset = offset correction factor, 0.0 by default

    References:

        OOI (2013). Data Product Specification for Seafloor Pressure from
        Sea-Bird SBE 26PLUS. Document Control Number 1341-00230.
        https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
        >> Controlled >> 1000 System Level >>
        1341-00230_Data_Product_SPEC_SFLPRES_OOI.pdf)
    """
    # replace type integer fill values with nans
    p_dec_tide = replace_fill_with_nan(None, p_dec_tide)

    psia = ne.evaluate('slope * ((p_dec_tide - b) / m) + offset')
    tide = ne.evaluate('0.689475728 * psia')
    return tide


def sfl_sflpres_wave(ptcn, p_dec_wave, u0, y1, y2, y3, c1, c2, c3, d1, d2,
                     t1, t2, t3, t4, poff, slope=1.0, offset=0.0):
    """
    Description:

        The OOI Level 1 Seafloor Pressure core data products, SFLPRES and
        sub-parameters SFLPRES-RTIME, SFLPRES-TIDE, and SFLPRES-WAVE, are
        created from the Sea-Bird Electronics SBE 26plus member of the Seafloor
        Pressure (SFL, PRESF) family of instruments by either a) polling, in
        real-time, for L0 ASCII text format data output and converting from
        psia to decibar units or b) converting, after instrument recovery, L0
        hexadecimal pressure data into decimal format and the resulting tide
        and wave pressure data in psia to decibar units.

        This code creates the SFLPRES-WAVE data product.

    Implemented by:

        2014-09-23: Christopher Wingard. Initial code
        2015-07-20: Russell Desiderio. Modified code to accept p_dec_wave as 2D array.
        2015-07-22: Russell Desiderio. Added call to replace_fill_with_nan.

    Usage:

        wave = sfl_sflpres_wave(ptcn, p_dec_wave, u0, y1, y2, y3, c1, c2, c3,
                                d1, d2, t1, t2, t3, t4, slope, offset)

            where

        wave = wave burst pressure (SFLPRES-WAVE_L1) (hydrostatic + atmospheric) [dbar]
        ptcn = pressure temperature compensation number
        p_dec_wave = wave burst pressure (SFLPRES-WAVE_L0) [].
        u0 = calibration coefficient.
        y1 = calibration coefficient.
        y2 = calibration coefficient.
        y3 = calibration coefficient.
        c1 = calibration coefficient.
        c2 = calibration coefficient.
        c3 = calibration coefficient.
        d1 = calibration coefficient.
        d2 = calibration coefficient.
        t1 = calibration coefficient.
        t2 = calibration coefficient.
        t3 = calibration coefficient.
        t4 = calibration coefficient.
        poff = pressure offset calibration coefficient
        slope = slope correction factor, 1.0 by default
        offset = offset correction factor, 0.0 by default

    References:

        OOI (2013). Data Product Specification for Seafloor Pressure from
        Sea-Bird SBE 26PLUS. Document Control Number 1341-00230.
        https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
        >> Controlled >> 1000 System Level >>
        1341-00230_Data_Product_SPEC_SFLPRES_OOI.pdf)
    """
    # replace type integer fill values with nans
    p_dec_wave, ptcn = replace_fill_with_nan(None, p_dec_wave, ptcn)

    # if p_dec_wave is a 1D array make it into a row vector
    p_dec_wave = np.atleast_2d(p_dec_wave)
    n_time_points, n_values_in_burst = p_dec_wave.shape

    # compute the pressure temperature compensation frequency (PTCF) and
    # pressure frequency (PF) from raw inputs
    PTCF = ne.evaluate('ptcn / 256.0')
    PF = ne.evaluate('p_dec_wave / 256.0')

    # use calibration coefficients to compute scale factors.
    U = ne.evaluate('((1.0 / PTCF) * 1000000) - u0')
    C = ne.evaluate('c1 + (c2 * U) + (c3 * U**2)')
    D = ne.evaluate('d1 + d2')
    T0 = ne.evaluate('(t1 + t2 * U + t3 * U**2 + t4 * U**3) / 1000000')
    # broadcast T0 to the shape of PF
    T0 = np.tile(T0, (n_values_in_burst, 1)).T
    W = ne.evaluate('1.0 - (T0**2 * PF**2)')
    # broadcast C, D, and poff to the shape of W
    C = np.tile(C, (n_values_in_burst, 1)).T
    D = np.tile(D, (n_values_in_burst, 1)).T
    poff = np.tile(poff, (n_values_in_burst, 1)).T
    # broadcast slope and offset to the shape of W if not a scalar
    if np.atleast_1d(slope).shape[0] != 1:
        slope = np.tile(slope, (n_values_in_burst, 1)).T
        offset = np.tile(offset, (n_values_in_burst, 1)).T

    # compute the wave pressure data in dbar
    psia = ne.evaluate('slope * ((C * W * (1.0 - D * W)) + poff) + offset')
    wave = ne.evaluate('0.689475728 * psia')
    return wave
