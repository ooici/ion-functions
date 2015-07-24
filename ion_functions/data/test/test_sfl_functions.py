#!/usr/bin/env python

"""
@package ion_functions.test.sfl_functions
@file ion_functions/test/sfl_functions.py
@author Christopher Wingard, Russell Desiderio
@brief Unit tests for sfl_functions module
"""

from nose.plugins.attrib import attr
from ion_functions.test.base_test import BaseUnitTestCase

import numpy as np
from ion_functions.data import sfl_functions as sflfunc

# CI type integer fill value
SYSTEM_FILLVALUE = -999999999
# CI type float fill value has been determined to be np.nan
#    note: initial code and unit tests were written long before the CI convention
#    for handling fill values was established. rather than changing all instances
#    of the variable fill_value in the unit tests to np.nan, this global variable
#    (local to this module) is employed instead.
fill_value = np.nan


@attr('UNIT', group='func')
class TestSFLFunctionsUnit(BaseUnitTestCase):

    def test_sfl_thsph_L2_products(self):
        """
        Test the 6 functions that calculate the THSPH L2 data products:
            sfl_thsph_hydrogen      : THSPHHC
            sfl_thsph_sulfide       : THSPHHS
            sfl_thsph_ph            : THSPHPH-PH
            sfl_thsph_ph_noref      : THSPHPH-PH-NOREF
            sfl_thsph_ph_acl        : THSPHPH-PH-ACL
            sfl_thsph_ph_noref_acl  : THSPHPH-TH-NOREF-ACL

        Test values based on a preliminary version of a spreadsheet supplied with
        preliminary versions of the DPSs: 1341-00190 (PH), 1341-00200 (HS),
        1341-00210 (HC). The L2 test values in the spreadsheets were not derived
        from L0 values, so that a back-calculation had to be made to figure out
        a set of L0 values that would be consistent with the L2 values.

        Implemented by:

            2014-07-09: Russell Desiderio. Initial Code
            2015-07-22: Russell Desiderio. Added test for replacing system fill values
                                           with nans.

        References:

        OOI (2014). Data Product Specification for Vent Fluid pH. Document Control
            Number 1341-00190. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI>> Controlled >> 1000 System Level >>
            1341-00190_Data_Product_Spec_THSPHPH_OOI.pdf)

        OOI (2014). Data Product Specification for Vent Fluid Hydrogen Sulfide
            Concentration. Document Control Number 1341-00200.
            https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI>> Controlled >> 1000 System Level >>
            1341-00200_Data_Product_Spec_THSPHHS_OOI.pdf)

        OOI (2014). Data Product Specification for Vent Fluid Hydrogen Concentration.
            Document Control Number 1341-00210. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI>> Controlled >> 1000 System Level >>
            1341-00210_Data_Product_Spec_THSPHHC_OOI.pdf)

        Prospective location(s) of data file with unit test values. Note that all L2
            data product test data will be contained in one excel file (presumably v5).

        OOI (2014). THSPH L2 data product unit test data. 1341-00190_THSPHPH DPS Artifact.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI >>
            >> REFERENCE >> Data Product Specification Artifacts >> 1341-00190_THSPHPH >>
            THSPHTE Test Data Set v5.xls).   Also see adjoining 1341-00200_THSPHHS and
            1341-00210_THSPHHC entries.
        """
        # for convenience
        sfill = SYSTEM_FILLVALUE

        # calibration polynomial coefficients:

        # electrode engineering to lab calibrated units
        e2l_h2 = np.array([0.0, 0.0, 0.0, 0.0, 1.0, -0.00375])
        e2l_hs = np.array([0.0, 0.0, 0.0, 0.0, 1.0, -0.00350])
        e2l_ysz = np.array([0.0, 0.0, 0.0, 0.0, 1.0, -0.00375])
        e2l_agcl = np.array([0.0, 0.0, 0.0, 0.0, 1.0, -0.00225])
        # electrode material response
        arr_hgo = np.array([0.0, 0.0, 4.38978E-10, -1.88519E-07, -1.88232E-04, 9.23720E-01])
        arr_agcl = np.array([0.0, -8.61134E-10, 9.21187E-07, -3.7455E-04, 6.6550E-02, -4.30086])
        # calculated theoretical reference electrode potential
        arr_agclref = np.array([0.0, 0.0, -2.5E-10, -2.5E-08, -2.5E-06, -9.025E-02])
        # for calculation of chl activity polynomial coefficients
        arr_tac = np.array([0.0, 0.0, -2.80979E-09, 2.21477E-06, -5.53586E-04, 5.723E-02])
        arr_tbc1 = np.array([0.0, 0.0, -6.59572E-08, 4.52831E-05, -1.204E-02, 1.70059])
        arr_tbc2 = np.array([0.0, 0.0, 8.49102E-08, -6.20293E-05, 1.485E-02, -1.41503])
        arr_tbc3 = np.array([-1.86747E-12, 2.32877E-09, -1.18318E-06, 3.04753E-04, -3.956E-02, 2.2047])
        # h2 and h2s fugacity/activity calculations
        arr_logkfh2g = np.array([0.0, 0.0, -1.51904000E-07, 1.16655E-04, -3.435E-02, 6.32102])
        arr_eh2sg = np.array([0.0, 0.0, 0.0, 0.0, -4.49477E-05, -1.228E-02])
        arr_yh2sg = np.array([2.3113E+01, -1.8780E+02, 5.9793E+02, -9.1512E+02, 6.7717E+02, -1.8638E+02])

        # SINGLE-VALUED TESTS

        # hydrogen concentration THSPHHC
        counts_h2 = 4907
        counts_ysz = 7807
        temperature = 320.0
        h2_xpctd = 0.02712

        h2_calc = sflfunc.sfl_thsph_hydrogen(counts_h2, counts_ysz, temperature, e2l_h2, e2l_ysz, arr_hgo,
                                             arr_logkfh2g)
        np.testing.assert_allclose(h2_calc, h2_xpctd, rtol=0.0, atol=0.0001)

        # sulfide concentration THSPHHS
        counts_hs = 3806
        counts_ysz = 7007
        temperature = 320.0
        h2s_xpctd = 0.95744

        h2s_calc = sflfunc.sfl_thsph_sulfide(counts_hs, counts_ysz, temperature, e2l_hs, e2l_ysz, arr_hgo,
                                             arr_logkfh2g, arr_eh2sg, arr_yh2sg)
        np.testing.assert_allclose(h2s_calc, h2s_xpctd, rtol=0.0, atol=0.0001)

        # pH: THSPHPH-PH  and  THSPHPH-PH-ACL
        counts_agcl = 7801
        counts_ysz = 6607
        temperature = 300.0
        trhphcc = 400.0     # [chloride] from trhphcc data product has units of mmol/kg
        # THSPHPH-PH
        pH_xpctd = 5.3857
        pH_calc = sflfunc.sfl_thsph_ph(counts_ysz, counts_agcl, temperature, e2l_ysz, e2l_agcl, arr_hgo,
                                       arr_agcl, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3, trhphcc)
        np.testing.assert_allclose(pH_calc, pH_xpctd, rtol=0.0, atol=0.001)
        # THSPHPH-PH-ACL (no chloride measurement supplied)
        pH_acl_xpctd = 5.2357
        pH_calc = sflfunc.sfl_thsph_ph_acl(counts_ysz, counts_agcl, temperature, e2l_ysz, e2l_agcl,
                                           arr_hgo, arr_agcl, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3)
        np.testing.assert_allclose(pH_calc, pH_acl_xpctd, rtol=0.0, atol=0.001)

        # pH: THSPHPH-PH-NOREF  and  THSPHPH-PH-NOREF-ACL
        counts_ysz = 6207
        temperature = 300.0
        # THSPHPH-PH-NOREF
        pH_noref_xpctd = 4.5064
        pH_calc = sflfunc.sfl_thsph_ph_noref(counts_ysz, temperature, arr_agclref, e2l_ysz, arr_hgo,
                                             arr_agcl, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3, trhphcc)
        np.testing.assert_allclose(pH_calc, pH_noref_xpctd, rtol=0.0, atol=0.001)
        # THSPHPH-PH-NOREF-ACL (no chloride measurement supplied)
        pH_noref_acl_xpctd = 4.3564
        pH_calc = sflfunc.sfl_thsph_ph_noref_acl(counts_ysz, temperature, arr_agclref, e2l_ysz, arr_hgo,
                                                 arr_agcl, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3)
        np.testing.assert_allclose(pH_calc, pH_noref_acl_xpctd, rtol=0.0, atol=0.001)

        # VECTORIZED TESTS

        # set up cal coeff arrays.
        # since i will use 7 sets of L0 values for the pH function inputs,
        # i'll set up all coeff arrays to be 7x6, and also run 7 values
        # for the hydrogen and sulfide functions.
        npackets = 7
        e2l_h2 = np.tile(e2l_h2, (npackets, 1))
        e2l_hs = np.tile(e2l_hs, (npackets, 1))
        e2l_ysz = np.tile(e2l_ysz, (npackets, 1))
        e2l_agcl = np.tile(e2l_agcl, (npackets, 1))
        # electrode material response
        arr_hgo = np.tile(arr_hgo, (npackets, 1))
        arr_agcl = np.tile(arr_agcl, (npackets, 1))
        # calculated theoretical reference electrode potential
        arr_agclref = np.tile(arr_agclref, (npackets, 1))
        # for calculation of chl activity polynomial coefficients
        arr_tac = np.tile(arr_tac, (npackets, 1))
        arr_tbc1 = np.tile(arr_tbc1, (npackets, 1))
        arr_tbc2 = np.tile(arr_tbc2, (npackets, 1))
        arr_tbc3 = np.tile(arr_tbc3, (npackets, 1))
        # h2 and h2s fugacity/activity calculations
        arr_logkfh2g = np.tile(arr_logkfh2g, (npackets, 1))
        arr_eh2sg = np.tile(arr_eh2sg, (npackets, 1))
        arr_yh2sg = np.tile(arr_yh2sg, (npackets, 1))

        # hydrogen concentration THSPHHC
        counts_h2 = np.array([4907, 4207, 4907, 4207, sfill, 4207, 4907])
        counts_ysz = np.array([7807, 7407, 7807, 7407, 7807, sfill, 7807])
        temperature = np.array([320.0, 250.0, 320.0, 250.0, 320.0, 250.0, 320.0])
        h2_xpctd = np.array([0.02712, 0.09265, 0.02712, 0.09265, fill_value, fill_value, 0.02712])

        h2_calc = sflfunc.sfl_thsph_hydrogen(counts_h2, counts_ysz, temperature, e2l_h2, e2l_ysz, arr_hgo,
                                             arr_logkfh2g)
        np.testing.assert_allclose(h2_calc, h2_xpctd, rtol=0.0, atol=0.0001)

        # sulfide concentration THSPHHS
        counts_hs = np.array([3806, sfill, 3806, 3166, 3806, 3166, 3806])
        counts_ysz = np.array([7007, 6607, 7007, 6607, sfill, 6607, 7007])
        temperature = np.array([320.0, 260.0, 320.0, 260.0, 320.0, 260.0, 320.0])
        h2s_xpctd = np.array([0.95744, fill_value, 0.95744, 3.70778, fill_value, 3.70778, 0.95744])

        h2s_calc = sflfunc.sfl_thsph_sulfide(counts_hs, counts_ysz, temperature, e2l_hs, e2l_ysz, arr_hgo,
                                             arr_logkfh2g, arr_eh2sg, arr_yh2sg)
        np.testing.assert_allclose(h2s_calc, h2s_xpctd, rtol=0.0, atol=0.0001)

        # pH: THSPHPH-PH  and  THSPHPH-PH-ACL
        counts_agcl = np.array([7801, 7001, sfill, 7801, 8201, 7801, 7801])
        counts_ysz = np.array([7407, 6207, 6607, 6207, 6207, 5407, 8207])
        temperature = np.tile(300, npackets)
        trhphcc = np.tile(400.0, npackets)  # [chloride] from trhphcc data product has units of mmol/kg
        # THSPHPH-PH
        pH_xpctd = np.array([fill_value, 6.26505, fill_value, 4.50641, 3.62709, fill_value, fill_value])
        pH_calc = sflfunc.sfl_thsph_ph(counts_ysz, counts_agcl, temperature, e2l_ysz, e2l_agcl, arr_hgo,
                                       arr_agcl, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3, trhphcc)
        np.testing.assert_allclose(pH_calc, pH_xpctd, rtol=0.0, atol=0.001)
        # THSPHPH-PH-ACL (no chloride measurement supplied)
        pH_acl_xpctd = np.array([6.99431, 6.11499, fill_value, 4.35635, 3.47704, fill_value, fill_value])
        pH_calc = sflfunc.sfl_thsph_ph_acl(counts_ysz, counts_agcl, temperature, e2l_ysz, e2l_agcl,
                                           arr_hgo, arr_agcl, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3)
        np.testing.assert_allclose(pH_calc, pH_acl_xpctd, rtol=0.0, atol=0.001)

        # pH: THSPHPH-PH-NOREF  and  THSPHPH-PH-NOREF-ACL
        counts_ysz = np.array([7807, 7407, 7007, 6607, 6207, 5807, 5407])
        temperature = np.tile(300, npackets)
        # THSPHPH-PH-NOREF
        pH_noref_xpctd = np.array([fill_value, fill_value, 6.26505, 5.38573, 4.50641, 3.62709, fill_value])
        pH_calc = sflfunc.sfl_thsph_ph_noref(counts_ysz, temperature, arr_agclref, e2l_ysz, arr_hgo,
                                             arr_agcl, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3, trhphcc)
        np.testing.assert_allclose(pH_calc, pH_noref_xpctd, rtol=0.0, atol=0.001)
        # THSPHPH-PH-NOREF-ACL (no chloride measurement supplied)
        pH_noref_acl_xpctd = np.array([fill_value, 6.99431, 6.11499, 5.23567, 4.35635, 3.47704, fill_value])
        pH_calc = sflfunc.sfl_thsph_ph_noref_acl(counts_ysz, temperature, arr_agclref, e2l_ysz, arr_hgo,
                                                 arr_agcl, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3)
        np.testing.assert_allclose(pH_calc, pH_noref_acl_xpctd, rtol=0.0, atol=0.001)

    def test_sfl_thsph_temp(self):
        """
        Test the 6 functions that calculate the THSPHTE L1 data products:
            sfl_thsph_temp_int : THSPHTE-INT
            sfl_thsph_temp_ref : THSPHTE-REF
            sfl_thsph_temp_tcl : THSPHTE-TCL
            sfl_thsph_temp_tl  : THSPHTE-TL
            sfl_thsph_temp_tch : THSPHTE-TCH
            sfl_thsph_temp_th  : THSPHTE-TH

        Values based on those described in DPS and DPS artifact as available on Alfresco:

        OOI (2014). Data Product Specification for Vent Fluid Temperature from
            TRHPH. Document Control Number 1341-00120.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00120_Data_Product_Spec_THSPHTE_OOI.pdf)

        OOI (2014). THSPHTE unit test data. 1341-00120_THSPHTE DPS Artifact.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI >>
            >> REFERENCE >> Data Product Specification Artifacts >> 1341-00120_THSPHTE >>
            THSPHTE Test Data Set v4.xls)

        Implemented by:

            2014-05-01: Russell Desiderio. Initial Code
            2014-07-02: Russell Desiderio. Incorporated changes in DPS calibration algorithms.
            2015-07-24: Russell Desiderio. Incorporated SYSTEM_FILLVALUE checks.
        """
        # calibration constants: b thermistor
        # engineering values to lab calibrated values
        e2l_b = np.array([0.0, 0.0, 0.0, 0.0, 1.04938, -275.5])
        # lab calibrated values to scientific values
        l2s_b = np.array([0.0, 0.0, 8.7755e-08, 0.0, 0.000234101, 0.001129306])

        # calibration constants: r thermistor
        # engineering values to lab calibrated values
        e2l_r = np.array([0.0, 0.0, 0.0, 0.0, 1.04938, -275.5])
        # lab calibrated values to scientific values
        l2s_r = np.array([0.0, 0.0, 8.7755e-08, 0.0, 0.000234101, 0.001129306])

        # calibration constants: L thermocouple
        # engineering values to lab calibrated values
        e2l_L = np.array([0.0, 0.0, 0.0, 0.0, 0.9964, -0.46112])
        # lab calibrated values to scientific values
        l2s_L = np.array([9.32483e-7, -0.000122268, 0.00702, -0.23532, 17.06172, 0.0])

        # calibration constants: H thermocouple
        # engineering values to lab calibrated values
        e2l_H = np.array([0.0, 0.0, 0.0, 0.0, 0.9979, -0.10287])
        # lab calibrated values to scientific values
        l2s_H = np.array([9.32483e-7, -0.000122268, 0.00702, -0.23532, 17.06172, 0.0])

        # calibration constants: convert 'r' thermistor scientific (temperature)
        # values to thermocouple equivalent voltage [mV]
        s2v_r = np.array([5.83124e-14, -4.09038e-11, -3.44498e-8, 5.14528e-5, 0.05841, 0.00209])

        # rawdata: aH200B200720C420A1108D3E8C22421FFC#
        # test inputs
        ts_rawdec_b = 8188
        ts_rawdec_r = 8770
        tc_rawdec_L = 16012
        tc_rawdec_H = 4237
        # set expected outputs
        int_xpctd = 24.53
        ref_xpctd = 21.25
        tcl_xpctd = 637.88
        tl_xpctd = 655.28
        tch_xpctd = 7.94
        th_xpctd = 28.91

        ## rawdata: aH2009200820C220A0108C3E8922361FF9#
        ## test inputs
        #ts_rawdec_b = 8185.0
        #ts_rawdec_r = 8758.0
        #tc_rawdec_L = 16009.0
        #tc_rawdec_H = 4236.0
        ## set expected outputs
        #int_xpctd = 24.55
        #ref_xpctd = 21.31
        #tcl_xpctd = 637.72
        #tl_xpctd = 655.17
        #tch_xpctd = 7.87
        #th_xpctd = 28.91

        # SINGLE VALUED TESTS

        # calculate the various temperature dataproducts and
        # compare to the expected outputs; lab accuracy is
        # between 0.2 and 1.0 degC, so check calculation to 0.01 degC.
        int_calc = sflfunc.sfl_thsph_temp_int(ts_rawdec_b, e2l_b, l2s_b)
        np.testing.assert_allclose(int_calc, int_xpctd, rtol=0.0, atol=0.01)

        ref_calc = sflfunc.sfl_thsph_temp_ref(ts_rawdec_r, e2l_r, l2s_r)
        np.testing.assert_allclose(ref_calc, ref_xpctd, rtol=0.0, atol=0.01)

        tcl_calc = sflfunc.sfl_thsph_temp_tcl(tc_rawdec_L, e2l_L, l2s_L)
        np.testing.assert_allclose(tcl_calc, tcl_xpctd, rtol=0.0, atol=0.01)

        tch_calc = sflfunc.sfl_thsph_temp_tch(tc_rawdec_H, e2l_H, l2s_H)
        np.testing.assert_allclose(tch_calc, tch_xpctd, rtol=0.0, atol=0.01)

        tl_calc = sflfunc.sfl_thsph_temp_tl(tc_rawdec_L, e2l_L, l2s_L,
                                            ts_rawdec_r, e2l_r, l2s_r, s2v_r)
        np.testing.assert_allclose(tl_calc, tl_xpctd, rtol=0.0, atol=0.01)

        th_calc = sflfunc.sfl_thsph_temp_th(tc_rawdec_H, e2l_H, l2s_H,
                                            ts_rawdec_r, e2l_r, l2s_r, s2v_r)
        np.testing.assert_allclose(th_calc, th_xpctd, rtol=0.0, atol=0.01)

        # trap-out-nans-and-inf test -> replace with fill_value
        # test inputs
        ts_rawdec_b = SYSTEM_FILLVALUE
        ts_rawdec_r = SYSTEM_FILLVALUE
        tc_rawdec_L = SYSTEM_FILLVALUE
        tc_rawdec_H = 4237.0
        # set expected outputs
        int_xpctd = fill_value
        ref_xpctd = fill_value
        tcl_xpctd = fill_value
        tl_xpctd = fill_value   # calc uses ref_xpctd
        tch_xpctd = 7.94
        th_xpctd = fill_value   # calc uses ref_xpctd

        int_calc = sflfunc.sfl_thsph_temp_int(ts_rawdec_b, e2l_b, l2s_b)
        np.testing.assert_allclose(int_calc, int_xpctd, rtol=0.0, atol=0.01)

        ref_calc = sflfunc.sfl_thsph_temp_ref(ts_rawdec_r, e2l_r, l2s_r)
        np.testing.assert_allclose(ref_calc, ref_xpctd, rtol=0.0, atol=0.01)

        tcl_calc = sflfunc.sfl_thsph_temp_tcl(tc_rawdec_L, e2l_L, l2s_L)
        np.testing.assert_allclose(tcl_calc, tcl_xpctd, rtol=0.0, atol=0.01)

        tch_calc = sflfunc.sfl_thsph_temp_tch(tc_rawdec_H, e2l_H, l2s_H)
        np.testing.assert_allclose(tch_calc, tch_xpctd, rtol=0.0, atol=0.01)

        tl_calc = sflfunc.sfl_thsph_temp_tl(tc_rawdec_L, e2l_L, l2s_L,
                                            ts_rawdec_r, e2l_r, l2s_r, s2v_r)
        np.testing.assert_allclose(tl_calc, tl_xpctd, rtol=0.0, atol=0.01)

        th_calc = sflfunc.sfl_thsph_temp_th(tc_rawdec_H, e2l_H, l2s_H,
                                            ts_rawdec_r, e2l_r, l2s_r, s2v_r)
        np.testing.assert_allclose(th_calc, th_xpctd, rtol=0.0, atol=0.01)

        # DOUBLE VALUED (VECTORIZED) TESTS
        # test inputs
        ts_rawdec_b = np.array([8188, 8185])
        ts_rawdec_r = np.array([8770, 8758])
        tc_rawdec_L = np.array([16012, 16009])
        tc_rawdec_H = np.array([4237, SYSTEM_FILLVALUE])
        # set expected outputs
        int_xpctd = np.array([24.53, 24.55])
        ref_xpctd = np.array([21.25, 21.31])
        tcl_xpctd = np.array([637.88, 637.72])
        tl_xpctd = np.array([655.28, 655.17])
        tch_xpctd = np.array([7.94, fill_value])
        th_xpctd = np.array([28.91, fill_value])
        # tile the cal coeff arrays
        e2l_b = np.tile(e2l_b, (2, 1))
        l2s_b = np.tile(l2s_b, (2, 1))
        e2l_r = np.tile(e2l_r, (2, 1))
        l2s_r = np.tile(l2s_r, (2, 1))
        e2l_L = np.tile(e2l_L, (2, 1))
        l2s_L = np.tile(l2s_L, (2, 1))
        e2l_H = np.tile(e2l_H, (2, 1))
        l2s_H = np.tile(l2s_H, (2, 1))
        s2v_r = np.tile(s2v_r, (2, 1))

        int_calc = sflfunc.sfl_thsph_temp_int(ts_rawdec_b, e2l_b, l2s_b)
        np.testing.assert_allclose(int_calc, int_xpctd, rtol=0.0, atol=0.01)

        ref_calc = sflfunc.sfl_thsph_temp_ref(ts_rawdec_r, e2l_r, l2s_r)
        np.testing.assert_allclose(ref_calc, ref_xpctd, rtol=0.0, atol=0.01)

        tcl_calc = sflfunc.sfl_thsph_temp_tcl(tc_rawdec_L, e2l_L, l2s_L)
        np.testing.assert_allclose(tcl_calc, tcl_xpctd, rtol=0.0, atol=0.01)

        tch_calc = sflfunc.sfl_thsph_temp_tch(tc_rawdec_H, e2l_H, l2s_H)
        np.testing.assert_allclose(tch_calc, tch_xpctd, rtol=0.0, atol=0.01)

        tl_calc = sflfunc.sfl_thsph_temp_tl(tc_rawdec_L, e2l_L, l2s_L,
                                            ts_rawdec_r, e2l_r, l2s_r, s2v_r)
        np.testing.assert_allclose(tl_calc, tl_xpctd, rtol=0.0, atol=0.01)

        th_calc = sflfunc.sfl_thsph_temp_th(tc_rawdec_H, e2l_H, l2s_H,
                                            ts_rawdec_r, e2l_r, l2s_r, s2v_r)
        np.testing.assert_allclose(th_calc, th_xpctd, rtol=0.0, atol=0.01)

    def test_sfl_trhph_vfltemp(self):
        """
        Test sfl_trhph_vfltemp function.

        Values based on those described in DPS as available on Alfresco:

        OOI (2012). Data Product Specification for Vent Fluid Temperature from
            TRHPH. Document Control Number 1341-00150.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00150_Data_Product_Spec_TRHPHTE_OOI.pdf)

        Implemented by Christopher Wingard, April 2013
        """

        # calibration constants
        tc_slope = np.ones(24) * 4.22e-5
        ts_slope = np.ones(24) * 0.003

        #   V_ts, V_tc, T_ts, T
        test_array = np.array([
            [1.506,	0.000,	12.01,	12.0],
            [1.479,	0.015,	12.67,	17.1],
            [1.926,	0.001,	2.47,	2.1],
            [1.932,	0.274,	2.34,	69.5],
            [1.927,	0.306,	2.45,	77.5],
            [1.930,	0.393,	2.38,	99.0],
            [1.929,	0.383,	2.40,	96.6],
            [1.930,	0.388,	2.38,	97.8],
            [1.930,	0.469,	2.38,	117.8],
            [1.931,	1.077,	2.36,	268.3],
            [1.926,	1.288,	2.47,	320.6],
            [1.926,	1.305,	2.47,	324.8],
            [1.928,	1.319,	2.43,	328.2],
            [1.929,	1.318,	2.40,	328.0],
            [1.601,	0.091,	9.75,	29.7],
            [1.958,	0.407,	1.78,	102.0],
            [1.962,	1.216,	1.70,	302.2],
            [1.441,	-0.300,	13.61,	13.6],
            [1.441,	1.216,  13.61,	312.1],
            [1.321,	1.216,	16.68,	315.0],
            [1.110,	1.216,	22.60,	320.8],
            [1.600,	1.216,	9.77,	308.1],
            [1.590,	1.216,	10.01,	308.9],
            [1.591,	1.216,	9.98,	308.3],
        ])

        # set inputs
        V_ts = test_array[:, 0]
        V_tc = test_array[:, 1]

        # set known outputs
        T = test_array[:, 3]

        # calculate the Vent Fluid Temperature
        Tout = sflfunc.sfl_trhph_vfltemp(V_ts, V_tc, tc_slope, ts_slope)

        # How'd we do?
        np.testing.assert_allclose(Tout, T, rtol=0.1, atol=0.0)

    def test_sfl_trhph_vfl_thermistor_temp(self):
        """
        Test sfl_trhph_vfl_thermistor_temp function.

        Values based on those described in DPS as available on Alfresco:

        OOI (2012). Data Product Specification for Vent Fluid Temperature from
            TRHPH. Document Control Number 1341-00150.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00150_Data_Product_Spec_TRHPHTE_OOI.pdf)

        Implemented by Russell Desiderio, 28-Feb-2014
        """

        #   V_ts, V_tc, T_ts, T
        test_array = np.array([
            [1.506,	0.000,	12.01,	12.0],
            [1.479,	0.015,	12.67,	17.1],
            [1.926,	0.001,	2.47,	2.1],
            [1.932,	0.274,	2.34,	69.5],
            [1.927,	0.306,	2.45,	77.5],
            [1.930,	0.393,	2.38,	99.0],
            [1.929,	0.383,	2.40,	96.6],
            [1.930,	0.388,	2.38,	97.8],
            [1.930,	0.469,	2.38,	117.8],
            [1.931,	1.077,	2.36,	268.3],
            [1.926,	1.288,	2.47,	320.6],
            [1.926,	1.305,	2.47,	324.8],
            [1.928,	1.319,	2.43,	328.2],
            [1.929,	1.318,	2.40,	328.0],
            [1.601,	0.091,	9.75,	29.7],
            [1.958,	0.407,	1.78,	102.0],
            [1.962,	1.216,	1.70,	302.2],
            [1.441,	-0.300,	13.61,	13.6],
            [1.441,	1.216,  13.61,	312.1],
            [1.321,	1.216,	16.68,	315.0],
            [1.110,	1.216,	22.60,	320.8],
            [1.600,	1.216,	9.77,	308.1],
            [1.590,	1.216,	10.01,	308.9],
            [1.591,	1.216,	9.98,	308.3],
        ])

        # set inputs
        V_ts = test_array[:, 0]

        # set known outputs
        T_ts = test_array[:, 2]

        # calculate the thermistor temperature
        Tout = sflfunc.sfl_trhph_vfl_thermistor_temp(V_ts)

        # How'd we do?
        np.testing.assert_allclose(Tout, T_ts, rtol=0.01, atol=0.0)

    def test_sfl_trhph_vflorp(self):
        """
        Test sfl_trhph_vflorp function.

        Values based on those described in DPS as available on Alfresco:

        OOI (2012). Data Product Specification for Vent Fluid Oxidation-
            Reduction Potential (ORP). Document Control Number 1341-00170.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00170_Data_Product_Spec_TRHPHEH_OOI.pdf)

        Implemented by Russell Desiderio, 28-Feb-2014
        """
        # the test data from the DPS do not compute correctly.
        # to make it work, the 3rd test result was changed from -49 to -48,
        # and the offset was changed from 2008.0 to 2004.0

        # test calibration constants
        offset = 2004.0
        gain = 4.0

        #   [  V      ORP[mV] ]
        test_array = np.array([
            [1.806,   -50.],
            [1.541,  -116.],
            [1.810,   -48.],
            [0.735,  -317.],
            [0.745,  -315.],
            [0.715,  -322.],
            [0.775,  -307.],
            [0.799,  -301.],
            [0.757,  -312.],
            [0.542,  -366.],
            [0.831,  -293.],
            [0.867,  -284.],
            [0.911,  -273.]
        ])

        # set inputs
        V = test_array[:, 0]

        # set known outputs
        ORP = test_array[:, -1]

        # calculate the oxidation-reduction potential
        ORP_out = sflfunc.sfl_trhph_vflorp(V, offset, gain)

        # How'd we do?
        np.testing.assert_allclose(ORP_out, ORP, rtol=0.001, atol=0.0)

    def test_sfl_trhph_chloride(self):
        """
        Test sfl_trhph_chloride function.

        Values based on those described in DPS as available on Alfresco:

        OOI (2012). Data Product Specification for Vent Fluid Chloride
            Concentration. Document Control Number 1341-00160.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00160_Data_Product_Spec_TRHPHCC_OOI.pdf)

        Implemented by Christopher Wingard, April 2013
        Modified by Russell Desiderio, Feb. 28, 2014.
            Extended unit tests; see comments below.
        """
        # original unit test data from DPS:
        #    does not test all v_r1 v_r2 v_r3 branches.
        #    does not unambiguously test resistivity value out of range.
        #        (all of the test result fill_values could arise
        #         solely from temperature out-of-range)
        #test_array = np.array([
        #    [0.906, 4.095, 4.095, 11.8, 4.530, 0.2208, fill_value, fill_value],
        #    [0.890, 4.095, 4.095, 15.9, 4.450, 0.2247, fill_value, fill_value],
        #    [0.891, 4.095, 4.095, 3.2, 4.455, 0.2245, fill_value, fill_value],
        #    [0.184, 0.915, 4.064, 67.7, 0.915, 1.0929, fill_value, fill_value],
        #    [0.198, 1.002, 4.095, 75.8, 1.002, 0.9980, fill_value, fill_value],
        #    [0.172, 0.857, 4.082, 97.5, 0.857, 1.1669, fill_value, fill_value],
        #    [0.183, 0.926, 4.076, 95.0, 0.926, 1.0799, fill_value, fill_value],
        #    [0.233, 1.182, 4.072, 96.2, 1.182, 0.8460, fill_value, fill_value],
        #    [0.146, 0.747, 3.634, 116.8, 0.727, 1.3759, 0.19507, 195],
        #    [0.134, 0.681, 3.405, 272.8, 0.681, 1.4684, 0.10893, 109],
        #    [0.131, 0.673, 3.293, 325.8, 0.659, 1.5184, 0.12813, 128],
        #    [0.133, 0.678, 3.396, 330.0, 0.679, 1.4723, 0.12664, 127],
        #    [0.135, 0.681, 3.409, 333.4, 0.682, 1.4667, 0.12878, 129],
        #    [0.135, 0.681, 3.426, 333.2, 0.685, 1.4594, 0.12802, 128]
        #])

        # the first fill_value tests "out-of-range" resistivity value v_r1 which
        #     results in a conductivity value that, with the corresponding
        #     temperature, does not result in a [chl] value because of the
        #     curvature of the empirical [chl] surface as a function of (Temp,Cond).
        # the second fill_value tests an out-of-range temperature value (84.6).
        # the third fill_value tests a resistivity value (v_r2=2.8) which
        #     results in an out-of-range conductivity value.
        #
        #     v_r1    v_r2    v_r3   temp    chl [mmol/kg]
        test_array = np.array([
            [0.440,  4.095,  4.095,  105.4,    59.0],
            [0.380,  4.095,  4.095,  241.9,   fill_value],
            [0.320,  4.095,  4.095,  374.2,    60.0],
            [0.184,  0.915,  4.064,  105.4,   175.0],
            [0.198,  1.002,  4.095,  241.9,    71.0],
            [0.172,  0.857,  4.082,  374.2,   132.0],
            [0.183,  0.926,  4.076,   84.6,   fill_value],
            [0.233,  2.800,  4.072,  250.2,   fill_value],
            [0.146,  0.747,  3.634,  116.8,   195.0],
            [0.134,  0.681,  3.405,  272.8,   109.0],
            [0.131,  0.673,  3.293,  325.8,   128.0],
            [0.133,  0.678,  3.396,  330.0,   127.0],
            [0.135,  0.681,  2.000,  333.4,   239.0],
            [0.135,  0.681,  1.000,  333.2,   501.0]
        ])

        # set inputs
        V_R1 = test_array[:, 0]
        V_R2 = test_array[:, 1]
        V_R3 = test_array[:, 2]
        T = test_array[:, 3]

        # set known output
        Cl = test_array[:, -1]

        # calculate the chloride concentration
        Clout = sflfunc.sfl_trhph_chloride(V_R1, V_R2, V_R3, T)

        ###########################################################################
        # How'd we do?
        np.testing.assert_allclose(Clout, Cl, rtol=0, atol=0)
        ###########################################################################

    def test_sfl_sflpres_rtime(self):
                """
        Test the sfl_sflpres_rtime function.

        Value based on that described in DPS as available on Alfresco:

        OOI (2013). Data Product Specification for Seafloor Pressure from
        Sea-Bird SBE 26PLUS. Document Control Number 1341-00230.
        https://alfresco.oceanobservatories.org/
        (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
        1341-00230_Data_Product_SPEC_SFLPRES_OOI.pdf)

        Implemented by Craig Risien, February 2014
        Added time-vectorized unit tests, Russell Desiderio, 2015-07-20
        """

    ### scalar test
    # test inputs
    pressure_output = np.array([14.868])
    # expected outputs
    pressure_expected = np.array([10.2511])
    # compute value
    pressure_calc = sflfunc.sfl_sflpres_rtime(pressure_output)
    # compare calculated results to expected
    np.testing.assert_allclose(pressure_calc, pressure_expected, rtol=0.0001, atol=0.0001)

    ### time-vectorized test
    # test inputs
    pressure_output = np.array([14.868, 14.868])
    # expected outputs
    pressure_expected = np.array([10.2511, 10.2511])
    # compute value
    pressure_calc = sflfunc.sfl_sflpres_rtime(pressure_output)
    # compare calculated results to expected
    np.testing.assert_allclose(pressure_calc, pressure_expected, rtol=0.0001, atol=0.0001)

    def test_sfl_sflpres_tide(self):
                """
        Test the sfl_sflpres_tide function.

        Value based on that described in DPS as available on Alfresco:

        OOI (2013). Data Product Specification for Seafloor Pressure from
        Sea-Bird SBE 26PLUS. Document Control Number 1341-00230.
        https://alfresco.oceanobservatories.org/
        (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
        1341-00230_Data_Product_SPEC_SFLPRES_OOI.pdf)

        Implemented by Craig Risien, February 2014
        Added time-vectorized unit tests, Russell Desiderio, 2015-07-20,22
        """

    ### scalar test
    # test inputs
    p_dec_tide = 4175754
    m = 279620.2
    b = 18641.3
    # expected output
    tide = 10.2504
    # compute values
    out = sflfunc.sfl_sflpres_tide(p_dec_tide, b, m)
    # compare calculated results to expected
    np.testing.assert_allclose(out, tide, rtol=0.0001, atol=0.0001)

    ### time-vectorized test with system fill value
    # test inputs
    p_dec_tide = np.array([4175754, SYSTEM_FILLVALUE, 4175754])
    m = np.array([279620.2, 279620.2, 279620.2])
    b = np.array([18641.3, 18641.3, 18641.3])
    # expected output
    tide = np.array([10.2504, np.nan, 10.2504])
    # compute values
    out = sflfunc.sfl_sflpres_tide(p_dec_tide, b, m)
    # compare calculated results to expected
    np.testing.assert_allclose(out, tide, rtol=0.0001, atol=0.0001)

    ### time-vectorized test with system fill value,
    ### with time-vectorized slope and offset input
    # test inputs
    p_dec_tide = np.array([4175754, SYSTEM_FILLVALUE, 4175754])
    m = np.array([279620.2, 279620.2, 279620.2])
    b = np.array([18641.3, 18641.3, 18641.3])
    slope = np.array([1.0, 1.0, 1.0])
    offset = np.array([0.0, 0.0, 0.0])
    # expected output
    tide = np.array([10.2504, np.nan, 10.2504])
    # compute values
    out = sflfunc.sfl_sflpres_tide(p_dec_tide, b, m, slope, offset)
    # compare calculated results to expected
    np.testing.assert_allclose(out, tide, rtol=0.0001, atol=0.0001)

    def test_sfl_sflpres_wave(self):
                """
        Test the sfl_sflpres_wave function.

        Value based on that described in DPS as available on Alfresco:

        OOI (2013). Data Product Specification for Seafloor Pressure from
        Sea-Bird SBE 26PLUS. Document Control Number 1341-00230.
        https://alfresco.oceanobservatories.org/
        (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
        1341-00230_Data_Product_SPEC_SFLPRES_OOI.pdf)

        Implemented by Craig Risien, February 2014
        Russell Desiderio, 2015-07-20,22: Corrected p_dec_wave input to be a 1D array for
                                       the scalar time case. Added time-vectorized case.
        """

    ### scalar time case
    # test inputs
    p_dec_wave = np.array([8900312, SYSTEM_FILLVALUE, 8900312])
    ptcn = 43746280
    u0 = 5.856409e+00
    y1 = -3.987838e+03
    y2 = -1.049603e+04
    y3 = 0.000000e+00
    c1 = 2.305367e+02
    c2 = 1.198422e+01
    c3 = -2.401512e+02
    d1 = 4.095400e-02
    d2 = 0.000000e+00
    t1 = 2.781994e+01
    t2 = 6.760780e-01
    t3 = 1.761829e+01
    t4 = 6.000932e+00
    poff = 0.0
    slope = 1.0
    offset = 0.0

    # expected outputs - CI requires a 2D row vector
    wave = np.array([[10.2511, np.nan, 10.2511]])

    # compute values
    out = sflfunc.sfl_sflpres_wave(ptcn, p_dec_wave, u0, y1, y2, y3, c1,
                                   c2, c3, d1, d2, t1, t2, t3, t4, poff, slope, offset)

    # compare calculated results to expected
    np.testing.assert_allclose(out, wave, rtol=0.0001, atol=0.0001)

    ### time-vectorized
    npts = 5
    args = [ptcn, u0, y1, y2, y3, c1, c2, c3, d1, d2, t1, t2, t3, t4, poff, slope, offset]
    args = [np.tile(x, npts) for x in args]
    [ptcn, u0, y1, y2, y3, c1, c2, c3, d1, d2, t1, t2, t3, t4, poff, slope, offset] = args

    p_dec_wave = np.tile(p_dec_wave, (npts, 1))
    wave = np.tile(wave, (npts, 1))

    # change the last ptcn value to a system fill value,
    ptcn[-1] = SYSTEM_FILLVALUE
    # which should change the last row of the expected wave output to nans.
    wave[-1, :] = np.nan

    # compute values
    out = sflfunc.sfl_sflpres_wave(ptcn, p_dec_wave, u0, y1, y2, y3, c1,
                                   c2, c3, d1, d2, t1, t2, t3, t4, poff, slope, offset)

    # compare calculated results to expected
    np.testing.assert_allclose(out, wave, rtol=0.0001, atol=0.0001)
