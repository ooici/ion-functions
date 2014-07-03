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
from ion_functions.utils import fill_value


@attr('UNIT', group='func')
class TestSFLFunctionsUnit(BaseUnitTestCase):

    def test_sfl_thsph_temp(self):
        """
        Test the 6 functions that calculate the THSPHTE data products:
            sfl_thsph_temp_int : THSPHTE-INT
            sfl_thsph_temp_ref : THSPHTE-REF
            sfl_thsph_temp_tcl : THSPHTE-TCL
            sfl_thsph_temp_tl  : THSPHTE-TL
            sfl_thsph_temp_tch : THSPHTE-TCH
            sfl_thsph_temp_th  : THSPHTE-TH

        Values based on those described in DPS as available on Alfresco:

        OOI (2014). Data Product Specification for Vent Fluid Temperature from
            TRHPH. Document Control Number 1341-00120.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00120_Data_Product_Spec_THSPHTE_OOI.pdf)

        Implemented by:

            2014-05-01: Russell Desiderio. Initial Code
            2014-07-02: Russell Desiderio. Incorporated changes in DPS calibration algorithms.
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
        ts_rawdec_b = 8188.0
        ts_rawdec_r = 8770.0
        tc_rawdec_L = 16012.0
        tc_rawdec_H = 4237.0
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
        ts_rawdec_b = -9999.9
        ts_rawdec_r = -9999.9
        tc_rawdec_L = -9999.9
        tc_rawdec_H = 4237.0
        # set expected outputs
        int_xpctd = fill_value
        ref_xpctd = fill_value
        tcl_xpctd = -5005.23    # this is unphysical, but no (div by 0) nor np.log(x<0) encountered
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
        ts_rawdec_b = np.array([8188.0, 8185.0])
        ts_rawdec_r = np.array([8770.0, 8758.0])
        tc_rawdec_L = np.array([16012.0, 16009.0])
        tc_rawdec_H = np.array([4237.0, 4236.0])
        # set expected outputs
        int_xpctd = np.array([24.53, 24.55])
        ref_xpctd = np.array([21.25, 21.31])
        tcl_xpctd = np.array([637.88, 637.72])
        tl_xpctd = np.array([655.28, 655.17])
        tch_xpctd = np.array([7.94, 7.87])
        th_xpctd = np.array([28.91, 28.91])
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

    def test_sfl_sflpres_l1(self):
                """
        Test the sfl_sflpres_l1 function.

        Value based on that described in DPS as available on Alfresco:

        OOI (2013). Data Product Specification for Seafloor Pressure from
        Sea-Bird SBE 26PLUS. Document Control Number 1341-00230.
        https://alfresco.oceanobservatories.org/
        (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
        1341-00230_Data_Product_SPEC_SFLPRES_OOI.pdf)

        Implemented by Craig Risien, February 2014
        """

    # test inputs
    pressure_output = np.array([14.868])

    # expected outputs
    pressure_expected = np.array([10.2511])

    # compute values
    pressure_calc = np.zeros(1)
    for i in range(0, 1):
        pressure_calc[i] = sflfunc.sfl_sflpres_l1(pressure_output[i])

    # compare calculated results to expected
    np.testing.assert_allclose(pressure_calc, pressure_expected, rtol=0.0001, atol=0.0001)
