#!/usr/bin/env python
from ion_functions.data.perf.test_performance import PerformanceTestCase, a_deca
from ion_functions.data.sfl_functions import (sfl_trhph_vfltemp,
                                              sfl_trhph_vflorp,
                                              sfl_trhph_chloride,
                                              sfl_sflpres_l1,
                                              sfl_thsph_temp_int,
                                              sfl_thsph_temp_ref,
                                              sfl_thsph_temp_tcl,
                                              sfl_thsph_temp_tch,
                                              sfl_thsph_temp_tl,
                                              sfl_thsph_temp_th)
import numpy as np


class TestSFLPerformance(PerformanceTestCase):

    def test_sfl_thsph_temp_int(self):
        stats = []
        # create 10000 data packets
        # input, raw decimal counts
        ts_rawdec_b = np.zeros((a_deca, 1)) + 8188.0
        # cal coeffs: engineering values to lab calibrated values
        c0_e2l_b = np.zeros((a_deca, 1)) + 0.05935
        c1_e2l_b = np.zeros((a_deca, 1)) + 0.00099151
        c2_e2l_b = np.zeros((a_deca, 1)) + 3.82028e-10
        c3_e2l_b = np.zeros((a_deca, 1)) + 4.54486e-13
        c4_e2l_b = np.zeros((a_deca, 1)) + 0.0
        # cal coeffs: lab calibrated values to scientific values
        c0_l2s_b = np.zeros((a_deca, 1)) + 79.12599
        c1_l2s_b = np.zeros((a_deca, 1)) + -9.58863
        c2_l2s_b = np.zeros((a_deca, 1)) + 0.53886
        c3_l2s_b = np.zeros((a_deca, 1)) + -0.01432
        c4_l2s_b = np.zeros((a_deca, 1)) + 1.38009e-4

        # timing test
        self.profile(stats, sfl_thsph_temp_int, ts_rawdec_b,
                     c0_e2l_b, c1_e2l_b, c2_e2l_b, c3_e2l_b, c4_e2l_b,
                     c0_l2s_b, c1_l2s_b, c2_l2s_b, c3_l2s_b, c4_l2s_b)

    def test_sfl_thsph_temp_ref(self):
        stats = []
        # create 10000 data packets
        # input, raw decimal counts
        ts_rawdec_r = np.zeros((a_deca, 1)) + 8770.0
        # cal coeffs: engineering values to lab calibrated values
        c0_e2l_r = np.zeros((a_deca, 1)) + 0.05935
        c1_e2l_r = np.zeros((a_deca, 1)) + 0.00099151
        c2_e2l_r = np.zeros((a_deca, 1)) + 3.82028e-10
        c3_e2l_r = np.zeros((a_deca, 1)) + 4.54486e-13
        c4_e2l_r = np.zeros((a_deca, 1)) + 0.0
        # cal coeffs: lab calibrated values to scientific values
        c0_l2s_r = np.zeros((a_deca, 1)) + 79.12599
        c1_l2s_r = np.zeros((a_deca, 1)) + -9.58863
        c2_l2s_r = np.zeros((a_deca, 1)) + 0.53886
        c3_l2s_r = np.zeros((a_deca, 1)) + -0.01432
        c4_l2s_r = np.zeros((a_deca, 1)) + 1.38009e-4

        # timing test
        self.profile(stats, sfl_thsph_temp_ref, ts_rawdec_r,
                     c0_e2l_r, c1_e2l_r, c2_e2l_r, c3_e2l_r, c4_e2l_r,
                     c0_l2s_r, c1_l2s_r, c2_l2s_r, c3_l2s_r, c4_l2s_r)

    def test_sfl_thsph_temp_tcl(self):
        stats = []
        # create 10000 data packets
        # input, raw decimal counts
        tc_rawdec_L = np.zeros((a_deca, 1)) + 16012.0
        # cal coeffs: engineering values to lab calibrated values
        c0_e2l_L = np.zeros((a_deca, 1)) - 0.00055
        c1_e2l_L = np.zeros((a_deca, 1)) + 1.0
        c2_e2l_L = np.zeros((a_deca, 1)) + 0.0
        c3_e2l_L = np.zeros((a_deca, 1)) + 0.0
        c4_e2l_L = np.zeros((a_deca, 1)) + 0.0
        # cal coeffs: lab calibrated values to scientific values
        c0_l2s_L = np.zeros((a_deca, 1)) - 0.00444
        c1_l2s_L = np.zeros((a_deca, 1)) + 17.06172
        c2_l2s_L = np.zeros((a_deca, 1)) - 0.23532
        c3_l2s_L = np.zeros((a_deca, 1)) + 0.00702
        c4_l2s_L = np.zeros((a_deca, 1)) - 0.000122268
        c5_l2s_L = np.zeros((a_deca, 1)) + 0.000000932483

        # timing test
        self.profile(stats, sfl_thsph_temp_tcl, tc_rawdec_L,
                     c0_e2l_L, c1_e2l_L, c2_e2l_L, c3_e2l_L, c4_e2l_L,
                     c0_l2s_L, c1_l2s_L, c2_l2s_L, c3_l2s_L, c4_l2s_L, c5_l2s_L)

    def test_sfl_thsph_temp_tch(self):
        stats = []
        # create 10000 data packets
        # input, raw decimal counts
        tc_rawdec_H = np.zeros((a_deca, 1)) + 4237.0
        # cal coeffs: engineering values to lab calibrated values
        c0_e2l_H = np.zeros((a_deca, 1)) - 0.00055
        c1_e2l_H = np.zeros((a_deca, 1)) + 1.0
        c2_e2l_H = np.zeros((a_deca, 1)) + 0.0
        c3_e2l_H = np.zeros((a_deca, 1)) + 0.0
        c4_e2l_H = np.zeros((a_deca, 1)) + 0.0
        # cal coeffs: lab calibrated values to scientific values
        c0_l2s_H = np.zeros((a_deca, 1)) - 0.00444
        c1_l2s_H = np.zeros((a_deca, 1)) + 17.06172
        c2_l2s_H = np.zeros((a_deca, 1)) - 0.23532
        c3_l2s_H = np.zeros((a_deca, 1)) + 0.00702
        c4_l2s_H = np.zeros((a_deca, 1)) - 0.000122268
        c5_l2s_H = np.zeros((a_deca, 1)) + 0.000000932483

        # timing test
        self.profile(stats, sfl_thsph_temp_tch, tc_rawdec_H,
                     c0_e2l_H, c1_e2l_H, c2_e2l_H, c3_e2l_H, c4_e2l_H,
                     c0_l2s_H, c1_l2s_H, c2_l2s_H, c3_l2s_H, c4_l2s_H, c5_l2s_H)

    def test_sfl_thsph_temp_tl(self):
        stats = []
        # create 10000 data packets
        # input, raw decimal counts:L
        tc_rawdec_L = np.zeros((a_deca, 1)) + 16012.0
        # cal coeffs: engineering values to lab calibrated values
        c0_e2l_L = np.zeros((a_deca, 1)) - 0.00055
        c1_e2l_L = np.zeros((a_deca, 1)) + 1.0
        c2_e2l_L = np.zeros((a_deca, 1)) + 0.0
        c3_e2l_L = np.zeros((a_deca, 1)) + 0.0
        c4_e2l_L = np.zeros((a_deca, 1)) + 0.0
        # cal coeffs: lab calibrated values to scientific values
        c0_l2s_L = np.zeros((a_deca, 1)) - 0.00444
        c1_l2s_L = np.zeros((a_deca, 1)) + 17.06172
        c2_l2s_L = np.zeros((a_deca, 1)) - 0.23532
        c3_l2s_L = np.zeros((a_deca, 1)) + 0.00702
        c4_l2s_L = np.zeros((a_deca, 1)) - 0.000122268
        c5_l2s_L = np.zeros((a_deca, 1)) + 0.000000932483
        # calibration constants: final linear calibrations
        c0_s2f_L = 1.68019
        c1_s2f_L = 0.95567
        # input, raw decimal counts:ref
        ts_rawdec_r = np.zeros((a_deca, 1)) + 8770.0
        # cal coeffs: engineering values to lab calibrated values
        c0_e2l_r = np.zeros((a_deca, 1)) + 0.05935
        c1_e2l_r = np.zeros((a_deca, 1)) + 0.00099151
        c2_e2l_r = np.zeros((a_deca, 1)) + 3.82028e-10
        c3_e2l_r = np.zeros((a_deca, 1)) + 4.54486e-13
        c4_e2l_r = np.zeros((a_deca, 1)) + 0.0
        # cal coeffs: lab calibrated values to scientific values
        c0_l2s_r = np.zeros((a_deca, 1)) + 79.12599
        c1_l2s_r = np.zeros((a_deca, 1)) + -9.58863
        c2_l2s_r = np.zeros((a_deca, 1)) + 0.53886
        c3_l2s_r = np.zeros((a_deca, 1)) + -0.01432
        c4_l2s_r = np.zeros((a_deca, 1)) + 1.38009e-4

        # timing test
        self.profile(stats, sfl_thsph_temp_tl, tc_rawdec_L,
                     c0_e2l_L, c1_e2l_L, c2_e2l_L, c3_e2l_L, c4_e2l_L,
                     c0_l2s_L, c1_l2s_L, c2_l2s_L, c3_l2s_L, c4_l2s_L, c5_l2s_L,
                     c0_s2f_L, c1_s2f_L,
                     ts_rawdec_r,
                     c0_e2l_r, c1_e2l_r, c2_e2l_r, c3_e2l_r, c4_e2l_r,
                     c0_l2s_r, c1_l2s_r, c2_l2s_r, c3_l2s_r, c4_l2s_r)

    def test_sfl_thsph_temp_th(self):
        stats = []
        # create 10000 data packets
        # input, raw decimal counts:L
        tc_rawdec_H = np.zeros((a_deca, 1)) + 4237.0
        # cal coeffs: engineering values to lab calibrated values
        c0_e2l_H = np.zeros((a_deca, 1)) - 0.00055
        c1_e2l_H = np.zeros((a_deca, 1)) + 1.0
        c2_e2l_H = np.zeros((a_deca, 1)) + 0.0
        c3_e2l_H = np.zeros((a_deca, 1)) + 0.0
        c4_e2l_H = np.zeros((a_deca, 1)) + 0.0
        # cal coeffs: lab calibrated values to scientific values
        c0_l2s_H = np.zeros((a_deca, 1)) - 0.00444
        c1_l2s_H = np.zeros((a_deca, 1)) + 17.06172
        c2_l2s_H = np.zeros((a_deca, 1)) - 0.23532
        c3_l2s_H = np.zeros((a_deca, 1)) + 0.00702
        c4_l2s_H = np.zeros((a_deca, 1)) - 0.000122268
        c5_l2s_H = np.zeros((a_deca, 1)) + 0.000000932483
        # calibration constants: final linear calibrations
        c0_s2f_H = 1.68019
        c1_s2f_H = 0.95567
        # input, raw decimal counts:ref
        ts_rawdec_r = np.zeros((a_deca, 1)) + 8770.0
        # cal coeffs: engineering values to lab calibrated values
        c0_e2l_r = np.zeros((a_deca, 1)) + 0.05935
        c1_e2l_r = np.zeros((a_deca, 1)) + 0.00099151
        c2_e2l_r = np.zeros((a_deca, 1)) + 3.82028e-10
        c3_e2l_r = np.zeros((a_deca, 1)) + 4.54486e-13
        c4_e2l_r = np.zeros((a_deca, 1)) + 0.0
        # cal coeffs: lab calibrated values to scientific values
        c0_l2s_r = np.zeros((a_deca, 1)) + 79.12599
        c1_l2s_r = np.zeros((a_deca, 1)) + -9.58863
        c2_l2s_r = np.zeros((a_deca, 1)) + 0.53886
        c3_l2s_r = np.zeros((a_deca, 1)) + -0.01432
        c4_l2s_r = np.zeros((a_deca, 1)) + 1.38009e-4

        # timing test
        self.profile(stats, sfl_thsph_temp_tl, tc_rawdec_H,
                     c0_e2l_H, c1_e2l_H, c2_e2l_H, c3_e2l_H, c4_e2l_H,
                     c0_l2s_H, c1_l2s_H, c2_l2s_H, c3_l2s_H, c4_l2s_H, c5_l2s_H,
                     c0_s2f_H, c1_s2f_H,
                     ts_rawdec_r,
                     c0_e2l_r, c1_e2l_r, c2_e2l_r, c3_e2l_r, c4_e2l_r,
                     c0_l2s_r, c1_l2s_r, c2_l2s_r, c3_l2s_r, c4_l2s_r)

    def test_sfl_trhph_vfltemp(self):
        stats = []
        # create 10000 data packets
        V_s = np.zeros((a_deca, 1)) + 1.931
        V_c = np.zeros((a_deca, 1)) + 1.077

        # calibration constants
        a = 1.98e-9
        b = -2.45e-6
        c = 9.28e-4
        d = -0.0888
        e = 0.731

        # timing test
        self.profile(stats, sfl_trhph_vfltemp, V_s, V_c, a, b, c, d, e)

    def test_sfl_trhph_vflorp(self):
        stats = []
        # create 10000 data packets
        V = np.zeros((a_deca, 1)) + 1.541

        # calibration coefficients
        offset = 2004.0
        gain = 4.0

        # timing test
        self.profile(stats, sfl_trhph_vflorp, V, offset, gain)

    def test_sfl_sflpres_l1(self):
        stats = []
        # create 10000 data packets
        P = np.zeros((a_deca, 1)) + 14.868

        # timing test
        self.profile(stats, sfl_sflpres_l1, P)

    def test_sfl_trhph_chloride(self):
        stats = []

        # set of data packets which will use the various resistivity branches
        # in the code and has 1 out of range temperature value
        test_array = np.array([
            [0.440,  4.095,  4.095,  105.4],
            [0.380,  4.095,  4.095,  241.9],
            [0.320,  4.095,  4.095,  374.2],
            [0.184,  0.915,  4.064,  105.4],
            [0.172,  0.857,  4.082,  374.2],
            [0.183,  0.926,  4.076,  222.0],
            [0.131,  0.673,  3.293,  325.8],
            [0.133,  0.678,  3.396,  999.9],
            [0.135,  0.681,  2.000,  333.4],
            [0.135,  0.681,  1.000,  333.2]
        ])

        # create 10000 data packets
        tile_value = np.round(a_deca/test_array.shape[0])
        test_array = np.tile(test_array, (tile_value, 1))
        V_R1 = test_array[:, 0]
        V_R2 = test_array[:, 1]
        V_R3 = test_array[:, 2]
        T = test_array[:, 3]

        # timing test
        self.profile(stats, sfl_trhph_chloride, V_R1, V_R2, V_R3, T)

