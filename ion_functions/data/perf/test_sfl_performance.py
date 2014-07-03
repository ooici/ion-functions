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
    # Performance tests for seafloor instruments THSPH and TRHPH

    def setUp(self):
        ### test inputs for THSPHTE products
        self.ts_rawdec_b = np.tile(8185.0, a_deca)
        self.ts_rawdec_r = np.tile(8758.0, a_deca)
        self.tc_rawdec_L = np.tile(16009.0, a_deca)
        self.tc_rawdec_H = np.tile(4236.0, a_deca)

        ### calibration arrays for THSPHTE products
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

        ### tiled calibration arrays for THSPHTE products
        tile_spec = (a_deca, 1)
        self.e2l_b = np.tile(e2l_b, tile_spec)
        self.l2s_b = np.tile(l2s_b, tile_spec)
        self.e2l_r = np.tile(e2l_r, tile_spec)
        self.l2s_r = np.tile(l2s_r, tile_spec)
        self.e2l_L = np.tile(e2l_L, tile_spec)
        self.l2s_L = np.tile(l2s_L, tile_spec)
        self.e2l_H = np.tile(e2l_H, tile_spec)
        self.l2s_H = np.tile(l2s_H, tile_spec)
        self.s2v_r = np.tile(s2v_r, tile_spec)

    # Performance tests for seafloor instruments THSPH
    def test_sfl_thsph_temp_int(self):
        stats = []
        # timing test
        self.profile(stats, sfl_thsph_temp_int, self.ts_rawdec_b, self.e2l_b, self.l2s_b)

    def test_sfl_thsph_temp_ref(self):
        stats = []
        # timing test
        self.profile(stats, sfl_thsph_temp_ref, self.ts_rawdec_r, self.e2l_r, self.l2s_r)

    def test_sfl_thsph_temp_tcl(self):
        stats = []
        # timing test
        self.profile(stats, sfl_thsph_temp_tcl, self.tc_rawdec_L, self.e2l_L, self.l2s_L)

    def test_sfl_thsph_temp_tch(self):
        stats = []
        # timing test
        self.profile(stats, sfl_thsph_temp_tch, self.tc_rawdec_H, self.e2l_H, self.l2s_H)

    def test_sfl_thsph_temp_tl(self):
        stats = []
        # timing test
        self.profile(stats, sfl_thsph_temp_tl, self.tc_rawdec_L, self.e2l_L, self.l2s_L,
                     self.ts_rawdec_r, self.e2l_r, self.l2s_r, self.s2v_r)

    def test_sfl_thsph_temp_th(self):
        stats = []
        # timing test
        self.profile(stats, sfl_thsph_temp_th, self.tc_rawdec_H, self.e2l_H, self.l2s_H,
                     self.ts_rawdec_r, self.e2l_r, self.l2s_r, self.s2v_r)

    # Performance tests for seafloor instruments TRHPH

    def test_sfl_trhph_vfltemp(self):
        stats = []
        # create 10000 data packets
        V_s = np.zeros(a_deca) + 1.931
        V_c = np.zeros(a_deca) + 1.077

        # calibration constants
        tc_slope = np.zeros(a_deca) + 4.22e-5
        ts_slope = np.zeros(a_deca) + 0.003
        #c3 = np.zeros(a_deca) + -1.00e-6
        #c2 = np.zeros(a_deca) + 7.00e-6
        #c1 = np.zeros(a_deca) + 0.0024
        #c0 = np.zeros(a_deca) + 0.015

        # timing test
        #self.profile(stats, sfl_trhph_vfltemp, V_s, V_c, tc_slope, ts_slope, c0, c1, c2, c3)
        self.profile(stats, sfl_trhph_vfltemp, V_s, V_c, tc_slope, ts_slope)

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

