#!/usr/bin/env python
from ion_functions.data.perf.test_performance import PerformanceTestCase, a_deca
from ion_functions.data.sfl_functions import (sfl_trhph_vfltemp,
                                              sfl_trhph_vflorp,
                                              sfl_trhph_chloride,
                                              sfl_thsph_temp_int,
                                              sfl_thsph_temp_ref,
                                              sfl_thsph_temp_tcl,
                                              sfl_thsph_temp_tch,
                                              sfl_thsph_temp_tl,
                                              sfl_thsph_temp_th,
                                              sfl_thsph_hydrogen,
                                              sfl_thsph_sulfide,
                                              sfl_thsph_ph,
                                              sfl_thsph_ph_acl,
                                              sfl_thsph_ph_noref,
                                              sfl_thsph_ph_noref_acl,
                                              sfl_sflpres_rtime,
                                              sfl_sflpres_tide,
                                              sfl_sflpres_wave,
                                              )
import numpy as np


class TestSFLPerformance(PerformanceTestCase):
    # Performance tests for seafloor instruments THSPH and TRHPH

    def setUp(self):

        ### test inputs for THSPH L2 data products
        self.counts_h2 = np.tile(4907, a_deca)
        self.counts_hs = np.tile(3806, a_deca)
        self.counts_ysz = np.tile(7807, a_deca)
        self.counts_agcl = np.tile(7801, a_deca)
        self.temperature = np.tile(300.0, a_deca)
        self.chloride = np.tile(400.0, a_deca)

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

        ### tiled calibration arrays for THSPH L2 products
        tile_spec = (a_deca, 1)
        self.e2l_h2 = np.tile(e2l_h2, tile_spec)
        self.e2l_hs = np.tile(e2l_hs, tile_spec)
        self.e2l_ysz = np.tile(e2l_ysz, tile_spec)
        self.e2l_agcl = np.tile(e2l_agcl, tile_spec)
        self.arr_hgo = np.tile(arr_hgo, tile_spec)
        self.arr_agcl = np.tile(arr_agcl, tile_spec)
        self.arr_agclref = np.tile(arr_agclref, tile_spec)
        self.arr_tac = np.tile(arr_tac, tile_spec)
        self.arr_tbc1 = np.tile(arr_tbc1, tile_spec)
        self.arr_tbc2 = np.tile(arr_tbc2, tile_spec)
        self.arr_tbc3 = np.tile(arr_tbc3, tile_spec)
        self.arr_logkfh2g = np.tile(arr_logkfh2g, tile_spec)
        self.arr_eh2sg = np.tile(arr_eh2sg, tile_spec)
        self.arr_yh2sg = np.tile(arr_yh2sg, tile_spec)

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

    # Performance tests for seafloor instruments THSPH, L2 data products

    def test_sfl_thsph_hydrogen(self):
        stats = []
        # timing test
        self.profile(stats, sfl_thsph_hydrogen, self.counts_h2, self.counts_ysz,
                     self.temperature, self.e2l_h2, self.e2l_ysz, self.arr_hgo,
                     self.arr_logkfh2g)

    def test_sfl_thsph_sulfide(self):
        stats = []
        # timing test
        self.profile(stats, sfl_thsph_sulfide, self.counts_hs, self.counts_ysz,
                     self.temperature, self.e2l_hs, self.e2l_ysz, self.arr_hgo,
                     self.arr_logkfh2g, self.arr_eh2sg, self.arr_yh2sg)

    def test_sfl_thsph_ph(self):
        stats = []
        # timing test
        self.profile(stats, sfl_thsph_ph, self.counts_ysz, self.counts_agcl,
                     self.temperature, self.e2l_ysz, self.e2l_agcl, self.arr_hgo,
                     self.arr_agcl, self.arr_tac, self.arr_tbc1, self.arr_tbc2,
                     self.arr_tbc3, self.chloride)

    def test_sfl_thsph_ph_acl(self):
        stats = []
        # timing test
        self.profile(stats, sfl_thsph_ph_acl, self.counts_ysz, self.counts_agcl,
                     self.temperature, self.e2l_ysz, self.e2l_agcl, self.arr_hgo,
                     self.arr_agcl, self.arr_tac, self.arr_tbc1, self.arr_tbc2,
                     self.arr_tbc3)

    def test_sfl_thsph_ph_noref(self):
        stats = []
        # timing test
        self.profile(stats, sfl_thsph_ph_noref, self.counts_ysz, self.temperature,
                     self.arr_agclref, self.e2l_ysz, self.arr_hgo, self.arr_agcl,
                     self.arr_tac, self.arr_tbc1, self.arr_tbc2, self.arr_tbc3,
                     self.chloride)

    def test_sfl_thsph_ph_noref_acl(self):
        stats = []
        # timing test
        self.profile(stats, sfl_thsph_ph_noref_acl, self.counts_ysz, self.temperature,
                     self.arr_agclref, self.e2l_ysz, self.arr_hgo, self.arr_agcl,
                     self.arr_tac, self.arr_tbc1, self.arr_tbc2, self.arr_tbc3)

    # Performance tests for seafloor instruments THSPH, 6 THSPHTE data products

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

    # Performance tests for seafloor instruments PRESF

    def test_sfl_sflpres_rtime(self):
        stats = []
        # create 10000 data packets
        P = np.zeros(a_deca) + 14.868

        # timing test
        self.profile(stats, sfl_sflpres_rtime, P)

    def test_sfl_sflpres_tide(self):
        stats = []
        # create 10000 data packets
        p_dec_tide = np.zeros(a_deca, dtype=int) + 4175754
        m = np.zeros(a_deca) + 279620.2
        b = np.zeros(a_deca) + 18641.3
        slope = np.zeros(a_deca) + 1.0
        offset = np.zeros(a_deca) + 0.0
        # timing test
        self.profile(stats, sfl_sflpres_tide, p_dec_tide, b, m, slope, offset)

    def test_sfl_sflpres_wave(self):
        stats = []
        # create 10000 data packets
        # i don't know how many p_dec_wave values there are per burst ...
        nval_per_burst = 20
        p_dec_wave = np.tile(8900312, (a_deca, nval_per_burst))
        #
        lots_of_zeros = np.zeros(a_deca, dtype=int)
        ptcn = 43746280 + lots_of_zeros
        u0 = 5.856409e+00 + lots_of_zeros
        y1 = -3.987838e+03 + lots_of_zeros
        y2 = -1.049603e+04 + lots_of_zeros
        y3 = 0.000000e+00 + lots_of_zeros
        c1 = 2.305367e+02 + lots_of_zeros
        c2 = 1.198422e+01 + lots_of_zeros
        c3 = -2.401512e+02 + lots_of_zeros
        d1 = 4.095400e-02 + lots_of_zeros
        d2 = 0.000000e+00 + lots_of_zeros
        t1 = 2.781994e+01 + lots_of_zeros
        t2 = 6.760780e-01 + lots_of_zeros
        t3 = 1.761829e+01 + lots_of_zeros
        t4 = 6.000932e+00 + lots_of_zeros
        poff = 0.0 + lots_of_zeros
        slope = 1.0 + lots_of_zeros
        offset = 0.0 + lots_of_zeros

        # timing test
        self.profile(stats, sfl_sflpres_wave, ptcn, p_dec_wave, u0, y1, y2, y3,
                     c1, c2, c3, d1, d2, t1, t2, t3, t4, poff, slope, offset)

