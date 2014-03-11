#!/usr/bin/env python
from ion_functions.data.perf.test_performance import PerformanceTestCase, a_deca
from ion_functions.data.sfl_functions import (sfl_trhph_vfltemp,
                                              sfl_trhph_vflorp,
                                              sfl_trhph_chloride,
                                              sfl_sflpres_l1)
import numpy as np


class TestSFLPerformance(PerformanceTestCase):

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

