#!/usr/bin/env python
from ion_functions.data.perf.test_performance import PerformanceTestCase, a_deca
from ion_functions.data.co2_functions import pco2_thermistor, pco2_blank, pco2_pco2wat, pco2_co2flux
import numpy as np


class TestCO2Performance(PerformanceTestCase):
    def setUp(self):
        ### test data for the PCO2W related calculcations
        s = '*7E2704CBACF1230081007C0B2B00BF080800DB01BF0390007C00790B3000C7080B00EA0C5606C80A'
        self.traw = int(s[75:79], 16)
        self.ea434 = 19706.
        self.ea620 = 34.
        self.eb434 = 3073.
        self.eb620 = 44327.
        self.calt = 16.5
        self.cala = 0.0459
        self.calb = 0.6257
        self.calc = -1.5406

        self.light = np.zeros(14, dtype=np.int)
        self.mtype = int(s[5:7], 16)
        self.traw = int(s[75:79], 16)
        strt = 15
        step = 4
        for j in range(14):
            self.light[j] = int(s[strt:strt+step], 16)
            strt += step

        self.a434braw = self.light[6]
        self.a434blnk = pco2_blank(self.a434braw)
        self.a620braw = self.light[7]
        self.a620blnk = pco2_blank(self.a620braw)

        ### test data for the PCO2A CO2FLUX calculations
        self.flux_data = np.array([
            [360, 390, 5, 0, 34, -2.063e-08],
            [360, 390, 5, 0, 35, -2.052e-08],
            [360, 390, 5, 10, 34, -1.942e-08],
            [360, 390, 5, 10, 35, -1.932e-08],
            [360, 390, 5, 20, 34, -1.869e-08],
            [360, 390, 5, 20, 35, -1.860e-08],
            [360, 390, 10, 0, 34, -8.250e-08],
            [360, 390, 10, 0, 35, -8.207e-08],
            [360, 390, 10, 10, 34, -7.767e-08],
            [360, 390, 10, 10, 35, -7.728e-08],
            [360, 390, 10, 20, 34, -7.475e-08],
            [360, 390, 10, 20, 35, -7.440e-08],
            [360, 390, 20, 0, 34, -3.300e-07],
            [360, 390, 20, 0, 35, -3.283e-07],
            [360, 390, 20, 10, 34, -3.107e-07],
            [360, 390, 20, 10, 35, -3.091e-07],
            [360, 390, 20, 20, 34, -2.990e-07],
            [360, 390, 20, 20, 35, -2.976e-07],
            [400, 390, 5, 0, 34, 6.875e-09],
            [400, 390, 5, 0, 35, 6.839e-09],
            [400, 390, 5, 10, 34, 6.472e-09],
            [400, 390, 5, 10, 35, 6.440e-09],
            [400, 390, 5, 20, 34, 6.229e-09],
            [400, 390, 5, 20, 35, 6.200e-09],
            [400, 390, 10, 0, 34, 2.750e-08],
            [400, 390, 10, 0, 35, 2.736e-08],
            [400, 390, 10, 10, 34, 2.589e-08],
            [400, 390, 10, 10, 35, 2.576e-08],
            [400, 390, 10, 20, 34, 2.492e-08],
            [400, 390, 10, 20, 35, 2.480e-08],
            [400, 390, 20, 0, 34, 1.100e-07],
            [400, 390, 20, 0, 35, 1.094e-07],
            [400, 390, 20, 10, 34, 1.036e-07],
            [400, 390, 20, 10, 35, 1.030e-07],
            [400, 390, 20, 20, 34, 9.966e-08],
            [400, 390, 20, 20, 35, 9.920e-08],
            [440, 390, 5, 0, 34, 3.438e-08],
            [440, 390, 5, 0, 35, 3.420e-08],
            [440, 390, 5, 10, 34, 3.236e-08],
            [440, 390, 5, 10, 35, 3.220e-08],
            [440, 390, 5, 20, 34, 3.114e-08],
            [440, 390, 5, 20, 35, 3.100e-08],
            [440, 390, 10, 0, 34, 1.375e-07],
            [440, 390, 10, 0, 35, 1.368e-07],
            [440, 390, 10, 10, 34, 1.294e-07],
            [440, 390, 10, 10, 35, 1.288e-07],
            [440, 390, 10, 20, 34, 1.246e-07],
            [440, 390, 10, 20, 35, 1.240e-07],
            [440, 390, 20, 0, 34, 5.500e-07],
            [440, 390, 20, 0, 35, 5.471e-07],
            [440, 390, 20, 10, 34, 5.178e-07],
            [440, 390, 20, 10, 35, 5.152e-07],
            [440, 390, 20, 20, 34, 4.983e-07],
            [440, 390, 20, 20, 35, 4.960e-07]
        ])

    def test_pco2_blanks(self):
        stats = []

        # create 10000 data points
        data = np.ones(a_deca, dtype='int16')
        braw = data * self.a434braw
        self.profile(stats, pco2_blank, braw)
        braw = data * self.a620braw
        self.profile(stats, pco2_blank, braw)

    def test_pco2_thermistor(self):
        stats = []

        # create 10000 data points
        data = np.ones(a_deca, dtype='int16')
        traw = data * self.traw
        self.profile(stats, pco2_thermistor, traw)

    def test_pco2_calc_pco2(self):
        stats = []

        # create 10000 data points
        data = np.ones(a_deca, dtype='int16')
        mtype = data * self.mtype
        traw = data * self.traw
        ea434 = data * self.ea434
        eb434 = data * self.eb434
        ea620 = data * self.ea620
        eb620 = data * self.eb620
        calt = data * self.calt
        cala = data * self.cala
        calb = data * self.calb
        calc = data * self.calc
        a434blnk = data * self.a434blnk
        a620blnk = data * self.a620blnk
        light = np.ones((a_deca, 14)) * self.light

        tout = pco2_thermistor(traw)

        self.profile(stats, pco2_pco2wat, mtype, light, tout, ea434, eb434,
                     ea620, eb620, calt, cala, calb, calc, a434blnk, a620blnk)

    def test_pco2_co2flux(self):
        stats = []

        # generate 10000 data points from input array
        nPts = np.round(a_deca / self.flux_data.shape[0])
        pco2w = np.repeat(self.flux_data[:, 0], nPts)
        pco2a = np.repeat(self.flux_data[:, 1], nPts)
        u10 = np.repeat(self.flux_data[:, 2], nPts)
        t = np.repeat(self.flux_data[:, 3], nPts)
        s = np.repeat(self.flux_data[:, 4], nPts)

        self.profile(stats, pco2_co2flux, pco2w, pco2a, u10, t, s)
