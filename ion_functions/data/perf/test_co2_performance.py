#!/usr/bin/env python
from ion_functions.data.perf.test_performance import PerformanceTestCase
from ion_functions.data.co2_functions import pco2_thermistor, pco2_abs434_blank, pco2_abs620_blank, pco2_pco2wat, pco2_co2flux
import numpy as np
from ion_functions.utils import fill_value


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
        self.a434blnk = fill_value
        self.a620blnk = fill_value

        # expected outputs
        self.therm = np.array([18.8526, 18.8765, 18.9245, 18.9485,
                               18.9485, 18.9485, 18.8765, 19.0686,
                               19.0686, 19.0446, 18.9725])
        self.pco2 = np.array([fill_value, 294.1720, 311.3361, 319.0101,
                              319.8925, 319.8950, 305.8104, 317.9661,
                              284.3676, 280.2324, 280.0354])

        self.light = np.zeros(14, dtype=np.int)
        self.mtype = int(s[5:7], 16)
        self.traw = int(s[75:79], 16)
        strt = 15
        step = 4
        for j in range(14):
            self.light[j] = int(s[strt:strt+step], 16)
            strt += step

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

    def test_pco2_thermistor(self):
        stats = []

        sample_data = np.empty(3600 * 24 * 365, dtype='int32')
        sample_data.fill(self.traw)
        self.profile(stats, pco2_thermistor, sample_data)

    def test_pco2_calc_pco2(self):
        stats = []

        light = self.light
        mtype = self.mtype
        traw = np.empty(3600 * 24 * 365, dtype=np.int)
        tout = pco2_thermistor(traw)
        a434blnk = pco2_abs434_blank(mtype, light, self.a434blnk)
        a620blnk = pco2_abs620_blank(mtype, light, self.a620blnk)

        self.profile(stats, pco2_pco2wat, mtype, light, tout, self.ea434, self.eb434,
                     self.ea620, self.eb620, self.calt, self.cala, self.calb,
                     self.calc, a434blnk, a620blnk)

    def test_pco2_co2flux(self):
        stats = []

        # generate a very large dataset (~62 years worth of data, sampled
        # hourly).
        pco2w = np.repeat(self.flux_data[:, 0], 10000)
        pco2a = np.repeat(self.flux_data[:, 1], 10000)
        u10 = np.repeat(self.flux_data[:, 2], 10000)
        t = np.repeat(self.flux_data[:, 3], 10000)
        s = np.repeat(self.flux_data[:, 4], 10000)

        self.profile(stats, pco2_co2flux, pco2w, pco2a, u10, t, s)
