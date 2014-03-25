#!/usr/bin/env python
from ion_functions.data.perf.test_performance import PerformanceTestCase, a_deca
from ion_functions.data import opt_functions as optfunc
import numpy as np


class TestOPTAAPerformance(PerformanceTestCase):
    def setUp(self):

        ### realistic values for ac-s data packets:
        n_wvl = 90  # number of wavelengths specified in DPS is incorrect
        wvl_tile = n_wvl/6  # test arrays have 6 values
        n_tbins = 35
        tbin_tile = n_tbins/7  # test array has 7 tbin values

        ### test data common to both OPTATTN and OPTABSN
        self.traw = 48355
        self.tcal = 20.0
        self.T = 12.0
        self.S = 35.0
        # tbin values are used in an interpolation algorithm; make sure
        # their values are monotonic
        self.tbins = np.array(range(n_tbins))

        ### test data for OPTATTN
        self.c_sig = np.tile([150., 225., 200., 350., 450., 495.], (1, wvl_tile))
        self.c_ref = np.tile([550., 540., 530., 520., 510., 500.], (1, wvl_tile))
        self.c_off = np.tile([1.35, 1.30, 1.25, 1.20, 1.15, 1.10], (1, wvl_tile))
        self.c_wvl = np.tile([510., 540., 580., 630., 670., 710.], (1, wvl_tile))
        self.tc_arr = np.tile([
            [0.0, -0.004929, -0.004611, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.004611, -0.004418, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.004418, -0.004355, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.004355, -0.004131, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.004131, -0.003422, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.003422, -0.002442, 0.0, 0.0, 0.0, 0.0]], (wvl_tile, tbin_tile))

        ### test data for OPTABSN
        self.a_sig = np.tile([250., 300., 210., 430., 470., 495.], (1, wvl_tile))
        self.a_ref = np.tile([450., 460., 470., 480., 490., 500.], (1, wvl_tile))
        self.a_off = np.tile([0.35, 0.30, 0.25, 0.20, 0.15, 0.10], (1, wvl_tile))
        self.a_wvl = np.tile([500., 550., 600., 650., 700., 715.], (1, wvl_tile))
        # note, even though here ta_arr and tc_arr are identical, in actual calibration
        #   data they will be different.
        self.ta_arr = np.tile([
            [0.0, -0.004929, -0.004611, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.004611, -0.004418, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.004418, -0.004355, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.004355, -0.004131, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.004131, -0.003422, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.003422, -0.002442, 0.0, 0.0, 0.0, 0.0]], (wvl_tile, tbin_tile))

        self.cpd_ts = np.tile([6.553771, 4.807914, 5.156010, 2.788715, 1.655607, 1.171965],
                             (1, wvl_tile))

        # test data for ocr-507 (spectir = downwelling irradiance)
        #     counts         offset          scale        mrsn
        self.ocr_507 = np.transpose(np.array([
            [2148370944, 2148377867.8, 2.09023117662E-07, 1.368],
            [2200000000, 2148218092.4, 2.06543624674E-07, 1.410],
            [2300000000, 2147607229.7, 2.12484770952E-07, 1.365],
            [2400000000, 2147789959.1, 2.07241106309E-07, 1.354],
            [2500000000, 2148047456.7, 1.99358530187E-07, 1.372],
            [2600000000, 2147335412.8, 2.06033896796E-07, 1.404],
            [2700000000, 2146998228.4, 2.14806273478E-07, 1.347]
        ]))

    def test_opt_beam_attenuation(self):
        stats = []
        # create 10000 data packets
        # common input variables: traw, Tcal, T, and PS
        traw = np.tile(self.traw, (a_deca, 1))
        tcal = np.tile(self.tcal, (a_deca, 1))
        T = np.tile(self.T, (a_deca, 1))
        PS = np.tile(self.S, (a_deca, 1))
        tbins = np.tile(self.tbins, (a_deca, 1))
        # variables unique to beam attenuation: 1D -> 2D
        c_sig = np.tile(self.c_sig, (a_deca, 1))
        c_ref = np.tile(self.c_ref, (a_deca, 1))
        c_off = np.tile(self.c_off, (a_deca, 1))
        c_wvl = np.tile(self.c_wvl, (a_deca, 1))
        # variables unique to beam attenuation: 2D -> 3D
        tc_arr = np.tile(self.tc_arr, (a_deca, 1, 1))
        # timing test
        self.profile(stats, optfunc.opt_beam_attenuation, c_ref, c_sig, traw, c_wvl, c_off, tcal,
                     tbins, tc_arr, T, PS)

    def test_opt_optical_absorption(self):
        stats = []
        # create 10000 data packets
        # common input variables: traw, Tcal, T, and PS
        traw = np.tile(self.traw, (a_deca, 1))
        tcal = np.tile(self.tcal, (a_deca, 1))
        T = np.tile(self.T, (a_deca, 1))
        PS = np.tile(self.S, (a_deca, 1))
        tbins = np.tile(self.tbins, (a_deca, 1))
        # variables unique to beam attenuation: 1D -> 2D
        a_sig = np.tile(self.a_sig, (a_deca, 1))
        a_ref = np.tile(self.a_ref, (a_deca, 1))
        a_off = np.tile(self.a_off, (a_deca, 1))
        a_wvl = np.tile(self.a_wvl, (a_deca, 1))
        # variables unique to beam attenuation: 2D -> 3D
        ta_arr = np.tile(self.ta_arr, (a_deca, 1, 1))
        cpd_ts = np.tile(self.cpd_ts, (a_deca, 1))
        c_wvl = np.tile(self.c_wvl, (a_deca, 1))
        # timing test
        self.profile(stats, optfunc.opt_optical_absorption, a_ref, a_sig, traw, a_wvl, a_off, tcal,
                     tbins, ta_arr, cpd_ts, c_wvl, T, PS)

    def test_opt_ocr507_irradiance(self):
        stats = []
        # create 10000 data packets
        counts = np.tile(self.ocr_507[0, :], (a_deca, 1))
        offset = np.tile(self.ocr_507[1, :], (a_deca, 1))
        scale = np.tile(self.ocr_507[2, :], (a_deca, 1))
        immersion_factor = np.tile(self.ocr_507[3, :], (a_deca, 1))
        # timing test
        self.profile(stats, optfunc.opt_ocr507_irradiance, counts, offset, scale, immersion_factor)
