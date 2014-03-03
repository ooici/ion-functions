#!/usr/bin/env python
from ion_functions.data.perf.test_performance import PerformanceTestCase, a_deca
from ion_functions.data.ctd_functions import (ctd_sbe52mp_tempwat,
                                              ctd_sbe52mp_condwat,
                                              ctd_sbe52mp_preswat,
                                              ctd_sbe37im_tempwat,
                                              ctd_sbe37im_condwat,
                                              ctd_sbe37im_preswat,
                                              ctd_sbe16plus_tempwat,
                                              ctd_sbe16plus_condwat,
                                              ctd_sbe16plus_preswat,
                                              ctd_sbe16digi_preswat,
                                              ctd_pracsal,
                                              ctd_density)
import numpy as np


class TestCTDPerformance(PerformanceTestCase):
    def setUp(self):
        ### test data for sbe52mp T,C,P calculations (CTDPF series C,K,L)
        self.sbe52mp_T0 = 200000
        self.sbe52mp_C0 = 305000
        self.sbe52mp_P0 = 201000

        ### test data for sbe37im T,C,P calculations (CTDMO all series)
        self.sbe37im_T0 = 340357
        self.sbe37im_C0 = 400000
        self.sbe37im_P0 = 2789
        self.sbe37im_Prange_psia = 1000

        ### test data for sbe16plus_tempwat calculations (CTDPF series A,B; CTDBP all series)
        self.sbe16plus_T0 = 248471
        self.sbe16plus_a0 = 1.281651e-3
        self.sbe16plus_a1 = 2.706002e-4
        self.sbe16plus_a2 = -1.027561e-6
        self.sbe16plus_a3 = 1.749446e-7

        ### test data for sbe16plus_condwat calculations (CTDPF series A,B; CTDBP all series)
        self.sbe16plus_C0 = 1673175
        self.sbe16plus_t1 = 22.544681
        self.sbe16plus_p1 = 27.282116
        self.sbe16plus_g = -9.72193700e-001
        self.sbe16plus_h = 1.38675900e-001
        self.sbe16plus_i = -1.08398500e-004
        self.sbe16plus_j = 2.63219300e-005
        self.sbe16plus_cpcor = -9.57000000e-008
        self.sbe16plus_ctcor = 3.2500e-006

        ### test data for sbe16plus_preswat calculations (CTDPF series A,B; CTDBP series C,D,E,F)
        self.sbe16plus_P0_strain = 528418
        self.sbe16plus_t0_strain = 24303
        self.sbe16plus_ptempa0 = -6.87701000e+001
        self.sbe16plus_ptempa1 = 5.05406200e+001
        self.sbe16plus_ptempa2 = -2.15672900e-001
        self.sbe16plus_ptca0 = 5.24965500e+005
        self.sbe16plus_ptca1 = 7.23620100e+000
        self.sbe16plus_ptca2 = -9.94485900e-002
        self.sbe16plus_ptcb0 = 2.51220000e+001
        self.sbe16plus_ptcb1 = -2.00000000e-004
        self.sbe16plus_ptcb2 = 0.00000000e+000
        self.sbe16plus_pa0 = 1.73472300e+000
        self.sbe16plus_pa1 = 1.57475000e-002
        self.sbe16plus_pa2 = -6.51927800e-010

        ### test data for sbe16digi_preswat calculations (CTDBP series N,O)
        self.sbe16plus_P0_digi = 8833629
        self.sbe16plus_t0_digi = 34336
        self.sbe16plus_C1 = 991.3651
        self.sbe16plus_C2 = 1.01360e-05
        self.sbe16plus_C3 = -1.18210e-04
        self.sbe16plus_D1 = 0.031072
        self.sbe16plus_D2 = 0.0
        self.sbe16plus_T1 = 27.67412
        self.sbe16plus_T2 = -1.08033e-04
        self.sbe16plus_T3 = 1.03670e-06
        self.sbe16plus_T4 = 1.68749e-09
        self.sbe16plus_T5 = 0.0

        ### test data for salinity (C,T,P) and density (S,T,P,lat,lon); ALL CTD
        self.C_L1 = 3.463402
        self.T_L1 = 6.0
        self.P_L1 = 800.0
        self.S_L2 = 34.9
        self.lat = 15.00
        self.lon = -55.00

    def test_ctd_sbe52mp_tempwat(self):
        stats = []
        # create 10000 data points
        data = np.zeros(a_deca, dtype='int16') + self.sbe52mp_T0
        self.profile(stats, ctd_sbe52mp_tempwat, data)

    def test_ctd_sbe52mp_condwat(self):
        stats = []
        # create 10000 data points
        data = np.zeros(a_deca, dtype='int16') + self.sbe52mp_C0
        self.profile(stats, ctd_sbe52mp_condwat, data)

    def test_ctd_sbe52mp_preswat(self):
        stats = []
        # create 10000 data points
        data = np.zeros(a_deca, dtype='int16') + self.sbe52mp_P0
        self.profile(stats, ctd_sbe52mp_preswat, data)

    def test_ctd_sbe37im_tempwat(self):
        stats = []
        # create 10000 data points
        data = np.zeros(a_deca, dtype='int16') + self.sbe37im_T0
        self.profile(stats, ctd_sbe37im_tempwat, data)

    def test_ctd_sbe37im_condwat(self):
        stats = []
        # create 10000 data points
        data = np.zeros(a_deca, dtype='int16') + self.sbe37im_C0
        self.profile(stats, ctd_sbe37im_condwat, data)

    def test_ctd_sbe37im_preswat(self):
        stats = []
        # create 10000 data points
        data = np.zeros(a_deca, dtype='int16') + self.sbe37im_P0
        Prange = np.zeros(a_deca, dtype='int16') + self.sbe37im_Prange_psia
        self.profile(stats, ctd_sbe37im_preswat, data, Prange)

    def test_ctd_sbe16plus_tempwat(self):
        stats = []
        # create 10000 data points
        data = np.zeros(a_deca, dtype='int16') + self.sbe16plus_T0
        a0 = np.zeros(a_deca) + self.sbe16plus_a0
        a1 = np.zeros(a_deca) + self.sbe16plus_a1
        a2 = np.zeros(a_deca) + self.sbe16plus_a2
        a3 = np.zeros(a_deca) + self.sbe16plus_a3
        self.profile(stats, ctd_sbe16plus_tempwat, data, a0, a1, a2, a3)

    def test_ctd_sbe16plus_condwat(self):
        stats = []
        # create 10000 data points
        data = np.zeros(a_deca, dtype='int16') + self.sbe16plus_C0
        t1 = np.zeros(a_deca) + self.sbe16plus_t1
        p1 = np.zeros(a_deca) + self.sbe16plus_p1
        g = np.zeros(a_deca) + self.sbe16plus_g
        h = np.zeros(a_deca) + self.sbe16plus_h
        i = np.zeros(a_deca) + self.sbe16plus_i
        j = np.zeros(a_deca) + self.sbe16plus_j
        cpcor = np.zeros(a_deca) + self.sbe16plus_cpcor
        ctcor = np.zeros(a_deca) + self.sbe16plus_ctcor
        self.profile(stats, ctd_sbe16plus_condwat, data, t1, p1, g, h, i, j, cpcor, ctcor)

    def test_ctd_sbe16plus_preswat(self):
        stats = []
        # create 10000 data points
        P0 = np.zeros(a_deca, dtype='int16') + self.sbe16plus_P0_strain
        t0 = np.zeros(a_deca, dtype='int16') + self.sbe16plus_t0_strain
        ptempa0 = np.zeros(a_deca) + self.sbe16plus_ptempa0
        ptempa1 = np.zeros(a_deca) + self.sbe16plus_ptempa1
        ptempa2 = np.zeros(a_deca) + self.sbe16plus_ptempa2
        ptca0 = np.zeros(a_deca) + self.sbe16plus_ptca0
        ptca1 = np.zeros(a_deca) + self.sbe16plus_ptca1
        ptca2 = np.zeros(a_deca) + self.sbe16plus_ptca2
        ptcb0 = np.zeros(a_deca) + self.sbe16plus_ptcb0
        ptcb1 = np.zeros(a_deca) + self.sbe16plus_ptcb1
        ptcb2 = np.zeros(a_deca) + self.sbe16plus_ptcb2
        pa0 = np.zeros(a_deca) + self.sbe16plus_pa0
        pa1 = np.zeros(a_deca) + self.sbe16plus_pa1
        pa2 = np.zeros(a_deca) + self.sbe16plus_pa2
        self.profile(stats, ctd_sbe16plus_preswat, P0, t0, ptempa0, ptempa1, ptempa2, ptca0, ptca1, ptca2,
                     ptcb0, ptcb1, ptcb0, pa0, pa1, pa2)

    def test_ctd_sbe16digi_preswat(self):
        stats = []
        # create 10000 data points
        P0 = np.zeros(a_deca, dtype='int16') + self.sbe16plus_P0_digi
        t0 = np.zeros(a_deca, dtype='int16') + self.sbe16plus_t0_digi
        C1 = np.zeros(a_deca) + self.sbe16plus_C1
        C2 = np.zeros(a_deca) + self.sbe16plus_C2
        C3 = np.zeros(a_deca) + self.sbe16plus_C3
        D1 = np.zeros(a_deca) + self.sbe16plus_D1
        D2 = np.zeros(a_deca) + self.sbe16plus_D2
        T1 = np.zeros(a_deca) + self.sbe16plus_T1
        T2 = np.zeros(a_deca) + self.sbe16plus_T2
        T3 = np.zeros(a_deca) + self.sbe16plus_T3
        T4 = np.zeros(a_deca) + self.sbe16plus_T4
        T5 = np.zeros(a_deca) + self.sbe16plus_T5
        self.profile(stats, ctd_sbe16digi_preswat, P0, t0, C1, C2, C3, D1, D2, T1, T2, T3, T4, T5)

    def test_ctd_pracsal(self):
        stats = []
        # create 10000 data points
        C = np.zeros(a_deca) + self.C_L1
        T = np.zeros(a_deca) + self.T_L1
        P = np.zeros(a_deca) + self.P_L1
        self.profile(stats, ctd_pracsal, C, T, P)

    def test_ctd_density(self):
        stats = []
        # create 10000 data points
        S = np.zeros(a_deca) + self.S_L2
        T = np.zeros(a_deca) + self.T_L1
        P = np.zeros(a_deca) + self.P_L1
        lat = np.zeros(a_deca) + self.lat
        lon = np.zeros(a_deca) + self.lon
        self.profile(stats, ctd_density, S, T, P, lat, lon)

 