#!/usr/bin/env python

from ion_functions.data.perf.test_performance import PerformanceTestCase, a_day
from ion_functions.data.vel_functions import nobska_mag_corr_east, nobska_mag_corr_north
from ion_functions.data.vel_functions import nortek_mag_corr_east, nortek_mag_corr_north
from ion_functions.data.generic_functions import magnetic_declination, magnetic_correction

import numpy as np


class TestVel3DBPerformance(PerformanceTestCase):
    def setUp(self):
        self.lat = 14.6846
        self.lon = -51.044
        self.ts = np.ones(a_day*2, dtype=np.int) * 3319563600
        self.ve = np.ones(a_day*2, dtype=np.float) * -3.2
        self.vn = np.ones(a_day*2, dtype=np.float) * 18.2
        self.vu = np.ones(a_day*2, dtype=np.float) * -1.1

    def test_nobska_corr_east(self):
        stats = []
        self.profile(stats, nobska_mag_corr_east, self.ve, self.vn, self.lat, self.lon, self.ts, 3)

    def test_nobska_corr_north(self):
        stats = []
        self.profile(stats, nobska_mag_corr_north, self.ve, self.vn, self.lat, self.lon, self.ts, 3)

    def test_nortek_corr_east(self):
        stats = []
        self.profile(stats, nortek_mag_corr_east, self.ve, self.vn, self.lat, self.lon, self.ts, 3)

    def test_nortek_corr_north(self):
        stats = []
        self.profile(stats, nortek_mag_corr_north, self.ve, self.vn, self.lat, self.lon, self.ts, 3)

    def test_magnetic_declination(self):
        stats = []
        self.profile(stats, magnetic_declination, self.lat, self.lon, self.ts, 3)

    def test_magnetic_correction(self):
        stats = []
        theta = magnetic_declination(self.lat, self.lon, self.ts, 3)
        magcor = np.vectorize(magnetic_correction)
        self.profile(stats, magcor, theta, self.ve, self.vn)
