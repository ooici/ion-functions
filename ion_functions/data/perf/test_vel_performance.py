#!/usr/bin/env python

from ion_functions.data.perf.test_performance import PerformanceTestCase, a_day
from ion_functions.data.vel_functions import nobska_mag_corr_east, nobska_mag_corr_north
from ion_functions.data.vel_functions import nortek_mag_corr_east, nortek_mag_corr_north
from ion_functions.data.vel_functions import vel3dk_mag_corr_east, vel3dk_mag_corr_north
from ion_functions.data.generic_functions import magnetic_declination, magnetic_correction

import numpy as np


class TestVelPerformance(PerformanceTestCase):
    """
    Performance tests for the vel family of functions to compare to the
    optimization benchmark of running the function for 10000 data
    records in under 10 seconds.

    Implemented by:
        2014-03-04: Stuart Pearce. Initial code.
    """
    def setUp(self):
        self.lat = np.ones(10000, dtype=np.float) * 14.6846
        self.lon = np.ones(10000, dtype=np.float) * -51.044
        self.ts = np.ones(10000, dtype=np.int) * 3319563600
        self.ve = np.ones(10000, dtype=np.float) * -3.2
        self.vn = np.ones(10000, dtype=np.float) * 18.2
        #self.vu = np.ones(10000, dtype=np.float) * -1.1
        self.ve_vel3dk = np.ones(10000, dtype=np.int) * 12345
        self.vn_vel3dk = np.ones(10000, dtype=np.int) * 15432
        self.vscale = -4

    def test_nobska_corr_east(self):
        """
        Performance test for the nobska_mag_corr_east function for the
        VEL3D-B (Nobska MAVS-4) instrument data processing algorithm.
        """
        stats = []
        self.profile(
            stats, nobska_mag_corr_east, self.ve, self.vn,
            self.lat, self.lon, self.ts, 3)

    def test_nobska_corr_north(self):
        """
        Performance test for the nobska_mag_corr_north function for the
        VEL3D-B (Nobska MAVS-4) instrument data processing algorithm.
        """
        stats = []
        self.profile(
            stats, nobska_mag_corr_north, self.ve, self.vn,
            self.lat, self.lon, self.ts, 3)

    def test_nortek_corr_east(self):
        """
        Performance test for the nortek_mag_corr_east function for the
        VEL3D-CD & VELPT (Nortek Vector and Aquadopp) instrument data
        processing algorithm.
        """
        stats = []
        self.profile(
            stats, nortek_mag_corr_east, self.ve, self.vn,
            self.lat, self.lon, self.ts, 3)

    def test_nortek_corr_north(self):
        """
        Performance test for the nortek_mag_corr_north function for the
        VEL3D-CD & VELPT (Nortek Vector and Aquadopp) instrument data
        processing algorithm.
        """
        stats = []
        self.profile(
            stats, nortek_mag_corr_north, self.ve, self.vn,
            self.lat, self.lon, self.ts, 3)

    def test_vel3dk_corr_east(self):
        """
        Performance test for the vel3dk_mag_corr_east function for the
        VEL3D-K (Nortek Aquadopp 2) instrument data processing
        algorithm.
        """
        stats = []
        self.profile(
            stats, vel3dk_mag_corr_east, self.ve_vel3dk, self.vn_vel3dk,
            self.lat, self.lon, self.ts, 3, self.vscale)

    def test_vel3dk_corr_north(self):
        """
        Performance test for the vel3dk_mag_corr_north function for the
        VEL3D-K (Nortek Aquadopp 2) instrument data processing
        algorithm.
        """
        stats = []
        self.profile(
            stats, vel3dk_mag_corr_north, self.ve_vel3dk, self.vn_vel3dk,
            self.lat, self.lon, self.ts, 3, self.vscale)

    def test_magnetic_declination(self):
        """
        Performance test for the magnetic_declination function for the
        calculating magnetic declination, used in multiple data product
        algorithms.
        """
        stats = []
        self.profile(
            stats, magnetic_declination, self.lat, self.lon, self.ts, 3)

    def test_magnetic_correction(self):
        """
        Performance test for the magnetic_correction function for the
        correcting vectors for magnetic declination, used in multiple
        data product algorithms.
        """
        stats = []
        theta = magnetic_declination(self.lat, self.lon, self.ts, 3)
        magcor = np.vectorize(magnetic_correction)
        self.profile(stats, magcor, theta, self.ve, self.vn)
