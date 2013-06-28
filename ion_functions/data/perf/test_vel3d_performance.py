#!/usr/bin/env python
from ion_functions.data.perf.test_performance import PerformanceTestCase, a_year, a_day
from ion_functions.data.vel_functions import nobska_mag_corr_east, nobska_mag_corr_north, magnetic_declination, adcp_magvar

#def nobska_mag_corr_east(uu,vv,lat,lon,timestamp,z=0):
#def nobska_mag_corr_north(uu,vv,lat,lon,timestamp,z=0):
import numpy as np

class TestVel3DBPerformance(PerformanceTestCase):
    def setUp(self):
        self.lat = 14.6846
        self.lon = -51.044
        self.ts = np.ones(a_day*2,dtype=np.int) * 3319563600
        self.ve = np.ones(a_day*2, dtype=np.float) * -3.2
        self.vn = np.ones(a_day*2, dtype=np.float) * 18.2
        self.vu = np.ones(a_day*2, dtype=np.float) * -1.1
    
    def test_mag_corr_east(self):
        stats = []

        self.profile(stats, nobska_mag_corr_east, self.ve, self.vn, self.lat, self.lon, self.ts, 6)

    def test_mag_decl(self):
        stats = []
        self.profile(stats, magnetic_declination, self.lat, self.lon, self.ts, 6)

    def test_magvar(self):
        stats = []
        theta = magnetic_declination(self.lat, self.lon, self.ts, 6)
        magvar = np.vectorize(adcp_magvar)
        self.profile(stats, magvar, theta, self.ve, self.vn)

