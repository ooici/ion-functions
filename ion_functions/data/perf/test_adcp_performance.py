"""
@package ion_functions.data.perf.test_adcp_performance
@file ion_functions/data/perf/test_adcp_performance.py
@author Christopher Wingard
@brief Performance tests for adcp_functions module
"""

import numpy as np
from nose.plugins.attrib import attr

from ion_functions.data.perf.test_performance import PerformanceTestCase
from ion_functions.data import adcp_functions as af

# Note, the VADCP related data products use the same internal functions as the
# family of beam wrapper functions (e.g. adcp_beam_eastward). Thus, those
# functions won't be added to this test. The only way to really speed this
# process up any further is to set the wrapper functions to return all the data
# products for an instrument at once rather than singly. That way the functions
# only need to be run once rather than 4 times for each instrument (a factor of
# four improvement in performance).


@attr('PERF', group='func')
class TestADCPPerformance(PerformanceTestCase):

    def setUp(self):
        # set test inputs -- values from DPS
        self.b1 = np.array([-0.0300, -0.2950, -0.5140, -0.2340, -0.1880,
                            0.2030, -0.3250,  0.3050, -0.2040, -0.2940]) * 1000
        self.b2 = np.array([0.1800, -0.1320,  0.2130,  0.3090,  0.2910,
                            0.0490,  0.1880,  0.3730, -0.0020, 0.1720]) * 1000
        self.b3 = np.array([-0.3980, -0.4360, -0.1310, -0.4730, -0.4430,
                            0.1880, -0.1680,  0.2910, -0.1790, 0.0080]) * 1000
        self.b4 = np.array([-0.2160, -0.6050, -0.0920, -0.0580,  0.4840,
                            -0.0050,  0.3380,  0.1750, -0.0800, -0.5490]) * 1000
        self.echo = np.array([0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250])
        self.sfactor = 0.45
        self.heading = 98.4100 / 100.
        self.pitch = 0.6900 / 100.
        self.roll = -2.5400 / 100.
        self.orient = 1
        self.lat = 50.0000
        self.lon = -145.0000
        self.depth = 0.0
        self.ntp = 3545769600.0    # May 12, 2012

        # set expected results -- velocity profiles in earth coordinates
        # (values in DPS)
        self.uu = np.array([0.2175, -0.2814, -0.1002, 0.4831, 1.2380,
                            -0.2455, 0.6218, -0.1807, 0.0992, -0.9063])
        self.vv = np.array([-0.3367, -0.1815, -1.0522, -0.8676, -0.8919,
                            0.2585, -0.8497, -0.0873, -0.3073, -0.5461])
        self.ww = np.array([0.1401,  0.3977,  0.1870,  0.1637,  0.0091,
                            -0.1290,  0.0334, -0.3017, 0.1384, 0.1966])

    def test_adcp_backscatter(self):
        stats = []

        echo = np.tile(self.echo, (10000, 1))
        sfactor = np.repeat(self.sfactor, 10000)
        self.profile(stats, af.adcp_backscatter, echo, sfactor)

    def test_adcp_beam_eastward(self):
        stats = []

        b1 = np.tile(self.b1, (10000, 1))
        b2 = np.tile(self.b2, (10000, 1))
        b3 = np.tile(self.b3, (10000, 1))
        b4 = np.tile(self.b4, (10000, 1))

        h = np.repeat(self.heading, 10000)
        p = np.repeat(self.pitch, 10000)
        r = np.repeat(self.roll, 10000)
        vf = np.repeat(self.orient, 10000)
        lat = np.repeat(self.lat, 10000)
        lon = np.repeat(self.lon, 10000)
        z = np.repeat(self.depth, 10000)
        dt = np.repeat(self.ntp, 10000)

        self.profile(stats, af.adcp_beam_eastward, b1, b2, b3, b4, h, p, r, vf, lat, lon, z, dt)

    def test_adcp_beam_northward(self):
        stats = []

        b1 = np.tile(self.b1, (10000, 1))
        b2 = np.tile(self.b2, (10000, 1))
        b3 = np.tile(self.b3, (10000, 1))
        b4 = np.tile(self.b4, (10000, 1))

        h = np.repeat(self.heading, 10000)
        p = np.repeat(self.pitch, 10000)
        r = np.repeat(self.roll, 10000)
        vf = np.repeat(self.orient, 10000)
        lat = np.repeat(self.lat, 10000)
        lon = np.repeat(self.lon, 10000)
        z = np.repeat(self.depth, 10000)
        dt = np.repeat(self.ntp, 10000)

        self.profile(stats, af.adcp_beam_northward, b1, b2, b3, b4, h, p, r, vf, lat, lon, z, dt)

    def test_adcp_beam_vertical(self):
        stats = []

        b1 = np.tile(self.b1, (10000, 1))
        b2 = np.tile(self.b2, (10000, 1))
        b3 = np.tile(self.b3, (10000, 1))
        b4 = np.tile(self.b4, (10000, 1))

        h = np.repeat(self.heading, 10000)
        p = np.repeat(self.pitch, 10000)
        r = np.repeat(self.roll, 10000)
        vf = np.repeat(self.orient, 10000)

        self.profile(stats, af.adcp_beam_vertical, b1, b2, b3, b4, h, p, r, vf)

    def test_adcp_beam_error(self):
        stats = []

        b1 = np.tile(self.b1, (10000, 1))
        b2 = np.tile(self.b2, (10000, 1))
        b3 = np.tile(self.b3, (10000, 1))
        b4 = np.tile(self.b4, (10000, 1))

        self.profile(stats, af.adcp_beam_error, b1, b2, b3, b4)

    def test_adcp_earth_eastward(self):
        stats = []

        u = np.tile(self.uu, (10000, 1))
        v = np.tile(self.vv, (10000, 1))

        lat = np.repeat(self.lat, 10000)
        lon = np.repeat(self.lon, 10000)
        z = np.repeat(self.depth, 10000)
        dt = np.repeat(self.ntp, 10000)

        self.profile(stats, af.adcp_earth_eastward, u, v, z, lat, lon, dt)

    def test_adcp_earth_northward(self):
        stats = []

        u = np.tile(self.uu, (10000, 1))
        v = np.tile(self.vv, (10000, 1))

        lat = np.repeat(self.lat, 10000)
        lon = np.repeat(self.lon, 10000)
        z = np.repeat(self.depth, 10000)
        dt = np.repeat(self.ntp, 10000)

        self.profile(stats, af.adcp_earth_northward, u, v, z, lat, lon, dt)

    def test_adcp_earth_vertical(self):
        stats = []

        w = np.tile(self.ww, (10000, 1))

        self.profile(stats, af.adcp_earth_vertical, w)
        # adcp_earth_error is the same transform, so this test applies to both
