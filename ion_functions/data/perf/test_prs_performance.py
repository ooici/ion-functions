#!/usr/bin/env python

"""
@package ion_functions.perf.test_prs_performance
@file ion_functions/perf/test_prs_performance.py
@author Russell Desiderio
@brief Performance tests for prs_functions module
"""

from nose.plugins.attrib import attr
from ion_functions.data.perf.test_performance import PerformanceTestCase
from ion_functions.data import prs_functions as pf
import numpy as np


@attr('UNIT', group='func')
class TestPRSPerformance(PerformanceTestCase):

    # test 10000 data packets; test data will be 10 data packets, so
    repeat_value = 1000

    # setup up a list with the test data used by the three test functions
    #      xtilt ,  ytilt , scmp, sernum , ccmp, tmag, angle, tdir
    lily = [
        [-90.000, 50.000, 0.15, 'N9676', 200, 103.0, -29.1, 139],
        [0.000, 0.000, 359.76, 'N9676', 200, 0.0, 0.0, 290],
        [10.000, -30.000, 144.21, 'N9656', 87, 31.6, -71.6, 249],
        [40.000, -60.000, 359.76, 'N9656', 166, 72.1, -56.3, 312],
        [50.000, -70.000, 0.15, 'N9652', 173, 86.0, -54.5, 317],
        [100.000, -120.000, 359.76, 'N9652', 173, 156.2, -50.2, 313],
        [110.000, -130.000, 0.15, 'N9655', 173, 170.3, -49.8, 313],
        [160.000, 120.000, 359.76, 'N9655', 173, 200.0, 36.9, 226],
        [-60.000, -220.000, 71.99, 'N9651', 170, 228.0, 74.7, 5],
        [-310.000, -10.000, 359.76, 'N9651', 183, 310.2, 1.8, 91]
    ]

    def test_prs_bottilt_ccmp(self):
        stats = []

        # create 10000 data points
        scmp = np.atleast_1d([row[2] for row in self.lily])
        scmp = np.tile(scmp, (self.repeat_value))
        snum = np.atleast_1d([row[3] for row in self.lily])
        snum = np.tile(snum, (self.repeat_value))

        # time it
        self.profile(stats, pf.prs_bottilt_ccmp, scmp, snum)

    def test_prs_bottilt_tmag(self):
        stats = []

        # create 10000 data points
        xtilt = np.atleast_1d([row[0] for row in self.lily])
        xtilt = np.tile(xtilt, (self.repeat_value))
        ytilt = np.atleast_1d([row[1] for row in self.lily])
        ytilt = np.tile(ytilt, (self.repeat_value))

        # time it
        self.profile(stats, pf.prs_bottilt_tmag, xtilt, ytilt)

    def test_prs_bottilt_tdir(self):
        stats = []

        # create 10000 data points
        xtilt = np.atleast_1d([row[0] for row in self.lily])
        xtilt = np.tile(xtilt, (self.repeat_value))
        ytilt = np.atleast_1d([row[1] for row in self.lily])
        ytilt = np.tile(ytilt, (self.repeat_value))
        ccmp = np.atleast_1d([row[4] for row in self.lily])
        ccmp = np.tile(ccmp, (self.repeat_value))

        # time it
        self.profile(stats, pf.prs_bottilt_tdir, xtilt, ytilt, ccmp)
