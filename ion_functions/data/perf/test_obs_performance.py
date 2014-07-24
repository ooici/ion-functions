"""
@package ion_functions.data.perf.test_obs_performance
@file ion_functions/data/perf/test_obs_performance.py
@author Christopher Wingard
@brief Performance tests for obs_functions module
"""

import numpy as np
from nose.plugins.attrib import attr

from ion_functions.data.perf.test_performance import PerformanceTestCase
from ion_functions.data import obs_functions as obs


@attr('PERF', group='func')
class TestOBSPerformance(PerformanceTestCase):

    def setUp(self):
        # setup the test values for all the obs tests
        self.raw = np.array([0, 512, 1024, 2048, 3072, 4096])

    def test_obs_bb_ground_velocity(self):
        stats = []

        raw = np.tile(self.raw, (10000000, 1))
        self.profile(stats, obs.obs_bb_ground_velocity, raw)

    def test_obs_bb_ground_acceleration(self):
        stats = []

        raw = np.tile(self.raw, (10000000, 1))
        self.profile(stats, obs.obs_bb_ground_acceleration, raw)

    def test_obs_sp_ground_velocity(self):
        stats = []

        raw = np.tile(self.raw, (10000000, 1))
        self.profile(stats, obs.obs_sp_ground_velocity, raw)
