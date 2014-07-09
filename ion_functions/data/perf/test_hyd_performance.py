"""
@package ion_functions.data.perf.test_hyd_performance
@file ion_functions/data/perf/test_hyd_performance.py
@author Christopher Wingard
@brief Performance tests for hyd_functions module
"""

import numpy as np
from nose.plugins.attrib import attr

from ion_functions.data.perf.test_performance import PerformanceTestCase
from ion_functions.data import hyd_functions as hy


@attr('PERF', group='func')
class TestHYDPerformance(PerformanceTestCase):

    def setUp(self):
        # setup the test values for HYDAPBB
        self.gain = 6.
        self.wav = np.array([-2.40000, -0.31200, 0.01110, 0.00442])

        # setup the test values for HYDAPLF
        self.raw = np.array([0, 1024, 2048, 3072, 4096])

    def test_hyd_bb_acoustic_pwaves(self):
        stats = []

        gain = np.repeat(self.gain, 10000000)
        wav = np.tile(self.wav, (10000000, 1))
        self.profile(stats, hy.hyd_bb_acoustic_pwaves, wav, gain)

    def test_hyd_lf_acoustic_pwaves(self):
        stats = []

        raw = np.tile(self.raw, (10000000, 1))
        self.profile(stats, hy.hyd_lf_acoustic_pwaves, raw)
