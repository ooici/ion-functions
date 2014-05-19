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
class TestADCPPerformance(PerformanceTestCase):

    def setUp(self):
        # set test inputs
        # setup the test values
        self.gain = 6.
        self.wav = np.array([-2.40000, -0.31200, 0.01110, 0.00442])

    def test_flo_bback_total(self):
        stats = []

        gain = np.repeat(self.gain, 10000000)
        wav = np.tile(self.wav, (10000000, 1))

        self.profile(stats, hy.hyd_acoustic_pwaves, wav, gain)
