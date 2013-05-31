from ion_functions.data.perf.test_performance import PerformanceTestCase, a_year
from ion_functions.qc.qc_functions import dataqc_globalrangetest as grt
from ion_functions.qc.qc_functions import dataqc_spiketest as spiketest

import numpy as np
import unittest


class TestQCPerformance(PerformanceTestCase):
    def test_globalrangetest(self):
        stats = []
        sample_set = np.empty(a_year, dtype=np.float)
        sample_set.fill(17.)
        indexes = [i for i in xrange(a_year) if not i%20]
        sample_set[indexes] = 40

        self.profile(stats, grt, sample_set, [10,20])

    @unittest.skip('spike takes too long to profile')
    def test_spiketest(self):
        stats = []

        sample_set = np.empty(a_year, dtype=np.float)
        sample_set.fill(3)
        indexes = [i for i in xrange(a_year) if not i%20]
        sample_set[indexes] = 40

        print "Profiling..."
        self.profile(stats, spiketest, sample_set, 0.1)


