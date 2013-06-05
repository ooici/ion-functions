from ion_functions.data.perf.test_performance import PerformanceTestCase, a_year, a_day
from ion_functions.qc.qc_functions import dataqc_globalrangetest as grt
from ion_functions.qc.qc_functions import dataqc_spiketest as spiketest
from ion_functions.qc.qc_functions import dataqc_stuckvaluetest as stuckvalue
from ion_functions.qc.qc_functions import dataqc_polytrendtest as trend
from ion_functions.qc.qc_functions import dataqc_gradienttest as grad

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

    def test_spiketest(self):
        stats = []

        sample_set = np.empty(a_year, dtype=np.float)
        sample_set.fill(3)
        indexes = [i for i in xrange(a_day * 2) if not i%20]
        sample_set[indexes] = 40

        self.profile(stats, spiketest, sample_set, 0.1)


    def test_stuckvalue(self):
        stats = []
        
        sample_set = np.arange(a_year, dtype=np.float)
        v = [4.83, 1.40, 3.33, 3.33, 3.33, 3.33, 4.09, 2.97, 2.85, 3.67]
        sample_set[0:len(v)] = v
        self.profile(stats, stuckvalue, sample_set, 0.001, 4)

    def test_trend(self):
        stats = []
        x = np.arange(a_year, dtype=np.float)
        sample_set = np.sin(np.pi * 2 * x/60.) * 6 + 3.
        self.profile(stats, trend, sample_set, np.arange(a_year, dtype=np.float))

    def test_gradient(self):
        stats = []
        sample_set = np.arange(a_day * 2, dtype=np.float)

        self.profile(stats, grad, sample_set, sample_set, [-50,50], .1, [], 5)



        
