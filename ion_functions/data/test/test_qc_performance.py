from ion_functions.data.test.test_performance import PerformanceTestCase
from ion_functions.qc.qc_functions import dataqc_globalrangetest as grt
import numpy as np


class TestQCPerformance(PerformanceTestCase):
    def test_globalrangetest(self):
        stats = []
        sample_set = np.empty(24 * 3600 * 365, dtype=np.float)
        sample_set.fill(17.)
        indexes = [i for i in xrange(24 * 3600 * 365) if not i%20]
        sample_set[indexes] = 40

        print 'starting profile'
        self.profile(stats, grt, sample_set, [10,20])



