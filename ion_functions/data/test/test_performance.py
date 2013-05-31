#!/usr/bin/env python
'''
@author Luke Campbell <luke.s.campbell at gmail dot-ish-thing com>
@file ion_functions/data/test/test_performance.py
@brief Profiling of Data Functions
'''

from nose.plugins.attrib import attr
from unittest import TestCase

from ion_functions.data.co2_functions import pco2_thermistor

import time
import numpy as np

class TimeIt(object):
    def __init__(self, results=[]):
        self.results = results
    def __enter__(self):
        self.i = time.time()
    def __exit__(self, type, value, traceback):
        v = time.time() - self.i
        self.results.append(v)

@attr('PERF')
class TestPerformance(TestCase):

    def profile(self, stats, func, *args, **kwargs):
        for i in xrange(10):
            with TimeIt(stats):
                func(*args, **kwargs)
            print 'Run %i: %s' %(i,stats[i])

        print 'Total: ', np.average(np.array(stats))


    def test_pco2_thermistor(self):

        stats = []
        s = '*7E2704CBACF1230081007C0B2B00BF080800DB01BF0390007C00790B3000C7080B00EA0C5606C80A'
        traw = int(s[75:79], 16)

        sample_data = np.empty(3600 * 24 * 365, dtype='int32')
        sample_data.fill(traw)
        self.profile(stats,pco2_thermistor,sample_data)


