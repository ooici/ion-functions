#!/usr/bin/env python
'''
@author Luke Campbell <luke.s.campbell at gmail dot-ish-thing com>
@file ion_functions/data/test/test_performance.py
@brief Profiling of Data Functions
'''

from nose.plugins.attrib import attr
from unittest import TestCase


import time
import numpy as np

a_year = 24 * 3600 * 365
a_week = 24 * 3600 * 7
a_month = 24 * 3600 * 30
a_day = 24 * 3600

class TimeIt(object):
    def __init__(self, results=[]):
        self.results = results
    def __enter__(self):
        self.i = time.time()
    def __exit__(self, type, value, traceback):
        v = time.time() - self.i
        self.results.append(v)

@attr('PERF')
class PerformanceTestCase(TestCase):
    def profile(self, stats, func, *args, **kwargs):
        for i in xrange(10):
            with TimeIt(stats):
                func(*args, **kwargs)
            print 'Run %i: %s' %(i,stats[i])
            if stats[i] > 10:
                raise AssertionError('Performance standard failed. Method exceeded 20 seconds')

        print 'Mean: ', np.average(np.array(stats))

