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

        print 'Mean: ', np.average(np.array(stats))

