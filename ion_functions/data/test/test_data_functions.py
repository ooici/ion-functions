#!/usr/bin/env python

"""
@package ion_functions.test.data_functions
@file ion_functions/test/data_functions.py
@author Christopher Mueller
@brief Unit tests for data_functions module
"""

from nose.plugins.attrib import attr
from ion_functions.test.base_test import BaseUnitTestCase

import numpy as np
from ion_functions.data import data_functions as dfunc


@attr('UNIT', group='func')
class TestDataFunctionsUnit(BaseUnitTestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_data_density(self):
        """
        Test data_density function.

        Values based on those in pygsw.test.test_pygsw.py
        """

        from pygsw.test import test_pygsw

        arrlen = 100

        sp = np.array([test_pygsw.sp] * arrlen)
        p = np.array([test_pygsw.p] * arrlen)
        t = np.array([test_pygsw.t] * arrlen)
        lat = np.array([test_pygsw.lat] * arrlen)
        lon = np.array([test_pygsw.lon] * arrlen)

        # This value slightly different from what's expected in test_pygsw (1026.4562376198473e0)
        # because of rounding of the intermediate values in test_pygsw
        out = np.array([1027.69107725] * arrlen)

        self.assertTrue(np.allclose(out, dfunc.data_density(sp, p, t, lat, lon)))