#!/usr/bin/env python
"""
@package ion_functions.test.obs_functions
@file ion_functions/test/test_obs_functions.py
@author Christopher Wingard
@brief Unit tests for obs_functions module
"""
from nose.plugins.attrib import attr
from ion_functions.test.base_test import BaseUnitTestCase

import numpy as np
from ion_functions.data import obs_functions as obs


@attr('UNIT', group='func')
class TestOBSFunctionsUnit(BaseUnitTestCase):

    def setUp(self):
        # setup the test values
        self.raw = np.array([0, 512, 1024, 2048, 3072, 4096])
        self.grndvel = np.array([0.00000E+00, 5.46133E-07, 1.09227E-06,
                                 2.18453E-06, 3.27680E-06, 4.36907E-06])
        self.grndacc = np.array([0.00000E+00, 1.61260E-03, 3.22520E-03,
                                 6.45039E-03, 9.67559E-03, 1.29008E-02])
        self.sgrdvel = np.array([0.00000E+00, 6.05867E-07, 1.21173E-06,
                                 2.42347E-06, 3.63520E-06, 4.84693E-06])

    def test_obs_bb_ground_velocity(self):
        """
        Test obs_bb_ground_velocity function.

        Values based on a test data set created using Excel by the test author
        (DPS does not provide any test data).

        Implemented by Christopher Wingard, May 2014
        """
        # calculate broadband ground velocity
        out = obs.obs_bb_ground_velocity(self.raw)

        # How'd we do?
        np.testing.assert_allclose(out, self.grndvel, rtol=1e-10, atol=1e-10)

    def test_obs_bb_ground_acceleration(self):
        """
        Test obs_bb_ground_acceleration function.

        Values based on a test data set created using Excel by the test author
        (DPS does not provide any test data).

        Implemented by Christopher Wingard, May 2014
        """
        # calculate broadband ground acceleration
        out = obs.obs_bb_ground_acceleration(self.raw)

        # How'd we do?
        np.testing.assert_allclose(out, self.grndacc, rtol=1e-6, atol=0)

    def test_obs_sp_ground_velocity(self):
        """
        Test obs_sp_ground_velocity function.

        Values based on a test data set created using Excel by the test author
        (DPS does not provide any test data).

        Implemented by Christopher Wingard, May 2014
        """
        # calculate short period ground velocity
        out = obs.obs_sp_ground_velocity(self.raw)

        # How'd we do?
        np.testing.assert_allclose(out, self.sgrdvel, rtol=1e-10, atol=1e-10)
