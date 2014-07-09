#!/usr/bin/env python
"""
@package ion_functions.test.hyd_functions
@file ion_functions/test/test_hyd_functions.py
@author Christopher Wingard
@brief Unit tests for hyd_functions module
"""
from nose.plugins.attrib import attr
from ion_functions.test.base_test import BaseUnitTestCase

import numpy as np
from ion_functions.data import obs_functions as obs


@attr('UNIT', group='func')
class TestHYDFunctionsUnit(BaseUnitTestCase):

    def setUp(self):
        # setup the test values
        self.raw = np.array([0, 1024, 2048, 3072, 4096])
        self.grndvel = np.array([])
        self.grndacc = np.array([])
        self.sgrdvel = np.array([])

    def test_obs_bb_ground_velocity(self):
        """
        Test obs_acoustic_pwaves function.

        Values based on those described in DPS as available on Alfresco:

        OOI (2013). Data Product Specification for Acoustic Pressure Waves.
            Document Control Number 1341-00820.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00820_Data_Product_SPEC_HYDAPBB_OOI.pdf)

        Implemented by Christopher Wingard, May 2014
        """
        # calculate broadband ground velocity
        out = obs.obs_bb_ground_velocity(self.raw)

        # How'd we do?
        np.testing.assert_allclose(out, self.grndvel, rtol=1e-6, atol=1e-6)

    def test_obs_bb_ground_acceleration(self):
        """
        Test obs_acoustic_pwaves function.

        Values based on those described in DPS as available on Alfresco:

        OOI (2013). Data Product Specification for Acoustic Pressure Waves.
            Document Control Number 1341-00820.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00820_Data_Product_SPEC_HYDAPBB_OOI.pdf)

        Implemented by Christopher Wingard, May 2014
        """
        # calculate broadband ground acceleration
        out = obs.obs_bb_ground_acceleration(self.raw)

        # How'd we do?
        np.testing.assert_allclose(out, self.grndacc, rtol=1e-6, atol=1e-6)

    def test_obs_sp_ground_velocity(self):
        """
        Test obs_acoustic_pwaves function.

        Values based on those described in DPS as available on Alfresco:

        OOI (2013). Data Product Specification for Acoustic Pressure Waves.
            Document Control Number 1341-00820.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00820_Data_Product_SPEC_HYDAPBB_OOI.pdf)

        Implemented by Christopher Wingard, May 2014
        """
        # calculate short period ground velocity
        out = obs.obs_sp_ground_velocity(self.raw)

        # How'd we do?
        np.testing.assert_allclose(out, self.sgrdvel, rtol=1e-6, atol=1e-6)
