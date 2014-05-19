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
from ion_functions.data import hyd_functions as hydfunc


@attr('UNIT', group='func')
class TestHYDFunctionsUnit(BaseUnitTestCase):

    def setUp(self):
        # setup the test values
        self.gain = 6.
        self.wav = np.array([-2.40000, -0.31200, 0.01110, 0.00442])
        self.tsv = np.atleast_2d(np.array([-3.60855, -0.46911, 0.01669, 0.00665]))

    def test_hyd_acoustic_pwaves(self):
        """
        Test hyd_acoustic_pwaves function.

        Values based on those described in DPS as available on Alfresco:

        OOI (2013). Data Product Specification for Acoustic Pressure Waves.
            Document Control Number 1341-00820.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00820_Data_Product_SPEC_HYDAPBB_OOI.pdf)

        Note, once again, the test data provided in the DPS is wrong. Had to
        reset tsv based on running the oooohhhh soooo complex algorithm through
        matlab.

        Implemented by Christopher Wingard, May 2014
        """
        # calculate the corrected compass direction
        out = hydfunc.hyd_acoustic_pwaves(self.wav, self.gain)

        # How'd we do?
        np.testing.assert_allclose(out, self.tsv, rtol=1e-5, atol=1e-5)
